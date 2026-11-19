# shellcheck shell=bash
# Shared code for the R2 release scripts. Source it; do not run it.
#
# The bucket holds, and nothing else writes it:
#   app/LABEL/CLAIM                  reserves the label; written first
#   app/LABEL/dist/...               the built web app (web/dist)
#   app/LABEL/COMMIT                 the commit it was built from
#   app/LABEL/MANIFEST.sha256        "SHA256  ./path" for dist/; written last (the seal)
#   pub/NAME/CLAIM                   reserves the name; written first
#   pub/NAME/report/...              a report tree (grm report output)
#   pub/NAME/SOURCE.sha256           its manifest; written last (the seal)
#   state/active, state/previous     one JSON line each: the pair the Worker
#                                    serves, and the rollback target
#   state/history/TIME--LABEL--NAME  every activation, never rewritten
#
# A label or name is never written twice: publish refuses one whose CLAIM or
# seal exists, so an interrupted upload leaves its claim and the next try
# takes a new name. Only a sealed release can be activated.
#
# Configuration, from the environment only:
#   GRM_SITE_TARGET       local (Miniflare state under GRM_SITE_PERSIST) or
#                         remote (the Cloudflare account); required
#   GRM_SITE_PERSIST      local state directory (local only)
#   GRM_SITE_TOKEN_FILE   remote only: a file of mode 0600 holding the API
#                         token. It is read into wrangler's environment only,
#                         never printed, logged or put on a command line
#   GRM_SITE_ACCOUNT_ID   remote only: the account id (not a secret)
#   GRM_SITE_BUCKET       bucket name (default grm-site)
#   GRM_SITE_LOG          absolute path of the operator log; every step of an
#                         --apply run appends a line (required)
#   GRM_SITE_URL          optional base URL that activate.sh checks after the
#                         switch, rolling back if the check fails
#   GRM_SITE_JOBS         parallel object reads when verifying (default 8, local 4)

# Tracing would print the token as wrangler's environment is set, so it stays off
# even under bash -x.
{ set +x; } 2>/dev/null
set -euo pipefail
trap ':' PIPE

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKER=$(cd -- "$HERE/../worker" && pwd)
TOOLKIT=$(cd -- "$HERE/../../../scripts/deployment" && pwd)
WRANGLER=$WORKER/node_modules/.bin/wrangler
BUCKET=${GRM_SITE_BUCKET:-grm-site}
JOBS=${GRM_SITE_JOBS:-$( [ "${GRM_SITE_TARGET:-}" = local ] && echo 4 || echo 8)}
APPLY=0
ARGS=()
export WRANGLER_SEND_METRICS=false

say() {
    local line
    line="$(date -u +%FT%TZ) [$GRM_SITE_TARGET $BUCKET] $*"
    printf '%s\n' "$line" 2>/dev/null || true
    if ((APPLY)); then printf '%s\n' "$line" >>"$GRM_SITE_LOG"; fi
}
die() {
    say "FAIL: $1" >&2
    exit "${2:-1}"
}
would() { printf 'dry-run: would %s\n' "$*"; }

parse_args() {
    local arg
    for arg in "$@"; do
        case $arg in
            --apply) APPLY=1 ;;
            --rollback) ROLLBACK=1 ;;
            --*) echo "unknown option $arg" >&2; exit 2 ;;
            *) ARGS+=("$arg") ;;
        esac
    done
}

valid_name() { [[ $1 =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] && [[ $1 != *--* ]]; }
check_name() { valid_name "$1" || die "invalid name '$1'" 2; }

case ${GRM_SITE_TARGET:-} in
    local)
        [ -n "${GRM_SITE_PERSIST:-}" ] || die "GRM_SITE_TARGET=local needs GRM_SITE_PERSIST" 2
        WHERE=(--local --persist-to "$GRM_SITE_PERSIST")
        ;;
    remote)
        [ -n "${GRM_SITE_ACCOUNT_ID:-}" ] || die "GRM_SITE_TARGET=remote needs GRM_SITE_ACCOUNT_ID" 2
        tf=${GRM_SITE_TOKEN_FILE:-}
        [ -n "$tf" ] && [ -f "$tf" ] && [ ! -L "$tf" ] || die "GRM_SITE_TARGET=remote needs GRM_SITE_TOKEN_FILE, a regular file" 2
        [ "$(stat -c %a "$tf")" = 600 ] && [ "$(stat -c %u "$tf")" = "$(id -u)" ] ||
            die "$tf must be mode 0600 and owned by $(id -un)" 2
        WHERE=(--remote)
        ;;
    *) echo "set GRM_SITE_TARGET to local or remote" >&2; exit 2 ;;
esac
case ${GRM_SITE_LOG:-} in
    /?*) ;;
    *) echo "set GRM_SITE_LOG to the absolute path of the operator log" >&2; exit 2 ;;
esac
[[ $JOBS =~ ^[1-9][0-9]?$ ]] || die "GRM_SITE_JOBS must be 1 to 99" 2
[ -x "$WRANGLER" ] || die "missing $WRANGLER: run npm ci in $WORKER" 2

# wrangler with the token in its environment only. The token never reaches
# this shell's variables, argv or output.
wr() {
    if [ "$GRM_SITE_TARGET" = remote ]; then
        CLOUDFLARE_ACCOUNT_ID=$GRM_SITE_ACCOUNT_ID CLOUDFLARE_API_TOKEN=$(<"$GRM_SITE_TOKEN_FILE") \
            "$WRANGLER" "$@" "${WHERE[@]}"
    else
        "$WRANGLER" "$@" "${WHERE[@]}"
    fi
}

# get KEY FILE: 0 and FILE written, 1 when KEY does not exist; anything else
# (auth, network) is tried 3 times and then dies, so a failed read is never
# taken for an absent key. Parallel local reads of one Miniflare state can
# fail transiently, which is why it retries.
get() {
    local err rc try
    err=$(mktemp)
    for try in 1 2 3; do
        rc=0
        wr r2 object get "$BUCKET/$1" --file "$2" >/dev/null 2>"$err" || rc=$?
        if ((rc == 0)); then rm -f "$err"; return 0; fi
        if grep -q 'The specified key does not exist' "$err"; then rm -f "$err" "$2"; return 1; fi
        echo "  read $BUCKET/$1 failed (try $try): $(grep -m1 -i 'error' "$err" | cut -c1-200)" >&2
        sleep "$try"
    done
    sed 's/^/  wrangler: /' "$err" >&2
    rm -f "$err"
    die "reading $BUCKET/$1 failed"
}
exists() {
    local tmp rc=0
    tmp=$(mktemp)
    get "$1" "$tmp" || rc=$?
    rm -f "$tmp"
    return $rc
}
put() {
    wr r2 object put "$BUCKET/$1" --file "$2" >/dev/null 2>&1 || die "writing $BUCKET/$1 failed"
}
# put_checked KEY FILE: write, read back, compare byte for byte.
put_checked() {
    local back
    put "$1" "$2"
    back=$(mktemp)
    get "$1" "$back" || die "$BUCKET/$1 is missing right after it was written"
    cmp -s "$2" "$back" || { rm -f "$back"; die "$BUCKET/$1 differs from what was written"; }
    rm -f "$back"
}

# upload_tree DIR MANIFEST PREFIX: every file MANIFEST lists, under PREFIX.
upload_tree() {
    local dir=$1 manifest=$2 prefix=$3 list
    list=$(mktemp)
    python3 -I - "$dir" "$manifest" "$prefix" >"$list" <<'PY'
import json, os, re, sys
root, manifest, prefix = sys.argv[1:]
entries = []
for line in open(manifest, encoding="ascii"):
    m = re.fullmatch(r"[0-9a-f]{64} [ *](?:\./)?(.+)\n", line)
    entries.append({"key": prefix + m.group(1), "file": os.path.join(root, m.group(1))})
json.dump(entries, sys.stdout)
PY
    wr r2 bulk put "$BUCKET" --filename "$list" --concurrency 8 >/dev/null 2>&1 ||
        { rm -f "$list"; die "uploading $dir to $BUCKET/$prefix failed"; }
    rm -f "$list"
}

# verify_objects PREFIX MANIFEST: every object MANIFEST lists exists under
# PREFIX with the listed SHA-256. Read back through wrangler, JOBS at a time.
verify_objects() {
    local prefix=$1 manifest=$2 work n bad
    work=$(mktemp -d)
    n=$(wc -l <"$manifest")
    sed -E 's/^([0-9a-f]{64}) [ *](\.\/)?(.*)$/\1 \3/' "$manifest" |
        GRM_V_WORK=$work GRM_V_PREFIX=$prefix xargs -P "$JOBS" -L 1 bash -c '
            want=$0 path=$1 f=$(mktemp -p "$GRM_V_WORK")
            if "$GRM_V_SELF" "$GRM_V_PREFIX$path" "$f" && [ "$(sha256sum <"$f" | cut -d" " -f1)" = "$want" ]
            then echo ok >>"$GRM_V_WORK/ok.$$"; else echo "$path" >>"$GRM_V_WORK/bad.$$"; fi
            rm -f "$f"' || true
    bad=$(cat "$work"/bad.* 2>/dev/null | LC_ALL=C sort || true)
    local ok
    ok=$(cat "$work"/ok.* 2>/dev/null | wc -l)
    rm -rf "$work"
    [ -z "$bad" ] && [ "$ok" -eq "$n" ] || die "$BUCKET/$prefix: $ok of $n objects match; differing or missing: $(echo $bad | head -c 400)"
    say "verified $BUCKET/$prefix: $ok/$n objects match"
}
export GRM_V_SELF=$HERE/get-one.sh

# The pair record that state/active and state/previous hold, read into ST_*.
read_pair() {
    python3 -I - "$1" <<'PY'
import json, re, sys
p = json.load(open(sys.argv[1]))
name = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
for k in ("app", "publication"):
    assert name.match(p[k]) and "--" not in p[k], k
for k in ("app_manifest", "report_manifest"):
    assert re.fullmatch(r"[0-9a-f]{64}", p[k]), k
print(p["app"], p["publication"], p["app_manifest"], p["report_manifest"], p["schema"], p["revision"], p["datasets"])
PY
}
