#!/usr/bin/env bash
# Point the site at a sealed app release and a sealed publication.
# usage: activate.sh LABEL NAME [--apply]
# Both must be sealed, their claims must name their seals' digests, and every
# object of both is read back and hashed first. Then state/previous gets the
# pair state/active names (if any), a history record is written, and
# state/active gets the new pair: one PUT, so the Worker serves either the old
# pair or the new one, never a mix. With GRM_SITE_URL set, smoke.py then
# checks the served site, and a failure restores both state objects (exit 1).
# rollback.sh runs this with --rollback, which also requires the pair to be
# the one state/previous names.
ROLLBACK=0
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
parse_args "$@"
set -- "${ARGS[@]}"
[ $# -eq 2 ] || die "usage: activate.sh LABEL NAME [--apply]" 2
label=$1 name=$2
check_name "$label"
check_name "$name"
A=app/$label P=pub/$name
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

get "$A/MANIFEST.sha256" "$work/app.sha256" || die "$BUCKET/$A is not a sealed app release"
get "$A/CLAIM" "$work/app.claim" || die "$BUCKET/$A has no claim"
get "$P/SOURCE.sha256" "$work/report.sha256" || die "$BUCKET/$P is not a sealed publication"
get "$P/CLAIM" "$work/pub.claim" || die "$BUCKET/$P has no claim"
app_digest=$(sha256sum <"$work/app.sha256" | cut -d' ' -f1)
report_digest=$(sha256sum <"$work/report.sha256" | cut -d' ' -f1)
read -r commit schema revision datasets < <(python3 -I - "$work" "$label" "$name" "$app_digest" "$report_digest" <<'PY'
import json, os, sys
work, label, name, app_digest, report_digest = sys.argv[1:]
a = json.load(open(os.path.join(work, "app.claim")))
p = json.load(open(os.path.join(work, "pub.claim")))
assert a["app"] == label and a["app_manifest"] == app_digest, "app claim does not match its seal"
assert p["publication"] == name and p["report_manifest"] == report_digest, "publication claim does not match its seal"
print(a["commit"], p["schema"], p["revision"], p["datasets"])
PY
) || die "claims and seals disagree"
[ -n "${datasets:-}" ] || die "claims and seals disagree"
grep -q '^[0-9a-f]\{64\}  \./index\.html$' "$work/app.sha256" || die "$BUCKET/$A has no index.html"
verify_objects "$A/dist/" "$work/app.sha256"
verify_objects "$P/report/" "$work/report.sha256"

cur=''
if get state/active "$work/active"; then
    cur=$(read_pair "$work/active") || die "state/active is not a valid pair record"
    set -- $cur
    [ "$1 $2" != "$label $name" ] || die "$label with $name is already active" 2
    say "active now: $1 with $2"
else
    say "nothing is active yet"
fi
had_previous=0
if get state/previous "$work/previous"; then had_previous=1; fi
if ((ROLLBACK)); then
    ((had_previous)) || die "state/previous is missing" 2
    set -- $(read_pair "$work/previous")
    [ "$1 $2" = "$label $name" ] || die "state/previous now names $1 with $2, not $label with $name" 2
fi

now=$(date -u +%FT%TZ)
printf '{"app":"%s","publication":"%s","commit":"%s","app_manifest":"%s","report_manifest":"%s","schema":%d,"revision":%d,"datasets":%d,"time":"%s"}\n' \
    "$label" "$name" "$commit" "$app_digest" "$report_digest" "$schema" "$revision" "$datasets" "$now" >"$work/new"
history=state/history/${now//:/}--$label--$name
if ! ((APPLY)); then
    [ -n "$cur" ] && would "write state/previous = the active pair ($(set -- $cur; echo "$1 with $2"))"
    would "write $history and state/active = $label with $name (revision $revision, $datasets datasets)"
    [ -n "${GRM_SITE_URL:-}" ] && would "check $GRM_SITE_URL and restore the old state if it fails"
    exit 0
fi

[ -n "$cur" ] && { put_checked state/previous "$work/active"; say "state/previous = $(set -- $cur; echo "$1 with $2")"; }
put_checked "$history" "$work/new"
put_checked state/active "$work/new"
say "state/active = $label with $name (commit $commit, revision $revision, $datasets datasets)"
[ -n "${GRM_SITE_URL:-}" ] || { say "activated $label with $name; GRM_SITE_URL unset, no served check"; exit 0; }

if python3 -I "$HERE/smoke.py" "$GRM_SITE_URL" "$work/app.sha256" "$work/report.sha256" "$schema" "$revision" "$datasets" >"$work/smoke" 2>&1; then
    say "$(tail -n 1 "$work/smoke")"
    say "activated $label with $name"
    exit 0
fi
while read -r line; do say "$line"; done <"$work/smoke"
if [ -n "$cur" ]; then
    put_checked state/active "$work/active"
    ((had_previous)) && put_checked state/previous "$work/previous"
    say "restored state/active = $(set -- $cur; echo "$1 with $2")$( ((had_previous)) && echo ' and state/previous')"
else
    wr r2 object delete "$BUCKET/state/active" >/dev/null 2>&1 || die "removing state/active failed; the failed pair is still active" 4
    say "removed state/active; nothing was active before"
fi
die "$label with $name failed the served check at $GRM_SITE_URL"
