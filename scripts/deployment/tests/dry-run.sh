#!/usr/bin/env bash
# Exercise the deployment scripts' checks against a throwaway root in a
# temporary directory. Nothing here opens a port, starts a process that
# outlives its command, contacts Cloudflare or touches a real deployment: the
# release's grm is a stub that records any run, the "live" port must have no
# listener, and every script except publish.sh (into the temporary root) runs
# without --apply. tests/apply.sh covers the --apply paths.
# usage: bash scripts/deployment/tests/dry-run.sh
set -euo pipefail
S=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
REPO=$(cd -- "$S/../.." && pwd)
. "$S/tests/tmproot.sh"
T=$(make_test_root) || exit 2
# Every assertion below must run: a command that fails under set -e or
# pipefail (a SIGPIPE in a pipeline, say) would otherwise end the suite early
# without a FAIL line, so the EXIT trap reports it, and the total is checked.
EXPECTED_ASSERTIONS=0  # set below, before the summary
finished=0
on_exit() {
    local rc=$?
    if ((!finished)); then
        echo "FAIL dry-run suite aborted (exit $rc) after $((pass + fail)) assertions; a command outside expect failed"
        rc=1
    fi
    chmod -R u+w "$T"
    rm -rf "$T"
    exit "$rc"
}
trap on_exit EXIT
D=$T/root
LIVE_PORT=1
STAGE_PORT=2
for port in $LIVE_PORT $STAGE_PORT; do
    [ -z "$(fuser -n tcp "$port" 2>/dev/null)" ] || { echo "port $port is in use; cannot run" >&2; exit 2; }
done
export GRM_DEPLOY_ROOT=$D GRM_PORT=$LIVE_PORT GRM_EXPECT_SCHEMA=2 GRM_EXPECT_REVISION=7 GRM_EXPECT_DATASETS=2
unset GRM_PLAYWRIGHT_CORE GRM_SKIP_BROWSER GRM_SOURCE_REPO GRM_TUNNEL_TOKEN_FILE
export GRM_TEST_MARKER=$T/grm-ran
pass=0
fail=0
# Every guard against a hang is "timeout -k 5 N": a script that traps TERM and
# is blocked in open() never runs its trap, so TERM alone could wait forever.

expect() {  # expect CODE PATTERN NAME -- command...: exit code and output match
    local code=$1 pattern=$2 name=$3 got out
    shift 4
    set +e
    out=$("$@" 2>&1)
    got=$?
    set -e
    if [ "$got" = "$code" ] && grep -Eq -- "$pattern" <<<"$out"; then
        pass=$((pass + 1))
        echo "ok   $name"
    else
        fail=$((fail + 1))
        printf 'FAIL %s (exit %s, wanted %s /%s/)\n%s\n' "$name" "$got" "$code" "$pattern" "$out" | sed '2,$s/^/     /'
    fi
}
snapshot() { (cd "$D" && find . -printf '%p %m %s %T@\n' | LC_ALL=C sort | sha256sum); }
# lines FIRST [LAST]: input lines FIRST to LAST (default the end). It reads its
# whole input, so the writer of a pipeline into it never gets SIGPIPE (head
# exits early, and under pipefail a SIGPIPE'd writer fails the pipeline).
lines() { sed -n "$1,${2:-\$}p"; }
digest() { sha256sum <"$1" | cut -d' ' -f1; }
manifest() { (cd "$1" && find . -type f ! -path ./MANIFEST.sha256 -print0 | LC_ALL=C sort -z | xargs -0 sha256sum); }
# The approval stage.sh writes for LABEL with PUB, then extra or replaced fields.
approval() {
    local label=$1 pub=$2
    shift 2
    printf 'stage=0 browser=skipped commit=%s publication=%s release_manifest=%s report_manifest=%s schema=2 revision=7 datasets=2 receipts=stage/x%s\n' \
        "$(cat "$D/release/$label/COMMIT")" "$pub" "$(digest "$D/release/$label/MANIFEST.sha256")" \
        "$(digest "$D/published/$pub/SOURCE.sha256")" "${*:+ $*}"
}
# The state record activate.sh writes for LABEL with PUB.
state() {
    printf 'label=%s\npublication=%s\ncommit=%s\nrelease_manifest=%s\nreport_manifest=%s\nschema=2\nrevision=7\ndatasets=%s\nat=2026-10-27T00:00:00Z\n' \
        "$1" "$2" "$(cat "$D/release/$1/COMMIT")" "$(digest "$D/release/$1/MANIFEST.sha256")" \
        "$(digest "$D/published/$2/SOURCE.sha256")" "${3:-2}"
}

# Fixture: two-dataset report and two releases.
src=$T/src-report
mkdir -p "$src/datasets"
printf '{"schema":2,"evaluation_revision":7,"datasets":[{"id":"sp__a"},{"id":"sp__b"}]}\n' >"$src/index.json"
echo '{"corpus":1}' >"$src/corpus.json"
echo '{"id":"sp__a"}' >"$src/datasets/sp__a.json"
echo '{"id":"sp__b"}' >"$src/datasets/sp__b.json"
(cd "$src" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum) >"$T/SOURCE.sha256"
mkdir -p "$D/release"
make_release() {
    local R=$D/release/$1
    mkdir -p "$R/dist/assets"
    printf '#!/bin/sh\necho ran >>"$GRM_TEST_MARKER"\nexit 99\n' >"$R/grm"
    chmod +x "$R/grm"
    echo '<!doctype html><title>GRM</title>' >"$R/dist/index.html"
    echo "console.log('$1')" >"$R/dist/assets/index-test.js"
    printf '%s\n' "$2" >"$R/COMMIT"
    manifest "$R" >"$R/MANIFEST.sha256"
}
c1=1111111111111111111111111111111111111111
c2=2222222222222222222222222222222222222222
make_release rel1 "$c1"
make_release rel2 "$c2"

# Each file separately: "bash -n a.sh b.sh" checks only a.sh.
for f in "$S"/*.sh "$S"/tests/*.sh; do
    expect 0 '' "bash -n ${f#"$REPO"/}" -- bash -n "$f"
done
printf 'echo ok\n' >"$T/good.sh"
printf 'if then\n' >"$T/bad.sh"
expect 0 '' "a multi-file bash -n misses a broken second file" -- bash -n "$T/good.sh" "$T/bad.sh"
expect 1 'cannot create a temporary directory' "a failed mktemp stops a suite instead of using the current directory" -- \
    env TMPDIR="$T/missing" bash -c '. "$1"; make_test_root' _ "$S/tests/tmproot.sh"
expect 1 'contains the repository' "a temporary root that contains the repository is refused" -- \
    env GRM_FAKE_TMP="$REPO" bash -c 'mktemp() { printf "%s\n" "$GRM_FAKE_TMP"; }; . "$1"; make_test_root' _ "$S/tests/tmproot.sh"
expect 123 'syntax error' "the documented per-file check catches it" -- \
    sh -c 'printf "%s\n" "$1" "$2" | xargs -n1 bash -n' _ "$T/good.sh" "$T/bad.sh"

# Publication: dry run writes nothing; --apply makes a read-only copy.
expect 0 'dry-run: would copy' "publish dry run" -- bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub1
expect 0 '' "publish dry run wrote nothing" -- test ! -e "$D/published"
expect 0 'published pub1 files=4' "publish --apply into temporary root" -- \
    bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub1 --apply
expect 2 'already exists' "publish refuses an existing name" -- \
    bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub1 --apply
expect 1 '' "publication is read-only" -- test -w "$D/published/pub1/report/index.json"
mkdir -p "$D/go"
approval rel1 pub1 >"$D/go/rel1--pub1"
variant() {  # copy the source report to $T/NAME and print its path
    rm -rf "${T:?}/$1"
    cp -a "$src" "$T/$1"
    printf '%s\n' "$T/$1"
}
v=$(variant extra)
echo '{"id":"sp__c"}' >"$v/datasets/sp__c.json"
expect 1 'files its manifest does not list' "publish refuses an unlisted file" -- \
    bash "$S/publish.sh" "$v" "$T/SOURCE.sha256" pub-extra
v=$(variant inner-manifest)
cp "$T/SOURCE.sha256" "$v/MANIFEST.sha256"
expect 1 'files its manifest does not list' "publish refuses an unlisted report/MANIFEST.sha256" -- \
    bash "$S/publish.sh" "$v" "$T/SOURCE.sha256" pub-inner
v=$(variant link-unlisted)
ln -s ../corpus.json "$v/datasets/link.json"
expect 1 'link or special file: datasets/link.json' "publish refuses an unlisted symlink" -- \
    bash "$S/publish.sh" "$v" "$T/SOURCE.sha256" pub-link
v=$(variant link-listed)
cp "$v/datasets/sp__b.json" "$T/sp__b.copy"
ln -sf "$T/sp__b.copy" "$v/datasets/sp__b.json"
expect 1 'link or special file: datasets/sp__b.json' "publish refuses a listed file that is a symlink with matching bytes" -- \
    bash "$S/publish.sh" "$v" "$T/SOURCE.sha256" pub-link2
v=$(variant fifo)
mkfifo "$v/datasets/pipe"
expect 1 'link or special file: datasets/pipe' "publish refuses a FIFO" -- \
    bash "$S/publish.sh" "$v" "$T/SOURCE.sha256" pub-fifo
v=$(variant outside)
{ cat "$T/SOURCE.sha256"; sha256sum "$T/sp__b.copy" | sed "s#  .*#  ../sp__b.copy#"; } >"$T/outside.sha256"
expect 1 'files its manifest does not list' "publish refuses a manifest entry outside the tree" -- \
    bash "$S/publish.sh" "$v" "$T/outside.sha256" pub-outside
expect 1 'datasets' "publish refuses an unexpected dataset count" -- \
    env GRM_EXPECT_DATASETS=3 bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub-count
expect 1 'evaluation_revision 7, expected 6' "publish refuses an unexpected revision" -- \
    env GRM_EXPECT_REVISION=6 bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub-rev
expect 2 'invalid name' "publish refuses a path-like name" -- \
    bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" ../escape

# Damaged and malformed releases and publications.
cp -a "$D/published/pub1" "$D/published/pub-bad"
chmod -R u+w "$D/published/pub-bad"
echo '{"id":"sp__a","changed":1}' >"$D/published/pub-bad/report/datasets/sp__a.json"
approval rel1 pub-bad >"$D/go/rel1--pub-bad"
cp -a "$D/release/rel1" "$D/release/rel-bad"
chmod -R u+w "$D/release/rel-bad"
echo tampered >>"$D/release/rel-bad/dist/index.html"
approval rel-bad pub1 >"$D/go/rel-bad--pub1"
broken_release() {  # copy rel1 to NAME, writable, with its approval
    cp -a "$D/release/rel1" "$D/release/$1"
    chmod -R u+w "$D/release/$1"
}
broken_release rel-link
mv "$D/release/rel-link/grm" "$D/release/rel-link/grm.real"
ln -s grm.real "$D/release/rel-link/grm"
manifest "$D/release/rel-link" | grep -v ' \./grm$' >"$D/release/rel-link/MANIFEST.sha256"
approval rel-link pub1 >"$D/go/rel-link--pub1"
broken_release rel-nogrm
rm "$D/release/rel-nogrm/grm"
manifest "$D/release/rel-nogrm" >"$D/release/rel-nogrm/MANIFEST.sha256"
approval rel-nogrm pub1 >"$D/go/rel-nogrm--pub1"
broken_release rel-dirlink
mv "$D/release/rel-dirlink/dist/assets" "$T/assets"
ln -s "$T/assets" "$D/release/rel-dirlink/dist/assets"
broken_release rel-extra
echo extra >"$D/release/rel-extra/extra"
ln -s rel1 "$D/release/rel-alias"
cp -a "$D/published/pub1" "$D/published/pub-unapproved"
echo "stage=0 schema1=0 browser=1042/1042 commit=$c2" >"$D/GO-rel2"

before=$(snapshot)
expect 0 'dry-run: would start rel1 with pub1' "restore dry run of an approved pair" -- \
    bash "$S/restore.sh" rel1 pub1
expect 1 'GO-rel2 is a legacy approval .*approve-legacy.sh rel2 pub1' "a legacy GO-LABEL alone approves nothing" -- \
    bash "$S/restore.sh" rel2 pub1
expect 1 'no state/active' "restore without a target fails" -- bash "$S/restore.sh"
expect 1 'does not match' "restore refuses a damaged publication" -- bash "$S/restore.sh" rel1 pub-bad
expect 1 'does not match' "restore refuses a damaged release" -- bash "$S/restore.sh" rel-bad pub1
expect 1 'link or special file: grm' "restore refuses an unlisted grm symlink" -- bash "$S/restore.sh" rel-link pub1
expect 1 'incomplete: its manifest does not list grm' "restore refuses a release without grm" -- \
    bash "$S/restore.sh" rel-nogrm pub1
expect 1 'link or special file: dist/assets' "restore refuses a linked directory in a release" -- \
    bash "$S/restore.sh" rel-dirlink pub1
expect 1 'files its manifest does not list' "restore refuses an unlisted release file" -- \
    bash "$S/restore.sh" rel-extra pub1
expect 1 'release rel-alias is missing' "restore refuses a release directory that is a link" -- \
    bash "$S/restore.sh" rel-alias pub1
expect 0 'GRM_EXPECT_\* is not used here' "restore checks the approval's dataset count, not GRM_EXPECT_DATASETS" -- \
    env GRM_EXPECT_DATASETS=237 bash "$S/restore.sh" rel1 pub1
expect 1 'no stage approval' "restore refuses an unapproved pair" -- bash "$S/restore.sh" rel1 pub-unapproved
expect 1 'missing directory' "restore refuses a missing publication" -- bash "$S/restore.sh" rel1 nope
expect 2 'nothing listens' "activate refuses without a live managed service" -- bash "$S/activate.sh" rel1 pub1
expect 1 'no stage approval' "activate refuses an unapproved pair" -- \
    bash "$S/activate.sh" rel1 pub-unapproved
expect 1 'no state/previous' "rollback without a recorded previous release fails" -- bash "$S/rollback.sh"
expect 2 'live port' "stage refuses the live port" -- \
    env GRM_SKIP_BROWSER=1 bash "$S/stage.sh" rel2 pub1 "$LIVE_PORT"
expect 1 'GRM_PLAYWRIGHT_CORE' "stage requires the browser check or an explicit skip" -- \
    bash "$S/stage.sh" rel2 pub1 "$STAGE_PORT"
expect 0 'dry-run: would serve rel2' "stage dry run" -- \
    env GRM_SKIP_BROWSER=1 bash "$S/stage.sh" rel2 pub1 "$STAGE_PORT"
expect 1 'link or special file' "stage refuses a linked grm" -- \
    env GRM_SKIP_BROWSER=1 bash "$S/stage.sh" rel-link pub1 "$STAGE_PORT"
expect 2 'unknown option' "scripts reject unknown options" -- bash "$S/restore.sh" rel1 pub1 --force
expect 2 'absolute' "a relative deployment root is refused" -- \
    env GRM_DEPLOY_ROOT=relative bash "$S/restore.sh" rel1 pub1
expect 0 '' "dry runs left the root unchanged" -- test "$before" = "$(snapshot)"

# Approvals bind the exact commit, build, publication, report copy and index.
cp -a "$D/published/pub1" "$D/published/pub2"
chmod u+w "$D/published/pub2" "$D/published/pub2/SOURCE.sha256"
LC_ALL=C sort -r "$T/SOURCE.sha256" >"$D/published/pub2/SOURCE.sha256"
cp -a "$D/published/pub1" "$D/published/pub.1"
tamper() {  # tamper NAME PATTERN SED-EXPRESSION|FIELDS...: a modified rel1--pub2 approval is refused
    local name=$1 pattern=$2
    shift 2
    approval rel1 pub2 >"$D/go/rel1--pub2"
    if [ "${1:-}" = -e ]; then sed -i -E "$2" "$D/go/rel1--pub2"; else approval rel1 pub2 "$@" >"$D/go/rel1--pub2"; fi
    expect 1 "$pattern" "$name" -- bash "$S/restore.sh" rel1 pub2
}
approval rel1 pub2 >"$D/go/rel1--pub2"
expect 0 'dry-run: would start rel1 with pub2' "an exact approval is accepted" -- bash "$S/restore.sh" rel1 pub2
tamper "approval for another publication is refused" 'another publication' -e 's/publication=pub2/publication=other/'
tamper "approval without report_manifest is refused" 'lacks report_manifest' -e 's/ report_manifest=[0-9a-f]+//'
tamper "approval of another copy of the publication is refused" 'another copy of pub2' \
    -e "s/report_manifest=[0-9a-f]+/report_manifest=$(digest "$D/published/pub1/SOURCE.sha256")/"
tamper "approval for another build is refused" 'another build' -e "s/release_manifest=[0-9a-f]+/release_manifest=$(printf 0%.0s {1..64})/"
tamper "approval of another commit is refused" 'does not approve commit' -e "s/commit=$c1/commit=$c2/"
tamper "approval without a commit is refused" 'lacks commit' -e "s/ commit=$c1//"
tamper "approval without release_manifest is refused" 'lacks release_manifest' -e 's/ release_manifest=[0-9a-f]+//'
tamper "approval without datasets is refused" 'lacks datasets' -e 's/ datasets=2//'
tamper "a failed stage is not an approval" 'does not approve: stage=1' -e 's/stage=0/stage=1/'
tamper "a failed browser check is not an approval" 'browser=fail' -e 's/browser=skipped/browser=fail/'
tamper "a repeated field is refused" 'unique key=value' publication=pub2
tamper "an unknown field is refused" 'unknown field note' note=x
tamper "a two-line approval is refused" 'one-line approval' -e 's/ receipts=/\nreceipts=/'
tamper "a malformed field is refused" 'unique key=value' 'publication pub2'
tamper "browser=legacy needs legacy_sha256" 'legacy migration without legacy_sha256' -e 's/browser=skipped/browser=legacy/'
tamper "an approval's dataset count is checked against the index" 'datasets' -e 's/datasets=2/datasets=3/'
approval rel1 pub2 | sed -E "s/browser=skipped/browser=legacy/; s# receipts=stage/x# legacy_sha256=$(printf a%.0s {1..64})#" >"$D/go/rel1--pub2"
expect 0 'dry-run: would start rel1 with pub2' "a migrated legacy approval is accepted" -- bash "$S/restore.sh" rel1 pub2
approval rel1 pub2 | sed 's/publication=pub2/publication=pubx1/' >"$D/go/rel1--pub.1"
expect 1 'another publication' "publication names are compared exactly, not as a regex" -- \
    bash "$S/restore.sh" rel1 pub.1
ln -s rel1--pub1 "$D/go/rel2--pub1"
expect 1 'one-line approval' "an approval that is a link is refused" -- bash "$S/restore.sh" rel2 pub1
rm "$D/go/rel2--pub1"

# Every byte of an approval is checked before any field is parsed, so a NUL
# cannot hide conflicting fields after it, and the grammar is one line of
# fields separated by single spaces.
nul_then() {  # nul_then NAME FIELDS: rel1--pub2's exact approval, then a NUL and FIELDS on its line
    { approval rel1 pub2 | tr -d '\n'; printf '\0%s\n' "$2"; } >"$D/go/rel1--pub2"
    expect 1 'byte 0x00' "an approval with a NUL then $1 is refused" -- bash "$S/restore.sh" rel1 pub2
}
nul_then "stage=1" ' stage=1'
nul_then "another commit" " commit=$c2"
nul_then "another report_manifest" " report_manifest=$(digest "$D/published/pub1/SOURCE.sha256")"
nul_then "another release_manifest" " release_manifest=$(printf 0%.0s {1..64})"
nul_then "another schema" ' schema=9'
nul_then "another revision" ' revision=6'
nul_then "another dataset count" ' datasets=237'
nul_then "another publication" ' publication=pub1'
nul_then "nothing" ''
raw_approval() {  # raw_approval NAME PATTERN FORMAT: printf FORMAT with the exact approval line as %s
    # shellcheck disable=SC2059
    printf "$3" "$(approval rel1 pub2)" >"$D/go/rel1--pub2"
    expect 1 "$2" "$1" -- bash "$S/restore.sh" rel1 pub2
}
raw_approval "an approval with a CRLF ending is refused" 'byte 0x0d' '%s\r\n'
raw_approval "an approval with a tab is refused" 'byte 0x09' '\t%s\n'
raw_approval "an approval with a DEL byte is refused" 'byte 0x7f' '%s\177\n'
raw_approval "an approval with a non-ASCII byte is refused" 'byte 0xff' '%s\377\n'
raw_approval "an approval without a final newline is refused" 'does not end with a newline' '%s'
raw_approval "an approval with a trailing blank line is refused" '2 lines, expected one' '%s\n\n'
raw_approval "an approval with a leading space is refused" 'empty field' ' %s\n'
raw_approval "an approval with a trailing space is refused" 'empty field' '%s \n'
printf '' >"$D/go/rel1--pub2"
expect 1 'is empty' "an empty approval is refused" -- bash "$S/restore.sh" rel1 pub2
approval rel1 pub2 | sed 's/ commit=/  commit=/' >"$D/go/rel1--pub2"
expect 1 'empty field' "an approval with a doubled space is refused" -- bash "$S/restore.sh" rel1 pub2
rm "$D/go/rel1--pub2"
mkfifo "$D/go/rel1--pub2"
expect 1 'is not a regular file' "an approval that is a FIFO is refused without blocking" -- \
    timeout -k 5 20 bash "$S/restore.sh" rel1 pub2
rm "$D/go/rel1--pub2"
mkdir "$D/go/rel1--pub2"
expect 1 'is not a regular file' "an approval that is a directory is refused" -- bash "$S/restore.sh" rel1 pub2
rmdir "$D/go/rel1--pub2"
approval rel1 pub2 >"$D/go/rel1--pub2"
expect 0 'dry-run: would start rel1 with pub2' "the exact approval is accepted again" -- bash "$S/restore.sh" rel1 pub2

# Legacy GO-LABEL files get the same whole-file check before any field is read.
cp "$D/GO-rel2" "$T/GO-rel2"
legacy() {  # legacy NAME CODE PATTERN FORMAT ARGS...: approve-legacy.sh rel2 pub1 with GO-rel2 = printf FORMAT ARGS
    local name=$1 code=$2 pattern=$3
    shift 3
    # shellcheck disable=SC2059
    printf "$@" >"$D/GO-rel2"
    expect "$code" "$pattern" "$name" -- bash "$S/approve-legacy.sh" rel2 pub1
}
legacy "a clean legacy file reaches the live-pair check" 2 'only the live pair' 'stage=0 commit=%s\n' "$c2"
legacy "legacy: NUL then stage=1 is refused" 1 'byte 0x00' 'stage=0 commit=%s\0 stage=1\n' "$c2"
legacy "legacy: NUL then another commit is refused" 1 'byte 0x00' 'stage=0 commit=%s\0 commit=%s\n' "$c2" "$c1"
legacy "legacy: NUL then another publication is refused" 1 'byte 0x00' 'stage=0 commit=%s\0 publication=other\n' "$c2"
legacy "legacy: NUL then another report, schema, revision and count is refused" 1 'byte 0x00' \
    'stage=0 commit=%s\0 report_manifest=%s schema=9 revision=6 datasets=237\n' "$c2" "$(printf 0%.0s {1..64})"
legacy "legacy: a NUL before every field is refused" 1 'byte 0x00' '\0stage=0 commit=%s\n' "$c2"
legacy "legacy: a repeated commit is refused" 1 'repeats field commit' 'stage=0 commit=%s commit=%s\n' "$c2" "$c2"
legacy "legacy: a repeated stage is refused" 1 'repeats field stage' 'stage=0 stage=1 commit=%s\n' "$c2"
legacy "legacy: two lines are refused" 1 '2 lines, expected one' 'stage=0\ncommit=%s\n' "$c2"
legacy "legacy: a CRLF ending is refused" 1 'byte 0x0d' 'stage=0 commit=%s\r\n' "$c2"
legacy "legacy: a missing stage is refused" 1 'says stage=<none>' 'commit=%s\n' "$c2"
cp "$T/GO-rel2" "$D/GO-rel2"

# Recorded state: rollback resolves state/previous and still needs a live service.
approval rel2 pub1 >"$D/go/rel2--pub1"
mkdir -p "$D/state"
state rel1 pub1 >"$D/state/previous"
state rel2 pub1 >"$D/state/active"
before=$(snapshot)
expect 2 'nothing listens' "rollback dry run resolves state/previous" -- bash "$S/rollback.sh"
expect 0 'dry-run: would start rel2 with pub1' "restore dry run uses state/active" -- bash "$S/restore.sh"
expect 0 '' "dry runs with recorded state left the root unchanged" -- test "$before" = "$(snapshot)"
state rel2 pub1 3 >"$D/state/active"
expect 1 'state/active does not match' "restore refuses state whose dataset count disagrees with the approval" -- \
    bash "$S/restore.sh"
printf 'label=rel2\npublication=pub1\n' >"$D/state/active"
expect 1 'state/active lacks commit' "restore refuses state without the recorded expectations" -- bash "$S/restore.sh"
state rel1 pub1 | sed 's/^label=rel1/label=rel1 x/' >"$D/state/previous"
expect 1 'malformed' "rollback refuses a malformed state file" -- bash "$S/rollback.sh"
state rel1 pub1 | sed 's/^label=rel1/label=..\/x/' >"$D/state/previous"
expect 2 'invalid name' "rollback refuses a path-like name in state" -- bash "$S/rollback.sh"

# State files get the same byte and grammar checks, and a link or special
# file is never opened.
bad_state() {  # bad_state NAME PATTERN: state/active as prepared is refused
    expect 1 "$2" "restore refuses state/active $1" -- timeout -k 5 30 bash "$S/restore.sh"
}
{ state rel2 pub1; printf '\0label=rel1\n'; } >"$D/state/active"
bad_state "with a NUL then another label" 'malformed: it holds byte 0x00'
{ state rel2 pub1 | lines 1 1 | tr -d '\n'; printf '\0x\n'; state rel2 pub1 | lines 2; } >"$D/state/active"
bad_state "with a NUL inside its first line" 'malformed: it holds byte 0x00'
{ state rel2 pub1 | lines 1 1; echo; state rel2 pub1 | lines 2; } >"$D/state/active"
bad_state "with a blank line" 'malformed: it has an empty field'
state rel2 pub1 | head -c -1 >"$D/state/active"
bad_state "without a final newline" 'does not end with a newline'
state rel2 pub1 | sed 's/$/\r/' >"$D/state/active"
bad_state "with CRLF endings" 'byte 0x0d'
{ state rel2 pub1; echo label=rel1; } >"$D/state/active"
bad_state "with a repeated label" 'repeats field label'
state rel2 pub1 | lines 1 4 >"$D/state/active"
bad_state "cut short by a partial write" 'lacks report_manifest'
# The same cut from a writer far larger than a pipe buffer: lines drains it, so
# under pipefail the pipeline succeeds and every later assertion still runs;
# head on the same input is killed by SIGPIPE (exit 141).
big_state() { state rel2 pub1; seq 1 400000 | sed 's/^/pad=/'; }
big=$(big_state | lines 1 9)
expect 0 '' "a drain-safe cut of a 2.7 MB pipeline succeeds under pipefail" -- test "$big" = "$(state rel2 pub1)"
expect 141 '' "head on the same pipeline fails under pipefail with SIGPIPE (why lines is used)" -- \
    bash -c 'set -o pipefail; seq 1 400000 | sed "s/^/pad=/" | head -n 9 >/dev/null'
big_state | lines 1 4 >"$D/state/active"
bad_state "cut short from a large pipeline" 'lacks report_manifest'
rm "$D/state/active"
mkfifo "$D/state/active"
bad_state "that is a FIFO, without blocking" 'is not a regular file'
rm "$D/state/active"
mkdir "$D/state/active"
bad_state "that is a directory" 'is not a regular file'
rmdir "$D/state/active"
state rel2 pub1 >"$T/active-elsewhere"
ln -s "$T/active-elsewhere" "$D/state/active"
bad_state "that is a link to a valid record" 'is a link'
rm "$D/state/active"
state rel2 pub1 >"$D/state/active"
state rel1 pub1 >"$D/state/previous"

# A journal of an interrupted start blocks activation and is reported; only
# restore.sh --apply resolves it.
gone=$(sh -c 'echo $$')
tok=00000000-0000-4000-8000-000000000000
mark="mark=$tok.abcdef mark_id=1:1"
printf 'token=%s owner=%s owner_start=1 label=rel2 publication=pub1 port=%s %s\n' \
    "$tok" "$gone" "$LIVE_PORT" "$mark" >"$D/serve-$LIVE_PORT.starting"
before=$(snapshot)
expect 0 'dry-run: would stop pids none unless serve-1.pid shows' "restore dry run reports an interrupted start" -- bash "$S/restore.sh"
expect 2 'interrupted start of rel2 with pub1' "activate refuses while an interrupted start is unresolved" -- \
    bash "$S/activate.sh" rel1 pub1
expect 2 'interrupted start of rel2 with pub1' "rollback refuses while an interrupted start is unresolved" -- bash "$S/rollback.sh"
expect 1 'interrupted start of rel2 with pub1 left pids: none' "status reports the interrupted start" -- bash "$S/status.sh"
expect 0 '' "dry runs with a journal left the root unchanged" -- test "$before" = "$(snapshot)"
read -r -a self_stat <<<"$(sed 's/.*) //' "/proc/$$/stat")"
printf 'token=%s owner=%s owner_start=%s label=rel2 publication=pub1 port=%s %s\n' \
    "$tok" "$$" "${self_stat[19]}" "$LIVE_PORT" "$mark" >"$D/serve-$LIVE_PORT.starting"
expect 2 "pid $$ is still starting rel2" "restore refuses while the journal's owner is alive" -- bash "$S/restore.sh"
bad_journal() {  # bad_journal NAME PATTERN FIELDS: a journal of a dead owner with FIELDS is refused
    printf '%s\n' "$3" >"$D/serve-$LIVE_PORT.starting"
    expect 2 "$2" "restore refuses a journal $1" -- bash "$S/restore.sh"
}
bad_journal "without a spawn marker (an older format)" 'lacks mark' \
    "token=$tok owner=$gone owner_start=1 label=rel2 publication=pub1 port=$LIVE_PORT"
bad_journal "whose marker names another token" 'malformed spawn marker' \
    "token=$tok owner=$gone owner_start=1 label=rel2 publication=pub1 port=$LIVE_PORT mark=11111111-0000-4000-8000-000000000000.abcdef mark_id=1:1"
bad_journal "with a child but no child start" 'malformed child registration' \
    "token=$tok owner=$gone owner_start=1 label=rel2 publication=pub1 port=$LIVE_PORT $mark child=5"
bad_journal "with an unknown field" 'unknown fields' \
    "token=$tok owner=$gone owner_start=1 label=rel2 publication=pub1 port=$LIVE_PORT $mark note=x"
printf 'token=x owner=1\0 port=1\n' >"$D/serve-$LIVE_PORT.starting"
expect 2 'serve-1.starting holds byte 0x00' "restore refuses a malformed journal" -- bash "$S/restore.sh"
rm "$D/serve-$LIVE_PORT.starting"
mkdir "$D/serve-$LIVE_PORT.starting"
expect 2 'serve-1.starting is not a regular file' "restore refuses a journal path that is a directory" -- bash "$S/restore.sh"
rmdir "$D/serve-$LIVE_PORT.starting"
for path in "serve-$LIVE_PORT.pid" "logs/serve-$LIVE_PORT.log"; do
    mkdir -p "$D/logs"
    mkfifo "$D/$path"
    expect 2 'link, directory or special file; nothing started' "restore dry run refuses $path as a FIFO, without blocking" -- \
        timeout -k 5 30 bash "$S/restore.sh"
    rm "$D/$path"
done
rmdir "$D/logs"

# Build: dry run checks the commit and that every output path is new.
mkdir -p "$T/bin"
printf '#!/bin/sh\necho stub >>"$GRM_TEST_MARKER"\nexit 99\n' >"$T/bin/cargo"
cp "$T/bin/cargo" "$T/bin/npm"
chmod +x "$T/bin/cargo" "$T/bin/npm"
head=$(git -C "$REPO" rev-parse HEAD)
expect 0 "would extract $head" "build dry run" -- \
    env GRM_EXTRA_PATH="$T/bin" GRM_MIN_FREE_GIB=0 bash "$S/build-release.sh" rel3 HEAD
expect 2 'already exists' "build refuses an existing release label" -- \
    env GRM_EXTRA_PATH="$T/bin" GRM_MIN_FREE_GIB=0 bash "$S/build-release.sh" rel1 HEAD
expect 2 'not a commit' "build refuses an unknown commit" -- \
    env GRM_EXTRA_PATH="$T/bin" GRM_MIN_FREE_GIB=0 bash "$S/build-release.sh" rel3 no-such-ref
expect 2 'GiB free' "build refuses when space is short" -- \
    env GRM_EXTRA_PATH="$T/bin" GRM_MIN_FREE_GIB=999999 bash "$S/build-release.sh" rel3 HEAD

# Connector: check validates the exact binary hash and token permissions only.
C=$T/connector
mkdir -m 700 "$C"
printf '#!/bin/sh\necho "cloudflared version test"\n' >"$C/cloudflared"
chmod 700 "$C/cloudflared"
(cd "$C" && sha256sum cloudflared >cloudflared.sha256)
secret=placeholder-not-a-credential-$RANDOM
printf '%s\n' "$secret" >"$C/tunnel.token"
chmod 600 "$C/tunnel.token"
connector() { env GRM_CONNECTOR_DIR="$C" bash "$S/connector.sh" "$@"; }
expect 0 'connector check pass' "connector check" -- connector check
set +e
out=$(connector check 2>&1)
set -e
expect 1 '' "connector check does not print the token" -- grep -qF "$secret" <<<"$out"
chmod 644 "$C/tunnel.token"
expect 2 'must not be group- or world-accessible' "connector refuses a readable token" -- connector check
chmod 600 "$C/tunnel.token"
expect 2 'missing, empty or a link' "connector refuses a missing token file" -- \
    env GRM_CONNECTOR_DIR="$C" GRM_TUNNEL_TOKEN_FILE="$C/none" bash "$S/connector.sh" check
# run changes into the connector directory before cloudflared opens the token,
# so relative paths, checked from $T here, would name other files there.
expect 2 'GRM_TUNNEL_TOKEN_FILE must be an absolute path' "connector refuses a relative token path" -- \
    env -C "$T" GRM_CONNECTOR_DIR="$C" GRM_TUNNEL_TOKEN_FILE=connector/tunnel.token bash "$S/connector.sh" check
expect 2 'GRM_CONNECTOR_DIR must be an absolute path' "connector refuses a relative connector directory" -- \
    env -C "$T" GRM_CONNECTOR_DIR=connector bash "$S/connector.sh" check
good=$(cat "$C/cloudflared.sha256")
sum=${good%% *}
cp "$C/cloudflared" "$C/other"
# cloudflared.sha256 must be exactly one complete LF-terminated entry and no other byte.
one_entry() {  # one_entry NAME FORMAT ARGS...: the manifest printf FORMAT ARGS... is refused
    local name=$1
    shift
    # shellcheck disable=SC2059
    printf "$@" >"$C/cloudflared.sha256"
    expect 2 'exactly one newline-terminated entry, for cloudflared' "connector refuses a manifest with $name" -- connector check
}
one_entry "another file" '%s  other\n' "$sum"
one_entry "an extra complete line naming another path" '%s\n%s  other\n' "$good" "$sum"
one_entry "the valid entry twice" '%s\n%s\n' "$good" "$good"
one_entry "a path-like entry" '%s  ../connector/cloudflared\n' "$sum"
one_entry "a ./.. path" '%s  ./../connector/cloudflared\n' "$sum"
one_entry "an absolute path" '%s  %s\n' "$sum" "$C/cloudflared"
one_entry "unterminated junk after the entry" '%s\njunk' "$good"
one_entry "an unterminated entry" '%s' "$good"
one_entry "a NUL after the entry" '%s\n\0' "$good"
one_entry "a NUL inside the entry line" '%s\0junk\n' "$good"
one_entry "a CRLF line ending (refused: LF only)" '%s\r\n' "$good"
one_entry "a leading blank line" '\n%s\n' "$good"
one_entry "a trailing blank line" '%s\n\n' "$good"
one_entry "a tab separator" '%s\tcloudflared\n' "$sum"
one_entry "an uppercase hash" '%s  cloudflared\n' "${sum^^}"
one_entry "nothing in it" ''
printf '%s *./cloudflared\n' "$sum" >"$C/cloudflared.sha256"
expect 0 'connector check pass' "connector accepts a binary-mode ./cloudflared entry" -- connector check
printf '%s\n' "$good" >"$C/cloudflared.sha256"
expect 0 'connector check pass' "connector accepts the valid one-entry manifest" -- connector check
mv "$C/cloudflared" "$C/cloudflared.real"
ln -s cloudflared.real "$C/cloudflared"
expect 2 'not executable or a link' "connector refuses a cloudflared symlink" -- connector check
rm "$C/cloudflared"
mv "$C/cloudflared.real" "$C/cloudflared"
echo '#' >>"$C/cloudflared"
expect 2 'does not match its SHA-256' "connector refuses a changed binary" -- connector check
expect 2 'usage' "connector needs a verb" -- connector
expect 0 'not running' "connector stop without a supervisor changes nothing" -- connector stop
echo junk >"$C/supervisor.pid"
expect 2 'malformed; no signal sent' "connector stop refuses a malformed pid file" -- connector stop
printf '%s 1\n%s 1\n' "$$" "$$" >"$C/supervisor.pid"
expect 2 'malformed; no signal sent' "connector stop refuses a pid file with a second line" -- connector stop
printf '%s 1\0 2\n' "$$" >"$C/supervisor.pid"
expect 2 'malformed; no signal sent' "connector stop refuses a pid file with a NUL" -- connector stop
rm "$C/supervisor.pid"
mkfifo "$C/supervisor.pid"
expect 2 'malformed; no signal sent' "connector stop refuses a FIFO pid file without blocking" -- timeout -k 5 20 env GRM_CONNECTOR_DIR="$C" bash "$S/connector.sh" stop
rm "$C/supervisor.pid"
gone=$(sh -c 'echo $$')
echo "$gone 1" >"$C/supervisor.pid"
expect 0 'stale.*no signal sent' "connector stop treats an exited pid as stale" -- connector stop
echo "$$ 1" >"$C/supervisor.pid"
expect 0 'stale.*no signal sent' "connector stop does not signal a reused pid" -- connector stop
read -r -a self_stat <<<"$(sed 's/.*) //' "/proc/$$/stat")"
echo "$$ ${self_stat[19]}" >"$C/supervisor.pid"
expect 2 'not the supervisor of this connector directory: .*connector.lock is missing or not a regular file; no signal sent' \
    "connector stop refuses a live process that is not its supervisor" -- connector stop
: >"$C/connector.lock"
expect 2 'no process holds the lock on .*connector.lock; no signal sent' \
    "connector stop refuses a live process while connector.lock is unlocked" -- connector stop
rm "$C/connector.lock"
mkfifo "$C/connector.lock"
expect 2 'not the supervisor of this connector directory.*no signal sent' \
    "connector stop refuses a FIFO connector.lock without blocking" -- timeout -k 5 20 env GRM_CONNECTOR_DIR="$C" bash "$S/connector.sh" stop
rm "$C/connector.lock" "$C/supervisor.pid"
mkfifo "$C/connector.lock"
expect 2 'connector.lock is a link, directory, FIFO' \
    "connector stop without a supervisor refuses a FIFO connector.lock without blocking" -- \
    timeout -k 5 20 env GRM_CONNECTOR_DIR="$C" bash "$S/connector.sh" stop
rm "$C/connector.lock"
mkfifo "$C/cloudflared.sha256.fifo"
mv "$C/cloudflared.sha256" "$T/cloudflared.sha256"
mv "$C/cloudflared.sha256.fifo" "$C/cloudflared.sha256"
expect 2 'cloudflared.sha256 is missing' "connector check refuses a FIFO manifest without blocking" -- \
    timeout -k 5 20 env GRM_CONNECTOR_DIR="$C" bash "$S/connector.sh" check
rm "$C/cloudflared.sha256"
mv "$T/cloudflared.sha256" "$C/cloudflared.sha256"

expect 1 '' "no stub grm, cargo or npm ever ran" -- test -e "$GRM_TEST_MARKER"
EXPECTED_ASSERTIONS=187
if [ "$((pass + fail))" != "$EXPECTED_ASSERTIONS" ]; then
    fail=$((fail + 1))
    echo "FAIL dry-run suite ran $((pass + fail - 1)) assertions, expected $EXPECTED_ASSERTIONS"
fi
finished=1
echo "dry-run tests: $pass passed, $fail failed"
if ((fail)); then exit 1; fi
