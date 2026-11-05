#!/usr/bin/env bash
# Exercise the deployment scripts' checks against a throwaway root in a
# temporary directory. Nothing here opens a port, starts a process that
# outlives its command, contacts Cloudflare or touches a real deployment: the
# release's grm is a stub that records any run, the "live" port must have no
# listener, and every script except publish.sh (into the temporary root) runs
# without --apply, but for stage.sh --apply runs that are refused before they
# start anything. tests/apply.sh covers the --apply paths. It needs a git
# checkout with a commit: the build cases archive HEAD (build-release.sh).
# usage: bash scripts/deployment/tests/dry-run.sh
set -euo pipefail
S=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
REPO=$(cd -- "$S/../.." && pwd)
git -C "$REPO" rev-parse -q --verify HEAD >/dev/null 2>&1 ||
    { echo "FAIL tests/dry-run.sh needs a git checkout with a commit (the build cases archive HEAD); $REPO is not one"; exit 2; }
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
skipped=0
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
mkdir "$T/not-empty"
: >"$T/not-empty/file"
expect 1 'not empty' "a temporary root that is not empty is refused" -- \
    env GRM_FAKE_TMP="$T/not-empty" bash -c 'mktemp() { printf "%s\n" "$GRM_FAKE_TMP"; }; . "$1"; make_test_root' _ "$S/tests/tmproot.sh"
# / is not empty either, so ls is stubbed: only the depth check may refuse it.
expect 1 'too close to /' "a temporary root at / is refused even when it looks empty" -- \
    bash -c 'mktemp() { printf "/\n"; }; ls() { :; }; . "$1"; make_test_root' _ "$S/tests/tmproot.sh"
expect 123 'syntax error' "the documented per-file check catches it" -- \
    sh -c 'printf "%s\n" "$1" "$2" | xargs -n1 bash -n' _ "$T/good.sh" "$T/bad.sh"

# Publication: dry run writes nothing; --apply makes a read-only copy.
expect 0 'dry-run: would copy' "publish dry run" -- bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub1
expect 0 '' "publish dry run wrote nothing" -- test ! -e "$D/published"
expect 0 'published pub1 files=4' "publish --apply into temporary root" -- \
    bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub1 --apply
expect 2 'already exists' "publish refuses an existing name" -- \
    bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub1 --apply
if [ "$(id -u)" != 0 ]; then
    expect 1 '' "publication is read-only" -- test -w "$D/published/pub1/report/index.json"
else
    echo "skip publication is read-only (root ignores the mode)"
    skipped=$((skipped + 1))
fi
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

# From revision 8, every dataset's own evaluation revision, in its index row
# and in its report, must be the expected one too: grm report stamps the index
# with its own revision, so a root left with reports of an earlier revision
# passes the index alone. Revision 7 keeps the index checks only (the report
# above records no per-dataset revision).
audit_json() { [ "$1" = - ] || printf ',"evaluation":{"revision":%s}' "$1"; }
rehash() { (cd "$T/$1" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum) >"$T/$1.sha256"; }
rev_report() {  # rev_report NAME REVISION A B [A-REPORT B-REPORT]: $T/NAME and its manifest, an index at REVISION whose rows give sp__a and sp__b revisions A and B (JSON, or - for none), as do their reports unless given apart
    local r=$T/$1
    rm -rf "$r"
    mkdir -p "$r/datasets"
    printf '{"schema":2,"evaluation_revision":%s,"datasets":[{"id":"sp__a"%s},{"id":"sp__b"%s}]}\n' \
        "$2" "$(audit_json "$3")" "$(audit_json "$4")" >"$r/index.json"
    echo '{"corpus":1}' >"$r/corpus.json"
    printf '{"id":"sp__a"%s}\n' "$(audit_json "${5:-$3}")" >"$r/datasets/sp__a.json"
    printf '{"id":"sp__b"%s}\n' "$(audit_json "${6:-$4}")" >"$r/datasets/sp__b.json"
    rehash "$1"
}
rev_publish() {  # rev_publish REVISION NAME PUB [--apply]: publish $T/NAME as PUB, expecting REVISION
    env GRM_EXPECT_REVISION="$1" bash "$S/publish.sh" "$T/$2" "$T/$2.sha256" "$3" "${@:4}"
}
rev_report r11 11 11 11
rev_report r10 10 10 10
rev_report mix-one 11 11 10          # one dataset left at 10, indexed as grm report does
rev_report mix-all 11 10 10          # every dataset at 10, indexed by a revision-11 grm
rev_report mix-report 11 11 11 11 10 # the index says 11 throughout; one report says 10
rev_report mix-row 11 11 10 11 11
rev_report no-audit 11 - -
rev_report bad-audit 11 '"11"' 11.0
rev_report not-json 11 11 11
echo '{"id":"sp__b",' >"$T/not-json/datasets/sp__b.json"
rehash not-json
rev_report r8-no-audit 8 - -
rev_report r7-mixed 7 7 6
expect 0 'dry-run: would copy' "publish dry run accepts revision 11 with every dataset at 11" -- rev_publish 11 r11 pub-r11
expect 0 'dry-run: would copy' "publish dry run accepts revision 10 with every dataset at 10" -- rev_publish 10 r10 pub-r10
expect 1 'evaluation_revision 10, expected 11' "publish refuses revision 10 when 11 is expected" -- rev_publish 11 r10 pub-r10
expect 1 'evaluation_revision 11, expected 10' "publish refuses revision 11 when 10 is expected" -- rev_publish 10 r11 pub-r11
expect 1 'datasets with index evaluation.revision 10, expected 11: sp__b$' \
    "publish dry run refuses a revision-11 index with one dataset at 10, naming it" -- rev_publish 11 mix-one pub-mix
expect 1 'datasets whose report has evaluation.revision 10, expected 11: sp__a, sp__b$' \
    "publish dry run refuses a revision-11 index over revision-10 datasets, naming them" -- rev_publish 11 mix-all pub-mix
expect 1 'datasets whose report has evaluation.revision 10, expected 11: sp__b$' \
    "publish --apply refuses a revision-11 index with one dataset at 10, naming it" -- rev_publish 11 mix-one pub-mix --apply
expect 1 'datasets with index evaluation.revision 10, expected 11: sp__a, sp__b$' \
    "publish --apply refuses a revision-11 index over revision-10 datasets, naming them" -- rev_publish 11 mix-all pub-mix --apply
expect 0 '' "the refused publish --apply runs published nothing" -- \
    test ! -e "$D/published/pub-mix" -a ! -e "$D/published/pub-mix.partial"
expect 1 'datasets whose report has evaluation.revision 10, expected 11: sp__b$' \
    "publish refuses a report at 10 under an index that says 11 throughout" -- rev_publish 11 mix-report pub-mix
expect 1 'datasets with index evaluation.revision 10, expected 11: sp__b$' \
    "publish refuses an index row at 10 over a report at 11" -- rev_publish 11 mix-row pub-mix
expect 1 'datasets whose report has evaluation.revision missing, expected 11: sp__a, sp__b$' \
    "publish refuses revision-11 datasets that record no revision" -- rev_publish 11 no-audit pub-mix
expect 1 "datasets with index evaluation.revision '11', expected 11: sp__a\$" \
    "publish refuses a dataset revision that is a string" -- rev_publish 11 bad-audit pub-mix
expect 1 'datasets whose report has evaluation.revision 11.0, expected 11: sp__b$' \
    "publish refuses a dataset revision that is not an integer" -- rev_publish 11 bad-audit pub-mix
expect 1 'datasets whose report is not readable JSON: sp__b$' \
    "publish refuses a dataset report that is not JSON" -- rev_publish 11 not-json pub-mix
expect 1 'datasets with index evaluation.revision missing, expected 8: sp__a, sp__b$' \
    "from revision 8, every dataset must record its revision" -- rev_publish 8 r8-no-audit pub-mix
expect 0 'dry-run: would copy' "revision 7 keeps the index checks only" -- rev_publish 7 r7-mixed pub-mix
expect 0 'published pub-r11 files=4' "publish --apply accepts revision 11 with every dataset at 11" -- \
    rev_publish 11 r11 pub-r11 --apply
expect 0 'published pub-r10 files=4' "publish --apply accepts revision 10 with every dataset at 10" -- \
    rev_publish 10 r10 pub-r10 --apply
# stage.sh checks the publication the same way; these are copied in by hand,
# past publish.sh. Its --apply runs stop at the check, before anything starts.
hand_publish() {  # hand_publish NAME PUB: $T/NAME as published/PUB
    mkdir -p "$D/published/$2"
    cp -a "$T/$1" "$D/published/$2/report"
    cp "$T/$1.sha256" "$D/published/$2/SOURCE.sha256"
}
hand_publish mix-one pub-mix-one
hand_publish mix-all pub-mix-all
hand_publish mix-report pub-mix-report
rev_stage() {  # rev_stage REVISION PUB [--apply]: stage rel2 with PUB on the spare port, expecting REVISION
    env GRM_SKIP_BROWSER=1 GRM_EXPECT_REVISION="$1" bash "$S/stage.sh" rel2 "$2" "$STAGE_PORT" "${@:3}"
}
expect 1 'evaluation.revision 10, expected 11: sp__b$' \
    "stage dry run refuses a revision-11 index with one dataset at 10, naming it" -- rev_stage 11 pub-mix-one
expect 1 'evaluation.revision 10, expected 11: sp__b$' \
    "stage --apply refuses a revision-11 index with one dataset at 10, naming it" -- rev_stage 11 pub-mix-one --apply
expect 1 'evaluation.revision 10, expected 11: sp__a, sp__b$' \
    "stage dry run refuses a revision-11 index over revision-10 datasets, naming them" -- rev_stage 11 pub-mix-all
expect 1 'evaluation.revision 10, expected 11: sp__a, sp__b$' \
    "stage --apply refuses a revision-11 index over revision-10 datasets, naming them" -- rev_stage 11 pub-mix-all --apply
expect 1 'datasets whose report has evaluation.revision 10, expected 11: sp__b$' \
    "stage --apply refuses a report at 10 under an index that says 11 throughout" -- rev_stage 11 pub-mix-report --apply
expect 0 '' "the refused stage --apply runs started and approved nothing" -- \
    test ! -e "$D/stage" -a ! -e "$D/spawn" -a ! -e "$D/serve-$STAGE_PORT.starting" -a ! -e "$D/go/rel2--pub-mix-one" -a ! -e "$D/go/rel2--pub-mix-all"
expect 0 'dry-run: would serve rel2 with pub-r11' "stage dry run accepts revision 11 with every dataset at 11" -- rev_stage 11 pub-r11
expect 0 'dry-run: would serve rel2 with pub-r10' "stage dry run accepts revision 10 with every dataset at 10" -- rev_stage 10 pub-r10
expect 1 'evaluation_revision 10, expected 11' "stage refuses revision 10 when 11 is expected" -- rev_stage 11 pub-r10
expect 1 'evaluation_revision 11, expected 10' "stage refuses revision 11 when 10 is expected" -- rev_stage 10 pub-r11
# With --apply, a publication whose datasets all match passes the check and is
# stopped at the next one, before anything starts: its approval path is a
# directory here.
mkdir "$D/go/rel2--pub-r11" "$D/go/rel2--pub-r10"
expect 2 'rel2--pub-r11, serve-2.starting.* is a link, directory or special file' \
    "stage --apply passes revision 11 with every dataset at 11" -- rev_stage 11 pub-r11 --apply
expect 2 'rel2--pub-r10, serve-2.starting.* is a link, directory or special file' \
    "stage --apply passes revision 10 with every dataset at 10" -- rev_stage 10 pub-r10 --apply
rmdir "$D/go/rel2--pub-r11" "$D/go/rel2--pub-r10"
# activate.sh, rollback.sh and restore.sh check it against the approval.
approval rel1 pub-mix-one | sed 's/ revision=7 / revision=11 /' >"$D/go/rel1--pub-mix-one"
approval rel1 pub-r11 | sed 's/ revision=7 / revision=11 /' >"$D/go/rel1--pub-r11"
expect 1 'evaluation.revision 10, expected 11: sp__b$' \
    "restore refuses a revision-11 approval of a publication with one dataset at 10" -- bash "$S/restore.sh" rel1 pub-mix-one
expect 0 'dry-run: would start rel1 with pub-r11' \
    "restore accepts a revision-11 approval of a publication with every dataset at 11" -- bash "$S/restore.sh" rel1 pub-r11

# GRM_EXPECT_REVISION and GRM_EXPECT_DATASETS have no default (R-D3):
# publish.sh, stage.sh and approve-legacy.sh refuse with exit 5, before they
# read anything, when either is unset or not a positive whole number. With
# the old defaults, revision 7's, a revision-11 root that happened to hold 237
# reports (a stopped run, say) passed on revision 7's count.
r=$T/r11-237
mkdir -p "$r/datasets"
rows=''
for i in $(seq -w 1 237); do
    printf '{"id":"sp__d%s","evaluation":{"revision":11}}\n' "$i" >"$r/datasets/sp__d$i.json"
    rows+="${rows:+,}{\"id\":\"sp__d$i\",\"evaluation\":{\"revision\":11}}"
done
printf '{"schema":2,"evaluation_revision":11,"datasets":[%s]}\n' "$rows" >"$r/index.json"
echo '{"corpus":1}' >"$r/corpus.json"
rehash r11-237
no_count() { env -u GRM_EXPECT_DATASETS GRM_EXPECT_REVISION=11 "$@"; }
expect 5 'GRM_EXPECT_DATASETS is unset or empty' \
    "publish dry run refuses revision 11 without GRM_EXPECT_DATASETS, though the root holds 237 reports" -- \
    no_count bash "$S/publish.sh" "$r" "$T/r11-237.sha256" pub-r11-default
expect 5 'GRM_EXPECT_DATASETS is unset or empty' "publish --apply refuses it too" -- \
    no_count bash "$S/publish.sh" "$r" "$T/r11-237.sha256" pub-r11-default --apply
expect 0 '' "the refused publish published nothing" -- \
    test ! -e "$D/published/pub-r11-default" -a ! -e "$D/published/pub-r11-default.partial"
expect 0 'published pub-r11-237 files=239' "publish --apply accepts that root with GRM_EXPECT_DATASETS=237 set" -- \
    env GRM_EXPECT_DATASETS=237 GRM_EXPECT_REVISION=11 bash "$S/publish.sh" "$r" "$T/r11-237.sha256" pub-r11-237 --apply
expect 5 'GRM_EXPECT_DATASETS is unset or empty' "stage refuses revision 11 without GRM_EXPECT_DATASETS" -- \
    no_count env GRM_SKIP_BROWSER=1 bash "$S/stage.sh" rel2 pub-r11-237 "$STAGE_PORT"
# stage.sh refuses the unset count before it reads anything (REVIEW-DEPLOY
# L3): with a publication that does not exist, the count is still what it
# reports, not the missing publication.
expect 5 'GRM_EXPECT_DATASETS is unset or empty' \
    "stage refuses the unset count before it reads the release or the publication (a missing one here)" -- \
    no_count env GRM_SKIP_BROWSER=1 bash "$S/stage.sh" rel2 pub-missing "$STAGE_PORT"
expect 0 'dry-run: would serve rel2 with pub-r11-237' "stage dry run accepts it with GRM_EXPECT_DATASETS=237 set" -- \
    env GRM_SKIP_BROWSER=1 GRM_EXPECT_DATASETS=237 GRM_EXPECT_REVISION=11 bash "$S/stage.sh" rel2 pub-r11-237 "$STAGE_PORT"
expect 5 'GRM_EXPECT_DATASETS is unset or empty' "approve-legacy refuses revision 11 without GRM_EXPECT_DATASETS" -- \
    no_count bash "$S/approve-legacy.sh" rel1 pub-r11-237
# Revision 7, the old default, is no exception, nor is a count without a
# revision; and each must be a positive whole number in ASCII digits.
# refused SCRIPT ENV-ARGS...: a dry run of SCRIPT on pub1 under env ENV-ARGS
refused() {
    local script=$1
    shift
    case $script in
        publish) env "$@" bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub-refused ;;
        stage) env "$@" GRM_SKIP_BROWSER=1 bash "$S/stage.sh" rel2 pub1 "$STAGE_PORT" ;;
        approve-legacy) env "$@" bash "$S/approve-legacy.sh" rel1 pub1 ;;
    esac
}
bad_values=(0 07 abc "$(printf '7\r')" "$(printf '\357\274\227')" 1000000000000000000)
refused_before=$(snapshot)
bad_names=(zero 'a leading zero' letters 'a trailing CR' 'a full-width digit' '19 digits')
for script in publish stage approve-legacy; do
    expect 5 'GRM_EXPECT_REVISION is unset or empty' "$script refuses with no GRM_EXPECT_* set (no revision-7 default)" -- \
        refused "$script" -u GRM_EXPECT_SCHEMA -u GRM_EXPECT_REVISION -u GRM_EXPECT_DATASETS
    expect 5 'GRM_EXPECT_DATASETS is unset or empty' "$script refuses GRM_EXPECT_REVISION=7 without a count (no 237 default)" -- \
        refused "$script" -u GRM_EXPECT_DATASETS GRM_EXPECT_REVISION=7
    expect 5 'GRM_EXPECT_REVISION is unset or empty' "$script refuses a count without a revision" -- \
        refused "$script" -u GRM_EXPECT_REVISION GRM_EXPECT_DATASETS=2
    expect 5 'GRM_EXPECT_REVISION is unset or empty' "$script refuses an empty GRM_EXPECT_REVISION" -- \
        refused "$script" GRM_EXPECT_REVISION= GRM_EXPECT_DATASETS=2
    for k in "${!bad_values[@]}"; do
        expect 5 'GRM_EXPECT_REVISION must be a positive whole number' "$script refuses GRM_EXPECT_REVISION with ${bad_names[k]}" -- \
            refused "$script" "GRM_EXPECT_REVISION=${bad_values[k]}" GRM_EXPECT_DATASETS=2
    done
    expect 5 'GRM_EXPECT_DATASETS must be a positive whole number' "$script refuses GRM_EXPECT_DATASETS=0" -- \
        refused "$script" GRM_EXPECT_REVISION=7 GRM_EXPECT_DATASETS=0
    expect 5 'GRM_EXPECT_DATASETS must be a positive whole number' "$script refuses GRM_EXPECT_DATASETS=-2" -- \
        refused "$script" GRM_EXPECT_REVISION=7 GRM_EXPECT_DATASETS=-2
done
expect 0 '' "the refused runs changed nothing under the root" -- test "$(snapshot)" = "$refused_before"
expect 0 'GRM_EXPECT_\* is not used here' "restore, which checks the approval, is not refused for revision 11 without a count" -- \
    no_count bash "$S/restore.sh" rel1 pub1

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
# A release that carries its commit's browser check, a stale check in src-rel2
# that stage must not run, and a node stub that would record any run.
broken_release rel-browser
echo '// browser check' >"$D/release/rel-browser/browser-check.cjs"
manifest "$D/release/rel-browser" >"$D/release/rel-browser/MANIFEST.sha256"
mkdir -p "$D/src-rel2/web/scripts" "$T/nodebin" "$T/playwright"
echo '// stale browser check' >"$D/src-rel2/web/scripts/browser-check.cjs"
printf '#!/bin/sh\necho node >>"$GRM_TEST_MARKER"\nexit 0\n' >"$T/nodebin/node"
chmod +x "$T/nodebin/node"

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
browser_stage() { env PATH="$T/nodebin:$PATH" GRM_PLAYWRIGHT_CORE="$T/playwright" bash "$S/stage.sh" "$1" pub1 "$STAGE_PORT"; }
expect 2 'rel2 has no browser-check.cjs in its manifest' \
    "stage refuses a browser check its release does not carry, even with one in src-LABEL" -- browser_stage rel2
expect 0 'browser=run' "stage dry run accepts the browser check its release carries" -- browser_stage rel-browser
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
# browser=pass binds the SHA-256 of the release's own browser check; a stage
# from before that check was bound wrote browser=pass without it.
tamper "a browser=pass approval without browser_check (an older stage) is refused" 'browser=pass without browser_check.*re-stage rel1 with pub2' \
    -e 's/browser=skipped/browser=pass/'
tamper "a browser=pass approval whose browser_check the release does not carry is refused" "browser_check does not match release rel1's" \
    -e "s/browser=skipped/browser=pass browser_check=$(printf a%.0s {1..64})/"
tamper "browser_check without browser=pass is refused" 'browser_check without browser=pass' \
    -e "s/browser=skipped/browser=skipped browser_check=$(printf a%.0s {1..64})/"
approval rel-browser pub1 |
    sed "s/browser=skipped/browser=pass browser_check=$(digest "$D/release/rel-browser/browser-check.cjs")/" >"$D/go/rel-browser--pub1"
expect 0 'dry-run: would start rel-browser with pub1' "a browser=pass approval binding the release's browser check is accepted" -- \
    bash "$S/restore.sh" rel-browser pub1
rm "$D/go/rel-browser--pub1"
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
# head on the same input is killed by SIGPIPE (exit 141). GitHub's runner starts
# steps with SIGPIPE ignored, which a shell can't undo, so the writers got EPIPE
# and sed exited 4; env --default-signal gives the control the default action.
big_state() { state rel2 pub1; seq 1 400000 | sed 's/^/pad=/'; }
big=$(big_state | lines 1 9)
expect 0 '' "a drain-safe cut of a 2.7 MB pipeline succeeds under pipefail" -- test "$big" = "$(state rel2 pub1)"
expect 141 '' "head on the same pipeline fails under pipefail with SIGPIPE (why lines is used)" -- \
    env --default-signal=PIPE bash -c 'set -o pipefail; seq 1 400000 | sed "s/^/pad=/" | head -n 9 >/dev/null'
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

# The Vast boot hook, ops/vast/boot-active.sh (the PR 64 review, M1 and L3),
# run with PATH shims and a stub restore.sh in a directory of its own. The rm
# shim logs every call and never passes /dev/null on to the real rm, so the
# old hook's rm of /dev/null is caught without removing anything. Each run
# prints its exit code, the rm calls on /dev/null, the restores of the active
# pair, the fallback restores, the BOOT-FALLBACK lines and the temporary files
# left behind.
BOOT_HOOK=$REPO/ops/vast/boot-active.sh
cat >"$T/boot-hook.sh" <<'SH'
set -u
hook=$1 dir=$2 mktemp_mode=$3 active=$4 dry=${5:-}
mkdir -p "$dir/bin" "$dir/scripts/deployment" "$dir/root/state" "$dir/tmp"
real_rm=$(command -v rm)
cat >"$dir/bin/rm" <<EOF2
#!/bin/sh
echo "rm \$*" >>"$dir/rm.log"
for a in "\$@"; do [ "\$a" != /dev/null ] || exit 0; done
exec $real_rm "\$@"
EOF2
if [ "$mktemp_mode" = fail ]; then printf '#!/bin/sh\nexit 1\n' >"$dir/bin/mktemp"; fi
cat >"$dir/scripts/deployment/restore.sh" <<EOF2
n=\$(cat "$dir/active.count" 2>/dev/null || echo 0)
case "\$*" in
    "" | --apply)
        n=\$((n + 1)); echo "\$n" >"$dir/active.count"
        case $active in
            pass) exit 0 ;;
            once) [ "\$n" -ge 2 ] && exit 0; echo "2026-10-30T00:00:00Z FAIL: stub transient failure"; exit 1 ;;
            *) echo "2026-10-30T00:00:00Z FAIL: stub failure"; exit 1 ;;
        esac ;;
    "ui-r8-690815a rev8-273-20261029"*) echo fallback >>"$dir/fallback.log"; exit 0 ;;
esac
exit 9
EOF2
chmod +x "$dir"/bin/*
printf 'label=ui-r9-test\npublication=rev9-test\n' >"$dir/root/state/active"
env PATH="$dir/bin:$PATH" TMPDIR="$dir/tmp" BOOT_DEPLOY_ROOT="$dir/root" BOOT_SCRIPTS="$dir/scripts/deployment" \
    BOOT_RETRY_WAIT=0 bash "$hook" $dry >"$dir/out" 2>&1
rc=$?
cnt() { if [ -f "$2" ]; then grep -c -- "$1" "$2"; else echo 0; fi; }
echo "exit=$rc rm_devnull=$(cnt ' /dev/null' "$dir/rm.log") active_tries=$(cat "$dir/active.count" 2>/dev/null || echo 0) fallback=$(cnt . "$dir/fallback.log") marker=$(cnt . "$dir/root/BOOT-FALLBACK") tmp_left=$(ls -A "$dir/tmp" | wc -l)"
SH
# The hook must be in the tree: both suites fail, rather than skip, when an
# assertion cannot run (README.md, REVIEW-DS2 I4).
expect 0 '' "the boot hook, ops/vast/boot-active.sh, is in this tree" -- test -f "$BOOT_HOOK"
expect 0 '^exit=3 rm_devnull=0 active_tries=2 fallback=1 marker=2 tmp_left=0$' \
    "boot hook: with mktemp failing it never removes /dev/null, and falls back after a retry" -- \
    timeout -k 5 60 bash "$T/boot-hook.sh" "$BOOT_HOOK" "$T/boot1" fail fail
expect 0 '^exit=0 rm_devnull=0 active_tries=1 fallback=0 marker=0 tmp_left=0$' \
    "boot hook: with mktemp failing and the active pair healthy it never removes /dev/null" -- \
    timeout -k 5 60 bash "$T/boot-hook.sh" "$BOOT_HOOK" "$T/boot2" fail pass
expect 0 '^exit=0 rm_devnull=0 active_tries=2 fallback=0 marker=0 tmp_left=0$' \
    "boot hook: a transient failure of the active pair is retried once, with no fallback and no marker" -- \
    timeout -k 5 60 bash "$T/boot-hook.sh" "$BOOT_HOOK" "$T/boot3" ok once
expect 0 '^exit=3 rm_devnull=0 active_tries=2 fallback=1 marker=2 tmp_left=0$' \
    "boot hook: an active pair that fails twice falls back, writes BOOT-FALLBACK and removes its temporary file" -- \
    timeout -k 5 60 bash "$T/boot-hook.sh" "$BOOT_HOOK" "$T/boot4" ok fail
expect 0 '^exit=3 rm_devnull=0 active_tries=2 fallback=1 marker=0 tmp_left=0$' \
    "boot hook: the dry run retries without waiting, falls back and writes no marker" -- \
    timeout -k 5 60 bash "$T/boot-hook.sh" "$BOOT_HOOK" "$T/boot5" ok fail --dry-run

# The live port as smoke.sh reads it (REVIEW-DS4 L1): a leading zero, which
# grm and curl read as decimal, passed the digits check and then failed every
# health check. It is refused before anything runs, and so is such a stage port.
expect 2 'GRM_PORT must be a port number' "publish refuses a GRM_PORT with a leading zero before reading anything" -- \
    env GRM_PORT="0$LIVE_PORT" bash "$S/publish.sh" "$src" "$T/SOURCE.sha256" pub-port
expect 2 'invalid port 0641' "stage refuses a spare port with a leading zero" -- \
    bash "$S/stage.sh" rel1 pub1 0641

expect 1 '' "no stub grm, cargo, npm or node ever ran" -- test -e "$GRM_TEST_MARKER"
EXPECTED_ASSERTIONS=$((282 - skipped))  # less the root-only skip
if [ "$((pass + fail))" != "$EXPECTED_ASSERTIONS" ]; then
    fail=$((fail + 1))
    echo "FAIL dry-run suite ran $((pass + fail - 1)) assertions, expected $EXPECTED_ASSERTIONS"
fi
finished=1
echo "dry-run tests: $pass passed, $fail failed"
if ((fail)); then exit 1; fi
