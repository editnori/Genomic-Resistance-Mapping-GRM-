#!/usr/bin/env bash
# Tests of publish, activate, rollback and status against local Miniflare
# state; no Cloudflare account. usage: release.sh WORK_DIR APP_DIST PUB_A PUB_B
# PUB_A and PUB_B are publication layouts ({report/,SOURCE.sha256}) of two
# synthetic trees at revision 11 with 35 and 30 datasets (synthetic-report.mts
# --shape full and wide). WORK_DIR must not exist. A Worker is started on the
# test state for the served-check cases. Prints one line per case and exits 1
# if any fails.
set -uo pipefail
W=$1 DIST=$2 A=$3 B=$4
S=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../scripts" && pwd)
WORKER=$(cd -- "$S/../worker" && pwd)
[ ! -e "$W" ] || { echo "$W exists" >&2; exit 2; }
mkdir -p "$W"
export GRM_SITE_TARGET=local GRM_SITE_PERSIST=$W/state GRM_SITE_LOG=$W/site.log
unset GRM_SITE_URL GRM_EXPECT_REVISION GRM_EXPECT_DATASETS GRM_EXPECT_SCHEMA
PASS=0 FAIL=0
C=0000000000000000000000000000000000000001
ok() { echo "pass $1"; PASS=$((PASS + 1)); }
bad() { echo "FAIL $1: $2"; FAIL=$((FAIL + 1)); }
# expect NAME RC CMD...: runs CMD, output to $W/NAME.out, and wants exit RC.
expect() {
    local name=$1 want=$2 rc=0
    shift 2
    "$@" >"$W/$name.out" 2>&1 || rc=$?
    if [ "$rc" = "$want" ]; then ok "$name (rc $rc)"; else bad "$name" "rc $rc, wanted $want; $(tail -n 2 "$W/$name.out" | tr '\n' ' ')"; fi
}
state() { "$S/../worker/node_modules/.bin/wrangler" r2 object get "grm-site/state/$1" --pipe --local --persist-to "$W/state" 2>/dev/null; }
pair() { state "$1" | python3 -I -c 'import json,sys; p=json.load(sys.stdin); print(p["app"], p["publication"])' 2>/dev/null || echo none; }
want_pairs() {
    local a p
    a=$(pair active) p=$(pair previous)
    if [ "$a" = "$2" ] && [ "$p" = "$3" ]; then ok "$1 (active $a, previous $p)"; else bad "$1" "active '$a' previous '$p', wanted '$2' '$3'"; fi
}
pub() { GRM_EXPECT_REVISION=11 GRM_EXPECT_DATASETS=$2 bash "$S/publish.sh" "$1/report" "$1/SOURCE.sha256" "${@:3}"; }

# Publishing
expect publish-needs-datasets 2 env GRM_EXPECT_REVISION=11 bash "$S/publish.sh" "$A/report" "$A/SOURCE.sha256" syn-a --apply
expect publish-needs-revision 2 env GRM_EXPECT_DATASETS=35 bash "$S/publish.sh" "$A/report" "$A/SOURCE.sha256" syn-a --apply
expect publish-wrong-count 1 pub "$A" 34 syn-a --apply
expect publish-wrong-revision 1 env GRM_EXPECT_REVISION=9 GRM_EXPECT_DATASETS=35 bash "$S/publish.sh" "$A/report" "$A/SOURCE.sha256" syn-a --apply
cp -a "$A" "$W/tampered"
chmod -R u+w "$W/tampered"
printf ' ' >>"$W/tampered/report/corpus.json"
expect publish-tampered-source 1 pub "$W/tampered" 35 syn-t --apply
expect publish-bad-name 2 pub "$A" 35 'syn--a' --apply
# A dry run logs nothing and claims nothing (publish-a below would refuse).
touch "$W/site.log"
before=$(wc -l <"$W/site.log")
expect publish-dry-run 0 pub "$A" 35 syn-a
[ "$(wc -l <"$W/site.log")" = "$before" ] && ok "publish-dry-run-logs-nothing" || bad "publish-dry-run-logs-nothing" "log grew"
expect publish-a 0 pub "$A" 35 syn-a --apply
expect publish-b 0 pub "$B" 30 syn-b --apply
expect republish-refused 2 pub "$A" 35 syn-a --apply
expect republish-other-tree-refused 2 pub "$B" 30 syn-a --apply
expect publish-app 0 bash "$S/publish-app.sh" "$DIST" ui-t "$C" --apply
expect publish-app-again-refused 2 bash "$S/publish-app.sh" "$DIST" ui-t "$C" --apply
expect publish-app-bad-commit 2 bash "$S/publish-app.sh" "$DIST" ui-u abc --apply

# A claimed but unsealed publication (an interrupted upload) and a sealed one
# whose object was changed afterwards, both refused by activate.
printf '{}\n' >"$W/claim"
"$WORKER/node_modules/.bin/wrangler" r2 object put grm-site/pub/syn-half/CLAIM --file "$W/claim" --local --persist-to "$W/state" >/dev/null 2>&1
expect activate-unsealed 1 bash "$S/activate.sh" ui-t syn-half --apply
expect publish-over-claim-refused 2 pub "$A" 35 syn-half --apply
expect publish-c 0 pub "$A" 35 syn-c --apply
printf '{"tampered":true}\n' >"$W/evil"
"$WORKER/node_modules/.bin/wrangler" r2 object put grm-site/pub/syn-c/report/corpus.json --file "$W/evil" --local --persist-to "$W/state" >/dev/null 2>&1
expect activate-tampered-object 1 bash "$S/activate.sh" ui-t syn-c --apply
expect activate-missing-app 1 bash "$S/activate.sh" ui-none syn-a --apply
expect rollback-without-previous 2 bash "$S/rollback.sh" --apply
want_pairs nothing-active-yet none none

# Activation and rollback
expect activate-dry-run 0 bash "$S/activate.sh" ui-t syn-a
want_pairs dry-run-changes-nothing none none
expect activate-a 0 bash "$S/activate.sh" ui-t syn-a --apply
want_pairs first-activation "ui-t syn-a" none
expect activate-same-refused 2 bash "$S/activate.sh" ui-t syn-a --apply
expect activate-b 0 bash "$S/activate.sh" ui-t syn-b --apply
want_pairs second-activation "ui-t syn-b" "ui-t syn-a"
expect rollback-dry-run 0 bash "$S/rollback.sh"
grep -q 'rollback target from state/previous: ui-t with syn-a' "$W/rollback-dry-run.out" && ok rollback-dry-run-names-previous || bad rollback-dry-run-names-previous "$(tail -n 1 "$W/rollback-dry-run.out")"
want_pairs rollback-dry-run-changes-nothing "ui-t syn-b" "ui-t syn-a"
expect rollback 0 bash "$S/rollback.sh" --apply
want_pairs after-rollback "ui-t syn-a" "ui-t syn-b"
expect rollback-forward 0 bash "$S/rollback.sh" --apply
want_pairs after-second-rollback "ui-t syn-b" "ui-t syn-a"
expect status 0 bash "$S/status.sh"
grep -q '^active   {"app":"ui-t","publication":"syn-b"' "$W/status.out" && grep -q '^previous {"app":"ui-t","publication":"syn-a"' "$W/status.out" &&
    ok status-prints-both || bad status-prints-both "$(head -n 2 "$W/status.out" | cut -c1-80)"

# The served check: a failed check restores both state objects byte for byte.
port=$((20000 + RANDOM % 20000))
(cd "$WORKER" && exec setsid "$WORKER/node_modules/.bin/wrangler" dev --local --persist-to "$W/state" --ip 127.0.0.1 --port "$port" --show-interactive-dev-session=false) >"$W/dev.log" 2>&1 </dev/null &
dev=$!
trap 'kill -- -$dev 2>/dev/null' EXIT
for _ in $(seq 1 60); do curl -s -o /dev/null "http://127.0.0.1:$port/" && break; sleep 1; done
state active >"$W/active.before"
state previous >"$W/previous.before"
expect activate-unpublished 1 env GRM_SITE_URL=http://127.0.0.1:9 bash "$S/activate.sh" ui-t syn-c2 --apply
expect publish-c2 0 pub "$A" 35 syn-c2 --apply
expect served-check-fails-restores 1 env GRM_SITE_URL=http://127.0.0.1:9 bash "$S/activate.sh" ui-t syn-c2 --apply
state active >"$W/active.after"
state previous >"$W/previous.after"
cmp -s "$W/active.before" "$W/active.after" && cmp -s "$W/previous.before" "$W/previous.after" &&
    ok "failed-check-restored-state-bytes" || bad "failed-check-restored-state-bytes" "state objects differ"
grep -q 'restored state/active = ui-t with syn-b and state/previous' "$W/served-check-fails-restores.out" &&
    ok failed-check-logged || bad failed-check-logged "$(tail -n 1 "$W/served-check-fails-restores.out")"
expect served-check-passes 0 env GRM_SITE_URL="http://127.0.0.1:$port" bash "$S/activate.sh" ui-t syn-c2 --apply
want_pairs after-served-activation "ui-t syn-c2" "ui-t syn-b"
grep -q 'smoke pass' "$W/served-check-passes.out" && ok served-check-ran || bad served-check-ran "no smoke pass line"
served=$(curl -s "http://127.0.0.1:$port/data/index.json" | sha256sum | cut -d' ' -f1)
[ "$served" = "$(sha256sum <"$A/report/index.json" | cut -d' ' -f1)" ] && ok worker-serves-active || bad worker-serves-active "$served"
expect rollback-served 0 env GRM_SITE_URL="http://127.0.0.1:$port" bash "$S/rollback.sh" --apply
served=$(curl -s "http://127.0.0.1:$port/data/index.json" | sha256sum | cut -d' ' -f1)
[ "$served" = "$(sha256sum <"$B/report/index.json" | cut -d' ' -f1)" ] && ok worker-serves-rolled-back || bad worker-serves-rolled-back "$served"

# The log: every --apply step appended, and dry runs appended nothing.
grep -c . "$W/site.log" | xargs -I{} echo "log lines {}"
grep -q 'activated ui-t with syn-b' "$W/site.log" && grep -q 'published syn-a files=' "$W/site.log" &&
    grep -q 'FAIL: grm-site/pub/syn-a already exists' "$W/site.log" && ok log-records-steps || bad log-records-steps "missing lines"
grep -q 'dry-run' "$W/site.log" && bad log-has-no-dry-runs "dry-run line in log" || ok log-has-no-dry-runs
echo "release-tests pass=$PASS fail=$FAIL"
[ "$FAIL" -eq 0 ]
