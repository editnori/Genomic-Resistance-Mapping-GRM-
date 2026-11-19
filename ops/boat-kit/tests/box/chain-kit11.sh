#!/usr/bin/env bash
# KIT11 test box chain, in one detached script with the kit's chain-lib.sh (the chain exits 1 when any step fails).
# The harness is KIT10's (make-harness.sh's k10.tar.gz, at /home/user/k10 so setup-mirror.sh and the path-pinned
# tests run unchanged), with KIT11's runner and drive as cand/, KIT10's as k10final/, and KIT11's tests and stubs over
# tests/. Steps: KIT11's two tests; each on KIT10's files as a control, which must fail with exactly the pinned
# DIFFERS ids; KIT10's runner, H-1 box, H-1 drive, L2 drive and stub regression tests on KIT11's files (regression);
# then KIT11's mutants.
set -u
H=/home/user/k10 V=/var/tmp/k11
CHAIN_LOG=/home/user/k11-logs/chain
[ -e /home/user/k11-logs/.started ] && { echo "already started"; exit 0; }
mkdir -p /home/user/k11-logs "$V"; touch /home/user/k11-logs/.started
. "$H/kit/chain-lib.sh"
cd "$H"
# KIT10's runner has no stage-local and no GRM_HOST: every check but L7 (wsl mode, unchanged) and T7 and T22 (usage
# refusals, rc 2 from KIT10's usage too) differs
K10_RUNNER="L1-linux-needs-grm-boat L2-unknown-host-2 L3-linux-preflight-reads-this-host L4-linux-memory-under-the-floor-76
 L5-bad-floor-2 L6-linux-push-data-stages-here L6-box-checked L8-linux-pull T1-stage-local-passes T1-exact-set-and-manifest
 T1-hard-links-not-copies T1-sent-no-tar T1-noted T2-wrong-pin-97-no-box-call T3-denied-segment-97-no-box-call
 T4-receipt-goes-through-push-data-97 T5-grm-combine-refused-97 T6-denied-source-folder-97
 T8-a-pass-leaves-the-source-as-it-was T9-file-link-in-source-86 T9-names-the-link T10-folder-link-in-source-86
 T11-source-folder-a-link-86 T12-changed-bytes-86-source-kept T13-missing-source-file-86 T14-dest-not-empty-86-kept
 T15-dest-a-link-86 T16-manifest-changed-on-the-way-86 T17-line-twice-97 T18-two-sets-from-one-source
 T19-whole-pinned-vocabulary-links T20-partial-vocabulary-97 T21-real-pins-refuse-other-bytes"
# KIT10's drive has no snapshot mode: push-data mode (D0) and its DATA_ROOT rule (D8) hold, everything else differs
K10_DRIVE="D1-snapshot-pilot-rc D1-three-stage-local-no-push-data D1-from-the-snapshot-folders D1-same-manifests-as-push-data
 D1-pins-passed D1-box-from-the-data-snapshot D1-vocabulary-pin-names-its-file D1-stage-env-pins-match-push-mode
 D1-done-closed D2-combine-still-pushes-its-parts D3-snapshot-mode-is-h1-only D4-snapshot-mode-needs-a-grm-tmpl-template
 D5-missing-vocab-manifest-2-no-box D6-other-vocab-manifest-2-no-box D7-snapshot-folders-under-home-user
 D9-stage-local-refusal-stops-and-closes D10-snapshot-mode-needs-no-data-or-vocab-root"
info() { uname -a; nproc; free -g | head -n 2; grep -m1 'model name' /proc/cpuinfo; df -h /home/user | tail -n 1; python3 -V; id; }
syntax() {
  local f n=0
  for f in cand/grm-boat.sh cand/drive.sh k10final/*.sh kit/*.sh tests/*.sh box/*.sh; do bash -n "$f" || { echo "FAIL $f"; return 1; }; n=$((n + 1)); done
  for f in kit/*.py box/*.py; do python3 -I -B -c 'import ast, sys; ast.parse(open(sys.argv[1]).read())' "$f" || { echo "FAIL $f"; return 1; }; n=$((n + 1)); done
  echo "syntax ok: $n files"
}
control() { # control "PINNED IDS" CMD...: CMD must fail, and its DIFFERS ids must be exactly the pinned ones
  local want got rc o; want=$(printf '%s\n' $1 | LC_ALL=C sort | paste -sd' '); shift
  o=$(mktemp); "$@" > "$o" 2>&1; rc=$?; cat "$o"
  got=$(grep '^RESULT .* DIFFERS' "$o" | awk '{print $2}' | LC_ALL=C sort | paste -sd' '); rm -f "$o"
  echo "control rc=$rc"; echo "control pinned DIFFERS: $want"; echo "control got DIFFERS:    $got"
  if [ "$rc" != 0 ] && [ "$got" = "$want" ]; then echo "CONTROL OK: failed with exactly the pinned ids"; return 0; fi
  echo "CONTROL WRONG: rc=$rc, or the DIFFERS ids aren't the pinned ones"; return 1
}
sums() { sha256sum cand/* k10final/* kit/*.sh kit/*.py tests/* box/*; }
must info info
must sums sums
must mirror bash box/setup-mirror.sh
step syntax syntax
step runner-k11 bash tests/test-runner-kit11.sh cand/grm-boat.sh "$V/runner-k11"
step runner-k11-k10-control control "$K10_RUNNER" bash tests/test-runner-kit11.sh k10final/grm-boat.sh "$V/runner-k11-k10"
step drive-k11 bash tests/test-drive-h1-kit11.sh "$V/drive-k11" cand/drive.sh kit
step drive-k11-k10-control control "$K10_DRIVE" bash tests/test-drive-h1-kit11.sh "$V/drive-k11-k10" k10final/drive.sh kit
step runner-k10 bash tests/test-runner-kit10.sh cand/grm-boat.sh "$V/runner-k10" kit/chain-lib.sh
step h1box bash tests/test-h1-box.sh "$V/h1box" kit
step drive-h1 bash tests/test-drive-h1.sh "$V/drive-h1" cand/drive.sh kit drafts
step drive-l2 bash tests/test-drive-l2-kit10.sh "$V/drive-l2" cand/drive.sh
step stub-regress bash tests/stub-regress-kit10.sh "$V/regress" cand/drive.sh
step mutants python3 -I -B box/mutants-kit11.py "$H" "$V/mut" /home/user/k11-logs/mutants.txt "$(( $(nproc) / 2 ))"
finish
