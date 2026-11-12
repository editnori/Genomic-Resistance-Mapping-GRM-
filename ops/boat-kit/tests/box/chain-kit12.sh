#!/usr/bin/env bash
# KIT12 test box chain, in one detached script with the kit's chain-lib.sh (the chain exits 1 when any step fails).
# The harness is KIT11's (KIT10's harness at /home/user/k10, so setup-mirror.sh and the path-pinned tests run
# unchanged), with KIT12's drive as cand/drive.sh and KIT12's h1-box.sh and r11-box.sh in kit/; KIT11's drive and
# runner are k11final/ and KIT11's kit is k11kit/. Steps: KIT12's two tests; each on KIT11's files as a control,
# which must fail with exactly the pinned DIFFERS ids; round 3: each on KIT12 round 2's files (4ae8d7a: its drive as
# k12r2/drive.sh, its h1-box.sh in k12r2kit/) as a control, which must fail with exactly round 3's new ids; KIT11's and
# KIT10's runner, H-1 box, H-1 drive, L2 drive and stub regression tests on KIT12's files (regression); then KIT12's
# mutants.
set -u
H=/home/user/k10 V=/var/tmp/k12
CHAIN_LOG=/home/user/k12-logs/chain
[ -e /home/user/k12-logs/.started ] && { echo "already started"; exit 0; }
mkdir -p /home/user/k12-logs "$V"; touch /home/user/k12-logs/.started
. "$H/kit/chain-lib.sh"
cd "$H"
# KIT11's kit has no addon role: its checkpoint (B4) and a run box's arguments (B11) hold, everything else differs.
# Round 3: it also lacks the species union and prune, so it differs on round 3's ids (R3_BOX) as round 2's kit does
R3_BOX="U1-union-is-the-spec-s-species-only U2-addon-union-is-the-spec-s-species-only U3-spec-without-species-2-no-run
 U4-species-not-a-slug-2-no-run U5-verify-leaves-no-scratch-work P1-prune-removes-the-root P1-prune-receipt
 P4-prune-of-a-set-2 P5-prune-outside-the-home-2 P9-prune-of-a-root-the-sets-are-not-linked-from-2"
K11_BOX="B1-addon-union-passes-the-reference B2-addon-verify-passes-the-reference B3-addon-run-passes-the-reference
 B5-reference-pin-mismatch-2-no-run B6-reference-file-changed-2 B7-a-second-part-2 B8-no-reference-name-2
 B9-reference-name-a-path-2 B10-unfinished-main-part-2 B12-run-rechecks-the-set D1-addon-verify-goes-to-h1-box
 D2-addon-start-runs-h1-box-run-to-done D3-addon-reads-h1-run-s-pause D4-bad-role-names-addon $R3_BOX"
# round 3's drive ids: KIT12 round 2's drive has no DATA_SNAPSHOT_PRUNE (A17, no prune unless set, holds there)
R3_DRIVE="A15-pilot-prunes-after-its-fasta-set A16-addon-prunes-after-its-fasta-set
 A18-prune-refused-81-closed-before-verify A19-bad-prune-settings-2-no-box A19-each-names-its-setting"
# KIT11's drive refuses kind addon with 2 before any box, so only A6's "2, no box" holds, and A17 as above
K11_DRIVE="A1-addon-waits-for-its-main-part A2-addon-runs A2-a-run-box-s-data-then-the-main-part
 A2-records-by-push-data-data-by-stage-local A2-the-main-part-s-records-are-sent A2-records-hold-the-main-part-s-work
 A2-pinned-in-notes-before-the-box A2-stage-env-names-the-part-and-pin A3-a-run-box-names-no-reference
 A4-addon-box-from-the-data-snapshot-with-the-binary-gate A5-combine-takes-the-addon-as-a-part A6-each-names-its-cell
 A7-other-flags-than-the-main-part-2-no-box A8-main-part-without-flags-2-no-box A9-addon-waits-for-the-reference-binary
 A10-records-refused-97-closed A11-push-data-mode-addon A12-main-part-done-but-not-pulled-79
 A13-main-part-pulled-but-not-done-79 A14-rerun-builds-the-set-afresh $R3_DRIVE"
info() { uname -a; nproc; free -g | head -n 2; grep -m1 'model name' /proc/cpuinfo; df -h /home/user | tail -n 1; python3 -V; id; }
syntax() {
  local f n=0
  for f in cand/grm-boat.sh cand/drive.sh k11final/*.sh k11kit/*.sh kit/*.sh tests/*.sh box/*.sh; do bash -n "$f" || { echo "FAIL $f"; return 1; }; n=$((n + 1)); done
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
sums() { sha256sum cand/* k11final/* k12r2/* k12r2kit/h1-box.sh kit/*.sh kit/*.py tests/* box/*; }
must info info
must sums sums
must mirror bash box/setup-mirror.sh
step syntax syntax
step box-k12 bash tests/test-h1-box-kit12.sh "$V/box-k12" kit
step box-k12-k11-control control "$K11_BOX" bash tests/test-h1-box-kit12.sh "$V/box-k12-k11" k11kit
step drive-k12 bash tests/test-drive-h1-kit12.sh "$V/drive-k12" cand/drive.sh kit
step drive-k12-k11-control control "$K11_DRIVE" bash tests/test-drive-h1-kit12.sh "$V/drive-k12-k11" k11final/drive.sh kit
step box-k12-r2-control control "$R3_BOX" bash tests/test-h1-box-kit12.sh "$V/box-k12-r2" k12r2kit
step drive-k12-r2-control control "$R3_DRIVE" bash tests/test-drive-h1-kit12.sh "$V/drive-k12-r2" k12r2/drive.sh kit
step runner-k11 bash tests/test-runner-kit11.sh cand/grm-boat.sh "$V/runner-k11"
step drive-k11 bash tests/test-drive-h1-kit11.sh "$V/drive-k11" cand/drive.sh kit
step runner-k10 bash tests/test-runner-kit10.sh cand/grm-boat.sh "$V/runner-k10" kit/chain-lib.sh
step h1box bash tests/test-h1-box.sh "$V/h1box" kit
step drive-h1 bash tests/test-drive-h1.sh "$V/drive-h1" cand/drive.sh kit drafts
step drive-l2 bash tests/test-drive-l2-kit10.sh "$V/drive-l2" cand/drive.sh
step stub-regress bash tests/stub-regress-kit10.sh "$V/regress" cand/drive.sh
step mutants python3 -I -B box/mutants-kit12.py "$H" "$V/mut" /home/user/k12-logs/mutants.txt "$(( $(nproc) / 2 ))"
finish
