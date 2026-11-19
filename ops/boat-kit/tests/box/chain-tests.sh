#!/usr/bin/env bash
# KIT10 round 2 test box chain: every gate but the H-1 dry run, in one detached script, with the kit's own
# chain-lib.sh (the chain exits 1 when any step fails). Round 2 (REVIEW-KIT10 L5): a control passes only when its test
# fails with exactly the DIFFERS ids pinned here, so a control that fails for some other reason doesn't count. The
# restore step and its live-v8 control are skipped this round, since no restore set goes to Boat; mtb-box.sh is checked
# byte-equal to the file run 5's restore receipts cover instead. Round 3 (REVIEW-KIT10-R2 L1): the install test gains
# I14-I18, so round 1's install control pins them too, and install-r2-control runs the test on round 2's install.
set -u
H=/home/user/k10 V=/var/tmp/k10
CHAIN_LOG=/home/user/k10-logs/chain
[ -e /home/user/k10-logs/.started ] && { echo "already started"; exit 0; }
mkdir -p /home/user/k10-logs "$V"; touch /home/user/k10-logs/.started
. "$H/kit/chain-lib.sh"
cd "$H"
MTB_PIN=f67db557c7e642efb01e4a01184cc1103424e97a16e787d2d4a997569ee5f072
# the pinned DIFFERS ids: KIT9's set and round 1's set must fail exactly these checks
K9_RUNNER="B1-fresh-low-beacon-not-read B2-fresh-good-beacon-no-fallback-76 B3-garbled-beacon-irrelevant
 E1-memread-exe-not-run E2-hung-exe-irrelevant G1-slow-start-then-healthy-completes G2-no-first-byte-is-a-stall
 H1-h1-run-outputs-and-master-copy-pass S1-parallel-push-data-same-lane S2-both-boxes-checked S3b-lane-folder-kept
 S4-two-calls-race-on-the-lane-folder U1-unreadable-size-stops U2-before-the-try-timeout V0-runner-pins-13-names
 V1-whole-pinned-vocabulary-passes V4-partial-refused W4-tables-sigdrug-name-refused
 W5-sigdrug-name-in-other-case-refused W7-drive-combine-root-to-grm-combine-passes X1-shorter-refused-unchanged
 X2-longer-set X2b-end-line-is-the-new-end X3-noted X4-no-archiveAfter-refused X5-bad-ttl-2-and-moved-back-86
 X6-extend-that-changed-nothing-86 X7-extend-600s-short-86 X8-extend-60s-short-within-120s-0 X9-moved-back-within-120s-86"
R1_RUNNER="C5-no-step-fails S3b-lane-folder-kept S4-two-calls-race-on-the-lane-folder
 W1-h1work-sigdrug-tsv-elsewhere-refused W2-h1work-under-a-prefix-elsewhere-refused W3-master-manifest-at-depth-3-refused
 W4-tables-sigdrug-name-refused W5-sigdrug-name-in-other-case-refused W6-sigdrug-name-in-a-part-refused-on-combine
 W8-same-root-to-another-destination-refused W9-master-manifest-at-the-top-elsewhere-refused
 W10-h1work-outside-parts-refused-on-combine W11-master-manifest-in-a-part-refused-on-combine X2b-end-line-is-the-new-end
 X6-extend-that-changed-nothing-86 X7-extend-600s-short-86"
R1_H1BOX="D4-guard-pause-reads-paused M3-marker-that-cannot-move-2 R4-run-at-usable-free V7-estimate-at-usable-free
 V8-free-below-the-pause-floor-4"
R1_DRIVE_H1="H4d-combine-parts-with-empty-flags-2-no-box H14b-deadline-is-the-end-the-box-got
 H14c-extend-without-end-86-closed H17-guard-pause-stops-81-with-its-reason H18-git-status-failing-2-no-box
 H19-stage-table-placeholder-2-no-box H21-drafted-stage-table-refused-2-no-box"
K9_DRIVE_L2="L1-stale-gone-from-copy2 L1-copy2-is-out L2-no-done-mark-without-outputs L3-read-only-resume-set-rc
 L4-gone-box-not-86"
R1_INSTALL="I2-install-md-rollback-restores-and-checks I3-install-after-the-rollback
 I7-read-only-target-dir-stops-before-any-change I8-step-2-failure-rolls-back I8b-rerun-after-the-rollback
 I9-failure-after-mktemp-leaves-no-temp I10-last-move-fails-all-seven-undone I12-a-command-naming-drive-sh-does-not-stop-it
 I14-term-in-step-2-rolls-back I15-hup-in-step-2-rolls-back I16-int-to-the-group-in-step-2-rolls-back
 I16b-rerun-after-an-interrupted-install I17-term-just-after-a-move-rolls-that-file-back
 I18-install-md-rollback-removes-left-temp-files"
R2_INSTALL="I14-term-in-step-2-rolls-back I15-hup-in-step-2-rolls-back I16-int-to-the-group-in-step-2-rolls-back
 I16b-rerun-after-an-interrupted-install I17-term-just-after-a-move-rolls-that-file-back
 I18-install-md-rollback-removes-left-temp-files"
info() { uname -a; nproc; free -g | head -n 2; grep -m1 'model name' /proc/cpuinfo; df -h /home/user | tail -n 1; python3 -V; rsync --version | head -n 1; ls -l /usr/bin/time; id; grep -E '^(SigIgn|CapAmb)' /proc/$$/status; }
syntax() {
  local f n=0
  for f in cand/grm-boat.sh cand/drive.sh kit/*.sh tests/*.sh tests/*.next box/*.sh handback/install-kit10.sh \
           handback/next/*.next handback/next/kit/*.sh.next drafts/run-h1.conf r1/cand/*.sh r1/kit/*.sh r1/handback/install-kit10.sh \
           r2/handback/install-kit10.sh; do
    bash -n "$f" || { echo "FAIL $f"; return 1; }; n=$((n + 1))
  done
  for f in kit/*.py box/*.py handback/next/kit/*.py.next r1/kit/*.py; do
    python3 -I -B -c 'import ast, sys; ast.parse(open(sys.argv[1]).read())' "$f" || { echo "FAIL $f"; return 1; }; n=$((n + 1))
  done
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
mtbpin() { # kit, handback and round 1 hold the same mtb-box.sh, the one run 5's restore receipts cover
  local f s ok=0
  for f in kit/mtb-box.sh handback/next/kit/mtb-box.sh.next r1/kit/mtb-box.sh; do
    s=$(sha256sum < "$f" | cut -d' ' -f1); echo "$s  $f"; [ "$s" = "$MTB_PIN" ] || ok=1
  done
  [ $ok = 0 ] && echo "mtb-box.sh byte-equal to $MTB_PIN"; return $ok
}
sums() { sha256sum cand/* kit/*.sh kit/*.py tests/* box/* drafts/* handback/install-kit10.sh handback/next/*.next handback/next/kit/*; }
noise() { # REVIEW-KIT10 L5, a receipt (always rc 0): G1 and P8 (20 kB/s against the 10 kB/s floor) failing in a runner
  # test whose mutant targets neither is load noise; the runner and chain mutants ran half the vCPUs at a time
  local id want n=0 r=0 o
  mkdir -p /home/user/k10-logs/mut
  for o in "$V"/mut/*/test.out; do id=${o%/test.out}; id=${id##*/}; cp "$o" "/home/user/k10-logs/mut/$id.out"; done
  while read -r id _ _ _ want _; do
    o=$V/mut/$id/test.out; grep -q '^RESULT G1-' "$o" 2> /dev/null || continue
    case $want in G1-*|P8-*) continue;; esac
    r=$((r + 1)); grep -E '^RESULT (G1|P8)-' "$o" | grep ' DIFFERS' | sed "s/^/$id: /" | grep . && n=$((n + 1))
  done < <(grep -E ' (KILLED|WRONG-ID|NO-LINE|SURVIVED|NOT-APPLIED) ' /home/user/k10-logs/mutants.txt)
  echo "G1/P8 noise: $n of $r runner-test mutants that target neither failed G1 or P8"
}
must info info
must sums sums
must mirror bash box/setup-mirror.sh
step syntax syntax
step mtb-box-pin mtbpin
step runner bash tests/test-runner-kit10.sh cand/grm-boat.sh "$V/runner" kit/chain-lib.sh
step runner-kit9-control control "$K9_RUNNER" bash tests/test-runner-kit10.sh kit9/grm-boat.sh.next "$V/runner-kit9" kit/chain-lib.sh
step runner-r1-control control "$R1_RUNNER" bash tests/test-runner-kit10.sh r1/cand/grm-boat.sh "$V/runner-r1" r1/kit/chain-lib.sh
step h1box bash tests/test-h1-box.sh "$V/h1box" kit
step h1box-r1-control control "$R1_H1BOX" bash tests/test-h1-box.sh "$V/h1box-r1" r1/kit
step drive-h1 bash tests/test-drive-h1.sh "$V/drive-h1" cand/drive.sh kit drafts
step drive-h1-r1-control control "$R1_DRIVE_H1" bash tests/test-drive-h1.sh "$V/drive-h1-r1" r1/cand/drive.sh r1/kit drafts
step drive-l2 bash tests/test-drive-l2-kit10.sh "$V/drive-l2" cand/drive.sh
step drive-l2-kit9-control control "$K9_DRIVE_L2" bash tests/test-drive-l2-kit10.sh "$V/drive-l2-kit9" kit9/drive.sh.next
step stub-regress bash tests/stub-regress-kit10.sh "$V/regress" cand/drive.sh
step install bash tests/test-install-kit10.sh "$V/install" handback "$H"
step install-r1-control control "$R1_INSTALL" bash tests/test-install-kit10.sh "$V/install-r1" r1/handback "$H"
step install-r2-control control "$R2_INSTALL" bash tests/test-install-kit10.sh "$V/install-r2" r2/handback "$H"
step prov bash -c 'cd prov && python3 -I -B -m unittest discover -s . -p test_stage_manifest.py -k test_the_provenance_file_is_never_staged_and_nothing_else_changes -k test_without_the_file_both_agree_and_refusals_are_unchanged -v && cmp stage-manifest.py.next ../kit/stage-manifest.py && echo "kit stage-manifest.py is PROV next"'
step mutants python3 -I -B box/mutants.py "$H" "$V/mut" /home/user/k10-logs/mutants.txt $(( $(nproc) / 2 ))
step g1-p8-noise noise
finish
