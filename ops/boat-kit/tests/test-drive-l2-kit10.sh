#!/usr/bin/env bash
# KIT9, DRIVE7 L2: a rerun of a stage starts from a clean out/ and keeps the earlier attempt under attempts/. Runs
# drive.sh against PORT's stub runner and stub boat (no box), with the harness setup taken from PORT's
# drv-gates-v4.sh (everything above its first case). One "RESULT id want got verdict detail" line per check; the last
# line is "FAILS n", and the exit code is n. KIT10 adds REVIEW-KIT9's L1 (COPY2 keeps an earlier attempt's file)
# and L2 (a carry stopped between the outputs and the done mark), then its own L3 (a read-only resume set) and the H-1
# dry run's stale still-up entry (L4). L2 needs a drive that can't read a mode-000 folder: on Boat, a command started
# through boat exec (so the whole chain) holds the ambient capability cap_dac_read_search, which reads it anyway, so
# L2's second call runs under setpriv with the ambient capabilities dropped, and checks that premise first.
# usage: test-drive-l2-kit10.sh WORK DRIVE_SH
set -u
WK=${1:?work dir} DRV_UNDER_TEST=${2:?drive.sh to test}
H=/home/layth/review/grm-port-r11-20261105/dryrun/drv-gates-v4.sh
mkdir -p "$WK"
sed -n '1,/^# ---- signals/p' "$H" | sed '$d' > "$WK/setup.sh"
grep -q '^conf()' "$WK/setup.sh" && grep -q '^drive()' "$WK/setup.sh" || { echo "FAILS 1 (no harness setup in $H)"; exit 1; }
set -- "$WK/h" "$DRV_UNDER_TEST"
. "$WK/setup.sh"
FAILS=0
chk() { res "$@"; [ "$2" = "$3" ] || FAILS=$((FAILS + 1)); }
yn() { if "$@"; then echo yes; else echo no; fi; }

# L2a: a stage whose earlier box left extra files in out/ and ckpt/ (no species finished)
d=$W/l2a; conf "$d" "REF_BOX="; SD=$d/stages/R2
mkdir -p "$SD/out/data/report/datasets" "$SD/out/run" "$SD/ckpt/c1"
echo '{"stale": true}' > "$SD/out/data/report/datasets/stale.json"; echo old > "$SD/out/run/old-status"; echo x > "$SD/ckpt/c1/leftover"
before=$(cd "$SD" && find out ckpt -type f | LC_ALL=C sort | xargs sha256sum)
drive "$d" R2 SR_STATUS=done; rc=$?
A=$(ls -d "$SD"/attempts/* 2> /dev/null | head -n 1)
chk L2a-rc 0 "$rc" ""
chk L2a-stale-gone-from-out no "$(yn test -e "$SD/out/data/report/datasets/stale.json")" ""
chk L2a-stale-not-in-sums 0 "$(grep -c 'stale.json\|old-status' "$SD/out/SHA256SUMS")" ""
chk L2a-stale-not-in-copy2 no "$(yn test -e "$d/copy2/R2/data/report/datasets/stale.json")" ""
chk L2a-earlier-attempt-kept-byte-equal yes "$(yn test "$before" = "$( [ -n "$A" ] && cd "$A" && find out ckpt -type f | LC_ALL=C sort | xargs sha256sum)")" "kept at ${A#$W/}"
chk L2a-new-out-verifies yes "$(yn bash -c "cd '$SD/out' && sha256sum -c --strict --quiet SHA256SUMS")" ""

# L2b: an r11 stage that finished a species, then had an extra file in out/, is rerun: the resume set still carries
# the species, the extra file does not come back
d=$W/l2b; conf "$d" "REF_BOX="; SD=$d/stages/R2
drive "$d" R2 SR_STATUS=running FORCE_TTL_NEAR_AFTER=1 SR_DONE_SPECIES=campylobacter_jejuni; r1=$?
mkdir -p "$SD/out/data/report/datasets"; echo '{"extra": true}' > "$SD/out/data/report/datasets/old.json"
drive "$d" R2 SR_STATUS=done SR_DONE_SPECIES=campylobacter_jejuni; r2=$?
chk L2b-rc "74,0" "$r1,$r2" ""
chk L2b-resume-set-carries-species 1 "$(grep -c 'resume set: species campylobacter_jejuni' "$d/NOTES.md")" ""
chk L2b-resume-manifest-lists-species-results 1 "$(grep -c 'species/campylobacter_jejuni/results/campylobacter_jejuni__x.json' "$SD/RESUME-MANIFEST.sha256")" ""
chk L2b-extra-gone-from-out no "$(yn test -e "$SD/out/data/report/datasets/old.json")" ""
F=$(ls "$SD/ckpt" | grep final | tail -n 1)
other=$(cd "$SD/out" && find . -type f ! -name SHA256SUMS | sed 's#^\./##' | while read -r f; do
  [ -e "$SD/ckpt/$F/$f" ] && continue
  case "$f" in data/species/campylobacter_jejuni/*|run/gates/sp-campylobacter_jejuni.done) continue;; esac; echo "$f"; done | wc -l)
chk L2b-out-is-final-checkpoint-plus-resumed-species 0 "$other" "(files in out/ from neither)"
chk L2b-earlier-attempt-kept 1 "$(ls "$SD/attempts" | wc -l)" "(the first call had no out/ to keep)"
chk L2b-kept-attempt-has-extra yes "$(yn test -e "$(ls -d "$SD"/attempts/* | head -n 1)/out/data/report/datasets/old.json")" ""

# L2c: the reference box is rerun after it fixed BINARY-REF-C, whose ref_reports_dir names out/data: the manifest of
# its g1a reports must still verify there
d=$W/l2c; conf "$d" "REF_BOX=S1" "CANARY_SPECIES=campylobacter_jejuni"; SD=$d/stages/S1
drive "$d" S1 SR_STATUS=running FORCE_TTL_NEAR_AFTER=1 SR_DONE_SPECIES=campylobacter_jejuni; r1=$?
REFK=$(ls "$d"/BINARY-REF-* 2> /dev/null | grep -v '\.sha$\|REPORTS' | head -n 1)
drive "$d" S1 SR_STATUS=done SR_DONE_SPECIES=campylobacter_jejuni; r2=$?
chk L2c-rc "74,0" "$r1,$r2" ""
if [ -n "$REFK" ]; then
  . "$REFK"
  chk L2c-ref-reports-verify-after-rerun yes "$(yn bash -c "cd '$ref_reports_dir' && sha256sum -c --strict --quiet '$ref_reports_manifest'")" "dir ${ref_reports_dir#$W/}"
else chk L2c-ref-reports-verify-after-rerun yes "no BINARY-REF file" ""; fi
# L1 (REVIEW-KIT9 A1): attempt 1 stops (81) and its COPY2 holds a file the rerun's out/ doesn't; after the rerun,
# COPY2 must be exactly out/
d=$W/l1; conf "$d" "REF_BOX="; SD=$d/stages/R2
drive "$d" R2 SR_STATUS=stopped; r1=$?
mkdir -p "$d/copy2/R2/data/report/datasets"; echo '{"stale": true}' > "$d/copy2/R2/data/report/datasets/stale.json"
drive "$d" R2 SR_STATUS=done; r2=$?
chk L1-rc "81,0" "$r1,$r2" ""
chk L1-stale-gone-from-copy2 no "$(yn test -e "$d/copy2/R2/data/report/datasets/stale.json")" ""
chk L1-copy2-is-out yes "$(yn bash -c "diff -rq '$SD/out' '$d/copy2/R2' > /dev/null")" ""   # diff's own lines stay out of RESULT

# L2 (REVIEW-KIT9 A2): the carry of a finished species fails after the mark would have gone first (its outputs
# can't be read); the next out/ must not hold the species' done mark
d=$W/l2d; conf "$d" "REF_BOX="; SD=$d/stages/R2
drive "$d" R2 SR_STATUS=running FORCE_TTL_NEAR_AFTER=1 SR_DONE_SPECIES=campylobacter_jejuni; r1=$?
chmod 000 "$SD/out/data/species/campylobacter_jejuni/results"
chk L2-premise-outputs-unreadable yes "$(yn bash -c "! setpriv --ambient-caps -all ls '$SD/out/data/species/campylobacter_jejuni/results' > /dev/null 2>&1")" \
  "$(grep -E '^Cap(Eff|Amb)' /proc/$$/status | tr '\t\n' '  ')"
env $(envs "$d") SR_STATUS=done SR_DONE_SPECIES=campylobacter_jejuni setpriv --ambient-caps -all timeout 120 bash "$DRV" R2 > "$d/out" 2>&1; r2=$?
chmod -R u+rwx "$SD" 2> /dev/null
chk L2-rc "74,2" "$r1,$r2" "(the second call stops at the carry)"
chk L2-no-done-mark-without-outputs no "$(yn test -e "$SD/out/run/gates/sp-campylobacter_jejuni.done")" ""

# L3 (KIT10): a finished species' outputs are read-only (as a pull of read-only files leaves them). cp -al keeps the
# dir modes, so the second call's resume set is read-only too; the third call must still clear it and resume
d=$W/l3; conf "$d" "REF_BOX="; SD=$d/stages/R2
drive "$d" R2 SR_STATUS=running FORCE_TTL_NEAR_AFTER=1 SR_DONE_SPECIES=campylobacter_jejuni; r1=$?
chmod -R a-w "$SD/out/data/species/campylobacter_jejuni"
drive "$d" R2 SR_STATUS=stopped SR_DONE_SPECIES=campylobacter_jejuni; r2=$?
drive "$d" R2 SR_STATUS=done SR_DONE_SPECIES=campylobacter_jejuni; r3=$?
chmod -R u+rwx "$SD" 2> /dev/null
chk L3-read-only-resume-set-rc "74,81,0" "$r1,$r2,$r3" "$(grep -o 'STOP rc=2: .*' "$d/NOTES.md" | head -n 1)"

# L4 (KIT10, the H-1 dry run): the first close (bail's) fails and the direct delete too, so the box is named still up;
# the exit trap's close then deletes it. The exit is the stage's own 80, not 86 for a box that is gone.
cat > "$W/runner-close-fails-once.sh" << EOF
#!/usr/bin/env bash
# PORT's stub runner, but the first close refuses (86, nothing deleted), as the live runner does when boat usage fails
a=\$1; [ "\$a" = --notes ] && a=\$3
if [ "\$a" = close ] && [ ! -e "\$SR_LOG.close-failed" ]; then
  echo "\$*" >> "\$SR_LOG"; : > "\$SR_LOG.close-failed"; echo "grm-boat: close: boat usage failed" >&2; exit 86
fi
exec bash $ME/drv-stub-runner-v4.sh "\$@"
EOF
chmod +x "$W/runner-close-fails-once.sh"   # drive.sh runs the runner itself, not through bash
d=$W/l4; conf "$d" "REF_BOX=" "RUNNER=$W/runner-close-fails-once.sh" "DELETE_CONFIRM_TRIES=1"
drive "$d" R2 SR_HOST=slow SR_DELETE_FAIL=1; rc=$?
chk L4-gone-box-not-86 "80,0" "$rc,$(grep -c . "$d/live")" "$(grep -c 'BOX STILL UP' "$d/NOTES.md") still-up note, $(grep -c '^close' "$d/calls.log") closes"
echo "FAILS $FAILS"
exit "$FAILS"
