#!/usr/bin/env bash
# KIT10 runner tests, no box (stub-boat-kit10.sh, a stub PowerShell and a stub wslpath; every byte synthetic).
# KIT9's 30 checks (the progress floor on push-data, push and pull, the source push's timeout, the memcheck retry,
# exit codes), with B1-B3 and E1-E2 rewritten for main's ruling (no beacon, no grm-memread.exe: neither is read now),
# then KIT10's: per-call staging (S), the floor's start-up (G, REVIEW-KIT9 L3) and unreadable sizes (U, L6), extend
# (X, REVIEW-ASSURE L4), the vocabulary class (V, G11), H-1 run outputs (H), chain-lib.sh's exit (C, REVIEW-DEPS L1)
# and run returning a chain's exit (Q). One "RESULT id want got verdict detail" line per check, then "FAILS n";
# exits 1 on any failure. KIT10 round 2 (REVIEW-KIT10): W (push-data's SigDrug names and H-1 classes by destination,
# M1), S3b and S4 (the lane folder stays; two calls in a loop, L1), X2b and X6-X9 (extend's asked end and its END
# line, L2; X9 keeps the moved-back check live inside the 120 s), C5 (a chain with no step fails), P8 and G1 at 20 kB/s against the 10 kB/s floor (L5), and the V checks
# on a runner with no VOCAB block run that runner unchanged (L5: KIT9's control had rc 127 there).
# usage: test-runner-kit10.sh RUNNER WORK [CHAIN_LIB]
set -u
RUNNER=${1:?runner} W=${2:?work} CHAINLIB=${3:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/kit/chain-lib.sh}
T=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CHAINLIB=$(readlink -f "$CHAINLIB")   # the chains below source it from other folders
rm -rf "$W"; mkdir -p "$W/bin"
cat > "$W/bin/wslpath" << 'EOF'
#!/usr/bin/env bash
[ "$1" = -u ] && shift; p=${1//\\//}; echo "${WSLPATH_ROOT:?}${p#C:}"
EOF
chmod +x "$W/bin/wslpath"; export PATH=$W/bin:$PATH WSLPATH_ROOT=$W/winc
export GRM_BOAT=$T/stub-boat-kit10.sh GRM_BOAT_ADOPT=1 GRM_BOAT_STATE=$W/state GRM_BT_TIMEOUT=30
export GRM_STAGE_WIN='C:\Users\layth888\.t3\scratch\2026-10-28-test-3a358d9a\boat-grm\kit10-local'
STAGEBASE=$W/winc/Users/layth888/.t3/scratch/2026-10-28-test-3a358d9a/boat-grm/kit10-local
export GRM_LANES_MD=$W/lanes.md GRM_MEM_RETRY_SECS=0
echo "GRM_CAP=64" > "$W/lanes.md"
# the floor at test scale: 10 kB/s over 3 s windows (30 kB per window)
export GRM_SCP_WINDOW=3 GRM_SCP_FLOOR_MBPS=0.01
FAILS=0
res() { local v=as-predicted; [ "$2" = "$3" ] || { v=DIFFERS; FAILS=$((FAILS + 1)); }; printf 'RESULT %s want=%s got=%s %s %s\n' "$1" "$2" "$3" "$v" "${4:-}"; }
ID=bx_stub1
D=$W/root; mkdir -p "$D/fasta" "$D/catalog" "$D/curated"
python3 -I -c 'import random, sys; r = random.Random(9); open(sys.argv[1], "w").write(">c\n" + "".join(r.choice("ACGT") for _ in range(200000)) + "\n")' "$D/fasta/7770.1.fna"
echo synthetic > "$D/catalog/version.txt"; echo '[]' > "$D/curated/datasets.json"
(cd "$D" && sha256sum fasta/7770.1.fna catalog/version.txt curated/datasets.json) > "$W/m.sha256"; MS=$(sha256sum "$W/m.sha256" | cut -d' ' -f1)
G=$W/wt; mkdir -p "$G/src"; echo 'fn main() {}' > "$G/src/main.rs"; git -C "$G" init -q; git -C "$G" add .; git -C "$G" -c user.name=t -c user.email=t@t commit -qm t
python3 -I -c 'import os, sys; open(sys.argv[1], "wb").write(os.urandom(150000))' "$G/blob.bin"; git -C "$G" add .; git -C "$G" -c user.name=t -c user.email=t@t commit -qm blob
case_() { B=$W/box-$1 N=$W/notes-$1.md; rm -rf "$B"; mkdir -p "$B"; export STUB_ROOT=$B GRM_NOTES=$N; }
pd() { "$RUNNER" push-data "$ID" "$D" --manifest "$W/m.sha256" --manifest-sha256 "$MS" --dest "${1:-/home/user/grm-data}" > "$W/out${2:-}" 2>&1; }
tries() { grep -c 'try [0-9] of [0-9] failed' "$N" 2> /dev/null; }
scps() { cat "$B/scp.count" 2> /dev/null || echo 0; }
left() { find "$B/tmp" -name 'grm-*.tar' 2> /dev/null | wc -l; }

# ---------------------------------------------------------------- KIT9's checks
case_ p1; pd; rc=$?
res P1-push-data-fast "0,0,1" "$rc,$(tries),$(scps)" "rc,failed tries,scp calls; $(tail -n 1 "$W/out" | cut -c1-80)"
res P1-box-checked yes "$([ -s "$B/home/user/grm-data/.grm-input-manifest.sha256" ] && echo yes || echo no)" ""
case_ p2; t0=$(date +%s); STUB_SCP_RATE=2000 pd; rc=$?; el=$(( $(date +%s) - t0 ))
res P2-push-data-slow-all-tries "86,2,2,0" "$rc,$(tries),$(scps),$(left)" "rc,failed tries,scp calls,partial tars left; ${el}s"
res P2-floor-named-per-try 2 "$(grep -c 'failed | progress floor: the destination grew' "$N")" "$(grep -m1 -o 'progress floor[^|]*' "$N" | cut -c1-150)"
res P2-nothing-extracted no "$([ -e "$B/home/user/grm-data" ] && echo yes || echo no)" ""
res P2-stopped-by-floor-not-timeout yes "$([ "$el" -lt 40 ] && echo yes || echo no)" "${el}s (timeout per try is 1200 s)"
case_ p3; STUB_SCP_RATE=2000 STUB_SCP_SLOW_TRIES=1 pd; rc=$?
res P3-push-data-slow-first-try-only "0,1,2" "$rc,$(tries),$(scps)" "$(tail -n 1 "$W/out" | cut -c1-80)"
res P3-box-checked-after-retry yes "$(grep -q 'push-data root | ok' "$N" && echo yes || echo no)" ""
case_ p8; STUB_SCP_RATE=20000 pd; rc=$?
res P8-slow-but-above-floor-completes "0,0,1" "$rc,$(tries),$(scps)" "(20 kB/s against a 10 kB/s floor)"
case_ p4; t0=$(date +%s); STUB_SCP_HANG=1 "$RUNNER" push "$ID" "$G" > "$W/out" 2>&1; rc=$?; el=$(( $(date +%s) - t0 ))
res P4-source-push-silent-scp "86,2,2" "$rc,$(tries),$(scps)" "${el}s; $(grep -m1 -o 'progress floor[^|]*' "$N" | cut -c1-120)"
res P4-stopped-well-before-the-timeout yes "$([ "$el" -lt 60 ] && echo yes || echo no)" "${el}s (the try timeout is 600 s)"
case_ p5; t0=$(date +%s); STUB_SCP_HANG=1 GRM_SCP_FLOOR_MBPS=0 GRM_PUSH_TIMEOUT=4 "$RUNNER" push "$ID" "$G" > "$W/out" 2>&1; rc=$?; el=$(( $(date +%s) - t0 ))
res P5-source-push-timeout-floor-off "86,2,2" "$rc,$(tries),$(scps)" "${el}s; $(grep -m1 -o 'boat scp rc [0-9]*' "$N")"
res P5-timeout-rc-124 yes "$(grep -q 'boat scp rc 124' "$N" && echo yes || echo no)" ""
res P5-bounded yes "$([ "$el" -lt 45 ] && echo yes || echo no)" "${el}s for two 4 s tries"
case_ p9; "$RUNNER" push "$ID" "$G" > "$W/out" 2>&1; rc=$?
res P9-source-push-fast "0,0,1" "$rc,$(tries),$(scps)" "$(tail -n 1 "$W/out" | cut -c1-80)"
mkp() { mkdir -p "$B/home/user/stuff"; cp "$D/fasta/7770.1.fna" "$B/home/user/stuff/a.txt"; }
case_ p6; mkp; STUB_SCP_RATE=2000 "$RUNNER" pull "$ID" /home/user/stuff "$W/pulled6" > "$W/out" 2>&1; rc=$?
res P6-pull-slow "86,2,2" "$rc,$(tries),$(scps)" "$(grep -m1 -o 'progress floor[^|]*' "$N" | cut -c1-120)"
res P6-no-partial-local-tar 0 "$(find "$STAGEBASE" -name 'pull-*.tar' 2> /dev/null | wc -l)" ""
case_ p7; mkp; "$RUNNER" pull "$ID" /home/user/stuff "$W/pulled7" > "$W/out" 2>&1; rc=$?
res P7-pull-fast "0,0,yes" "$rc,$(tries),$(cmp -s "$D/fasta/7770.1.fna" "$W/pulled7/stuff/a.txt" && echo yes || echo no)" ""
case_ r1; bad=$W/bad.sha256; (cd "$D" && sha256sum fasta/7770.1.fna) | sed 's#  fasta/#  sealed/fasta/#' > "$bad"
"$RUNNER" push-data "$ID" "$D" --manifest "$bad" --manifest-sha256 "$(sha256sum "$bad" | cut -d' ' -f1)" > "$W/out" 2>&1; rc=$?
res R1-refusal-still-97-nothing-sent "97,0" "$rc,$(scps)" ""
"$RUNNER" > /dev/null 2>&1; res R2-usage-still-2 2 $? ""
PS=$W/pwsh.sh
cat > "$PS" << 'EOF'
#!/usr/bin/env bash
n=$(( $(cat "$PS_COUNT" 2> /dev/null || echo 0) + 1 )); echo $n > "$PS_COUNT"
case "$PS_MODE" in
  empty-then-ok) [ $n = 1 ] || echo "18.21 80.44";;
  empty) :;;
  garbled-then-ok) [ $n = 1 ] && echo "error: CIM" || echo "18.21 80.44";;
  low) echo "1.01 105.90";;
  ok) echo "18.21 80.44";;
esac
EOF
chmod +x "$PS"; export GRM_POWERSHELL=$PS
for m in "empty-then-ok 0 2" "garbled-then-ok 0 2" "empty 76 2" "low 76 1" "ok 0 1"; do
  set -- $m; case_ "m-$1"; export PS_MODE=$1 PS_COUNT=$W/ps-$1.count
  "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
  res "M-$1" "$2,$3" "$rc,$(cat "$PS_COUNT" 2> /dev/null || echo 0)" "rc,PowerShell reads; $(tail -n 1 "$W/out" | cut -c1-110)"
done
pscount() { cat "$PS_COUNT" 2> /dev/null || echo 0; }
# main's ruling (7 Nov): no beacon and no grm-memread.exe. Both are planted where KIT9 read them; neither is read.
export GRM_MEMDIR=$W/memdir; MD=$W/memdir; mkdir -p "$MD"
case_ b1; export PS_COUNT=$W/ps-b1.count PS_MODE=ok; echo "$(date +%s) 1.00 200.00" > "$MD/mem-beacon.txt"
GRM_MEM_BEACON=$MD/mem-beacon.txt "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res B1-fresh-low-beacon-not-read "0,1,0" "$rc,$(pscount),$(grep -c beacon "$W/out")" "$(tail -n 1 "$W/out" | cut -c1-110)"
case_ b2; export PS_COUNT=$W/ps-b2.count PS_MODE=empty; echo "$(date +%s) 18.21 80.44" > "$MD/mem-beacon.txt"
GRM_MEM_BEACON=$MD/mem-beacon.txt "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res B2-fresh-good-beacon-no-fallback-76 "76,2" "$rc,$(pscount)" "$(tail -n 1 "$W/out" | cut -c1-160)"
case_ b3; export PS_COUNT=$W/ps-b3.count PS_MODE=ok; echo "$(date +%s) n/a n/a" > "$MD/mem-beacon.txt"
GRM_MEM_BEACON=$MD/mem-beacon.txt "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res B3-garbled-beacon-irrelevant "0,1" "$rc,$(pscount)" "(REVIEW-KIT9 L5's case)"
rm -f "$MD/mem-beacon.txt"
printf '#!/usr/bin/env bash\ntouch "%s/exe-ran"; echo "17.50 81.00"\n' "$W" > "$MD/grm-memread.exe"; chmod +x "$MD/grm-memread.exe"
case_ e1; export PS_COUNT=$W/ps-e1.count PS_MODE=empty
GRM_MEMREAD=$MD/grm-memread.exe "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res E1-memread-exe-not-run "76,2,no" "$rc,$(pscount),$([ -e "$W/exe-ran" ] && echo yes || echo no)" ""
printf '#!/usr/bin/env bash\nexec sleep 1000\n' > "$MD/grm-memread.exe"
case_ e2; export PS_COUNT=$W/ps-e2.count PS_MODE=ok; t0=$(date +%s)
GRM_MEMREAD=$MD/grm-memread.exe GRM_MEM_TIMEOUT=2 "$RUNNER" preflight > "$W/out" 2>&1; rc=$?; el=$(( $(date +%s) - t0 ))
res E2-hung-exe-irrelevant "0,1,yes" "$rc,$(pscount),$([ $el -lt 10 ] && echo yes || echo no)" "${el}s (REVIEW-KIT9 L4's case)"
rm -f "$MD/grm-memread.exe"; export PS_MODE=ok PS_COUNT=$W/ps-rest.count

# ---------------------------------------------------------------- S: per-call staging (UI6's race)
# call a starts its scp 4 s late; call b, in the same lane, finishes first. With one shared stage-<lane>, b's exit
# removes a's tar before a reads it; here each call has its own folder.
case_ s1; t0=$(date +%s)
( STUB_SCP_DELAY=4 GRM_SCP_FLOOR_MBPS=0 pd /home/user/grm-a a ) & pa=$!
sleep 1; GRM_SCP_FLOOR_MBPS=0 pd /home/user/grm-b b; rb=$?; wait $pa; ra=$?
res S1-parallel-push-data-same-lane "0,0" "$ra,$rb" "$(tail -n 1 "$W/outa" | cut -c1-90)"
res S2-both-boxes-checked "yes,yes" "$([ -s "$B/home/user/grm-a/.grm-input-manifest.sha256" ] && echo yes || echo no),$([ -s "$B/home/user/grm-b/.grm-input-manifest.sha256" ] && echo yes || echo no)" ""
res S3-staging-cleaned 0 "$(find "$STAGEBASE" -mindepth 2 2> /dev/null | wc -l)" "(each call removed only its own folder)"
res S3b-lane-folder-kept 1 "$(find "$STAGEBASE" -mindepth 1 -maxdepth 1 -type d -name 'stage-*' 2> /dev/null | wc -l)" "(REVIEW-KIT10 L1: removing it when empty raced another call's mktemp)"
# S4 (REVIEW-KIT10 L1, as its scratch/stage-race/race.sh): the runner's own stage_dir and unstage, two calls in one
# lane at once, call b staging and unstaging while call a stages; a's mktemp must never find the lane folder gone.
# Round 1's file failed 84 of 2000 rounds on the controller (the review's race.sh), then 7 of 400 and 20 of 1500 in
# round 2's runs. S4_ROUNDS rounds (default 3000, about 10 s) leave a round-1 file a chance below 1e-9 of passing at 1%.
SR=$W/race; mkdir -p "$SR/bin"
printf '#!/usr/bin/env bash\n[ "$1" = -u ] && shift; echo "%s/winc"\n' "$SR" > "$SR/bin/wslpath"; chmod +x "$SR/bin/wslpath"
fn() { awk -v n="$1() {" 'index($0, n) == 1 { f = 1 } f { print } f && (/^}/ || (index($0, n) == 1 && /}$/)) { exit }' "$RUNNER"; }
{ fn stage_dir; fn unstage; } > "$SR/funcs.sh"   # the two functions only, each one line or up to its closing brace
cat > "$SR/call.sh" << 'EOS'
set -u
PATH=$1:$PATH STAGE_WIN='C:\x' LANE=k10 ID=bx_x ROUTE_RC=86
fail() { echo "FAIL rc=$1 $4" >&2; exit "$1"; }
. "$2"
stage_dir "$3"
[ "$3" = b ] && unstage
exit 0
EOS
if [ "$(grep -c 'mktemp -d' "$SR/funcs.sh")" = 1 ] && grep -q '^unstage() {' "$SR/funcs.sh"; then
  afail=0
  for i in $(seq "${S4_ROUNDS:-3000}"); do
    bash "$SR/call.sh" "$SR/bin" "$SR/funcs.sh" b 2> /dev/null & pb=$!
    bash "$SR/call.sh" "$SR/bin" "$SR/funcs.sh" a 2>> "$SR/a.err" & pa=$!
    wait $pb; wait $pa || afail=$((afail + 1))
    rm -rf "$SR/winc"
  done
  res S4-two-calls-race-on-the-lane-folder 0 "$afail" "(a's failures in ${S4_ROUNDS:-3000} rounds; first: $(head -n 1 "$SR/a.err" 2> /dev/null | cut -c1-90))"
else res S4-two-calls-race-on-the-lane-folder 0 no-stage_dir "(the runner has no per-call stage_dir and unstage)"; fi

# ---------------------------------------------------------------- G: start-up (REVIEW-KIT9 L3), U: unreadable sizes (L6)
case_ g1; STUB_SCP_DELAY=5 STUB_SCP_RATE=20000 GRM_SCP_START=12 pd; rc=$?
res G1-slow-start-then-healthy-completes "0,0,1" "$rc,$(tries),$(scps)" "(5 s silent, then 20 kB/s; KIT9's first window counted the silence)"
case_ g2; t0=$(date +%s); STUB_SCP_DELAY=30 GRM_SCP_START=6 pd; rc=$?; el=$(( $(date +%s) - t0 ))
res G2-no-first-byte-is-a-stall "86,2" "$rc,$(grep -c 'failed | progress floor: no byte reached the destination in 6 s' "$N")" "${el}s"
case_ u1; t0=$(date +%s); STUB_SCP_HANG=1 STUB_STAT_FAIL=1 GRM_PUSH_TIMEOUT=40 "$RUNNER" push "$ID" "$G" > "$W/out" 2>&1; rc=$?; el=$(( $(date +%s) - t0 ))
res U1-unreadable-size-stops "86,2" "$rc,$(grep -c "failed | progress floor: the destination's size could not be read 3 times in a row" "$N")" "${el}s"
res U2-before-the-try-timeout yes "$([ "$el" -lt 60 ] && echo yes || echo no)" "${el}s against two 40 s tries"

# ---------------------------------------------------------------- X: extend (REVIEW-ASSURE L4)
arch() { cat "$B/archive"; }
case_ x1; date -u -d '+3600 sec' +%FT%T.000Z > "$B/archive"; a0=$(arch)
"$RUNNER" extend "$ID" --ttl 1800 > "$W/out" 2>&1; rc=$?
res X1-shorter-refused-unchanged "2,yes,0" "$rc,$([ "$(arch)" = "$a0" ] && echo yes || echo no),$(grep -c '^extend' "$B/boat-calls.log")" "$(tail -n 1 "$W/out" | cut -c1-150)"
"$RUNNER" extend "$ID" --ttl 7200 > "$W/out" 2>&1; rc=$?; d=$(( $(date -u -d "$(arch)" +%s) - $(date +%s) ))
res X2-longer-set "0,yes,0" "$rc,$([ $d -ge 7150 ] && [ $d -le 7200 ] && echo yes || echo no),$(grep -c -- '--no-auto-stop' "$B/boat-calls.log")" "${d}s left; $(grep '^extend' "$B/boat-calls.log" | tail -n 1)"
res X3-noted 1 "$(grep -c 'extend --ttl 7200 | ok | archiveAfter' "$N")" ""
res X2b-end-line-is-the-new-end yes "$([ "$(sed -n 's/^END //p' "$W/out")" = "$(date -u -d "$(arch)" +%s)" ] && echo yes || echo no)" "$(grep '^END' "$W/out")"
case_ x4; echo null > "$B/archive"; rm -f "$B/archive"
"$RUNNER" extend "$ID" --ttl 7200 > "$W/out" 2>&1; rc=$?
res X4-no-archiveAfter-refused "86,0" "$rc,$(grep -c '^extend' "$B/boat-calls.log")" "$(tail -n 1 "$W/out" | cut -c1-120)"
case_ x5; date -u -d '+600 sec' +%FT%T.000Z > "$B/archive"
"$RUNNER" extend "$ID" --ttl abc > /dev/null 2>&1; r1=$?; STUB_ARCHIVE_BACK=1 "$RUNNER" extend "$ID" --ttl 3600 > "$W/out" 2>&1; r2=$?
res X5-bad-ttl-2-and-moved-back-86 "2,86" "$r1,$r2" "$(tail -n 1 "$W/out" | cut -c1-120)"
# REVIEW-KIT10 L2: boat extend returns 0 but leaves the end, or moves it less than asked
case_ x6; date -u -d '+600 sec' +%FT%T.000Z > "$B/archive"
STUB_ARCHIVE_SAME=1 "$RUNNER" extend "$ID" --ttl 3600 > "$W/out" 2>&1; rc=$?
res X6-extend-that-changed-nothing-86 "86,0" "$rc,$(grep -c '^END' "$W/out")" "$(tail -n 1 "$W/out" | cut -c1-140)"
case_ x7; date -u -d '+600 sec' +%FT%T.000Z > "$B/archive"
STUB_ARCHIVE_SHORT=600 "$RUNNER" extend "$ID" --ttl 3600 > "$W/out" 2>&1; rc=$?
res X7-extend-600s-short-86 "86,0" "$rc,$(grep -c '^END' "$W/out")" "$(tail -n 1 "$W/out" | cut -c1-140)"
case_ x8; date -u -d '+600 sec' +%FT%T.000Z > "$B/archive"
STUB_ARCHIVE_SHORT=60 "$RUNNER" extend "$ID" --ttl 3600 > "$W/out" 2>&1; rc=$?
res X8-extend-60s-short-within-120s-0 0 "$rc" "$(tail -n 1 "$W/out" | cut -c1-80)"
# an end moved back 60 s on an ask just past it: inside the 120 s, so only the a1 >= a0 check stops it
case_ x9; date -u -d '+3600 sec' +%FT%T.000Z > "$B/archive"
STUB_ARCHIVE_SHORT=70 "$RUNNER" extend "$ID" --ttl 3610 > "$W/out" 2>&1; rc=$?
res X9-moved-back-within-120s-86 "86,0" "$rc,$(grep -c '^END' "$W/out")" "$(tail -n 1 "$W/out" | cut -c1-140)"

# ---------------------------------------------------------------- V: the vocabulary (G11), H: H-1 run outputs
# A copy of the runner whose VOCAB pins are 13 synthetic files under the same names (the real pins are checked
# against grm-r9-curate/vocabulary on the laptop: REPORT.md); the unmodified runner refuses the synthetic bytes.
# The 13 names of grm-r9-curate/vocabulary (G11), pinned here so a runner without them is still tested on them.
V=$W/voc; mkdir -p "$V/vocabulary"
names=$(printf '%s\n' brand_generic_folds.tsv brand_name_spelling_folds.tsv closed_set_review.tsv display_format_rules.tsv \
  drug_name_typo_folds.tsv form_phrases.tsv ingredient_abbrev_folds.tsv preferred_base_names.tsv preferred_ingredients.tsv \
  release_phrases.tsv route_phrases.tsv salt_base_ingredients.tsv strength_grammar_review.tsv)
rnames=$(sed -n '/^VOCAB = {/,/^}/p' "$RUNNER" | grep -o '"[a-z_]*\.tsv"' | tr -d '"' | LC_ALL=C sort)
res V0-runner-pins-13-names yes "$([ "$rnames" = "$(LC_ALL=C sort <<< "$names")" ] && echo yes || echo no)" "($(grep -c . <<< "$rnames") names in the runner's VOCAB)"
for n in $names; do printf 'synthetic\t%s\n' "$n" > "$V/vocabulary/$n"; done
VR=$W/runner-synth-vocab.sh
# a copy of the runner with the 13 synthetic files' sha256s as its pins; a runner with no VOCAB block (KIT9) is
# copied as it is, so the V checks show how it treats a vocabulary push (REVIEW-KIT10 L5)
python3 -I - "$RUNNER" "$VR" "$V/vocabulary" << 'EOF'
import hashlib, os, re, sys
src, dst, v = sys.argv[1:]
t = open(src).read()
if "VOCAB = {\n" in t:
    a = t.index("VOCAB = {\n"); b = t.index("\n}\n", a)
    body = "".join('    "%s": "%s",\n' % (n, hashlib.sha256(open(os.path.join(v, n), "rb").read()).hexdigest()) for n in sorted(os.listdir(v)))
    t = t[:a] + "VOCAB = {\n" + body + t[b + 1:]
else:
    print("NOTE the runner has no VOCAB block; the V checks run it unchanged")
open(dst, "w").write(t); os.chmod(dst, 0o755)
EOF
vm() { (cd "$V" && find vocabulary -type f | LC_ALL=C sort | xargs sha256sum) > "$W/v.sha256"; VS=$(sha256sum "$W/v.sha256" | cut -d' ' -f1); }
vpd() { "${2:-$VR}" push-data "$ID" "$V" --manifest "$W/v.sha256" --manifest-sha256 "$VS" --dest "${1:-/home/user/grm-vocab}" > "$W/out" 2>&1; }
case_ v1; vm; vpd; rc=$?
res V1-whole-pinned-vocabulary-passes "0,13" "$rc,$(find "$B/home/user/grm-vocab/vocabulary" -type f 2> /dev/null | wc -l)" "$(tail -n 1 "$W/out" | cut -c1-80)"
case_ v2; vm; vpd "" "$RUNNER"; rc=$?
res V2-real-pins-refuse-other-bytes "97,0" "$rc,$(scps)" "$(grep -m1 'vocabulary/' "$W/out" | cut -c1-110)"
case_ v3; first=$(head -n 1 <<< "$names"); cp "$V/vocabulary/$first" "$W/keep"; echo x >> "$V/vocabulary/$first"; vm; vpd; rc=$?; cp "$W/keep" "$V/vocabulary/$first"
res V3-changed-file-refused "97,0" "$rc,$(scps)" ""
case_ v4; mv "$V/vocabulary/$first" "$W/keep"; vm; vpd; rc=$?; mv "$W/keep" "$V/vocabulary/$first"
res V4-partial-refused "97,1" "$rc,$(grep -q 'whole or not at all' "$W/out" && echo 1 || echo 0)" ""
case_ v5; echo extra > "$V/vocabulary/extra.tsv"; vm; vpd; rc=$?; rm -f "$V/vocabulary/extra.tsv"
res V5-unpinned-name-refused 97 "$rc" ""
case_ v6; mkdir -p "$V/vocabulary/sub"; cp "$V/vocabulary/$first" "$V/vocabulary/sub/$first"; vm; vpd; rc=$?; rm -rf "$V/vocabulary/sub"
res V6-subfolder-refused 97 "$rc" ""
case_ v7; mkdir -p "$V/x/vocabulary"; cp "$V/vocabulary/$first" "$V/x/vocabulary/"; (cd "$V" && sha256sum "x/vocabulary/$first") > "$W/v.sha256"; VS=$(sha256sum "$W/v.sha256" | cut -d' ' -f1); vpd; rc=$?; rm -rf "$V/x"
res V7-not-at-the-top-refused 97 "$rc" ""
P=$W/parts; mkdir -p "$P/parts/p1/h1work/arms/default/salmonella_enterica__x" "$P/parts/p1/h1work/j" "$P/parts/p1/h1work/report"
for f in spec.json inputs.json summary.json species-groups.json raw-manifest.txt j/contrasts.json report/ablation-report.md arms/default/salmonella_enterica__x/predictions.jsonl arms/default/salmonella_enterica__x/resume.json; do echo '{}' > "$P/parts/p1/h1work/$f"; done
echo "$(printf 'a%.0s' {1..64})  bvbrc/genomes.jsonl" > "$P/MASTER-MANIFEST-h1.sha256"
hm() { (cd "$P" && find . -type f | sed 's#^\./##' | LC_ALL=C sort | xargs sha256sum) > "$W/h.sha256"; HS=$(sha256sum "$W/h.sha256" | cut -d' ' -f1)
  "$RUNNER" push-data "$ID" "$P" --manifest "$W/h.sha256" --manifest-sha256 "$HS" --dest /home/user/grm-combine > "$W/out" 2>&1; }
case_ h1; hm; rc=$?
res H1-h1-run-outputs-and-master-copy-pass "0,10" "$rc,$(find "$B/home/user/grm-combine" -type f ! -name .grm-input-manifest.sha256 | wc -l)" "$(tail -n 1 "$W/out" | cut -c1-80)"
case_ h2; head -c 2000 /dev/urandom > "$P/parts/p1/h1work/arms/default/salmonella_enterica__x/model.bin"; hm; rc=$?; rm -f "$P/parts/p1/h1work/arms/default/salmonella_enterica__x/model.bin"
res H2-binary-member-refused 97 "$rc" "$(grep -m1 'model.bin' "$W/out" | cut -c1-100)"
case_ h3; head -c 3000 /dev/urandom > "$P/parts/p1/h1work/arms/default/salmonella_enterica__x/x.json"; hm; rc=$?; rm -f "$P/parts/p1/h1work/arms/default/salmonella_enterica__x/x.json"
res H3-not-text-refused 97 "$rc" "$(grep -m1 'UTF-8' "$W/out" | cut -c1-100)"
case_ h2b; echo '{}' > "$P/parts/p1/h1work/arms/default/salmonella_enterica__x/model.pkl"; hm; rc=$?; rm -f "$P/parts/p1/h1work/arms/default/salmonella_enterica__x/model.pkl"
res H2b-text-member-with-another-suffix-refused 97 "$rc" "$(grep -m1 'model.pkl' "$W/out" | cut -c1-100)"   # text, so only the name class can refuse it

# ---------------------------------------------------------------- W: SigDrug names and the H-1 classes by destination
# REVIEW-KIT10 M1: h1work/ and MASTER-MANIFEST*.sha256 passed for every destination, under any prefix, and a SigDrug
# vocabulary name passed under h1work/ or tables/ with any bytes. W1-W4 are the review's first four pd rows.
craft() { # NAME DEST PATH...: a root holding each PATH (synthetic text), pushed to DEST with its manifest
  local n=$1 dest=$2 p; shift 2; local CR=$W/craft-$n; rm -rf "$CR"; mkdir -p "$CR"
  for p in "$@"; do mkdir -p "$CR/$(dirname "$p")"; printf 'synthetic\t%s\n' "$p" > "$CR/$p"; done
  (cd "$CR" && printf '%s\n' "$@" | LC_ALL=C sort | xargs sha256sum) > "$W/c-$n.sha256"
  "$RUNNER" push-data "$ID" "$CR" --manifest "$W/c-$n.sha256" --manifest-sha256 "$(sha256sum "$W/c-$n.sha256" | cut -d' ' -f1)" --dest "$dest" > "$W/out" 2>&1; }
first_refusal() { grep -m1 -o 'line [0-9]*: .*' "$W/out" | cut -c1-120; }
case_ w1; craft w1 /home/user/grm-x h1work/x/brand_generic_folds.tsv; rc=$?
res W1-h1work-sigdrug-tsv-elsewhere-refused "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w2; craft w2 /home/user/grm-x p/q/r/h1work/a/b/c/d/e/notes.md; rc=$?
res W2-h1work-under-a-prefix-elsewhere-refused "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w3; craft w3 /home/user/grm-x p/q/r/MASTER-MANIFEST-anything.sha256; rc=$?
res W3-master-manifest-at-depth-3-refused "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w4; craft w4 /home/user/grm-x tables/brand_generic_folds.tsv; rc=$?
res W4-tables-sigdrug-name-refused "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w5; craft w5 /home/user/grm-x tables/Brand_Generic_Folds.tsv; rc=$?
res W5-sigdrug-name-in-other-case-refused "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w6; craft w6 /home/user/grm-combine parts/p1/h1work/brand_generic_folds.tsv; rc=$?
res W6-sigdrug-name-in-a-part-refused-on-combine "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w9; craft w9 /home/user/grm-x MASTER-MANIFEST-h1.sha256; rc=$?
res W9-master-manifest-at-the-top-elsewhere-refused "97,0" "$rc,$(scps)" "$(first_refusal)"
# W10, W11: on the combine destination too, h1work/ only under parts/<name>/ and the master manifest only at the top
case_ w10; craft w10 /home/user/grm-combine x/h1work/a.json; rc=$?
res W10-h1work-outside-parts-refused-on-combine "97,0" "$rc,$(scps)" "$(first_refusal)"
case_ w11; craft w11 /home/user/grm-combine parts/p1/MASTER-MANIFEST-h1.sha256; rc=$?
res W11-master-manifest-in-a-part-refused-on-combine "97,0" "$rc,$(scps)" "$(first_refusal)"
# W7, W8: the drive's own combine root, as h1-manifest.py combine builds it from two pulled WORKs (roots/ and logs/
# left out by it), to the combine destination and to another
H1MAN=$(dirname "$CHAINLIB")/h1-manifest.py
for pt in p1 s2; do
  pw=$W/pull-$pt/h1work; mkdir -p "$pw/arms/default/salmonella_enterica__x" "$pw/j" "$pw/k" "$pw/report" "$pw/roots/default/species" "$pw/logs"
  for f in spec.json inputs.json summary.json species-groups.json raw-manifest.txt j/contrasts.json k/manifest.json report/ablation-report.md \
    arms/default/salmonella_enterica__x/predictions.jsonl arms/default/salmonella_enterica__x/resume.json roots/default/species/x.json logs/run.log; do
    echo "{\"part\": \"$pt\"}" > "$pw/$f"; done
done
echo "$(printf 'a%.0s' {1..64})  bvbrc/genomes.jsonl" > "$W/MASTER-MANIFEST-h1.sha256"
rm -rf "$W/comb"; python3 -I "$H1MAN" combine "$W/comb" "$W/MASTER-MANIFEST-h1.sha256" "$(sha256sum "$W/MASTER-MANIFEST-h1.sha256" | cut -d' ' -f1)" \
  p1="$W/pull-p1/h1work" s2="$W/pull-s2/h1work" > "$W/comb.sha256"; crc=$?
cpush() { "$RUNNER" push-data "$ID" "$W/comb" --manifest "$W/comb.sha256" --manifest-sha256 "$(sha256sum "$W/comb.sha256" | cut -d' ' -f1)" --dest "$1" > "$W/out" 2>&1; }
case_ w7; cpush /home/user/grm-combine; rc=$?
res W7-drive-combine-root-to-grm-combine-passes "0,0,21,yes" "$crc,$rc,$(find "$B/home/user/grm-combine" -type f ! -name .grm-input-manifest.sha256 2> /dev/null | wc -l),$(grep -q '  MASTER-MANIFEST-h1.sha256$' "$W/comb.sha256" && echo yes || echo no)" "$(tail -n 1 "$W/out" | cut -c1-80)"
case_ w8; cpush /home/user/grm-data; rc=$?
res W8-same-root-to-another-destination-refused "97,0" "$rc,$(scps)" "$(first_refusal)"

# ---------------------------------------------------------------- C: chain-lib.sh (REVIEW-DEPS L1), Q: run returns it
mkchain() { # FILE: a chain with steps given as "name:rc" words after it
  local f=$1; shift
  { echo "CHAIN_LOG=$W/chainlog-$(basename "$f" .sh); . $CHAINLIB"
    for s in "$@"; do case "$s" in must:*) s=${s#must:}; echo "must ${s%%:*} bash -c 'exit ${s##*:}'";; *) echo "step ${s%%:*} bash -c 'exit ${s##*:}'";; esac; done
    echo finish; } > "$f"; }
mkchain "$W/c1.sh" a:0 b:3 c:0; bash "$W/c1.sh"; rc=$?
res C1-failed-middle-step-fails-the-chain "1,1" "$rc,$(grep -c 'CHAIN FAIL: b' "$W/chainlog-c1/STEPS")" "(the last step passing doesn't hide it)"
mkchain "$W/c2.sh" a:0 b:0; bash "$W/c2.sh"; rc=$?
res C2-all-pass-0 "0,1" "$rc,$(grep -c 'CHAIN PASS 2 steps' "$W/chainlog-c2/STEPS")" ""
mkchain "$W/c3.sh" must:a:4 b:0; bash "$W/c3.sh"; rc=$?
res C3-must-stops-the-chain "1,0" "$rc,$(grep -c ' b start' "$W/chainlog-c3/STEPS")" ""
printf 'CHAIN_LOG=%s/chainlog-c4\n. %s\nstep a bash -c "exit 2"\nstep b true\n' "$W" "$CHAINLIB" > "$W/c4.sh"; bash "$W/c4.sh"; rc=$?
res C4-control-no-finish-is-the-deps-bug 0 "$rc" "(a chain without finish exits with its last step's rc: why finish exists)"
printf 'CHAIN_LOG=%s/chainlog-c5\n. %s\nfinish\n' "$W" "$CHAINLIB" > "$W/c5.sh"; bash "$W/c5.sh"; rc=$?
res C5-no-step-fails "1,1" "$rc,$(grep -c 'CHAIN FAIL: no step ran' "$W/chainlog-c5/STEPS" 2> /dev/null)" "(REVIEW-KIT10 note: finish passed with 0 steps)"
case_ q1; mkdir -p "$B/home/user"; mkchain "$B/home/user/q1.sh" a:0 b:5 c:0
"$RUNNER" run "$ID" --log q1 --cwd /home/user -- bash /home/user/q1.sh > "$W/out" 2>&1; rc=$?
res Q1-run-returns-a-failed-chain 1 "$rc" "$(grep -o 'q1 ended rc=[0-9]*' "$W/out")"
mkchain "$B/home/user/q2.sh" a:0 b:0
"$RUNNER" run "$ID" --log q2 --cwd /home/user -- bash /home/user/q2.sh > "$W/out" 2>&1; rc=$?
res Q2-run-returns-a-passing-chain 0 "$rc" "$(grep -o 'q2 ended rc=[0-9]*' "$W/out")"
echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
