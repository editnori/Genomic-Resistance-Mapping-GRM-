#!/usr/bin/env bash
# KIT11 runner tests, no box (stub-boat-kit10.sh; every byte synthetic): the Linux host mode (L) and stage-local (T).
# L: GRM_HOST=linux needs GRM_BOAT, reads this host's memory and cgroup room against GRM_MIN_AVAIL_MB, and stages
# under GRM_STAGE_DIR with no wslpath call (a wslpath that fails is first on PATH); GRM_HOST=wsl stages as before.
# T: stage-local checks the manifest here (pin, lines, classes, deny, vocabulary; receipts and grm-combine refused)
# before any box call, sends only the manifest, re-hashes it in the box, links each file without following a link,
# and runs push-data's box checks; every box failure leaves no destination, no manifest copy and the source intact.
# One "RESULT id want got verdict detail" line per check, then "FAILS n"; exits 1 on any failure.
# usage: test-runner-kit11.sh RUNNER WORK
set -u
RUNNER=${1:?runner} W=${2:?work}
T=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
rm -rf "$W"; mkdir -p "$W/bin" "$W/nowsl"
# wsl mode's stub wslpath (as in test-runner-kit10.sh); for linux mode a wslpath that logs any call and fails
cat > "$W/bin/wslpath" << 'EOF'
#!/usr/bin/env bash
[ "$1" = -u ] && shift; p=${1//\\//}; echo "${WSLPATH_ROOT:?}${p#C:}"
EOF
cat > "$W/nowsl/wslpath" << 'EOF'
#!/usr/bin/env bash
echo "wslpath $*" >> "${WSL_CALLS:?}"; exit 1
EOF
chmod +x "$W/bin/wslpath" "$W/nowsl/wslpath"
PATH0=$PATH
export PATH=$W/bin:$PATH0 WSLPATH_ROOT=$W/winc WSL_CALLS=$W/wslpath-calls.log
export GRM_BOAT=$T/stub-boat-kit10.sh GRM_BOAT_ADOPT=1 GRM_BOAT_STATE=$W/state GRM_BT_TIMEOUT=30
export GRM_STAGE_WIN='C:\Users\layth888\.t3\scratch\2026-10-28-test-3a358d9a\boat-grm\kit11-local'
WSTAGE=$W/winc/Users/layth888/.t3/scratch/2026-10-28-test-3a358d9a/boat-grm/kit11-local
export GRM_LANES_MD=$W/lanes.md GRM_MEM_RETRY_SECS=0
echo "GRM_CAP=64" > "$W/lanes.md"
export GRM_SCP_WINDOW=3 GRM_SCP_FLOOR_MBPS=0.01
FAILS=0
res() { local v=as-predicted; [ "$2" = "$3" ] || { v=DIFFERS; FAILS=$((FAILS + 1)); }; printf 'RESULT %s want=%s got=%s %s %s\n' "$1" "$2" "$3" "$v" "${4:-}"; }
yn() { if "$@"; then echo yes; else echo no; fi; }
ID=bx_stub1
LSTAGE=$W/lstage
lin() { env PATH="$W/nowsl:$PATH0" GRM_HOST=linux GRM_STAGE_DIR="$LSTAGE" "$@"; }
case_() { B=$W/box-$1 N=$W/notes-$1.md; rm -rf "$B"; mkdir -p "$B/tmp" "$B/home/user"; export STUB_ROOT=$B GRM_NOTES=$N; }
scps() { cat "$B/scp.count" 2> /dev/null || echo 0; }
wcalls() { local n=0; [ -s "$WSL_CALLS" ] && n=$(grep -c . "$WSL_CALLS"); echo "$n"; }   # one line, even for an empty log
left() { find "$B/tmp" -name 'grm-stage-*' 2> /dev/null | wc -l; }
first() { grep -m1 -E "$1" "$W/out" | cut -c1-120; }
# the data: a public-class set (FASTA, catalog, curated) and its manifest, as the snapshot holds it on the box
D=$W/root; mkdir -p "$D/fasta" "$D/catalog" "$D/curated"
python3 -I -c 'import random, sys; r = random.Random(11); open(sys.argv[1], "w").write(">c\n" + "".join(r.choice("ACGT") for _ in range(50000)) + "\n")' "$D/fasta/7770.1.fna"
echo synthetic > "$D/catalog/version.txt"; echo '[]' > "$D/curated/datasets.json"
(cd "$D" && sha256sum catalog/version.txt curated/datasets.json fasta/7770.1.fna) > "$W/m.sha256"; MS=$(sha256sum "$W/m.sha256" | cut -d' ' -f1)
SRCB=/home/user/grm-h1raw
mkbox() { mkdir -p "$B$SRCB"; cp -a "$D/." "$B$SRCB/"; echo unlisted > "$B$SRCB/catalog/extra.txt"; }
before() { (cd "$B$SRCB" && find . -printf '%P %y %s %i\n' | LC_ALL=C sort) | sha256sum | cut -c1-16; }
pinof() { sha256sum "$1" | cut -d' ' -f1; }
sl() { # [runner] -- stage-local with FROM, MAN, PIN, DEST (defaults: the snapshot set to /home/user/grm-meta)
  lin "${SLR:-$RUNNER}" stage-local "$ID" --from "${FROM:-$SRCB}" --manifest "${MAN:-$W/m.sha256}" --manifest-sha256 "${PIN:-$MS}" \
    --dest "${DEST:-/home/user/grm-meta}" > "$W/out" 2>&1; }

# ---------------------------------------------------------------- L: the Linux host mode
case_ l1; env -u GRM_BOAT GRM_HOST=linux "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res L1-linux-needs-grm-boat "2,1" "$rc,$(grep -c 'needs GRM_BOAT' "$W/out")" "$(head -n 1 "$W/out" | cut -c1-90)"
case_ l2; GRM_HOST=mac "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res L2-unknown-host-2 2 "$rc" "$(head -n 1 "$W/out" | cut -c1-90)"
case_ l3; : > "$WSL_CALLS"; lin env GRM_MIN_AVAIL_MB=0 "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res L3-linux-preflight-reads-this-host "0,1,0" "$rc,$(grep -c 'MemAvailable [0-9]* MiB, cgroup room -\{0,1\}[0-9]* MiB' "$W/out"),$(wcalls)" "$(grep -o 'host:.*' "$W/out" | cut -c1-140)"
case_ l4; lin env GRM_MIN_AVAIL_MB=999999999 "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res L4-linux-memory-under-the-floor-76 "76,1" "$rc,$(grep -c 'each needs 999999999 or more' "$W/out")" "$(grep -o 'memcheck.*' "$W/out" | cut -c1-140)"
case_ l5; lin env GRM_MIN_AVAIL_MB=1e3 "$RUNNER" preflight > "$W/out" 2>&1; rc=$?
res L5-bad-floor-2 2 "$rc" ""
case_ l6; : > "$WSL_CALLS"; rm -rf "$LSTAGE"
lin "$RUNNER" push-data "$ID" "$D" --manifest "$W/m.sha256" --manifest-sha256 "$MS" --dest /home/user/grm-data > "$W/out" 2>&1; rc=$?
res L6-linux-push-data-stages-here "0,yes,0,0" "$rc,$(yn test -d "$LSTAGE/stage-kit10-stub"),$(find "$LSTAGE" -mindepth 2 | wc -l),$(wcalls)" "rc, lane folder, call folders left, wslpath calls; $(tail -n 1 "$W/out" | cut -c1-70)"
res L6-box-checked yes "$(yn test -s "$B/home/user/grm-data/.grm-input-manifest.sha256")" ""
case_ l7; rm -rf "$WSTAGE"
"$RUNNER" push-data "$ID" "$D" --manifest "$W/m.sha256" --manifest-sha256 "$MS" --dest /home/user/grm-data > "$W/out" 2>&1; rc=$?
res L7-wsl-mode-stages-as-before "0,yes" "$rc,$(yn test -d "$WSTAGE/stage-kit10-stub")" ""
case_ l8; mkdir -p "$B/home/user/pullme"; echo pulled > "$B/home/user/pullme/f.txt"; rm -rf "$W/pulled"; : > "$WSL_CALLS"
lin "$RUNNER" pull "$ID" /home/user/pullme "$W/pulled" > "$W/out" 2>&1; rc=$?
res L8-linux-pull "0,pulled,0" "$rc,$(cat "$W/pulled/pullme/f.txt" 2> /dev/null),$(wcalls)" "$(tail -n 1 "$W/out" | cut -c1-80)"

# ---------------------------------------------------------------- T: stage-local
case_ t1; mkbox; : > "$WSL_CALLS"; sl; rc=$?
res T1-stage-local-passes "0,1,0,0" "$rc,$(scps),$(left),$(wcalls)" "rc, scp calls (the manifest only), manifest copies left, wslpath calls; $(tail -n 1 "$W/out" | cut -c1-80)"
res T1-exact-set-and-manifest "4,yes" "$(find "$B/home/user/grm-meta" -type f | wc -l),$(yn cmp -s "$B/home/user/grm-meta/.grm-input-manifest.sha256" "$W/m.sha256")" "(3 files + the manifest; the unlisted extra.txt stays out)"
res T1-hard-links-not-copies yes "$(yn test "$(stat -c %i "$B/home/user/grm-meta/fasta/7770.1.fna")" = "$(stat -c %i "$B$SRCB/fasta/7770.1.fna")")" ""
res T1-sent-no-tar "0,$(stat -c %s "$W/m.sha256")" "$(grep -c '^scp .*\.tar' "$B/boat-calls.log"),$(stat -c %s "$B/home/user/grm-meta/.grm-input-manifest.sha256")" "(tar scp calls, bytes of the manifest that arrived)"
res T1-noted "1,0" "$(grep -c "stage-local $SRCB | ok | 3 files linked" "$N"),$(find "$LSTAGE/stage-kit10-stub" -mindepth 1 | wc -l)" "(notes line, call folders left)"
case_ t2; mkbox; PIN=$(printf 'f%.0s' {1..64}) sl; rc=$?
res T2-wrong-pin-97-no-box-call "97,0,no" "$rc,$(scps),$(yn test -e "$B/home/user/grm-meta")" "$(first 'pinned|refused' )"
case_ t3; mkbox; mkdir -p "$B$SRCB/curated"; echo '[]' > "$B$SRCB/curated/heldout.json"
{ cat "$W/m.sha256"; (cd "$B$SRCB" && sha256sum curated/heldout.json); } > "$W/t3.sha256"
MAN=$W/t3.sha256 PIN=$(pinof "$W/t3.sha256") sl; rc=$?
res T3-denied-segment-97-no-box-call "97,0" "$rc,$(scps)" "$(first 'refused segment')"
case_ t4; mkbox; mkdir -p "$B$SRCB/audit"; echo '{}' > "$B$SRCB/audit/run.json"
{ cat "$W/m.sha256"; (cd "$B$SRCB" && sha256sum audit/run.json); } > "$W/t4.sha256"
MAN=$W/t4.sha256 PIN=$(pinof "$W/t4.sha256") sl; rc=$?
res T4-receipt-goes-through-push-data-97 "97,0,yes" "$rc,$(scps),$(yn grep -q 'go through push-data' "$W/out")" ""
case_ t5; mkbox; DEST=/home/user/grm-combine sl; rc=$?
res T5-grm-combine-refused-97 "97,0,yes" "$rc,$(scps),$(yn grep -q 'does not fill /home/user/grm-combine' "$W/out")" ""
case_ t6; FROM=/home/user/grm-heldout-raw sl; r1=$?; FROM=/home/user/sigdrug-raw sl; r2=$?; FROM=/home/user/KSD-raw sl; r3=$?
res T6-denied-source-folder-97 "97,97,97,0" "$r1,$r2,$r3,$(scps)" ""
case_ t7; mkbox
FROM=/tmp/raw sl; a=$?; DEST=/tmp/meta sl; b=$?; DEST=$SRCB/sub sl; c=$?; DEST=$SRCB sl; d=$?; FROM=/home/user/a/../grm-h1raw sl; e=$?
DEST=/home/user/grm-meta/./x sl; f=$?; FROM=/home/user/grm-meta/raw DEST=/home/user/grm-meta sl; g=$?
res T7-folders-must-be-plain-and-apart "2,2,2,2,2,2,2,0" "$a,$b,$c,$d,$e,$f,$g,$(scps)" "(/tmp source, /tmp dest, dest inside source, dest = source, .., ., source inside dest; scp calls)"
case_ t8; mkbox; b0=$(before); sl > /dev/null; rc0=$?; rm -rf "$B/home/user/grm-meta"; b1=$(before)
res T8-a-pass-leaves-the-source-as-it-was "0,$b0" "$rc0,$b1" "(listing hash of the source: names, types, sizes, inodes)"
case_ t9; mkbox; mkdir -p "$B/elsewhere"; mv "$B$SRCB/catalog/version.txt" "$B/elsewhere/version.txt"; ln -s "$B/elsewhere/version.txt" "$B$SRCB/catalog/version.txt"
b0=$(before); sl; rc=$?
res T9-file-link-in-source-86 "86,no,0,$b0" "$rc,$(yn test -e "$B/home/user/grm-meta"),$(left),$(before)" "$(first 'LINKERR')"
res T9-names-the-link 1 "$(grep -c 'LINKERR not a regular file in the source folder: catalog/version.txt' "$W/out")" ""
case_ t10; mkbox; mkdir -p "$B/elsewhere"; mv "$B$SRCB/catalog" "$B/elsewhere/catalog"; ln -s "$B/elsewhere/catalog" "$B$SRCB/catalog"; sl; rc=$?
res T10-folder-link-in-source-86 "86,no,1" "$rc,$(yn test -e "$B/home/user/grm-meta"),$(grep -c 'LINKERR catalog/version.txt' "$W/out")" "$(first 'LINKERR')"
case_ t11; mkbox; ln -s "$B$SRCB" "$B/home/user/grm-h1link"; FROM=/home/user/grm-h1link sl; rc=$?
res T11-source-folder-a-link-86 "86,no,1" "$rc,$(yn test -e "$B/home/user/grm-meta"),$(grep -c 'source folder is missing or a link' "$W/out")" ""
case_ t12; mkbox; printf 'G' >> "$B$SRCB/fasta/7770.1.fna"; sl; rc=$?
res T12-changed-bytes-86-source-kept "86,no,yes,0" "$rc,$(yn test -e "$B/home/user/grm-meta"),$(yn test -s "$B$SRCB/fasta/7770.1.fna"),$(left)" "$(first 'REFUSED')"
case_ t13; mkbox; rm "$B$SRCB/curated/datasets.json"; sl; rc=$?
res T13-missing-source-file-86 "86,no,1" "$rc,$(yn test -e "$B/home/user/grm-meta"),$(grep -c 'LINKERR curated/datasets.json' "$W/out")" ""
case_ t14; mkbox; mkdir -p "$B/home/user/grm-meta"; echo keep > "$B/home/user/grm-meta/keep.txt"; sl; rc=$?
res T14-dest-not-empty-86-kept "86,keep,0" "$rc,$(cat "$B/home/user/grm-meta/keep.txt" 2> /dev/null),$(left)" ""
case_ t15; mkbox; mkdir -p "$B/elsewhere/empty"; ln -s "$B/elsewhere/empty" "$B/home/user/grm-meta"; sl; rc=$?
res T15-dest-a-link-86 "86,0,1" "$rc,$(find "$B/elsewhere/empty" -mindepth 1 | wc -l),$(grep -c 'destination is a link' "$W/out")" ""
case_ t16; mkbox; STUB_SCP_TAMPER=1 sl; rc=$?
res T16-manifest-changed-on-the-way-86 "86,no,0,1" "$rc,$(yn test -e "$B/home/user/grm-meta"),$(left),$(grep -c 'manifest copy is not the pinned one' "$W/out")" ""
case_ t17; mkbox; { cat "$W/m.sha256"; head -n 1 "$W/m.sha256"; } > "$W/t17.sha256"; MAN=$W/t17.sha256 PIN=$(pinof "$W/t17.sha256") sl; rc=$?
res T17-line-twice-97 "97,0" "$rc,$(scps)" "$(first 'twice')"
case_ t18; mkbox; sl > /dev/null; DEST=/home/user/grm-meta2 sl; rc=$?
res T18-two-sets-from-one-source "0,yes" "$rc,$(yn test "$(stat -c %i "$B/home/user/grm-meta2/fasta/7770.1.fna")" = "$(stat -c %i "$B/home/user/grm-meta/fasta/7770.1.fna")")" ""
# the vocabulary class, on a runner copy whose 13 pins are synthetic files' (as test-runner-kit10.sh's V checks)
V=$W/voc; mkdir -p "$V/vocabulary"
names=$(sed -n '/^VOCAB = {/,/^}/p' "$RUNNER" | grep -o '"[a-z_]*\.tsv"' | tr -d '"' | LC_ALL=C sort)
for n in $names; do printf 'synthetic\t%s\n' "$n" > "$V/vocabulary/$n"; done
VR=$W/runner-synth-vocab.sh
python3 -I - "$RUNNER" "$VR" "$V/vocabulary" << 'EOF'
import hashlib, os, sys
src, dst, v = sys.argv[1:]
t = open(src).read()
a = t.index("VOCAB = {\n"); b = t.index("\n}\n", a)
body = "".join('    "%s": "%s",\n' % (n, hashlib.sha256(open(os.path.join(v, n), "rb").read()).hexdigest()) for n in sorted(os.listdir(v)))
open(dst, "w").write(t[:a] + "VOCAB = {\n" + body + t[b + 1:]); os.chmod(dst, 0o755)
EOF
(cd "$V" && find vocabulary -type f | LC_ALL=C sort | xargs sha256sum) > "$W/v.sha256"
vbox() { mkdir -p "$B/home/user/grm-h1vocab"; cp -a "$V/vocabulary" "$B/home/user/grm-h1vocab/"; }
case_ t19; vbox; SLR=$VR FROM=/home/user/grm-h1vocab MAN=$W/v.sha256 PIN=$(pinof "$W/v.sha256") DEST=/home/user/grm-vocab sl; rc=$?
res T19-whole-pinned-vocabulary-links "0,13" "$rc,$(find "$B/home/user/grm-vocab/vocabulary" -type f 2> /dev/null | wc -l)" "$(tail -n 1 "$W/out" | cut -c1-80)"
case_ t20; vbox; sed 1d "$W/v.sha256" > "$W/v12.sha256"; SLR=$VR FROM=/home/user/grm-h1vocab MAN=$W/v12.sha256 PIN=$(pinof "$W/v12.sha256") DEST=/home/user/grm-vocab sl; rc=$?
res T20-partial-vocabulary-97 "97,0,yes" "$rc,$(scps),$(yn grep -q 'whole or not at all' "$W/out")" ""
case_ t21; vbox; FROM=/home/user/grm-h1vocab MAN=$W/v.sha256 PIN=$(pinof "$W/v.sha256") DEST=/home/user/grm-vocab sl; rc=$?
res T21-real-pins-refuse-other-bytes "97,0" "$rc,$(scps)" "$(first 'vocabulary/')"
case_ t22; mkbox; lin "$RUNNER" stage-local "$ID" "$SRCB" --manifest "$W/m.sha256" --manifest-sha256 "$MS" > "$W/out" 2>&1; a=$?
lin "$RUNNER" stage-local "$ID" --from "$SRCB" --manifest-sha256 "$MS" > "$W/out" 2>&1; b=$?
lin "$RUNNER" stage-local "$ID" --from "$SRCB" --manifest "$W/nofile" --manifest-sha256 "$MS" > "$W/out" 2>&1; c=$?
res T22-usage-2 "2,2,2,0" "$a,$b,$c,$(scps)" "(a bare folder, no --manifest, a missing manifest)"

echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
