#!/usr/bin/env bash
# KIT10: drive.sh's KIND=h1 on a stub runner and stub boat (no box; every byte synthetic). One "RESULT id want got
# verdict detail" line per check, then "FAILS n"; exits 1 on any failure. KIT10 round 2 (REVIEW-KIT10): H4d (empty
# H1-FLAGS), H14b and H14c (the deadline is the runner's END, L2), H17 (a guard pause stops with its reason, L3), H18
# (a failed git status is not clean), H19 (a stage-table placeholder) and H20, H21 (the drafted run-h1.conf and
# stages-h1.tsv are refused before any box while they hold placeholders).
# usage: test-drive-h1.sh WORK DRIVE_SH KIT_DIR [DRAFTS_DIR]
set -u
W=${1:?work} DRV=${2:?drive.sh} KIT=${3:?kit dir}
T=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DRAFTS=$(readlink -f "${4:-$T/../drafts}")
rm -rf "$W"; mkdir -p "$W/bin"
cat > "$W/bin/boat" << 'EOF'
#!/usr/bin/env bash
echo "boat $*" >> "$SR_LOG.boat"
case "$1" in
  usage) echo '{"seconds": 2000}';;
  delete) i=$(printf '%s\n' "$@" | grep -m1 '^bx_'); sed -i "/^$i /d" "$SR_LIVE"; echo "deleted $i" >> "$SR_CLOSED";;
  info) i=$(printf '%s\n' "$@" | grep -m1 '^bx_'); grep -q "^$i " "$SR_LIVE" 2> /dev/null && { echo '{"sandbox": {"state": "ready"}}'; exit 0; }; echo "not found (404)"; exit 1;;
  list) echo '{"sandboxes": []}';;
  *) exit 1;;
esac
EOF
chmod +x "$W/bin/boat"
FAILS=0
res() { local v=as-predicted; [ "$2" = "$3" ] || { v=DIFFERS; FAILS=$((FAILS + 1)); }; printf 'RESULT %s want=%s got=%s %s %s\n' "$1" "$2" "$3" "$v" "${4:-}"; }
yn() { if "$@"; then echo yes; else echo no; fi; }
cnt() { [ -e "$1" ] && grep -c . "$1" || echo 0; }
# a synthetic data root, its master manifest, a vocabulary root and a source with the specs
DR=$W/root; mkdir -p "$DR/bvbrc/features" "$DR/fasta" "$DR/catalog" "$DR/species-check" "$DR/curated" "$DR/rev9"
echo '{}' > "$DR/bvbrc/genomes.jsonl"; echo '{}' > "$DR/bvbrc/genome_amr.jsonl"; echo '[]' > "$DR/bvbrc/features/1.1.json"
echo v > "$DR/catalog/version.txt"; echo '{}' > "$DR/species-check/species_check.json"; echo '[]' > "$DR/curated/datasets.json"
for i in 1.1 2.1 3.1; do echo ">$i" | gzip -n > "$DR/fasta/$i.fna.gz"; done
(cd "$DR" && find . -type f | sed 's#^\./##' | LC_ALL=C sort | xargs sha256sum) > "$W/MASTER-MANIFEST-test.sha256"
MSHA=$(sha256sum "$W/MASTER-MANIFEST-test.sha256" | cut -d' ' -f1)
VR=$W/vocab; mkdir -p "$VR/vocabulary"; printf 'a\tb\n' > "$VR/vocabulary/a.tsv"; printf 'c\td\n' > "$VR/vocabulary/b.tsv"
VSHA=$(python3 -I "$KIT/h1-manifest.py" vocab "$VR" | sha256sum | cut -d' ' -f1)
SRC=$W/src; mkdir -p "$SRC/src" "$SRC/scripts/input-comparison/h1"; echo 'pub const EVALUATION_REVISION: u32 = 12;' > "$SRC/src/analysis.rs"
echo '{"species": "salmonella_enterica"}' > "$SRC/scripts/input-comparison/h1/a.json"
git -C "$SRC" init -q; git -C "$SRC" add .; git -C "$SRC" -c user.name=t -c user.email=t@t commit -qm one; SC=$(git -C "$SRC" rev-parse HEAD)
printf '1.1\n2.1\n' > "$W/union"; printf '1.1\n9.9\n' > "$W/union-missing"
S=scripts/input-comparison/h1/a.json
{ printf 'P1\tpilot\txlarge\t0.67\t0.93\t1\t-\t%s\t1.5\t-\t1\n' "$S"
  printf 'S2\ts1\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\t-\t1\n' "$S"
  printf 'S3\ts1\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\t-\t1\n' "$S"
  printf 'SL\ts1\txlarge\t0.3\t0.9\t1\t-\t%s\t0.5\t-\t1\n' "$S"
  printf 'C1\tcombine\tdefault\t0.25\t0.4\t1\t$RUN_DIR/stages/P1/DONE,$RUN_DIR/stages/S2/DONE\t%s\t1.0\tP1,S2\t2\n' "$S"
  printf 'PH\ts1\txlarge\t0.3\t@@SLOW_ESTIMATE_H@@\t0\t-\t%s\t1.0\t-\t1\n' "$S"; } > "$W/stages.tsv"
conf() { # dir extra...
  local d=$1; shift; mkdir -p "$d"
  cat > "$d/run.conf" << EOF
RUNNER=$T/drv-stub-runner-kit10.sh
KIND=h1
STAGES=$W/stages.tsv
SRC=$SRC
SRC_COMMIT=$SC
REVISION=
KIT=$KIT KIT_REPO=$KIT
DATA_ROOT=$DR
MASTER_MANIFEST=$W/MASTER-MANIFEST-test.sha256
MASTER_SHA256=$MSHA
VOCAB_ROOT=$VR
VOCAB_MANIFEST_SHA256=$VSHA
REF_BOX=P1
H1_TRACE_BOXES=P1
LANE_PREFIX=k10t-
TEMPLATE=grm-tmpl-rust1981-c90ebaf
RAM_NEED_GIB=13.5
POLL=1 CKPT_SECS=1800
BUDGET_STOP=14 BUDGET_CAP=20
GRM_SCP_FLOOR_MBPS=5
RUN_DIR=$d
NOTES=$d/NOTES.md
COPY2=$d/copy2
EOF
  for x in "$@"; do printf '%s\n' "$x" >> "$d/run.conf"; done
}
drive() { local d=$1 b=$2; shift 2
  env GRM_BOAT_STATE="$d/state" SR_LOG="$d/calls.log" SR_CREATED="$d/created" SR_LIVE="$d/live" SR_CLOSED="$d/closed" \
    GRM_BOAT="$W/bin/boat" RUN_CONF="$d/run.conf" SR_UNION="$W/union" "$@" timeout 120 bash "$DRV" "$b" > "$d/out.$b" 2>&1; }
before=$( (cd "$W" && find root vocab -printf '%P %s %T@ %m\n' | LC_ALL=C sort) | sha256sum)

# P1 and S2 are both traced here, so C1 combines parts that ran with the same flags (H4b and H4c: they differ)
d=$W/a; conf "$d" "H1_TRACE_BOXES=P1,S2"; drive "$d" P1; rc=$?
res H1-pilot-box-rc 0 "$rc" "$(tail -n 1 "$d/out.P1")"
res H1-reference-binary-fixed yes "$(yn test -s "$d/BINARY-REF-$SC.sha")" ""
res H1-pushes-one-at-a-time-in-order "grm-meta,grm-vocab,grm-fasta" "$(grep -o 'push-data [^|]*--dest /home/user/grm-[a-z]*' "$d/calls.log" | grep -o 'grm-[a-z]*$' | paste -sd, -)" ""
res H1-meta-set-no-features-curated "4,0" "$(cnt "$d/calls.log.pd-grm-meta"),$(grep -c 'features/\|curated/\|rev9/\|fasta/' "$d/calls.log.pd-grm-meta")" ""
res H1-fasta-set-is-the-union-verbatim yes "$(yn bash -c "grep -E '  fasta/(1\.1|2\.1)\.fna\.gz$' '$W/MASTER-MANIFEST-test.sha256' | cmp -s - '$d/calls.log.pd-grm-fasta'")" ""
res H1-vocab-set-pinned yes "$(yn test "$(sha256sum < "$d/calls.log.pd-grm-vocab" | cut -d' ' -f1)" = "$VSHA")" ""
res H1-pins-in-notes "1,1,1" "$(grep -c 'vocabulary manifest pinned' "$d/NOTES.md"),$(grep -c 'FASTA manifest pinned' "$d/NOTES.md"),$(grep -c 'stage input manifest pinned' "$d/NOTES.md")" ""
res H1-stage-env "1,1,1,1" "$(grep -c "SPEC=$S" "$d/calls.log.env"),$(grep -c 'H1_ROLE=run' "$d/calls.log.env"),$(grep -c 'H1_TRACE=1' "$d/calls.log.env"),$(grep -c 'H1_FASTA_SHA256=' "$d/calls.log.env")" ""
res H1-floor-5-reaches-the-runner 0 "$(grep 'push-data' "$d/calls.log" | grep -vc 'floor=5$')" "$(grep -m1 'push-data' "$d/calls.log" | grep -o 'floor=.*')"
res H1-ttl-from-column-and-template yes "$(yn grep -q -- '--ttl 5400 --from-template grm-tmpl-rust1981-c90ebaf' "$d/created")" "$(cut -d' ' -f3- "$d/created")"
res H1-done-closed "yes,1" "$(yn test -e "$d/stages/P1/DONE"),$(cnt "$d/closed")" ""
res H1-measured-disk-reported 1 "$(grep -c 'H-1 measured disk' "$d/NOTES.md")" ""
drive "$d" S2 SR_SHA="$(printf 'b%.0s' {1..64})"; rc=$?
res H2-other-binary-stops-82 "82,yes" "$rc,$(yn test -e "$d/BINARY-MISMATCH-$SC")" ""
rm -f "$d/BINARY-MISMATCH-$SC"; drive "$d" S2; rc=$?
res H3-equal-binary-runs 0 "$rc" "$(grep -c 'equals source' "$d/NOTES.md") equal-gate lines"
nb=$(grep -c 'r11-box.sh host' "$d/calls.log")
drive "$d" C1; rc=$?
res H4-combine-box-rc 0 "$rc" "$(tail -n 1 "$d/out.C1")"
res H4-combine-default-no-template yes "$(yn bash -c "grep 'c1-' '$d/created' | grep -q -- '--type default --ttl 3600$'")" "$(grep 'c1-' "$d/created" | cut -d' ' -f3-)"
res H4-combine-parts-and-master "yes,yes,yes" "$(yn test -s "$d/stages/C1/combine-in/parts/p1/h1work/spec.json"),$(yn test -s "$d/stages/C1/combine-in/parts/s2/h1work/arms/default/ds1/runs.json"),$(yn cmp -s "$d/stages/C1/combine-in/MASTER-MANIFEST-test.sha256" "$W/MASTER-MANIFEST-test.sha256")" ""
res H4-combine-pushed-once-to-grm-combine "1,0" "$(grep -c 'push-data.*--dest /home/user/grm-combine' "$d/calls.log"),$(( $(grep -c 'r11-box.sh host' "$d/calls.log") - nb ))" "(no build on a combine box)"
res H4-combine-env "1,1" "$(grep -c 'H1_ROLE=combine' "$d/calls.log.env"),$(grep -c 'H1_PARTS=p1,s2' "$d/calls.log.env")" ""
# DC1 (the H-1 dry run): a traced and an untraced part refuse before the combine box, and so do parts with no H1-FLAGS
d=$W/m; conf "$d"; drive "$d" P1; drive "$d" S2; drive "$d" C1; rc=$?
res H4b-combine-parts-with-other-flags-2-no-box "2,0" "$rc,$(grep -c 'c1-' "$d/created")" "$(grep -o 'the parts ran with.*' "$d/out.C1" | cut -c1-110)"
rm -f "$d/stages/P1/H1-FLAGS" "$d/stages/S2/H1-FLAGS"; drive "$d" C1; rc=$?
res H4c-combine-parts-without-flags-2-no-box "2,0" "$rc,$(grep -c 'c1-' "$d/created")" ""
# KIT10 round 2: empty H1-FLAGS files read as one shared empty flag set, so round 1 combined them
: > "$d/stages/P1/H1-FLAGS"; : > "$d/stages/S2/H1-FLAGS"; drive "$d" C1; rc=$?
res H4d-combine-parts-with-empty-flags-2-no-box "2,0" "$rc,$(grep -c 'c1-' "$d/created")" "$(grep -o 'the parts ran with.*' "$d/out.C1" | cut -c1-120)"

d=$W/b; conf "$d"; drive "$d" S2; rc=$?
res H5-run-box-before-the-reference-79 "79,0" "$rc,$(cnt "$d/created")" ""
drive "$d" C1; rc=$?
res H6-combine-before-its-parts-79 "79,0" "$rc,$(cnt "$d/created")" ""
d=$W/c; conf "$d" "VOCAB_MANIFEST_SHA256=$(printf 'c%.0s' {1..64})"; drive "$d" P1; rc=$?
res H7-vocabulary-pin-mismatch-2-no-box "2,0" "$rc,$(cnt "$d/created")" "$(grep -o 'the vocabulary in.*' "$d/out.P1" | cut -c1-90)"
d=$W/e; conf "$d"; drive "$d" P1 SR_UNION="$W/union-missing"; rc=$?
res H8-union-genome-without-fasta-2-box-closed-reason-named "2,1,1,yes" "$rc,$(cnt "$d/created"),$(cnt "$d/closed"),$(yn grep -q "STOP rc=2: FASTA manifest: 1 of the union's 2 genomes have no FASTA line" "$d/out.P1")" "$(grep -o 'FASTA manifest:.*' "$d/out.P1" | head -n 1 | cut -c1-100)"
d=$W/f; conf "$d" "SRC_COMMIT=$(printf '1%.0s' {1..40})"; drive "$d" P1; rc=$?
res H9-other-source-commit-2-no-box "2,0" "$rc,$(cnt "$d/created")" ""
d=$W/g; conf "$d" "COPY2=@@COPY2@@"; drive "$d" P1; rc=$?
res H10-placeholder-2 2 "$rc" ""
d=$W/h; conf "$d"; drive "$d" P1 SR_STATUS=running FORCE_TTL_NEAR_AFTER=1; rc=$?
res H11-ttl-near-is-83-no-resume "83,yes" "$rc,$(yn test -e "$d/stages/P1/TTL-EXHAUSTED")" ""
d=$W/i; conf "$d"; drive "$d" P1 SR_HOST=slow; rc=$?
res H12-slow-host-slow-ok-in-ttl-goes-on "0,0" "$rc,$(grep -c '^extend' "$d/calls.log")" "(P1: TTL 1.5 h = slow + 0.5 h, so no extend)"
drive "$d" S2 SR_HOST=slow; rc=$?
res H13-slow-host-not-slow-ok-80 80 "$rc" ""
# H14 has its own run folder: after H12 and H13 the ledger of $W/i ends in two slow hosts, and drive.sh's budget rule
# (two slow hosts in a row stop the run) refuses SL with 78 before it makes a box. Here P1 (fast) fixes the reference.
d=$W/k; conf "$d"; drive "$d" P1; drive "$d" SL SR_HOST=slow; rc=$?
res H14-extend-through-the-runner "0,1,0" "$rc,$(grep -c '^extend bx_[a-z0-9]* --ttl [0-9]*' "$d/calls.log"),$(grep -c 'extend' "$d/calls.log.boat")" "$(grep -m1 '^extend' "$d/calls.log" | cut -c1-60)"
# REVIEW-KIT10 L2: the box got an end about 140 s away (the stub's END is 4900 s short of the asked); the deadline is
# that end, so the first running poll is TTL-near (83, its checkpoint pulled). Round 1 kept T0 + want and ran on.
d=$W/k2; conf "$d"; drive "$d" P1; drive "$d" SL SR_HOST=slow SR_EXTEND_SHORT=4900 SR_RUNNING_POLLS=3; rc=$?
res H14b-deadline-is-the-end-the-box-got "83,yes,1" "$rc,$(yn test -e "$d/stages/SL/TTL-EXHAUSTED"),$(grep -c 'TTL-near' "$d/NOTES.md")" "$(grep -o 'the box ends [^ ]*' "$d/NOTES.md" | head -n 1)"
d=$W/k3; conf "$d"; drive "$d" P1; drive "$d" SL SR_HOST=slow SR_EXTEND_NOEND=1; rc=$?
res H14c-extend-without-end-86-closed "86,2" "$rc,$(cnt "$d/closed")" "$(grep -o 'extend gave no END.*' "$d/NOTES.md" | head -n 1 | cut -c1-90)"
d=$W/j; conf "$d"; drive "$d" P1 SR_PUSHDATA_FAIL=fasta; rc=$?
res H15-refused-push-data-stops-closes "97,1" "$rc,$(cnt "$d/closed")" ""
# REVIEW-KIT10 L3: r11-box.sh's status says paused when the H-1 stage's guard holds it stopped
d=$W/p; conf "$d"; drive "$d" P1 SR_STATUS=paused; rc=$?
res H17-guard-pause-stops-81-with-its-reason "81,1,1,no" "$rc,$(grep -c 'STOP rc=81: P1 was paused by its guard (guard disk: ' "$d/NOTES.md"),$(grep -c 'log stop' "$d/calls.log"),$(yn test -e "$d/stages/P1/TTL-EXHAUSTED")" "$(grep -o 'STOP rc=81.*' "$d/NOTES.md" | cut -c1-120)"
# REVIEW-KIT10 note: a git status that fails is not a clean worktree
G=$W/gitfail; mkdir -p "$G"; printf '#!/usr/bin/env bash\nfor a in "$@"; do [ "$a" = --untracked-files=all ] && exit 128; done\nexec %s "$@"\n' "$(command -v git)" > "$G/git"; chmod +x "$G/git"
d=$W/q; conf "$d"; drive "$d" P1 PATH="$G:$PATH"; rc=$?
res H18-git-status-failing-2-no-box "2,0" "$rc,$(cnt "$d/created")" "$(grep -o 'git status failed.*' "$d/out.P1" | cut -c1-90)"
# a placeholder in the stage table's row (here the slow estimate, which round 1 read as 0 and went on)
d=$W/r; conf "$d"; drive "$d" P1; drive "$d" PH; rc=$?
res H19-stage-table-placeholder-2-no-box "2,1" "$rc,$(cnt "$d/created")" "$(grep -o 'placeholder left in.*' "$d/out.PH" | cut -c1-110)"
# H20, H21: the drafted run file and stage table (next/run-h1.conf, next/stages-h1.tsv). As drafted, the run file's
# placeholders refuse first; with every one of them given a test value, the stage table's do. No box either way.
if [ -f "$DRAFTS/run-h1.conf" ] && [ -f "$DRAFTS/stages-h1.tsv" ]; then
  d=$W/s; mkdir -p "$d"
  env GRM_BOAT_STATE="$d/state" SR_LOG="$d/calls.log" SR_CREATED="$d/created" SR_LIVE="$d/live" SR_CLOSED="$d/closed" \
    GRM_BOAT="$W/bin/boat" RUN_CONF="$DRAFTS/run-h1.conf" timeout 60 bash "$DRV" P1 > "$d/out.P1" 2>&1; rc=$?
  ph=$(grep -o 'placeholder left in: .*' "$d/out.P1")
  res H20-drafted-run-file-refused-2-no-box "2,yes,0,0" "$rc,$(yn test -n "$ph"),$(cnt "$d/created"),$(cnt "$d/calls.log")" "${ph:0:200}"
  d=$W/t; mkdir -p "$d"; cp "$DRAFTS/run-h1.conf" "$d/run.conf"
  # the test's own values for the stub run, then one for each placeholder the draft names
  { echo "RUNNER=$T/drv-stub-runner-kit10.sh"; echo "STAGES=$DRAFTS/stages-h1.tsv"; echo "SRC=$SRC"; echo "SRC_COMMIT=$SC"
    echo "KIT=$KIT KIT_REPO=$KIT DATA_ROOT=$DR MASTER_MANIFEST=$W/MASTER-MANIFEST-test.sha256 MASTER_SHA256=$MSHA"
    echo "VOCAB_ROOT=$VR VOCAB_MANIFEST_SHA256=$VSHA RUN_DIR=$d NOTES=$d/NOTES.md COPY2=$d/copy2"; } >> "$d/run.conf"
  unknown=""
  for v in $(sed -n 's/^drive: placeholder left in: //p' "$W/s/out.P1"); do
    case "$v" in
      RUNNER|STAGES|SRC|SRC_COMMIT|KIT|KIT_REPO|DATA_ROOT|MASTER_MANIFEST|MASTER_SHA256|VOCAB_ROOT|VOCAB_MANIFEST_SHA256|RUN_DIR|NOTES|COPY2) ;;
      H1_TRACE_BOXES) echo "H1_TRACE_BOXES=P1" >> "$d/run.conf";;
      BUDGET_STOP) echo "BUDGET_STOP=14" >> "$d/run.conf";;
      TEMPLATE) echo "TEMPLATE=grm-tmpl-test" >> "$d/run.conf";;
      RAM_NEED_GIB) echo "RAM_NEED_GIB=13.5" >> "$d/run.conf";;
      *) unknown="$unknown $v";;
    esac
  done
  env GRM_BOAT_STATE="$d/state" SR_LOG="$d/calls.log" SR_CREATED="$d/created" SR_LIVE="$d/live" SR_CLOSED="$d/closed" \
    GRM_BOAT="$W/bin/boat" RUN_CONF="$d/run.conf" SR_UNION="$W/union" timeout 60 bash "$DRV" P1 > "$d/out.P1" 2>&1; rc=$?
  res H21-drafted-stage-table-refused-2-no-box "2,yes,0,none" "$rc,$(yn grep -q "placeholder left in $DRAFTS/stages-h1.tsv: @@" "$d/out.P1"),$(cnt "$d/created"),${unknown:-none}" "$(grep -o 'STOP rc=2: .*' "$d/out.P1" | cut -c1-160)"
else
  res H20-drafted-run-file-refused-2-no-box present missing "(no $DRAFTS/run-h1.conf or stages-h1.tsv)"
fi
after=$( (cd "$W" && find root vocab -printf '%P %s %T@ %m\n' | LC_ALL=C sort) | sha256sum)
res H16-data-and-vocabulary-roots-untouched yes "$(yn test "$before" = "$after")" ""
echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
