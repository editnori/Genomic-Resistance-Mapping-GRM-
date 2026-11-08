#!/usr/bin/env bash
# KIT12: drive.sh's k add-on box (kind addon, KIND=h1; G2(b)) on KIT10's stub runner and a stub boat (no box; every
# byte synthetic). An add-on waits for the reference binary and for its main part's DONE and pulled WORK; refuses,
# before any box, a column 10 that isn't one run box of the table and H-1 flags that aren't the main part's; is a run
# box in every step, from the data snapshot or by push-data; sends the main part's records, pinned here, through
# push-data to /home/user/grm-combine after the vocabulary; and names the part and pin in stage.env. A combine box
# takes the add-on as a part. One "RESULT id want got verdict detail" line per check, then "FAILS n"; exits 1 on any
# failure. The stub boat and synthetic inputs are test-drive-h1-kit11.sh's, copied.
# usage: test-drive-h1-kit12.sh WORK DRIVE_SH KIT_DIR
set -u
W=${1:?work} DRV=${2:?drive.sh} KIT=${3:?kit dir}
T=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
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
cnt() { local n=0; [ -s "$1" ] && n=$(grep -c . "$1"); echo "$n"; }   # one line, even for an empty file
# a synthetic data root (read only by the push-data run that the snapshot run is compared with), its master
# manifest, a vocabulary root and its manifest, and a source with the spec
DR=$W/root; mkdir -p "$DR/bvbrc/features" "$DR/fasta" "$DR/catalog" "$DR/species-check" "$DR/curated"
echo '{}' > "$DR/bvbrc/genomes.jsonl"; echo '{}' > "$DR/bvbrc/genome_amr.jsonl"; echo '[]' > "$DR/bvbrc/features/1.1.json"
echo v > "$DR/catalog/version.txt"; echo '{}' > "$DR/species-check/species_check.json"; echo '[]' > "$DR/curated/datasets.json"
for i in 1.1 2.1 3.1; do echo ">$i" | gzip -n > "$DR/fasta/$i.fna.gz"; done
(cd "$DR" && find . -type f | sed 's#^\./##' | LC_ALL=C sort | xargs sha256sum) > "$W/MASTER-MANIFEST-test.sha256"
MSHA=$(sha256sum "$W/MASTER-MANIFEST-test.sha256" | cut -d' ' -f1)
VR=$W/vocab; mkdir -p "$VR/vocabulary"; printf 'a\tb\n' > "$VR/vocabulary/a.tsv"; printf 'c\td\n' > "$VR/vocabulary/b.tsv"
python3 -I "$KIT/h1-manifest.py" vocab "$VR" > "$W/VOCAB-MANIFEST.sha256"; VSHA=$(sha256sum "$W/VOCAB-MANIFEST.sha256" | cut -d' ' -f1)
SRC=$W/src; mkdir -p "$SRC/src" "$SRC/scripts/input-comparison/h1"; echo 'pub const EVALUATION_REVISION: u32 = 12;' > "$SRC/src/analysis.rs"
echo '{"species": "salmonella_enterica"}' > "$SRC/scripts/input-comparison/h1/a.json"
git -C "$SRC" init -q; git -C "$SRC" add .; git -C "$SRC" -c user.name=t -c user.email=t@t commit -qm one; SC=$(git -C "$SRC" rev-parse HEAD)
printf '1.1\n2.1\n' > "$W/union"
S=scripts/input-comparison/h1/a.json
{ printf 'P1\tpilot\txlarge\t0.67\t0.93\t1\t-\t%s\t1.5\t-\t1\n' "$S"
  printf 'S2\ts1\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\t-\t1\n' "$S"
  printf 'K2\taddon\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\tS2\t1\n' "$S"
  printf 'KD\taddon\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\t-\t1\n' "$S"
  printf 'KT\taddon\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\tS2,P1\t1\n' "$S"
  printf 'KN\taddon\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\tNOPE\t1\n' "$S"
  printf 'KA\taddon\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\tK2\t1\n' "$S"
  printf 'KC\taddon\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\tC1\t1\n' "$S"
  printf 'C1\tcombine\tdefault\t0.25\t0.4\t1\t-\t%s\t1.0\tS2,K2\t2\n' "$S"; } > "$W/stages.tsv"
SNAP=grm-tmpl-h1data-7dd00779
conf() { # dir mode(snap|push) extra...
  local d=$1 m=$2; shift 2; mkdir -p "$d"
  cat > "$d/run.conf" << EOF2
RUNNER=$T/drv-stub-runner-kit10.sh
KIND=h1
STAGES=$W/stages.tsv
SRC=$SRC
SRC_COMMIT=$SC
REVISION=
KIT=$KIT KIT_REPO=$KIT
MASTER_MANIFEST=$W/MASTER-MANIFEST-test.sha256
MASTER_SHA256=$MSHA
VOCAB_MANIFEST_SHA256=$VSHA
REF_BOX=P1
H1_TRACE_BOXES=P1,S2,K2,KA
LANE_PREFIX=k12t-
RAM_NEED_GIB=13.5
POLL=1 CKPT_SECS=1800
BUDGET_STOP=14 BUDGET_CAP=20
GRM_SCP_FLOOR_MBPS=5
RUN_DIR=$d
NOTES=$d/NOTES.md
EOF2
  if [ "$m" = snap ]; then
    printf '%s\n' "DATA_SNAPSHOT_DIR=/home/user/grm-h1raw" "VOCAB_SNAPSHOT_DIR=/home/user/grm-h1vocab" \
      "VOCAB_MANIFEST=$W/VOCAB-MANIFEST.sha256" "TEMPLATE=$SNAP" >> "$d/run.conf"
  else
    printf '%s\n' "DATA_ROOT=$DR" "VOCAB_ROOT=$VR" "TEMPLATE=grm-tmpl-rust1981-c90ebaf" >> "$d/run.conf"
  fi
  for x in "$@"; do printf '%s\n' "$x" >> "$d/run.conf"; done
}
drive() { local d=$1 b=$2; shift 2
  env GRM_BOAT_STATE="$d/state" SR_LOG="$d/calls.log" SR_CREATED="$d/created" SR_LIVE="$d/live" SR_CLOSED="$d/closed" \
    GRM_BOAT="$W/bin/boat" RUN_CONF="$d/run.conf" SR_UNION="$W/union" "$@" timeout 120 bash "$DRV" "$b" > "$d/out.$b" 2>&1; }
boxes() { cnt "$1/created"; }
mark() { cnt "$1/calls.log"; }                     # calls so far, for since()
since() { tail -n +"$(( $2 + 1 ))" "$1/calls.log"; } # the calls after mark
envof() { grep -F "STAGE=$2" "$1/calls.log.env" | tail -n 1; }   # the stage-env line drive.sh wrote for box $2

# the snapshot run: P1 fixes the binary; K2 before S2 waits; then S2 and K2
d=$W/a; conf "$d" snap; drive "$d" P1; r0=$?; drive "$d" K2; rc=$?
res A1-addon-waits-for-its-main-part "0,79,1,1" "$r0,$rc,$(boxes "$d"),$(grep -c 'main part S2 has no finished, pulled WORK' "$d/out.K2")" "$(tail -n 1 "$d/out.K2" | cut -c1-90)"
drive "$d" S2; r1=$?; m=$(mark "$d"); drive "$d" K2; rc=$?; K=$(since "$d" "$m")
res A2-addon-runs "0,0,yes,1" "$r1,$rc,$(yn test -e "$d/stages/K2/DONE"),$(grep -c '^close ' <<< "$K")" "$(tail -n 1 "$d/out.K2" | cut -c1-90)"
res A2-a-run-box-s-data-then-the-main-part "grm-meta,grm-vocab,grm-combine,grm-fasta" "$(grep -oE '^(stage-local|push-data) .*--dest /home/user/grm-[a-z]*' <<< "$K" | grep -o 'grm-[a-z]*$' | paste -sd, -)" ""
res A2-records-by-push-data-data-by-stage-local "1,3" "$(grep -c '^push-data .*--dest /home/user/grm-combine' <<< "$K"),$(grep -c '^stage-local ' <<< "$K")" ""
RM=$d/stages/K2/REF-MANIFEST.sha256 RS=$(sha256sum "$d/stages/K2/REF-MANIFEST.sha256" 2> /dev/null | cut -d' ' -f1)
cp "$d/calls.log.pd-grm-combine" "$W/k2-sent.sha256" 2> /dev/null
res A2-the-main-part-s-records-are-sent "yes,yes,0" "$(yn cmp -s "$RM" "$W/k2-sent.sha256"),$(yn grep -q -- "--manifest-sha256 $RS " <<< "$(grep '^push-data .*grm-combine' <<< "$K")"),$(cut -c67- "$RM" | grep -vc -e '^parts/s2/h1work/' -e '^MASTER-MANIFEST-test\.sha256$')" "$(wc -l < "$RM" 2> /dev/null) lines"
res A2-records-hold-the-main-part-s-work "yes,yes" "$(yn grep -q '  parts/s2/h1work/summary.json$' "$RM"),$(yn grep -q '  parts/s2/h1work/arms/default/ds1/runs.json$' "$RM")" ""
res A2-pinned-in-notes-before-the-box 1 "$(grep -c "main part S2's records pinned: $RM sha256 $RS" "$d/NOTES.md")" ""
E=$(envof "$d" K2)
res A2-stage-env-names-the-part-and-pin "H1_ROLE=addon,H1_REFERENCE=s2,H1_REFERENCE_SHA256=$RS,1" "$(grep -o 'H1_ROLE=[a-z]*' <<< "$E"),$(grep -o 'H1_REFERENCE=[a-z0-9]*' <<< "$E"),$(grep -o 'H1_REFERENCE_SHA256=[0-9a-f]*' <<< "$E"),$(grep -c 'H1_META_SHA256=[0-9a-f]\{64\}' <<< "$E")" ""
E2=$(envof "$d" S2)
res A3-a-run-box-names-no-reference "H1_ROLE=run,H1_REFERENCE=,H1_REFERENCE_SHA256=" "$(grep -o 'H1_ROLE=[a-z]*' <<< "$E2"),$(grep -o 'H1_REFERENCE=[a-z0-9]*' <<< "$E2"),$(grep -o 'H1_REFERENCE_SHA256=[0-9a-f]*' <<< "$E2")" ""
res A4-addon-box-from-the-data-snapshot-with-the-binary-gate "yes,1" "$(yn bash -c "grep k12t-k2 '$d/created' | grep -q -- '--from-template $SNAP'"),$(grep -c 'K2 | binary gate: grm [0-9a-f]* equals source' "$d/NOTES.md")" ""
m=$(mark "$d"); drive "$d" C1; rc=$?
res A5-combine-takes-the-addon-as-a-part "0,yes,yes" "$rc,$(yn grep -q '  parts/k2/h1work/' "$d/calls.log.pd-grm-combine"),$(yn grep -q '  parts/s2/h1work/' "$d/calls.log.pd-grm-combine")" "$(tail -n 1 "$d/out.C1" | cut -c1-90)"
# a column 10 that isn't one run box of the table: none, two, an unknown box, an add-on, a combine box. KA traces as
# its main part K2 does (H1_TRACE_BOXES above), so only the kind check stops it, not the flags check after it
rcs=$(for b in KD KT KN KA KC; do drive "$d" "$b"; echo $?; done | paste -sd, -)
res A6-bad-main-part-cells-2-no-box "2,2,2,2,2,0" "$rcs,$(grep -cE ' grm-k12t-k[dtnac]-[0-9a-f]{6} ' "$d/created")" ""
res A6-each-names-its-cell "1,1,1,1,1" "$(grep -c 'names its one main part in column 10, not -$' "$d/out.KD"),$(grep -c 'names its one main part in column 10, not S2,P1$' "$d/out.KT"),$(grep -c "its main part NOPE is not a box of $W/stages.tsv\$" "$d/out.KN"),$(grep -c 'its main part K2 is of kind addon, not a run box$' "$d/out.KA"),$(grep -c 'its main part C1 is of kind combine, not a run box$' "$d/out.KC")" ""
# the add-on's H-1 flags must be the main part's (run.py refuses a reference that ran with other flags)
d=$W/b; conf "$d" snap "H1_TRACE_BOXES=P1,S2"; drive "$d" P1; drive "$d" S2; drive "$d" K2; rc=$?
res A7-other-flags-than-the-main-part-2-no-box "2,2,1" "$rc,$(boxes "$d"),$(grep -c "this add-on's H-1 flags (H1_TRACE=0 H1_EXTRA=) are not main part S2's (H1_TRACE=1 H1_EXTRA=)" "$d/out.K2")" "$(tail -n 1 "$d/out.K2" | cut -c1-90)"
d=$W/c; conf "$d" snap; drive "$d" P1; drive "$d" S2; rm -f "$d/stages/S2/H1-FLAGS"; drive "$d" K2; rc=$?
res A8-main-part-without-flags-2-no-box "2,2,1" "$rc,$(boxes "$d"),$(grep -c 'main part S2 records no H-1 flags' "$d/out.K2")" ""
d=$W/e; conf "$d" snap; drive "$d" K2; rc=$?
res A9-addon-waits-for-the-reference-binary "79,0,1" "$rc,$(boxes "$d"),$(grep -c 'the H-1 reference box P1 has not fixed' "$d/out.K2")" ""
d=$W/f; conf "$d" snap; drive "$d" P1; drive "$d" S2; drive "$d" K2 SR_PUSHDATA_FAIL=combine; rc=$?
res A10-records-refused-97-closed "97,3,3,yes" "$rc,$(boxes "$d"),$(cnt "$d/closed"),$(yn grep -q "push-data of main part S2's records refused or failed" "$d/NOTES.md")" ""
d=$W/g; conf "$d" push; drive "$d" P1; drive "$d" S2; m=$(mark "$d"); drive "$d" K2; rc=$?; K=$(since "$d" "$m")
res A11-push-data-mode-addon "0,grm-meta,grm-vocab,grm-combine,grm-fasta,0" "$rc,$(grep -o '^push-data .*--dest /home/user/grm-[a-z]*' <<< "$K" | grep -o 'grm-[a-z]*$' | paste -sd, -),$(grep -c '^stage-local ' <<< "$K")" ""
d=$W/h; conf "$d" snap; drive "$d" P1; drive "$d" S2; rm -rf "$d/stages/S2/out/h1work"; drive "$d" K2; rc=$?
res A12-main-part-done-but-not-pulled-79 "79,2" "$rc,$(boxes "$d")" "$(tail -n 1 "$d/out.K2" | cut -c1-90)"
d=$W/i; conf "$d" snap; drive "$d" P1; drive "$d" S2; rm -f "$d/stages/S2/DONE"; drive "$d" K2; rc=$?
res A13-main-part-pulled-but-not-done-79 "79,2" "$rc,$(boxes "$d")" "$(tail -n 1 "$d/out.K2" | cut -c1-90)"
# a rerun of the add-on (its DONE removed) builds the main part's set afresh: the same manifest, sent again
d=$W/j; conf "$d" snap; drive "$d" P1; drive "$d" S2; drive "$d" K2; r1=$?; cp "$d/stages/K2/REF-MANIFEST.sha256" "$W/j-ref1"
rm -f "$d/stages/K2/DONE"; m=$(mark "$d"); drive "$d" K2; rc=$?; K=$(since "$d" "$m")
res A14-rerun-builds-the-set-afresh "0,0,yes,1" "$r1,$rc,$(yn cmp -s "$W/j-ref1" "$d/stages/K2/REF-MANIFEST.sha256"),$(grep -c '^push-data .*--dest /home/user/grm-combine' <<< "$K")" "$(tail -n 1 "$d/out.K2" | cut -c1-90)"

echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
