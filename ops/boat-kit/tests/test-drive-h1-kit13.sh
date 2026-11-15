#!/usr/bin/env bash
# KIT13: drive.sh's two-slow-hosts stop (step 1) on KIT10's stub runner and a stub boat (no box; every byte
# synthetic). The stop stays the default: with the budget ledger's last two rows on slow hosts, a drive with
# SLOW_PAIR_PASS unset or 0 exits 78 before any box, with the text ops' restart scripts read. SLOW_PAIR_PASS=1, from
# the environment (a restart script) or the run file, lets a SLOW_OK stage pass with a note and run to DONE; a stage
# that isn't SLOW_OK still stops, the 35 cap still refuses, and the ledger still records the new box's host. Any
# value but 0 or 1 is refused with 2 before any box. One "RESULT id want got verdict detail" line per check, then
# "FAILS n"; exits 1 on any failure. The stub boat and synthetic inputs are test-drive-h1-kit12.sh's, cut down.
# usage: test-drive-h1-kit13.sh WORK DRIVE_SH KIT_DIR
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
# a synthetic data root, its master manifest, a vocabulary root and a source with the spec
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
# P1 is the reference box and SLOW_OK; S2 is a main part that isn't
{ printf 'P1\tpilot\txlarge\t0.67\t0.93\t1\t-\t%s\t1.5\t-\t1\n' "$S"
  printf 'S2\ts1\txlarge\t0.3\t0.5\t0\t-\t%s\t1.0\t-\t1\n' "$S"; } > "$W/stages.tsv"
conf() { # dir extra...: a push-data run file
  local d=$1; shift; mkdir -p "$d"
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
LANE_PREFIX=k13t-
RAM_NEED_GIB=13.5
POLL=1 CKPT_SECS=1800
BUDGET_STOP=14 BUDGET_CAP=20
GRM_SCP_FLOOR_MBPS=5
RUN_DIR=$d
NOTES=$d/NOTES.md
DATA_ROOT=$DR
VOCAB_ROOT=$VR
TEMPLATE=grm-tmpl-rust1981-c90ebaf
EOF2
  for x in "$@"; do printf '%s\n' "$x" >> "$d/run.conf"; done
}
ledger() { # dir xlarge_hours host...: a budget ledger whose rows end with these hosts, in order
  local d=$1 xl=$2 h; shift 2
  printf 'utc\tbox\tid\tcpu\thost\txlarge_hours\tvm_seconds\tresult\n' > "$d/budget.tsv"
  for h in "$@"; do printf '2026-11-11T07:00:00Z\tX\tbx_old\tStub CPU\t%s\t%s\t1\tdone\n' "$h" "$xl" >> "$d/budget.tsv"; done
}
# SLOW_PAIR_PASS comes from the environment unless the case sets it; env -u clears any the caller had
drive() { local d=$1 b=$2; shift 2
  env -u SLOW_PAIR_PASS GRM_BOAT_STATE="$d/state" SR_LOG="$d/calls.log" SR_CREATED="$d/created" SR_LIVE="$d/live" \
    SR_CLOSED="$d/closed" GRM_BOAT="$W/bin/boat" RUN_CONF="$d/run.conf" SR_UNION="$W/union" "$@" \
    timeout 120 bash "$DRV" "$b" > "$d/out.$b" 2>&1; }
boxes() { cnt "$1/created"; }
STOP='the last two boxes had slow hosts; stop and report'   # the text finish.sh (LEDGER C16) reads
NOTE='the last two boxes had slow hosts; SLOW_PAIR_PASS=1 and stage P1 is SLOW_OK, so this drive goes on'

# the default: the stop, on a SLOW_OK stage too, before any box; the exact text a restart script reads
d=$W/a; conf "$d"; ledger "$d" 0.1 fast slow slow; drive "$d" P1; rc=$?
res S1-default-stops-78-no-box "78,0,1" "$rc,$(boxes "$d"),$(grep -c "STOP rc=78: budget: $STOP\$" "$d/NOTES.md")" "$(tail -n 1 "$d/out.P1" | cut -c1-90)"
d=$W/b; conf "$d"; ledger "$d" 0.1 slow slow; drive "$d" P1 SLOW_PAIR_PASS=0; rc=$?
res S2-pass-0-stops-78-no-box "78,0,1" "$rc,$(boxes "$d"),$(grep -c "$STOP" "$d/out.P1")" ""
d=$W/c; conf "$d" SLOW_PAIR_PASS=0; ledger "$d" 0.1 slow slow; drive "$d" P1 SLOW_PAIR_PASS=1; rc=$?
res S3-run-file-0-overrides-the-environment "78,0" "$rc,$(boxes "$d")" ""
# the pass, from a restart script's environment: a SLOW_OK stage runs to DONE, with the note, and its row is booked
d=$W/e; conf "$d"; ledger "$d" 0.1 slow slow; drive "$d" P1 SLOW_PAIR_PASS=1; rc=$?
res S4-pass-1-slow-ok-stage-runs "0,1,yes,1" "$rc,$(boxes "$d"),$(yn test -e "$d/stages/P1/DONE"),$(cnt "$d/closed")" "$(tail -n 1 "$d/out.P1" | cut -c1-90)"
res S4-the-note-is-in-notes-not-the-stop "1,0" "$(grep -c "P1 | budget: $NOTE\$" "$d/NOTES.md"),$(grep -c "$STOP" "$d/NOTES.md")" ""
res S4-the-ledger-books-the-new-box "4,P1,fast" "$(cnt "$d/budget.tsv"),$(tail -n 1 "$d/budget.tsv" | cut -f2),$(tail -n 1 "$d/budget.tsv" | cut -f5)" ""
# the pass from the run file
d=$W/f; conf "$d" SLOW_PAIR_PASS=1; ledger "$d" 0.1 slow slow; drive "$d" P1; rc=$?
res S5-pass-1-from-the-run-file "0,1,1" "$rc,$(boxes "$d"),$(grep -c "$NOTE" "$d/NOTES.md")" "$(tail -n 1 "$d/out.P1" | cut -c1-90)"
# a stage that isn't SLOW_OK still stops (P1 fixes the reference binary first, on a ledger with no slow pair)
d=$W/g; conf "$d"; ledger "$d" 0.1 fast; drive "$d" P1; r0=$?
for _ in 1 2; do printf '2026-11-11T07:00:00Z\tX\tbx_old\tStub CPU\tslow\t0.1\t1\tdone\n' >> "$d/budget.tsv"; done
drive "$d" S2 SLOW_PAIR_PASS=1; rc=$?
res S6-pass-1-not-slow-ok-still-78 "0,78,1,1" "$r0,$rc,$(boxes "$d"),$(grep -c "S2 | STOP rc=78: budget: $STOP\$" "$d/NOTES.md")" "$(tail -n 1 "$d/out.S2" | cut -c1-90)"
# the cap still refuses a SLOW_OK stage under the pass: 19.5 used + 0.93 slow estimate passes 20
d=$W/h; conf "$d"; ledger "$d" 9.75 slow slow; drive "$d" P1 SLOW_PAIR_PASS=1; rc=$?
res S7-pass-1-the-cap-still-refuses "78,0,1,0" "$rc,$(boxes "$d"),$(grep -c 'passes the 20 cap; refused' "$d/NOTES.md"),$(grep -c "$NOTE" "$d/NOTES.md")" ""
# no slow pair: the pass changes nothing and writes no note
d=$W/i; conf "$d"; ledger "$d" 0.1 slow fast; drive "$d" P1 SLOW_PAIR_PASS=1; rc=$?
res S8-no-slow-pair-no-note "0,1,0" "$rc,$(boxes "$d"),$(grep -c "$NOTE" "$d/NOTES.md")" ""
# a slow host of its own on a passed SLOW_OK stage runs, and the ledger records it slow
d=$W/j; conf "$d"; ledger "$d" 0.1 slow slow; drive "$d" P1 SLOW_PAIR_PASS=1 SR_HOST=slow; rc=$?
res S9-a-slow-host-of-its-own-is-booked-slow "0,yes,slow" "$rc,$(yn test -e "$d/stages/P1/DONE"),$(tail -n 1 "$d/budget.tsv" | cut -f5)" ""
# any other value is refused with 2 before any box, from the environment or the run file
d=$W/k1; conf "$d"; ledger "$d" 0.1 slow slow; drive "$d" P1 SLOW_PAIR_PASS=yes; r1=$?
d=$W/k2; conf "$d" SLOW_PAIR_PASS=2; ledger "$d" 0.1 fast; drive "$d" P1; r2=$?
d=$W/k3; conf "$d"; ledger "$d" 0.1 fast; drive "$d" P1 "SLOW_PAIR_PASS= 1"; r3=$?
res S10-bad-values-2-no-box "2,2,2,0" "$r1,$r2,$r3,$(( $(boxes "$W/k1") + $(boxes "$W/k2") + $(boxes "$W/k3") ))" ""
res S10-each-names-its-value "1,1,1" "$(grep -c 'SLOW_PAIR_PASS is 0 or 1, not yes$' "$W/k1/out.P1"),$(grep -c 'SLOW_PAIR_PASS is 0 or 1, not 2$' "$W/k2/out.P1"),$(grep -c 'SLOW_PAIR_PASS is 0 or 1, not  1$' "$W/k3/out.P1")" ""

echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
