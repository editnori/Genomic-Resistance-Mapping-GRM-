#!/usr/bin/env bash
# KIT11: drive.sh's KIND=h1 with the data on a snapshot (DATA_SNAPSHOT_DIR), on the stub runner and a stub boat (no
# box; every byte synthetic). The same manifests and pins as push-data mode reach the runner as three stage-local
# calls, from the snapshot's folders, to the same destinations; the box starts from the data snapshot; DATA_ROOT and
# VOCAB_ROOT aren't read; the mode's settings are refused before any box when they are wrong; a combine box still
# gets its parts through push-data. One "RESULT id want got verdict detail" line per check, then "FAILS n"; exits 1
# on any failure.
# usage: test-drive-h1-kit11.sh WORK DRIVE_SH KIT_DIR
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
  printf 'C1\tcombine\tdefault\t0.25\t0.4\t1\t$RUN_DIR/stages/P1/DONE,$RUN_DIR/stages/S2/DONE\t%s\t1.0\tP1,S2\t2\n' "$S"; } > "$W/stages.tsv"
SNAP=grm-tmpl-h1data-7dd00779
conf() { # dir mode(snap|push) extra...
  local d=$1 m=$2; shift 2; mkdir -p "$d"
  cat > "$d/run.conf" << EOF
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
H1_TRACE_BOXES=P1,S2
LANE_PREFIX=k11t-
RAM_NEED_GIB=13.5
POLL=1 CKPT_SECS=1800
BUDGET_STOP=14 BUDGET_CAP=20
GRM_SCP_FLOOR_MBPS=5
RUN_DIR=$d
NOTES=$d/NOTES.md
COPY2=$d/copy2
EOF
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

# the push-data run, for comparison: its three manifests are what the snapshot run must hand stage-local
p=$W/push; conf "$p" push; drive "$p" P1; rc=$?
res D0-push-data-mode-unchanged "0,3,0" "$rc,$(grep -c '^push-data ' "$p/calls.log"),$(grep -c '^stage-local ' "$p/calls.log")" "$(tail -n 1 "$p/out.P1")"
# DATA_ROOT and VOCAB_ROOT point nowhere here, so a read of either would fail
d=$W/a; conf "$d" snap "DATA_ROOT=$W/no-such-root" "VOCAB_ROOT=$W/no-such-vocab"; drive "$d" P1; rc=$?
res D1-snapshot-pilot-rc 0 "$rc" "$(tail -n 1 "$d/out.P1")"
res D1-three-stage-local-no-push-data "grm-meta,grm-vocab,grm-fasta,0" "$(grep -o '^stage-local [^|]*--dest /home/user/grm-[a-z]*' "$d/calls.log" | grep -o 'grm-[a-z]*$' | paste -sd, -),$(grep -c '^push-data ' "$d/calls.log")" ""
res D1-from-the-snapshot-folders "/home/user/grm-h1raw,/home/user/grm-h1vocab,/home/user/grm-h1raw" "$(cat "$d/calls.log.slfrom-grm-meta"),$(cat "$d/calls.log.slfrom-grm-vocab"),$(cat "$d/calls.log.slfrom-grm-fasta")" ""
res D1-same-manifests-as-push-data "yes,yes,yes" "$(yn cmp -s "$d/calls.log.sl-grm-meta" "$p/calls.log.pd-grm-meta"),$(yn cmp -s "$d/calls.log.sl-grm-vocab" "$p/calls.log.pd-grm-vocab"),$(yn cmp -s "$d/calls.log.sl-grm-fasta" "$p/calls.log.pd-grm-fasta")" ""
res D1-pins-passed "3,1,1" "$(grep '^stage-local ' "$d/calls.log" | grep -c -- '--manifest-sha256 [0-9a-f]\{64\} '),$(grep '^stage-local .*grm-vocab' "$d/calls.log" | grep -c -- "--manifest-sha256 $VSHA "),$(grep -c 'FASTA manifest pinned' "$d/NOTES.md")" ""
res D1-box-from-the-data-snapshot yes "$(yn grep -q -- "--from-template $SNAP" "$d/created")" "$(cut -d' ' -f3- "$d/created")"
res D1-vocabulary-pin-names-its-file 1 "$(grep -c "vocabulary manifest pinned: .* from $W/VOCAB-MANIFEST.sha256 (the files at /home/user/grm-h1vocab in $SNAP)" "$d/NOTES.md")" ""
res D1-stage-env-pins-match-push-mode "yes" "$(yn bash -c "grep -o 'H1_\(META\|VOCAB\)_SHA256=[0-9a-f]*' '$d/calls.log.env' | sort | cmp -s - <(grep -o 'H1_\(META\|VOCAB\)_SHA256=[0-9a-f]*' '$p/calls.log.env' | sort)")" ""
res D1-done-closed "yes,1" "$(yn test -e "$d/stages/P1/DONE"),$(cnt "$d/closed")" ""
drive "$d" S2; rc=$?
drive "$d" C1; rc2=$?
res D2-combine-still-pushes-its-parts "0,0,1,0" "$rc,$rc2,$(grep -c '^push-data .*--dest /home/user/grm-combine' "$d/calls.log"),$(grep 'c1-' "$d/created" | grep -c -- '--from-template')" "(S2 rc, C1 rc, combine push-data calls, combine box templates)"
d=$W/b; conf "$d" snap "KIND=r11"; drive "$d" P1; rc=$?
res D3-snapshot-mode-is-h1-only "2,0,1" "$rc,$(boxes "$d"),$(grep -c 'is for KIND=h1 only' "$d/out.P1")" ""
d=$W/c; conf "$d" snap "TEMPLATE="; drive "$d" P1; r1=$?; c1=$(boxes "$d")
d=$W/c2; conf "$d" snap "TEMPLATE=some-other-snapshot"; drive "$d" P1; r2=$?; c2=$(boxes "$d")
res D4-snapshot-mode-needs-a-grm-tmpl-template "2,0,2,0,2" "$r1,$c1,$r2,$c2,$(cat "$W/c/out.P1" "$W/c2/out.P1" | grep -c 'needs TEMPLATE, the grm-tmpl-')" ""
d=$W/e; conf "$d" snap "VOCAB_MANIFEST=$W/no-such-manifest"; drive "$d" P1; rc=$?
res D5-missing-vocab-manifest-2-no-box "2,0,1" "$rc,$(boxes "$d"),$(grep -c 'needs VOCAB_MANIFEST' "$d/out.P1")" "$(tail -n 1 "$d/out.P1" | cut -c1-90)"
printf '%s  vocabulary/a.tsv\n' "$(printf 'e%.0s' {1..64})" > "$W/other-vocab.sha256"
d=$W/f; conf "$d" snap "VOCAB_MANIFEST=$W/other-vocab.sha256"; drive "$d" P1; rc=$?
res D6-other-vocab-manifest-2-no-box "2,0,1" "$rc,$(boxes "$d"),$(grep -c "not the pinned $VSHA" "$d/NOTES.md")" ""
d=$W/g; conf "$d" snap "DATA_SNAPSHOT_DIR=/tmp/raw"; drive "$d" P1; r1=$?
d=$W/g2; conf "$d" snap "VOCAB_SNAPSHOT_DIR=relative/vocab"; drive "$d" P1; r2=$?
res D7-snapshot-folders-under-home-user "2,2,2" "$r1,$r2,$(cat "$W/g/out.P1" "$W/g2/out.P1" | grep -c 'must be a box folder under /home/user')" ""
d=$W/h; conf "$d" push "DATA_ROOT="; drive "$d" P1; rc=$?
res D8-push-mode-still-needs-data-root "2,0,1" "$rc,$(boxes "$d"),$(grep -c 'DATA_ROOT is unset' "$d/out.P1")" ""
d=$W/i; conf "$d" snap; drive "$d" P1 SR_STAGELOCAL_FAIL=fasta; rc=$?
res D9-stage-local-refusal-stops-and-closes "97,1,yes" "$rc,$(cnt "$d/closed"),$(yn grep -q 'stage-local of the FASTAs refused or failed' "$d/NOTES.md")" ""
d=$W/j; conf "$d" snap; drive "$d" P1; rc=$?
res D10-snapshot-mode-needs-no-data-or-vocab-root "0,3" "$rc,$(grep -c '^stage-local ' "$d/calls.log")" "(no DATA_ROOT or VOCAB_ROOT line in the run file)"

echo "FAILS $FAILS"
[ "$FAILS" = 0 ]
