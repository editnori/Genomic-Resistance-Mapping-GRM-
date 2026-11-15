#!/usr/bin/env bash
# Builds one species' count, structure and matrix in the scratch data root with grm at 852ebc9,
# timing each stage with /usr/bin/time -v, then removes the ksets (as `grm pipeline --drop-ksets` does).
# usage: build-species.sh SLUG THREADS
set -u -o pipefail
S=$1 T=$2
SC=/home/user/scratch/lassoscan
GRM=$SC/target/release/grm D=$SC/data L=$SC/logs/build-$S
mkdir -p "$L"
export RAYON_NUM_THREADS=$T
{ echo "start $(date -u +%FT%TZ) species=$S threads=$T nproc=$(nproc) load=$(cut -d' ' -f1-3 /proc/loadavg)"; sha256sum "$GRM"; } > "$L/stages.log"
rc=0
for stage in count structure matrix; do
  echo "$stage start $(date -u +%FT%TZ)" >> "$L/stages.log"
  /usr/bin/time -v -o "$L/$stage.time" "$GRM" --data "$D" "$stage" --species "$S" > "$L/$stage.out" 2>&1
  rc=$?
  echo "$stage end $(date -u +%FT%TZ) rc=$rc maxrss_kb=$(awk -F': ' '/Maximum resident/ {print $2}' "$L/$stage.time")" >> "$L/stages.log"
  [ $rc = 0 ] || break
done
if [ $rc = 0 ]; then
  du -sb "$D/species/$S/ksets" >> "$L/stages.log" 2>&1
  rm -rf "$D/species/$S/ksets" && echo "ksets removed $(date -u +%FT%TZ)" >> "$L/stages.log"
  du -sb "$D/species/$S" >> "$L/stages.log"
fi
echo "rc=$rc $(date -u +%FT%TZ)" > "$L/done"
