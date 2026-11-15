#!/bin/bash
# Lane clinda: scratch data roots and grm's own stage commands (PLAN.md section 5; addendum A2, A3, A4).
#
#   arms.sh make-root-with
#   arms.sh stages  <root> <species> <threads>     grm count, structure, matrix (defaults: k 31, min support 2)
#   arms.sh learn   <root> <species> <threads>     grm learn --only <species>__clindamycin
#   arms.sh make-root-without <name> <species> <removal ids>
#   arms.sh guard   <species>                      sha256 of every built file of root-with/species/<species>
#
# Every grm call runs under `unshare -rn` (no network). The public root is only read.
set -euo pipefail
S=/home/user/scratch/clinda
PUB=/home/user/grm-h1raw
GRM=$S/target/release/grm
OUT=/home/user/grm-ops/program-20261114/clinda
W=$S/root-with
export RUSTUP_TOOLCHAIN=1.98.1
SPECIES="campylobacter_coli campylobacter_jejuni enterococcus_faecium staphylococcus_aureus streptococcus_agalactiae streptococcus_pneumoniae"

cmd=$1; shift
case $cmd in
make-root-with)
  mkdir -p "$W"
  for e in bindings bvbrc catalog fasta species-check; do cp -rs "$PUB/$e" "$W/$e"; done
  cp -r "$PUB/curated" "$W/curated"
  for sp in $SPECIES; do mkdir -p "$W/species/$sp"; cp "$PUB/species/$sp/genomes.txt" "$W/species/$sp/"; done
  ;;
stages)
  R=$1; sp=$2; th=$3
  export RAYON_NUM_THREADS=$th
  for st in count structure matrix; do
    date -u +"$sp $st start %FT%TZ threads $th"
    unshare -rn "$GRM" --data "$R" "$st" --species "$sp"
    date -u +"$sp $st end %FT%TZ"
  done
  ;;
learn)
  R=$1; sp=$2; th=$3
  date -u +"$sp learn $R start %FT%TZ threads $th"
  RAYON_NUM_THREADS=$th unshare -rn "$GRM" --data "$R" learn --species "$sp" --only "${sp}__clindamycin"
  date -u +"$sp learn $R end %FT%TZ"
  ;;
make-root-without)
  name=$1; sp=$2; ids=$3
  R=$S/roots/$name
  mkdir -p "$R/curated" "$R/species/$sp/results"
  for e in bindings bvbrc catalog fasta species-check; do ln -s "$W/$e" "$R/$e"; done
  for f in "$W"/curated/*; do b=$(basename "$f"); [ "$b" = phenotypes.jsonl ] || cp "$f" "$R/curated/$b"; done
  python3 -I "$OUT/scripts/filter_phenotypes.py" "$W/curated/phenotypes.jsonl" "$ids" "$R/curated/phenotypes.jsonl"
  for f in "$W/species/$sp"/*; do b=$(basename "$f"); [ "$b" = results ] || ln -s "$f" "$R/species/$sp/$b"; done
  ;;
guard)
  sp=$1
  (cd "$W/species/$sp" && find . -type f -not -path './results/*' -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d' ' -f1)
  ;;
*)
  echo "unknown command $cmd" >&2; exit 2
  ;;
esac
