#!/usr/bin/env bash
# unflagged.sh MAIN_GRM BRANCH_GRM PILOT TEST_DNA OUT
#
# The check that grm with no C4 flag writes main's bytes (docs/screen-arms.md,
# The unflagged bytes). MAIN_GRM is grm built at main, BRANCH_GRM grm built
# at this branch; PILOT is the OUT folder of a finished run_species.sh run
# (its data root, its species and datasets in lists.txt, its reviewed study
# plans); TEST_DNA is the AMRFinderPlus program's bundled test_dna.fa. Each
# binary runs, on its own copy of a data root and with no network:
#   grm learn --only DATASETS --force --seed 42 (the reports),
#   the same with --exp-pick DIR --exp-pick-inner 0 (the default's results),
#   grm learn --fold-file on each study plan with --exp-pick (the study runs),
# on two roots:
#   A. the pilot's data root (its ksets linked, not copied);
#   B. a root of SYNTHETIC genomes, one per id of the species' list, made by
#      unflagged.py fasta from test_dna.fa's contigs and drawn filler, with
#      the pilot's catalog, labels and BV-BRC records and an empty feature
#      list per genome; grm count, structure and matrix build it, then the
#      same fits.
# unflagged.py compare then holds every file of the two copies equal, the
# run clocks of JSON files apart. Prints counts only; exits 1 if a file
# differs or stands in one copy only. Logs stay in OUT.
set -euo pipefail
MAIN=$(realpath "${1:?grm built at main}"); BR=$(realpath "${2:?grm built at this branch}")
P=$(realpath "${3:?a finished pilot OUT}"); TD=$(realpath "${4:?test_dna.fa}"); E=${5:?out folder}
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SP=$(awk -F'\t' '$1 == "species" { print $2 }' "$P/lists.txt")
D=$(awk -F'\t' '$1 == "datasets" { print $2 }' "$P/lists.txt")
S=$(paste -sd, "$P/work/plans/study.txt")
[ -n "$SP" ] && [ -n "$D" ] || { echo "unflagged.sh: no species or datasets in $P/lists.txt" >&2; exit 2; }
mkdir -p "$E"; E=$(realpath "$E"); cd "$E"
say() { echo "$(date -u +%FT%TZ) $*"; }
say "main grm $(sha256sum "$MAIN" | cut -c1-16); branch grm $(sha256sum "$BR" | cut -c1-16); $(echo "$D" | tr ',' '\n' | wc -l) datasets, $(echo "$S" | tr ',' '\n' | grep -c .) study plans"

fits() {  # fits <bin> <root> <out>
  local bin=$1 root=$2 out=$3 d
  mkdir -p "$out/picks" "$out/study"
  unshare -rn "$bin" --data "$root" learn --species "$SP" --only "$D" --force --seed 42 > "$out/plain.log" 2>&1
  cp -a "$root/species/$SP/results" "$out/reports-plain"
  unshare -rn "$bin" --data "$root" learn --species "$SP" --only "$D" --force --seed 42 \
    --exp-pick "$out/picks" --exp-pick-inner 0 > "$out/pick.log" 2>&1
  for d in $(echo "$S" | tr ',' ' '); do
    unshare -rn "$bin" --data "$root" learn --species "$SP" --fold-file "$P/work/plans/reviewed/folds/$d.study-whole.json" \
      --fold-out "$out/study" --exp-pick "$out/picks" --exp-pick-inner 0 --force > "$out/study-$d.log" 2>&1
  done
}

for who in main branch; do
  bin=$MAIN; [ "$who" = branch ] && bin=$BR
  rm -rf "A-$who"; mkdir -p "A-$who/root/species/$SP"
  for x in bvbrc catalog curated fasta; do cp -a "$P/work/root/$x" "A-$who/root/$x"; done
  for x in "$P/work/root/species/$SP"/*; do
    if [ "$(basename "$x")" = ksets ]; then cp -al "$x" "A-$who/root/species/$SP/ksets"; else cp -a "$x" "A-$who/root/species/$SP/"; fi
  done
  [ ! -f "$P/work/root/species/$SP.json" ] || cp -a "$P/work/root/species/$SP.json" "A-$who/root/species/"
  say "A $who: fits"
  fits "$bin" "A-$who/root" "A-$who/out"
done
say "A: compare"
rcA=0; python3 -I "$HERE/unflagged.py" compare A-main A-branch "root/species/$SP/ksets" || rcA=$?

rm -rf B-template; mkdir -p "B-template/root/species/$SP"
for x in bvbrc catalog curated; do cp -a "$P/work/root/$x" "B-template/root/$x"; done
cp -a "$P/work/root/species/$SP/genomes.txt" "B-template/root/species/$SP/"
[ ! -f "$P/work/root/species/$SP.json" ] || cp -a "$P/work/root/species/$SP.json" "B-template/root/species/"
python3 -I "$HERE/unflagged.py" fasta "$TD" "B-template/root/species/$SP/genomes.txt" B-template/root/fasta
# The genomes are SYNTHETIC, so no BV-BRC annotation is theirs: an empty
# feature list for each, and no download (the fits run with no network).
rm -rf B-template/root/bvbrc/features; mkdir -p B-template/root/bvbrc/features
while read -r g; do [ -z "$g" ] || echo "[]" > "B-template/root/bvbrc/features/$g.json"; done < "B-template/root/species/$SP/genomes.txt"
for who in main branch; do
  bin=$MAIN; [ "$who" = branch ] && bin=$BR
  rm -rf "B-$who"; cp -a B-template "B-$who"
  say "B $who: count, structure, matrix"
  for c in count structure matrix; do unshare -rn "$bin" --data "B-$who/root" "$c" --species "$SP" > "B-$who/$c.log" 2>&1; done
  say "B $who: fits"
  fits "$bin" "B-$who/root" "B-$who/out"
done
say "B: compare"
rcB=0; python3 -I "$HERE/unflagged.py" compare B-main B-branch || rcB=$?
say "reports: A $(find A-branch/out/reports-plain -name '*.json' | wc -l), B $(find B-branch/out/reports-plain -name '*.json' | wc -l); result files: A $(find A-branch/out/picks -type f | wc -l), B $(find B-branch/out/picks -type f | wc -l)"
say "done A=$rcA B=$rcB"
[ "$rcA" = 0 ] && [ "$rcB" = 0 ]
