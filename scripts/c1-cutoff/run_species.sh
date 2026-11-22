#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS] [LINES] [DRAWS] [TWIN]
#
# The small-cohort run of one species (plan block C1 part B,
# docs/cohort-draw.md) on a fresh box with no agent attached: no prompt,
# nothing read from another machine but the public sources below, one log,
# and a restart that skips every finished step and every finished result.
# Run it from a clean checkout of this branch:
#
#   setsid nohup bash scripts/c1-cutoff/run_species.sh KIT SPECIES DATASETS OUT 16 > /dev/null 2>&1 < /dev/null &
#
# DATASETS is a comma list of the species' dataset ids. LINES (default
# 10,15,20,25,50) are the genomes per class, DRAWS (default 5) the draws per
# dataset and line, TWIN (default 1) fits the r 25 cohorts again at seed 45.
# The timing pilot ran with LINES 25, DRAWS 1, TWIN 0. KIT holds these files,
# checked against the sha256 pinned here before anything else runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   MASTER-MANIFEST.sha256      the data root's list of files and sha256
#   curated.tar, curated.sha256 the rebuilt revision 14 table's curated files
#   fold-plans.tar.xz, study-plans.sha256
#                               the reviewed study plans per dataset and the
#                               sha256 of each, taken before any fit
#
# Steps, each once: the release; grm built from this checkout, with its
# commit and sha256; the data root (the curated table; the catalog's three
# files from NCBI, each equal to the master list; the genomes from BV-BRC,
# whose list equals the master list and each of which has the release's
# length, contig count and N50, since a refetched file can differ in its gzip
# bytes; then count, structure and matrix with no network); the study-whole
# plans of the datasets that have one; the draws (draw.py); the fits; the
# scorer.
#
# Fits: grm learn --exp-cohort at seed 42 over every cohort file, and with
# TWIN 1 at seed 45 over the r 25 cohort files. grm learn keeps a result that
# exists, so a restart fits only the missing ones. After each run every
# result it should have written must exist. A restart refuses to go on if the
# checkout, the binary, the dataset list or the drawn cohorts moved since the
# first start.
#
# OUT gets logs/run.log, logs/<step>.time and .log, logs/disk.tsv, flags/,
# work/ (release, root, build, plans, cohorts: the cohort files list genome
# ids and stay on the box), results/cohorts (one result per cohort and seed,
# aggregates with no genome id), results/draws.tsv, score/ (aggregates),
# binary.txt, lists.txt and expect.tsv. The end writes SHA256SUMS and
# OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> SHA256SUMS <sha256>
# The binary's sha256 changes with the checkout's path, which its debug
# information records; the code sha256 is that of its loaded sections less
# the build-id note, which is the same from any path.
# Any failure writes OUT/FAILED naming the step and exits 1. It stops
# nothing but itself.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?datasets}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}; LINES=${6:-10,15,20,25,50}; DRAWS=${7:-5}; TWIN=${8:-1}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; PLANS=$W/plans; COH=$W/cohorts; RES=$OUT/results/cohorts
export RUSTUP_TOOLCHAIN=1.98.1 CARGO_BUILD_JOBS=$JOBS RAYON_NUM_THREADS=$JOBS
export PATH=$HOME/.cargo/bin:/usr/local/cargo/bin:$PATH
# One start at a time on an OUT: a second one would share its steps and files.
exec 9> "$OUT/run.lock"
flock -n 9 || { echo "$(date -u +%FT%TZ) another start runs on this OUT; this one ends"; exit 1; }
STEP=start
say() { echo "$(date -u +%FT%TZ) $*"; }
disk() { echo -e "$(date -u +%FT%TZ)\t$1\t$(du -sb "$OUT" | cut -f1)\t$(df -B1 --output=avail "$OUT" | tail -1 | tr -d ' ')" >> "$OUT/logs/disk.tsv"; }
trap 'say "FAILED at $STEP"; echo "FAILED $SP at $STEP $(date -u +%FT%TZ)" > "$OUT/FAILED"; exit 1' ERR
trap 'say "stopped at $STEP"; echo "FAILED $SP stopped at $STEP $(date -u +%FT%TZ)" > "$OUT/FAILED"; exit 1' TERM INT
rm -f "$OUT/FAILED" "$OUT/DONE"
say "start $SP, $JOBS jobs, repo $(git -C "$REPO" rev-parse HEAD)"
[ -s "$OUT/logs/disk.tsv" ] || echo -e "time\tafter\tout_bytes\tfree_bytes" > "$OUT/logs/disk.tsv"
disk start

pinned() {  # pinned <file> <sha256>
  [ "$(sha256sum "$1" | cut -c1-64)" = "$2" ] || { say "$1: sha256 differs from the pinned value"; return 1; }
}
master() {  # master <path in the data root> <file>
  local want; want=$(grep -m1 "  $1\$" "$M" | cut -c1-64)
  [ -n "$want" ] && [ "$(sha256sum "$2" | cut -c1-64)" = "$want" ]
}
# Genome ids stay on the box. grm names a genome in a download's retry line
# (even when the retry succeeds), a cohort file check and a count failure,
# and logs/ is copied back. So each step's and fit's log moves to
# work/raw-logs, and logs/ gets it back with every genome id of the species'
# list written as [genome id]; a failure's last 30 lines in logs/run.log come
# from that copy.
scrub() {  # scrub <log>
  local raw; raw=$W/raw-logs/$(basename "$1")
  mkdir -p "$W/raw-logs"; mv -f "$1" "$raw"
  python3 -I -c '
import json, os, re, sys
ids = set()
if os.path.isfile(sys.argv[1]):
    with open(sys.argv[1]) as f:
        ids = {g["id"].encode() for g in json.load(f)["genomes"]}
token = re.compile(rb"(?<![\w.])\d+\.\d+(?!\d)")
for line in sys.stdin.buffer:
    sys.stdout.buffer.write(token.sub(lambda m: b"[genome id]" if m.group(0) in ids else m.group(0), line))
' "$R/species/$SP.json" < "$raw" > "$1"
}
step() {  # step <name> <function>: once, timed, in a shell that stops at the first failure
  STEP=$1
  if [ -e "$OUT/flags/$1.done" ]; then say "$1: done before, skipped"; return 0; fi
  say "$1: start"
  /usr/bin/time -f "$1 %e s wall, %M KB peak, %U s user, %S s sys" -o "$OUT/logs/$1.time" \
    bash -euo pipefail -c "$(declare -f say pinned master "$2"); $2" < /dev/null > "$OUT/logs/$1.log" 2>&1 \
    || { scrub "$OUT/logs/$1.log"; tail -n 30 "$OUT/logs/$1.log"; false; }
  scrub "$OUT/logs/$1.log"
  touch "$OUT/flags/$1.done"
  disk "$1"
  say "$1: end, $(cat "$OUT/logs/$1.time")"
}
export KIT SP DATASETS OUT JOBS LINES DRAWS TWIN REPO W ROOT R M GRM PLANS COH RES

STEP=kit
# grm matrix keeps one kset open per genome, up to 2,000 a species; a
# service's soft open-file limit can be 1,024. Raise it to the hard limit.
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
[ "$(ulimit -n)" = unlimited ] || [ "$(ulimit -n)" -ge 4096 ] \
  || { say "the open-file limit is $(ulimit -n), below the 4,096 grm matrix needs here"; false; }
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the binary must be its commit"; false; }
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/curated.tar" d77e2c6f75f8d986268d0c075cb2ea364fba9972fa22a55139495ac0eb376d97
pinned "$KIT/curated.sha256" 8b0e43dd98ca1b3fc7609b7e25b1b52d4bb4dd0c4c555cf5e771622d24fd9c6e
pinned "$KIT/fold-plans.tar.xz" b13cad7dbf6318fd175a1d4a1fb1715e7763f19030b48338cbd2cd3f14ca71ec
pinned "$KIT/study-plans.sha256" c9c902f2c9e9b808b8867aff129af1ea736e52dedc93ae1d9e3e631cf39654a1
for d in $(echo "$DATASETS" | tr ',' '\n'); do
  case $d in "$SP"__*) ;; *) say "$d is not a dataset of $SP"; false ;; esac
done
[ "$(echo "$DATASETS" | tr ',' '\n' | sort | uniq -d | wc -l)" = 0 ] || { say "a dataset is listed twice"; false; }
case $TWIN in 0|1) ;; *) say "TWIN is 0 or 1"; false ;; esac
printf 'species\t%s\ndatasets\t%s\nlines\t%s\ndraws\t%s\ntwin\t%s\n' "$SP" "$DATASETS" "$LINES" "$DRAWS" "$TWIN" > "$OUT/lists.txt.part"
if [ -s "$OUT/lists.txt" ] && ! cmp -s "$OUT/lists.txt" "$OUT/lists.txt.part"; then
  say "the dataset list or the lines differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/lists.txt.part" "$OUT/lists.txt"

release() {
  rm -rf "$W/release"; mkdir -p "$W/release"
  tar --zstd -xf "$KIT/rev14-274-20261112.tar.zst" -C "$W/release"
  test -f "$R/species/$SP.json"
  for d in $(echo "$DATASETS" | tr ',' '\n'); do
    test -f "$R/datasets/$d.json" || { echo "the release has no report $d"; exit 1; }
  done
  echo "reports of $SP in the release: $(ls "$R"/datasets/"$SP"__*.json | wc -l)"
}
step release release

build() {
  cd "$REPO" && CARGO_TARGET_DIR=$W/target cargo build --release --locked --bin grm
  echo "commit $(git -C "$REPO" rev-parse HEAD)" > "$OUT/binary.txt"
  echo "grm $(sha256sum "$GRM" | cut -c1-64)" >> "$OUT/binary.txt"
  objcopy -O binary --remove-section=.note.gnu.build-id "$GRM" "$W/grm.image"
  echo "code $(sha256sum "$W/grm.image" | cut -c1-64)" >> "$OUT/binary.txt"
  rm -f "$W/grm.image"
  cat "$OUT/binary.txt"
}
step build build
STEP=build-check
[ "$(sed -n 's/^commit //p' "$OUT/binary.txt")" = "$(git -C "$REPO" rev-parse HEAD)" ] \
  || { say "the checkout moved since the build; use a fresh OUT"; false; }
[ "$(sed -n 's/^grm //p' "$OUT/binary.txt")" = "$(sha256sum "$GRM" | cut -c1-64)" ] \
  || { say "the binary changed since the build"; false; }

curated() {
  mkdir -p "$ROOT/curated" "$ROOT/bvbrc"
  tar -xf "$KIT/curated.tar" -C "$ROOT/curated"
  (cd "$ROOT/curated" && sha256sum --quiet -c "$KIT/curated.sha256")
  # Genome metadata feeds only a report's descriptive fields and is not in
  # the input fingerprint; the cohort experiment writes no report.
  : > "$ROOT/bvbrc/genomes.jsonl"
}
step curated curated

catalog() {
  mkdir -p "$ROOT/catalog"
  base=https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/2026-08-07.1
  for f in AMR_CDS.fa ReferenceGeneCatalog.txt version.txt; do
    [ -s "$ROOT/catalog/$f" ] || { curl --fail --silent --show-error --location --retry 3 -A curl -o "$ROOT/catalog/$f.part" "$base/$f"; mv "$ROOT/catalog/$f.part" "$ROOT/catalog/$f"; }
    master "catalog/$f" "$ROOT/catalog/$f" || { echo "catalog/$f differs from the master list"; exit 1; }
  done
  echo "catalog: 3 of 3 files equal the master list"
}
step catalog catalog

genomes() {
  mkdir -p "$ROOT/species/$SP" "$ROOT/fasta"
  python3 -I -c 'import json,sys; [print(g["id"]) for g in json.load(open(sys.argv[1]))["genomes"]]' \
    "$R/species/$SP.json" > "$ROOT/species/$SP/genomes.txt"
  master "species/$SP/genomes.txt" "$ROOT/species/$SP/genomes.txt"
  "$GRM" --data "$ROOT" fetch fasta --species "$SP" --workers 4
  n=0; bad=0
  while read -r id; do
    n=$((n + 1)); master "fasta/$id.fna.gz" "$ROOT/fasta/$id.fna.gz" || bad=$((bad + 1))
  done < "$ROOT/species/$SP/genomes.txt"
  # A refetched file can differ in its bytes from the master list. Every
  # genome must still have the release's length, contig count and N50.
  echo "genomes: $n, bytes differ from the master list: $bad"
  python3 -I -c '
import gzip, json, sys
gs = json.load(open(sys.argv[1]))["genomes"]; bad = 0
for g in gs:
    lens = []
    for line in gzip.open(sys.argv[2] + "/" + g["id"] + ".fna.gz", "rt"):
        if line.startswith(">"):
            lens.append(0)
        else:
            lens[-1] += len(line.strip())
    tot, acc, n50 = sum(lens), 0, -1
    for x in sorted(lens, reverse=True):
        acc += x
        if 2 * acc >= tot:
            n50 = x
            break
    bad += (tot, len(lens), n50) != (g["length"], g["contigs"], g["n50"])
print("assembly stats differ from the release:", bad, "of", len(gs))
sys.exit(bad != 0)
' "$R/species/$SP.json" "$ROOT/fasta"
}
step genomes genomes

count() { unshare -rn "$GRM" --data "$ROOT" count --species "$SP"; }
step count count
structure() { unshare -rn "$GRM" --data "$ROOT" structure --species "$SP"; }
step structure structure
matrix() { unshare -rn "$GRM" --data "$ROOT" matrix --species "$SP"; }
step matrix matrix

# The study-whole plan of each listed dataset that has a reviewed one.
plans() {
  rm -rf "$PLANS"; mkdir -p "$PLANS/reviewed"
  names=()
  for d in $(echo "$DATASETS" | tr ',' '\n'); do
    if grep -q "/folds/$d.study-whole.json\$\|  folds/$d.study-whole.json\$" "$KIT/study-plans.sha256"; then
      names+=("folds/$d.study-whole.json")
    fi
  done
  if [ "${#names[@]}" = 0 ]; then mkdir -p "$PLANS/reviewed/folds"; echo "study plans: none of the datasets has one"; return 0; fi
  tar -xJf "$KIT/fold-plans.tar.xz" -C "$PLANS/reviewed" "${names[@]}"
  for f in "${names[@]}"; do
    want=$(grep -m1 "/$f\$\|  $f\$" "$KIT/study-plans.sha256" | cut -c1-64)
    [ -n "$want" ] && [ "$(sha256sum "$PLANS/reviewed/$f" | cut -c1-64)" = "$want" ] \
      || { echo "$f differs from the reviewed plan's sha256"; exit 1; }
  done
  echo "study plans: ${#names[@]}, each equal to its sha256"
}
step plans plans

draws() {
  rm -rf "$COH"; mkdir -p "$COH"
  python3 -I -B "$REPO/scripts/c1-cutoff/draw.py" "$R" "$DATASETS" "$LINES" "$DRAWS" "$COH" "$PLANS/reviewed/folds"
  cp "$COH/draws.tsv" "$OUT/results/draws.tsv"
  # The twin's cohorts: the r 25 files, by link.
  mkdir -p "$COH/twin"
  for f in "$COH"/*.r25.d*.cohort.json; do [ -e "$f" ] && ln -s "../$(basename "$f")" "$COH/twin/"; done
  true
}
step draws draws
STEP=expect
python3 -I -B "$REPO/scripts/c1-cutoff/expect.py" "$SP" "$OUT/results/draws.tsv" "$TWIN" > "$OUT/expect.tsv.part"
if [ -s "$OUT/expect.tsv" ] && ! cmp -s "$OUT/expect.tsv" "$OUT/expect.tsv.part"; then
  say "the expected results differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"
say "expect: $(($(wc -l < "$OUT/expect.tsv") - 1)) results"
# Each cohort file must still be the one drawn.
while IFS=$'\t' read -r d _ _ _ _ _ cohort sha _; do
  [ "$d" = dataset ] && continue
  [ "$(sha256sum "$COH/$cohort" | cut -c1-64)" = "$sha" ] || { say "a cohort file of $d changed since the draw"; false; }
done < "$OUT/expect.tsv"

# One fit run: fit <name> <grm learn arguments...>, timed into logs/<name>.time.
fit() {
  local name=$1; shift
  STEP=$name
  say "$name: start"
  /usr/bin/time -v -o "$OUT/logs/$name.time" "$GRM" --data "$ROOT" learn --species "$SP" "$@" \
    < /dev/null > "$OUT/logs/$name.log" 2>&1 || { scrub "$OUT/logs/$name.log"; tail -n 30 "$OUT/logs/$name.log"; false; }
  scrub "$OUT/logs/$name.log"
  say "$name: end, $(grep -E 'Elapsed|Maximum resident' "$OUT/logs/$name.time" | sed 's/^\s*//' | tr '\n' ';')"
}
# The datasets with a missing result at a seed, as a comma list.
missing() { awk -F'\t' -v s="$1" 'NR > 1 && $5 == s { print $1 "\t" $9 }' "$OUT/expect.tsv" \
  | while IFS=$'\t' read -r d file; do [ -e "$RES/$file" ] || echo "$d"; done | sort -u | paste -sd, -; }

mkdir -p "$RES"
for seed in 42 45; do
  STEP=fit-s$seed
  if [ "$seed" = 45 ] && [ "$TWIN" = 0 ]; then continue; fi
  todo=$(missing "$seed")
  [ -n "$todo" ] || { say "fit-s$seed: all results present, skipped"; continue; }
  if [ "$seed" = 42 ]; then dir=$COH; else dir=$COH/twin; fi
  fit "fit-s$seed" --only "$todo" --seed "$seed" --exp-cohort "$dir" --exp-cohort-out "$RES"
  left=$(missing "$seed")
  [ -z "$left" ] || { say "fit-s$seed: no result written for $left"; false; }
  disk "fit-s$seed"
done

STEP=score
say "score: start"
python3 -I -B "$REPO/scripts/c1-cutoff/score.py" "$OUT/expect.tsv" "$COH" "$RES" "$R" "$OUT/score" > "$OUT/logs/score.log" 2>&1 \
  || { cat "$OUT/logs/score.log"; false; }
cat "$OUT/logs/score.log"

STEP=done
cd "$OUT"
{ find results score -type f \( -name '*.json' -o -name '*.tsv' \); echo binary.txt; echo expect.tsv; echo lists.txt; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
