#!/usr/bin/env bash
# run_species.sh KIT SPECIES LINEAGE STUDY OUT [JOBS]
#
# The MIC models of one species on the lineage and study splits
# (docs/mic-plan.md), on a fresh box with no agent attached: no prompt,
# nothing read from another machine but the public sources below, one log,
# and a restart that skips every finished step and every finished result.
# Run it from a clean checkout of this branch:
#
#   setsid nohup bash scripts/c6-mic/run_species.sh KIT SPECIES LINEAGE STUDY OUT 16 > /dev/null 2>&1 < /dev/null &
#
# LINEAGE and STUDY are comma lists of the species' dataset ids, or - for
# none. KIT holds these files, checked against the sha256 pinned here before
# anything else runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   MASTER-MANIFEST.sha256      the data root's list of files and sha256
#   curated.tar, curated.sha256 the rebuilt revision 14 table's curated files
#   fold-plans.tar.xz, study-plans.sha256 (only with a STUDY list)
#                               the five reviewed study plans per dataset and
#                               the sha256 of each, taken before any fit
#
# Steps, each once: the release; grm built from this checkout, with its
# commit and sha256; the data root (the curated table; the catalog's three
# files from NCBI, each equal to the master list; the genomes from BV-BRC,
# whose list equals the master list and each of which has the release's
# length, contig count and N50, since a refetched file can differ in its gzip
# bytes; then count, structure and matrix with no network); the study plan for
# seed 47 (scripts/c6-mic/study_plans.py); the fits; the probe; the scorer.
#
# The probe: the rebuilt table can give a dataset a breakpoint fitted on one
# more MIC pair than the release's, which moves the input fingerprint and
# nothing the fits read. For each such dataset (release_specs.py) seed 42 is
# refitted on a copy of the data root whose breakpoints are the release's,
# into results/probe; it must give the release's fingerprint, and the scorer
# needs it to equal the run's seed-42 result in all else.
#
# Fits: grm learn --exp-mic on the lineage split at seeds 42 to 47 (released
# and ordered at 42 to 44, released at 45 to 47), one run per seed over the
# datasets whose result is missing; and with --fold-file on the six study
# plans, one run per dataset and plan. After each run every result it should
# have written must exist: grm learn writes no file for a dataset it skips
# and still exits 0. A restart fits only the missing results, and refuses to
# go on if the checkout has moved since the build.
#
# OUT gets logs/run.log, logs/<step>.time and .log, logs/disk.tsv, flags/,
# work/ (release, root, build, plans), results/lineage and results/study (one
# result per dataset and plan, with a row per tested genome and no genome id;
# they stay on the box), score/checks.tsv and score/scores.tsv (aggregates),
# binary.txt and expect.tsv. The end writes SHA256SUMS and OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> SHA256SUMS <sha256>
# The binary's sha256 changes with the checkout's path, which its debug
# information records; the code sha256 is that of its loaded sections less
# the build-id note, which is the same from any path.
# Any failure writes OUT/FAILED naming the step and exits 1. It stops
# nothing but itself.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; LINEAGE=${3:?lineage datasets or -}
STUDY=${4:?study datasets or -}; OUT=${5:?out folder}; JOBS=${6:-$(nproc)}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; PLANS=$W/plans
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
# (even when the retry succeeds), a fold file check, a copy refusal and a
# count failure, and logs/ is copied back. So each step's and fit's log moves
# to work/raw-logs, and logs/ gets it back with every genome id of the
# species' list written as [genome id]; a failure's last 30 lines in
# logs/run.log come from that copy.
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
export KIT SP LINEAGE STUDY OUT JOBS REPO W ROOT R M GRM PLANS

STEP=kit
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the binary must be its commit"; false; }
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/curated.tar" d77e2c6f75f8d986268d0c075cb2ea364fba9972fa22a55139495ac0eb376d97
pinned "$KIT/curated.sha256" 8b0e43dd98ca1b3fc7609b7e25b1b52d4bb4dd0c4c555cf5e771622d24fd9c6e
if [ "$STUDY" != - ]; then
  pinned "$KIT/fold-plans.tar.xz" b13cad7dbf6318fd175a1d4a1fb1715e7763f19030b48338cbd2cd3f14ca71ec
  pinned "$KIT/study-plans.sha256" c9c902f2c9e9b808b8867aff129af1ea736e52dedc93ae1d9e3e631cf39654a1
fi
python3 -I -B "$REPO/scripts/c6-mic/expect.py" "$SP" "$LINEAGE" "$STUDY" > "$OUT/expect.tsv.part"
if [ -s "$OUT/expect.tsv" ] && ! cmp -s "$OUT/expect.tsv" "$OUT/expect.tsv.part"; then
  say "the dataset lists differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"
say "expect: $(($(wc -l < "$OUT/expect.tsv") - 1)) results"

release() {
  rm -rf "$W/release"; mkdir -p "$W/release"
  tar --zstd -xf "$KIT/rev14-274-20261112.tar.zst" -C "$W/release"
  test -f "$R/species/$SP.json"
  for d in $(echo "$LINEAGE,$STUDY" | tr ',' '\n' | grep -vx -- '-' | sort -u); do
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
  # the input fingerprint; the MIC experiment writes no report.
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

plans() {
  rm -rf "$PLANS"; mkdir -p "$PLANS/reviewed" "$PLANS/s47"
  echo "$STUDY" | tr ',' '\n' > "$PLANS/study.txt"
  names=()
  while read -r d; do
    for p in study-whole study-d1 study-d2 study-d3 study-d4; do names+=("folds/$d.$p.json"); done
  done < "$PLANS/study.txt"
  tar -xJf "$KIT/fold-plans.tar.xz" -C "$PLANS/reviewed" "${names[@]}"
  for f in "${names[@]}"; do
    want=$(grep -m1 "/$f\$" "$KIT/study-plans.sha256" | cut -c1-64)
    [ -n "$want" ] && [ "$(sha256sum "$PLANS/reviewed/$f" | cut -c1-64)" = "$want" ] \
      || { echo "$f differs from the reviewed plan's sha256"; exit 1; }
  done
  echo "reviewed study plans: ${#names[@]}, each equal to its sha256"
  python3 -I -B "$REPO/scripts/c6-mic/study_plans.py" "$R" "$PLANS/reviewed/folds" "$PLANS/study.txt" "$PLANS/s47"
  cp "$PLANS/s47/study-plans.tsv" "$OUT/results/study-plans.tsv"
}
if [ "$STUDY" != - ]; then step plans plans; fi

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
# The expected results of a split and seed: file<TAB>dataset<TAB>arms, one a line.
expected() { awk -F'\t' -v s="$1" -v e="$2" 'NR > 1 && $3 == s && $5 == e { print $7 "\t" $1 "\t" $6 }' "$OUT/expect.tsv"; }
missing_of() { while IFS=$'\t' read -r file d arms; do [ -e "$OUT/results/$file" ] || echo "$d"; done; }

for seed in 42 43 44 45 46 47; do
  STEP=lineage-s$seed
  todo=$(expected lineage "$seed" | missing_of | paste -sd, -)
  [ -n "$todo" ] || { say "lineage-s$seed: all results present, skipped"; continue; }
  arms=$(expected lineage "$seed" | awk -F'\t' 'NR == 1 { print $3 }')
  fit "lineage-s$seed" --only "$todo" --seed "$seed" --exp-mic "$arms" --exp-mic-out "$OUT/results/lineage"
  left=$(expected lineage "$seed" | missing_of | paste -sd, -)
  [ -z "$left" ] || { say "lineage-s$seed: no result written for $left"; false; }
  disk "lineage-s$seed"
done

for seed in 42 43 44 45 46 47; do
  while IFS=$'\t' read -r file d arms; do
    plan=$(basename "$file" .mic.json); plan=${plan#"$d."}
    STEP=study-$d-$plan
    if [ -e "$OUT/results/$file" ]; then say "$d $plan: result present, skipped"; continue; fi
    if [ "$plan" = c6mic-s47 ]; then ff=$PLANS/s47/$d.$plan.json; else ff=$PLANS/reviewed/folds/$d.$plan.json; fi
    fit "study-$d-$plan" --fold-file "$ff" --fold-out "$OUT/results/study" --exp-mic "$arms"
    [ -e "$OUT/results/$file" ] || { say "$d $plan: no result written"; false; }
  done < <(expected study "$seed")
done
disk fits

probe() {
  rm -rf "$W/probe" "$OUT/results/probe"; mkdir -p "$W/probe/root/curated" "$OUT/results/probe"
  all=$(echo "$LINEAGE,$STUDY" | tr ',' '\n' | grep -vx -- '-' | sort -u | paste -sd, -)
  python3 -I -B "$REPO/scripts/c6-mic/release_specs.py" "$ROOT/curated/datasets.json" "$R" "$all" "$W/probe"
  cp "$W/probe/release-specs.tsv" "$OUT/results/release-specs.tsv"
  todo=$(cat "$W/probe/differ.txt")
  if [ -z "$todo" ]; then echo "probe: every breakpoint is the release's, none needed"; return 0; fi
  for e in "$ROOT"/*; do [ "$(basename "$e")" = curated ] || ln -s "$e" "$W/probe/root/"; done
  for f in "$ROOT"/curated/*; do [ "$(basename "$f")" = datasets.json ] || ln -s "$f" "$W/probe/root/curated/"; done
  cp "$W/probe/datasets.json" "$W/probe/root/curated/datasets.json"
  arms=$(awk -F'\t' 'NR > 1 && $3 == "lineage" && $5 == 42 { print $6; exit }' "$OUT/expect.tsv")
  "$GRM" --data "$W/probe/root" learn --species "$SP" --only "$todo" --seed 42 --exp-mic "${arms:-released,ordered}" --exp-mic-out "$OUT/results/probe"
  python3 -I -c '
import json, os, sys
probe, release, todo = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
bad = 0
for d in todo:
    path = os.path.join(probe, d + ".lineage-s42.mic.json")
    with open(os.path.join(release, "datasets", d + ".json")) as f:
        rep = json.load(f)
    ok = os.path.isfile(path)
    if ok:
        with open(path) as f:
            doc = json.load(f)
        ok = doc["input_fingerprint"] == rep["input_fingerprint"] and doc["eval_config"] == rep["eval_config"]
    if not ok:
        print(d + ": the probe does not give the release fingerprint, so the inputs differ in more than the breakpoint")
    bad += not ok
print("probe:", len(todo) - bad, "of", len(todo), "give the release fingerprint")
sys.exit(bad != 0)
' "$OUT/results/probe" "$R" "$todo"
}
step probe probe

STEP=score
say "score: start"
python3 -I -B "$REPO/scripts/c6-mic/score.py" "$OUT/expect.tsv" "$OUT/results" "$R" "$OUT/score" "$OUT/results/probe" > "$OUT/logs/score.log" 2>&1 \
  || { cat "$OUT/logs/score.log"; false; }
cat "$OUT/logs/score.log"

STEP=done
cd "$OUT"
{ find results score -type f -name '*.json' -o -type f -name '*.tsv'; echo binary.txt; echo expect.tsv; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
