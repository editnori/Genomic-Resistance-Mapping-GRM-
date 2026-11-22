#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARMS]
#
# The input-arm run of one species (plan block C3, docs/input-arms.md) on a
# fresh box with no agent attached: no prompt, nothing read from another
# machine but the public sources below, one log, and a restart that skips
# every finished step and every finished result. Run it from a clean checkout
# of this branch:
#
#   setsid nohup bash scripts/c3-inputs/run_species.sh KIT SPECIES DATASETS OUT 16 > /dev/null 2>&1 < /dev/null &
#
# DATASETS is a comma list of the species' dataset ids that are not marked
# M. The program arms run on those of them in A (scripts/c3-inputs/lists.tsv,
# in_A). SEEDS (default 42,43,44) are the draws of the default (k31) and of
# each arm; TWIN (default 45,46,47, or - for none) the twin's seeds, the
# default fitted again. ARMS (default program,k31-program, the arms built so
# far) are the arms fitted beside the default. Every fit runs with
# --exp-pick-inner 0: each outer fold's calls per learner, no inner pick.
# SEEDS must hold 42, the release's seed: the default's seed-42 reports are
# checked against the release's before any arm runs. KIT holds these files,
# checked against the sha256 pinned here before anything else runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   MASTER-MANIFEST.sha256      the data root's list of files and sha256
#   curated.tar, curated.sha256 the rebuilt revision 14 table's curated files
#   fold-plans.tar.xz, study-plans.sha256
#                               the reviewed study plans per dataset and the
#                               sha256 of each, taken before any fit
#   marks.tsv                   the label count's table: a dataset is M when
#                               pairs_S+I or pairs_I+R is above 0
#   amrkit.tar, tool-tables.sha256
#                               the AMRFinderPlus kit (AMRFinderPlus 4.2.7,
#                               database 2026-08-07.1) and the sha256 of each
#                               program table the release recorded
#
# Steps, each once: the release; grm built from this checkout, with its
# commit, sha256 and code sha256; the data root (the curated table with the
# release's breakpoints; the catalog's three files from NCBI, each equal to
# the master list; the genomes from BV-BRC, whose list equals the master list
# and each of which has the release's length, contig count and N50; then
# count, structure and matrix with no network); the study-whole plans of the
# datasets that have one; then, when a program arm runs, the program's
# environment (micromamba, the pinned environment, the database, the kit's
# own tests and its fixture tables) and the program on every genome of the
# species with a recorded table, each table equal to its recorded sha256;
# expect.tsv (expect.py).
#
# Fits, per seed of SEEDS: the default, grm learn --seed s --exp-pick
# results/picks/k31, then each arm, grm learn --seed s --exp-inputs ARM
# --exp-program-tables <the tables> --exp-pick results/picks/ARM, each over
# the datasets whose result at that seed is missing (--only, --force). After
# the default's seed-42 run, scripts/refit/refit_check.py compares each
# seed-42 report with the release's; any difference stops the run before
# another fit. Then the twin at each seed of TWIN into results/picks/twin;
# then, at the release's seeds, grm learn --fold-file for the default and
# each arm on each dataset with a study plan whose result is missing. After
# each run every result it should have written must exist. A restart on a
# later checkout goes on when the binary is unchanged and the commits differ
# only under scripts/c3-inputs/, docs/, tests/ and .github/; it writes a
# restart line in binary.txt.
#
# OUT gets logs/run.log, logs/<step>.time and .log (each with genome ids
# written as [genome id]; the raw logs stay in work/raw-logs), logs/disk.tsv,
# logs/timing.tsv (per step and fit: wall seconds and peak memory), flags/,
# work/ (release, root with the reports, build, plans, the program's tables,
# the study results: these hold genome ids and stay on the box),
# results/picks/<arm>/ (one result per dataset, split and seed and one inputs
# file per dataset: counts and calls, no genome id), check/ (refit_check.py's
# tables), score/ (score.py's tables), binary.txt, lists.txt and expect.tsv.
# The end checks that no file under results/, check/, score/ or logs/ names a
# genome of the species, then writes SHA256SUMS and OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> results <n> SHA256SUMS <sha256>
# Any failure writes OUT/FAILED naming the step and exits 1. It stops
# nothing but itself. A restart that the kit, list or pin checks refuse
# leaves a finished OUT as it was: DONE stays and no FAILED is written.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?datasets}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}; SEEDS=${6:-42,43,44}; TWIN=${7:-45,46,47}; ARMS=${8:-program,k31-program}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/check "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/c3-inputs
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; PLANS=$W/plans; PICKS=$OUT/results/picks; STUDY=$W/study; AMR=$W/amr
TABLES=$AMR/tool/$SP; REPORTS=$ROOT/species/$SP/results
export RUSTUP_TOOLCHAIN=1.98.1 CARGO_BUILD_JOBS=$JOBS RAYON_NUM_THREADS=$JOBS
export PATH=$HOME/.cargo/bin:/usr/local/cargo/bin:$PATH
# One start at a time on an OUT: a second one would share its steps and files.
exec 9> "$OUT/run.lock"
flock -n 9 || { echo "$(date -u +%FT%TZ) another start runs on this OUT; this one ends"; exit 1; }
STEP=start
say() { echo "$(date -u +%FT%TZ) $*"; }
disk() { echo -e "$(date -u +%FT%TZ)\t$1\t$(du -sb "$OUT" | cut -f1)\t$(df -B1 --output=avail "$OUT" | tail -1 | tr -d ' ')" >> "$OUT/logs/disk.tsv"; }
fail() {  # fail <how>: at <step> or stopped at <step>
  if [ -e "$OUT/DONE" ]; then say "refused $1 $STEP; the finished run's DONE is kept"
  else say "FAILED $1 $STEP"; echo "FAILED $SP $1 $STEP $(date -u +%FT%TZ)" > "$OUT/FAILED"; fi
  exit 1
}
trap 'fail at' ERR
trap 'fail "stopped at"' TERM INT
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
# Genome ids stay on the box. grm and the program's runner name a genome in
# some log lines, and logs/ is copied back. So each step's and fit's log
# moves to work/raw-logs, and logs/ gets it back with every genome id of the
# species' list written as [genome id].
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
export KIT SP DATASETS OUT JOBS SEEDS TWIN ARMS REPO HERE W ROOT R M GRM PLANS PICKS STUDY AMR TABLES REPORTS

STEP=kit
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
[ "$(ulimit -n)" = unlimited ] || [ "$(ulimit -n)" -ge 4096 ] \
  || { say "the open-file limit is $(ulimit -n), below the 4,096 grm matrix needs here"; false; }
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the binary must be its commit"; false; }
# Every tool a later step needs, refused here rather than after hours of fits:
# the scorer imports numpy, the build needs cargo and objcopy, the matrix
# steps unshare, the program's step curl, tar, xz and zstd.
python3 -I -c 'import numpy' 2> /dev/null || { say "python3 has no numpy, which score.py needs; install python3-numpy"; false; }
for tool in cargo objcopy unshare flock /usr/bin/time curl tar xz zstd sha256sum; do
  command -v "$tool" > /dev/null || { say "$tool is not installed, which a later step needs"; false; }
done
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/curated.tar" d77e2c6f75f8d986268d0c075cb2ea364fba9972fa22a55139495ac0eb376d97
pinned "$KIT/curated.sha256" 8b0e43dd98ca1b3fc7609b7e25b1b52d4bb4dd0c4c555cf5e771622d24fd9c6e
pinned "$KIT/fold-plans.tar.xz" b13cad7dbf6318fd175a1d4a1fb1715e7763f19030b48338cbd2cd3f14ca71ec
pinned "$KIT/study-plans.sha256" c9c902f2c9e9b808b8867aff129af1ea736e52dedc93ae1d9e3e631cf39654a1
pinned "$KIT/marks.tsv" ba662467864edce9e8cc617651b81598feac5f32bb1219d7409467ff094aa313
pinned "$KIT/amrkit.tar" da07f277a7d271f1ae006919d2c2df4185737526348302f665a9bcb3a083d0e6
pinned "$KIT/tool-tables.sha256" 9291b83f4c7d01a0384f35446e6de8f986dff5f659439ee3d673bccdd09427ba
for d in $(echo "$DATASETS" | tr ',' '\n'); do
  case $d in "$SP"__*) ;; *) say "$d is not a dataset of $SP"; false ;; esac
done
[ "$(echo "$DATASETS" | tr ',' '\n' | sort | uniq -d | wc -l)" = 0 ] || { say "a dataset is listed twice"; false; }
# Only datasets that are not M are fitted (docs/assumption-plan.md, Rules, Labels).
python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    marks = {r["dataset"]: int(r["pairs_S+I"]) > 0 or int(r["pairs_I+R"]) > 0 for r in csv.DictReader(f, delimiter="\t")}
listed = sys.argv[2].split(",")
gone = [d for d in listed if d not in marks]
m = [d for d in listed if marks.get(d)]
if gone or m:
    print(f"{len(gone)} listed datasets are not in marks.tsv, {len(m)} are marked M: " + ",".join(gone + m))
    sys.exit(1)
' "$KIT/marks.tsv" "$DATASETS" || { say "the list must hold datasets that are not M only"; false; }
echo "$SEEDS" | grep -Eq '^[0-9]+(,[0-9]+)*$' || { say "SEEDS is a comma list of seeds"; false; }
[ "$TWIN" = - ] || echo "$TWIN" | grep -Eq '^[0-9]+(,[0-9]+)*$' || { say "TWIN is a comma list of seeds or -"; false; }
echo ",$SEEDS," | grep -q ',42,' || { say "SEEDS must hold 42, the release's seed"; false; }
echo "$ARMS" | grep -Eq '^(program|k31-program)(,(program|k31-program))*$' \
  || { say "ARMS is a comma list of the built arms, program and k31-program; the others have no fit code yet"; false; }
# The program arms run on the listed datasets in A.
A_LIST=$(python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    a = {r["dataset"] for r in csv.DictReader(f, delimiter="\t") if r["in_A"] == "1"}
print(",".join(d for d in sys.argv[2].split(",") if d in a))
' "$HERE/lists.tsv" "$DATASETS")
export A_LIST
printf 'species\t%s\ndatasets\t%s\na\t%s\nseeds\t%s\ntwin\t%s\narms\t%s\n' "$SP" "$DATASETS" "$A_LIST" "$SEEDS" "$TWIN" "$ARMS" > "$OUT/lists.txt.part"
if [ -s "$OUT/lists.txt" ] && ! cmp -s "$OUT/lists.txt" "$OUT/lists.txt.part"; then
  say "the dataset list, the seeds or the arms differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/lists.txt.part" "$OUT/lists.txt"
rm -f "$OUT/FAILED" "$OUT/DONE"
say "datasets $(echo "$DATASETS" | tr ',' '\n' | wc -l), of them in A $(echo "$A_LIST" | tr ',' '\n' | grep -c . || true)"

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
[ "$(sed -n 's/^grm //p' "$OUT/binary.txt")" = "$(sha256sum "$GRM" | cut -c1-64)" ] \
  || { say "the binary changed since the build"; false; }
built=$(sed -n 's/^commit //p' "$OUT/binary.txt"); head=$(git -C "$REPO" rev-parse HEAD)
if [ "$built" != "$head" ]; then
  git -C "$REPO" cat-file -e "$built^{commit}" 2> /dev/null \
    || { say "the checkout moved since the build and the built commit is not in it; use a fresh OUT"; false; }
  moved=$(git -C "$REPO" diff --name-only "$built" HEAD -- . ':(exclude)scripts/c3-inputs/' ':(exclude)docs/' \
    ':(exclude)tests/' ':(exclude).github/')
  [ -z "$moved" ] || { say "the checkout moved since the build, $(echo "$moved" | wc -l) changed paths outside scripts/c3-inputs/, docs/, tests/ and .github/, the first $(echo "$moved" | head -1); use a fresh OUT"; false; }
  echo "restart $head $(date -u +%FT%TZ)" >> "$OUT/binary.txt"
  say "restart on $head, built at $built: $(git -C "$REPO" diff --name-only "$built" HEAD | wc -l) changed paths, all under scripts/c3-inputs/, docs/, tests/ or .github/"
fi

curated() {
  mkdir -p "$ROOT/curated" "$ROOT/bvbrc"
  tar -xf "$KIT/curated.tar" -C "$ROOT/curated"
  (cd "$ROOT/curated" && sha256sum --quiet -c "$KIT/curated.sha256")
  # The listed datasets take the release's breakpoints, so the seed-42 check
  # compares like with like (scripts/c9-pick/run_species.sh, curated).
  python3 -I -B "$REPO/scripts/c6-mic/release_specs.py" "$ROOT/curated/datasets.json" "$R" "$DATASETS" "$W/specs"
  cp "$W/specs/datasets.json" "$ROOT/curated/datasets.json"
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
  "$GRM" --data "$ROOT" fetch fasta --species "$SP" --workers 8
  n=0; bad=0
  while read -r id; do
    n=$((n + 1)); master "fasta/$id.fna.gz" "$ROOT/fasta/$id.fna.gz" || bad=$((bad + 1))
  done < "$ROOT/species/$SP/genomes.txt"
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
  rm -rf "$PLANS"; mkdir -p "$PLANS/reviewed/folds"
  names=()
  for d in $(echo "$DATASETS" | tr ',' '\n'); do
    if grep -q "/folds/$d.study-whole.json\$\|  folds/$d.study-whole.json\$" "$KIT/study-plans.sha256"; then
      names+=("folds/$d.study-whole.json")
    fi
  done
  if [ "${#names[@]}" = 0 ]; then echo "study plans: none of the datasets has one"; : > "$PLANS/study.txt"; return 0; fi
  tar -xJf "$KIT/fold-plans.tar.xz" -C "$PLANS/reviewed" "${names[@]}"
  for f in "${names[@]}"; do
    want=$(grep -m1 "/$f\$\|  $f\$" "$KIT/study-plans.sha256" | cut -c1-64)
    [ -n "$want" ] && [ "$(sha256sum "$PLANS/reviewed/$f" | cut -c1-64)" = "$want" ] \
      || { echo "$f differs from the reviewed plan's sha256"; exit 1; }
  done
  printf '%s\n' "${names[@]}" | sed 's|^folds/||; s|\.study-whole\.json$||' | LC_ALL=C sort > "$PLANS/study.txt"
  echo "study plans: ${#names[@]}, each equal to its sha256"
}
step plans plans

# The program, from the AMRFinderPlus kit, as scripts/refit/run_species.sh
# installs and proves it.
program_env() {
  mkdir -p "$AMR" && cd "$AMR"
  rm -rf kit && mkdir kit && tar -xf "$KIT/amrkit.tar" -C kit
  curl -fsSL --retry 3 -A curl -o micromamba \
    "https://github.com/mamba-org/micromamba-releases/releases/download/2.10.0-0/micromamba-linux-64"
  pinned micromamba 701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd
  chmod +x micromamba
  rm -rf "$AMR/env"
  MAMBA_ROOT_PREFIX=$AMR/mamba ./micromamba create -y -p "$AMR/env" --file kit/spec/env-amrfinder.explicit.txt
  export PATH="$AMR/env/bin:$PATH"
  amrfinder --version | grep -qx '4.2.7'
  cd kit
  rm -rf "$AMR/db/2026-08-07.1"
  bash scripts/fetch_database.sh 2026-08-07.1 spec/db-2026-08-07.1.sha256 "$AMR/db/2026-08-07.1"
  python3 -I -B tests/test_contract.py
  python3 -I -B tests/test_runner.py
  mkdir -p "$AMR/testdna"
  curl -fsSL --retry 3 -A curl -o "$AMR/testdna/test_dna.fa" \
    "https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_dna.fa"
  pinned "$AMR/testdna/test_dna.fa" 2c330ed5f979e18ba915ed3aea36c3a0779829e8421d1e1c57c6df842aa2b394
  rm -rf "$AMR/fx" "$AMR/fx-out"
  python3 -I -B scripts/make_fixture.py --test-dna "$AMR/testdna/test_dna.fa" --out "$AMR/fx"
  bash scripts/run_amrfinder.sh "$AMR/fx/cohort.tsv" "$AMR/fx/genomes" "$AMR/fx-out" "$AMR/fx/run-manifest.json" \
    "$AMR/db/2026-08-07.1" 2 1
  diff <(cd "$AMR/fx-out" && sha256sum */*.tsv) <(cd fixtures/e2e/tool-out && sha256sum */*.tsv)
  echo "install proved: the fixture's tool tables equal the kit's record"
}
# The program on every genome of the species with a table the release
# recorded; each table must equal its recorded sha256, so the calls are the
# release's (the cohort comes from tool-tables.sha256's rows).
program_run() {
  export PATH="$AMR/env/bin:$PATH"; cd "$AMR/kit"
  grep "  $SP/" "$KIT/tool-tables.sha256" | sort -k2 > "$AMR/tables.recorded"
  sed -E "s#^[0-9a-f]{64}  $SP/(.+)\.tsv\$#\1#" "$AMR/tables.recorded" | LC_ALL=C sort > "$AMR/cohort.ids"
  outside=$(LC_ALL=C comm -23 "$AMR/cohort.ids" <(LC_ALL=C sort "$ROOT/species/$SP/genomes.txt") | wc -l)
  [ "$outside" -eq 0 ] || { echo "recorded tables: $outside are not of a genome in the species list"; exit 1; }
  echo "cohort: $(wc -l < "$AMR/cohort.ids") genomes with a recorded table, of $(wc -l < "$ROOT/species/$SP/genomes.txt") in the species list"
  { printf 'species_slug\tgenome_id\n'; sed "s/^/$SP\t/" "$AMR/cohort.ids"; } > "$AMR/cohort.tsv"
  rm -rf "$AMR/tool"
  bash scripts/run_amrfinder.sh "$AMR/cohort.tsv" "$ROOT/fasta" "$AMR/tool" spec/run-manifest.rev14.json \
    "$AMR/db/2026-08-07.1" "$JOBS" 1 4000
  ( cd "$AMR/tool" && find "$SP" -name '*.tsv' | sort | xargs sha256sum ) | sort -k2 > "$AMR/tables.now"
  echo "tool tables: $(wc -l < "$AMR/tables.now") made, $(wc -l < "$AMR/tables.recorded") recorded"
  diff "$AMR/tables.recorded" "$AMR/tables.now"
  echo "every tool table equals its recorded sha256"
}
if [ -n "$A_LIST" ]; then
  step program_env program_env
  step program_run program_run
else
  say "program: no listed dataset is in A, so no program arm runs"
fi

STEP=expect
python3 -I -B "$HERE/expect.py" "$SP" "$DATASETS" "${A_LIST:--}" "$(paste -sd, "$PLANS/study.txt")" "$SEEDS" "$TWIN" "$ARMS" \
  > "$OUT/expect.tsv.part"
if [ -s "$OUT/expect.tsv" ] && ! cmp -s "$OUT/expect.tsv" "$OUT/expect.tsv.part"; then
  say "the expected results differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"
say "expect: $(($(wc -l < "$OUT/expect.tsv") - 1)) results"

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
# The datasets with a missing result for an arm, split and seed, as a comma
# list. An arm's result is missing until its pick file and its inputs file
# both exist.
missing() { awk -F'\t' -v a="$1" -v p="$2" -v s="$3" 'NR > 1 && $3 == a && $4 == p && $6 == s { print $1 "\t" $8 "\t" $9 }' "$OUT/expect.tsv" \
  | while IFS=$'\t' read -r d file inputs; do
      { [ -e "$PICKS/$file" ] && { [ "$inputs" = - ] || [ -e "$PICKS/$inputs" ]; }; } || echo "$d"
    done | sort -u | paste -sd, -; }
check42() {
  STEP=check
  say "check: start"
  python3 -I -B "$REPO/scripts/refit/refit_check.py" "$REPORTS" "$R" "$DATASETS" "$OUT/check" > "$OUT/logs/check.log" 2>&1 \
    || { cat "$OUT/logs/check.log"; false; }
  cat "$OUT/logs/check.log"
}
arm_args() {  # arm_args <arm>: the flags of an arm's fit, beside --exp-pick
  [ "$1" = k31 ] || [ "$1" = twin ] || echo "--exp-inputs $1 --exp-program-tables $TABLES"
}
lineage_fit() {  # lineage_fit <arm> <seed>
  local arm=$1 seed=$2 todo left
  STEP=fit-$arm-s$seed
  todo=$(missing "$arm" lineage "$seed")
  if [ -z "$todo" ]; then say "fit-$arm-s$seed: all results present, skipped"; return 0; fi
  mkdir -p "$PICKS/$arm"
  # shellcheck disable=SC2046
  fit "fit-$arm-s$seed" --only "$todo" --force --seed "$seed" --exp-pick "$PICKS/$arm" --exp-pick-inner 0 $(arm_args "$arm")
  left=$(missing "$arm" lineage "$seed")
  [ -z "$left" ] || { say "fit-$arm-s$seed: no result written for $left"; false; }
  disk "fit-$arm-s$seed"
}

for s in $(echo "$SEEDS" | tr ',' '\n'); do
  lineage_fit k31 "$s"
  # The gate: no arm runs unless the default's seed-42 reports are the release's.
  [ "$s" != 42 ] || check42
  for arm in $(echo "$ARMS" | tr ',' '\n'); do
    [ -z "$A_LIST" ] || lineage_fit "$arm" "$s"
  done
done
[ "$TWIN" = - ] || for s in $(echo "$TWIN" | tr ',' '\n'); do lineage_fit twin "$s"; done
study_have() {  # study_have <arm> <dataset>: its study result, and an arm's inputs file, exist
  [ -e "$PICKS/$1/$2.study-whole-s42.pick.json" ] && { [ "$1" = k31 ] || [ -e "$PICKS/$1/$2.study-whole.inputs.json" ]; }
}
while read -r d; do
  [ -n "$d" ] || continue
  for arm in k31 $(echo "$ARMS" | tr ',' '\n'); do
    case $arm in k31) ;; *) echo ",$A_LIST," | grep -q ",$d," || continue ;; esac
    STEP=study-$arm-$d
    if study_have "$arm" "$d"; then say "study-$arm-$d: result present, skipped"; continue; fi
    mkdir -p "$PICKS/$arm" "$STUDY/$arm"
    # shellcheck disable=SC2046
    fit "study-$arm-$d" --fold-file "$PLANS/reviewed/folds/$d.study-whole.json" --fold-out "$STUDY/$arm" \
      --exp-pick "$PICKS/$arm" --exp-pick-inner 0 --force $(arm_args "$arm")
    study_have "$arm" "$d" || { say "study-$arm-$d: no result written"; false; }
    disk "study-$arm-$d"
  done
done < "$PLANS/study.txt"

STEP=timing
{ printf 'step\twall_s\tpeak_kb\n'
  for f in "$OUT"/logs/*.time; do
    name=$(basename "$f" .time)
    if grep -q 'Elapsed (wall clock)' "$f"; then
      wall=$(sed -n 's/.*Elapsed (wall clock) time (h:mm:ss or m:ss): //p' "$f" \
        | awk -F: '{ s = 0; for (i = 1; i <= NF; i++) s = s * 60 + $i; printf "%.2f", s }')
      peak=$(sed -n 's/.*Maximum resident set size (kbytes): //p' "$f")
    else
      wall=$(sed -E 's/^.* ([0-9.]+) s wall.*$/\1/' "$f"); peak=$(sed -E 's/^.* ([0-9]+) KB peak.*$/\1/' "$f")
    fi
    printf '%s\t%s\t%s\n' "$name" "$wall" "$peak"
  done; } > "$OUT/logs/timing.tsv"

STEP=score
say "score: start"
# score.py prints the check counts only; its tables stay in score/.
python3 -I -B "$HERE/score.py" "$OUT/expect.tsv" "$PICKS" "$R" "$OUT/score" > "$OUT/logs/score.log" 2>&1 \
  || { scrub "$OUT/logs/score.log"; cat "$OUT/logs/score.log"; false; }
scrub "$OUT/logs/score.log"
cat "$OUT/logs/score.log"

STEP=done
cd "$OUT"
python3 -I -c '
import json, os, re, sys
with open(sys.argv[1]) as f:
    ids = {g["id"] for g in json.load(f)["genomes"]}
token = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")
hit = 0
for top in sys.argv[2:]:
    for d, _, files in os.walk(top):
        for name in files:
            with open(os.path.join(d, name), errors="replace") as f:
                hit += sum(1 for m in token.finditer(f.read()) if m.group(0) in ids)
print("genome ids in results, check, score and logs:", hit)
sys.exit(hit != 0)
' "$R/species/$SP.json" results check score logs
{ find results check score -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.txt' \); echo binary.txt; echo expect.tsv; echo lists.txt; echo logs/timing.tsv; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) results $(find results/picks -name '*.pick.json' | wc -l) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
