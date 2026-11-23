#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARMS] [TABLES] [GRAPH]
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
# default fitted again. ARMS (default all five: program, k31-program,
# k31-counts, k31-aa and unitigs) are the arms fitted beside the default;
# k31-counts, k31-aa and unitigs run on every listed dataset. TABLES (default none) is a
# folder of the program's tables made before, as TABLES/SPECIES/<genome>.tsv:
# when it is given the program is neither installed nor run, and each table
# of the cohort is copied from it and held to its recorded sha256. Every fit
# runs with --exp-pick-inner 0: each outer fold's calls per learner, no
# inner pick. GRAPH (default grm) is who builds the unitig graph: grm, GRM's
# own compaction, which the lead's ruling allows only once its unitigs equal
# an outside build's as a set on a pilot species (they did on Klebsiella
# aerogenes, 93 genomes, and on Escherichia coli, 2,000: docs/input-arms.md);
# bcalm2, the outside program, whose unitigs GRM reads; or check, both, with
# GRM's graph used and the two compared, the equality counts written in the
# build's record, and any unitig or 31-mer in one build only stopping the
# species.
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
#   prodigal.linux              with k31-aa: Prodigal 2.6.3, the v2.6.3
#                               release's Linux binary (hyattpd/Prodigal)
#   bcalm-binaries-v2.2.3-Linux.tar.gz
#                               with unitigs and GRAPH bcalm2 or check:
#                               bcalm2 2.2.3's release binaries (GATB/bcalm)
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
# species with a recorded table, each table equal to its recorded sha256,
# or with TABLES, those tables copied in and held to the same sha256 (a
# table missing or unequal stops the run with the counts); expect.tsv
# (expect.py). Every start, a restart too, holds the tables to
# their recorded sha256 again before any fit.
#
# The matrix arms' builds, after matrix, each once, timed and reading no
# label, phenotype, fold or seed (docs/input-arms.md): tools (Prodigal and
# bcalm2 from the kit, their version and sha256 written in binary.txt);
# inputs_counts (grm exp-input-build counts); prodigal (Prodigal 2.6.3 on
# every genome, -p single -g 11, a separate program) and inputs_aa (grm
# exp-input-build aa, which holds every protein to GRM's own table 11
# translation); bcalm (bcalm2 over every genome, k 31, abundance 1) and
# inputs_unitigs (grm exp-input-build unitigs). Each writes its matrix under
# work/inputs/<kind>, the sha256 of its files in logs/inputs-<kind>.sha256
# and its counts-only record in results/builds/<kind>.build.json, and
# removes what it no longer needs (the proteins, bcalm2's unitigs). Every
# start holds each built matrix to its recorded sha256 again before any fit.
# A build the kernel kills for want of memory fails its step as any failure
# does: OUT/FAILED names it, nothing after it runs, and a restart with a
# smaller JOBS (fewer threads, so fewer genomes and buckets held at once; 32 GB
# is the largest box) runs that step again and skips the finished ones.
#
# Fits, per seed of SEEDS: the default, grm learn --seed s --exp-pick
# results/picks/k31, then each arm, grm learn --seed s --exp-inputs ARM,
# with --exp-program-tables <the tables> for the program arms and
# --exp-input-matrix work/inputs/<kind> for the others, and --exp-pick
# results/picks/ARM, each over
# the datasets whose result at that seed is missing (--only, --force). After
# the default's seed-42 run, scripts/refit/refit_check.py compares each
# seed-42 report with the release's; any difference stops the run before
# another fit. Then the twin at each seed of TWIN into results/picks/twin;
# then, at the release's seeds, grm learn --fold-file for the default and
# each arm on each dataset with a study plan whose result is missing. After
# each run every result it should have written must exist. When grm learn
# refuses an arm's study fit because an outer fold's training genomes give
# no candidate feature, and the default's study result at that dataset
# exists, refused.py writes results/picks/ARM/<dataset>.study-whole-s42.refused.json
# (the reason and the fold's counts) in place of the result and the run goes
# on: score.py scores it as the constant model on the default's folds (the
# lead's ruling 16 of 23 Nov 2026). Any other failure, and that refusal of
# the default, stops the run. The drawn splits never refuse for it: a fold
# with no candidate gets the release's fits on no feature. A restart on a
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
# tables), score/ (score.py's tables and its output, score.log, kept out of
# logs/ and run.log), binary.txt, lists.txt and expect.tsv.
# The end checks that no file under results/, check/, score/ or logs/ names a
# genome of the species, then writes SHA256SUMS and OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> results <n> SHA256SUMS <sha256>
# Any failure writes OUT/FAILED naming the step and exits 1. It stops
# nothing but itself. A restart that the kit, list or pin checks refuse
# leaves a finished OUT as it was: DONE stays and no FAILED is written.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?datasets}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}; SEEDS=${6:-42,43,44}; TWIN=${7:-45,46,47}; ARMS=${8:-program,k31-program,k31-counts,k31-aa,unitigs}
TABLES_FROM=${9:-}; GRAPH=${10:-grm}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/check "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/c3-inputs
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; PLANS=$W/plans; PICKS=$OUT/results/picks; STUDY=$W/study; AMR=$W/amr
TABLES=$AMR/tool/$SP; REPORTS=$ROOT/species/$SP/results
INPUTS=$W/inputs; PRODIGAL=$W/tools/prodigal; BCALM=$W/tools/bcalm-binaries-v2.2.3-Linux/bin/bcalm
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
# keep_build <kind>: the sha256 of the built matrix, held again at every
# start, and the build's counts-only record among the results. A step that
# builds one runs it in its own shell, so step() hands it over.
keep_build() {
  ( cd "$INPUTS" && find "$1" -type f | LC_ALL=C sort | xargs sha256sum ) > "$OUT/logs/inputs-$1.sha256"
  mkdir -p "$OUT/results/builds"; cp "$INPUTS/$1/build.json" "$OUT/results/builds/$1.build.json"
  python3 -I -c '
import json, sys
b = json.load(open(sys.argv[1]))
print(" ".join(f"{k} {v}" for k, v in b.items() if isinstance(v, int) and not isinstance(v, bool)))
for k, v in (b.get("equality") or {}).items():
    print(f"equality {k} {v}")
' "$INPUTS/$1/build.json"
}
# killed_for_memory <time file>: whether GNU time's record says the kernel's
# SIGKILL ended the program, as it ends one when memory runs out. A program
# time runs itself is "terminated by signal 9"; a step runs in bash, whose exit
# is then 137, which time writes as "non-zero status 137".
killed_for_memory() {
  grep -qE 'terminated by signal 9|non-zero status 137' "$1"
}
step() {  # step <name> <function>: once, timed, in a shell that stops at the first failure
  STEP=$1
  if [ -e "$OUT/flags/$1.done" ]; then say "$1: done before, skipped"; return 0; fi
  say "$1: start"
  /usr/bin/time -f "$1 %e s wall, %M KB peak, %U s user, %S s sys" -o "$OUT/logs/$1.time" \
    bash -euo pipefail -c "$(declare -f say pinned master keep_build "$2"); $2" < /dev/null > "$OUT/logs/$1.log" 2>&1 \
    || { scrub "$OUT/logs/$1.log"; tail -n 30 "$OUT/logs/$1.log"
         ! killed_for_memory "$OUT/logs/$1.time" \
           || say "$1: killed by signal 9, as the kernel kills a program when memory runs out; restart with a smaller JOBS, which holds fewer genomes and buckets at once"
         false; }
  scrub "$OUT/logs/$1.log"
  touch "$OUT/flags/$1.done"
  disk "$1"
  say "$1: end, $(cat "$OUT/logs/$1.time")"
}
export KIT SP DATASETS OUT JOBS SEEDS TWIN ARMS TABLES_FROM GRAPH REPO HERE W ROOT R M GRM PLANS PICKS STUDY AMR TABLES REPORTS
export INPUTS PRODIGAL BCALM
has_arm() { echo ",$ARMS," | grep -q ",$1,"; }

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
ARM_NAMES='program|k31-program|k31-counts|k31-aa|unitigs'
echo "$ARMS" | grep -Eq "^($ARM_NAMES)(,($ARM_NAMES))*\$" \
  || { say "ARMS is a comma list of the arms: program, k31-program, k31-counts, k31-aa, unitigs"; false; }
case $GRAPH in grm|bcalm2|check) ;; *) say "GRAPH is grm, bcalm2 or check"; false ;; esac
# The outside programs of the matrix arms, refused here when the kit lacks them.
if has_arm k31-aa; then
  [ -f "$KIT/prodigal.linux" ] || { say "k31-aa needs Prodigal 2.6.3, the kit's prodigal.linux, which is absent"; false; }
  pinned "$KIT/prodigal.linux" 80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246
fi
if has_arm unitigs && [ "$GRAPH" != grm ]; then
  [ -f "$KIT/bcalm-binaries-v2.2.3-Linux.tar.gz" ] \
    || { say "unitigs with GRAPH $GRAPH needs bcalm2 2.2.3, the kit's bcalm-binaries-v2.2.3-Linux.tar.gz, which is absent"; false; }
  pinned "$KIT/bcalm-binaries-v2.2.3-Linux.tar.gz" 20479729a9447980a53164f9b62bc04ce8e08246fffb24af9d49e9d07858e2ff
fi
if [ -n "$TABLES_FROM" ]; then
  [ -d "$TABLES_FROM/$SP" ] || { say "TABLES has no folder $SP"; false; }
  TABLES_FROM=$(cd "$TABLES_FROM" && pwd)
fi
# The program arms run on the listed datasets in A.
A_LIST=$(python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    a = {r["dataset"] for r in csv.DictReader(f, delimiter="\t") if r["in_A"] == "1"}
print(",".join(d for d in sys.argv[2].split(",") if d in a))
' "$HERE/lists.tsv" "$DATASETS")
export A_LIST
printf 'species\t%s\ndatasets\t%s\na\t%s\nseeds\t%s\ntwin\t%s\narms\t%s\ngraph\t%s\n' "$SP" "$DATASETS" "$A_LIST" "$SEEDS" "$TWIN" "$ARMS" "$GRAPH" > "$OUT/lists.txt.part"
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

# The matrix arms' label-free builds, once per species over all its genomes.
# Each reads the genomes, the species' ksets or an outside program's output;
# none reads a label, a phenotype, a fold or a seed.
tools() {
  rm -rf "$W/tools"; mkdir -p "$W/tools"
  case ",$ARMS," in *,k31-aa,*)
    install -m 755 "$KIT/prodigal.linux" "$PRODIGAL"
    "$PRODIGAL" -v 2>&1 | grep -q 'Prodigal V2\.6\.3' || { echo "the kit's Prodigal is not 2.6.3"; exit 1; }
    echo "prodigal 2.6.3 $(sha256sum "$PRODIGAL" | cut -c1-64)" >> "$OUT/binary.txt" ;;
  esac
  case ",$ARMS,:$GRAPH" in *,unitigs,*:bcalm2 | *,unitigs,*:check)
    tar -xzf "$KIT/bcalm-binaries-v2.2.3-Linux.tar.gz" -C "$W/tools"
    "$BCALM" -version 2>&1 | grep -q 'version v2\.2\.3' || { echo "the kit's bcalm is not 2.2.3"; exit 1; }
    echo "bcalm 2.2.3 $(sha256sum "$BCALM" | cut -c1-64)" >> "$OUT/binary.txt" ;;
  esac
  grep -E '^(prodigal|bcalm) ' "$OUT/binary.txt" || echo "no outside program needed"
}
inputs_counts() {
  rm -rf "$INPUTS/counts"
  unshare -rn "$GRM" --data "$ROOT" exp-input-build counts --species "$SP" --out "$INPUTS/counts"
  keep_build counts
}
prodigal() {
  rm -rf "$W/prodigal"; mkdir -p "$W/prodigal"
  one() {
    zcat "$ROOT/fasta/$1.fna.gz" > "$W/prodigal/$1.fna"
    "$PRODIGAL" -i "$W/prodigal/$1.fna" -a "$W/prodigal/$1.faa.part" -p single -g 11 -q -o /dev/null
    mv "$W/prodigal/$1.faa.part" "$W/prodigal/$1.faa"; rm -f "$W/prodigal/$1.fna"
  }
  export -f one
  xargs -P "$JOBS" -I{} bash -euo pipefail -c 'one "$1"' _ {} < "$ROOT/species/$SP/genomes.txt"
  n=$(find "$W/prodigal" -name '*.faa' | wc -l); want=$(wc -l < "$ROOT/species/$SP/genomes.txt")
  echo "proteins: Prodigal's file for $n of $want genomes, $(du -sb "$W/prodigal" | cut -f1) bytes"
  [ "$n" -eq "$want" ]
}
inputs_aa() {
  rm -rf "$INPUTS/aa"
  unshare -rn "$GRM" --data "$ROOT" exp-input-build aa --species "$SP" --proteins "$W/prodigal" --out "$INPUTS/aa"
  keep_build aa
  rm -rf "$W/prodigal"
}
bcalm() {
  rm -rf "$W/bcalm"; mkdir -p "$W/bcalm/tmp"; cd "$W/bcalm"
  sed "s#^#$ROOT/fasta/#; s#\$#.fna.gz#" "$ROOT/species/$SP/genomes.txt" > list.txt
  mem=$(awk '/^MemTotal/ { print int($2 / 2048) }' /proc/meminfo)
  "$BCALM" -in list.txt -kmer-size 31 -abundance-min 1 -nb-cores "$JOBS" -max-memory "$mem" -out graph -out-tmp tmp
  rm -rf tmp graph.h5 graph.unitigs.fa.glue* list.txt
  echo "bcalm2: $(grep -c '^>' graph.unitigs.fa) unitigs, $(stat -c %s graph.unitigs.fa) bytes"
}
inputs_unitigs() {
  rm -rf "$INPUTS/unitigs"
  case $GRAPH in
    grm) unshare -rn "$GRM" --data "$ROOT" exp-input-build unitigs --species "$SP" --graph grm --out "$INPUTS/unitigs" ;;
    bcalm2) unshare -rn "$GRM" --data "$ROOT" exp-input-build unitigs --species "$SP" --graph outside \
              --outside-unitigs "$W/bcalm/graph.unitigs.fa" --out "$INPUTS/unitigs" ;;
    check) unshare -rn "$GRM" --data "$ROOT" exp-input-build unitigs --species "$SP" --graph grm \
             --outside-unitigs "$W/bcalm/graph.unitigs.fa" --out "$INPUTS/unitigs" ;;
  esac
  keep_build unitigs
  # Check mode: by ruling 1 GRM's graph stands only while its unitigs equal
  # bcalm2's as a set, so any unitig or 31-mer in one build only stops the
  # species here; it then runs again in a fresh OUT with GRAPH bcalm2.
  if [ "$GRAPH" = check ]; then
    python3 -I -c 'import json, sys
e = json.load(open(sys.argv[1]))["equality"]
print("graph check: " + ", ".join(f"{k} {v}" for k, v in sorted(e.items())))
sys.exit(1 if any(v for k, v in e.items() if k.endswith("_only")) else 0)' "$INPUTS/unitigs/build.json" \
      || { say "GRM's unitigs differ from bcalm2's as a set: by ruling 1 the species runs again in a fresh OUT with GRAPH bcalm2"; false; }
  fi
  rm -rf "$W/bcalm"
}
if has_arm k31-aa || { has_arm unitigs && [ "$GRAPH" != grm ]; }; then step tools tools; fi
if has_arm k31-counts; then step inputs_counts inputs_counts; fi
if has_arm k31-aa; then step prodigal prodigal; step inputs_aa inputs_aa; fi
if has_arm unitigs; then
  [ "$GRAPH" = grm ] || step bcalm bcalm
  step inputs_unitigs inputs_unitigs
fi
# Every start holds each built matrix to the sha256 its step recorded: a
# restart skips the step that built it, and the fits read it.
for kind in counts aa unitigs; do
  [ -e "$OUT/flags/inputs_$kind.done" ] || continue
  STEP=inputs-check-$kind
  ( cd "$INPUTS" && sha256sum --quiet -c "$OUT/logs/inputs-$kind.sha256" ) \
    || { say "the $kind matrix differs from the sha256 its build recorded"; false; }
  say "inputs: the $kind matrix equals its recorded sha256, $(wc -l < "$OUT/logs/inputs-$kind.sha256") files"
done

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
# With TABLES: the same cohort's tables, made before, copied in and held to
# the same recorded sha256; no program is installed or run.
program_tables() {
  mkdir -p "$AMR"
  grep "  $SP/" "$KIT/tool-tables.sha256" | sort -k2 > "$AMR/tables.recorded"
  sed -E "s#^[0-9a-f]{64}  $SP/(.+)\.tsv\$#\1#" "$AMR/tables.recorded" | LC_ALL=C sort > "$AMR/cohort.ids"
  outside=$(LC_ALL=C comm -23 "$AMR/cohort.ids" <(LC_ALL=C sort "$ROOT/species/$SP/genomes.txt") | wc -l)
  [ "$outside" -eq 0 ] || { echo "recorded tables: $outside are not of a genome in the species list"; exit 1; }
  echo "cohort: $(wc -l < "$AMR/cohort.ids") genomes with a recorded table, of $(wc -l < "$ROOT/species/$SP/genomes.txt") in the species list"
  rm -rf "$AMR/tool"; mkdir -p "$TABLES"
  missing=0
  while read -r id; do
    if [ -f "$TABLES_FROM/$SP/$id.tsv" ]; then cp "$TABLES_FROM/$SP/$id.tsv" "$TABLES/$id.tsv"; else missing=$((missing + 1)); fi
  done < "$AMR/cohort.ids"
  ( cd "$AMR/tool" && find "$SP" -name '*.tsv' | sort | xargs -r sha256sum ) | sort -k2 > "$AMR/tables.now"
  unequal=$(LC_ALL=C join -1 2 -2 2 <(LC_ALL=C sort -k2,2 "$AMR/tables.recorded") <(LC_ALL=C sort -k2,2 "$AMR/tables.now") \
    | awk '$2 != $3' | wc -l)
  echo "tool tables: $(wc -l < "$AMR/tables.recorded") recorded, $(wc -l < "$AMR/tables.now") copied from TABLES, $missing missing there, $unequal unequal to their recorded sha256; $(find "$TABLES_FROM/$SP" -name '*.tsv' | wc -l) in TABLES"
  if [ "$missing" -ne 0 ] || [ "$unequal" -ne 0 ] || ! cmp -s "$AMR/tables.recorded" "$AMR/tables.now"; then
    echo "the tables in TABLES are not the recorded ones"; exit 1
  fi
  echo "every tool table equals its recorded sha256; the program was not run"
}
PROGRAM_RUNS=0
if [ -n "$A_LIST" ] && { has_arm program || has_arm k31-program; }; then PROGRAM_RUNS=1; fi
if [ "$PROGRAM_RUNS" = 0 ]; then
  say "program: no listed dataset is in A or ARMS holds no program arm, so no program arm runs"
elif [ -e "$OUT/flags/program_run.done" ] || { [ -z "$TABLES_FROM" ] && [ ! -e "$OUT/flags/program_tables.done" ]; }; then
  step program_env program_env
  step program_run program_run
else
  [ -z "$TABLES_FROM" ] || say "program: the tables come from TABLES; the program is neither installed nor run"
  step program_tables program_tables
fi
# Every start holds the tables to their recorded sha256 again before any
# fit: a restart skips the step that made them, and the fits read them.
if [ "$PROGRAM_RUNS" = 1 ]; then
  STEP=tables-check
  grep "  $SP/" "$KIT/tool-tables.sha256" | sort -k2 > "$AMR/tables.recorded.again"
  ( cd "$AMR/tool" && find "$SP" -name '*.tsv' | sort | xargs -r sha256sum ) | sort -k2 > "$AMR/tables.again"
  if ! cmp -s "$AMR/tables.recorded.again" "$AMR/tables.again"; then
    say "the program's tables differ from their recorded sha256: $(wc -l < "$AMR/tables.again") on the box, $(wc -l < "$AMR/tables.recorded.again") recorded, $(diff "$AMR/tables.recorded.again" "$AMR/tables.again" | grep -c '^[<>]' || true) lines differ"
    false
  fi
  say "tables: $(wc -l < "$AMR/tables.again") of $(wc -l < "$AMR/tables.recorded.again") recorded, each equal to its sha256"
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
    < /dev/null > "$OUT/logs/$name.log" 2>&1 || { scrub "$OUT/logs/$name.log"; tail -n 30 "$OUT/logs/$name.log"
      ! killed_for_memory "$OUT/logs/$name.time" \
        || say "$name: killed by signal 9, as the kernel kills a program when memory runs out; restart with a smaller JOBS"
      false; }
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
  case $1 in
    k31 | twin) ;;
    program | k31-program) echo "--exp-inputs $1 --exp-program-tables $TABLES" ;;
    k31-counts) echo "--exp-inputs $1 --exp-input-matrix $INPUTS/counts" ;;
    k31-aa) echo "--exp-inputs $1 --exp-input-matrix $INPUTS/aa" ;;
    unitigs) echo "--exp-inputs $1 --exp-input-matrix $INPUTS/unitigs" ;;
  esac
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
  # The program arms run on the A datasets, the others on every listed one.
  for arm in $(echo "$ARMS" | tr ',' '\n'); do
    case $arm in program | k31-program) [ "$PROGRAM_RUNS" = 0 ] || lineage_fit "$arm" "$s" ;; *) lineage_fit "$arm" "$s" ;; esac
  done
done
[ "$TWIN" = - ] || for s in $(echo "$TWIN" | tr ',' '\n'); do lineage_fit twin "$s"; done
study_have() {  # study_have <arm> <dataset>: its study result (or an arm's refusal record) and an arm's inputs file exist
  { [ -e "$PICKS/$1/$2.study-whole-s42.pick.json" ] \
    || { [ "$1" != k31 ] && [ -e "$PICKS/$1/$2.study-whole-s42.refused.json" ]; }; } \
    && { [ "$1" = k31 ] || [ -e "$PICKS/$1/$2.study-whole.inputs.json" ]; }
}
# One study fit, as fit, except that an arm's refusal for no candidate,
# where the default's study result exists, becomes a refusal record and the
# run goes on (ruling 16). The default's refusal stops the run.
study_fit() {  # study_fit <arm> <dataset>
  local arm=$1 d=$2 name=study-$1-$2 rc=0
  STEP=$name
  say "$name: start"
  # shellcheck disable=SC2046
  /usr/bin/time -v -o "$OUT/logs/$name.time" "$GRM" --data "$ROOT" learn --species "$SP" \
    --fold-file "$PLANS/reviewed/folds/$d.study-whole.json" --fold-out "$STUDY/$arm" \
    --exp-pick "$PICKS/$arm" --exp-pick-inner 0 --force $(arm_args "$arm") \
    < /dev/null > "$OUT/logs/$name.log" 2>&1 || rc=$?
  scrub "$OUT/logs/$name.log"
  if [ "$rc" != 0 ]; then
    if [ "$arm" != k31 ] && [ -e "$PICKS/k31/$d.study-whole-s42.pick.json" ] \
      && python3 -I -B "$HERE/refused.py" "$OUT/logs/$name.log" "$d" "$arm" "$PICKS/$arm/$d.study-whole-s42.refused.json"; then
      say "$name: $arm/$d.study-whole-s42.refused.json written, to be scored as the constant model"
    else
      tail -n 30 "$OUT/logs/$name.log"; false
    fi
  fi
  say "$name: end, $(grep -E 'Elapsed|Maximum resident' "$OUT/logs/$name.time" | sed 's/^\s*//' | tr '\n' ';')"
}
while read -r d; do
  [ -n "$d" ] || continue
  for arm in k31 $(echo "$ARMS" | tr ',' '\n'); do
    case $arm in program | k31-program) echo ",$A_LIST," | grep -q ",$d," || continue ;; esac
    STEP=study-$arm-$d
    if study_have "$arm" "$d"; then say "study-$arm-$d: result present, skipped"; continue; fi
    mkdir -p "$PICKS/$arm" "$STUDY/$arm"
    study_fit "$arm" "$d"
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
# score.py prints the check counts, or why it stopped. Its output stays with
# its tables in score/, which no one opens before every species has run: not
# in logs/ and not in run.log, which are read while the run goes on.
python3 -I -B "$HERE/score.py" "$OUT/expect.tsv" "$PICKS" "$R" "$OUT/score" > "$OUT/score/score.log" 2>&1 \
  || { scrub "$OUT/score/score.log"; say "score: stopped; its output is in score/score.log"; false; }
scrub "$OUT/score/score.log"
say "score: end"

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
{ find results check score -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.txt' \); echo binary.txt; echo expect.tsv; echo lists.txt; echo logs/timing.tsv
  find logs -name 'inputs-*.sha256'; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
say "arm results refused for no candidate: $(find results/picks -name '*.refused.json' | wc -l)"
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) results $(find results/picks -name '*.pick.json' | wc -l) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
