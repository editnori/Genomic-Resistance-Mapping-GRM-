#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS] [DRAWS] [TWIN] [ARMS] [SPLITS] [SCORE] [AT_ONCE]
#
# The cap-arm run of one species (plan block C2, docs/cap-arms.md) on a
# fresh box with no agent attached: no prompt, nothing read from another
# machine but the public sources below, one log, and a restart that skips
# every finished step and every finished fit. Run it from a clean checkout
# of this branch:
#
#   setsid nohup bash scripts/c2-cap/run_species.sh KIT SPECIES DATASETS OUT 16 > /dev/null 2>&1 < /dev/null &
#
# DATASETS is a comma list of the species' first-run datasets
# (scripts/c2-cap/lists.tsv: CAP, not M). DRAWS (default 3) are the lineage
# draws, grm arm's keyed slots 0 to DRAWS-1; TWIN (default 3, or 0 for none)
# the twin's slots, the default's next TWIN slots. ARMS (default cap-all)
# are the wider arms; cap-4000 is withdrawn for the run (docs/cap-arms.md)
# and builds no root unless named here. SPLITS (default lineage,study) are
# the splits fitted; the study split runs on the listed datasets with a
# reviewed study-whole plan. SCORE (default yes) runs score.py at the end;
# with no, nothing reads a score: after each dataset's fits their
# predictions.jsonl, runs.json and repeats.json move to work/unread, which
# nothing opens (the pilot). AT_ONCE (default 1, at most 4) is how many
# datasets' lineage fits run at once; see Fits below. KIT holds these
# files, checked against the sha256 pinned here before anything else runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   MASTER-MANIFEST.sha256      the data root's list of files and sha256
#   curated.tar, curated.sha256 the rebuilt revision 14 table's curated files
#   fold-plans.tar.xz, study-plans.sha256
#                               the reviewed study plans per dataset and the
#                               sha256 of each, taken before any fit
#   marks.tsv                   the label count's table: a dataset is M when
#                               pairs_S+I or pairs_I+R is above 0
#   master-genome-meta.jsonl.gz the master root's genome rows, read for the
#                               BioProject of each added genome
#
# Steps, each once: the release; grm built from this checkout, with its
# commit, sha256 and code sha256; the curated table with the listed
# datasets' breakpoints set to the release's and the catalog's three files
# from NCBI, each equal to the master list (shared by the roots); the
# cohorts (cohorts.py: the release's 2,000, and the top 4,000 and 4,096 by
# the release's rule on the rebuilt table, and per dataset and wider arm
# the added genomes with a call, by class); the BioProject of each added
# genome from BV-BRC's genome table by id, two fields (bioprojects.py,
# only when a listed dataset has a study plan); the fetch of every cap-all
# genome from BV-BRC (grm fetch fasta), each release genome held to the
# master list and the release's length, contig count and N50
# (fetch_check.py, which also times each genome and drops a genome BV-BRC no
# longer serves from the wider cohorts); one data root per cohort; count in
# the cap-all root, then each other root takes its ksets as hard links
# (roots.py, the sharing of scripts/input-comparison/run.py) and counts
# what is left, which is nothing; structure and matrix in each root, with
# no network; the study plans (the reviewed study-whole plan of each
# dataset that has one, held to its sha256) and each wider arm's study fold
# file (study.py); expect.tsv (expect.py). A wider arm whose cohort equals
# another's (Enterococcus faecium's cap-4000 and cap-all) shares its root
# and its fits.
#
# Fits, per dataset: the default, grm arm --arm cap-2000 on its root with
# DRAWS + TWIN keyed lineage slots and the species grouping of the cap-all
# root (--groups); then each wider arm, --pairing
# common_outer_test_set_balanced against the default, DRAWS slots; each
# with --trace-training; then check_folds.py holds every fit's folds to the
# design (the scored genomes the default's, no added genome scored or
# trained beside a scored genome of its lineage or duplicate group) and
# writes the counts. Then on each study dataset grm learn --fold-file at
# seed 42, the default on the reviewed plan, each wider arm on its own file,
# with --exp-pick at no inner fold. A wider arm that adds no training genome
# to a study dataset is not fitted there, as expect.tsv records, and a
# wider arm with no added genome labelled on a dataset (added-labelled.tsv,
# none_added 1) is fitted on neither split there.
#
# Lineage fits at once: each dataset's chain (the default's fit, each wider
# arm's, the fold check) is one job. The datasets run largest first, by the
# cap-all cohort's labelled genomes; the first runs alone, and its peak
# memory sets N, the chains at once after it: the smaller of AT_ONCE and
# (the box's memory less 4 GiB) / (the largest fit peak so far), at least 1,
# read again before each chain starts. N is written to results/at-once.txt.
# A fit killed for memory fails the run once the running chains end; the
# next start runs at half the N in force then. Study fits run one at a time.
#
# OUT gets logs/run.log, logs/<step>.time and .log (each with genome ids
# written as [genome id]; the raw logs stay in work/raw-logs), logs/disk.tsv,
# logs/timing.tsv (per step and fit: wall seconds, CPU seconds and peak
# memory), flags/, work/ (release, cohorts, fasta, roots, fits, plans,
# study, picks: these hold genome ids or scores and stay on the box),
# results/ (counts only: cohorts.json with each cohort list's line count and
# sha256, added-labelled.tsv, bioprojects.json, fetch.json, roots.json,
# folds tables, study-files tables, at-once.txt), score/ (score.py's tables and its output, score.log,
# kept out of logs/ and run.log), binary.txt, lists.txt and expect.tsv. The
# end checks that no file under results/, score/ or logs/ names a genome of
# the species' cohorts, then writes SHA256SUMS and OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> fits <n> SHA256SUMS <sha256>
# Any failure writes OUT/FAILED naming the step and exits 1. It stops
# nothing but itself.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?datasets}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}; DRAWS=${6:-3}; TWIN=${7:-3}; ARMS=${8:-cap-all}
SPLITS=${9:-lineage,study}; SCORE=${10:-yes}; AT_ONCE=${11:-1}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/c2-cap
W=$OUT/work; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; COH=$W/cohorts; ROOTS=$W/roots; FITS=$W/fits; UNREAD=$W/unread
PLANS=$W/plans; STUDY=$W/study; PICKS=$W/picks; GROUPING=$W/roots/cap-all/species/$SP/structure.json
export RUSTUP_TOOLCHAIN=1.98.1 CARGO_BUILD_JOBS=$JOBS RAYON_NUM_THREADS=$JOBS
export PATH=$HOME/.cargo/bin:/usr/local/cargo/bin:$PATH
exec 9> "$OUT/run.lock"
flock -n 9 || { echo "$(date -u +%FT%TZ) another start runs on this OUT; this one ends"; exit 1; }
STEP=start
say() { echo "$(date -u +%FT%TZ) $*"; }
disk() { echo -e "$(date -u +%FT%TZ)\t$1\t$(du -sb "$OUT" | cut -f1)\t$(df -B1 --output=avail "$OUT" | tail -1 | tr -d ' ')" >> "$OUT/logs/disk.tsv"; }
fail() {
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
# Genome ids stay on the box. grm names a genome in some log lines, and
# logs/ is copied back. So each step's and fit's log moves to work/raw-logs,
# and logs/ gets it back with every genome id of the release's list and of
# the cohorts written as [genome id] (scrub.py).
scrub() {  # scrub <log>
  local raw; raw=$W/raw-logs/$(basename "$1")
  mkdir -p "$W/raw-logs"; mv -f "$1" "$raw"
  python3 -I -B "$HERE/scrub.py" "$R/species/$SP.json" "$COH" < "$raw" > "$1"
}
killed_for_memory() { grep -qE 'terminated by signal 9|non-zero status 137' "$1"; }
TIME_FORMAT="%e s wall, %M KB peak, %U s user, %S s sys"
step() {  # step <name> <function>: once, timed, in a shell that stops at the first failure
  STEP=$1
  if [ -e "$OUT/flags/$1.done" ]; then say "$1: done before, skipped"; return 0; fi
  say "$1: start"
  /usr/bin/time -f "$1 $TIME_FORMAT" -o "$OUT/logs/$1.time" \
    bash -euo pipefail -c "$(declare -f say pinned master make_root has study_listed "$2"); $2" < /dev/null > "$OUT/logs/$1.log" 2>&1 \
    || { scrub "$OUT/logs/$1.log"; tail -n 30 "$OUT/logs/$1.log"
         ! killed_for_memory "$OUT/logs/$1.time" \
           || say "$1: killed by signal 9, as the kernel kills a program when memory runs out; restart with a smaller JOBS"
         false; }
  scrub "$OUT/logs/$1.log"
  touch "$OUT/flags/$1.done"
  disk "$1"
  say "$1: end, $(cat "$OUT/logs/$1.time")"
}
has() { echo ",$1," | grep -q ",$2,"; }
study_listed() {  # the listed datasets with a reviewed study plan
  local d
  for d in $(echo "$DATASETS" | tr ',' '\n'); do
    if grep -q "/folds/$d.study-whole.json\$\|  folds/$d.study-whole.json\$" "$KIT/study-plans.sha256"; then echo "$d"; fi
  done
}
# make_root <arm>: a data root holding the arm's cohort, the shared curated
# table, catalog and genomes linked in.
make_root() {
  local r=$ROOTS/$1
  mkdir -p "$r/species/$SP" "$r/bvbrc"
  for part in curated catalog fasta; do [ -e "$r/$part" ] || ln -s "../../$part" "$r/$part"; done
  [ -e "$r/bvbrc/genomes.jsonl" ] || : > "$r/bvbrc/genomes.jsonl"
  cp "$COH/$1.txt" "$r/species/$SP/genomes.txt"
}
export KIT SP DATASETS OUT JOBS DRAWS TWIN ARMS SPLITS SCORE AT_ONCE REPO HERE W R M GRM COH ROOTS FITS UNREAD PLANS STUDY PICKS GROUPING

STEP=kit
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
[ "$(ulimit -n)" = unlimited ] || [ "$(ulimit -n)" -ge 4096 ] \
  || { say "the open-file limit is $(ulimit -n), below the 4,096 grm matrix needs here"; false; }
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the binary must be its commit"; false; }
python3 -I -c 'import numpy' 2> /dev/null || { say "python3 has no numpy, which score.py needs; install python3-numpy"; false; }
for tool in cargo objcopy unshare flock /usr/bin/time curl tar xz zstd sha256sum stat; do
  command -v "$tool" > /dev/null || { say "$tool is not installed, which a later step needs"; false; }
done
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/curated.tar" d77e2c6f75f8d986268d0c075cb2ea364fba9972fa22a55139495ac0eb376d97
pinned "$KIT/curated.sha256" 8b0e43dd98ca1b3fc7609b7e25b1b52d4bb4dd0c4c555cf5e771622d24fd9c6e
pinned "$KIT/fold-plans.tar.xz" b13cad7dbf6318fd175a1d4a1fb1715e7763f19030b48338cbd2cd3f14ca71ec
pinned "$KIT/study-plans.sha256" c9c902f2c9e9b808b8867aff129af1ea736e52dedc93ae1d9e3e631cf39654a1
pinned "$KIT/marks.tsv" ba662467864edce9e8cc617651b81598feac5f32bb1219d7409467ff094aa313
pinned "$KIT/master-genome-meta.jsonl.gz" 84ae57a777cb234aa0a39102ab7f9ca9448e52d1474654179f26165644c25f64
echo "$DRAWS" | grep -Eq '^[1-9][0-9]*$' || { say "DRAWS is a count of draws, 1 or more"; false; }
echo "$TWIN" | grep -Eq '^[0-9]+$' || { say "TWIN is a count of twin slots, 0 or more"; false; }
echo "$ARMS" | grep -Eq '^(cap-4000|cap-all)(,(cap-4000|cap-all))?$' || { say "ARMS is cap-4000, cap-all or both"; false; }
echo "$SPLITS" | grep -Eq '^(lineage|study)(,(lineage|study))?$' || { say "SPLITS is lineage, study or both"; false; }
case $SCORE in yes|no) ;; *) say "SCORE is yes or no"; false ;; esac
echo "$AT_ONCE" | grep -Eq '^[1-4]$' || { say "AT_ONCE is 1 to 4 datasets at once"; false; }
# The listed datasets: the species' own, each once, each in the first run
# (CAP, not M, by lists.tsv and by the label count's marks.tsv).
python3 -I -B "$HERE/lists_check.py" "$HERE/lists.tsv" "$KIT/marks.tsv" "$SP" "$DATASETS" \
  || { say "the list must hold first-run datasets of $SP only, each once"; false; }
printf 'species\t%s\ndatasets\t%s\ndraws\t%s\ntwin\t%s\narms\t%s\nsplits\t%s\nscore\t%s\n' \
  "$SP" "$DATASETS" "$DRAWS" "$TWIN" "$ARMS" "$SPLITS" "$SCORE" > "$OUT/lists.txt.part"
if [ -s "$OUT/lists.txt" ] && ! cmp -s "$OUT/lists.txt" "$OUT/lists.txt.part"; then
  say "the dataset list, the draws or the arms differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/lists.txt.part" "$OUT/lists.txt"
rm -f "$OUT/FAILED" "$OUT/DONE"
say "datasets $(echo "$DATASETS" | tr ',' '\n' | wc -l)"

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
  moved=$(git -C "$REPO" diff --name-only "$built" HEAD -- . ':(exclude)scripts/c2-cap/' ':(exclude)docs/' \
    ':(exclude)tests/' ':(exclude).github/')
  [ -z "$moved" ] || { say "the checkout moved since the build, $(echo "$moved" | wc -l) changed paths outside scripts/c2-cap/, docs/, tests/ and .github/; use a fresh OUT"; false; }
  echo "restart $head $(date -u +%FT%TZ)" >> "$OUT/binary.txt"
  say "restart on $head, built at $built"
fi

curated() {
  rm -rf "$W/curated"; mkdir -p "$W/curated"
  tar -xf "$KIT/curated.tar" -C "$W/curated"
  (cd "$W/curated" && sha256sum --quiet -c "$KIT/curated.sha256")
  python3 -I -B "$REPO/scripts/c6-mic/release_specs.py" "$W/curated/datasets.json" "$R" "$DATASETS" "$W/specs"
  cp "$W/specs/datasets.json" "$W/curated/datasets.json"
}
step curated curated

catalog() {
  mkdir -p "$W/catalog"
  base=https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/2026-08-07.1
  for f in AMR_CDS.fa ReferenceGeneCatalog.txt version.txt; do
    [ -s "$W/catalog/$f" ] || { curl --fail --silent --show-error --location --retry 3 -A curl -o "$W/catalog/$f.part" "$base/$f"; mv "$W/catalog/$f.part" "$W/catalog/$f"; }
    master "catalog/$f" "$W/catalog/$f" || { echo "catalog/$f differs from the master list"; exit 1; }
  done
  echo "catalog: 3 of 3 files equal the master list"
}
step catalog catalog

cohorts() {
  rm -rf "$COH"
  python3 -I -B "$HERE/cohorts.py" "$W/curated" "$R/species/$SP.json" "$SP" "$COH" 4000 4096
  cp "$COH/cohorts.json" "$OUT/results/cohorts.json"
  cp "$COH/added-labelled.tsv" "$OUT/results/added-labelled.tsv"
}
step cohorts cohorts

# The added genomes' BioProjects, by id, from BV-BRC's genome table: two
# fields, genome_id and bioproject_accession. The answers stay in work/;
# results/bioprojects.json holds the endpoint, the date and the counts.
bioprojects() {
  if ! has "$SPLITS" study || [ -z "$(study_listed)" ]; then
    echo "bioprojects: no listed dataset has a study plan; nothing asked"; return 0
  fi
  python3 -I -B "$HERE/bioprojects.py" "$COH" "$KIT/master-genome-meta.jsonl.gz" "$W/bioprojects.part.jsonl.gz" \
    "$OUT/results/bioprojects.json"
  mv "$W/bioprojects.part.jsonl.gz" "$W/bioprojects.jsonl.gz"
}
step bioprojects bioprojects
META=$KIT/master-genome-meta.jsonl.gz
[ ! -e "$W/bioprojects.jsonl.gz" ] || META=$META,$W/bioprojects.jsonl.gz
export META

# Every cap-all genome from BV-BRC, into the shared fasta folder.
fetch() {
  mkdir -p "$W/fasta"; make_root cap-all
  "$GRM" --data "$ROOTS/cap-all" fetch fasta --species "$SP" --workers 8
  python3 -I -B "$HERE/fetch_check.py" "$R/species/$SP.json" "$COH" "$W/fasta" "$M" "$KIT/master-genome-meta.jsonl.gz" \
    "$ROOTS/cap-all/species/$SP/fasta_failures.tsv" "$OUT/results/fetch.json" "$W/fetched.sha256"
  make_root cap-all
}
step fetch fetch

# The roots: cap-all counts its ksets, each other root takes them as links.
# A wider arm whose cohort equals another's shares its root and its fits
# (roots.json: roots, the arms fitted and the root each fits on, same_as).
STEP=roots
python3 -I -B "$HERE/roots.py" plan "$COH" "$ARMS" "$OUT/results/roots.json"
OWN_ROOTS=$(python3 -I -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["roots"]))' "$OUT/results/roots.json")
FIT_ARMS=$(python3 -I -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["fits"]))' "$OUT/results/roots.json")
export OWN_ROOTS FIT_ARMS
root_of() { python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["fits"][sys.argv[2]])' "$OUT/results/roots.json" "$1"; }
count_cap_all() { unshare -rn "$GRM" --data "$ROOTS/cap-all" count --species "$SP"; }
step count-cap-all count_cap_all
for arm in $OWN_ROOTS; do
  [ "$arm" != cap-all ] || continue
  fn=${arm//-/_}
  eval "count_$fn() { make_root $arm
    python3 -I -B \"\$HERE/roots.py\" share \"\$ROOTS/cap-all\" \"\$ROOTS/$arm\" \"\$SP\" \"\$W/linked-$arm.txt\"
    unshare -rn \"\$GRM\" --data \"\$ROOTS/$arm\" count --species \"\$SP\"
    python3 -I -B \"\$HERE/roots.py\" record \"\$ROOTS/cap-all\" \"\$ROOTS/$arm\" \"\$SP\" \"\$W/linked-$arm.txt\"; }"
  step "count-$arm" "count_$fn"
done
for arm in $OWN_ROOTS; do
  fn=${arm//-/_}
  eval "structure_$fn() { unshare -rn \"\$GRM\" --data \"\$ROOTS/$arm\" structure --species \"\$SP\"; }"
  step "structure-$arm" "structure_$fn"
  eval "matrix_$fn() { unshare -rn \"\$GRM\" --data \"\$ROOTS/$arm\" matrix --species \"\$SP\"; }"
  step "matrix-$arm" "matrix_$fn"
done
STEP=roots-check
python3 -I -B "$HERE/roots.py" counts "$ROOTS" "$SP" "$OUT/results/roots.json"

plans() {
  rm -rf "$PLANS"; mkdir -p "$PLANS/reviewed/folds"
  : > "$PLANS/study.txt"
  has "$SPLITS" study || { echo "study plans: the study split is not run"; return 0; }
  names=()
  for d in $(study_listed); do names+=("folds/$d.study-whole.json"); done
  if [ "${#names[@]}" = 0 ]; then echo "study plans: none of the datasets has one"; return 0; fi
  tar -xJf "$KIT/fold-plans.tar.xz" -C "$PLANS/reviewed" "${names[@]}"
  for f in "${names[@]}"; do
    want=$(grep -m1 "/$f\$\|  $f\$" "$KIT/study-plans.sha256" | cut -c1-64)
    [ -n "$want" ] && [ "$(sha256sum "$PLANS/reviewed/$f" | cut -c1-64)" = "$want" ] \
      || { echo "$f differs from the reviewed plan's sha256"; exit 1; }
  done
  printf '%s\n' "${names[@]}" | sed 's|^folds/||; s|\.study-whole\.json$||' | LC_ALL=C sort > "$PLANS/study.txt"
  echo "study plans: ${#names[@]}, each equal to its sha256"
  for arm in $FIT_ARMS; do
    python3 -I -B "$HERE/study.py" "$PLANS/reviewed/folds" "$PLANS/study.txt" "$COH" "$arm" "$W/curated" \
      "$META" "$GROUPING" "$PLANS/$arm" "$OUT/results/study-files-$arm.tsv"
  done
}
step plans plans

STEP=expect
python3 -I -B "$HERE/expect.py" "$SP" "$DATASETS" "$PLANS" "$DRAWS" "$TWIN" "$ARMS" "$SPLITS" "$OUT/results/roots.json" \
  "$OUT/results/added-labelled.tsv" > "$OUT/expect.tsv.part"
if [ -s "$OUT/expect.tsv" ] && ! cmp -s "$OUT/expect.tsv" "$OUT/expect.tsv.part"; then
  say "the expected results differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"
say "expect: $(($(wc -l < "$OUT/expect.tsv") - 1)) rows"

# One fit: fit <name> <grm arguments...>, timed into logs/<name>.time, done once.
fit() {
  local name=$1; shift
  STEP=$name
  if [ -e "$OUT/flags/$name.done" ]; then say "$name: done before, skipped"; return 0; fi
  say "$name: start"
  /usr/bin/time -f "$name $TIME_FORMAT" -o "$OUT/logs/$name.time" "$GRM" "$@" \
    < /dev/null > "$OUT/logs/$name.log" 2>&1 || { scrub "$OUT/logs/$name.log"; tail -n 30 "$OUT/logs/$name.log"
      ! killed_for_memory "$OUT/logs/$name.time" \
        || { touch "$OUT/flags/memory-killed"
             say "$name: killed by signal 9, as the kernel kills a program when memory runs out; the next start runs half as many fits at once"; }
      false; }
  scrub "$OUT/logs/$name.log"
  touch "$OUT/flags/$name.done"
  disk "$name"
  say "$name: end, $(cat "$OUT/logs/$name.time")"
}
# Moves a fit's score files where nothing reads them (SCORE no).
unread() {  # unread <arm> <dataset>
  [ "$SCORE" = no ] || return 0
  mkdir -p "$UNREAD/$1/$2"
  for f in predictions.jsonl runs.json repeats.json; do
    [ ! -e "$FITS/$1/$2/$f" ] || mv -f "$FITS/$1/$2/$f" "$UNREAD/$1/$2/$f"
  done
}
# none_added <arm> <dataset>: the arm adds no genome labelled on the dataset (added-labelled.tsv).
none_added() {
  awk -F'\t' -v a="$1" -v d="$2" 'NR > 1 && $1 == d && $4 == a && $10 == 1 { f = 1 } END { exit !f }' \
    "$OUT/results/added-labelled.tsv"
}
# chain <dataset>: the dataset's lineage fits, one after another, then its fold check.
chain() {
  local d=$1 arm name
  name=fit-cap-2000-$d
  [ -e "$OUT/flags/$name.done" ] || rm -rf "${FITS:?}/cap-2000/${d:?}"
  fit "$name" --data "$ROOTS/cap-2000" arm --species "$SP" --dataset "$d" --arm cap-2000 --blocks K \
    --repeats $((DRAWS + TWIN)) --groups "$GROUPING" --trace-training --out "$FITS/cap-2000/$d"
  for arm in $FIT_ARMS; do
    name=fit-$arm-$d
    if none_added "$arm" "$d"; then say "$name: no added genome labelled; none-added, not fitted"; continue; fi
    [ -e "$OUT/flags/$name.done" ] || rm -rf "${FITS:?}/${arm:?}/${d:?}"
    fit "$name" --data "$ROOTS/$(root_of "$arm")" arm --species "$SP" --dataset "$d" --arm "$arm" --blocks K \
      --repeats "$DRAWS" --pairing common_outer_test_set_balanced --reference "$FITS/cap-2000/$d" \
      --groups "$GROUPING" --trace-training --out "$FITS/$arm/$d"
  done
  STEP=folds-$d
  python3 -I -B "$HERE/check_folds.py" lineage "$FITS" "$d" "$COH" "$GROUPING" "$OUT/results/roots.json" \
    "$DRAWS" "$TWIN" "$OUT/results/folds-lineage-$d.tsv" "$OUT/results/added-labelled.tsv" > "$OUT/logs/folds-$d.log" 2>&1 \
    || { scrub "$OUT/logs/folds-$d.log"; cat "$OUT/logs/folds-$d.log"; false; }
  scrub "$OUT/logs/folds-$d.log"; cat "$OUT/logs/folds-$d.log"
  for arm in cap-2000 $FIT_ARMS; do unread "$arm" "$d"; done
}
# at_once: the chains to run at once now, from the largest fit peak so far.
at_once() {
  local total peak n
  total=$(awk '/^MemTotal:/ { print $2 }' /proc/meminfo)
  peak=$(cat "$OUT"/logs/fit-*.time 2> /dev/null | sed -n 's/.* s wall, \([0-9]*\) KB peak.*/\1/p' | sort -n | tail -n 1)
  [ -n "$peak" ] && [ "$peak" -gt 0 ] || { echo 1; return 0; }
  n=$(( (total - 4194304) / peak ))
  [ "$n" -ge 1 ] || n=1
  [ "$n" -le "$LIMIT" ] || n=$LIMIT
  echo "$n"
}
if has "$SPLITS" lineage; then
  STEP=at-once
  LIMIT=$AT_ONCE
  [ ! -s "$OUT/results/at-once.txt" ] || [ "$(cat "$OUT/results/at-once.txt")" -ge "$LIMIT" ] \
    || LIMIT=$(cat "$OUT/results/at-once.txt")
  if [ -e "$OUT/flags/memory-killed" ]; then
    LIMIT=$(( LIMIT / 2 )); [ "$LIMIT" -ge 1 ] || LIMIT=1
    rm -f "$OUT/flags/memory-killed"
    say "a fit was killed for memory at the last start; this start runs at most $LIMIT at once"
  fi
  echo "$LIMIT" > "$OUT/results/at-once.txt"
  # Largest first, by the cap-all cohort's labelled genomes on the dataset.
  ORDER=$(for d in $(echo "$DATASETS" | tr ',' '\n'); do
    echo "$(awk -F'\t' -v d="$d" 'NR > 1 && $1 == d && $5 > m { m = $5 } END { print m + 0 }' \
      "$OUT/results/added-labelled.tsv") $d"
  done | sort -k1,1nr -k2,2 | cut -d' ' -f2)
  set -- $ORDER
  chain "$1"; shift
  failed=0
  for d in "$@"; do
    while [ "$(jobs -rp | wc -l)" -ge "$(at_once)" ]; do wait -n || failed=1; done
    [ "$failed" = 0 ] || break
    n=$(at_once)
    [ "$(cat "$OUT/results/at-once.txt")" = "$n" ] || { echo "$n" > "$OUT/results/at-once.txt"; say "lineage fits: $n at once"; }
    chain "$d" &
  done
  while [ "$(jobs -rp | wc -l)" -gt 0 ]; do wait -n || failed=1; done
  STEP=lineage-fits
  [ "$failed" = 0 ] || { say "a dataset's lineage fits failed; the chains running then were let finish"; false; }
fi

if has "$SPLITS" study; then
  STEP=folds-study
  python3 -I -B "$HERE/check_folds.py" study "$PLANS" "$COH" "$META" "$GROUPING" "$FIT_ARMS" \
    "$OUT/results/folds-study.tsv" > "$OUT/logs/folds-study.log" 2>&1 \
    || { scrub "$OUT/logs/folds-study.log"; cat "$OUT/logs/folds-study.log"; false; }
  scrub "$OUT/logs/folds-study.log"; cat "$OUT/logs/folds-study.log"
  while read -r d; do
    [ -n "$d" ] || continue
    mkdir -p "$STUDY/cap-2000" "$PICKS/cap-2000"
    fit "study-cap-2000-$d" --data "$ROOTS/cap-2000" learn --species "$SP" \
      --fold-file "$PLANS/reviewed/folds/$d.study-whole.json" --fold-out "$STUDY/cap-2000" \
      --exp-pick "$PICKS/cap-2000" --exp-pick-inner 0 --force
    # One at a time: the proof that fits at once write the same files ran on grm arm only.
    for arm in $FIT_ARMS; do
      [ -e "$PLANS/$arm/$d.study-whole.json" ] || { say "study-$arm-$d: the arm adds no training genome, not fitted"; continue; }
      root=$(root_of "$arm")
      mkdir -p "$STUDY/$arm" "$PICKS/$arm"
      fit "study-$arm-$d" --data "$ROOTS/$root" learn --species "$SP" \
        --fold-file "$PLANS/$arm/$d.study-whole.json" --fold-out "$STUDY/$arm" \
        --exp-pick "$PICKS/$arm" --exp-pick-inner 0 --force
    done
  done < "$PLANS/study.txt"
fi

STEP=timing
python3 -I -B "$HERE/timing.py" "$OUT/logs" > "$OUT/logs/timing.tsv"

if [ "$SCORE" = yes ]; then
  STEP=score
  say "score: start"
  # score.py's output stays with its tables in score/, which no one opens
  # before every species has run: not in logs/ and not in run.log.
  python3 -I -B "$HERE/score.py" "$OUT/expect.tsv" "$FITS" "$PICKS" "$OUT/score" > "$OUT/score/score.log" 2>&1 \
    || { scrub "$OUT/score/score.log"; say "score: stopped; its output is in score/score.log"; false; }
  scrub "$OUT/score/score.log"
  say "score: end"
fi

STEP=done
cd "$OUT"
python3 -I -B "$HERE/scrub.py" --scan "$R/species/$SP.json" "$COH" results score logs
{ find results score -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.txt' \); echo binary.txt; echo expect.tsv; echo lists.txt; echo logs/timing.tsv; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) fits $(find flags -name 'fit-*.done' -o -name 'study-*.done' | wc -l) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
