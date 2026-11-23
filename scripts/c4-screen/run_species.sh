#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARMS]
#
# The screen-arm run of one species (plan block C4, docs/screen-arms.md) on
# a fresh box with no agent attached: no prompt, nothing read from another
# machine but the public sources scripts/c3-inputs/run_species.sh names,
# one log, and a restart that skips every finished step and every finished
# result. Run it from a clean checkout of this branch:
#
#   setsid nohup bash scripts/c4-screen/run_species.sh KIT SPECIES DATASETS OUT 16 > /dev/null 2>&1 < /dev/null &
#
# DATASETS is a comma list of the species' dataset ids that are not marked
# M (scripts/c4-screen/lists.tsv, first_run). SEEDS (default 42,43,44) are
# the draws of the default (k31) and of each arm; TWIN (default 45,46,47,
# or - for none) the twin's seeds, the default fitted again. ARMS (default
# cmh-4000,cmh-16000,floor-3) are the arms fitted beside the default:
#   cmh-4000   --exp-screen cmh
#   cmh-16000  --exp-screen cmh --screened-features 16000
#   floor-3    --exp-floor 3
# Every fit runs with --exp-pick-inner 0 (each outer fold's calls per
# learner, no inner pick) and every arm with --force. SEEDS must hold 42,
# the release's seed: the default's seed-42 reports are checked against the
# release's before any arm runs. KIT holds these files, each checked
# against the sha256 pinned here before anything else runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   MASTER-MANIFEST.sha256      the data root's list of files and sha256
#   curated.tar, curated.sha256 the rebuilt revision 14 table's curated files
#   fold-plans.tar.xz, study-plans.sha256
#                               the reviewed study plans per dataset and the
#                               sha256 of each, taken before any fit
#   marks.tsv                   the label count's table: a dataset is M when
#                               pairs_S+I or pairs_I+R is above 0
#
# The steps up to the fits are C3's, run from C3's own script: this one
# reads the functions pinned, master, scrub, step, release, build, curated,
# catalog, genomes, count, structure, matrix, plans, fit, missing, check42
# and lineage_fit out of scripts/c3-inputs/run_species.sh and copies none of
# them. So the release, grm built from this checkout, the data root (the
# curated table with the release's breakpoints, the catalog's three files
# from NCBI, the genomes from BV-BRC, each held to the master list and the
# release's assembly figures, then count, structure and matrix with no
# network) and the study-whole plans are made as C3's run makes them, and
# no AMRFinderPlus program is installed or run.
#
# Fits, per seed of SEEDS: the default, grm learn --seed s --exp-pick
# results/picks/k31, then each arm into results/picks/ARM, each over the
# datasets whose result or inputs file at that seed is missing (--only,
# --force). After the default's seed-42 run, scripts/refit/refit_check.py
# compares each seed-42 report with the release's; any difference stops the
# run before another fit. Then the twin at each seed of TWIN into
# results/picks/twin; then, at the release's seeds, grm learn --fold-file
# for the default and each arm on each dataset with a study plan whose
# result is missing. After each run every result it should have written
# must exist. An arm's study fit that grm learn refuses because an outer
# fold's training genomes give no candidate, where the default's study
# result exists, becomes scripts/c3-inputs/refused.py's record and the run
# goes on (ruling 16, kept as C3 has it); it cannot arise here, since the
# screen arms take the default's candidates and floor 3 gives no fewer. Any
# other failure stops the run. A restart on a later checkout goes on when
# the binary is unchanged and the commits differ only under
# scripts/c4-screen/, docs/, tests/ and .github/; it writes a restart line
# in binary.txt.
#
# OUT gets logs/run.log, logs/<step>.time and .log (each with genome ids
# written as [genome id]; the raw logs stay in work/raw-logs), logs/disk.tsv,
# logs/timing.tsv (per step and fit: wall seconds and peak memory),
# logs/fits.tsv (per arm, split, fold and learner from the inputs files:
# wall and CPU seconds and peak memory), flags/, work/ (release, root with
# the reports, build, plans, the study results: these hold genome ids and
# stay on the box), results/picks/<arm>/ (one result per dataset, split and
# seed, and per arm one inputs file per dataset, split and seed: counts,
# calls and fit times, no genome id), check/ (refit_check.py's tables),
# score/ (score.py's tables and its output, score.log, kept out of logs/
# and run.log), binary.txt, lists.txt and expect.tsv. The end
# checks that no file under results/, check/, score/ or logs/ names a genome
# of the species, then writes SHA256SUMS and OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> results <n> SHA256SUMS <sha256>
# Any failure writes OUT/FAILED naming the step and exits 1. It stops
# nothing but itself. A restart that the kit, list or pin checks refuse
# leaves a finished OUT as it was: DONE stays and no FAILED is written.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?datasets}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}; SEEDS=${6:-42,43,44}; TWIN=${7:-45,46,47}
ARMS=${8:-cmh-4000,cmh-16000,floor-3}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/check "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/c4-screen; C3=$REPO/scripts/c3-inputs
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; PLANS=$W/plans; PICKS=$OUT/results/picks; STUDY=$W/study
REPORTS=$ROOT/species/$SP/results; TABLES_FROM=; AMR=; TABLES=
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

# C3's step functions, read from its script (see the head of this file).
SHARED="pinned master scrub step release build curated catalog genomes count structure matrix plans fit missing check42 lineage_fit"
eval "$(sed -n -e '/^\(pinned\|master\|scrub\|step\|release\|build\|curated\|catalog\|genomes\|plans\|fit\|check42\|lineage_fit\)() {/,/^}$/p' \
  -e '/^\(count\|structure\|matrix\)() {.*; }$/p' -e '/^missing() {/,/; }$/p' "$C3/run_species.sh")"
for f in $SHARED; do
  declare -F "$f" > /dev/null || { say "scripts/c3-inputs/run_species.sh has no function $f"; false; }
done
export KIT SP DATASETS OUT JOBS SEEDS TWIN ARMS TABLES_FROM REPO HERE C3 W ROOT R M GRM PLANS PICKS STUDY AMR TABLES REPORTS

STEP=kit
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
[ "$(ulimit -n)" = unlimited ] || [ "$(ulimit -n)" -ge 4096 ] \
  || { say "the open-file limit is $(ulimit -n), below the 4,096 grm matrix needs here"; false; }
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the binary must be its commit"; false; }
# Every tool a later step needs, refused here rather than after hours of fits.
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
echo "$ARMS" | grep -Eq '^(cmh-4000|cmh-16000|floor-3)(,(cmh-4000|cmh-16000|floor-3))*$' \
  || { say "ARMS is a comma list of cmh-4000, cmh-16000 and floor-3"; false; }
printf 'species\t%s\ndatasets\t%s\nseeds\t%s\ntwin\t%s\narms\t%s\n' "$SP" "$DATASETS" "$SEEDS" "$TWIN" "$ARMS" > "$OUT/lists.txt.part"
if [ -s "$OUT/lists.txt" ] && ! cmp -s "$OUT/lists.txt" "$OUT/lists.txt.part"; then
  say "the dataset list, the seeds or the arms differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/lists.txt.part" "$OUT/lists.txt"
rm -f "$OUT/FAILED" "$OUT/DONE"
say "datasets $(echo "$DATASETS" | tr ',' '\n' | wc -l), arms $ARMS"

step release release
step build build
STEP=build-check
[ "$(sed -n 's/^grm //p' "$OUT/binary.txt")" = "$(sha256sum "$GRM" | cut -c1-64)" ] \
  || { say "the binary changed since the build"; false; }
built=$(sed -n 's/^commit //p' "$OUT/binary.txt"); head=$(git -C "$REPO" rev-parse HEAD)
if [ "$built" != "$head" ]; then
  git -C "$REPO" cat-file -e "$built^{commit}" 2> /dev/null \
    || { say "the checkout moved since the build and the built commit is not in it; use a fresh OUT"; false; }
  moved=$(git -C "$REPO" diff --name-only "$built" HEAD -- . ':(exclude)scripts/c4-screen/' ':(exclude)docs/' \
    ':(exclude)tests/' ':(exclude).github/')
  [ -z "$moved" ] || { say "the checkout moved since the build, $(echo "$moved" | wc -l) changed paths outside scripts/c4-screen/, docs/, tests/ and .github/, the first $(echo "$moved" | head -1); use a fresh OUT"; false; }
  echo "restart $head $(date -u +%FT%TZ)" >> "$OUT/binary.txt"
  say "restart on $head, built at $built: $(git -C "$REPO" diff --name-only "$built" HEAD | wc -l) changed paths, all under scripts/c4-screen/, docs/, tests/ or .github/"
fi
step curated curated
step catalog catalog
step genomes genomes
step count count
step structure structure
step matrix matrix
step plans plans

STEP=expect
python3 -I -B "$HERE/expect.py" "$SP" "$DATASETS" "$(paste -sd, "$PLANS/study.txt")" "$SEEDS" "$TWIN" "$ARMS" \
  > "$OUT/expect.tsv.part"
if [ -s "$OUT/expect.tsv" ] && ! cmp -s "$OUT/expect.tsv" "$OUT/expect.tsv.part"; then
  say "the expected results differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"
say "expect: $(($(wc -l < "$OUT/expect.tsv") - 1)) results"

arm_args() {  # arm_args <arm>: the flags of an arm's fit beside --exp-pick (lineage_fit gives --force)
  case $1 in
    k31 | twin) ;;
    cmh-4000) echo "--exp-screen cmh" ;;
    cmh-16000) echo "--exp-screen cmh --screened-features 16000" ;;
    floor-3) echo "--exp-floor 3" ;;
    *) say "no such arm $1"; return 1 ;;
  esac
}

for s in $(echo "$SEEDS" | tr ',' '\n'); do
  lineage_fit k31 "$s"
  # The gate: no arm runs unless the default's seed-42 reports are the release's.
  [ "$s" != 42 ] || check42
  for arm in $(echo "$ARMS" | tr ',' '\n'); do lineage_fit "$arm" "$s"; done
done
[ "$TWIN" = - ] || for s in $(echo "$TWIN" | tr ',' '\n'); do lineage_fit twin "$s"; done

study_have() {  # study_have <arm> <dataset>: its study result (or an arm's refusal record) and an arm's inputs file exist
  { [ -e "$PICKS/$1/$2.study-whole-s42.pick.json" ] \
    || { [ "$1" != k31 ] && [ -e "$PICKS/$1/$2.study-whole-s42.refused.json" ]; }; } \
    && { [ "$1" = k31 ] || [ -e "$PICKS/$1/$2.study-whole-s42.inputs.json" ] \
      || [ -e "$PICKS/$1/$2.study-whole-s42.refused.json" ]; }
}
# One study fit, as fit, except that an arm's refusal for no candidate,
# where the default's study result exists, becomes C3's refusal record and
# the run goes on (ruling 16). The default's refusal stops the run.
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
      && python3 -I -B "$C3/refused.py" "$OUT/logs/$name.log" "$d" "$arm" "$PICKS/$arm/$d.study-whole-s42.refused.json"; then
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
# Per arm, split and learner, each outer fit's wall and CPU seconds and the
# process's peak after it, from the inputs files: times only.
python3 -I -B "$HERE/fit_times.py" "$PICKS" "$OUT/logs/fits.tsv"

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
{ find results check score -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.txt' \); echo binary.txt; echo expect.tsv; echo lists.txt; echo logs/timing.tsv; echo logs/fits.tsv; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
say "arm results refused for no candidate: $(find results/picks -name '*.refused.json' | wc -l)"
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) results $(find results/picks -name '*.pick.json' | wc -l) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
