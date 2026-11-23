#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARM_SEEDS] [ENET_JOBS] [ENET]
#
# The outside-learner run of one species (plan block C7, docs/outside-learners.md)
# with no agent attached: no prompt, one log, and a restart that skips every
# finished step and every finished result. Run it from a clean checkout of
# this branch:
#
#   setsid nohup bash scripts/c7-outside/run_species.sh KIT SPECIES DATASETS OUT 16 > /dev/null 2>&1 < /dev/null &
#
# DATASETS is a comma list of the species' dataset ids in C7's first run
# (scripts/c7-outside/lists.tsv, first_run 1, so none is M). KIT is C3's kit
# (scripts/c3-inputs/run_species.sh lists it; the unitig arm with GRAPH grm
# needs no outside program). SEEDS (default 42,43,44) are the default's
# draws, TWIN (default 45,46,47) the twin's, ARM_SEEDS (default SEEDS) the
# library arms' (the Rules fix three draws). ENET_JOBS is how many elastic net
# fits run at once (default: as many as three quarters of the box's memory
# holds, from 1 to JOBS; see the fits step); the boosters run JOBS at once.
# ENET is run (the default) or leave-out: left out, no elastic net fit is
# listed or run, and the count names the arm withdrawn (score.py --count
# --withdrawn elasticnet).
#
# The default (k31) and the twin are a C3 run's. On the species' C3 machine
# OUT/c3 is a link to that run's folder (the lead's ruling on question 11):
# nothing of stage 1 runs again, and the run refuses unless that run ended
# DONE, built at 8bfd4ca (binary.txt's commit), with its SHA256SUMS that
# DONE names and every line of it outside score/ holding; nothing under its
# score/ is opened. On a new box (OUT/c3 not a link) stage 1 is C3's run as
# it is (question 12): scripts/c3-inputs/run_species.sh KIT SPECIES DATASETS
# OUT/c3 JOBS SEEDS TWIN unitigs "" grm, which builds grm and the data root,
# the unitig matrix with GRM's graph, the default and the twin by grm learn
# --exp-pick at no inner fold, and writes OUT/c3/DONE. A restart skips it
# when OUT/c3/DONE exists, and otherwise restarts it.
#
# Then, each once, in OUT: libraries (a venv of the box's python3 with the
# wheels of requirements.txt, each held to its sha256 on pypi.org, and the
# versions read back from the imports); build (grm from this checkout into
# work/target for the export alone, since the C3 run's binary has no
# --seed; export-binary.txt records it); export (grm export --seed s of the
# listed datasets with the C3 run's data root, no network, into
# work/export-s<s> for each seed of ARM_SEEDS, every listed dataset's folder
# present); expect.tsv (expect.py, the study datasets those with a reviewed
# study plan in stage 1); plans (plans.py: each export's lineage plan at seed
# s equals the default's folds at seed s by class, and at seed 42 the
# export's input fingerprint equals the release report's, or the run stops);
# the fits, one outside.py fit per arm, dataset, split and seed whose result
# is missing, JOBS at once, one thread each, timed by GNU time; score
# (score.py, its output in score/score.log, kept out of logs/); the genome-id
# scan (scan.py) of results/, score/ and logs/; SHA256SUMS and OUT/DONE:
#   DONE <species> <time> grm <sha256> code <sha256> export-grm <sha256> elasticnet <run|leave-out> results <n> SHA256SUMS <sha256>
# A fit's log goes to work/raw-logs and logs/fits/ gets it back with each
# genome id of the species written as [genome id]. Any failure writes
# OUT/FAILED naming the step and exits 1; a fit that fails stops the run after
# the fits running beside it end. It stops nothing but itself.
#
# OUT gets logs/run.log, logs/<step>.log and .time, logs/fits/ (each fit's
# log, GNU time record and timing file: wall and CPU seconds and peak memory),
# logs/timing.tsv, results/outside/<arm>/ (one report per dataset, split and
# seed: counts and calls, no genome id), score/, libraries.txt, lists.txt,
# expect.tsv, SHA256SUMS, DONE, and work/ (the exports, the venv, the raw
# logs: these hold genome ids and stay on the box).
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?datasets}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}; SEEDS=${6:-42,43,44}; TWIN=${7:-45,46,47}; ARM_SEEDS=${8:-$SEEDS}
ENET_JOBS=${9:-}; ENET=${10:-run}
# The commit every C3 run was built at; a C3 run read on its own machine must be of it.
C3_COMMIT=8bfd4ca8a64678c69431a6bc22b869af6d5cb81b
mkdir -p "$OUT"/logs/fits "$OUT"/flags "$OUT"/work "$OUT"/results/outside "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/c7-outside
C3=$OUT/c3; W=$OUT/work; GRM=$C3/work/target/release/grm; ROOT=$C3/work/root; R=$C3/work/release/rev14-274-20261112
PLANS=$C3/work/plans; VENV=$W/venv; PY=$W/venv/bin/python; EXPORT_GRM=$W/target/release/grm
UNITIGS=$C3/work/inputs/unitigs/matrix
exec 9> "$OUT/run.lock"
flock -n 9 || { echo "$(date -u +%FT%TZ) another start runs on this OUT; this one ends"; exit 1; }
STEP=start
say() { echo "$(date -u +%FT%TZ) $*"; }
fail() {
  if [ -e "$OUT/DONE" ]; then say "refused $1 $STEP; the finished run's DONE is kept"
  else say "FAILED $1 $STEP"; echo "FAILED $SP $1 $STEP $(date -u +%FT%TZ)" > "$OUT/FAILED"; fi
  exit 1
}
trap 'fail at' ERR
trap 'fail "stopped at"' TERM INT
say "start $SP, $JOBS jobs, repo $(git -C "$REPO" rev-parse HEAD)"

STEP=check
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the run must be its commit"; false; }
for tool in python3 flock /usr/bin/time unshare sha256sum cargo git; do
  command -v "$tool" > /dev/null || { say "$tool is not installed, which a later step needs"; false; }
done
python3 -I -c 'import ensurepip, venv' 2> /dev/null || { say "python3 cannot make a venv; install python3-venv"; false; }
for list in "$SEEDS" "$TWIN" "$ARM_SEEDS"; do
  echo "$list" | grep -Eq '^[0-9]+(,[0-9]+)*$' || { say "SEEDS, TWIN and ARM_SEEDS are comma lists of seeds"; false; }
done
case "$ENET" in
  run) ARMS=xgboost,lightgbm,elasticnet ;;
  leave-out) ARMS=xgboost,lightgbm ;;
  *) say "ENET is run or leave-out, not $ENET"; false ;;
esac
echo "$JOBS" | grep -Eq '^[1-9][0-9]*$' || { say "JOBS is a whole number from 1"; false; }
[ -z "$ENET_JOBS" ] || echo "$ENET_JOBS" | grep -Eq '^[1-9][0-9]*$' || { say "ENET_JOBS is empty or a whole number from 1"; false; }
python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    first = {r["dataset"] for r in csv.DictReader(f, delimiter="\t") if r["first_run"] == "1"}
listed = sys.argv[2].split(",")
bad = [d for d in listed if d not in first or not d.startswith(sys.argv[3] + "__")]
if bad or len(set(listed)) != len(listed):
    print(f"{len(bad)} listed datasets are not first-run datasets of the species, or one is listed twice: " + ",".join(bad))
    sys.exit(1)
' "$HERE/lists.tsv" "$DATASETS" "$SP" || { say "the list must hold first-run datasets of $SP, each once"; false; }
printf 'species\t%s\ndatasets\t%s\nseeds\t%s\ntwin\t%s\narm_seeds\t%s\narms\t%s\n' "$SP" "$DATASETS" "$SEEDS" "$TWIN" \
  "$ARM_SEEDS" "$ARMS" > "$OUT/lists.txt.part"
if [ -s "$OUT/lists.txt" ] && ! cmp -s "$OUT/lists.txt" "$OUT/lists.txt.part"; then
  say "the dataset list or the seeds differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/lists.txt.part" "$OUT/lists.txt"
rm -f "$OUT/FAILED" "$OUT/DONE"

# scrub <log>: the log moves to work/raw-logs and comes back with every genome id of the species as [genome id].
scrub() {
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
    bash -euo pipefail -c "$(declare -f "$2"); $2" < /dev/null > "$OUT/logs/$1.log" 2>&1 \
    || { scrub "$OUT/logs/$1.log"; tail -n 30 "$OUT/logs/$1.log"; false; }
  scrub "$OUT/logs/$1.log"
  touch "$OUT/flags/$1.done"
  say "$1: end, $(cat "$OUT/logs/$1.time")"
}
export KIT SP DATASETS OUT JOBS SEEDS TWIN ARM_SEEDS ENET_JOBS ENET REPO HERE C3 W GRM EXPORT_GRM ROOT R PLANS VENV PY
export ARMS UNITIGS C3_COMMIT

# c3_stage: the C3 run C7 reads. On the C3 machine (OUT/c3 a link) it is read and checked, never run; on a new box
# it is run when it has no DONE. Either way it must end DONE with its SHA256SUMS holding outside score/.
c3_stage() {
  if [ -L "$C3" ]; then
    [ -e "$C3/DONE" ] || { say "c3: the C3 run at $(readlink "$C3") has no DONE; on its machine stage 1 never runs"; false; }
    [ "$(sed -n 's/^commit //p' "$C3/binary.txt")" = "$C3_COMMIT" ] \
      || { say "c3: the C3 run was built at $(sed -n 's/^commit //p' "$C3/binary.txt"), not $C3_COMMIT"; false; }
    say "c3: the C3 run on this machine, read and not run: $(cat "$C3/DONE")"
  elif [ -e "$C3/DONE" ]; then
    say "c3: done before, skipped: $(cat "$C3/DONE")"
  else
    say "c3: start, scripts/c3-inputs/run_species.sh into $C3 (its log is $C3/logs/run.log)"
    bash "$REPO/scripts/c3-inputs/run_species.sh" "$KIT" "$SP" "$DATASETS" "$C3" "$JOBS" "$SEEDS" "$TWIN" unitigs "" grm \
      < /dev/null > /dev/null 2>&1 || { say "c3: stopped: $(cat "$C3/FAILED" 2> /dev/null || echo 'no FAILED file')"; false; }
    [ -e "$C3/DONE" ] || { say "c3: ended with no DONE"; false; }
    say "c3: end, $(cat "$C3/DONE")"
  fi
  [ "$(sed -n 's/.* SHA256SUMS \([0-9a-f]\{64\}\)$/\1/p' "$C3/DONE")" = "$(sha256sum "$C3/SHA256SUMS" | cut -c1-64)" ] \
    || { say "c3: the C3 run's SHA256SUMS is not the one its DONE names"; false; }
  # Lines of files under score/ are left out: C7 opens nothing there.
  (cd "$C3" && grep -v '^[0-9a-f]\{64\}  score/' SHA256SUMS | sha256sum -c --quiet --strict -) \
    || { say "c3: a file of the C3 run outside score/ is not the one its SHA256SUMS names"; false; }
  [ "$(sed -n 's/^grm //p' "$C3/binary.txt")" = "$(sha256sum "$GRM" | cut -c1-64)" ] \
    || { say "c3: the C3 run's binary changed since its build"; false; }
  say "c3: checked, $(grep -vc '^[0-9a-f]\{64\}  score/' "$C3/SHA256SUMS") files outside score/ held to SHA256SUMS"
}
STEP=c3
c3_stage

libraries() {
  rm -rf "$VENV"
  python3 -m venv "$VENV"
  "$VENV/bin/pip" install --quiet --no-deps --require-hashes -r "$HERE/requirements.txt"
  env -i PATH=/usr/bin:/bin "$PY" -I -c '
import sys, numpy, scipy, sklearn, xgboost, lightgbm
print(f"python {sys.version.split()[0]}")
for m in (numpy, scipy, sklearn, xgboost, lightgbm):
    print(f"{m.__name__} {m.__version__}")
' 2> /dev/null | tee "$OUT/libraries.txt"
}
step libraries libraries

build() {
  cd "$REPO" && CARGO_TARGET_DIR=$W/target cargo build --release --locked --bin grm
  printf 'commit %s\ngrm %s\n' "$(git -C "$REPO" rev-parse HEAD)" "$(sha256sum "$EXPORT_GRM" | cut -c1-64)" \
    | tee "$OUT/export-binary.txt"
}
step build build
[ "$(sed -n 's/^grm //p' "$OUT/export-binary.txt")" = "$(sha256sum "$EXPORT_GRM" | cut -c1-64)" ] \
  || { say "the export's binary changed since its build"; false; }

export_step() {
  for s in $(echo "$ARM_SEEDS" | tr ',' ' '); do
    rm -rf "$W/export-s$s"
    unshare -rn "$EXPORT_GRM" --data "$ROOT" export --species "$SP" --only "$DATASETS" --seed "$s" --out "$W/export-s$s" \
      --force
    for d in $(echo "$DATASETS" | tr ',' ' '); do
      test -f "$W/export-s$s/$d/dataset.json" || { echo "the export has no $d"; exit 1; }
    done
    echo "export seed $s: $(echo "$DATASETS" | tr ',' '\n' | wc -l) datasets"
  done
}
step export export_step

STEP=expect
python3 -I -B "$HERE/expect.py" "$SP" "$DATASETS" "$(paste -sd, "$PLANS/study.txt")" "$SEEDS" "$TWIN" "$ARM_SEEDS" \
  "$ARMS" > "$OUT/expect.tsv.part"
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"
say "expect: $(($(wc -l < "$OUT/expect.tsv") - 1)) results, of them the library arms' $(grep -c 'results/outside/' "$OUT/expect.tsv")"

STEP=plans
"$PY" -I -B "$HERE/plans.py" "$OUT/expect.tsv" "$OUT" "$R" "$W/export-s{seed}" > "$OUT/logs/plans.log" 2>&1 \
  || { scrub "$OUT/logs/plans.log"; tail -n 20 "$OUT/logs/plans.log"; say "plans: an export is not the C3 run's folds"; false; }
scrub "$OUT/logs/plans.log"
say "$(tail -n 1 "$OUT/logs/plans.log")"

# fit <arm> <dataset> <scheme> <seed>: one fit, skipped when its report exists.
fit() {
  local arm=$1 d=$2 scheme=$3 s=$4 name out extra=()
  name=$arm-$d-$scheme-s$s; out=$OUT/results/outside/$arm/$d.$scheme-s$s.outside.json
  [ -e "$out" ] && return 0
  mkdir -p "$OUT/results/outside/$arm"
  [ "$arm" = elasticnet ] && extra+=(--unitigs "$UNITIGS")
  [ "$scheme" = study-whole ] && extra+=(--fold-file "$PLANS/reviewed/folds/$d.study-whole.json")
  if /usr/bin/time -v -o "$OUT/logs/fits/$name.time" env -i PATH=/usr/bin:/bin HOME="$HOME" OMP_NUM_THREADS=1 \
      OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$PY" -I -B "$HERE/outside.py" fit --arm "$arm" \
      --export "$W/export-s$s/$d" --scheme "$scheme" --seed "$s" --out "$out" --timing "$OUT/logs/fits/$name.json" \
      "${extra[@]}" < /dev/null > "$OUT/logs/fits/$name.log" 2>&1; then rc=0; else rc=$?; fi
  scrub "$OUT/logs/fits/$name.log"
  echo "$(date -u +%FT%TZ) fit $name rc $rc"
  return "$rc"
}
export -f fit scrub
STEP=fits
awk -F'\t' 'NR > 1 && $8 ~ /^results\/outside\// { print $3, $1, $5, $6 }' "$OUT/expect.tsv" > "$W/fits.txt"
# An elastic net fit holds its cohort's projection of the unitig build, which peaks at about 5.2 times the build's
# patterns.bin: 4,364,420 KB for a SYNTHETIC build of Escherichia coli's size, 3,351,087 patterns over 2,000 genomes
# (docs/outside-learners.md). Unless ENET_JOBS is given, those fits run as many at once as three quarters of the
# box's memory holds at 6 times patterns.bin each, from 1 to JOBS. The boosters read the export's screened rows
# only and run JOBS at once, first.
if [ "$ENET" = leave-out ]; then
  ENET_JOBS=1
elif [ -z "$ENET_JOBS" ]; then
  build_kb=$(( $(stat -c %s "$UNITIGS/patterns.bin") / 1024 )); mem_kb=$(awk '/^MemTotal/ { print $2 }' /proc/meminfo)
  ENET_JOBS=$(( mem_kb * 3 / 4 / (6 * build_kb + 1) ))
fi
[ "$ENET_JOBS" -ge 1 ] || ENET_JOBS=1
[ "$ENET_JOBS" -le "$JOBS" ] || ENET_JOBS=$JOBS
if [ "$ENET" = leave-out ]; then say "fits: start, the boosters $JOBS at once, the elastic net left out"
else say "fits: start, the boosters $JOBS at once, the elastic net $ENET_JOBS at once"; fi
for group in boosters elasticnet; do
  if [ "$group" = boosters ]; then n=$JOBS; pick='$1 != "elasticnet"'; else n=$ENET_JOBS; pick='$1 == "elasticnet"'; fi
  awk "$pick" "$W/fits.txt" | xargs -r -P "$n" -L 1 bash -c 'fit "$@"' fit \
    || { say "fits: a fit failed; its log is under logs/fits"; false; }
done
say "fits: end, $(find "$OUT/results/outside" -name '*.outside.json' | wc -l) reports"

STEP=timing
{ printf 'fit\twall_s\tcpu_s\tpeak_kb\n'
  for f in "$OUT"/logs/fits/*.json; do
    python3 -I -c '
import json, sys
t = json.load(open(sys.argv[1]))
print(f"{sys.argv[2]}\t{t[\"wall_s\"]}\t{t[\"cpu_s\"]}\t{t[\"peak_kb\"]}")
' "$f" "$(basename "$f" .json)"
  done; } > "$OUT/logs/timing.tsv"

STEP=score
say "score: start"
# score.py prints the check counts, or why it stopped. Its output stays with its tables in score/, which no one
# opens before every species has run: not in logs/ and not in run.log.
"$PY" -I -B "$HERE/score.py" "$OUT/expect.tsv" "$OUT" "$R" "$OUT/score" > "$OUT/score/score.log" 2>&1 \
  || { scrub "$OUT/score/score.log"; say "score: stopped; its output is in score/score.log"; false; }
scrub "$OUT/score/score.log"
say "score: end"

STEP=done
cd "$OUT"
python3 -I -B "$HERE/scan.py" "$R/species/$SP.json" results score logs
{ find results score -type f \( -name '*.json' -o -name '*.tsv' \); echo libraries.txt; echo expect.tsv; echo lists.txt
  echo export-binary.txt; echo logs/timing.tsv; } | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
line="DONE $SP $(date -u +%FT%TZ) grm $(sed -n 's/^grm //p' "$C3/binary.txt") code $(sed -n 's/^code //p' "$C3/binary.txt") export-grm $(sed -n 's/^grm //p' export-binary.txt) elasticnet $ENET results $(find results/outside -name '*.outside.json' | wc -l) SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
