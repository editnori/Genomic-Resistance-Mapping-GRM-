#!/usr/bin/env bash
# run_m.sh KIT SPECIES OUT [JOBS] [TABLES]
#
# Plan block C3's run on its datasets of M (docs/input-arms.md, The run on
# the datasets of M; docs/assumption-plan-corrections.md, entry M run) for
# one species, on a fresh box with no agent attached, from a clean checkout
# of this branch:
#
#   setsid nohup bash scripts/c3-inputs/run_m.sh KIT SPECIES OUT 16 > /dev/null 2>&1 < /dev/null &
#
# It is the one driver that admits a dataset of M: run_species.sh, which
# blocks C4 and C7 read, still refuses one. It runs the species' datasets of
# M on scripts/c3-inputs/lists.tsv (in_M), at the first run's seeds (42, 43
# and 44), twin (45, 46 and 47), arms (all five; the program arms on the
# datasets in A) and graph (grm, as every first-run species used), and the
# study split on each of them in S. TABLES is run_species.sh's: a folder of
# the program's tables made before. The data steps, the builds, the program
# and the fit and study functions are run_species.sh's own, read from it.
#
# KIT holds run_species.sh's files, each held to the sha256 it pins, and:
#   genome_amr.jsonl  the 1,285,111 lab records the pinned table (curated.tar)
#                     was rebuilt from on 21 Nov 2026, held to its own
#                     sha256 (ruling R26): the table holds no record's call,
#                     so the list of genomes with an intermediate call is
#                     read from these
#   m-checked.txt     on every species but the check species: the check
#                     species' record that the three checks passed
#
# Steps, in order:
# 0. Before any build: GRM_COMMIT's diff from FIRST_COMMIT must touch the
#    CHANGED files only.
# 1. Same build. FIRST_COMMIT (8bfd4ca, the first run's commit) is built from
#    git archive with rustc 1.98.1, --locked and cargo's home at
#    /home/user/.cargo, the first run's (the registry paths it holds are
#    built into the binary), and its code sha256 (the binary less its build
#    id) must be the species' row of reference/input-arms/page-binaries.tsv.
# 2. The changed binary: GRM_COMMIT, the reviewed commit of the branch
#    run-m/input-arms-grm, FIRST_COMMIT with one change (grm learn
#    --exp-scores), built the same way. Every later step runs it.
# 3. The data steps, the builds, the study plans and the program, as
#    run_species.sh runs them. Then m_list.py, before any fit: the list of
#    genomes with an intermediate call (work/m/list, on this machine only)
#    and the floor table, results/m/floor.tsv, counts only.
# 4. The check species (m_checks.py species: the fewest first-run datasets
#    among those with a dataset of M, Neisseria gonorrhoeae on today's
#    list) fits its first-run datasets with the changed binary. From them
#    score.py's datasets.tsv must give back that species' rows of the
#    committed reference/input-arms/page-datasets.tsv at 12 decimals (Same
#    results), and each fold's counts must come back from the genomes'
#    scores (Same counts). check/m-checked.txt then records the binary's
#    code sha256. Any other species holds KIT/m-checked.txt to its own
#    binary's code sha256 before its first fit.
# 5. The fits of M: the default, each arm and the twin at their seeds, each
#    lineage fit with --exp-scores work/scores/<arm>, each report of the
#    default and the twin copied to work/reports/<arm>/ as it is written; the
#    default's seed-42 reports held to the release's before any arm runs;
#    then the study fits, their fold results kept in work/study/<arm>/.
# 6. Same counts over every fit of M; score.py (the first reading); second.py
#    (the second reading, with no refit); the end's genome-id count; pack.py.
#
# On a difference in Same build, Same results or Same counts the run names
# what differs, writes OUT/FAILED and stops: the lead rules on the entry's
# page before another build, machine or genome set is used.
#
# OUT gets run_species.sh's folders, plus work/scores, work/reports and
# work/m (genome ids: they stay on the box), results/m/floor.tsv,
# check/same-build.txt, check/same-results.txt, check/same-counts.txt,
# check/m-checked.txt on the check species, score/second/, binary-first.txt
# and expect-first.tsv on the check species. The end writes SHA256SUMS,
# OUT/pack.tar (pack.py: every file but work/ and flags/, refused if it would
# carry a score file, a report, a fold result, a list, a compressed or
# non-UTF-8 file, or a genome id in a path or a file) and OUT/DONE:
#   DONE <species> <time> binary <sha256> code <sha256> results <n> pack <sha256>
# Only pack.tar and DONE leave the machine. The pack holds both readings, so
# it goes unopened to the one machine that takes the count over all the
# species, and to no other: it is no page file.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; OUT=${3:?out folder}
JOBS=${4:-$(nproc)}; TABLES_FROM=${5:-}
SEEDS=42,43,44; TWIN=45,46,47; ARMS=program,k31-program,k31-counts,k31-aa,unitigs; GRAPH=grm
FIRST_COMMIT=8bfd4ca8a64678c69431a6bc22b869af6d5cb81b
GRM_COMMIT=f9affda749c3a77d06bb6921c43e22247fe36e09
CHANGED="docs/ci.md src/analysis.rs src/learn/held_scores.rs src/learn/mod.rs src/main.rs tests/rust-test-names.txt"
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results "$OUT"/check "$OUT"/score; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/c3-inputs
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; PLANS=$W/plans; PICKS=$OUT/results/picks; STUDY=$W/study; AMR=$W/amr
TABLES=$AMR/tool/$SP; REPORTS=$ROOT/species/$SP/results
INPUTS=$W/inputs; PRODIGAL=$W/tools/prodigal; BCALM=$W/tools/bcalm-binaries-v2.2.3-Linux/bin/bcalm
SCORES=$W/scores; KEPT=$W/reports; LIST=$W/m/list; FLOOR=$OUT/results/m/floor.tsv
export RUSTUP_TOOLCHAIN=1.98.1 CARGO_BUILD_JOBS=$JOBS RAYON_NUM_THREADS=$JOBS
# Cargo's home is the first run's: the registry paths under it are built into the binary, so the code sha256
# moves with it. Both builds use it, and binary-first.txt and binary.txt name it.
export CARGO_HOME=/home/user/.cargo
export PATH=$HOME/.cargo/bin:/usr/local/cargo/bin:$PATH
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
say "start $SP (the run on M), $JOBS jobs, repo $(git -C "$REPO" rev-parse HEAD)"
[ -s "$OUT/logs/disk.tsv" ] || echo -e "time\tafter\tout_bytes\tfree_bytes" > "$OUT/logs/disk.tsv"
disk start

# run_species.sh's step functions, read from it as block C4's script reads them.
SHARED="pinned master scrub keep_build killed_for_memory step release curated catalog genomes count structure matrix tools inputs_counts prodigal inputs_aa bcalm inputs_unitigs plans program_env program_run program_tables fit arm_args study_have study_fit"
eval "$(sed -n -e '/^\(pinned\|master\|scrub\|keep_build\|killed_for_memory\|step\|release\|curated\|catalog\|genomes\|tools\|inputs_counts\|prodigal\|inputs_aa\|bcalm\|inputs_unitigs\|plans\|program_env\|program_run\|program_tables\|fit\|arm_args\|study_have\|study_fit\)() {/,/^}$/p' \
  -e '/^\(count\|structure\|matrix\)() {.*; }$/p' "$HERE/run_species.sh")"
for f in $SHARED; do
  declare -F "$f" > /dev/null || { say "scripts/c3-inputs/run_species.sh has no function $f"; false; }
done
has_arm() { echo ",$ARMS," | grep -q ",$1,"; }
# scrub reads an id as pack.py does (after "_" or "." too, and inside a longer run of points and digits), so a log
# keeps no id that pack.py's check would refuse. run_species.sh's own scrub is unchanged.
scrub() {  # scrub <log>
  local raw; raw=$W/raw-logs/$(basename "$1")
  mkdir -p "$W/raw-logs"; mv -f "$1" "$raw"
  python3 -I -B "$HERE/pack.py" mask "$R/species/$SP.json" < "$raw" > "$1"
}
# lists.tsv's datasets of a species: lists <column> [<column> ...], each 1.
lists() { python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    rows = [r for r in csv.DictReader(f, delimiter="\t") if r["species_slug"] == sys.argv[2]]
print(",".join(sorted(r["dataset"] for r in rows if all(r[c] == "1" for c in sys.argv[3:]))))
' "$HERE/lists.tsv" "$SP" "$@"; }

STEP=kit
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
[ "$(ulimit -n)" = unlimited ] || [ "$(ulimit -n)" -ge 4096 ] \
  || { say "the open-file limit is $(ulimit -n), below the 4,096 grm matrix needs here"; false; }
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the scripts must be its commit"; false; }
python3 -I -c 'import numpy' 2> /dev/null || { say "python3 has no numpy, which score.py needs; install python3-numpy"; false; }
for tool in cargo objcopy unshare flock /usr/bin/time curl tar xz zstd sha256sum git; do
  command -v "$tool" > /dev/null || { say "$tool is not installed, which a later step needs"; false; }
done
{ mkdir -p "$CARGO_HOME" && [ -w "$CARGO_HOME" ]; } 2> /dev/null \
  || { say "cargo's home $CARGO_HOME, the first run's, cannot be written; both builds need it"; false; }
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/curated.tar" d77e2c6f75f8d986268d0c075cb2ea364fba9972fa22a55139495ac0eb376d97
pinned "$KIT/curated.sha256" 8b0e43dd98ca1b3fc7609b7e25b1b52d4bb4dd0c4c555cf5e771622d24fd9c6e
pinned "$KIT/fold-plans.tar.xz" b13cad7dbf6318fd175a1d4a1fb1715e7763f19030b48338cbd2cd3f14ca71ec
pinned "$KIT/study-plans.sha256" c9c902f2c9e9b808b8867aff129af1ea736e52dedc93ae1d9e3e631cf39654a1
pinned "$KIT/marks.tsv" ba662467864edce9e8cc617651b81598feac5f32bb1219d7409467ff094aa313
pinned "$KIT/amrkit.tar" da07f277a7d271f1ae006919d2c2df4185737526348302f665a9bcb3a083d0e6
pinned "$KIT/tool-tables.sha256" 9291b83f4c7d01a0384f35446e6de8f986dff5f659439ee3d673bccdd09427ba
pinned "$KIT/prodigal.linux" 80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246
# The records the list is read from (m_list.py): the 1,285,111 lab records the pinned table was rebuilt from on
# 21 Nov 2026 (reference/release-summary/README.md), held to their own sha256 (ruling R26), not to a master list line.
[ -f "$KIT/genome_amr.jsonl" ] || { say "the kit has no genome_amr.jsonl, the records m_list.py reads the intermediate calls from"; false; }
pinned "$KIT/genome_amr.jsonl" d231fa768a31fc4ce4c8a6bb7024d05d632621846070688adb892fbfe340fd78
git -C "$REPO" cat-file -e "$FIRST_COMMIT^{commit}" 2> /dev/null || { say "the checkout does not hold the first run's commit $FIRST_COMMIT"; false; }
git -C "$REPO" cat-file -e "$GRM_COMMIT^{commit}" 2> /dev/null || { say "the checkout does not hold the changed commit $GRM_COMMIT"; false; }
# The changed commit's diff from the first run's touches the CHANGED files alone: checked here, before any build.
changed_files() {
  local moved; moved=$(git -C "$REPO" diff --name-only "$FIRST_COMMIT" "$GRM_COMMIT" | LC_ALL=C sort | paste -sd' ' -)
  [ "$moved" = "$(echo "$CHANGED" | tr ' ' '\n' | LC_ALL=C sort | paste -sd' ' -)" ] \
    || { say "the changed commit's diff from the first run's touches: $moved"; return 1; }
}
changed_files
M_DATASETS=$(lists in_M)
[ -n "$M_DATASETS" ] || { say "$SP holds no dataset of M on lists.tsv"; false; }
M_A=$(lists in_M in_A); M_S=$(lists in_M in_S)
# The kit's label count must mark each of them M, as lists.tsv does.
python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    marks = {r["dataset"]: int(r["pairs_S+I"]) > 0 or int(r["pairs_I+R"]) > 0 for r in csv.DictReader(f, delimiter="\t")}
off = [d for d in sys.argv[2].split(",") if not marks.get(d)]
if off:
    print(f"{len(off)} datasets of M on lists.tsv are not marked M in marks.tsv: " + ",".join(off))
    sys.exit(1)
' "$KIT/marks.tsv" "$M_DATASETS" || { say "lists.tsv and the kit's marks.tsv differ on M"; false; }
CHECK_SPECIES=$(python3 -I -B "$HERE/m_checks.py" species "$HERE/lists.tsv" 2> /dev/null)
if [ "$SP" = "$CHECK_SPECIES" ]; then
  FIRST_DATASETS=$(lists first_run); FIRST_A=$(lists first_run_A)
else
  FIRST_DATASETS=; FIRST_A=
  [ -f "$KIT/m-checked.txt" ] || { say "$SP is not the check species $CHECK_SPECIES: the kit needs its m-checked.txt"; false; }
fi
DATASETS=$(echo "$M_DATASETS,$FIRST_DATASETS" | tr ',' '\n' | grep . | LC_ALL=C sort | paste -sd, -)
A_LIST=$(echo "$M_A,$FIRST_A" | tr ',' '\n' | grep . | LC_ALL=C sort | paste -sd, - || true)
if [ -n "$TABLES_FROM" ]; then
  [ -d "$TABLES_FROM/$SP" ] || { say "TABLES has no folder $SP"; false; }
  TABLES_FROM=$(cd "$TABLES_FROM" && pwd)
fi
export KIT SP DATASETS OUT JOBS SEEDS TWIN ARMS TABLES_FROM GRAPH REPO HERE W ROOT R M GRM PLANS PICKS STUDY AMR TABLES REPORTS
export INPUTS PRODIGAL BCALM A_LIST SCORES KEPT LIST FLOOR FIRST_COMMIT GRM_COMMIT CHANGED M_DATASETS
printf 'species\t%s\nm\t%s\nm_a\t%s\nm_s\t%s\nfirst\t%s\nfirst_a\t%s\nseeds\t%s\ntwin\t%s\narms\t%s\ngraph\t%s\nfirst_commit\t%s\ngrm_commit\t%s\n' \
  "$SP" "$M_DATASETS" "$M_A" "$M_S" "$FIRST_DATASETS" "$FIRST_A" "$SEEDS" "$TWIN" "$ARMS" "$GRAPH" "$FIRST_COMMIT" "$GRM_COMMIT" \
  > "$OUT/lists.txt.part"
if [ -s "$OUT/lists.txt" ] && ! cmp -s "$OUT/lists.txt" "$OUT/lists.txt.part"; then
  say "the datasets, seeds, arms or commits differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/lists.txt.part" "$OUT/lists.txt"
rm -f "$OUT/FAILED" "$OUT/DONE"
say "datasets of M $(echo "$M_DATASETS" | tr ',' '\n' | wc -l), in A $(echo "$M_A" | tr ',' '\n' | grep -c . || true), in S $(echo "$M_S" | tr ',' '\n' | grep -c . || true); check species $CHECK_SPECIES"

step release release

# grm from a commit, built from git archive into a target folder of its own. step runs each function in a
# shell of its own that knows say, pinned, master and keep_build only, so each build is written out whole.
first_build() {
  rm -rf "$W/src-first" "$W/target-first"; mkdir -p "$W/src-first"
  git -C "$REPO" archive "$FIRST_COMMIT" | tar -x -C "$W/src-first"
  ( cd "$W/src-first" && CARGO_TARGET_DIR=$W/target-first cargo build --release --locked --bin grm )
  objcopy -O binary --remove-section=.note.gnu.build-id "$W/target-first/release/grm" "$W/grm-first.image"
  { echo "commit $FIRST_COMMIT"; echo "grm $(sha256sum "$W/target-first/release/grm" | cut -c1-64)"
    echo "code $(sha256sum "$W/grm-first.image" | cut -c1-64)"; echo "cargo_home $CARGO_HOME"; } > "$OUT/binary-first.txt"
  rm -rf "$W/grm-first.image" "$W/src-first" "$W/target-first"
  cat "$OUT/binary-first.txt"
}
changed_build() {
  rm -rf "$W/src-changed" "$W/target"; mkdir -p "$W/src-changed"
  git -C "$REPO" archive "$GRM_COMMIT" | tar -x -C "$W/src-changed"
  ( cd "$W/src-changed" && CARGO_TARGET_DIR=$W/target cargo build --release --locked --bin grm )
  objcopy -O binary --remove-section=.note.gnu.build-id "$GRM" "$W/grm.image"
  { echo "commit $GRM_COMMIT"; echo "grm $(sha256sum "$GRM" | cut -c1-64)"
    echo "code $(sha256sum "$W/grm.image" | cut -c1-64)"; echo "cargo_home $CARGO_HOME"; } > "$OUT/binary.txt"
  rm -rf "$W/grm.image" "$W/src-changed"
  cat "$OUT/binary.txt"
}
step same_build first_build
STEP=same-build
python3 -I -B "$HERE/m_checks.py" build "$REPO/reference/input-arms/page-binaries.tsv" "$SP" "$FIRST_COMMIT" \
  "$(sed -n 's/^code //p' "$OUT/binary-first.txt")" "$(sed -n 's/^cargo_home //p' "$OUT/binary-first.txt")" \
  > "$OUT/check/same-build.txt" \
  || { cat "$OUT/check/same-build.txt"; say "Same build differs (cargo's home $CARGO_HOME): the run on M does not start; the lead rules"; false; }
say "$(cat "$OUT/check/same-build.txt")"
step build changed_build
STEP=build-check
[ "$(sed -n 's/^grm //p' "$OUT/binary.txt")" = "$(sha256sum "$GRM" | cut -c1-64)" ] || { say "the binary changed since the build"; false; }
[ "$(sed -n 's/^commit //p' "$OUT/binary.txt")" = "$GRM_COMMIT" ] || { say "binary.txt names another commit"; false; }
[ "$(sed -n 's/^cargo_home //p' "$OUT/binary.txt")" = "$CARGO_HOME" ] || { say "binary.txt names another cargo home"; false; }
CODE=$(sed -n 's/^code //p' "$OUT/binary.txt")
if [ "$SP" != "$CHECK_SPECIES" ]; then
  python3 -I -c '
import sys
got = dict(line.split(" ", 1) for line in open(sys.argv[1]).read().splitlines() if " " in line)
want = {"species": sys.argv[2], "code": sys.argv[3], "commit": sys.argv[4]}
bad = [k for k in want if got.get(k) != want[k]]
if bad:
    print("m-checked.txt differs in " + ", ".join(bad))
    sys.exit(1)
' "$KIT/m-checked.txt" "$CHECK_SPECIES" "$CODE" "$GRM_COMMIT" \
    || { say "the kit's m-checked.txt is not the check species' record for this binary's code sha256; no fit of M runs"; false; }
  say "m-checked: $CHECK_SPECIES passed the three checks with code $CODE"
fi

step curated curated
step catalog catalog
step genomes genomes
step count count
step structure structure
step matrix matrix
step tools tools
step inputs_counts inputs_counts
step prodigal prodigal; step inputs_aa inputs_aa
step inputs_unitigs inputs_unitigs
for kind in counts aa unitigs; do
  STEP=inputs-check-$kind
  ( cd "$INPUTS" && sha256sum --quiet -c "$OUT/logs/inputs-$kind.sha256" ) \
    || { say "the $kind matrix differs from the sha256 its build recorded"; false; }
done
step plans plans
# Each dataset of M in S runs the study split (the entry M run, Arms by block): its reviewed plan must be in the kit.
STEP=study-plans
for d in $(echo "$M_S" | tr ',' '\n'); do
  grep -qx "$d" "$PLANS/study.txt" || { say "$d is in S and the kit holds no reviewed study plan for it"; false; }
done
if [ -e "$OUT/flags/program_run.done" ] || { [ -z "$TABLES_FROM" ] && [ ! -e "$OUT/flags/program_tables.done" ]; }; then
  if [ -n "$A_LIST" ]; then step program_env program_env; step program_run program_run; fi
elif [ -n "$A_LIST" ]; then
  step program_tables program_tables
fi
if [ -n "$A_LIST" ]; then
  STEP=tables-check
  grep "  $SP/" "$KIT/tool-tables.sha256" | sort -k2 > "$AMR/tables.recorded.again"
  ( cd "$AMR/tool" && find "$SP" -name '*.tsv' | sort | xargs -r sha256sum ) | sort -k2 > "$AMR/tables.again"
  cmp -s "$AMR/tables.recorded.again" "$AMR/tables.again" || { say "the program's tables differ from their recorded sha256"; false; }
fi

# The list and the floor, from the label table, before any fit and any score.
m_list() {
  python3 -I -B "$HERE/m_list.py" "$ROOT" "$KIT/genome_amr.jsonl" "$HERE/lists.tsv" "$SP" "$M_DATASETS" "$LIST" "$FLOOR"
}
step m_list m_list

STEP=expect
study_of() { comm -12 <(echo "$1" | tr ',' '\n' | grep . | LC_ALL=C sort) <(LC_ALL=C sort "$PLANS/study.txt") | paste -sd, -; }
python3 -I -B "$HERE/expect.py" "$SP" "$M_DATASETS" "${M_A:--}" "$(study_of "$M_DATASETS")" "$SEEDS" "$TWIN" "$ARMS" > "$OUT/expect.tsv.part"
if [ -n "$FIRST_DATASETS" ]; then
  python3 -I -B "$HERE/expect.py" "$SP" "$FIRST_DATASETS" "${FIRST_A:--}" "$(study_of "$FIRST_DATASETS")" "$SEEDS" "$TWIN" "$ARMS" \
    > "$OUT/expect-first.tsv.part"
  mv "$OUT/expect-first.tsv.part" "$OUT/expect-first.tsv"
fi
if [ -s "$OUT/expect.tsv" ] && ! cmp -s "$OUT/expect.tsv" "$OUT/expect.tsv.part"; then
  say "the expected results differ from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.tsv.part" "$OUT/expect.tsv"

# The datasets of <expect> whose result for an arm, split and seed is missing, as a comma list:
# its pick file, an arm's inputs file, a lineage fit's score file and the default's or the twin's kept report.
m_missing() {  # m_missing <expect> <arm> <split> <seed>
  awk -F'\t' -v a="$2" -v p="$3" -v s="$4" 'NR > 1 && $3 == a && $4 == p && $6 == s { print $1 "\t" $8 "\t" $9 }' "$1" \
    | while IFS=$'\t' read -r d file inputs; do
        { [ -e "$PICKS/$file" ] && { [ "$inputs" = - ] || [ -e "$PICKS/$inputs" ]; } \
          && [ -e "$SCORES/$2/$d.lineage-s$4.scores.json" ] \
          && { [ "$2" != k31 ] && [ "$2" != twin ] || [ -e "$KEPT/$2/$d.lineage-s$4.json" ]; }; } || echo "$d"
      done | sort -u | paste -sd, -; }
# Each report of the default and the twin, copied as its fit writes it: keep_reports <arm> <seed> <datasets>.
keep_reports() {
  local from=$REPORTS
  [ "$2" = 42 ] || from=$ROOT/species/$SP/draws/seed-$2
  mkdir -p "$KEPT/$1"
  for d in $(echo "$3" | tr ',' '\n'); do
    cp "$from/$d.json" "$KEPT/$1/$d.lineage-s$2.json.part" && mv "$KEPT/$1/$d.lineage-s$2.json.part" "$KEPT/$1/$d.lineage-s$2.json"
  done
}
m_lineage_fit() {  # m_lineage_fit <expect> <arm> <seed> <name>: first (the check species' first run) or m
  local arm=$2 seed=$3 todo left
  STEP=fit-$arm-s$seed-$4
  todo=$(m_missing "$1" "$arm" lineage "$seed")
  if [ -z "$todo" ]; then say "fit-$arm-s$seed-$4: all results present, skipped"; return 0; fi
  mkdir -p "$PICKS/$arm" "$SCORES/$arm"
  # shellcheck disable=SC2046
  fit "fit-$arm-s$seed-$4" --only "$todo" --force --seed "$seed" --exp-pick "$PICKS/$arm" \
    --exp-pick-inner 0 --exp-scores "$SCORES/$arm" $(arm_args "$arm")
  case $arm in k31 | twin) keep_reports "$arm" "$seed" "$todo" ;; esac
  left=$(m_missing "$1" "$arm" lineage "$seed")
  [ -z "$left" ] || { say "fit-$arm-s$seed-$4: no result written for $left"; false; }
  disk "fit-$arm-s$seed-$4"
}
m_check42() {  # m_check42 <datasets> <name>: the default's seed-42 reports against the release's
  STEP=check-$2
  python3 -I -B "$REPO/scripts/refit/refit_check.py" "$REPORTS" "$R" "$1" "$OUT/check/$2" > "$OUT/logs/check-$2.log" 2>&1 \
    || { cat "$OUT/logs/check-$2.log"; false; }
  cat "$OUT/logs/check-$2.log"
}
m_study() {  # m_study <datasets> <A datasets>: every study fit, its fold result kept
  while read -r d; do
    [ -n "$d" ] || continue
    for arm in k31 $(echo "$ARMS" | tr ',' '\n'); do
      case $arm in program | k31-program) echo ",$2," | grep -q ",$d," || continue ;; esac
      STEP=study-$arm-$d
      if study_have "$arm" "$d" && { [ -e "$STUDY/$arm/$d.study-whole.json" ] || [ -e "$PICKS/$arm/$d.study-whole-s42.refused.json" ]; }; then
        say "study-$arm-$d: result present, skipped"; continue
      fi
      mkdir -p "$PICKS/$arm" "$STUDY/$arm"
      study_fit "$arm" "$d"
      study_have "$arm" "$d" || { say "study-$arm-$d: no result written"; false; }
      disk "study-$arm-$d"
    done
  done < <(study_of "$1" | tr ',' '\n')
}
m_fits() {  # m_fits <expect> <datasets> <A datasets> <name>
  for s in $(echo "$SEEDS" | tr ',' '\n'); do
    m_lineage_fit "$1" k31 "$s" "$4"
    [ "$s" != 42 ] || m_check42 "$2" "$4"
    for arm in $(echo "$ARMS" | tr ',' '\n'); do
      case $arm in program | k31-program) [ -n "$3" ] || continue ;; esac
      m_lineage_fit "$1" "$arm" "$s" "$4"
    done
  done
  for s in $(echo "$TWIN" | tr ',' '\n'); do m_lineage_fit "$1" twin "$s" "$4"; done
  m_study "$2" "$3"
}
same_counts() {  # same_counts <expect> <name>
  STEP=same-counts-$2
  python3 -I -B "$HERE/m_checks.py" counts "$1" "$PICKS" "$SCORES" "$KEPT" "$STUDY" > "$OUT/check/same-counts-$2.txt" 2>&1 \
    || { scrub "$OUT/check/same-counts-$2.txt"; cat "$OUT/check/same-counts-$2.txt"; say "Same counts differs: the lead rules"; false; }
  scrub "$OUT/check/same-counts-$2.txt"
  say "$(cat "$OUT/check/same-counts-$2.txt")"
}

# The three checks, on the check species, before any fit of M.
if [ "$SP" = "$CHECK_SPECIES" ] && [ ! -e "$OUT/flags/m_checked.done" ]; then
  m_fits "$OUT/expect-first.tsv" "$FIRST_DATASETS" "$FIRST_A" first
  STEP=same-results
  mkdir -p "$OUT/check/first-score"
  python3 -I -B "$HERE/score.py" "$OUT/expect-first.tsv" "$PICKS" "$R" "$OUT/check/first-score" > "$OUT/check/first-score/score.log" 2>&1 \
    || { scrub "$OUT/check/first-score/score.log"; cat "$OUT/check/first-score/score.log"; false; }
  scrub "$OUT/check/first-score/score.log"
  python3 -I -B "$HERE/m_checks.py" results "$REPO/reference/input-arms" "$SP" "$OUT/check/first-score/datasets.tsv" \
    > "$OUT/check/same-results.txt" || { cat "$OUT/check/same-results.txt"; say "Same results differs: the run on M does not start; the lead rules"; false; }
  say "$(cat "$OUT/check/same-results.txt")"
  same_counts "$OUT/expect-first.tsv" first
  printf 'species %s\ncode %s\ncommit %s\nchecked %s\n' "$SP" "$CODE" "$GRM_COMMIT" "$(date -u +%FT%TZ)" > "$OUT/check/m-checked.txt"
  touch "$OUT/flags/m_checked.done"
  say "the three checks passed on $SP: check/m-checked.txt goes in every other species' kit"
fi
[ "$SP" != "$CHECK_SPECIES" ] || [ -e "$OUT/flags/m_checked.done" ] || { say "the check species has not passed its checks"; false; }

# The fits of M.
m_fits "$OUT/expect.tsv" "$M_DATASETS" "$M_A" m
same_counts "$OUT/expect.tsv" m

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

# The two readings. Their output stays with their tables in score/, which no one opens before every species has run.
STEP=score
python3 -I -B "$HERE/score.py" "$OUT/expect.tsv" "$PICKS" "$R" "$OUT/score" > "$OUT/score/score.log" 2>&1 \
  || { scrub "$OUT/score/score.log"; say "score: stopped; its output is in score/score.log"; false; }
scrub "$OUT/score/score.log"
STEP=second
mkdir -p "$OUT/score/second"
python3 -I -B "$HERE/second.py" "$OUT/expect.tsv" "$PICKS" "$SCORES" "$STUDY" "$R" "$FLOOR" "$LIST" "$OUT/score/second" \
  > "$OUT/score/second/second.log" 2>&1 \
  || { scrub "$OUT/score/second/second.log"; say "second: stopped; its output is in score/second/second.log"; false; }
scrub "$OUT/score/second/second.log"
say "score: both readings written"

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
{ find results check score -type f \( -name '*.json' -o -name '*.tsv' -o -name '*.txt' -o -name '*.log' \)
  for f in binary.txt binary-first.txt expect.tsv expect-first.tsv lists.txt logs/timing.tsv; do [ ! -e "$f" ] || echo "$f"; done
  find logs -name 'inputs-*.sha256'; } | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
python3 -I -B "$HERE/pack.py" check "$OUT" "$R/species/$SP.json"
pack=$(python3 -I -B "$HERE/pack.py" tar "$OUT" "$R/species/$SP.json" "$OUT/pack.tar" | sed -n 's/.*sha256 //p')
[ -n "$pack" ] || { say "pack.py wrote no pack"; false; }
disk done
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $CODE results $(find results/picks -name '*.pick.json' | wc -l) pack $pack"
echo "$line" > DONE
say "$line"
