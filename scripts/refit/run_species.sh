#!/usr/bin/env bash
# run_species.sh KIT SPECIES DATASETS OUT [JOBS]
#
# The refit of the release's rule that saves fold models, for one species,
# then plan diagnostic D0 and block C8 on it (docs/assumption-refit.md), on a
# fresh box with no agent attached: no prompt, nothing read from another
# machine but the public sources below, one log, and a restart that skips
# every finished step. Run it from a clean checkout of this branch:
#
#   setsid nohup bash scripts/refit/run_species.sh KIT SPECIES DATASETS OUT 16 > /dev/null 2>&1 < /dev/null &
#
# DATASETS is a comma list of the species' dataset ids that are not marked M.
# KIT holds these files, checked against the sha256 pinned here before
# anything else runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   MASTER-MANIFEST.sha256      the data root's list of files and sha256
#   curated.tar, curated.sha256 the rebuilt revision 14 table's curated files
#   amrkit.tar, tool-tables.sha256
#                               the AMRFinderPlus kit and the sha256 of each
#                               of the program's per-genome tables it made
#                               for the release
#   marks.tsv                   the label count's table: a dataset is M when
#                               pairs_S+I or pairs_I+R is above 0
#
# The steps up to the matrix are scripts/c6-mic/run_species.sh's, in the same
# text. Then:
#   learn   grm learn at the release's seed and settings over the datasets
#           whose report is missing or has no fold_models, into
#           ROOT/species/SPECIES/results; every listed report must then hold
#           them;
#   check   refit_check.py: each report equals the release's in its input
#           fingerprint and every fitted field, with fold_models complete;
#           any difference stops the run before D0 or C8 reads a report;
#   C8      c8.py over the species' refit reports;
#   D0      for the species' datasets in the plan's MISS set that are listed
#           (scripts/assumption-d0-c8/miss_sets.py): the program installed
#           from the kit and proved on the tool's own test file, run on the
#           species' genomes that have a table recorded for the release,
#           each table equal to its recorded sha256, then per dataset grm
#           export on the same data root and d0.py dataset on the refit
#           report; a species with none skips D0 and says so in its DONE
#           line.
#
# OUT gets logs/ (each step's log with genome ids written as [genome id]; the
# raw logs stay in work/raw-logs), flags/, work/ (release, root with the
# refit reports, build, program, exports), binary.txt, expect.txt, d0.txt,
# check/ (counts and field names), results/d0 and results/c8 (aggregates:
# counts, rates and dataset names), and at the end SHA256SUMS and DONE:
#   DONE <species> <time> binary <sha256> code <sha256> d0 <n> datasets SHA256SUMS <sha256>
# or `d0 skipped: no MISS dataset among the listed` in place of `d0 <n> datasets`.
# work/ holds genome ids and stays on the box. Any failure writes OUT/FAILED
# naming the step and exits 1. It stops nothing but itself.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; DATASETS=${3:?dataset ids}
OUT=${4:?out folder}; JOBS=${5:-$(nproc)}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/check "$OUT"/results; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/refit; DC=$REPO/scripts/assumption-d0-c8
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; M=$KIT/MASTER-MANIFEST.sha256
GRM=$W/target/release/grm; AMR=$W/amr; RES=$W/root/species/$SP/results
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
export KIT SP DATASETS OUT JOBS REPO HERE DC W ROOT R M GRM AMR RES

STEP=kit
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the binary must be its commit"; false; }
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/curated.tar" d77e2c6f75f8d986268d0c075cb2ea364fba9972fa22a55139495ac0eb376d97
pinned "$KIT/curated.sha256" 8b0e43dd98ca1b3fc7609b7e25b1b52d4bb4dd0c4c555cf5e771622d24fd9c6e
pinned "$KIT/amrkit.tar" da07f277a7d271f1ae006919d2c2df4185737526348302f665a9bcb3a083d0e6
pinned "$KIT/tool-tables.sha256" 9291b83f4c7d01a0384f35446e6de8f986dff5f659439ee3d673bccdd09427ba
pinned "$KIT/marks.tsv" ba662467864edce9e8cc617651b81598feac5f32bb1219d7409467ff094aa313
echo "$DATASETS" | tr ',' '\n' | grep -v '^$' | LC_ALL=C sort -u > "$OUT/expect.txt.part"
grep -q . "$OUT/expect.txt.part" || { say "no dataset listed"; false; }
if grep -vq "^${SP}__" "$OUT/expect.txt.part"; then say "a listed dataset is not of $SP"; false; fi
# Only datasets that are not M are fitted (docs/assumption-plan.md, Rules, Labels).
python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    marks = {r["dataset"]: int(r["pairs_S+I"]) > 0 or int(r["pairs_I+R"]) > 0 for r in csv.DictReader(f, delimiter="\t")}
listed = open(sys.argv[2]).read().split()
gone = [d for d in listed if d not in marks]
m = [d for d in listed if marks.get(d)]
if gone or m:
    print(f"{len(gone)} listed datasets are not in marks.tsv, {len(m)} are marked M: " + ",".join(gone + m))
    sys.exit(1)
' "$KIT/marks.tsv" "$OUT/expect.txt.part" || { say "the list must hold datasets that are not M only"; false; }
if [ -s "$OUT/expect.txt" ] && ! cmp -s "$OUT/expect.txt" "$OUT/expect.txt.part"; then
  say "the dataset list differs from the first start's; use a fresh OUT"; false
fi
mv "$OUT/expect.txt.part" "$OUT/expect.txt"
# D0's datasets: the plan's MISS set among the listed ones.
python3 -I -B "$DC/miss_sets.py" --species "$SP" --list miss | LC_ALL=C sort -u | LC_ALL=C comm -12 - "$OUT/expect.txt" > "$OUT/d0.txt.part"
if [ -e "$OUT/d0.txt" ] && ! cmp -s "$OUT/d0.txt" "$OUT/d0.txt.part"; then
  say "D0's dataset list differs from the first start's; use a fresh OUT"; false
fi
mv "$OUT/d0.txt.part" "$OUT/d0.txt"
say "expect: $(wc -l < "$OUT/expect.txt") reports, D0 on $(wc -l < "$OUT/d0.txt")"

release() {
  rm -rf "$W/release"; mkdir -p "$W/release"
  tar --zstd -xf "$KIT/rev14-274-20261112.tar.zst" -C "$W/release"
  test -f "$R/species/$SP.json"
  while read -r d; do test -f "$R/datasets/$d.json" || { echo "the release has no report $d"; exit 1; }; done < "$OUT/expect.txt"
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
  mkdir -p "$ROOT/curated" "$ROOT/bvbrc" "$W/specs"
  tar -xf "$KIT/curated.tar" -C "$ROOT/curated"
  (cd "$ROOT/curated" && sha256sum --quiet -c "$KIT/curated.sha256")
  # Genome metadata feeds only a sample's country, year, host and mlst, and
  # is not in the input fingerprint. The check reports those four apart.
  : > "$ROOT/bvbrc/genomes.jsonl"
  # The release's breakpoint for every listed dataset (c6's release_specs.py),
  # so the fingerprint can equal the release's.
  python3 -I -B "$REPO/scripts/c6-mic/release_specs.py" "$ROOT/curated/datasets.json" "$R" \
    "$(paste -sd, "$OUT/expect.txt")" "$W/specs"
  cp "$W/specs/release-specs.tsv" "$OUT/check/release-specs.tsv"
  cp "$ROOT/curated/datasets.json" "$W/specs/curated-datasets.json"
  cp "$W/specs/datasets.json" "$ROOT/curated/datasets.json"
  echo "breakpoints set to the release's: $(tr ',' '\n' < "$W/specs/differ.txt" | grep -c . || true)"
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

missing() {  # listed datasets whose report is missing or has no fold_models
  python3 -I -c '
import json, os, sys
for d in open(sys.argv[1]).read().split():
    p = os.path.join(sys.argv[2], d + ".json")
    ok = os.path.isfile(p)
    if ok:
        with open(p) as f:
            ok = "fold_models" in json.load(f)
    if not ok:
        print(d)
' "$OUT/expect.txt" "$RES"
}
STEP=learn
todo=$(missing | paste -sd, -)
if [ -z "$todo" ]; then
  say "learn: every report present with fold_models, skipped"
else
  say "learn: start, $todo"
  # With the network: a feature card's gene and product come from the
  # BV-BRC CDS list of its genome, cached in ROOT/bvbrc/features.
  /usr/bin/time -v -o "$OUT/logs/learn.time" "$GRM" --data "$ROOT" learn --species "$SP" --only "$todo" --force \
    < /dev/null > "$OUT/logs/learn.log" 2>&1 || { scrub "$OUT/logs/learn.log"; tail -n 30 "$OUT/logs/learn.log"; false; }
  scrub "$OUT/logs/learn.log"
  say "learn: end, $(grep -E 'Elapsed|Maximum resident' "$OUT/logs/learn.time" | sed 's/^\s*//' | tr '\n' ';')"
  left=$(missing | paste -sd, -)
  [ -z "$left" ] || { say "learn: no report with fold_models for $left"; false; }
fi
n=0; same=0
for f in "$ROOT"/bvbrc/features/*.json; do
  [ -e "$f" ] || continue
  n=$((n + 1)); master "bvbrc/features/$(basename "$f")" "$f" && same=$((same + 1))
done
say "feature lists fetched: $n, equal to the master list: $same"
echo -e "fetched\tequal_to_master\n$n\t$same" > "$OUT/check/feature-lists.tsv"
disk learn

# The gate: no report reaches D0 or C8 unless it is the release's rule.
STEP=check
say "check: start"
python3 -I -B "$HERE/refit_check.py" "$RES" "$R" "$(paste -sd, "$OUT/expect.txt")" "$OUT/check" > "$OUT/logs/check.log" 2>&1 \
  || { cat "$OUT/logs/check.log"; false; }
cat "$OUT/logs/check.log"

# C8 reads the checked refit reports only, so it runs before the program
# and D0: a failure there (the network, a tool table) still leaves C8's
# tables in results/c8 for the copy back.
c8() {
  rm -rf "$OUT/results/c8"
  python3 -I -B -W error "$DC/c8.py" --reports "$RES" --only "$OUT/expect.txt" --marks "$KIT/marks.tsv" --out "$OUT/results/c8"
}
step c8 c8

if [ -s "$OUT/d0.txt" ]; then
  # The program, from the kit, proved on the tool's own test file.
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
    # The kit's fixture genomes are made from the tool's own bundled test file.
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
  step program_env program_env

  # The program on every genome of the species that has a recorded table; each
  # table must equal the sha256 recorded for the release, so the calls are
  # the release's. The release ran the program on the genomes of its rule
  # datasets only (the kit's RERUN-SPEC.md, 23,751 genomes), so a genome of
  # the species list that is in no rule dataset has no recorded table and is
  # not run: Staphylococcus aureus lists 2,000 genomes and 1,994 have one.
  # Each recorded genome must be in the species list, and the adapter stops
  # on a D0 report genome with no table.
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
    # The calls, read on the refit reports of D0's datasets.
    rm -rf "$AMR/reports" && mkdir -p "$AMR/reports"
    while read -r d; do ln -s "$RES/$d.json" "$AMR/reports/$d.json"; done < "$OUT/d0.txt"
    python3 -I -B scripts/amrfinder_adapter.py --reports "$AMR/reports" --tool-dir "$AMR/tool" \
      --run-manifest spec/run-manifest.rev14.json --calls "$AMR/calls.tsv" --summary "$AMR/adapter-summary.json"
  }
  step program_run program_run

  # D0 per dataset: the export of the same data root, read beside the refit
  # report's fold models. An export holds the dataset's matrix and is
  # removed once read.
  while read -r id; do
    d0_one() {
      rm -rf "$W/export/$id" "$W/d0/$id"
      unshare -rn "$GRM" --data "$ROOT" export --species "$SP" --only "$id" --out "$W/export"
      python3 -I -B -W error "$DC/d0.py" dataset --report "$RES/$id.json" --export "$W/export/$id" \
        --calls "$AMR/calls.tsv" --reference "$REPO/reference/amrfinderplus-comparison/callers.tsv" --out "$W/d0/$id"
      rm -rf "$W/export/$id"
    }
    export id
    step "d0-$id" d0_one
  done < "$OUT/d0.txt"

  d0_results() {
    dest=$OUT/results/d0; rm -rf "$dest"; mkdir -p "$dest"
    dirs=()
    while read -r id; do dirs+=("$W/d0/$id"); done < "$OUT/d0.txt"
    # The miss set's datasets and pairs among D0's datasets, from the record,
    # so a dataset table left out or read twice stops the step.
    pairs=$(python3 -I -c '
import sys
sys.path.insert(0, sys.argv[1])
import miss_sets
want = set(open(sys.argv[2]).read().split())
print(sum(d in want for v in miss_sets.sets()[1].values() for d in v))
' "$DC" "$OUT/d0.txt")
    python3 -I -B -W error "$DC/d0.py" combine --marks "$KIT/marks.tsv" \
      --expect-datasets "$(wc -l < "$OUT/d0.txt")" --expect-pairs "$pairs" --out "$dest" "${dirs[@]}"
    python3 -I -c '
import csv, sys
out = sys.argv[1]
for name in ("d0-datasets.tsv", "d0-folds.tsv"):
    rows, header = [], None
    for d in sys.argv[2:]:
        with open(d + "/" + name, encoding="utf-8", newline="") as f:
            r = list(csv.reader(f, delimiter="\t"))
        header = header or r[0]
        assert r[0] == header, name
        rows += r[1:]
    with open(out + "/" + name, "w", encoding="utf-8", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows([header] + rows)
' "$dest" "${dirs[@]}"
    cp "$AMR/adapter-summary.json" "$dest/program-adapter-summary.json"
  }
  step d0-results d0_results
  d0line="d0 $(wc -l < "$OUT/d0.txt") datasets"
else
  say "D0 skipped: no dataset of the plan's MISS set is listed for $SP"
  d0line="d0 skipped: no MISS dataset among the listed"
fi

STEP=done
cd "$OUT"
{ find check results -type f \( -name '*.tsv' -o -name '*.txt' -o -name '*.json' \); echo binary.txt; echo expect.txt; echo d0.txt; } \
  | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
disk done
line="DONE $SP $(date -u +%FT%TZ) binary $(sed -n 's/^grm //p' binary.txt) code $(sed -n 's/^code //p' binary.txt) $d0line SHA256SUMS $(sha256sum SHA256SUMS | cut -c1-64)"
echo "$line" > DONE
say "$line"
