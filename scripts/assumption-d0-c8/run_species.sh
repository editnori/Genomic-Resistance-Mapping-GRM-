#!/usr/bin/env bash
# run_species.sh KIT SPECIES OUT [JOBS]
#
# D0 for one species of the miss set on a fresh box, with no agent attached:
# no prompt, nothing read from the controller, one log, and a restart that
# skips every finished step. Run it from a checkout of this branch:
#
#   setsid nohup bash scripts/assumption-d0-c8/run_species.sh KIT SPECIES OUT 16 > /dev/null 2>&1 < /dev/null &
#
# KIT holds the five files below, checked against the sha256 pinned here
# before anything else runs. OUT gets logs/run.log, flags/<step>.done, the
# work folders and results/<species>/. The last step writes OUT/DONE: a
# `DONE` line with the sha256 of results/<species>/SHA256SUMS. Any failed
# check stops the script with rc 1 and writes OUT/FAILED naming the step.
#
# What it reads: the release's reports and species lists, the species'
# public genomes (fetched from BV-BRC by grm's own step, or copied from
# DATA/fasta when DATA names a data root holding them), the pinned
# AMRFinderPlus kit, tool, database and bundled test file, and the
# repository's comparison record. Nothing is fitted or drawn.
#
# The released reports hold no fold models, so D0 here reads the causes
# Absent and Removed and leaves the rest open (scripts/assumption-d0-c8/d0.py,
# `release`). With EXPORT=1 and DATA naming a data root with the release's
# curated tables, it exports each dataset and also recomputes boosting's and
# the lasso's screens from the plan's holdouts (`dataset`); that path was not
# run in the pilot.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; OUT=${3:?out folder}; JOBS=${4:-$(nproc)}
DATA=${DATA:-}; EXPORT=${EXPORT:-0}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work; OUT=$(cd "$OUT" && pwd)
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
W=$OUT/work; ROOT=$W/root; R=$W/release/rev14-274-20261112; AMR=$W/amr
M=$KIT/MASTER-MANIFEST.sha256
export RUSTUP_TOOLCHAIN=1.98.1 CARGO_BUILD_JOBS=$JOBS RAYON_NUM_THREADS=$JOBS PATH=$HOME/.cargo/bin:$PATH
STEP=start
say() { echo "$(date -u +%FT%TZ) $*"; }
trap 'say "FAILED at $STEP"; echo "FAILED $SP at $STEP $(date -u +%FT%TZ)" > "$OUT/FAILED"; exit 1' ERR
trap 'say "stopped at $STEP"; echo "FAILED $SP stopped at $STEP $(date -u +%FT%TZ)" > "$OUT/FAILED"; exit 1' TERM INT
rm -f "$OUT/FAILED"
say "start $SP, $JOBS jobs, repo $(git -C "$REPO" rev-parse HEAD)"

pinned() {  # pinned <file> <sha256>
  [ "$(sha256sum "$1" | cut -c1-64)" = "$2" ] || { say "$1: sha256 differs from the pinned value"; return 1; }
}
master() {  # master <path in the data root> <file>
  local want; want=$(grep -m1 "  $1\$" "$M" | cut -c1-64)
  [ -n "$want" ] && [ "$(sha256sum "$2" | cut -c1-64)" = "$want" ]
}
step() {  # step <name> <function>: once, timed, in a shell that stops at the first failure
  STEP=$1
  if [ -e "$OUT/flags/$1.done" ]; then say "$1: done before, skipped"; return 0; fi
  say "$1: start"
  /usr/bin/time -f "$1 %e s wall, %M KB peak, %U s user, %S s sys" -o "$OUT/logs/$1.time" \
    bash -euo pipefail -c "$(declare -f say pinned master "$2"); $2" < /dev/null > "$OUT/logs/$1.log" 2>&1 \
    || { tail -n 30 "$OUT/logs/$1.log"; false; }
  touch "$OUT/flags/$1.done"
  say "$1: end, $(cat "$OUT/logs/$1.time")"
}
export KIT SP OUT JOBS DATA EXPORT REPO W ROOT R AMR M
STEP=kit
pinned "$KIT/repo.bundle" 4ef909f09cd5b28f7566379509326c1dbf4ead600a2b8eb72d0b828e47c25337
pinned "$KIT/amrkit.tar" da07f277a7d271f1ae006919d2c2df4185737526348302f665a9bcb3a083d0e6
pinned "$M" 7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887
pinned "$KIT/tool-tables.sha256" 9291b83f4c7d01a0384f35446e6de8f986dff5f659439ee3d673bccdd09427ba
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
git -C "$REPO" merge-base --is-ancestor 934292c HEAD || { say "the checkout does not hold main at 934292c"; false; }
python3 -I -B "$REPO/scripts/assumption-d0-c8/miss_sets.py" --species "$SP" --list miss > "$OUT/logs/miss.txt"
grep -q . "$OUT/logs/miss.txt" || { say "$SP has no dataset in the miss set"; false; }

release() {
  rm -rf "$W/release"; mkdir -p "$W/release"
  tar --zstd -xf "$KIT/rev14-274-20261112.tar.zst" -C "$W/release"
  test -f "$R/species/$SP.json"
  ls "$R"/datasets/"$SP"__*.json | wc -l
}
step release release

build() {
  cd "$REPO" && CARGO_TARGET_DIR=$W/target cargo build --release --locked
  sha256sum "$W/target/release/grm"
}
step build build
GRM=$W/target/release/grm; export GRM

genomes() {
  mkdir -p "$ROOT/species/$SP" "$ROOT/fasta"
  python3 -I -c 'import json,sys; [print(g["id"]) for g in json.load(open(sys.argv[1]))["genomes"]]' \
    "$R/species/$SP.json" > "$ROOT/species/$SP/genomes.txt"
  master "species/$SP/genomes.txt" "$ROOT/species/$SP/genomes.txt"
  if [ -n "$DATA" ] && [ -d "$DATA/fasta" ]; then
    while read -r id; do cp -n "$DATA/fasta/$id.fna.gz" "$ROOT/fasta/"; done < "$ROOT/species/$SP/genomes.txt"
  fi
  "$GRM" --data "$ROOT" fetch fasta --species "$SP" --workers 4
  n=0; bad=0
  while read -r id; do
    n=$((n + 1)); master "fasta/$id.fna.gz" "$ROOT/fasta/$id.fna.gz" || bad=$((bad + 1))
  done < "$ROOT/species/$SP/genomes.txt"
  # A refetched file can differ in its bytes from the master list (28 of 93
  # did in the pilot). Every genome must still have the release's length,
  # contig count and N50, and the program's tables are checked after it runs.
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
  echo "install proved: the fixture's 13 tool tables equal the kit's record"
}
step program_env program_env

program_run() {
  export PATH="$AMR/env/bin:$PATH"; cd "$AMR/kit"
  { printf 'species_slug\tgenome_id\n'; sed "s/^/$SP\t/" "$ROOT/species/$SP/genomes.txt"; } > "$AMR/cohort.tsv"
  bash scripts/run_amrfinder.sh "$AMR/cohort.tsv" "$ROOT/fasta" "$AMR/tool" spec/run-manifest.rev14.json \
    "$AMR/db/2026-08-07.1" "$JOBS" 1 4000
  grep "  $SP/" "$KIT/tool-tables.sha256" | sort -k2 > "$AMR/tables.recorded"
  ( cd "$AMR/tool" && find "$SP" -name '*.tsv' | sort | xargs sha256sum ) | sort -k2 > "$AMR/tables.now"
  echo "tool tables: $(wc -l < "$AMR/tables.now") made, $(wc -l < "$AMR/tables.recorded") recorded"
  diff "$AMR/tables.recorded" "$AMR/tables.now"
  echo "every tool table equals its recorded sha256"
  rm -rf "$AMR/reports" && mkdir -p "$AMR/reports"
  for f in "$R"/datasets/"$SP"__*.json; do ln -s "$f" "$AMR/reports/"; done
  python3 -I -B scripts/amrfinder_adapter.py --reports "$AMR/reports" --tool-dir "$AMR/tool" \
    --run-manifest spec/run-manifest.rev14.json --calls "$AMR/calls.tsv" --summary "$AMR/adapter-summary.json"
  python3 -I -B scripts/matched_eval.py --reports "$AMR/reports" --calls "$AMR/calls.tsv" \
    --by-dataset "$AMR/by-dataset.tsv" > "$AMR/matched-eval.json"
}
step program_run program_run

catalog() {
  mkdir -p "$ROOT/catalog"
  for f in AMR_CDS.fa ReferenceGeneCatalog.txt version.txt; do
    cp "$AMR/db/2026-08-07.1/$f" "$ROOT/catalog/$f"
    master "catalog/$f" "$ROOT/catalog/$f"
  done
  echo "catalog: 3 of 3 files equal the master list"
}
step catalog catalog

count() { unshare -rn "$GRM" --data "$ROOT" count --species "$SP"; }
structure() { unshare -rn "$GRM" --data "$ROOT" structure --species "$SP"; }
matrix() { unshare -rn "$GRM" --data "$ROOT" matrix --species "$SP"; }
step count count
step structure structure
step matrix matrix

curated() {
  # Only with EXPORT=1: the release's curated tables and record files, each
  # equal to the master list.
  for f in curated/datasets.json curated/drug_names.json curated/funnel.json curated/phenotypes.jsonl \
           curated/qc.json curated/species.json bvbrc/genomes.jsonl bvbrc/genome_amr.jsonl \
           species-check/species_check.json; do
    master "$f" "$DATA/$f"
    mkdir -p "$(dirname "$ROOT/$f")"; cp "$DATA/$f" "$ROOT/$f"
  done
  echo "curated tables: 9 of 9 equal the master list"
}
if [ "$EXPORT" = 1 ]; then
  [ -n "$DATA" ] || { say "EXPORT=1 needs DATA, a data root with the curated tables"; false; }
  step curated curated
fi

D0=$REPO/scripts/assumption-d0-c8/d0.py
REF=$REPO/reference/amrfinderplus-comparison/callers.tsv
export D0 REF
while read -r id; do
  d0_one() {
    rep=$R/datasets/$id.json
    if [ "$EXPORT" = 1 ]; then
      rm -rf "$W/export/$id"
      unshare -rn "$GRM" --data "$ROOT" export --species "$SP" --only "$id" --out "$W/export"
      python3 -I -B -W error "$D0" dataset --report "$rep" --export "$W/export/$id" --calls "$AMR/calls.tsv" \
        --reference "$REF" --out "$W/d0/$id"
      rm -rf "$W/export/$id"
    else
      python3 -I -B -W error "$D0" release --report "$rep" --species-matrix "$ROOT/species/$SP/matrix" \
        --calls "$AMR/calls.tsv" --reference "$REF" --out "$W/d0/$id"
    fi
  }
  export id
  step "d0-$id" d0_one
done < "$OUT/logs/miss.txt"

results() {
  dest=$OUT/results/$SP; rm -rf "$dest"; mkdir -p "$dest"
  dirs=()
  while read -r id; do dirs+=("$W/d0/$id"); done < "$OUT/logs/miss.txt"
  python3 -I -B -W error "$D0" combine --out "$dest" "${dirs[@]}"
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
  cat "$OUT"/logs/*.time > "$dest/cost.txt"
  for d in release root amr target; do echo "disk $d $(du -sb "$W/$d" | cut -f1) bytes"; done >> "$dest/cost.txt"
  ( cd "$dest" && sha256sum d0-*.tsv d0-ordering.txt program-adapter-summary.json cost.txt > SHA256SUMS )
}
step results results
line="DONE $SP $(date -u +%FT%TZ) results/$SP/SHA256SUMS $(sha256sum "$OUT/results/$SP/SHA256SUMS" | cut -c1-64)"
echo "$line" > "$OUT/DONE"
say "$line"
