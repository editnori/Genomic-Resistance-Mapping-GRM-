#!/usr/bin/env bash
# run_species.sh KIT SPECIES OUT META [TABLES] [DB] [FOLDS] [FASTA] [JOBS]
#
# The mechanism split of one species (docs/mechanism-split.md): block M1,
# recall by mechanism (strata.py), and block M2, the mechanism mix over
# collection years (history.py), on a box with no agent attached: no
# prompt, one log, and a restart that skips every finished step. Run it from
# a clean checkout of this branch:
#
#   setsid nohup bash scripts/mechanism/run_species.sh KIT SPECIES OUT META TABLES DB > /dev/null 2>&1 < /dev/null &
#
# KIT holds these files, held to the sha256 pinned here before anything runs:
#   rev14-274-20261112.tar.zst  the public release data-rev14-274-20261112
#   amrkit.tar                  the AMRFinderPlus kit (AMRFinderPlus 4.2.7,
#                               database 2026-08-07.1, the adapter)
#   tool-tables.sha256          the sha256 of each program table the release
#                               recorded, as <species slug>/<genome id>.tsv
# META is BV-BRC's genome table the study and time splits read
# (master-genome-meta.jsonl.gz), held to its sha256 in history.py.
# TABLES (or -) is a folder of the program's tables made before, as
# TABLES/SPECIES/<genome id>.tsv: each table of the species recorded in
# tool-tables.sha256 is copied from it and held to its recorded sha256, and
# the program is neither installed nor run. With TABLES -, the program's
# environment is installed from the kit, proved on the tool's own test file,
# and run on every genome of the species with a recorded table, read from
# FASTA (a folder as the kit's run_amrfinder.sh reads it); every table made
# must equal its recorded sha256.
# DB (or -) is a copy of the database 2026-08-07.1, held to the kit's
# spec/db-2026-08-07.1.sha256; with -, the kit's fetch_database.sh fetches
# it from NCBI and holds it to the same list.
# FOLDS (or -) is the time split's results folder (design.sha256 and
# fold-files/); with it, history.py writes the time-split link for the
# species' datasets that have a past arm.
#
# Steps, each once: kit (pins and tools), release (the species' reports),
# amrkit, db, tables (copied, or the program run), adapter (the kit's
# amrfinder_adapter.py, unchanged, on the species' reports of A), rules
# (each report's catalog.rule equal to scripts/mechanism/drug_rules.tsv),
# strata, history, scan, end. strata.py stops unless its own program call
# equals the adapter's on every genome.
#
# OUT gets logs/run.log, logs/<step>.log and .time (genome ids written as
# [genome id]; the raw logs stay in work/raw-logs), logs/timing.tsv (per
# step: wall seconds and peak memory), flags/, work/ (the release, the kit,
# the database, the tables, the adapter's calls: these name genomes and stay
# on the box), and results/m1, results/m2 (aggregate tables, no genome id).
# The scan step fails when any file under results/ or logs/ names a genome
# of the species or a BioProject accession, holds a sequence or a position,
# writes a rate under its floor or names a catalogue family or gene carried
# by under 5 genomes. The end writes
# results/SHA256SUMS and OUT/DONE:
#   DONE <species> <time> commit <sha> results <n> SHA256SUMS <sha256>
# Any failure writes OUT/FAILED naming the step and exits 1.
set -Eeuo pipefail
KIT=$(cd "${1:?kit folder}" && pwd); SP=${2:?species slug}; OUT=${3:?out folder}; META=${4:?BV-BRC genome table}
TABLES_FROM=${5:--}; DB_FROM=${6:--}; FOLDS=${7:--}; FASTA=${8:--}; JOBS=${9:-$(nproc)}
mkdir -p "$OUT"/flags "$OUT"/logs "$OUT"/work "$OUT"/results; OUT=$(cd "$OUT" && pwd)
META=$(cd "$(dirname "$META")" && pwd)/$(basename "$META")
exec >> "$OUT/logs/run.log" 2>&1
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HERE=$REPO/scripts/mechanism
W=$OUT/work; R=$W/release/rev14-274-20261112; AMR=$W/amr; TABLES=$AMR/tool/$SP
DB=$AMR/db/2026-08-07.1; RA=$W/reports-a; RS=$W/reports-species
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
[ -s "$OUT/logs/timing.tsv" ] || printf 'step\twall_seconds\tpeak_kb\n' > "$OUT/logs/timing.tsv"

pinned() {  # pinned <file> <sha256>
  [ "$(sha256sum "$1" | cut -c1-64)" = "$2" ] || { say "$1: sha256 differs from the pinned value"; return 1; }
}
# Genome ids stay on the box: each step's log moves to work/raw-logs and
# logs/ gets it back with every genome id of the species written as
# [genome id].
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
  /usr/bin/time -f "%e %M" -o "$OUT/logs/$1.time" \
    bash -euo pipefail -c "$(declare -f say pinned "$2"); $2" < /dev/null > "$OUT/logs/$1.log" 2>&1 \
    || { scrub "$OUT/logs/$1.log"; tail -n 30 "$OUT/logs/$1.log"; false; }
  scrub "$OUT/logs/$1.log"
  read -r wall peak < <(tail -n 1 "$OUT/logs/$1.time")
  printf '%s\t%s\t%s\n' "$1" "$wall" "$peak" >> "$OUT/logs/timing.tsv"
  touch "$OUT/flags/$1.done"
  say "$1: end, $wall s wall, $peak KB peak; $(tail -n 1 "$OUT/logs/$1.log")"
}
export KIT SP OUT META TABLES_FROM DB_FROM FOLDS FASTA JOBS REPO HERE W R AMR TABLES DB RA RS

STEP=kit
[ -z "$(git -C "$REPO" status --porcelain)" ] || { say "the checkout has changes; the run must be its commit"; false; }
for tool in python3 tar zstd sha256sum flock /usr/bin/time; do
  command -v "$tool" > /dev/null || { say "$tool is not installed, which a later step needs"; false; }
done
pinned "$KIT/rev14-274-20261112.tar.zst" c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de
pinned "$KIT/amrkit.tar" da07f277a7d271f1ae006919d2c2df4185737526348302f665a9bcb3a083d0e6
pinned "$KIT/tool-tables.sha256" 9291b83f4c7d01a0384f35446e6de8f986dff5f659439ee3d673bccdd09427ba
pinned "$META" 84ae57a777cb234aa0a39102ab7f9ca9448e52d1474654179f26165644c25f64
grep -qx "$SP" <(cut -f2 "$HERE/../c3-inputs/lists.tsv") || { say "$SP has no dataset in scripts/c3-inputs/lists.tsv"; false; }
if [ "$TABLES_FROM" != - ]; then
  [ -d "$TABLES_FROM/$SP" ] || { say "TABLES has no folder $SP"; false; }
  TABLES_FROM=$(cd "$TABLES_FROM" && pwd)
elif [ "$FASTA" = - ]; then
  say "with no TABLES the program runs, and it needs FASTA"; false
fi
[ "$DB_FROM" = - ] || DB_FROM=$(cd "$DB_FROM" && pwd)
if [ "$FOLDS" != - ]; then
  [ -f "$FOLDS/design.sha256" ] && [ -d "$FOLDS/fold-files" ] || { say "FOLDS needs design.sha256 and fold-files/"; false; }
  FOLDS=$(cd "$FOLDS" && pwd)
fi
export TABLES_FROM DB_FROM FOLDS
rm -f "$OUT/FAILED" "$OUT/DONE"

release() {
  rm -rf "$W/release" "$RA" "$RS"; mkdir -p "$W/release" "$RA" "$RS"
  tar --zstd -xf "$KIT/rev14-274-20261112.tar.zst" -C "$W/release"
  test -f "$R/species/$SP.json"
  for f in "$R"/datasets/"$SP"__*.json; do ln -s "$f" "$RS/$(basename "$f")"; done
  python3 -I -c '
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["species_slug"] == sys.argv[2] and r["in_A"] == "1":
            print(r["dataset"])
' "$HERE/../c3-inputs/lists.tsv" "$SP" > "$W/a.txt"
  while read -r d; do
    test -f "$R/datasets/$d.json" || { echo "the release has no report $d"; exit 1; }
    ln -s "$R/datasets/$d.json" "$RA/$d.json"
  done < "$W/a.txt"
  echo "reports of $SP in the release: $(ls "$RS" | wc -l); of them in A: $(wc -l < "$W/a.txt")"
}
step release release

amrkit() {
  rm -rf "$AMR/kit"; mkdir -p "$AMR/kit"
  tar -xf "$KIT/amrkit.tar" -C "$AMR/kit"
  echo "kit adapter sha256 $(sha256sum "$AMR/kit/scripts/amrfinder_adapter.py" | cut -c1-64)"
  python3 -I -c '
import sys; sys.path.insert(0, sys.argv[1]); import amr_rule
want = {"amrfinder_adapter.py": amr_rule.ADAPTER_SHA256, "grm_eval.py": amr_rule.GRM_EVAL_SHA256}
import hashlib, os
for name, h in want.items():
    got = hashlib.sha256(open(os.path.join(sys.argv[2], name), "rb").read()).hexdigest()
    if got != h:
        sys.exit(f"{name}: the kit holds another version than amr_rule.py ports")
print("the kit holds the adapter and grm_eval amr_rule.py ports")
' "$HERE" "$AMR/kit/scripts"
}
step amrkit amrkit

db() {
  rm -rf "$DB"; mkdir -p "$(dirname "$DB")"
  if [ "$DB_FROM" != - ]; then
    cp -r "$DB_FROM" "$DB"
  else
    ( cd "$AMR/kit" && bash scripts/fetch_database.sh 2026-08-07.1 spec/db-2026-08-07.1.sha256 "$DB" )
  fi
  ( cd "$DB" && sha256sum --check --quiet "$AMR/kit/spec/db-2026-08-07.1.sha256" )
  echo "database 2026-08-07.1: $(wc -l < "$AMR/kit/spec/db-2026-08-07.1.sha256") files equal the kit's list"
}
step db db

tables() {
  grep "  $SP/" "$KIT/tool-tables.sha256" | LC_ALL=C sort -k2 > "$AMR/tables.recorded"
  [ -s "$AMR/tables.recorded" ] || { echo "no table of $SP is recorded"; exit 1; }
  rm -rf "$AMR/tool"; mkdir -p "$TABLES"
  if [ "$TABLES_FROM" != - ]; then
    missing=0
    while read -r _ path; do
      if [ -f "$TABLES_FROM/$path" ]; then cp "$TABLES_FROM/$path" "$AMR/tool/$path"; else missing=$((missing + 1)); fi
    done < "$AMR/tables.recorded"
    [ "$missing" -eq 0 ] || { echo "tables: $missing recorded tables are not in TABLES"; exit 1; }
    echo "tables copied from TABLES: $(wc -l < "$AMR/tables.recorded")"
  else
    export PATH="$AMR/env/bin:$PATH"; cd "$AMR"
    curl -fsSL --retry 3 -A curl -o micromamba \
      "https://github.com/mamba-org/micromamba-releases/releases/download/2.10.0-0/micromamba-linux-64"
    pinned micromamba 701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd
    chmod +x micromamba; rm -rf "$AMR/env"
    MAMBA_ROOT_PREFIX=$AMR/mamba ./micromamba create -y -p "$AMR/env" --file kit/spec/env-amrfinder.explicit.txt
    amrfinder --version | grep -qx '4.2.7'
    cd kit
    python3 -I -B tests/test_contract.py
    python3 -I -B tests/test_runner.py
    mkdir -p "$AMR/testdna"
    curl -fsSL --retry 3 -A curl -o "$AMR/testdna/test_dna.fa" \
      "https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_dna.fa"
    pinned "$AMR/testdna/test_dna.fa" 2c330ed5f979e18ba915ed3aea36c3a0779829e8421d1e1c57c6df842aa2b394
    rm -rf "$AMR/fx" "$AMR/fx-out"
    python3 -I -B scripts/make_fixture.py --test-dna "$AMR/testdna/test_dna.fa" --out "$AMR/fx"
    bash scripts/run_amrfinder.sh "$AMR/fx/cohort.tsv" "$AMR/fx/genomes" "$AMR/fx-out" "$AMR/fx/run-manifest.json" "$DB" 2 1
    diff <(cd "$AMR/fx-out" && sha256sum */*.tsv) <(cd fixtures/e2e/tool-out && sha256sum */*.tsv)
    echo "install proved: the fixture's tool tables equal the kit's record"
    sed -E "s#^[0-9a-f]{64}  $SP/(.+)\.tsv\$#\1#" "$AMR/tables.recorded" > "$AMR/cohort.ids"
    { printf 'species_slug\tgenome_id\n'; sed "s/^/$SP\t/" "$AMR/cohort.ids"; } > "$AMR/cohort.tsv"
    rm -rf "$AMR/tool"
    bash scripts/run_amrfinder.sh "$AMR/cohort.tsv" "$FASTA" "$AMR/tool" spec/run-manifest.rev14.json "$DB" "$JOBS" 1 4000
  fi
  ( cd "$AMR/tool" && find "$SP" -name '*.tsv' | LC_ALL=C sort | xargs sha256sum ) | LC_ALL=C sort -k2 > "$AMR/tables.now"
  diff -q "$AMR/tables.recorded" "$AMR/tables.now" > /dev/null \
    || { echo "tables: $(wc -l < "$AMR/tables.now") here, $(wc -l < "$AMR/tables.recorded") recorded; they differ"; exit 1; }
  echo "tables: $(wc -l < "$AMR/tables.now") of $SP, each equal to its recorded sha256"
}
step tables tables

adapter() {
  cd "$AMR/kit"
  python3 -I -B scripts/amrfinder_adapter.py --reports "$RA" --tool-dir "$AMR/tool" \
    --run-manifest spec/run-manifest.rev14.json --calls "$AMR/calls.tsv" --summary "$AMR/adapter-summary.json"
  echo "adapter: $(($(wc -l < "$AMR/calls.tsv") - 1)) call rows over $(ls "$RA" | wc -l) reports"
}
step adapter adapter

rules() {
  python3 -I "$HERE/drug_rules.py" --reports "$RA" --lists "$HERE/../c3-inputs/lists.tsv" --species "$SP" \
    --check "$HERE/drug_rules.tsv"
}
step rules rules

m1_strata() {
  rm -rf "$OUT/results/m1"
  python3 -I "$HERE/strata.py" --species "$SP" --reports "$RA" --lists "$HERE/../c3-inputs/lists.tsv" \
    --drug-rules "$HERE/drug_rules.tsv" --tables "$TABLES" --calls "$AMR/calls.tsv" \
    --callers "$REPO/reference/amrfinderplus-comparison/callers.tsv" --out "$OUT/results/m1"
}
step strata m1_strata

m2_history() {
  rm -rf "$OUT/results/m2"
  folds=()
  [ "$FOLDS" = - ] || folds=(--folds "$FOLDS/fold-files" --design "$FOLDS/design.sha256")
  python3 -I "$HERE/history.py" --species "$SP" --reports "$RS" --lists "$HERE/../c3-inputs/lists.tsv" \
    --drug-rules "$HERE/drug_rules.tsv" --tables "$TABLES" --db "$DB" --meta "$META" --out "$OUT/results/m2" \
    "${folds[@]}"
}
step history m2_history

STEP=scan
rm -f "$OUT/flags/scan.done"
python3 -I "$HERE/scan.py" --ids "$R/species/$SP.json" "$OUT/results" "$OUT/logs" \
  || { say "a file under results/ or logs/ names a genome of $SP, a BioProject accession, a sequence or a position, a rate under its floor or a name carried by under 5 genomes"; false; }
touch "$OUT/flags/scan.done"

STEP=end
( cd "$OUT/results" && find . -type f ! -name SHA256SUMS | LC_ALL=C sort | xargs sha256sum ) > "$OUT/results/SHA256SUMS"
n=$(wc -l < "$OUT/results/SHA256SUMS")
echo "DONE $SP $(date -u +%FT%TZ) commit $(git -C "$REPO" rev-parse HEAD) results $n SHA256SUMS $(sha256sum "$OUT/results/SHA256SUMS" | cut -c1-64)" > "$OUT/DONE"
say "$(cat "$OUT/DONE")"
