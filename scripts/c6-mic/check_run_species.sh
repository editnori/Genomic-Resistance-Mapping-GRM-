#!/usr/bin/env bash
# check_run_species.sh KIT RELEASE WORK
#
# Checks run_species.sh's refusals on the box, before a run: KIT is the run's
# kit, RELEASE a folder holding the release extracted (rev14-274-20261112/),
# WORK an empty folder for scratch clones and OUTs. Each case runs this
# checkout's run_species.sh from a scratch clone on a fresh OUT, with the
# steps that build the data root flagged done and a stub in place of grm. The
# stub writes SYNTHETIC results ({"SYNTHETIC": true}, no row) and reads no
# data; two cases make it print the first genome id of the species' list on
# stderr, as grm does in a download retry or a fold file check. The scorer
# refuses SYNTHETIC results, so a run with no fault ends FAILED at score with
# no result missing. Prints one line per case and exits 1 when any case is
# not as wanted. Takes about three minutes.
set -uo pipefail
KIT=$(cd "${1:?kit}" && pwd); REL=$(cd "${2:?release folder}" && pwd); H=${3:?work folder}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SP=campylobacter_coli
LIN=campylobacter_coli__azithromycin,campylobacter_coli__ciprofloxacin
G=(git -c user.name=c6-mic-check -c user.email=c6-mic-check@localhost)
test -f "$REL/rev14-274-20261112/species/$SP.json" || { echo "no release under $REL"; exit 2; }
[ -z "$(git -C "$REPO" status --porcelain)" ] || { echo "the checkout has changes; commit them first"; exit 2; }
mkdir -p "$H"; H=$(cd "$H" && pwd)
DSJ=$(tar -tf "$KIT/curated.tar" | grep -m1 'datasets.json$')

cat > "$H/stub-grm" <<'STUB'
#!/usr/bin/env bash
# SYNTHETIC stand-in for grm learn. Reads no data.
seed=42; only=; out=; arms=released,ordered
while [ $# -gt 0 ]; do
  case $1 in
    --only) only=$2; shift ;; --seed) seed=$2; shift ;; --exp-mic-out|--fold-out) out=$2; shift ;;
    --exp-mic) arms=$2; shift ;;
  esac
  shift
done
first_id() { python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["genomes"][0]["id"])' "$R/species/$SP.json"; }
case "${STUB_MODE:-}:$seed" in
  fail:43) echo "SYNTHETIC failure at seed 43" >&2; exit 1 ;;
  skip:43) only=${only#*,} ;;
  id:44) i=$(first_id); echo "  $i: attempt 1 failed: SYNTHETIC" >&2; echo "SYNTHETIC: test fold 0 tests genome $i" >&2; exit 1 ;;
  idok:44) i=$(first_id); echo "  $i: attempt 1 failed: SYNTHETIC, then fetched" >&2 ;;
  slow:42) sleep 20 ;;
esac
mkdir -p "$out"
for d in ${only//,/ }; do
  python3 -I -c '
import json, sys
d, seed, arms, out = sys.argv[1:]
doc = {"SYNTHETIC": True, "format": "grm-mic-plan-result-v1", "id": d, "plan": {"kind": "lineage", "seed": int(seed)},
       "arms": sorted(arms.split(",")), "tested": [], "folds": [], "n_mic": 0}
json.dump(doc, open(f"{out}/{d}.lineage-s{seed}.mic.json", "w"))
' "$d" "$seed" "$arms" "$out"
done
STUB
chmod +x "$H/stub-grm"

# case <name>: a scratch clone of this checkout and a fresh OUT with the
# release linked, the stub as the build and the data root's steps flagged done.
case_() {
  local c=$1
  rm -rf "${H:?}/repo-$c" "${H:?}/out-$c"
  git clone -q "$REPO" "$H/repo-$c" && git -C "$H/repo-$c" checkout -q "$(git -C "$REPO" rev-parse HEAD)"
  local out=$H/out-$c
  mkdir -p "$out/flags" "$out/work/target/release" "$out/work/root/curated"
  ln -s "$REL" "$out/work/release"
  cp "$H/stub-grm" "$out/work/target/release/grm"
  tar -xf "$KIT/curated.tar" -C "$out/work/root/curated" "$DSJ"
  echo "commit $(git -C "$H/repo-$c" rev-parse HEAD)" > "$out/binary.txt"
  echo "grm $(sha256sum "$out/work/target/release/grm" | cut -c1-64)" >> "$out/binary.txt"
  for s in release build curated catalog genomes count structure matrix; do touch "$out/flags/$s.done"; done
}
go() {  # go <name> [kit] [lineage]
  timeout 600 bash "$H/repo-$1/scripts/c6-mic/run_species.sh" "${2:-$KIT}" "$SP" "${3:-$LIN}" - "$H/out-$1" 2 \
    > /dev/null 2>&1 < /dev/null
}
failed_at() {
  local f=$H/out-$1/FAILED
  if [ -f "$f" ]; then sed -E 's/^FAILED [a-z_]+ (stopped at|at) ([^ ]+) .*/\1 \2/' "$f"; else echo "no FAILED file"; fi
}
ids_out() {  # genome ids of the species' list in what RUN.md copies back
  python3 -I -c '
import json, os, re, sys
out, sp = sys.argv[1], sys.argv[2]
ids = {g["id"] for g in json.load(open(os.path.join(out, "work/release/rev14-274-20261112/species", sp + ".json")))["genomes"]}
tok = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")
files = ["binary.txt", "expect.tsv", "SHA256SUMS", "DONE", "FAILED", "results/study-plans.tsv", "results/release-specs.tsv"]
for d in ["score", "logs"]:
    if os.path.isdir(os.path.join(out, d)):
        files += [os.path.join(d, f) for f in os.listdir(os.path.join(out, d))]
n = 0
for f in files:
    p = os.path.join(out, f)
    if os.path.isfile(p):
        with open(p, errors="replace") as fh:
            n += sum(m.group(0) in ids for m in tok.finditer(fh.read()))
print(n)
' "$H/out-$1" "$SP"
}
bad=0
report() {  # report <case> <wanted> <seen> <ok 0/1>
  printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$([ "$4" = 1 ] && echo ok || echo NOT AS WANTED)"
  [ "$4" = 1 ] || bad=$((bad + 1))
}
has() { [[ $1 == *"$2"* ]] && echo 1 || echo 0; }

case_ none; STUB_MODE= go none
seen="$(failed_at none), $(grep -ho '[0-9]* missing' "$H/out-none/logs/run.log" | tail -1)"
report none "at score, 0 missing" "$seen" "$(has "$seen" "at score, 0 missing")"

case_ fit-fails; STUB_MODE=fail go fit-fails; seen=$(failed_at fit-fails)
report fit-fails "a FAILED file naming lineage-s43" "$seen" "$(has "$seen" "at lineage-s43")"

case_ result-missing; STUB_MODE=skip go result-missing; seen=$(failed_at result-missing)
report result-missing "at lineage-s43, no result written" "$seen" \
  "$([[ $seen == *"at lineage-s43"* ]] && grep -q "no result written" "$H/out-result-missing/logs/run.log" && echo 1 || echo 0)"

case_ checkout-moved; "${G[@]}" -C "$H/repo-checkout-moved" commit -q --allow-empty -m "SYNTHETIC later commit"
go checkout-moved; seen=$(failed_at checkout-moved)
report checkout-moved "at build-check" "$seen" "$(has "$seen" "at build-check")"

case_ checkout-dirty; echo SYNTHETIC > "$H/repo-checkout-dirty/untracked-SYNTHETIC.txt"
go checkout-dirty; seen=$(failed_at checkout-dirty)
report checkout-dirty "at kit" "$seen" "$(has "$seen" "at kit")"

case_ list-changed; go list-changed; go list-changed "$KIT" campylobacter_coli__azithromycin; seen=$(failed_at list-changed)
report list-changed "at kit on the restart" "$seen" "$(has "$seen" "at kit")"

rm -rf "${H:?}/kit-changed"; mkdir -p "$H/kit-changed"
for f in "$KIT"/*; do ln -s "$f" "$H/kit-changed/$(basename "$f")"; done
rm "$H/kit-changed/curated.sha256"; { cat "$KIT/curated.sha256"; echo; } > "$H/kit-changed/curated.sha256"
case_ kit-changed; go kit-changed "$H/kit-changed"; seen=$(failed_at kit-changed)
report kit-changed "at kit (curated.sha256 one byte longer)" "$seen" "$(has "$seen" "at kit")"

case_ id-failed; STUB_MODE=id go id-failed; n=$(ids_out id-failed)
report id-failed "at lineage-s44, 0 genome ids in the copy-back" "$(failed_at id-failed), $n" \
  "$([ "$n" = 0 ] && [[ $(failed_at id-failed) == *"at lineage-s44"* ]] && echo 1 || echo 0)"

case_ id-retried; STUB_MODE=idok go id-retried; n=$(ids_out id-retried)
report id-retried "0 genome ids in the copy-back" "$(failed_at id-retried), $n" "$([ "$n" = 0 ] && echo 1 || echo 0)"

case_ second-start
( STUB_MODE=slow go second-start ) & first=$!
sleep 8; STUB_MODE=slow go second-start; wait $first
starts=$(grep -c " start $SP," "$H/out-second-start/logs/run.log")
refused=$(grep -c "another start runs on this OUT" "$H/out-second-start/logs/run.log")
report second-start "one start logged, one refusal" "starts $starts, refusals $refused" \
  "$([ "$starts" = 1 ] && [ "$refused" = 1 ] && echo 1 || echo 0)"

echo "cases not as wanted: $bad"
[ "$bad" = 0 ]
