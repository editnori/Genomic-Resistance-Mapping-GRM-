#!/usr/bin/env bash
# check_run_species.sh KIT RELEASE WORK
#
# Checks run_species.sh's refusals and failure paths on the box, before a
# run: KIT is the run's kit, RELEASE a folder holding the release extracted
# (rev14-274-20261112/), WORK an empty folder for scratch clones and OUTs.
# Each case runs this checkout's run_species.sh from a scratch clone on a
# fresh OUT, with the steps that build the data root flagged done and a stub
# in place of grm. The stub answers `learn` only and reads no genome. It
# writes SYNTHETIC reports ({"SYNTHETIC": true} with an empty fold_models),
# so the check refuses them and a run with no fault ends FAILED at check
# with every report written; in mode copy it writes the release's report
# with a fold_models of null models, which the check takes, so the run goes
# on through C8 to DONE. Two cases make it print the first genome id of the
# species' list on stderr, as grm does in a download retry. Restart cases
# commit a SYNTHETIC file in the scratch clone after the build: under
# scripts/refit/ the restart goes on, under src/ or scripts/assumption-d0-c8/
# it is refused. Two cases hide numpy with a python3 that starts without its
# site packages (the box needs python3-numpy 1.26.4 for the others). Prints
# one line per case and the cases as wanted of cases, and exits 1 when any
# case is not as wanted. Takes about four minutes.
set -uo pipefail
KIT=$(cd "${1:?kit}" && pwd); REL=$(cd "${2:?release folder}" && pwd); H=${3:?work folder}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SP=klebsiella_aerogenes
LIST=klebsiella_aerogenes__ertapenem,klebsiella_aerogenes__imipenem,klebsiella_aerogenes__meropenem
# No dataset of this species is in the plan's MISS set, so D0 is skipped.
SP0=klebsiella_quasipneumoniae
LIST0=klebsiella_quasipneumoniae__cefotaxime,klebsiella_quasipneumoniae__ceftazidime,klebsiella_quasipneumoniae__sulfamethoxazole_trimethoprim
# Acinetobacter baumannii with amikacin is marked M; with cefepime it is not.
SPM=acinetobacter_baumannii
LISTM=acinetobacter_baumannii__cefepime,acinetobacter_baumannii__amikacin
G=(git -c user.name=refit-check -c user.email=refit-check@localhost)
PY3=$(command -v python3)
for s in "$SP" "$SP0"; do
  test -f "$REL/rev14-274-20261112/species/$s.json" || { echo "no release under $REL"; exit 2; }
done
[ -z "$(git -C "$REPO" status --porcelain)" ] || { echo "the checkout has changes; commit them first"; exit 2; }
mkdir -p "$H"; H=$(cd "$H" && pwd)

cat > "$H/stub-grm" <<'STUB'
#!/usr/bin/env bash
# SYNTHETIC stand-in for grm learn. Reads no genome.
data=; sp=; only=
while [ $# -gt 0 ]; do
  case $1 in
    --data) data=$2; shift ;; --species) sp=$2; shift ;; --only) only=$2; shift ;;
  esac
  shift
done
first_id() { python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["genomes"][0]["id"])' "$R/species/$sp.json"; }
case "${STUB_MODE:-}" in
  fail) echo "SYNTHETIC failure in learn" >&2; exit 1 ;;
  skip) only=${only#*,} ;;
  id) i=$(first_id); echo "  $i: attempt 1 failed: SYNTHETIC" >&2; exit 1 ;;
  idok) i=$(first_id); echo "  $i: attempt 1 failed: SYNTHETIC, then fetched" >&2 ;;
  slow) sleep 20 ;;
esac
out=$data/species/$sp/results; mkdir -p "$out"
for d in ${only//,/ }; do
  python3 -I -c '
import json, sys
from collections import OrderedDict
d, out, rel, mode = sys.argv[1:]
if mode == "copy":
    with open(f"{rel}/datasets/{d}.json") as f:
        doc = json.load(f, object_pairs_hook=OrderedDict)
    splits = [{"scheme": s["scheme"], **{l: [None] * s["folds"] for l in ("scm", "cart", "rf", "gbt", "lasso")}}
              for s in doc["splits"] if s.get("learners")]
    doc["fold_models"] = {"features": [], "splits": splits}
elif mode == "nofold":
    doc = {"SYNTHETIC": True, "id": d}
else:
    doc = {"SYNTHETIC": True, "id": d, "fold_models": {"features": [], "splits": []}}
with open(f"{out}/{d}.json", "w") as f:
    json.dump(doc, f)
' "$d" "$out" "$R" "${STUB_MODE:-}"
done
STUB
chmod +x "$H/stub-grm"
# A python3 that does not see the box's numpy: no site packages.
mkdir -p "$H/no-numpy"; printf '#!/bin/sh\nexec %s -S "$@"\n' "$PY3" > "$H/no-numpy/python3"; chmod +x "$H/no-numpy/python3"

# case_ <name>: a scratch clone of this checkout and a fresh OUT with the
# release linked, the stub as the build and the data root's steps flagged done.
case_() {
  local c=$1
  rm -rf "${H:?}/repo-$c" "${H:?}/out-$c"
  git clone -q "$REPO" "$H/repo-$c" && git -C "$H/repo-$c" checkout -q "$(git -C "$REPO" rev-parse HEAD)"
  local out=$H/out-$c
  mkdir -p "$out/flags" "$out/work/target/release" "$out/work/root"
  ln -s "$REL" "$out/work/release"
  cp "$H/stub-grm" "$out/work/target/release/grm"
  echo "commit $(git -C "$H/repo-$c" rev-parse HEAD)" > "$out/binary.txt"
  echo "grm $(sha256sum "$out/work/target/release/grm" | cut -c1-64)" >> "$out/binary.txt"
  echo "code SYNTHETIC" >> "$out/binary.txt"
  for s in release build curated catalog genomes count structure matrix; do touch "$out/flags/$s.done"; done
}
go() {  # go <name> [kit] [species] [list]
  timeout 600 bash "$H/repo-$1/scripts/refit/run_species.sh" "${2:-$KIT}" "${3:-$SP}" "${4:-$LIST}" "$H/out-$1" 2 \
    > /dev/null 2>&1 < /dev/null
}
failed_at() {
  local f=$H/out-$1/FAILED
  if [ -f "$f" ]; then sed -E 's/^FAILED [a-z_]+ (stopped at|at) ([^ ]+) .*/\1 \2/' "$f"; else echo "no FAILED file"; fi
}
logged() { grep -c -- "$2" "$H/out-$1/logs/run.log" 2> /dev/null || true; }
ids() {  # ids <case> copy|raw: genome ids of the release's lists in the copy list, or in the raw logs
  python3 -I -c '
import glob, json, os, re, sys
out, which = sys.argv[1], sys.argv[2]
ids = set()
for p in glob.glob(os.path.join(out, "work/release/rev14-274-20261112/species/*.json")):
    with open(p) as fh:
        ids.update(str(g["id"]) for g in json.load(fh)["genomes"])
tok = re.compile(r"(?<![\d.])\d+\.\d+(?!\d)")
if which == "copy":
    files = ["binary.txt", "expect.txt", "d0.txt", "SHA256SUMS", "DONE", "FAILED"]
    dirs = ["logs", "check", "results"]
else:
    files, dirs = [], ["work/raw-logs"]
for d in dirs:
    for root, _, names in os.walk(os.path.join(out, d)):
        files += [os.path.relpath(os.path.join(root, n), out) for n in names]
n = 0
for f in files:
    p = os.path.join(out, f)
    if os.path.isfile(p):
        with open(p, errors="replace") as fh:
            n += sum(m.group(0) in ids for m in tok.finditer(fh.read()))
print(n)
' "$H/out-$1" "$2"
}
bad=0; cases=0
report() {  # report <case> <wanted> <seen> <ok 0/1>
  printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$([ "$4" = 1 ] && echo ok || echo NOT AS WANTED)"
  cases=$((cases + 1)); [ "$4" = 1 ] || bad=$((bad + 1))
}
later() {  # later <case> <path>: a SYNTHETIC commit in the case's clone that adds <path>
  mkdir -p "$(dirname "$H/repo-$1/$2")"; echo SYNTHETIC > "$H/repo-$1/$2"
  "${G[@]}" -C "$H/repo-$1" add "$2" && "${G[@]}" -C "$H/repo-$1" commit -q -m "SYNTHETIC later commit"
}
restarts() { grep -c '^restart ' "$H/out-$1/binary.txt" 2> /dev/null || true; }
has() { [[ $1 == *"$2"* ]] && echo 1 || echo 0; }

case_ none; STUB_MODE= go none
seen="$(failed_at none), learn ends $(logged none 'learn: end'), reports missing $(logged none 'no report with fold_models')"
report none "at check, learn ends 1, reports missing 0" "$seen" "$(has "$seen" "at check, learn ends 1, reports missing 0")"

case_ done-no-miss; STUB_MODE=copy go done-no-miss "$KIT" "$SP0" "$LIST0"
dl=$(cat "$H/out-done-no-miss/DONE" 2> /dev/null || echo "no DONE file"); n=$(ids done-no-miss copy)
seen="$(failed_at done-no-miss); $(echo "$dl" | grep -o 'd0 skipped: no MISS dataset among the listed' || echo "$dl" | cut -c1-80); c8 rows $(($(wc -l < "$H/out-done-no-miss/results/c8/c8-datasets.tsv" 2> /dev/null || echo 1) - 1)); ids $n"
report done-no-miss "no FAILED file; d0 skipped: no MISS dataset among the listed; c8 rows 15; ids 0" "$seen" \
  "$(has "$seen" "no FAILED file; d0 skipped: no MISS dataset among the listed; c8 rows 15; ids 0")"

case_ fit-fails; STUB_MODE=fail go fit-fails; seen=$(failed_at fit-fails)
report fit-fails "a FAILED file naming learn" "$seen" "$(has "$seen" "at learn")"

case_ report-missing; STUB_MODE=skip go report-missing
seen="$(failed_at report-missing), $(grep -o 'no report with fold_models for [a-z_]*' "$H/out-report-missing/logs/run.log" | head -1)"
report report-missing "at learn, no report with fold_models for klebsiella_aerogenes__ertapenem" "$seen" \
  "$(has "$seen" "at learn, no report with fold_models for klebsiella_aerogenes__ertapenem")"

case_ no-fold-models; STUB_MODE=nofold go no-fold-models
seen="$(failed_at no-fold-models), reports without fold_models $(grep -o 'no report with fold_models for [a-z_,]*' "$H/out-no-fold-models/logs/run.log" | head -1 | tr ',' '\n' | grep -c __)"
report no-fold-models "at learn, reports without fold_models 3" "$seen" "$(has "$seen" "at learn, reports without fold_models 3")"

# A restart on a later checkout that changed scripts/refit/ only: the fits
# are kept (the stub fails if asked to fit), the check runs again, and the
# restart is written in binary.txt.
case_ restart-refit; STUB_MODE=copy go restart-refit "$KIT" "$SP0" "$LIST0"; rm -f "$H/out-restart-refit/DONE"
later restart-refit scripts/refit/SYNTHETIC-later.txt
STUB_MODE=fail go restart-refit "$KIT" "$SP0" "$LIST0"
seen="$(failed_at restart-refit); DONE $(ls "$H/out-restart-refit/DONE" 2> /dev/null | wc -l); learn ends $(logged restart-refit 'learn: end'), learn skipped $(logged restart-refit 'learn: every report present'); check starts $(logged restart-refit 'check: start'); restart lines $(restarts restart-refit); ids $(ids restart-refit copy)"
report restart-refit "no FAILED file; DONE 1; learn ends 1, learn skipped 1; check starts 2; restart lines 1; ids 0" "$seen" \
  "$(has "$seen" "no FAILED file; DONE 1; learn ends 1, learn skipped 1; check starts 2; restart lines 1; ids 0")"

for c in restart-src:src/SYNTHETIC-later.txt restart-d0c8:scripts/assumption-d0-c8/SYNTHETIC-later.txt; do
  n=${c%%:*}; case_ "$n"; later "$n" "${c#*:}"; go "$n"
  seen="$(failed_at "$n"), refused $(logged "$n" 'changed paths outside scripts/refit/'), restart lines $(restarts "$n")"
  report "$n" "at build-check, refused 1, restart lines 0" "$seen" "$(has "$seen" "at build-check, refused 1, restart lines 0")"
done

case_ binary-changed; later binary-changed scripts/refit/SYNTHETIC-later.txt; echo >> "$H/out-binary-changed/work/target/release/grm"
go binary-changed; seen="$(failed_at binary-changed), refused $(logged binary-changed 'the binary changed since the build')"
report binary-changed "at build-check, refused 1" "$seen" "$(has "$seen" "at build-check, refused 1")"

case_ no-numpy-d0; PATH=$H/no-numpy:$PATH go no-numpy-d0
seen="$(failed_at no-numpy-d0), refused $(logged no-numpy-d0 'D0 needs numpy 1.26.4'), release step $(ls "$H/out-no-numpy-d0/logs/release.log" 2> /dev/null | wc -l)"
report no-numpy-d0 "at kit, refused 1, release step 0" "$seen" "$(has "$seen" "at kit, refused 1, release step 0")"

case_ no-numpy-no-d0; STUB_MODE=copy PATH=$H/no-numpy:$PATH go no-numpy-no-d0 "$KIT" "$SP0" "$LIST0"
dl=$(cat "$H/out-no-numpy-no-d0/DONE" 2> /dev/null || echo "no DONE file")
seen="$(failed_at no-numpy-no-d0); $(echo "$dl" | grep -o 'd0 skipped: no MISS dataset among the listed' || echo "$dl" | cut -c1-80)"
report no-numpy-no-d0 "no FAILED file; d0 skipped: no MISS dataset among the listed" "$seen" \
  "$(has "$seen" "no FAILED file; d0 skipped: no MISS dataset among the listed")"

case_ checkout-dirty; echo SYNTHETIC > "$H/repo-checkout-dirty/untracked-SYNTHETIC.txt"
go checkout-dirty; seen=$(failed_at checkout-dirty)
report checkout-dirty "at kit" "$seen" "$(has "$seen" "at kit")"

case_ list-changed; go list-changed; go list-changed "$KIT" "$SP" klebsiella_aerogenes__ertapenem
seen="$(failed_at list-changed), refused $(logged list-changed 'the dataset list differs')"
report list-changed "at kit on the restart, refused 1" "$seen" "$(has "$seen" "at kit, refused 1")"

case_ m-listed; go m-listed "$KIT" "$SPM" "$LISTM"
seen="$(failed_at m-listed), $(grep -o '[0-9]* are marked M' "$H/out-m-listed/logs/run.log" | head -1)"
report m-listed "at kit, 1 are marked M" "$seen" "$(has "$seen" "at kit, 1 are marked M")"

rm -rf "${H:?}/kit-changed"; mkdir -p "$H/kit-changed"
for f in "$KIT"/*; do ln -s "$f" "$H/kit-changed/$(basename "$f")"; done
rm "$H/kit-changed/curated.sha256"; { cat "$KIT/curated.sha256"; echo; } > "$H/kit-changed/curated.sha256"
case_ kit-changed; go kit-changed "$H/kit-changed"; seen=$(failed_at kit-changed)
report kit-changed "at kit (curated.sha256 one byte longer)" "$seen" "$(has "$seen" "at kit")"

case_ id-failed; STUB_MODE=id go id-failed; n=$(ids id-failed copy); raw=$(ids id-failed raw)
seen="$(failed_at id-failed), copy list $n, raw logs $raw"
report id-failed "at learn, copy list 0, raw logs 1" "$seen" "$(has "$seen" "at learn, copy list 0, raw logs 1")"

case_ id-retried; STUB_MODE=idok go id-retried; n=$(ids id-retried copy); raw=$(ids id-retried raw)
seen="$(failed_at id-retried), copy list $n, raw logs $raw"
report id-retried "at check, copy list 0, raw logs 1" "$seen" "$(has "$seen" "at check, copy list 0, raw logs 1")"

case_ second-start
( STUB_MODE=slow go second-start ) & first=$!
sleep 8; STUB_MODE=slow go second-start; wait $first
seen="starts $(logged second-start " start $SP,"), refusals $(logged second-start "another start runs on this OUT")"
report second-start "starts 1, refusals 1" "$seen" "$(has "$seen" "starts 1, refusals 1")"

echo "cases as wanted: $((cases - bad)) of $cases"
[ "$bad" = 0 ]
