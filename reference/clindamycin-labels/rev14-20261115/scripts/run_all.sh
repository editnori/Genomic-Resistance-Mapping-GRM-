#!/bin/bash
# Lane clinda: every fit of step two, in order (PLAN.md section 5; addendum A2 to A4).
# Per refit species: grm count, structure, matrix in root-with; arm "with"; arm "without". Then the placebo arms,
# in descending order of D(d,T) (only campylobacter_coli has placebos under the plan's rule), none launched after
# 21:30Z on 15 Nov. Writes the flag /home/user/scratch/clinda/flags/run_all.done at its end.
set -uo pipefail
S=/home/user/scratch/clinda
OUT=/home/user/grm-ops/program-20261114/clinda
A=$OUT/scripts/arms.sh
mkdir -p "$S/reports" "$S/flags"
touch "$S/flags/run_all.start"

threads() {  # addendum A2: 5 when the machine is busy with both other lanes, else 8
  local l; l=$(cut -d' ' -f1 /proc/loadavg)
  if awk -v l="$l" 'BEGIN{exit !(l>6)}'; then echo 5; else echo 8; fi
}
step() { echo "STEP $* $(date -u +%FT%TZ) load $(cut -d' ' -f1-3 /proc/loadavg)"; }
fail() { echo "FAILED $* $(date -u +%FT%TZ)"; echo "rc=1 $*" > "$S/flags/run_all.done"; exit 1; }

[ -d "$S/root-with" ] || { step make-root-with; bash "$A" make-root-with || fail make-root-with; }
for sp in campylobacter_coli staphylococcus_aureus streptococcus_pneumoniae; do
  d=${sp}__clindamycin
  t=$(threads); step stages "$sp" threads "$t"; df -h /home/user/scratch | tail -1
  bash "$A" stages "$S/root-with" "$sp" "$t" || fail stages "$sp"
  bash "$A" guard "$sp" > "$S/reports/$sp.guard-before"
  t=$(threads); step learn-with "$d" threads "$t"
  bash "$A" learn "$S/root-with" "$sp" "$t" || fail learn-with "$d"
  cp "$S/root-with/species/$sp/results/$d.json" "$S/reports/$d.with.json" || fail copy-with "$d"
  step make-without "$d"
  bash "$A" make-root-without "without-$d" "$sp" "$OUT/results/removal/$d.ids" || fail make-without "$d"
  t=$(threads); step learn-without "$d" threads "$t"
  bash "$A" learn "$S/roots/without-$d" "$sp" "$t" || fail learn-without "$d"
  cp "$S/roots/without-$d/species/$sp/results/$d.json" "$S/reports/$d.without.json" || fail copy-without "$d"
  bash "$A" guard "$sp" > "$S/reports/$sp.guard-after-main"
done

for d in campylobacter_coli__clindamycin; do
  sp=${d%__clindamycin}
  for s in $(seq 1 10); do
    if [ "$(date -u +%Y%m%d%H%M)" -ge 202611152130 ]; then step "placebo launch limit reached before seed $s"; break; fi
    step make-placebo "$d" "$s"
    bash "$A" make-root-without "placebo-$s-$d" "$sp" "$OUT/results/placebo/$d/placebo-$s.ids" || fail make-placebo "$d" "$s"
    t=$(threads); step learn-placebo "$d" "$s" threads "$t"
    bash "$A" learn "$S/roots/placebo-$s-$d" "$sp" "$t" || fail learn-placebo "$d" "$s"
    cp "$S/roots/placebo-$s-$d/species/$sp/results/$d.json" "$S/reports/$d.placebo-$s.json" || fail copy-placebo "$d" "$s"
  done
  bash "$A" guard "$sp" > "$S/reports/$sp.guard-after-placebo"
done

step "public root written since start (must be empty):"
find /home/user/grm-h1raw -newer "$S/flags/run_all.start" | head -20
step done
echo "rc=0" > "$S/flags/run_all.done"
