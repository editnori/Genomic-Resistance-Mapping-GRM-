#!/usr/bin/env bash
# Item 5, run 2: the 30 independent path solves in 15 parallel workers (one
# thread each), then item5.py, which resumes from the saved solutions and does
# the table, the held-out calls and the time-boxed second solver.
set -u
SC=/home/user/scratch/lassoscan; O=/home/user/grm-ops/program-20261114/lassoscan
cd "$SC" || exit 1
D=dump/mycobacterium_tuberculosis__delamanid__random__0.json
mkdir -p item5/ind logs/item5-workers
{ echo "start $(date -u +%FT%TZ)"; cat /proc/loadavg; } > item5/run2.machine
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
pids=()
for w in $(seq 0 14); do
  venv/bin/python -I "$O/scripts/item5-solve.py" "$D" item5 "$w" "$((29 - w))" > "logs/item5-workers/w$w.log" 2>&1 &
  pids+=("$!")
done
wrc=0
for p in "${pids[@]}"; do wait "$p" || wrc=1; done
echo "workers done $(date -u +%FT%TZ) rc=$wrc" >> item5/run2.machine
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
venv/bin/python -I "$O/scripts/item5.py" "$D" item5 --saga-seconds 1800 > logs/item5.log 2>&1
rc=$?
echo "end $(date -u +%FT%TZ)" >> item5/run2.machine
echo "rc=$rc workers_rc=$wrc" > item5/run.rc
