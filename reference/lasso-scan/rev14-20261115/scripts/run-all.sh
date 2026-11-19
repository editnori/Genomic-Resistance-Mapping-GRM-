#!/usr/bin/env bash
# Pass 1 (timed, unmodified lasso-scan), the second delamanid timing, then pass 2 (instrumented, untimed).
set -u
SC=/home/user/scratch/lassoscan; O=/home/user/grm-ops/program-20261114/lassoscan
REL=/home/user/grm-ops/owner2/work/rev14-full/rev14-274-20261112; LIST=$SC/src/scripts/lasso-scan/scan-list.json
while [ ! -e "$SC/logs/builds-all.done" ]; do sleep 10; done
{ echo "pass1 start $(date -u +%FT%TZ)"; nproc; grep -m1 'model name' /proc/cpuinfo; cat /proc/loadavg; free -g | sed -n 2p; } > "$SC/run-all.machine"
python3 -I "$O/scripts/scan14.py" --binary "$SC/target/release/examples/lasso-scan" --data "$SC/data" --release "$REL" \
  --list "$LIST" --out "$SC/scan" --threads 16 > "$SC/logs/pass1.log" 2>&1
{ echo "pass1 end $(date -u +%FT%TZ)"; cat /proc/loadavg; } >> "$SC/run-all.machine"
bash "$O/scripts/time-delamanid.sh" > "$SC/logs/time2.log" 2>&1
echo "pass2 start $(date -u +%FT%TZ)" >> "$SC/run-all.machine"
python3 -I "$O/scripts/scan14.py" --binary "$SC/target-dump/release/examples/lasso-scan" --data "$SC/data" --release "$REL" \
  --list "$LIST" --out "$SC/scan-dump" --threads 16 --dump "$SC/dump" \
  --dump-full "mycobacterium_tuberculosis__delamanid|random|0" > "$SC/logs/pass2.log" 2>&1
echo "pass2 end $(date -u +%FT%TZ)" >> "$SC/run-all.machine"
date -u +%FT%TZ > "$SC/run-all.done"
