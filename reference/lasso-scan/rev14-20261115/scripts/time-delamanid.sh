#!/usr/bin/env bash
# Second timed run of mycobacterium_tuberculosis__delamanid alone (item 5), with the unmodified lasso-scan binary.
set -u
SC=/home/user/scratch/lassoscan; O=/home/user/grm-ops/program-20261114/lassoscan
{ date -u +%FT%TZ; nproc; grep -m1 'model name' /proc/cpuinfo; cat /proc/loadavg; } > "$SC/scan-time2.machine"
python3 -I "$O/scripts/scan14.py" --binary "$SC/target/release/examples/lasso-scan" --data "$SC/data" \
  --release /home/user/grm-ops/owner2/work/rev14-full/rev14-274-20261112 --list "$SC/src/scripts/lasso-scan/scan-list.json" \
  --out "$SC/scan-time2" --threads 16 --only mycobacterium_tuberculosis__delamanid
{ date -u +%FT%TZ; cat /proc/loadavg; } >> "$SC/scan-time2.machine"
