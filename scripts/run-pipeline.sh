#!/usr/bin/env bash
# Run from the checkout. Use nohup/setsid outside this script to detach it.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."

export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-64}"
mkdir -p data/logs
exec 9>data/logs/pipeline.lock
flock -n 9 || {
    echo 'Another pipeline owns data/logs/pipeline.lock.' >&2
    exit 1
}

printf '%s\n' "$$" >data/logs/pipeline.pid
rm -f -- data/logs/pipeline.exit
finish() {
    local result=$?
    printf '%s\n' "$result" >data/logs/pipeline.exit
    rm -f -- data/logs/pipeline.pid
    printf '%s pipeline finished with exit code %s\n' "$(date -u +%FT%TZ)" "$result"
}
trap finish EXIT
ulimit -n 8192 || {
    echo 'GRM needs an open-file limit of at least 8192 for the 2000-genome merge.' >&2
    exit 1
}
printf '%s pipeline started; threads=%s; open-files=%s\n' \
    "$(date -u +%FT%TZ)" "$RAYON_NUM_THREADS" "$(ulimit -n)"
nice -n 10 ./target/release/grm pipeline --workers 12 --drop-ksets "$@"
