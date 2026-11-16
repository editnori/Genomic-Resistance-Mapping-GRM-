#!/usr/bin/env bash
# check.sh
#
# Whether the lasso's column sums keep their running total in a register
# however the crate is split (docs/lasso-layout.md). It builds `grm` from
# this working tree three times, with the release profile, with one codegen
# unit and with 16, each into its own target folder, and checks in each build
# that `column_sum` (src/learn/lasso.rs) is a function of its own and that
# none of its instructions moves a floating-point register to or from the
# stack. Needs nm and objdump (binutils) and an x86-64 build. Exit 1 when a
# build fails the check.
#
#   scripts/lasso-layout/check.sh [DATA_ROOT SPECIES DATASET [ROUNDS]]
#
# With DATA_ROOT, SPECIES and DATASET it then also learns DATASET with each
# build in turn, ROUNDS times (default 3), at RAYON_NUM_THREADS (default 5),
# and prints the lasso's seconds of every run (the sum over both splits of
# the report's lasso `seconds`) and each build's median. DATA_ROOT is a data
# root where SPECIES is counted and built (count, structure, matrix). Its
# reports are not touched: each run writes to a new folder that links every
# input of DATA_ROOT. Builds, logs and reports go under GRM_LAYOUT_WORK, or a
# new folder from mktemp, which is kept. Timing is no CI gate: run it on a
# quiet machine, with nothing else running.
set -euo pipefail

case $# in
    0 | 3 | 4) ;;
    *)
        echo "usage: check.sh [DATA_ROOT SPECIES DATASET [ROUNDS]]" >&2
        exit 2
        ;;
esac
root=$(cd "$(dirname "$0")/../.." && pwd)
work=${GRM_LAYOUT_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/grm-lasso-layout.XXXXXX")}
mkdir -p "$work"
symbol=grm::learn::lasso::column_sum

# Codegen units per layout; empty keeps the release profile's.
layouts=(profile one-unit sixteen-units)
declare -A units=([profile]="" [one-unit]=1 [sixteen-units]=16)

failed=0
# One result per layout, printed after the three builds' output.
results=()
for layout in "${layouts[@]}"; do
    if [ -n "${units[$layout]}" ]; then
        (cd "$root" && CARGO_PROFILE_RELEASE_CODEGEN_UNITS=${units[$layout]} \
            CARGO_TARGET_DIR="$work/target-$layout" cargo build --release --locked --bin grm)
    else
        (cd "$root" && env -u CARGO_PROFILE_RELEASE_CODEGEN_UNITS \
            CARGO_TARGET_DIR="$work/target-$layout" cargo build --release --locked --bin grm)
    fi
    bin=$work/target-$layout/release/grm
    addr=
    size=
    read -r addr size < <(nm -C --print-size "$bin" |
        awk -v s="$symbol" '$3 ~ /^[tT]$/ && $4 == s { print $1, $2; exit }') || true
    if [ -z "$addr" ]; then
        results+=("$layout: FAIL, $symbol is not a function of its own (inlined or renamed)")
        failed=1
        continue
    fi
    start=$((16#$addr))
    stop=$((start + 16#$size))
    objdump -d --no-show-raw-insn --start-address="$start" --stop-address="$stop" "$bin" \
        > "$work/column_sum-$layout.txt"
    # Each instruction line of objdump -d starts with its address, a colon and a tab.
    count=$(grep -cE $'^ *[0-9a-f]+:\t' "$work/column_sum-$layout.txt" || true)
    spills=$(grep -E '%xmm[0-9]+.*\(%r[sb]p\)|\(%r[sb]p\).*%xmm[0-9]+' "$work/column_sum-$layout.txt" || true)
    if [ "$count" -eq 0 ]; then
        results+=("$layout: FAIL, objdump read no instruction of $symbol")
        failed=1
    elif [ -n "$spills" ]; then
        results+=("$layout: FAIL, $symbol moves a floating-point register to or from the stack:" "$spills")
        failed=1
    else
        results+=("$layout: ok, $symbol is $((16#$size)) bytes, $count instructions read, none moves a floating-point register to or from the stack")
    fi
done
echo "machine code of each build's $symbol in $work/column_sum-*.txt"
printf '%s\n' "${results[@]}"
if [ "$#" -eq 0 ] || [ "$failed" -ne 0 ]; then
    exit "$failed"
fi

data=$(cd "$1" && pwd)
species=$2
dataset=$3
rounds=${4:-3}
if [ ! -d "$data/species/$species/matrix" ]; then
    echo "check.sh: $data/species/$species/matrix is missing; build it with grm count, structure and matrix" >&2
    exit 1
fi
export RAYON_NUM_THREADS=${RAYON_NUM_THREADS:-5}

# A data root of links: every entry of DATA_ROOT but species/, and SPECIES
# with every entry but its results.
run_root=$work/root
rm -rf "$run_root"
mkdir -p "$run_root/species/$species"
for entry in "$data"/*; do
    [ "$(basename "$entry")" = species ] || ln -s "$entry" "$run_root/"
done
for entry in "$data/species/$species"/*; do
    [ "$(basename "$entry")" = results ] || ln -s "$entry" "$run_root/species/$species/"
done

seconds() {
    python3 -I -B -c '
import json, sys
r = json.load(open(sys.argv[1]))
print("%.2f" % sum(l["seconds"] for s in r["splits"] for l in s["learners"] if l["learner"] == "lasso"))
' "$1"
}

declare -A runs
for round in $(seq 1 "$rounds"); do
    for layout in "${layouts[@]}"; do
        "$work/target-$layout/release/grm" --data "$run_root" learn --species "$species" \
            --only "$dataset" --force > "$work/learn-$layout-$round.log" 2>&1
        s=$(seconds "$run_root/species/$species/results/$dataset.json")
        runs[$layout]="${runs[$layout]:-} $s"
        echo "round $round $layout: lasso $s s"
    done
done
for layout in "${layouts[@]}"; do
    median=$(printf '%s\n' ${runs[$layout]} | sort -g |
        awk '{ v[NR] = $1 } END { print (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2 }')
    echo "$layout: median lasso $median s over $rounds runs:${runs[$layout]}"
done
echo "builds, logs and reports in $work"
