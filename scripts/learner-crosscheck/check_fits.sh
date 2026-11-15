#!/usr/bin/env bash
# check_fits.sh
#
# Whether the committed result under reference/learner-crosscheck is of the
# code under test. It builds the Rust harness from this working tree, runs it
# on the 16 saved panels with 1 and with 8 threads, and sets the sha256 of
# every file each run writes against reference/learner-crosscheck/grm.sha256,
# the list run.sh wrote with the committed result. It needs Rust and none of
# the reference libraries, and it reads no genome data.
#
# It exits 0 when both runs wrote exactly the listed files, byte for byte.
# Otherwise it names each file that has other bytes, is missing or is not
# listed, says to run run.sh again, and exits 1.
set -euo pipefail

if [ "$#" -ne 0 ]; then
    echo "usage: check_fits.sh" >&2
    exit 2
fi
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
fixtures=$root/tests/data/lasso-stopping-fixtures.txt
list=reference/learner-crosscheck/grm.sha256

# The sha256 of every file under a folder, sorted by name, in the form
# sha256sum -c reads from that folder. run.sh writes grm.sha256 the same way.
listing() {
    (cd "$1" && find . -type f -printf '%P\n' | LC_ALL=C sort | xargs -r -d '\n' sha256sum --)
}

# One line for each file of a run that is not as the list has it.
differences() {
    listing "$1" | LC_ALL=C awk '
        FILENAME == ARGV[1] { want[substr($0, 67)] = substr($0, 1, 64); next }
        { got[substr($0, 67)] = substr($0, 1, 64) }
        END {
            for (name in want) if (!(name in got)) print name ": listed, not written"
            for (name in got)
                if (!(name in want)) print name ": written, not listed"
                else if (got[name] != want[name]) print name ": other bytes than listed, sha256 " got[name]
        }' "$list" - | LC_ALL=C sort
}

cd "$root"
if [ ! -f "$list" ]; then
    echo "check_fits.sh: $list is missing; run.sh writes it (scripts/learner-crosscheck/README.md)" >&2
    exit 1
fi
cargo build --release --locked --example learner-crosscheck
harness=${CARGO_TARGET_DIR:-$root/target}/release/examples/learner-crosscheck

out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT
failed=false
for threads in 1 8; do
    if ! RAYON_NUM_THREADS=$threads "$harness" "$fixtures" "$out/$threads" 2> "$out/$threads.log"; then
        cat "$out/$threads.log" >&2
        echo "check_fits.sh: the harness stopped with RAYON_NUM_THREADS=$threads" >&2
        exit 1
    fi
    found=$(differences "$out/$threads")
    if [ -n "$found" ]; then
        echo "check_fits.sh: with RAYON_NUM_THREADS=$threads, the harness built from this tree did not write the files $list lists:" >&2
        sed 's/^/  /' <<< "$found" >&2
        failed=true
    fi
done
if $failed; then
    echo "check_fits.sh: so the committed result under reference/learner-crosscheck is not of this tree's code." \
        "Run scripts/learner-crosscheck/run.sh ENV OUT again with the reference libraries and commit its files" \
        "(scripts/learner-crosscheck/README.md)." >&2
    exit 1
fi
echo "check_fits.sh: with 1 and with 8 threads the harness built from this tree wrote the $(wc -l < "$list") files of $list, byte for byte"
