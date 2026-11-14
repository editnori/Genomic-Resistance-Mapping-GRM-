#!/usr/bin/env bash
# run.sh ENV OUT
#
# The whole learner cross-check, from this checkout. ENV is a Python
# environment with the pinned reference libraries (README.md says how to make
# it). OUT is a new folder for the results.
#
# It builds the Rust harness, runs it on the 16 saved panels with 1 and with 8
# threads and requires the same bytes from both, sets the fits against the
# references (compare.py), and runs the exploratory counts (exploratory.py).
# It writes OUT/run.json: the commit, the tree of src/, the compiler and the
# times. It uses at most 8 threads and reads no genome data.
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: run.sh ENV OUT" >&2
    exit 2
fi
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
python=$1/bin/python
out=$2
fixtures=$root/tests/data/lasso-stopping-fixtures.txt
jobs=${CROSSCHECK_JOBS:-8}

if [ -e "$out" ]; then
    echo "run.sh: $out exists; give a new folder" >&2
    exit 2
fi
mkdir -p "$out"
started=$(date -u +%Y-%m-%dT%H:%M:%SZ)

"$python" -I "$here/check_env.py" > "$out/environment.tsv"

cd "$root"
cargo build --release --locked --example learner-crosscheck
harness=${CARGO_TARGET_DIR:-$root/target}/release/examples/learner-crosscheck

RAYON_NUM_THREADS=1 "$harness" "$fixtures" "$out/grm"
RAYON_NUM_THREADS=8 "$harness" "$fixtures" "$out/grm-8-threads"
diff -r "$out/grm" "$out/grm-8-threads"
rm -r "$out/grm-8-threads"

"$python" -I "$here/compare.py" "$fixtures" "$out/grm" "$out" --jobs "$jobs"
"$python" -I "$here/exploratory.py" "$fixtures" "$out/grm" "$out" > /dev/null

if git diff --quiet HEAD -- src Cargo.lock; then unchanged=true; else unchanged=false; fi
cat > "$out/run.json" <<EOF
{
 "head": "$(git rev-parse HEAD)",
 "src_tree": "$(git rev-parse HEAD:src)",
 "src_and_lockfile_as_committed": $unchanged,
 "rustc": "$(rustc --version)",
 "harness_identical_at_1_and_8_threads": true,
 "started": "$started",
 "ended": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF
echo "run.sh: results in $out"
