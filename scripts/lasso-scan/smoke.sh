#!/usr/bin/env bash
# End-to-end check of the lasso scan harness on a synthetic data root, with no
# real data: builds the lasso-scan example, writes the 150-genome fixture root
# with an ignored analysis test, and runs test_scan.py, whose RealBinary tests
# drive the real binary under the registered rule.
#
# Usage: scripts/lasso-scan/smoke.sh
#
# Writes only a new directory from mktemp (kept for inspection, about 1 MB)
# and the cargo target directory ($CARGO_TARGET_DIR or target/).
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$here/../.."

work=$(mktemp -d "${TMPDIR:-/tmp}/grm-lasso-scan-smoke.XXXXXX")
target=${CARGO_TARGET_DIR:-$PWD/target}
cargo build --locked --offline --example lasso-scan
GRM_LASSO_SCAN_FIXTURE="$work/root" cargo test --locked --offline --lib -- \
    --ignored --exact analysis::tests::write_lasso_scan_fixture_root
GRM_LASSO_SCAN_BIN="$target/debug/examples/lasso-scan" GRM_LASSO_SCAN_FIXTURE="$work/root" \
    python3 -I -B scripts/lasso-scan/test_scan.py
printf 'fixture root kept in %s\n' "$work"
