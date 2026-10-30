#!/usr/bin/env bash
# Run the full synthetic assurance tests offline and compare every output file
# with the hashes recorded by the 2026-10-27 receipts.
#
# Usage: scripts/synthetic-assurance/run-full.sh [--debug] [OUT_DIR]
#
# OUT_DIR must be new or empty; by default a fresh directory under ${TMPDIR:-/tmp}
# is used. Release writes about 0.95 GB and takes 1-4 minutes after the build.
# --debug runs the debug profile (5-16 minutes) without the k = 13, 14 check,
# matching the recorded debug run. Outputs are kept for inspection or the
# offline Bend workflow in docs/synthetic-assurance.md; delete them afterwards.
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

profile=release
if [[ ${1-} == --debug ]]; then
    profile=debug
    shift
fi
if (($# > 1)); then
    echo 'Usage: run-full.sh [--debug] [OUT_DIR]' >&2
    exit 2
fi
if (($# == 1)); then
    out=$1
    mkdir -p -- "$out"
    if [[ -n $(ls -A -- "$out") ]]; then
        echo "Refusing to write into non-empty $out." >&2
        exit 1
    fi
else
    out=$(mktemp -d "${TMPDIR:-/tmp}/grm-synthetic-assurance.XXXXXX")
fi
out=$(cd -- "$out" && pwd)
manifest="$here/expected-$profile.sha256"
cd -- "$here/../.."

flags=(--locked --offline)
kmer_env=(GRM_KMER_OUT="$out/kmer")
kmer_args=(--nocapture --test-threads=1)
if [[ $profile == release ]]; then
    flags+=(--release)
    kmer_env+=(GRM_KMER_LARGE_OUT="$out/kmer-large.json")
    kmer_args+=(--include-ignored)
fi

printf '%s %s run; output in %s\n' "$(date -u +%FT%TZ)" "$profile" "$out"
env GRM_SYNTH_OUT="$out/foldmask" cargo test "${flags[@]}" --test foldmask_extended \
    -- --ignored --exact extended_foldmask_vectors --nocapture --test-threads=1 \
    2>&1 | tee "$out/foldmask.log"
env "${kmer_env[@]}" cargo test "${flags[@]}" --test kmer_space_extended \
    -- "${kmer_args[@]}" 2>&1 | tee "$out/kmer.log"

# Every recorded file must match, and no unrecorded output may exist.
(cd -- "$out" && sha256sum --check --strict --quiet "$manifest")
diff <(cut -c67- "$manifest") \
    <(cd -- "$out" && find foldmask kmer kmer-large.json -type f 2>/dev/null | LC_ALL=C sort)
printf '%s %s: %s output files match %s\n' \
    "$(date -u +%FT%TZ)" "$profile" "$(wc -l <"$manifest")" "${manifest#"$PWD"/}"
