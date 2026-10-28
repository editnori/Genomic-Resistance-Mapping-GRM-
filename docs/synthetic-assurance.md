# Synthetic assurance tests

Two integration test targets check the fold-mask and k-mer code against specifications computed inside the tests, using seeded synthetic integers and byte strings only. They call the public `grm` functions directly (`grm::learn::{fold_masks, training_folds, validate_folds}`, `grm::bits`, `grm::dna::KmerSpace`) and write every input and output as evidence for an independent oracle and for the offline Bend2 checks below. No genome, phenotype or report data is read.

The full runs were recorded on 2026-10-27 against source commit 9f193bf (archive sha256 `2c32a295…a691a2`) with Rust 1.98.1, `--locked --offline`, in `/root/grm-synthetic-assurance-20261027-cont-yAuKsr` on Vast; its `receipts/VERDICT.md`, `final-results.json` and `output-files.sha256` are the authority for the numbers here. `src/` has not changed since that commit.

## Tests

| Test | Runs by default | Size | Debug | Release |
|---|---|---|---|---|
| `foldmask_extended::extended_foldmask_smoke` | yes | 11,815 case IDs, temporary directory | (10–14 s) | (1 s) |
| `foldmask_extended::extended_foldmask_vectors` | no, `#[ignore]` | 2,669,501 case IDs, about 0.9 GB | 317 s (940 s) | 36 s (112 s) |
| `kmer_space_extended::multi_k_space_vectors` | yes | k = 5..=32, about 30 MB | 4 s (13–17 s) | 12 s with the large test (91 s) |
| `kmer_space_extended::multi_k_mixer_bijection_large` | no, `#[ignore]` | k = 13, 14, 32 MiB bitmap | not run | (included above) |

Times are from the recorded Vast run, with a 12-core WSL workstation in parentheses.

The full fold-mask run assigns 2,462,269 exhaustive case IDs (every assignment over `0..k`, `k` and `usize::MAX` with every eligibility mask, for k = 2, 3 at n = 6, k = 4, 5 at n ≤ 5, k = 6..=8 at n ≤ 4), 100,000 random cases up to n = 4,097 and k = 300 with sentinel fold IDs up to `usize::MAX`, 33,732 ext45 Bend cases (k = 4, 5, n ≤ 4), 1,500 Bend random cases, 2,000 out-of-contract k ∈ {0, 1} cases, 50,000 `validate_folds` cases and 20,000 `training_folds` cases. The test derives the expected total from its counts and asserts both that and the recorded 2,669,501. The smoke test runs every tier and every assertion with the same generators at smaller counts (exhaustive n ≤ 2..4, ext45 n ≤ 2, 1,000 random, 150 Bend random, 200 out-of-contract, 2,000 validate, 1,000 training); it is a regression check, not recorded evidence, and its output is not hashed.

The k-mer vectors run in full by default: 250 sequences, 4,012 values and 1,000 texts per k (7,000 sequences, 576,765 canonical k-mers, 112,336 values, 28,000 texts), plus an exhaustive mix/unmix bijection for every k ≤ 12 (up to 2^24 values) with full mix-table digests for k ≤ 10. The large test extends the bijection to k = 13 (2^26 values) and k = 14 (2^28); it holds a 32 MiB bitmap and mixes and unmixes 3.4 × 10^8 values; it has only been run in release, so run it with `--release`.

## Running

```bash
cargo test --locked --offline                       # routine: smoke and k-mer vectors
scripts/synthetic-assurance/run-full.sh             # full release run, checked against the receipts
scripts/synthetic-assurance/run-full.sh --debug     # full debug run, as recorded (no k = 13, 14); 15+ min
```

`run-full.sh [--debug] [OUT_DIR]` is offline (`--locked --offline`), works from any directory, refuses a non-empty `OUT_DIR`, and by default writes to a new `mktemp -d` directory under `${TMPDIR:-/tmp}`, which it leaves in place for inspection; delete it afterwards. It fails unless every output file matches `scripts/synthetic-assurance/expected-<profile>.sha256` and no unrecorded file exists.

The equivalent manual commands keep the outputs where the variables point:

```bash
GRM_SYNTH_OUT=$OUT/foldmask cargo test --locked --offline --release \
    --test foldmask_extended -- --ignored --exact extended_foldmask_vectors --nocapture
GRM_KMER_OUT=$OUT/kmer GRM_KMER_LARGE_OUT=$OUT/kmer-large.json cargo test --locked --offline \
    --release --test kmer_space_extended -- --include-ignored --nocapture
(cd "$OUT" && sha256sum --check --strict "$OLDPWD/scripts/synthetic-assurance/expected-release.sha256")
```

| Variable | Used by | Meaning |
|---|---|---|
| `GRM_SYNTH_OUT` | `extended_foldmask_vectors` | Empty directory to keep the full fold-mask output; a temporary directory otherwise |
| `GRM_KMER_OUT` | `multi_k_space_vectors` | Empty directory to keep the k-mer output; a temporary directory otherwise |
| `GRM_KMER_LARGE_OUT` | `multi_k_mixer_bijection_large` | New file for the k = 13, 14 result; printed only otherwise |

A named directory is created if missing and must be empty, and `GRM_KMER_LARGE_OUT` must not exist, so nothing is overwritten. The recorded commands set these variables without `--ignored`; the full tests are now ignored by default, so those commands run only the smoke and k-mer vector tests and leave `GRM_SYNTH_OUT` and `GRM_KMER_LARGE_OUT` untouched. The smoke test never writes to `GRM_SYNTH_OUT`.

## Output hashes

`expected-release.sha256` (60 files) and `expected-debug.sha256` (59 files, no k = 13, 14 result) are the receipt's `output-files.sha256` entries for `out/{foldmask,kmer}-{release,debug}-r2` and `out/kmer-large-release-r2.json`, with paths rewritten to `foldmask/`, `kmer/` and `kmer-large.json`. Debug and release outputs are byte-identical except `foldmask/summary.json`, which records `debug_assertions` (`4034ffd7…` release, `240c2de2…` debug). The large outputs are:

| File | Bytes | sha256 |
|---|---|---|
| `foldmask/exhaustive.txt` | 309,079,244 | `23a24376044601ceb6fd8c379efbe65f83b88bbd5edf04c03322839ed890004e` |
| `foldmask/random.jsonl` | 521,791,514 | `0793988c1f671b0671ac236a753c2b610e042ab911bfaec23970c5d32593be3c` |
| `foldmask/bend-cases.jsonl` | 24,256,068 | `4d9774ffe9e87d41f53636e321774058a388eda386ea0b6b10f54d6df0c89d3a` |
| `kmer/sequences.jsonl` | 16,537,966 | `ee6b8a5507556308be14bcc6fed3e714d3993fc9977b6e47a750e7e3f477ef27` |

The hashes depend on the locked `rand`, `rand_chacha` and `serde_json` versions and a 64-bit `usize`. Do not commit generated outputs; keep them outside the checkout or delete them.

## Bend2 and kernel checks (offline, separate)

Bend and Lean are not cargo dependencies and nothing in `cargo test` runs them. The full fold-mask run writes two Bend books: `foldmask/bend/ext45/` (33,732 cases in 33 shards) and `foldmask/bend/rand/` (1,500 Bend random plus 40 out-of-contract cases, 1,540 in 13 shards). Each `cases_NNNN.bend` imports `../source/MODEL.bend` and `./CHECK.bend`, which are not in this repository. They and the toolchain come from the Bend2 fold-mask lane on Vast:

| Input | sha256 |
|---|---|
| `source/MODEL.bend` | `563a6865ad1d3c97630c3605233681c9d92c1b32bd0738f3c4b01b9371393ec9` |
| `source/LAWS.bend` | `9cf4789ea6d8c9353eb2b59f8d21781d37070f5ceb54066f53221ddf8038a83e` |
| `source/PROOF.bend` | `654ea84369230f8f557bcf7fabf7d3571a4b2bf65fa9e5e23d5eb60ae00d2d8d` |
| `CHECK.bend` | `803cd790dcfb757603c3df3766279697313f6ccd9d6c9d3f20ace2e6858ac616` |
| Bend 2.0.39 `bend-2.0.39-linux-x64.tar.gz` | `6cf462fd8f70b094a4950ebb709a9c4ffd8f1d1189d9cedf1fd63cd5d041c7ab` |
| Lean 4.34.0 `lean-4.34.0-linux.zip` | `5f14e0f3762058a529b1598f4355e6787c5f155be1bcce72ae288b9c35e3adb4` |
| `bendtt` kernel | `e943466a914ee7eb0ca53c8d720ff57ed07c55c4026a861362fe17faaba3e05f` |

The model files are in `/root/grm-bend2-foldmask-20261026/{source,vectors}` and the toolchain in `/root/grm-bend2-pilot-20261026`. To check a run, verify those hashes, copy each book's `VECTORS.bend` and `cases_*.bend` into a fresh `bendcheck/BOOK/` beside a copy of `CHECK.bend`, put `MODEL.bend` in `bendcheck/source/`, make the copies read-only, and run each step for `PROOF.bend` and for `bendcheck/{ext45,rand}/VECTORS.bend`:

```bash
bend FILE                                  # check
bend FILE -o FILE.bendtt                   # emit the kernel term
bend FILE --verdict
LEAN_STACK_SIZE_KB=4194304 bendtt FILE.bendtt
```

Each must exit 0 and report `ALL PROOFS CHECK`. On 2026-10-27 all 12 steps passed (vector steps 35–115 s each). The receipt's `bendcheck-staged.sha256` lists the staged files; the Bend book hashes in `expected-*.sha256` must match them.

## Trust boundaries

- The specification the Rust tests check against is written in the same test files; it is not an independent implementation. No independent oracle has yet re-derived the emitted JSONL (random, validate, training, k-mer files).
- `validate_tier` records each `validate_folds` outcome without asserting it. `training_tier` asserts the fold-mask precondition, group integrity and independence from samples outside `train`, and records `validate_folds` results.
- The k-mer test checks mix/unmix inverses, the exhaustive bijection for k ≤ 14, reverse-complement involution and canonical ordering; `canonical()`, `decode()` and the window extraction of `for_each_canonical` are recorded, not checked against an independent extractor.
- Only 35,272 of the 2,669,501 fold-mask cases reach Bend. The universal laws in `PROOF.bend` are about `MODEL.bend`, not about the Rust code; the Rust-to-model abstraction (fold IDs at or above k become `Outside`) and the Bend-to-`bendtt` translation are trusted, not proved.
- The recorded runs used one toolchain, platform and seed set. The smoke test covers the same code paths at small counts and does not reach the full exhaustive domains or the largest random sizes.
