# Streaming set operations and bounded matrix export (`grm kcount`)

This page documents streaming k-mer set algebra (`grm kcount op`) and bounded matrix export (`grm kcount matrix-export`) implemented in `src/kcount/ops.rs`. Existing `grm count`, `grm matrix`, `.kset` v1 (`GRMKSET1`), and Revision 14 learners are unchanged.

## 1. Streaming set operations (`grm kcount op`)

`grm kcount op` runs bucket-by-bucket k-way merge joins across 2 to `4,096` (`REVIEWED_COLUMNS`) input tables in any mix of `.kset` v1 (`GRMKSET1`), `.grk` (`GRMKGR1\0`), and `.ktab` (`GRMKTB1\0`) formats for `k = 5..=64`:

- **`union`:** Keys present in at least one input. Counts are accumulated with checked `u64` arithmetic. Any summed count exceeding `4,294,967,295` (`2^32 - 1`) fails with an error by default; passing `--allow-saturation` clamps overflowing counts to `2^32 - 1` and records both the exact `saturated_counts` count and the saturation flag in the written `.ktab` header and JSON summary.
- **`intersect`:** Keys present in every input, keeping the minimum count across inputs.
- **`diff`:** Keys present in the first input and absent from all subsequent inputs, keeping the first input's count.
- **`symdiff`:** Keys present in an odd number of inputs (exact symmetric difference for 2 inputs), summing counts across carrying inputs with the same checked `u64` / `--allow-saturation` ceiling capped at `2^32 - 1`.

### Guardrails and metadata contracts

- **Column cap (`REVIEWED_COLUMNS = 4,096`):** Enforced before opening any input file so passing more than 4,096 paths fails immediately without touching the filesystem.
- **Input integrity and equal `k`:** Stored SHA-256 header digests on `.grk` and `.ktab` inputs are verified before merging, and inputs with unequal `k` are rejected.
- **Output extension check:** Output paths must end in `.ktab` by default, or `.grk` when `--presence-only` is set; mismatched extensions are refused.
- **Modality propagation:** Output headers record `reads` if every input records `reads`, `assembly` if every input records `assembly`, and `unknown` otherwise (including any `.kset` v1 input or mixed modalities).
- **Default sample naming:** When `--name` is omitted, output headers and JSON summaries record `"sample"`. Input file stems and paths are never copied into outputs.

## 2. Bounded matrix export (`grm kcount matrix-export`)

`grm kcount matrix-export` streams 1 to `4,096` (`REVIEWED_COLUMNS`) verified tables bucket by bucket into four deterministic files under `--out-dir`:

- `features.tsv`: `feature_idx`, `mixed_key_hex`, `kmer_seq`, `support`.
- `matrix.tsv`:
  - `--layout sparse` (default): `feature_idx`, `sample_idx`, `value` for present non-zero entries (`1` for `--mode pa`, exact `u32` count for `--mode counts`).
  - `--layout dense`: `feature_idx` followed by one column per sample (`0`/`1` for `--mode pa`, exact `u32` count for `--mode counts`).
- `missing_mask.tsv`: `sample_idx`, `sample`, `missing_core_coverage`, `dosage_status`, `d_core`, `core_breadth`.
- `manifest.json`: schema version, `k`, `mode`, `layout`, sample count, feature count, missing-sample count, per-sample mask entries, and SHA-256 digests of `features.tsv`, `matrix.tsv`, and `missing_mask.tsv`.

### Filtering and anonymized naming

- **Support and prevalence bounds:** `--min-support`, `--max-support`, `--min-prevalence` (`[0.0, 1.0]`), and `--max-prevalence` (`[0.0, 1.0]`) filter emitted rows by `support` and `support / n_samples`.
- **Default column labels:** Unless `--sample-names` is explicitly provided (with one unique non-empty name per input), columns are labeled `"s000001"`, `"s000002"`, ..., never copying input filenames, stems, or embedded `.kset` header names into `matrix.tsv`, `missing_mask.tsv`, or `manifest.json`.

## 3. CLI usage

```text
grm kcount op union     --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only] [--allow-saturation]
grm kcount op intersect --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only]
grm kcount op diff      --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only]
grm kcount op symdiff   --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only] [--allow-saturation]

grm kcount matrix-export --input <T1> [--input <T2> ...] --out-dir <DIR> \
  [--mode pa|counts] [--layout sparse|dense] \
  [--min-support <N>] [--max-support <N>] \
  [--min-prevalence <F>] [--max-prevalence <F>] \
  [--sample-names <S1,S2,...>]
```
