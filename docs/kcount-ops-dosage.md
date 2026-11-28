# Streaming set operations, bounded matrix export, and core-panel dosage (`grm kcount`)

This page documents streaming k-mer set algebra (`grm kcount op`), bounded matrix export (`grm kcount matrix-export`), training-fold core-panel construction (`grm kcount build-core-panel`), and modality-verified core-panel dosage extraction (`grm kcount dosage`), implemented in `src/kcount/ops.rs` and `src/kcount/dosage.rs`. Existing `grm count`, `grm matrix`, `.kset` v1 (`GRMKSET1`), and Revision 14 learners are unchanged.

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
  - `--layout sparse` (default): `feature_idx`, `sample_idx`, `value` for valid present non-zero entries (`1` for `--mode pa`, exact `u32` count for `--mode counts`, `count / d_core` to 6 decimal places for `--mode dosage`). In `--mode dosage --layout sparse`, samples with `insufficient_core_coverage` are recorded only in `missing_mask.tsv` (`missing_core_coverage = 1`) and `manifest.json`; no per-feature `NA` rows are emitted into `matrix.tsv`.
  - `--layout dense`: `feature_idx` followed by one column per sample (`0`/`1` for `--mode pa`, exact `u32` count for `--mode counts`, `count / d_core` to 6 decimal places for passing samples in `--mode dosage`, and explicit `NA` in every row of a column whose sample has `insufficient_core_coverage`).
- `missing_mask.tsv`: `sample_idx`, `sample`, `missing_core_coverage`, `dosage_status`, `d_core`, `core_breadth`.
- `manifest.json`: schema version, `k`, `mode`, `layout`, sample count, feature count, missing-sample count, per-sample mask entries, and SHA-256 digests of `features.tsv`, `matrix.tsv`, and `missing_mask.tsv`.

### Filtering and anonymized naming

- **Support and prevalence bounds:** `--min-support`, `--max-support`, `--min-prevalence` (`[0.0, 1.0]`), and `--max-prevalence` (`[0.0, 1.0]`) filter emitted rows by `support` and `support / n_samples`.
- **Default column labels:** Unless `--sample-names` is explicitly provided (with one unique non-empty name per input), columns are labeled `"s000001"`, `"s000002"`, ..., never copying input filenames, stems, or embedded table header names into `matrix.tsv`, `missing_mask.tsv`, or `manifest.json`.

## 3. Core-panel construction and modality-verified dosage (`src/kcount/dosage.rs`)

### What the code verifies vs what the caller must supply

`grm kcount dosage` and `grm kcount matrix-export --mode dosage` verify:

1. **Core-panel schema and canonical digest:** `CorePanel` (`schema_version = 1`, `k` in `5..=64`, non-empty unique canonical k-mers) requires a non-empty `panel_sha256` matching `SHA-256("GRM-CORE-PANEL-V1\n" || "schema_version=...\n" || "k=...\n" || "input_type=...\n" || sorted_canonical_kmers)`. Missing, empty, or mismatched digests are refused unconditionally.
2. **Stored table modality and integrity:** `.ktab` SHA-256 integrity is verified and the table header's stored `InputModality` (`reads` or `assembly`) must match `panel.input_type`. Tables with `InputModality::Unknown` or cross-modality headers are refused even if `--input-type` on the CLI agrees with the panel. When `--input-type` is supplied, both the CLI flag and the stored `.ktab` header modality must match `panel.input_type`.
3. **Optional fold identifier:** When `--expected-fold-id <ID>` is passed, `panel.fold_id` must equal `<ID>`.
4. **Exact `k` matching:** Sample table `k`, core panel `k`, and optional target panel `k` must all agree.
5. **Core-coverage thresholds and explicit missingness:** `d_core` is the exact median count across core k-mers with count `>= min_count` (average of the two middle values for even length), and `core_breadth` is `core_detected / core_panel_size`. When `d_core < min_core_depth` (default `10.0`) or `core_breadth < min_core_breadth` (default `0.80`), `dosage_status` is `"insufficient_core_coverage"` and every normalized dosage field (`normalized_dosage`, `detected_normalized_dosage`, and per-k-mer `dosage`) serializes as JSON `null` (`None`).
6. **Optional QC vector (`DosageQcVector`):** Every QC field (`total_reads`, `total_bases`, `valid_bases`, `masked_bases`, `gc_fraction`, `usable_read_len_median`, `usable_read_len_mean`, `usable_read_len_p10`, `usable_read_len_p90`, `n50`, `singleton_fraction`) is `Option<...>`. When `--qc-json` is omitted, every QC field serializes as JSON `null` (not `0` or `0.0`). When `--qc-json` is provided, either a `CountSummary`-shaped object (`{"qc": {...}}`) or a standalone QC object is parsed with `Some(...)` for present values and `null` for absent fields.
7. **Default naming:** `grm kcount dosage` defaults `--name` to `"sample"` and never copies input filenames, paths, or embedded table names into `dosage.json`.

**Caller responsibility:** When constructing a core panel for cross-validation with `grm kcount build-core-panel` (`build_training_fold_core_panel(tables: &[PathBuf], cfg: &CorePanelBuildConfig)`), the caller must pass only training-fold tables and exclude held-out test genomes; the binary cannot determine split membership from table files alone.

### Building a training-fold core panel (`grm kcount build-core-panel`)

`build_training_fold_core_panel(tables: &[PathBuf], cfg: &CorePanelBuildConfig)` enforces `REVIEWED_COLUMNS` (`4,096`) before opening files, verifies table SHA-256 digests, requires every table to match `cfg.k` and `cfg.input_type`, computes each training table's exact median k-mer count from a sorted count histogram, and retains up to `cfg.max_kmers` k-mers that satisfy:

- prevalence across training tables `>= cfg.min_prevalence` (`(0.0, 1.0]`, default `0.95`), and
- count `<= sample_median * cfg.max_copy_factor` (`>= 1.0`, default `1.5`) in every training table carrying that k-mer (excluding high-copy plasmid or mobile-element k-mers).

The resulting `CorePanel` JSON records its canonical `panel_sha256` and contains no input file paths or stems.

## 4. CLI usage

```text
grm kcount op union     --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only] [--allow-saturation]
grm kcount op intersect --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only]
grm kcount op diff      --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only]
grm kcount op symdiff   --input <T1> --input <T2> [...] --out <OUT.ktab|OUT.grk> [--name <NAME>] [--presence-only] [--allow-saturation]

grm kcount matrix-export --input <T1> [--input <T2> ...] --out-dir <DIR> \
  [--mode pa|counts|dosage] [--layout sparse|dense] \
  [--min-support <N>] [--max-support <N>] \
  [--min-prevalence <F>] [--max-prevalence <F>] \
  [--sample-names <S1,S2,...>] \
  [--core-panel <PANEL.json>] [--input-type reads|assembly] [--expected-fold-id <ID>] \
  [--min-count <N>] [--min-core-depth <F>] [--min-core-breadth <F>]

grm kcount build-core-panel --input <T1> [--input <T2> ...] --k <5..=64> \
  --input-type reads|assembly [--min-prevalence <F>] [--max-copy-factor <F>] \
  [--max-kmers <N>] [--fold-id <ID>] [--out <PANEL.json>]

grm kcount dosage --table <SAMPLE.ktab> --core-panel <PANEL.json> \
  [--name <NAME>] [--input-type reads|assembly] [--expected-fold-id <ID>] \
  [--targets <TARGETS.json>] [--qc-json <QC.json>] \
  [--min-count <N>] [--min-core-depth <F>] [--min-core-breadth <F>] \
  [--out <DOSAGE.json>]
```
