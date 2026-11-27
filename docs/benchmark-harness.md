# Unified published AMR tool benchmark harness

This page specifies Track A (`scripts/bench/`): a revision-14-safe, offline-capable benchmark harness for comparing published antimicrobial resistance (AMR) prediction tools against GRM's learners on identical isolates, splits, and scoring rules.

The harness produces deterministic, zero-ID result summaries (`grm-bench-result-v1`) and separates end-to-end delivered pipeline comparisons from fixed-learner representation comparisons. No external validation cohort may enter tuning or model selection, and no unsupported endpoint, `no_call`, `I` (intermediate), or runtime failure is ever mapped to Susceptible (`S`).

## Two-mode adapter protocol

Every tool under `scripts/bench/adapters/<tool>/` declares its mode and comparison scope in `manifest.json` (`grm-bench-adapter-v1`), its drug name mapping in `drug_map.tsv`, and its pre-committed species-drug support tiers in `species_support.tsv`.

| Tool | Mode | Comparison kind | Bootstrap cluster unit | Emits `I` |
| --- | --- | --- | --- | --- |
| `amrfinderplus` | `pretrained` | `end_to_end_pipeline` | `species` | No |
| `resfinder_pointfinder` | `pretrained` | `end_to_end_pipeline` | `species` | No |
| `card_rgi` | `pretrained` | `end_to_end_pipeline` | `species` | No |
| `staramr` | `pretrained` | `end_to_end_pipeline` | `species` | No |
| `kleborate` | `pretrained` | `end_to_end_pipeline` | `lineage` | No |
| `tbprofiler` | `pretrained` | `end_to_end_pipeline` | `lineage` | No |
| `mykrobe` | `pretrained` | `end_to_end_pipeline` | `lineage` | Yes (`r` mapped to `I`) |
| `kover` | `trainable` | `end_to_end_pipeline` | `species` | No |
| `xgboost` | `trainable` | `fixed_learner_representation` | `species` | No |
| `lightgbm` | `trainable` | `fixed_learner_representation` | `species` | No |
| `elasticnet` | `trainable` | `fixed_learner_representation` | `species` | No |

### Pretrained mode (`scripts/bench/pretrained.py`)

Pretrained callers run once per unique genome FASTA rather than once per genome-drug row. Each tool emits canonical `genome_calls.jsonl` records in local scratch storage (never committed):

```json
{"call":"R","genome_id":"<id>","reading":"who_assoc","species_slug":"mycobacterium_tuberculosis","tool_drug":"rifampicin"}
```

Required keys are `genome_id`, `species_slug`, `reading`, `tool_drug`, and `call` (`"R"`, `"S"`, `"I"`, `"no_call"`, or `"failure"`). Optional keys are `score` (finite float in `[0, 1]`) and `reason` (mandatory non-empty string whenever `call` is `"no_call"`, `"failure"`, or `"I"`).

When a genome-level runtime error occurs, the harness writes explicit `"failure"` records rather than omitting the genome or defaulting to `"S"`. When projecting multi-drug proxy mappings (such as TBProfiler `levofloxacin,moxifloxacin` mapped to `ciprofloxacin`), calls combine deterministically by precedence `R > failure > no_call > I > S`, so a detected resistance determinant stands over incomplete coverage while failures, no-calls, and `I` never collapse to `"S"`.

Genome-level calls project onto three split types:

1. GRM revision 14 drawn splits (`random`, `lineage`) from `grm export` `dataset.json` or revision 14 compact reports.
2. Supplied fold files (`grm-fold-file-v1`, matching `src/learn/foldfile.rs`) for `study-whole`, `time`, and other reviewed partitions, scoring only `test` folds and recording the raw-byte SHA-256 of the fold file.
3. External deduplicated cohorts (`grm-external-cohort-v1`), which must declare `deduplicated: true` and `selection_role: "evaluation_only"`. Any attempt to perform tuning, threshold selection, or tool/model selection on an `evaluation_only` cohort is refused.

### Trainable mode (`scripts/bench/trainable.py`)

Trainable tools run per `(dataset, split, outer_fold)` consuming `grm export` (`dataset.json` + `matrix.bin`) or a `grm-fold-file-v1` partition:

- **Fold parity.** Drawn splits verify exact outer-fold parity against GRM's own reports via `scripts/c7-outside/plans.py`: fold counts by class (`export_shape` vs. `result_shape`), outer-fold membership by sorted genome-ID SHA-256 (`export_members` vs. `report_members`), and seed-42 `input_fingerprint`.
- **Training-only preprocessing and tuning.** Every learned step (candidate support floor, chi-square feature screening, feature standardization, inner cross-validation, and boosting early-stopping holdouts) executes strictly on outer training folds. Outer test labels are withheld from the `FoldFitContext` passed to every fitter.
- **Per-fold accounting.** Each outer fold emits `{fold, n, n_train, n_train_resistant, n_candidates, n_screened, tp, fp, tn, fn, no_call_r, no_call_s, i_r, i_s, failure_r, failure_s}` plus optional `loss_sum` (clipped to `[1e-4, 1 - 1e-4]` as in `src/cutoff.rs`).
- **Local-only out-of-fold vectors.** Optional out-of-fold score vectors for selective-calling interop (`grm-bench-oof-scratch-v1`) may be written only to an explicit scratch directory verified to sit outside `reference/` and outside the result output directory.

### Canonical drug vocabulary and `I` handling

`scripts/bench/vocab.py` parses `GRM_VOCABULARY`, `INHIBITORS`, `OVERRIDES`, and `single_class` directly from `src/drugs.rs`. Every `canonical_drug` in every adapter's `drug_map.tsv` and `species_support.tsv` must satisfy `src/drugs.rs` canonicalization (lowercase principal ingredients sorted alphabetically, beta-lactamase inhibitors sorted alphabetically last, joined by `/`).

Every adapter manifest declares:

- `i_handling`: whether the tool emits intermediate calls and its policy (`"separate_count_never_map_to_s"`). Tools such as Mykrobe that emit minor-allele or intermediate tokens (`"r"`) map them to canonical `"I"`, counted separately in `i_r` and `i_s`, and never mapped to `"S"`.
- `genotype_to_drug_rule`: the deterministic rule translating tool output into `"R"`, `"S"`, `"I"`, `"no_call"`, or `"failure"`.
- `comparison_kind`: `"end_to_end_pipeline"` for delivered FASTA-to-call pipelines (`amrfinderplus`, `resfinder_pointfinder`, `card_rgi`, `staramr`, `kleborate`, `tbprofiler`, `mykrobe`, `kover`) vs. `"fixed_learner_representation"` for external learners fit on GRM's exported feature representation (`xgboost`, `lightgbm`, `elasticnet`).

## Four-tier pre-committed support registry

Before running a tool, every `(species_slug, canonical_drug)` endpoint in the 274 revision-14 GRM datasets is assigned one of four tiers in `scripts/bench/adapters/<tool>/species_support.tsv`:

| Tier | Pooling rule | Denominator rule |
| --- | --- | --- |
| `supported_primary` | Eligible for `summary_primary` when common-support counts reach at least 25 R and 25 S | Counted in both primary and full-cohort denominators |
| `supported_ungraded` | Reported in `summary_ungraded`; never pooled into `summary_primary` | Counted in full-cohort denominator (e.g., TBProfiler cycloserine and para-aminosalicylic acid with zero WHO v2 graded entries) |
| `proxy_class` | Reported in `summary_proxy_class`; never pooled into `summary_primary` | Counted in full-cohort denominator (class-level or related-drug projections) |
| `unsupported` | Excluded from performance medians; never mapped to `S` | Counted in `full_denominator_coverage` with an explicit `unsupported_reason` |

Any `(species_slug, canonical_drug)` pair absent from `species_support.tsv` (such as an unmapped external-cohort endpoint) defaults to `unsupported` (`not_in_precommitted_support_registry`) and is never mapped to `"S"`.

## Time-split provenance rule

Time-split evaluations must carry a validated `time_provenance` block enforced by `protocol.validate_time_provenance`:

- **`historical_prospective`** is accepted only when `fold_cutoff_date`, `software_release_date`, `training_labels_max_date`, and (for pretrained tools) `database_release_date` and any `weights_release_date` are all explicitly provided as `YYYY-MM-DD` and dated on or before `fold_cutoff_date`.
- **`retrospective_modern_db`** is mandatory whenever a pretrained tool with a modern or post-cutoff database or software release is run on a BV-BRC insertion-date time split (`is_bvbrc_insertion_time_split: true`). Attempting to stamp such a run as `historical_prospective` raises `BenchError`.
- **`retrospective_same_db`** and **`external_retrospective`** label non-prospective splits within BV-BRC or across external cohorts.

## Unified scoring core (`grm-bench-result-v1`)

`scripts/bench/score.py` aggregates projected calls or trainable fold counts into a deterministic `grm-bench-result-v1` JSON document using only Python's standard library and NumPy:

- **Pairwise common-support and full-denominator coverage.** For each comparator learner, tool and learner are restricted to the identical set of isolates where both emit a binary `"R"` or `"S"` call, requiring `called_r >= 25` and `called_s >= 25` for eligibility (`eligible_common_support`). Paired discordance (`tool_only_right`, `learner_only_right`, `both_right`, `both_wrong`, `calls_agree`) and exact two-sided McNemar `mcnemar_p` are reported alongside `full_denominator_coverage` across all four support tiers.
- **Binomial confidence intervals.** Every endpoint rate (`sensitivity`, `specificity`, `vme`, `me`) carries both exact Clopper-Pearson 95% CIs (regularized incomplete beta continued fraction with bisection) and Wilson score 95% CIs.
- **No-call and failure sensitivity bounds.** Every endpoint reports worst-case and best-case bounds for sensitivity, specificity, VME, ME, and balanced accuracy under `no_call + failure` (`bounds_no_call_failure`) and under all uncalled categories (`bounds_all_uncalled`).
- **Clustered bootstraps (4,000 draws, deterministic SHA-256 seed).**
  - Multi-species tools (`bootstrap_cluster_unit: "species"`) resample sorted species with replacement (4,000 draws, 95% interval at the 101st and 3,900th sorted draws) and set `coarse_cluster_warning: true` when fewer than 6 species enter the set.
  - Single-species and species-complex tools (`tbprofiler`, `mykrobe`, `kleborate`; `bootstrap_cluster_unit: "lineage"`) join lineage clusters across endpoints by union-find ordered by smallest genome ID and draw 4,000 lineage-cluster resamples.
- **Zero-ID hygiene and determinism.** `grm-bench-result-v1` contains no genome IDs, no per-genome predictions, no raw features or sequences, no wall clock, and no peak RSS. `check_zero_id_hygiene` enforces a strict top-level key allow-list, a recursive forbidden-key check, accession/sequence regex scans, and `scripts/c7-outside/scan.py`'s exact genome-ID token scan.

## Lambda execution runner (`scripts/bench/run_lambda.py`)

`scripts/bench/run_lambda.py` provides the execution harness for Lambda hosts:

- **Pinned `run-manifest.json`.** Records tool version, database pin, deterministic SHA-256 of the local database directory/file, SHA-256 of `manifest.json`, `drug_map.tsv`, and `species_support.tsv`, Python/NumPy versions, optional environment lock SHA-256, `offline` flag, `split_scheme`, `selection_role`, and validated `time_provenance`.
- **Offline mode (`--offline`).** Strips proxy environment variables, exports offline environment flags (`GRM_BENCH_OFFLINE=1`, `HF_HUB_OFFLINE=1`, `CONDA_OFFLINE=1`, `PIP_NO_INDEX=1`, `NO_NETWORK=1`), requires `--db-path` for pretrained adapters, and verifies `--expected-db-sha256` before starting work.
- **Atomic checkpointing and timing separation.** Completed work units write atomic `checkpoints/<unit_id>.DONE` sentinels keyed to the `run-manifest.json` SHA-256; failures write `checkpoints/<unit_id>.FAILED` and run-level `FAILED` with known genome IDs scrubbed. Resuming an existing output directory skips completed units and refuses any configuration or pin change. Wall seconds, CPU seconds, and peak RSS are written exclusively to `logs/timing.tsv`. The run-level `DONE` file (listing `SHA256SUMS` of `run-manifest.json` and every result file) is written only after every result JSON passes `check_zero_id_hygiene`.

## Verification

Run the synthetic unit suite in Python 3.12 with NumPy 1.26.4 (no external ML or bioinformatics libraries required):

```sh
for seed in 0 1 4242; do
  env -i PATH=/usr/bin:/bin HOME="$HOME" PYTHONHASHSEED="$seed" \
    python3 -s -P -B -W error -m unittest discover -s tests/bench -t . -v
done
python3 -I -B -W error -m unittest discover -s tests/bench -t . -v
```

Preflight a Lambda manifest in offline dry-run mode:

```sh
python3 -I scripts/bench/run_lambda.py \
  --tool tbprofiler \
  --out /tmp/bench-preflight \
  --split-scheme lineage \
  --time-provenance retrospective_same_db \
  --offline \
  --db-path /path/to/local/tbdb \
  --dry-run
```
