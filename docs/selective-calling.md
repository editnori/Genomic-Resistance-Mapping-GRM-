# Nested selective calling and abstention

An experiment, off by default, for evaluating whether calibrated abstention reduces false-susceptible calls (very major errors, VME) under lineage and study shifts, and whether nearest-training genomic distance adds value beyond classifier score margin alone.

`grm learn --exp-selective <DIR>` turns the experiment on for that run. Without the flag, `grm learn` is unchanged byte for byte, `EVALUATION_REVISION` remains 14, `REPORT_SCHEMA` remains 2, `LearnerResult.conformal` remains empty, and no report or fold-file output changes.

Every split here (`random`, `lineage`, `study-whole`, `time`) is a retrospective split inside one public database. A retrospective split by lineage, study, or date is not external or clinical validation.

## Design principles

1. **Zero outer-test leakage.** Every calibration threshold `(tau_blind, t_lo, t_hi, tau_margin, tau_distance, tau_joint)`, every class-conditional Mondrian conformal quantile `(q_r, q_s)`, and every parameter of the joint margin+distance risk model is fit strictly inside outer training fold `T_f`. Outer test labels `y[E_f]`, outer test features `X[E_f]`, and outer test distances `D[E_f, *]` are never read during calibration.
2. **Fold-safe Mash distances.** Pairwise Mash distances are computed over active cohort genomes from cached MinHash sketches (`species/<slug>/sketches.bin` if present, otherwise read-only cached `species-check/sketches/<genome_id>.json` at `k = 21, s = 1000`) via `cluster::distance_matrix`, never from `structure.json`'s species-wide `nearest` (which leaks across outer test folds):
   - Calibration sample `i in T_f`: `d_inner(i) = min_{j in T_f, inner_fold[j] != inner_fold[i]} D[i, j]`, strictly excluding outer test fold `E_f` and `i`'s own inner fold (and entire lineage group on `lineage` splits).
   - Outer test sample `u in E_f`: `d_outer(u) = min_{j in T_f} D[u, j]`, strictly against outer training references `T_f` only.
   - When sketches are absent, score-only policies (`outcome_blind`, `score_margin`, and Mondrian conformal) still run while distance-dependent policies (`mash_distance`, `margin_plus_distance`) are explicitly marked unavailable with a zero-ID reason.
3. **Nested out-of-inner-fold score collection.** Inside each outer training fold `T_f`, a nested inner-fold collector (`collect_inner_oof` in `src/learn/selective.rs`, reusing `inner_pick::inner_partition` and `inner_pick::inner_seed`) draws `k` group-separated inner folds using `tuning_groups(scheme, groups, duplicates)` (lineages for `lineage`; k-mer-identical copy groups for `random`, `study-whole`, and `time`), fits each of the five release learners (`scm`, `cart`, `rf`, `gbt`, `lasso`) strictly on `T_f \ inner_j`, and returns raw out-of-inner-fold scores `s_inner[l][i]` and null-fit flags for every `i in T_f`.

## Pre-committed constants

Frozen in `src/learn/selective.rs` and recorded in every `grm-selective-v1` file:

| Constant | Value | Role |
|---|---|---|
| `PRIMARY_COVERAGE` | `0.85` | Pre-committed primary target coverage for cross-policy contrasts |
| `COVERAGE_LADDER` | `[1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50]` | Full coverage ladder evaluated per policy |
| `CONFORMAL_ALPHAS` | `[0.01, 0.02, 0.05, 0.10, 0.15, 0.20]` | Class-conditional Mondrian conformal error levels |
| `SPARSE_DENOMINATOR_MIN` | `25` | Minimum class-conditional called denominator for primary medians |
| `CI_ALPHA` | `0.05` | Two-sided 95% exact Clopper-Pearson intervals |
| `JOINT_L2_LAMBDA` | `0.05` | L2 penalty for the strictly convex monotone joint risk model |
| `TIE_RULE` | `fnv1a64(selective-tie/<seed>/<scheme>/<fold>/<sample_pos>)` | Deterministic, label-blind tie-breaking rank |

## Four auditable policies and Mondrian conformal calibration

For every learner (`scm`, `cart`, `rf`, `gbt`, `lasso`) on every outer fold `f`:

1. **`outcome_blind`:** Deterministic label-blind random thinning derived solely from `tie_rank`. Serves as the baseline comparator testing whether abstention beats random deferral.
2. **`score_margin`:** Margin `|s - 0.5|` calibrated on `s_inner[T_f]` (`t_lo = 0.5 - tau_margin`, `t_hi = 0.5 + tau_margin`), alongside class-conditional Mondrian conformal prediction sets on `s_inner[T_f]` with nonconformity `1 - s` for R and `s` for S. Conformal sets explicitly distinguish `{R, S}` (`conformal_both`, which covers the true label) from `{}` (`conformal_empty`, which does not cover the true label) and report `set_true_label_coverage`. Because `scm` produces binary `{0, 1}` scores (`|s - 0.5| == 0.5` everywhere), ties in `score_margin` fall through to the deterministic label-blind tie rule while `mash_distance` and `margin_plus_distance` provide non-degenerate abstention for `scm`.
3. **`mash_distance`:** Genomic-novelty gate scoring confidence as `-d_outer(u)` with distance threshold `tau_d` calibrated on `d_inner[T_f]`.
4. **`margin_plus_distance`:** Monotone L2-regularized logistic error-risk model `P(error | m, d) = sigmoid(b0 - w_m * z_m + w_d * z_d)` with `w_m >= 0, w_d >= 0` fit strictly on `(err_inner[i], |s_inner[i] - 0.5|, d_inner[i])` for `i in T_f` (standardization parameters `mean_m, sd_m, mean_d, sd_d` computed strictly on `T_f`; deterministic equal-weight fallback when `T_f` has zero inner errors or zero weight). Tests whether genomic distance adds value beyond score margin alone.

## Dual evaluation layers

Every policy reports two strictly separated layers in `grm-selective-v1` JSON:

- **`calibrated_operating_point`:** Applies thresholds `(t_lo, t_hi, tau_d, tau_joint)` fit strictly on `T_f` directly to outer test fold `E_f`. Under lineage or study shift, achieved test coverage can drift away from the nominal target `c`; this layer exposes that calibration shift directly. Outer null fits for `gbt` and `lasso` abstain with reason `null_fit` under score-dependent policies.
- **`matched_coverage`:** Retains exactly `floor(c * n_f)` of the same `n_f` outer-fold rows on every outer fold `f` for every policy using only that policy's outer confidence score and the frozen deterministic, label-blind tie rule. This removes coverage drift so policies are compared at identical call volume on identical folds.

## Clinical accounting and zero-ID hygiene

Each `(dataset, scheme, learner, policy, target_coverage)` cell and each Mondrian conformal point emits aggregate counts only (zero genome IDs, zero per-genome prediction arrays):

- **Raw counts:** `n`, `n_r`, `n_s`, `called`, `called_r`, `called_s`, `abstain_r`, `abstain_s`, `tp`, `fp`, `tn`, `fn`, plus `abstain_reasons` (`null_fit`, `ood_mash_distance`, `uncertain_score`, `conformal_both`, `conformal_empty`, `random_thinning`), summing to `abstain_r + abstain_s`.
- **Coverage:** `coverage = called / n`, `coverage_r = called_r / n_r`, `coverage_s = called_s / n_s`, and `set_true_label_coverage` for conformal sets.
- **Both denominators and predictive values with exact 95% Clopper-Pearson CIs:**
  - Among called: `vme_called = fn / (tp + fn)`, `me_called = fp / (tn + fp)`, `npv = tn / (tn + fn)`, `for_s = fn / (tn + fn)`, `ppv = tp / (tp + fp)`, `balanced_accuracy_called`.
  - Over all eligible: `vme_all = fn / n_r`, `me_all = fp / n_s`, plus per-100 clinical accounting (`caught`/`missed`/`deferred` per 100 R; `cleared`/`overcalled`/`deferred` per 100 S).
- **Sparse denominators:** Any rate with denominator `< 25` carries `sparse_denominator: true` and is excluded from primary species medians while raw counts stay recorded.

## CLI usage

Drawn `random` and `lineage` splits:

```bash
grm learn --species <slug> --exp-selective <out_dir>
```

Supplied `--fold-file` splits (`study-whole`, `time`, etc.):

```bash
grm learn --species <slug> --fold-file <plan.json> --fold-out <fold_out> --exp-selective <out_dir>
```

Writes `<out_dir>/<dataset>.<scheme>-s<seed>.selective.json`.

## Species-level cross-split aggregation

`scripts/selective/summarize.py` validates `grm-selective-v1` outputs and computes the pre-committed primary contrasts at `c = 0.85` (`score_margin` vs `outcome_blind` and `margin_plus_distance` vs `score_margin` on incremental reduction in `vme_called`):

```bash
python3 -I scripts/selective/summarize.py <selective_dir> --out <summary_dir>
```

- Excludes datasets with `called_r < 25` on either policy from primary species medians while retaining their raw counts in `datasets.tsv`.
- Aggregates paired differences to one value per species (median across eligible datasets in that species) and emits absolute species differences (`species_differences.tsv`), wins/ties/losses, and the one-sided exact binomial sign test on non-tied species (`sign_tests.tsv`). The sign test tests direction of paired differences across species, not mean magnitude.
- Requires passing on **both** `lineage` and `study-whole` splits (`joint_claims.tsv`) for any joint cross-split claim.
- Passing `--scan-species-json <species.json>` verifies zero-ID hygiene across inputs and outputs using the token rule of `scripts/c7-outside/scan.py`.
