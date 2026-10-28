# Lasso stopping fix (2026-10-27)

`src/learn/lasso.rs` stopped its inner coordinate-descent loop after 300 sweeps, without reporting it. On ill-conditioned panels the cap ended the solve well short of the optimum. This change raises that cap to 100,000 sweeps (glmnet's default `maxit`) and reports every penalty point that stops at a cap. Nothing else changes: the objective, unpenalized intercept, standardization, class weights, penalty grid, tolerances (1e-6 sweep, 1e-5 IRLS), the 25-step IRLS cap, the 1e-5 IRLS weight floor, warm starts, `max_nonzero` and refit semantics are all as before.

**This is an evaluation behavior change.** Lasso scores and final lasso models in every saved report were computed with the old cap, so it needs a new `EVALUATION_REVISION`. This lane did not bump the revision or rerun the pipeline. See [Revision and migration](#revision-and-migration).

The defect was confirmed in `data/research/grm-learner-lasso-gbt-resolve-20261027-d0a401022d9d/VERDICT.md`. The verification for this fix is in the Vast workspace `/root/grm-research/grm-lasso-stopping-fix-20261027-r8kq3v`.

## What changed in `lasso.rs`

- `LIMITS` sets 25 IRLS steps and 100,000 sweeps per step (previously the literals 25 and 300).
- `Convergence { irls_steps, sweep_cap_hits, converged }` records how each penalty point stopped. `converged` is true only when the IRLS change fell below its tolerance *and* that step's coordinate descent settled. A point that ran out of IRLS steps, or whose last inner solve stopped at the sweep cap, is `converged: false`.
- `fit_path_with_status` and `refit_path_with_status` return `LassoFit { model, convergence }` for each point.
- `fit_path` and `refit_path` keep their signatures and return the same models. They print one stderr line for each unconverged point (`lasso fit_path: point i (lambda …) stopped at a cap unconverged: …`), so a capped fit can no longer pass silently.
- `LassoModel` and its serialized form are unchanged. The evaluator (`eval.rs`) still calls `fit_path`/`refit_path`, so **reports do not yet record the flag**; see [Remaining gaps](#remaining-gaps).

## Tests

The saved cross-check panels are in `tests/data/lasso-stopping-fixtures.txt` (sha256 `7f354884…9a1f99550`). That file holds the 16 fixtures of `data/research/grm-learner-crosscheck-rev7-20261027-zz1lfkza/fixtures`, each hash-checked against that workspace's `SHA256SUMS` when the file was written. It also holds each saved candidate path, plus the best objective any solver reached at each point and the 300-sweep trace, both from the resolution lane's `results/lasso-resolution.csv`.

| Test (`learn::lasso::tests`) | Default | Checks |
|---|---|---|
| `a_cap_that_ends_the_solve_is_flagged` | yes | A 1-step/1-sweep limit is flagged (`converged: false`, one IRLS step, one sweep-cap hit); the default limits converge; `fit_path` returns the same models as `fit_path_with_status`. |
| `saved_panels_converge_flag_their_caps_and_keep_clean_paths` | yes | Every check below, on the 14 panels other than C. jejuni (77 points). |
| `all_91_saved_points` | `#[ignore]` | All 16 panels, 91 points. |
| `campylobacter_failures_reach_the_reference_optimum` | `#[ignore]` | The three failing C. jejuni ciprofloxacin points (random 5 and 6, lineage 6, zero-based). |

Both panel tests run each fixture's path with the old limits (25 IRLS steps, 300 sweeps) and the new ones, and check every saved point:

- With the old limits the path reproduces the saved candidate path bit for bit, with the same IRLS steps and sweep-cap hits as the saved trace. Every point that ran out of IRLS steps is flagged, and every flagged point stopped at a cap.
- With the new limits, every two-class point is `converged`, its KKT residual is at most 1e-6 and its objective is within 1e-6 of the best known. KKT and objective are recomputed in the test independently of the solver: the standardized-scale stationarity violation over the intercept and every fittable coefficient.
- `one_class_train` (one training class) is flagged unconverged at all 7 points under both limits. Its intercept diverges, so no cap can converge it; the evaluator rejects one-class folds before fitting.
- A point is unchanged when no earlier point on its path, and not the point itself, stopped at a cap under the old limits (warm starts carry each point into the next). Every unchanged point has a bit-identical model and the same status under both limits.

Pinned counts over all 91 points:

| | Count |
|---|---|
| Old limits: sweep cap reached in some IRLS step | 46 |
| Old limits: ran out of IRLS steps | 15 (7 are `one_class_train`) |
| Old limits: any cap reached | 53 |
| Old limits: flagged unconverged | 26: those 15, plus 11 where the IRLS change met its tolerance while the last inner solve stopped at 300 sweeps |
| New limits: flagged unconverged | 7, all `one_class_train` |
| New limits: sweep cap still reached in an intermediate IRLS step | 4 C. jejuni points (random 5 and 6: 6 hits each; lineage 5: 1; lineage 6: 3), all settling in a later step |
| Unchanged, bit for bit | 36 |

At the three prior failures, the old limits stop at 25 IRLS steps with 25 sweep-cap hits, KKT above 1e-4 and excess objective above 1e-4. The new limits converge in 7, 9 and 18 IRLS steps, with KKT ≤ 1e-6 and objective within 1e-6 of the reference optimum. A refit at the same penalties on the same mask returns the identical model and status.

These figures correct two statements in the earlier verdict. "The inner cap binds at 53 of 91 points" counts any cap; the sweep cap alone binds at 46. And "all 91 converge" holds for the 84 two-class points only; the 7 one-class points cannot converge.

The C. jejuni panels need about 2.3 million sweeps under the new limits: seconds in release, but many minutes unoptimized. Their tests are therefore ignored by default. Run everything with:

```bash
cargo test --locked --offline --release --lib learn::lasso -- --include-ignored
```

## Verification

All runs used Rust 1.98.1 with `--locked --offline`. On Vast, the source was a scratch copy of the candidate tree (`/root/grm-extra-species-analysis-20261027/source`) with only `src/learn/lasso.rs` replaced. The candidate's own `lasso.rs` (sha256 `12dd174c…`) was kept beside it as `lasso_old.rs`, and a copy of the fixed file with the sweep cap set to 300 was used only to read the old flags. The candidate source and the earlier cross-check evidence hashes were checked before and after every Vast step, and none changed.

### Tests

The final `lasso.rs` is sha256 `8394fa1f…b4b75c0683`. The deployment and cross-check runs below used an earlier revision of the file whose solver code is identical; only a test helper changed afterwards, for clippy.

- Local WSL, shared checkout: `cargo fmt --check` and `cargo clippy --all-targets -- -D warnings` are clean. `cargo test --locked --offline`: 205 + 19 + 1 + 1 passed, 0 failed (lib, bin, and the two integration targets). 8 tests are ignored: the two lasso tests above and 6 existing benchmark or full-run tests. `learn::eval::tests`: 21/21. The lasso tests take 44 s in debug. The checkout also had another lane's uncommitted `src/matrix.rs` (sha256 `26a65319…`), which this lane did not touch.
- Vast, scratch source, release: `learn::lasso` with `--include-ignored` gives 7/7 (79 s), and `learn::eval::tests` gives 21/21.

### Original 16-case cross-check, rerun on the fixed source

`xcheck/` reruns the original cross-check (`run.py`, `lasso_adjudicate.py`, `audit.py`, probe and plan) unchanged. Only two lines differ: the probe now links the fixed scratch source, and `references.py` points at the new workspace.

- 362 comparison rows. All 278 non-lasso rows (SCM, CART, forest, boosting formula checks, partitions, metrics) are identical to the original run. All 84 lasso rows differ, as expected.
- Candidate objective more than 2e-5 above the long SciPy reference: **3 → 0** of 84 penalty points. Candidate worse than any reference: 3 → 0.
- Objective outside tolerance in either direction: 33 → 30. The remaining 30 are scikit-learn reference fits with convergence warnings (30 before and after). In those rows the reference is the worse side.
- The adjudicator's largest candidate KKT is 1.1e-3, at `marker_dup_constant` point 6. That KKT also counts columns that are constant on the fitting samples, which the solver clears by design. At that point the candidate's objective is 5.8e-4 *below* the SciPy reference. The tests' KKT excludes those columns, and it is at most 1e-6 at every two-class point.

### Deployment-scale penalty selection and refit

`lasso_deploy` rebuilds each dataset exactly as `analysis::analyze_dataset` does: cohort, `Dataset::project` at `min_support = max(n/200, 3)`, lineage groups, and both fold schemes. It then replays `eval::fit_lasso` for each of the 5 outer folds and for the final model, once with the old solver and once with the fixed one. Each replay uses `training_candidates`, the 0.2 tuning holdout, 4,000 chi-square-screened features, 30 penalties, `max_nonzero` 400, minimum-validation-log-loss selection, and a refit on the outer training set. That makes 44 fits over four datasets:

| Dataset | Samples (R) | Features | Data root | Max held-out score shift |
|---|---|---|---|---|
| C. jejuni ciprofloxacin | 508 (113) | 266,769 | revision-7 preflight | 2.9e-4 |
| E. cloacae tobramycin | 99 (45) | 438,451 | revision-7 preflight | 2.9e-4 |
| S. haemolyticus fusidic acid | 114 (39) | 37,848 | all-23 candidate | 1.7e-5 |
| S. agalactiae clindamycin | 141 (115) | 8,014 | all-23 candidate | 2.0e-6 |

The two preflight reports are byte-identical to the candidate's (sha256 checked).

- **Replay fidelity.** With the old solver, the replay reproduces every saved lasso held-out score bit for bit (0 mismatches over 40 fold fits) and all 4 saved final models. The 300-sweep copy matches the candidate file on all 44 fits.
- **Convergence.** Under the old solver, 300 of 2,107 path and refit points stop at a cap unconverged. Up to 18 of 30 path points fail in a single fit, and up to 16 of 30 refit points. Under the fix, 0 of 2,107 are unconverged and none reaches the sweep cap.
- **Selection and calls.** 34 of 44 models change numerically. No selected penalty changes (0/44). There are 0 held-out call changes and 0 no-call changes. One CV refit gains a weight (C. jejuni random fold 3, 16 → 17 non-zero). All four final models keep the same feature sets.
- **Runtime.** Total lasso fitting time goes from 15.0 s to 28.8 s (1.9×). The largest single increases are 1.2 → 3.5 s (C. jejuni final model) and 0.23 → 0.93 s (S. haemolyticus lineage fold 3).

### Boosted trees on a separate early-stopping holdout

XGBoost cannot be imported in any Vast Python (`gwas`, `tbp6`, system) or in local WSL `python3`, and no XGBoost wheel or library exists on the Vast disk. Nothing was downloaded, so **XGBoost remains unverified**.

The substitute check uses scikit-learn 1.3.2's `TreeGrower`, as in the resolution lane. It runs on deployment-scale inputs: outer fold 0 of each dataset and scheme, 4,000 screened features, the evaluator's separate 0.15 early-stopping holdout, 400 rounds, patience 25, and `colsample = 1.0`. With every feature offered, no column-sampling random numbers are drawn. That gives 8 fits. The loop, gradients and early-stopping rule are re-implemented in the script; only tree growth is external. Columns are in matrix-row order, so that both implementations break exact ties by the same index.

- In all 8 fits, the external grower agrees tree for tree, to within 1.3e-8 in leaf values, up to its first divergence. Every first divergence is at a split whose two candidate gains are equal under correctly rounded float64 sums (relative gap 0 to 1.0e-11), which float32 histograms cannot separate. One fit (S. agalactiae lineage) has no divergence.
- Kept-tree counts match in 4 of 8 fits and differ by 1 to 6 in the others.
- There are **0 call differences** on held-out test samples and on the early-stopping holdout, in all 8 fits.

Before the columns were ordered by row, the two implementations broke ties differently, and 4 test calls and 1 holdout call differed. Those were artefacts of the harness, not of either implementation.

## Revision and migration

**Nothing in this lane bumps the revision, rewrites a report, or touches the candidate or published data.** What follows is the plan for a separate, independently reviewed re-evaluation.

**Why a revision is required.** `input_fingerprint` (`src/analysis.rs`) hashes inputs and settings, not code: "Code changes need a new EVALUATION_REVISION". The fix leaves every `EvalConfig` field unchanged, so fingerprints do not move. Without a bump, `grm learn` would treat every saved report as current and keep its old lasso results. The re-evaluation lane must set `EVALUATION_REVISION` to 8. `stale_reason` will then report "evaluation audit revision changed" for every saved report.

**Affected reports: all of them.** Every report has a lasso row in both schemes and a final lasso model, so this covers:

- all 275 reports of the all-23 candidate (`/root/grm-extra-species-analysis-20261027/data`, backup `data/backups/revision7-all23-candidate-20261027`);
- the 237 published revision-7 reports. The candidate's 13 original species are byte-identical to that release (`MANIFEST.json`), so **the published revision-7 release carries the same 300-sweep lasso results**. It stays as published until a reviewed revision 8 replaces it.

The 275-report candidate is **not publication-ready** after this fix.

**What changes in a report.**

- Changes: in each lasso learner row, scores, metrics, ROC, fold balanced accuracy, fold complexity and timing; `models.lasso`; the `features` list and its annotations, through `model_features`, which takes the top 30 lasso weights; any `best_learner` that the lasso wins or loses; `evaluation.revision`; and timing fields.
- Should not change: SCM, CART, forest, boosting, MIC (boosting-based), catalog and probe outputs. Their code is untouched and seeded, and the cross-check rerun reproduced every non-lasso comparison. The re-evaluation should *verify* this: diff each revision-8 report against its candidate, and require that only lasso-derived fields, the revision and timings differ.

**Expected size.** The deployment exercise covers 4 small datasets (99 to 508 samples). There, held-out lasso scores move by at most 2.9e-4, with no call, no-call, penalty-selection or final-feature-set change. Larger datasets were not exercised, so treat calls, metrics and best-learner labels as able to change until the rerun shows otherwise. Expect the lasso part of `learn` to take roughly twice as long.

**Suggested migration.**

1. Review this change independently.
2. In a new lane, bump `EVALUATION_REVISION` to 8. Decide whether the lasso status should be persisted in reports (see gaps); that is a schema decision, and it belongs before the run.
3. Run `grm learn` for all 23 species into a **new** data root, copied from the candidate inputs. Leave the candidate root, both local backups and the published release untouched.
4. Audit every report: the field-level diff above, `check-results.py`, completion logs, and a scan of the logs for `lasso … stopped at a cap unconverged` lines. Any such line on a two-class fold is a finding.
5. Rebuild the report and web bundle from revision 8. Do not publish before review.

**Downstream research that reads saved lasso outputs or `best_learner` must be rechecked against revision 8 and rerun where its inputs changed.** Lanes under `data/research` that mention the lasso:

- `grm-research-temporal-study-splits-rev7cand-r1`
- `grm-research-provenance-dup-all23-r1`
- `grm-research-allmic-cohort-r1`
- `grm-all23-selective-20261027-kimi`
- `grm-all23-saved-intervals-20261027-kimi`
- `grm-research-selective-calling-r1`
- `grm-research-mtb-scaleup-rev7cand-r1`
- `grm-research-deescalation-r1`
- `grm-research-clindamycin-imlsb-r1`
- `grm-clindamycin-label-recheck-20261027-tuwf2end`
- `grm-all23-mic-baseline-20261027-kimi`
- `coresistance-time-all23-20261027T101316Z`

Older revision-5 lanes, such as `grm-research-temporal-study-splits-r1`, are historical and need no rerun. Which lanes actually consume lasso scores, rather than only mentioning the learner, was not audited here. Lanes that refit their own lasso through `grm` code paths (selective and temporal reruns, for example) need rerunning on revision-8 source in any case.

## Remaining gaps

- **The flag does not reach reports.** `eval.rs` was outside this lane, so evaluation still calls `fit_path`/`refit_path`. A capped fit is reported on stderr, but `LearnerResult` and `FinalModels` do not record it. Before revision 8, switch `fit_lasso` to the `*_with_status` functions. Then either record the status (a report schema change) or refuse an unconverged two-class fit.
- **Scale coverage.** Deployment-scale behavior was exercised on 4 datasets of 99 to 508 samples only. No large species was run (M. tuberculosis, E. coli, K. pneumoniae), and there was no lasso-only scan of all 275 datasets. `lasso_deploy` accepts any dataset ids and would do that scan, at roughly the cost of the lasso part of `learn`.
- **The sweep cap can still bind.** At 4 C. jejuni points, 100,000 sweeps are reached in an early IRLS step, and a later step settles. A larger or worse-conditioned panel could reach it in the final step. That would now be flagged, not hidden, but this lane has no further remedy.
- **Boosting.** XGBoost is unverified. The deployment default `colsample = 0.5` (a ChaCha8 stream no external library replays) was not compared externally. The boosting loop in the comparison is re-implemented; only tree growth is external.
- **Default test coverage.** The two C. jejuni tests are `#[ignore]`d in routine `cargo test` and need the release command above.
- **Evidence location.** Results are in the Vast workspace `/root/grm-research/grm-lasso-stopping-fix-20261027-r8kq3v`: `results/deploy*`, `results/gbt-holdout.json`, `xcheck/results`, `logs/` and the harness in `scripts/`. They are not copied under `data/research`, which this lane did not write.
