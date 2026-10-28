# Lasso stopping fix and revision-8 migration prerequisites

Status: the solver change is approved and the code/test/checker prerequisites are implemented. **The large-species scan and the full 275-report revision-8 evaluation have not run.** No biological evaluation, report rebuild, deployment, or publication is authorized by this prerequisite lane.

The preserved independent verdict is [`data/research/grm-lasso-stopping-fix-20261027-r8kq3v/MIGRATION-REVIEW.md`](../data/research/grm-lasso-stopping-fix-20261027-r8kq3v/MIGRATION-REVIEW.md). It approved the solver change and required evaluator enforcement, a revision bump, schema/checker migration, edge guards, and stronger tests before any revision-8 run.

## Approved solver change

`src/learn/lasso.rs` formerly stopped coordinate descent after 300 sweeps without preserving a usable status. On ill-conditioned panels this stopped short of the optimum. The approved change retains 25 IRLS steps and raises the coordinate-descent limit to 100,000 sweeps per IRLS step (glmnet's default `maxit`). The objective, unpenalized intercept, standardization, class weights, penalty grid, tolerances (1e-6 sweep and 1e-5 IRLS), 1e-5 weight floor, warm starts, `max_nonzero`, and refit semantics are unchanged.

`Convergence { irls_steps, sweep_cap_hits, converged }` records how each penalty point stopped. A point converges only when:

- the IRLS change is below 1e-5;
- that step's coordinate descent settled rather than reaching the sweep cap; and
- the intercept and every coefficient are finite.

`fit_path_with_status` and `refit_path_with_status` return `LassoFit { model, convergence }`. The compatibility wrappers `fit_path` and `refit_path` retain their model-only signatures and report unconverged points on stderr; the evaluator no longer relies on those warnings.

### One-class and non-finite edges

Both public status-returning entry points return an **empty path by design** when the fitting mask contains only one class. One-class logistic regression has no finite intercept optimum; a long divergent path must not eventually claim convergence because of floating-point behavior. The legacy model-only wrappers inherit the same empty-path behavior. The evaluator's existing two-class outer-fold, inner-fold, and tuning-holdout gates remain intact, and one-class input is not converted into an ordinary fitted no-call.

A deliberate empty path on a valid two-class fit remains different: it means no penalty admits a feature and is the existing null-fit/no-call case. A non-empty path with no converged candidate is a solver failure, not a null fit.

The finite-parameter guard closes the separate NaN issue: `f64::max` can ignore NaN, so loop deltas alone are not sufficient evidence of convergence. NaN, positive infinity, and negative infinity in the intercept or any coefficient force `converged = false`.

## Evaluator enforcement in revision 8

`eval::fit_lasso` now uses `fit_path_with_status` and `refit_path_with_status`.

1. Every unconverged tuning-path point is excluded before validation loss and penalty selection.
2. A deliberately empty two-class path remains a null fit. A non-empty path with no converged candidate returns an error.
3. The selected point is defensively checked again for convergence.
4. The refit at the selected penalty is checked before any score or model is returned. An empty or unconverged selected refit is an error.
5. Path-selection and refit errors name the dataset, split scheme, zero-based outer fold, and stage. Final-model errors say `final model on all samples (no outer fold)` rather than pretending to be an outer fold.
6. Solver failure is never converted into a no-call, and an unconverged selected/final model is never scored or saved.

Converged candidates retain the existing minimum-validation-log-loss selection and refit behavior. Non-lasso learners are unchanged.

Any recorded unconverged two-class lasso path point makes `evaluation.valid` false, even when a different converged candidate could be selected. Missing per-fold/final solver evidence also makes a newly generated report invalid. Thus an excluded, non-selected failure remains visible and blocks a valid revision-8 report rather than being silently accepted.

## Report schema and independent checker

`REPORT_SCHEMA` remains **2**. `EVALUATION_REVISION` is now **8**. The optional Rust deserialization defaults preserve the ability to read older reports, but revision-8 output must contain the new evidence.

Each revision-8 lasso `LearnerResult` has `lasso_fits`, one record per outer fold. `FinalModels` has `models.lasso_fit` for the final fit. Each record has:

```json
{
  "path_points": 30,
  "path_unconverged": 0,
  "selected_point": 2,
  "selected_convergence": {
    "irls_steps": 3,
    "sweep_cap_hits": 0,
    "converged": true
  },
  "refit_convergence": {
    "irls_steps": 4,
    "sweep_cap_hits": 0,
    "converged": true
  }
}
```

For a deliberate empty path, `path_points` and `path_unconverged` are zero and `selected_point`, `selected_convergence`, and `refit_convergence` are null. Otherwise all three selected/refit fields are present and the selected index must identify a path point. `models.lasso` remains the saved model (or null for an empty path); the separate `models.lasso_fit` is its solver evidence.

`scripts/check-results.py` now:

- accepts revisions 3 through 8 and rejects future revisions;
- keeps the existing revision 3–7 schema, no-call, metric, split, and log semantics;
- requires per-fold and final lasso status on every available revision-8 split, including reports whose lasso rows have null fits;
- rejects missing fields, wrong types, invalid counts/indices, inconsistent selected/refit presence, and a final model inconsistent with an empty-path record; and
- requires zero unconverged two-class lasso path points and true selected/final convergence before accepting `valid: true`.

The checker tests include a valid revision-8 report, a valid empty-path final fit with a null model, future-revision rejection, and negative missing/malformed/unconverged cases, including a saved model whose record has no refit. Revision-7 reports do not need to be retrofitted with these fields.

### Stale and reuse behavior

`input_fingerprint` deliberately hashes inputs and settings, not code. The revision bump is therefore part of the same change as evaluator enforcement. A saved revision-7 report with an otherwise identical fingerprint is stale with reason `evaluation audit revision changed`; it cannot be reused inside a revision-8 report. Tests cover both `stale_reason` and the `analyze_species` reuse path. No forced partial real-data rebuild was run.

## Test coverage and required revision-8 gate

The saved cross-check panels are in `tests/data/lasso-stopping-fixtures.txt`. They preserve the 16 fixtures and candidate paths from the earlier cross-check, the best known objective at each point, and the 300-sweep trace. Bit-exact fixture checks remain pinned to the Linux platform because `powf` and `exp` can differ by one ULP on other platforms.

The lasso suite now has ten tests: the original seven plus public-wrapper, long-one-class, and finite-parameter regressions. In particular:

- `public_wrapper_uses_the_raised_cap_on_a_saved_panel` is a default, non-ignored test through `fit_path_with_status` on `streptococcus_agalactiae__clindamycin_lineage_fold0`. Its old 300-sweep path has unconverged points; reverting the inner cap to 300 fails this default test. Wrapper equivalence checks path length, penalties/models with `same()`, and statuses for both `fit_path` and `refit_path`.
- `one_class_long_paths_are_empty_by_design` exercises both public wrappers with 4,001 penalties for all-resistant and all-susceptible fitting masks.
- `non_finite_parameters_never_report_convergence` covers NaN and both infinities and exercises the solver with a NaN working weight.
- Evaluator tests prove that a lower-loss unconverged candidate cannot be selected, a path with no converged candidate cannot be refit, and an unconverged refit returns an error with dataset/scheme/fold/stage context. An unconverged or empty final-model refit names the final model instead of a fold.
- Analysis tests prove that unconverged or missing status cannot produce `valid: true` and that revision-7 results cannot be reused at revision 8.

The required release command, including both intentionally ignored C. jejuni tests, is:

```bash
cargo test --locked --offline --release --lib learn::lasso -- --include-ignored
```

The complete prerequisite gate is:

```bash
cargo fmt --check
cargo clippy --locked --offline --all-targets -- -D warnings
cargo test --locked --offline
cargo test --locked --offline --release --lib learn::lasso -- --include-ignored
(cd scripts/metric-audit/foldgen && cargo build --locked --offline --release)
python3 -I scripts/metric-audit/tests/test_check_results.py
python3 -I scripts/metric-audit/tests/test_audit.py
python3 -I scripts/metric-audit/tests/test_scripts.py
```

The Python suites call `scripts/metric-audit/foldgen/target/release/foldgen`, so build it first. No real or candidate data root is used by these tests.

The shared checkout was deleted at 04:41 UTC on 2026-10-28 while this work was in progress. Its edits were recovered into commits `a728d11` and `c87c4b5` and finished in the following commit. The pre-change files and hashes are under `data/research/grm-lasso-revision8-migration-20261028`. The completed gate's logs, mutation checks and validation summary are under `data/research/grm-lasso-revision8-finish-20261028`.

### Fixture accounting: 46 is not 38

Across all 91 saved points under the old limits:

| Measure | Count |
|---|---:|
| Sweep cap hit in some IRLS step | **46** |
| Sweep-cap only (not also out of IRLS steps) | **38** |
| Ran out of IRLS steps | 15 (7 are `one_class_train`) |
| Any cap | 53 |
| Flagged unconverged by the approved status rule | 26 |
| Never capped | 38 |
| Bit-identical when no earlier cap affected the warm start | 36 |

Thus **46** means the sweep cap was reached in at least one IRLS step; eight of those also exhausted IRLS steps. **38** is the sweep-cap-only count. Calling 46 “inner-only” or “sweep-cap only” is incorrect. The 26 flagged points comprise 15 that exhausted IRLS steps and 11 whose last inner solve hit the sweep cap.

Under the new limits, all 84 two-class saved points converge, with independently recomputed KKT residual at most 1e-6 and objective within 1e-6 of the best known. Four C. jejuni points still hit the 100,000-sweep cap in an intermediate IRLS step and settle later. The seven raw one-class fixture points remain unconverged in the internal trace used for accounting; public one-class calls now return empty paths instead of running that trace.

At the three prior C. jejuni failures, the old limits stop at 25 IRLS steps with 25 sweep-cap hits, KKT above 1e-4, and excess objective above 1e-4. The new limits converge in 7, 9, and 18 IRLS steps, which the ignored release test now pins. Refitting at the same penalties on the same mask returns the identical model and status.

## What revision 8 can change

The small deployment exercise (four datasets, 99–508 samples) reproduced every old-solver score/model on replay and found zero unconverged points among 2,107 after the fix. Thirty-four of 44 models changed numerically without a selected-penalty or call change, and lasso time increased from 15.0 s to 28.8 s (1.92×). This is evidence for the revision bump, not a prediction of full-run results or runtime.

A full revision-8 field-level diff must treat the following as liable to change:

- lasso scores, metrics, ROC, fold balanced accuracy, fold complexity, no-call counts, and timing;
- `models.lasso` and the new solver-status fields;
- `features` and annotations, because the list includes the top 30 lasso weights (newly included loci need annotation downloads);
- best-learner labels and corresponding log wording;
- `evaluation.revision`, report validity evidence, timing fields, and the report index.

SCM, CART, forest, boosting, MIC, catalog, and probe outputs should remain identical. The earlier cross-check reproduced all 278 non-lasso rows; the full-run diff must enforce that rather than assume it. Boosting is unaffected by this solver fix. XGBoost and the default `colsample = 0.5` remain a separate verification gap and do not block this revision-8 prerequisite work.

## Runtime stop rule and still-unrun biological gates

The new limits permit 25 × 100,000 sweeps per point, or 75 million sweeps for a 30-point path: 333× the old theoretical worst case. The observed 1.92× runtime covers only four small datasets, and intermediate cap hits persist on the motivating C. jejuni panels. There is no tested large-species runtime upper bound.

Before the full revision-8 run, a separate authorized lane must use the existing `lasso_deploy` harness for a lasso-only scan over the largest datasets, explicitly including **M. tuberculosis, E. coli, and K. pneumoniae**, and record sweeps, seconds, and status for every fit.

That harness is not in the repository. It is `scripts/lasso_deploy.rs` in the `grm-lasso-stopping-fix-20261027-r8kq3v` evidence workspace on Vast, which reads the candidate data root there. When these prerequisites were finished, the local copy of that evidence had not been restored and no local data root held the large-species inputs the scan needs. The harness replays the revision-7 `fit_lasso` against old and fixed solver modules. For the scan it has to be rebuilt against revision-8 source so that each fit goes through `fit_path_with_status` and `refit_path_with_status` and stops on the same conditions as the evaluator.

The scan must pre-register its stop rule in its evidence before launch: a per-fit wall-clock budget, a total wall-clock budget, and immediate stop on an unconverged selected/refit two-class fit. Reaching a budget or solver-failure condition stops the scan, preserves the partial logs, and blocks the full run pending review. No numeric time budget is approved by this prerequisite lane; it must be chosen from the operational constraints rather than extrapolated uncritically from the 1.92× small-panel figure.

The following remain separate, gated steps:

1. Run and review the large-species lasso-only scan with that stop rule.
2. Only after it passes, run the full 275-report evaluation with revision-8 source into a **new data root**, leaving revision-7 candidate, backup, and published roots untouched.
3. Require every report to pass the revision-8 checker with zero unconverged two-class lasso fits, then perform the field-level migration diff and completion-log audit.
4. Recheck/rerun downstream consumers as listed below, and mechanically audit which lanes actually consume lasso output before publication.
5. Rebuild and review the report/web bundle. Publication is a further explicit decision.

## Deployment settings and downstream work

The current published release remains revision 7. `scripts/deployment/` retains its revision-7 defaults and 237-dataset published-release expectations; this lane deliberately did not edit it while its independent review is active.

At a future, authorized revision-8 staging/publishing operation, the operator must explicitly set:

```bash
GRM_EXPECT_REVISION=8
GRM_EXPECT_DATASETS=275
```

Those are operation-time requirements, not changes to the current defaults and not evidence that staging or publication has occurred.

### Lanes that refit lasso through `grm` code

These need revision-8-source reruns when their inputs changed:

- `grm-research-selective-calling-r1`
- `grm-research-temporal-study-splits-rev7cand-r1`
- `grm-research-provenance-dup-all23-r1`
- `grm-research-allmic-cohort-r1`
- `grm-research-mtb-scaleup-rev7cand-r1`

### Lanes that read saved scores or best-learner labels

These need rechecks against revision-8 reports:

- `grm-all23-saved-intervals-20261027-kimi`
- `grm-all23-selective-20261027-kimi`
- `grm-all23-mic-baseline-20261027-kimi`
- `grm-research-clindamycin-imlsb-r1`
- `grm-clindamycin-label-recheck-20261027-tuwf2end`
- `coresistance-time-all23-20261027T101316Z`
- `grm-research-deescalation-r1`
- `deescalation-all23-saved-20261027-vozjjhqn`

The last lane uses GBT as its pre-specified learner, so its required lasso-specific recheck is the 544-row lasso no-call audit table. Which lanes consume lasso output rather than merely mention it still requires a mechanical audit before publication; that audit does not block the full revision-8 run itself.

Older revision-5 lanes remain historical. AMRFinder/GWAS work reads the revision-5 snapshot and is not on this migration list; this does not authorize restarting hard-stopped work.

**Completion of these prerequisites does not authorize the large-species scan, the full revision-8 run, downstream reruns, staging, or publication.**
