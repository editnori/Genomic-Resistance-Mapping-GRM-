# Lasso stopping fix and revision-8 migration prerequisites

Status: the solver change is approved, the code/test/checker prerequisites are implemented and independently reviewed, and the review's findings are closed. The large-species scan harness is prepared and its stop rule is registered below. The first registered scan stopped on its first fit, and the coordinate-descent solver it ran has since been replaced by an exact active-set solver (below), which still needs its independent review. **The rerun of the large-species scan and the full 275-report revision-8 evaluation have not run.** No scan, biological evaluation, report rebuild, deployment, or publication is authorized by this document.

The blocking review that set these prerequisites wrote its verdict to `data/research/grm-lasso-stopping-fix-20261027-r8kq3v/MIGRATION-REVIEW.md`. That file was lost in the 2026-10-28 deletion and is not on Vast either, so its findings survive only in the review's thread and in this summary: it approved the solver change and required evaluator enforcement, a revision bump, schema/checker migration, edge guards, and stronger tests before any revision-8 run. The guards were then reviewed independently; that review, its mutation logs and the revision-7 timing tables it measured are in [`data/research/grm-lasso-r8-review-20261028/`](../data/research/grm-lasso-r8-review-20261028/README.md).

## Approved solver change

`src/learn/lasso.rs` formerly stopped coordinate descent after 300 sweeps without preserving a usable status. On ill-conditioned panels this stopped short of the optimum. Revision 8 first kept 25 IRLS steps and raised the coordinate-descent limit to 100,000 sweeps per IRLS step (glmnet's default `maxit`), with a 1e-6 coefficient-change test in every step.

### Why coordinate descent was replaced

The registered scan ran that solver and stopped on its first fit on 2026-10-28: `mycobacterium_tuberculosis__delamanid random fold 0 passed its 842 s budget` at 845 s. The diagnosis is in [`data/research/grm-lasso-r8-slowdown-20261028/`](../data/research/grm-lasso-r8-slowdown-20261028/README.md).

On fold 0's 798 fitting samples, 2,820 of the 4,000 screened k-mer columns have a twin that differs in at most one sample. Moving weight between such twins barely changes the objective, so coordinate descent kept making small coefficient moves for tens of thousands of sweeps after the objective had stopped improving. At path point 9 it took 330,134 sweeps and 76.5 s. Points 20 to 26 used all 100,000 sweeps in every IRLS step and were flagged unconverged although they were within 2e-8 of the certified optimum with KKT residuals of about 1e-7, so the scan would have stopped even without a time budget. Revision 8 did at least 108 times revision 7's coordinate-descent work on that fold.

### The current solver

IRLS is unchanged: the objective (class-weighted mean log loss plus lambda times the standard-deviation-weighted L1 norm), the unpenalized intercept, standardization, class weights, the penalty grid, the 1e-5 floor on the IRLS weight p(1 - p), warm starts down the path, the evaluator's selection and refit semantics. What changed is how each IRLS step's weighted least-squares problem is solved and how a point is judged converged.

Each IRLS step is solved exactly by feature-sign search (Lee et al. 2007), as the diagnosis prototyped in numpy (`scripts/asnewton.py` there):

1. The non-zero coefficients form the active set. With each one's sign fixed, the least-squares objective is a quadratic, and its minimizer is one Cholesky solve on the active columns' weighted Gram matrix, intercept included.
2. A line search from the current point toward that minimizer evaluates the true least-squares objective at every point where a non-zero coefficient reaches zero, and at the full step, and moves to the lowest. A coefficient that reaches zero leaves the active set.
3. Once every active coefficient is stationary, the zero coefficients whose condition |g_k| <= lambda sd_k fails by more than 1e-9 (on the standardized scale) enter at zero with the sign of their gradient, worst first, in batches of at most one more than a tenth of the active set's size. If such a batch's step does not lower the objective, only its worst column is kept: from a stationary point, one entering column's step always does.
4. A column that is, to rounding, a linear combination of the intercept and the active columns cannot join the factor (its Cholesky pivot is below 1e-10 of its weighted sum of squares). k-mer columns can be exactly dependent, for example one column the union of two disjoint others. The solver then moves that column's coefficient and the active ones along the direction that leaves the fit unchanged, which lowers the objective at the rate pen_k - |g_k|, until an active coefficient reaches zero and leaves; the new column then joins in its place.
5. The solve ends when no zero coefficient violates its condition by more than 1e-9.

The Gram matrix and its Cholesky factor change by one row as a column enters (one forward solve) or leaves (Givens rotations on the rows after it), so a solve costs one factorization of the warm-started active set per IRLS step plus quadratic-cost updates, not a refactorization per change. The weighted Gram products are dense dot products over the fitting samples.

After the solve, the IRLS step is halved, up to 30 times, while it raises the penalized objective; near separation a full Newton step can overshoot. A step that no halving makes lower, or a solve that returns the current point, ends the point unconverged.

### Revision 9: a step that rounding hides

The revision-8 run stopped twice on the same solver weakness, each time with an invalid report: `vibrio_cholerae__nalidixic_acid` random fold 0, tuning path point 13 (100 IRLS steps, 99 of them with a solve at its 5,000-iteration guard, KKT 6.1e-4), and `salmonella_enterica__azithromycin` random fold 3, point 19 (3 IRLS steps, 2 solves stopped short, then no halving lowered the objective, KKT 2.8e-3). Every other point of the run converged.

Both started with an active coefficient within rounding of zero. In V. cholerae, near-copies 2209, 2389 and 2401 (they differ in 2 of the 57 fitting samples) reached zero almost together in one step; only 2389 was set to zero, and 2209 was left at 2.4e-15. Every later sign-fixed step crossed 2209 at t ≈ 9e-14, where the objective's fall is below rounding, so every candidate step showed a change of zero or more (+1.8e-18) and none was taken. The solve then took the point for stationary and checked the zero columns there. They violated their conditions by up to 3.0e-4 because the rest of the active set still had far to go. In V. cholerae the worst was 2985, the exact complement of 2209 on the fitting samples, so `exchange` swapped the two on a rate of ±1e-18 (pure rounding), back and forth until the guard. In S. enterica (a crossing at t ≈ 1.6e-15, change +3.0e-19) each violator entered, its step was blocked by the same crossing, and it was blocked in turn, 31 of them, many near-copies of active columns.

Exact complements were not the cause. Fold 0's screened set has 10,397 exact complement pairs (the independent 275 review's count), and every other point converges with them. Removing 2985 and 3174 from it (which is what folding complements into `representatives` would amount to) leaves point 13 failing exactly as before. Refusing the exchange, or detecting the 2-cycle, stops the swapping but leaves the solve stuck at the same point.

The fix is in the line search. Every non-zero active coefficient carries its fixed sign. When each entering column's step moves in the direction of its own sign, the objective along the sign-fixed step equals the sign-fixed quadratic up to the first crossing, and that quadratic falls all the way to its minimizer at t = 1. So the first crossing always lowers the objective in exact arithmetic. When no candidate step shows a change below zero but there is a crossing before t = 1, the solve now steps to the first crossing anyway, and that coefficient leaves. Three smaller rules come with it, as the admin ruled:

- An exchange must strictly improve: its rate has to be below zero by more than `INNER_TOL` (1e-9 on the standardized scale), the margin a zero column needs to enter at all. An exact complement of an active column has a rate of exactly zero, so rounding alone never decides it.
- Tied violators go to the lowest column. Revision 8 already ordered them so, through a stable sort of a list built in column order; the comparison now states it, so the order no longer rests on how the list is built.
- Every active coefficient that reaches zero in an exchange, to rounding, leaves, ties included, so no coefficient stays active at zero or past it.

On the fitting sets saved from the two stops, V. cholerae point 13 now converges in 3 IRLS steps at KKT 2.2e-9 and S. enterica point 19 in 3 steps at KKT 8.2e-8, with no solve stopping short. Over all 22 fits of the two datasets (660 points) the new step fired at 4 points and the margin refusal never fired. Every point whose path never reached a firing is bit-identical to revision 8, and so is every point of the 16 saved fixture panels, where nothing fired. Downstream of a firing point the objective moves by at most 6.3e-14 and fitted probabilities on the fitting samples by at most 1.6e-6. Held-out scores can move more, up to 0.098 at V. cholerae fold 0 points 17 to 29: there the optimum puts weight on columns that are exact complements on the fitting samples but differ outside them, any split of that weight is optimal, and revision 8's split depended on its history. Revision 9's split is deterministic, the same on every run and at any thread count, but not canonical: like revision 8's, it follows the path the solve takes, not a rule such as the collapse of exact copies.

### Convergence test

A point converges when its KKT residual is at most 1e-6, on a two-class fit, with a finite intercept and coefficients. The KKT residual is the largest violation of the penalized objective's optimality conditions on the standardized scale: the intercept's gradient, and for each column with a carrier, with g its gradient over its standard deviation, |g + lambda sign(b)| when its coefficient b is non-zero and max(|g| - lambda, 0) when it is zero. This is the residual the fixture tests have always computed independently, at the same 1e-6 acceptance. It is checked before each IRLS step, so a warm start that already meets it takes 0 steps (point 0, at the largest penalty, always does). A NaN anywhere in it counts as a failure, since `f64::max` would drop it.

The two-class condition matters only for internal calls: one-class logistic regression has no finite optimum and its gradient vanishes only as the intercept diverges, so a small residual certifies nothing. The public entry points return an empty one-class path before solving (below).

Two guards bound a point: 100 IRLS steps, and 5,000 active-set iterations per least-squares solve. Neither is the stopping rule. A point that reaches the IRLS guard is unconverged. An IRLS step whose solve reaches its iteration guard, or where rounding leaves a column no descent, is counted in `sweep_cap_hits`; the point still converges only if its KKT residual later meets the tolerance. On the saved fixture panels and on MTB delamanid fold 0 neither guard is reached and no solve stops short.

`Convergence { irls_steps, sweep_cap_hits, converged }` keeps its fields and serialized names, so the report schema and the checker are unchanged. `sweep_cap_hits` now counts IRLS steps whose active-set solve stopped short; the name dates from the coordinate-descent sweep cap. `fit_path_with_status` and `refit_path_with_status` return `LassoFit { model, convergence }` and take an observer that sees each point as soon as it is solved (`PointProgress`: point, lambda, IRLS steps, active coefficients, KKT residual, seconds, converged). The compatibility wrappers `fit_path` and `refit_path` retain their model-only signatures and report unconverged points on stderr; the evaluator no longer relies on those warnings.

### Exact duplicate columns

Columns identical on the fitting samples are collapsed: only the lowest-index column of each group can enter the active set, so the group shares one coefficient. Copies have the same standard deviation and so the same penalty factor, so every same-sign split of that weight among them has the same loss and penalty, and the optimal objective does not change. Coordinate descent split it arbitrarily, depending on its history. The KKT test still covers every copy. The Cholesky pivot test would also refuse a copy of an active column, but only through rounding; the collapse makes the lowest index a rule and skips the copies in each gradient pass (410 of fold 0's 4,000 columns).

The copies can differ on samples outside the fitting mask, so which copy carries the weight can change held-out scores: on MTB delamanid fold 0 the 410 redundant copies are identical on the 798 fitting samples and differ only on the 199 tuning-holdout samples. Keeping the weight on the lowest index makes that choice deterministic. `training_candidates` already removes patterns duplicated on the outer-training samples, so the refit and the final model have no exact copies; only the tuning path does.

### max_nonzero

The path still stops after the first point with more than `max_nonzero` (400) non-zero coefficients, and the count is still the saved model's weight count. Because exact copies share one coefficient, that count is now of distinct columns. Coordinate descent counted whichever copies its history had given weight; counting every copy of an active column instead would stop the path on columns the model does not use. On the diagnosis's final MTB delamanid fit, the prototype's distinct count stopped the path after 19 points where revision 7 stopped after 16, and both selected point 3.

### Costs

On MTB delamanid random fold 0, the fit that stopped the scan, the 30-point tuning path took 13.4 s and 11.5 s in two release-build runs on the diagnosis machine (i5-10400, 12 threads, with other lanes holding the load average at 13 to 21), and the refit 0.7 s and 1.0 s. Each point took 0 to 4 IRLS steps and each least-squares solve at most 44 active-set iterations. On the same machine revision 7's capped coordinate descent took 51 to 153 s for this path depending on load, the numpy prototype 35 to 67 s, and revision 8's coordinate descent was stopped at its 2-hour limit after 7,166 s at point 27.

Every point is within 8.4e-9 above the diagnosis's certified optimum (point 28, KKT 8.6e-7). Revision 8's coordinate descent was within 1.7e-8 at the points it finished, and revision 7 was up to 1.9e-4 above. The selection, refit and calls match revision 7: point 3, 53 refit weights, and none of the 251 test calls changed (largest score change 1.7e-5, test log loss 0.68618 for both).

| Point | Revision 7 above optimum | Revision 8 above optimum | Active-set solver above optimum | Active-set KKT |
|---:|---:|---:|---:|---:|
| 3 | 3.5e-10 | 2.0e-12 | -2.2e-16 | 1.1e-9 |
| 9 | 1.3e-5 | 1.4e-11 | 3.2e-13 | 1.7e-7 |
| 20 | 1.9e-4 | 7.0e-9, unconverged | 7.0e-11 | 1.1e-7 |
| 26 | not reached | 8.1e-9, unconverged | 4.3e-10 | 3.7e-7 |
| 28 | not reached | not reached | 8.4e-9 | 8.6e-7 |

Revision 7's and revision 8's columns are from the diagnosis's per-point table; revision 7's path stopped after point 20 and revision 8's run was killed during point 27.

Each IRLS step costs a pass over the screened columns' carriers for the KKT residual (2.36 million on fold 0), the warm-started active set's Gram matrix and factor (about n m²/2 and m³/6 multiply-adds for n fitting samples and m active columns: 53 million and 8 million at fold 0's largest set, 366 columns), and per active-set change one Gram row (n m), a forward solve (m²/2) and a gradient pass for the next check. Cost grows with the active set and the number of IRLS steps. Near twins cost no more than distinct columns, because each solve is exact whatever their correlation. Only fold 0 has been measured; the scan measures the rest.

### One-class and non-finite edges

Both public status-returning entry points return an **empty path by design** when the fitting mask contains only one class. One-class logistic regression has no finite intercept optimum; a long divergent path must not eventually claim convergence because of floating-point behavior. The legacy model-only wrappers inherit the same empty-path behavior.

The evaluator never hands the lasso one-class input: its two-class outer-fold, inner-fold, and tuning-holdout gates (`tuning_holdout` in `src/learn/mod.rs`) refuse such a partition first. Those gates are what keep one-class input from becoming an ordinary no-call: inside the evaluator's lasso code, an empty path from the one-class guard would look the same as a two-class null fit.

A deliberate empty path on a valid two-class fit remains different: it means no penalty admits a feature and is the existing null-fit/no-call case. A non-empty path with no converged candidate is a solver failure, not a null fit.

The finite-parameter guard closes the separate NaN issue: `f64::max` can ignore NaN, so a maximum over violations is not evidence of convergence by itself. The KKT residual is NaN whenever any of its terms is, which fails the tolerance, and NaN, positive infinity, and negative infinity in the intercept or any coefficient force `converged = false`. Weights that do not sum to a positive finite number stop the point at once, unconverged.

## Evaluator enforcement in revision 8

`eval::fit_lasso` now uses `fit_path_with_status` and `refit_path_with_status`.

1. Every unconverged tuning-path point is excluded before validation loss and penalty selection.
2. A deliberately empty two-class path remains a null fit. A non-empty path with no converged candidate returns an error.
3. The selected point is defensively checked again for convergence.
4. The refit at the selected penalty is checked before any score or model is returned. An empty or unconverged selected refit is an error.
5. Path-selection and refit errors name the dataset, split scheme, zero-based outer fold, and stage. Final-model errors say `final model on all samples (no outer fold)` rather than pretending to be an outer fold.
6. Solver failure is never converted into a no-call, and an unconverged selected/final model is never scored or saved.

Converged candidates retain the existing minimum-validation-log-loss selection and refit behavior. Non-lasso learners are unchanged.

Any recorded unconverged two-class lasso path point makes `evaluation.valid` false, even when a different converged candidate could be selected. Missing per-fold/final solver evidence, or the wrong number of per-fold records on a split, also makes a newly generated report invalid. Thus an excluded, non-selected failure remains visible and blocks a valid revision-8 report rather than being silently accepted.

An invalid report is written for review but is not a finished dataset. `grm learn` logs no `done` line for it: it stops the species with an error naming the dataset and the audit's reasons, as it does on a lasso error. A resumed run recomputes it rather than reusing it (see below). The decision is made where the report is written and where reuse is decided, not left to the checker.

## Report schema and independent checker

`REPORT_SCHEMA` remains **2**. `EVALUATION_REVISION` is now **10**. The optional Rust deserialization defaults preserve the ability to read older reports, but revision-8, revision-9 and revision-10 output must contain the new evidence.

Revision 10 keeps this lasso and its evidence unchanged. It adds the species identity check ([species-check.md](species-check.md)) and keeps k-mer-identical duplicates whole in one outer fold, listed in the new `kset_duplicates` field ([duplicate-isolation.md](duplicate-isolation.md)). A revision-10 report meets the revision-8 rules below, and a same-fingerprint revision-9 report is stale against revision 10.

Revision 9 changes only the lasso's active-set solve (above). It adds no report field and changes no field's meaning: `sweep_cap_hits` still counts IRLS steps whose solve stopped short, and the new step to a first crossing is part of an exact solve, not a short one. A revision-9 report meets the revision-8 rules below. The bump keeps revision-8 reports, whose lasso fits the old solve could leave unconverged where the new one converges, from being reused: a same-fingerprint revision-8 report is stale against revision 9.

The revision exists to keep reports from different evaluation code apart: a saved report is stale when its revision differs, because `input_fingerprint` hashes inputs and settings, not code. `EVALUATION_REVISION` stayed 8 when the active-set solver replaced coordinate descent, because no revision-8 report existed then: the only revision-8 run until that point, the large-species scan, wrote scan records, not reports, and stopped on its first fit. The 275-report revision-8 run that followed did write reports, so the revision-9 fix to the active-set solve bumps the revision, and a same-fingerprint revision-8 report is stale against revision 9.

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

- accepts revisions 3 through 10 and rejects future revisions, applying the revision-8 lasso rules to revisions 9 and 10;
- refuses a root whose reports have more than one evaluation revision, such as revision-8 reports left beside revision-9 ones by a run stopped part-way, so a revision-10 root must hold only revision-10 reports;
- keeps the existing revision 3–7 schema, no-call, metric, split, and log semantics;
- requires per-fold and final lasso status on every available revision-8 split, including reports whose lasso rows have null fits;
- rejects missing fields, wrong types, negative or out-of-range counts and indices, a wrong number of fold records, a non-empty path with no selected point (or the reverse), inconsistent selected/refit presence, and a final model inconsistent with an empty-path record;
- checks each fold's record against that fold's scores: a lasso fold is called exactly when its `fold_complexity` (the refit's weight count) is positive, and a fold whose record has no selected point has no calls; and
- requires zero unconverged two-class lasso path points and true selected/final convergence before accepting `valid: true`.

The checker tests cover a valid revision-8 report, a valid revision-9 report and the revision-8 rules applied to it, a valid revision-10 report with an isolated duplicate pair, a root mixing two revisions (9 with 10, 8 with 9, 7 with 9, 7 with 8) beside one holding two revision-9 reports, a valid empty-path final fit with a null model, a valid empty-path record on an uncalled fold, and future-revision rejection (revision 11). `metric-audit/audit.py` has no revision list: `--revision 9` audits a revision-9 run and `--revision 10` a revision-10 one, and tests audit both clean.

Negative cases cover missing, malformed and unconverged records; a saved model whose record has no refit; a dropped or extra fold record; a 30-point path with no selected point; an empty-path record on a called fold; a weight count that disagrees with the fold's calls; and a short `fold_complexity`.

Each check's removal fails at least one test (mutations P9, P11 and N4a–c in the close-out lane). Revision-7 reports do not need to be retrofitted with these fields, and the record-against-score check applies only from revision 8.

### Stale and reuse behavior

`input_fingerprint` deliberately hashes inputs and settings, not code. The revision bump is therefore part of the same change as evaluator enforcement. A saved revision-7 report with an otherwise identical fingerprint is stale with reason `evaluation audit revision changed`; it cannot be reused inside a revision-8 report.

A current-revision report whose `evaluation.valid` is not `true` is stale with reason `saved report's evaluation audit is not valid`, so an invalid report is recomputed on resume, never passed off as done. Recomputing it on the same code and inputs gives the same invalid report and stops the species again. Tests cover both rules in `stale_reason` and in the `analyze_species` reuse path, and `done_line` refuses an invalid audit. No forced partial real-data rebuild was run.

## Test coverage and required revision-8 gate

The saved cross-check panels are in `tests/data/lasso-stopping-fixtures.txt`. They preserve the 16 fixtures and candidate paths from the earlier cross-check, the best known objective at each point, and the 300-sweep trace. The tests no longer rerun the 300-sweep solver, which is gone; its saved trace and models stay as the record of what revision 7 did. Bit-exact checks (a refit or wrapper repeating a path) remain pinned to the Linux platform because `powf` and `exp` can differ by one ULP on other platforms.

The lasso suite has 28 default tests and one ignored release test. All 91 saved points now run by default: the C. jejuni panels that once took 2.3 million sweeps take seconds unoptimized. In particular:

- Revision 9's two real fitting sets are committed fixtures: `tests/data/lasso-vibrio-fold0-fixture.txt` (57 samples, 4,000 screened columns) and `tests/data/lasso-salmonella-fold3-fixture.txt.gz` (867 samples, 4,000 columns, gzipped). `vibrio_cholerae_fold0_converges_at_every_point` and `salmonella_enterica_fold3_converges_at_every_point` require every point of each 30-point path converged with no short solve and an independently recomputed KKT residual at most 1e-6; both fail on revision 8 at the point that stopped it. `the_saved_real_paths_are_the_same_at_any_thread_count` runs both paths in rayon pools of 1 and 64 threads, the second twice, and requires the same models and statuses, bit for bit.
- `a_coefficient_within_rounding_of_zero_does_not_stall_the_solve` builds the stalled state directly, a least-squares solve with column 0 at plus or minus 1e-200 and the rest of the active set past its optimum, with four kinds of partner for column 0: an exact complement, an exact copy (collapse turned off), two disjoint columns whose union it is, and near-copies differing from it in one and two samples. It checks that the state is the stalled one (no candidate step shows a fall, column 0 crosses first at t below 1e-100, the complement violates its condition and is dependent), then requires the solve to finish exactly at the brute-force optimum. Revision 8 fails all four.
- `a_hidden_fall_steps_to_the_first_crossing_only_while_every_sign_holds` tests the step decision directly: the line search's point when it lowers the objective; otherwise the first crossing, also while a column is entering in the direction of its sign; and no step when an entering column moves against its sign or not at all, or nothing crosses before t = 1. `the_first_crossing_is_the_earliest_before_the_full_step` pins the crossing itself, the intercept and a crossing at exactly t = 1 left out. No real or synthetic fit has reached a hidden fall while a column was entering, so this branch is tested only here.
- `an_exchange_must_lower_the_objective_by_more_than_rounding` sets the standard deviations so that an exchange's rate is -1e-14 (refused, nothing changed) or -5e-4 (taken, with both coefficients that reach zero together leaving, also when they start one ulp apart), and `violators_tie_to_the_lowest_column` orders 64 tied violations.

- `saved_panels_converge_and_agree_with_the_old_solver_where_it_settled` requires every two-class saved point to converge with an independently recomputed KKT residual at most 1e-6 and an objective within 1e-6 of the best known, no solve to stop short, and the raw one-class points to stay unconverged. At the 38 points where the 300-sweep solver settled without a cap, the two solutions must agree within 1e-9 in objective and 1e-4 in any fitted probability.
- `campylobacter_failures_reach_the_reference_optimum` pins the three former C. jejuni failures: the saved 300-sweep models have KKT and excess objective above 1e-4, and the current solver converges there in 6, 6 and 5 IRLS steps with no short solve.
- `public_wrapper_certifies_a_saved_panel_the_old_cap_left_unconverged` goes through `fit_path_with_status` on `streptococcus_agalactiae__clindamycin_lineage_fold0`, whose saved point 6 hit the sweep cap in all 16 IRLS steps, and requires every point certified.
- `public_refit_and_model_only_wrappers_repeat_the_converged_path` refits that panel's penalties through `refit_path_with_status` and requires the same models and statuses back, bit for bit, as well as from the model-only `fit_path` and `refit_path`.
- `a_guard_that_ends_the_solve_is_flagged` shows that both guards, the active-set iteration guard and the IRLS guard, leave a point unconverged and are counted, and checks wrapper equivalence on a synthetic panel.
- `near_twin_columns_are_certified_in_few_irls_steps` fits twelve columns with three one-sample twins each, at most 10 IRLS steps per point.
- `exact_copies_share_one_coefficient_and_leave_the_objective_unchanged` compares a panel with copies to one without them, checks the collapse map, and checks that a path capped at any `max_nonzero` stops right after its first point with more distinct non-zero columns.
- `a_column_dependent_on_the_active_set_is_exchanged_in` uses a column that is the union of two disjoint active ones, so the Gram matrix with all three is singular. `only_an_exactly_dependent_column_is_refused_by_the_factor` checks that the pivot test refuses that column and accepts one that differs from it in a single floor-weight sample.
- `a_coefficient_warm_started_with_the_wrong_sign_crosses_zero` starts a coefficient at the wrong sign. It and the dependent-column test compare the active-set solve with a brute-force minimum over every sign pattern.
- `the_line_search_stops_where_a_coefficient_reaches_zero` and `an_irls_step_that_raises_the_objective_is_halved` test the line search and the step halving directly.
- `the_kkt_residual_covers_the_intercept_and_every_condition` checks each term of the convergence test on hand-computed values, and `every_point_is_reported_as_it_is_solved` checks the observer, including that the solver's own KKT residual matches the independent one to 1e-12.
- `one_class_long_paths_are_empty_by_design` exercises both public wrappers with 4,001 penalties for all-resistant and all-susceptible fitting masks.
- `non_finite_parameters_never_report_convergence` covers NaN and both infinities, a NaN in the KKT residual, and NaN or infinite working weights in the solver.
- The ignored `mtb_delamanid_fold0_reaches_the_certified_optima` runs the fit that stopped the scan, from the diagnosis folder (gitignored, so it runs only where that folder exists). It requires all 30 tuning-path points converged with no short solve, each within 1e-6 above the diagnosis's certified optimum (KKT 1e-9) and no more than 1e-9 below it, selection of point 3, a converged refit with 53 weights, and no call changed against revision 7 on the fold's 251 test samples.
- Evaluator tests prove that a lower-loss unconverged candidate cannot be selected, a path with no converged candidate cannot be refit, and an unconverged refit returns an error with dataset/scheme/fold/stage context. An unconverged or empty final-model refit names the final model instead of a fold.
- Analysis tests prove that each of these makes `valid` false: an unconverged path point, selected point or refit; a missing final record; a split with no per-fold records or no lasso row; and a split with too few or too many records. An unavailable split needs none.
- Further analysis tests prove that revision-7 results cannot be reused at revision 8, and that an invalid current report is neither logged as done nor reused.
- `lasso_only_repeats_the_evaluator_lasso_fits` and `the_lasso_scan_fits_the_report_inputs_as_the_evaluator_does` pin the scan harness to the evaluator (see the scan section), and check that every solved point is reported once, in order, with its fit's status.

The required release command, including the ignored MTB delamanid test, is:

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
python3 -I -B scripts/lasso-scan/test_scan.py
scripts/lasso-scan/smoke.sh
```

The metric-audit suites call `scripts/metric-audit/foldgen/target/release/foldgen`, so build it first. `smoke.sh` builds the scan binary, writes a 150-genome synthetic data root, and drives the real binary under the registered rule. No real or candidate data root is used by these tests.

The shared checkout was deleted at 04:41 UTC on 2026-10-28 while this work was in progress. Its edits were recovered into commits `5302848` and `fee9f9f` and finished in `fb7ba09`. The pre-change files and hashes are under `data/research/grm-lasso-revision8-migration-20261028`. The finishing lane's gate logs, mutation checks and validation summary are under `data/research/grm-lasso-revision8-finish-20261028`, and the independent review's are under `data/research/grm-lasso-r8-review-20261028`.

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

These counts come from the saved trace. The 26 flagged points needed the 300-sweep solver's status rule to count, and the revision-8 tests that verified them ran that solver; the current tests check the 46 and the 15 from the saved trace alone.

Under the coordinate-descent solver with 100,000 sweeps, all 84 two-class saved points converged, but four C. jejuni points hit the sweep cap in an intermediate IRLS step, and the three prior C. jejuni failures took 7, 9 and 18 IRLS steps. Under the active-set solver all 84 converge with independently recomputed KKT residual at most 1e-6 and objective within 1e-6 of the best known, no solve stops short, and the three prior failures take 6, 6 and 5 IRLS steps. The seven raw one-class fixture points remain unconverged in the internal trace used for accounting; public one-class calls return empty paths instead of running that trace. Refitting at the same penalties on the same mask returns the identical model and status.

## What revision 8 can change

The small deployment exercise (four datasets, 99–508 samples) reproduced every old-solver score/model on replay and found zero unconverged points among 2,107 after the fix. Thirty-four of 44 models changed numerically without a selected-penalty or call change, and lasso time increased from 15.0 s to 28.8 s (1.92×). This is evidence for the revision bump, not a prediction of full-run results or runtime. Its evidence directory, `grm-lasso-stopping-fix-20261027-r8kq3v`, is not available locally, so these figures cannot be rechecked here; the blocking review checked them on 2026-10-28. These figures are from coordinate descent. The active-set solver reaches the same optima within the KKT tolerance, so its models differ from those only within that tolerance, but it was not rerun on those four datasets.

A full revision-8 field-level diff must treat the following as liable to change:

- lasso scores, metrics, ROC, fold balanced accuracy, fold complexity, no-call counts, and timing;
- `models.lasso` and the new solver-status fields;
- `features` and annotations, because the list includes the top 30 lasso weights (newly included loci need annotation downloads);
- best-learner labels and corresponding log wording;
- `evaluation.revision`, report validity evidence, timing fields, and the report index.

SCM, CART, forest, boosting, MIC, catalog, and probe outputs should remain identical. The earlier cross-check reproduced all 278 non-lasso rows; the full-run diff must enforce that rather than assume it. Boosting is unaffected by this solver fix. XGBoost and the default `colsample = 0.5` remain a separate verification gap and do not block this revision-8 prerequisite work.

## Large-species lasso scan

The scan measures large-species lasso cost under the rule registered here before any 275-report run. Its first run, with coordinate descent, stopped on its first fit (see "Why coordinate descent was replaced"). The active-set solver's only large measurement so far is that same fit, MTB delamanid random fold 0 (see "Costs"); the other ten fits of that dataset and the other 17 datasets have not been run with it. The rerun uses the same list, rule and budgets. **Running it needs separate user approval.**

### The harness

`scripts/lasso-scan/` runs the shipped revision-8 code, not a copy. The earlier `lasso_deploy.rs` harness, which replayed a copied revision-7 `fit_lasso`, is not used.

- `main.rs`, built as the Cargo example `lasso-scan`, calls `grm::analysis::lasso_scan`. That function loads a species with `load_inputs` and builds each dataset's cohort, projection and lineages with `eval_inputs`, the same function `analyze_dataset` uses. It reports the dataset's input fingerprint, then runs `eval::lasso_only`.
- `eval::lasso_only` runs only the lasso fits of `evaluate`, in its order: each available split's outer folds, then the final model. It records each penalty point's `Convergence`, for the tuning path and the refit, and reports each point as soon as it is solved. The binary writes that report to stderr as one line per point, such as `mycobacterium_tuberculosis__delamanid random fold 0 path point 3: 3 IRLS steps, 74 active, KKT 1.1e-9, 0.33 s`, so a fit that is stopped shows how far it got.
- `scan.py` is the driver. It checks the list against each dataset's revision-7 report, runs the binary on one dataset at a time, and enforces the rule below. It writes only its output directory: `plan.json`, `events.jsonl`, one stderr log per dataset (holding the progress lines), `summary.json`, and `stop.json`. It resolves `--binary` to an absolute path first: `Path("./lasso-scan")` prints as `lasso-scan`, which `subprocess` would look up on `PATH`.
- `scan-list.json` is the registered scan list, and `make_scan_list.py` regenerates it from the 275 revision-7 candidate reports.
- `test_scan.py` drives every stop condition with a scripted fake binary, and `smoke.sh` runs the real binary on a synthetic data root.

`lasso_only` shares `evaluate`'s input checks, outer partitions, split plan, per-fold candidates, fold seeds, final-model plan and `fit_lasso` code, and times each fit the way `evaluate` times a learner. Two default tests require its records and models to equal `evaluate`'s, one on synthetic data and one through a synthetic data root. Mutations that change its seed, its candidates or its lineage groups fail them.

The binary is a Cargo example, so `cargo build --release` and the deployment build leave it out.

### Scan set

The set is chosen by measured revision-7 lasso cost, not by species name. The five species with the most revision-7 lasso seconds are M. tuberculosis (3,814 s), K. pneumoniae (2,389 s), S. enterica (1,701 s), E. coli (1,559 s) and S. aureus (1,273 s). Each contributes its three costliest datasets, plus its datasets with the most samples and the most features.

That gives 18 datasets. They cover cost ranks 1–4, 9, 11, 13 and 14 of all 275, and the largest matrices: E. coli ampicillin, with 3.1 million features, and S. enterica ampicillin, with 1,997 samples. Together they hold 4,005 s of revision-7 lasso work (1.11 h) once the final fits are counted. The scan runs them most expensive first:

| Rank | Dataset | Rev-7 lasso s | Fit budget s | Dataset budget s |
|---:|---|---:|---:|---:|
| 1 | `mycobacterium_tuberculosis__delamanid` | 420.9 | 842 | 2,315 |
| 2 | `mycobacterium_tuberculosis__linezolid` | 328.6 | 657 | 1,807 |
| 3 | `salmonella_enterica__ampicillin` | 315.3 | 631 | 1,734 |
| 4 | `mycobacterium_tuberculosis__clofazimine` | 310.9 | 622 | 1,710 |
| 9 | `mycobacterium_tuberculosis__rifampin` | 229.1 | 458 | 1,260 |
| 11 | `staphylococcus_aureus__erythromycin` | 210.2 | 420 | 1,156 |
| 13 | `salmonella_enterica__tetracycline` | 203.0 | 406 | 1,117 |
| 14 | `klebsiella_pneumoniae__ciprofloxacin` | 200.5 | 401 | 1,103 |
| 16 | `staphylococcus_aureus__clindamycin` | 172.9 | 346 | 951 |
| 17 | `klebsiella_pneumoniae__piperacillin_tazobactam` | 167.3 | 335 | 920 |
| 20 | `klebsiella_pneumoniae__sulfamethoxazole_trimethoprim` | 161.3 | 323 | 887 |
| 21 | `salmonella_enterica__chloramphenicol` | 159.8 | 320 | 879 |
| 25 | `escherichia_coli__ceftazidime` | 147.2 | 294 | 810 |
| 26 | `escherichia_coli__ampicillin` | 147.1 | 294 | 809 |
| 29 | `staphylococcus_aureus__penicillin` | 139.9 | 280 | 770 |
| 37 | `escherichia_coli__amoxicillin_clavulanate` | 122.7 | 245 | 675 |
| 38 | `klebsiella_pneumoniae__gentamicin` | 121.3 | 243 | 667 |
| 65 | `klebsiella_pneumoniae__amikacin` | 83.1 | 166 | 457 |

`scan-list.json` also records each dataset's samples, features, outer fits, report seconds, revision-7 input fingerprint and why it was chosen.

### Registered scan stop rule

`r7_fit` is a dataset's revision-7 lasso seconds divided by its outer fits, which number 10 for every listed dataset. `h` scales every budget: it is the scan host's time for a fixed piece of work divided by the revision-7 host's. It is 1 on the Vast host class that ran revision 7, where the scan is to run.

| Limit | Rule | Registered value at h = 1 |
|---|---|---|
| Per fit (tuning path and refit) | max(60 s, 20 × `r7_fit` × h) | 166–842 s |
| Per dataset (fit time of the 10 outer fits and the final fit) | 5 × 11 × `r7_fit` × h | 457–2,315 s |
| Whole scan (wall-clock from launch, loading included) | the sum of the dataset budgets | 20,026 s (5.56 h) |
| Disk | free space on the output filesystem, minus what the scan still has to write, stays at or above 30 GiB | estimate below |

The 1.92× in "What revision 8 can change" was coordinate descent on four small datasets and does not predict the rerun. On MTB delamanid fold 0, the active-set path took less time on the laptop than revision 7's path did there (see "Costs").

The scan stops at once on any of the following. It kills the binary, keeps every log and partial output, and runs no further dataset:

- a lasso error: no converged candidate, an unconverged selected point, or an empty or unconverged refit;
- any unconverged path point, selected point or refit. That is also exactly what makes a report's `evaluation.valid` false;
- any budget overrun, whether per fit, per dataset or for the whole scan;
- free disk that would cross the floor before the next dataset, which counts as a budget overrun;
- inputs that are not the revision-7 report's: a different input fingerprint, or a report that does not match the list;
- a binary failure, unexpected output, or a signal that stops the driver.

It never retries, with raised limits or otherwise. A stop blocks the 275-report run pending review.

Without stopping, it records each penalty point's `irls_steps` and `sweep_cap_hits`, each fit's seconds, and each dataset's ratio. The binary runs with `RAYON_NUM_THREADS=64`, as revision 7 did. If the driver itself is killed, the kernel kills the binary too, so no fit outlives its budget.

**Go/no-go for the 275-report run.** For each dataset, r is its scan seconds over the 10 outer fits divided by (its revision-7 lasso seconds × h). The final fit is timed and budgeted but left out of r, because revision 7 did not time it. Go only if all of these hold:

- the worst r is at most 5;
- the median r is at most 3;
- the projected 275-report run is at most 16 h, which fits one overnight Vast session;
- at most 5% of the solved penalty points, tuning and refit together, hit the sweep cap in some IRLS step. Above that, a review comes first. With the active-set solver this counts points with an IRLS step whose solve stopped short (`sweep_cap_hits` above 0, the same field `scan.py` reads); the threshold is unchanged. The saved fixtures' rate is 0 of 84 with the active-set solver and was 4 of 84 with coordinate descent.

The projection is (25,088 s + (median r − 1) × 1.1 × 14,989 s) ÷ 3,600, where 25,088 s and 14,989 s are the 275 revision-7 reports' total report and lasso seconds. It gives 11.2 h at 1.92×, 16.1 h at 3× and 25.3 h at 5×.

`scan.py` exits 0 for go, 1 for no-go, 2 when the rule stops the scan (`stop.json` says why) and 3 when it refuses to start.

### Disk

The scan copies no input and writes no temporary file. It reads about 4.4 GiB in place under the input root: the five species' matrices, counts and structures, the curated phenotypes and the BV-BRC genome list.

Each fit line it writes holds at most 30 tuning and 30 refit statuses and one weight per screened feature, so it stays under 256 KiB. The stderr log holds one progress line of under 140 bytes per solved point, under 100 KiB for a dataset's 11 fits. Add a 1 MiB stderr allowance per dataset and 4 MiB of fixed files, and the 18 datasets need at most **75,249,664 bytes (71.8 MiB)**. For scale, the synthetic smoke run wrote 106,674 bytes for one dataset, about half of it the 607 progress lines. Staging the release binary (7,725,344 bytes as built for `znver2` from this lane's source, needing glibc 2.34 or later), `scan.py` and the list adds 7,759,021 bytes (7.4 MiB). The scan therefore needs at most **83,008,685 bytes (79.2 MiB)** of new disk in all; the inputs-only check writes a few kilobytes more. Building on Vast instead would add a release target directory of about 443 MB, which is why the binary is built elsewhere.

`scan.py` prints its estimate before launch. It refuses to start if free space minus the estimate would fall below 30 GiB, and it checks again before each dataset. On 2026-10-28 Vast had 32.35 GiB free against that floor, and other lanes there were capped at 1.6 GiB of new disk combined; the scan's need is well inside that cap.

### Inputs on Vast

The revision-7 candidate root, `/root/grm-extra-species-analysis-20261027/data`, holds no matrices for the five scan species. They are in `/root/grm-revision7-preflight/data`. The candidate copied that root's counts, structures and reports for the 13 original species hash-equal, and its reports carry the same input fingerprints (M. tuberculosis delamanid: `275e461dec002050` in both). The scan reads that root.

Its `catalog` is a symlink to `/root/grm/data/catalog`. If that catalog has changed since revision 7, loading fails or the fingerprints differ. The inputs-only check below catches either before any fit. It also catches a matrix the current code cannot load, since the matrix commits `7cc39ed` and `0df8eb9` came after revision 7.

### Launch procedure

The launch needs separate user approval. Build the binary from a clean export of the frozen commit, on Ubuntu 24.04 with glibc 2.39 as on Vast, so nothing is compiled on Vast.

`.cargo/config.toml` builds for the build machine's CPU (`target-cpu=native`), and revision 7 was built on Vast for its AMD EPYC 7K62, which LLVM calls `znver2`. So set `RUSTFLAGS` to build for that CPU instead. Rust does not contract or reorder floating-point operations for a target CPU, so the fits are the same; the timing is comparable with revision 7's.

Then stage the binary in a new directory on Vast, check the plan and inputs without fitting, and launch detached:

```bash
# Build machine, clean export of the frozen commit
RUSTFLAGS="-C target-cpu=znver2" cargo build --locked --offline --release --example lasso-scan
# -> target/release/examples/lasso-scan

# Vast: stage lasso-scan, scripts/lasso-scan/scan.py and scripts/lasso-scan/scan-list.json in $SCAN
SCAN=/root/grm-lasso-scan-r8-<date>; ROOT=/root/grm-revision7-preflight/data; COMMIT=<frozen sha>
cd "$SCAN"
python3 -I -B scan.py --plan --input-root "$ROOT" --out "$SCAN/inputs" --binary ./lasso-scan --commit "$COMMIT"
python3 -I -B scan.py --inputs-only --input-root "$ROOT" --out "$SCAN/inputs" --binary ./lasso-scan --commit "$COMMIT"
nohup nice -n 10 python3 -I -B scan.py --input-root "$ROOT" --out "$SCAN/scan" --binary ./lasso-scan \
    --commit "$COMMIT" > "$SCAN/driver.log" 2>&1 &
```

### Still gated

1. Before the scan: an independent review of the active-set solver (its math, termination and duplicate handling, `EVALUATION_REVISION` staying 8, and this document). Then rerun the release lasso tests, including the MTB delamanid fold-0 test, and the r8kq3v `xcheck.sh` on the frozen commit, with a separate Python interpreter rather than the hard-stopped GWAS lane's. Then run and review the scan under this rule.
2. Decide how to respond if a lasso error or an invalid report stops `grm learn` for a species, since no downstream consumer handles fewer than 275 reports.
3. Only after the scan passes, run the full 275-report evaluation with revision-8 source into a **new data root**, leaving revision-7 candidate, backup, and published roots untouched.
4. Require every report to pass the revision-8 checker with zero unconverged two-class lasso fits, and run `metric-audit --revision 8`. Then do the field-level migration diff and the completion-log audit.
5. Rerun the downstream consumers listed below.
6. Rebuild and review the report/web bundle. Publication is a further explicit decision.

`models.lasso_fit` sits between `models.lasso` and `rf_importance`, so whole revision-8 files will not be byte-identical to revision-7 ones. The migration diff has to compare non-lasso fields field by field.

## Deployment settings and downstream work

The current published release remains revision 7. `scripts/deployment/` retains its revision-7 defaults and 237-dataset published-release expectations; this lane deliberately did not edit it while its independent review is active.

At a future, authorized revision-10 staging/publishing operation, the operator must explicitly set:

```bash
GRM_EXPECT_REVISION=10
GRM_EXPECT_DATASETS=<reports in the revision-10 run>
```

The revision-9 run wrote 275 reports. The species check reselects some cohorts, so take the revision-10 count from its run, not from this page.

Those are operation-time requirements, not changes to the current defaults and not evidence that staging or publication has occurred.

### Downstream consumers

The consumer map in [`data/research/grm-lasso-consumers-20261028/`](../data/research/grm-lasso-consumers-20261028/README.md) traced every reader of lasso output from code and saved artifacts, locally and on Vast. Follow its minimum rerun set once the 275-report revision-8 candidate exists.

Two lanes refit lasso and must be ported before they rerun: the temporal/study held-out lane (`grm-research-temporal-study-splits-rev7cand-r1`) and the duplicate-isolated lane (`grm-research-dup-leakage-all23-r1`). Each runs a probe that carries its own copy of the revision-7 `fit_lasso`, and that copy still selects among all path points. Rebuilding a probe on revision-8 source is therefore not enough: its unconverged points could still win, silently.

So first port each probe's `fit_lasso` to the revision-8 rules: use `fit_path_with_status`, select only converged points, refuse an unconverged selected or refit point, and record the status. Then re-point the probe at revision-8 source and rerun every job. The other consumers:

- **Saved-output recomputes:**
  - Phase A intervals (`grm-all23-saved-intervals-20261027-kimi`);
  - the all-23 selective saved analysis (`grm-all23-selective-20261027-kimi`);
  - the de-escalation all-learner no-call audit table, but not the GBT replay (`deescalation-all23-saved-20261027-vozjjhqn`);
  - the MTB revision-7-candidate scale-up's `audit.py` and `tables.py` section 4 (`grm-research-mtb-scaleup-rev7cand-r1`).
- **Second-order:**
  - the clindamycin label recheck (`grm-clindamycin-label-recheck-20261027-tuwf2end`), after the candidate and the temporal rerun;
  - the duplicate-isolated lane's lasso-excluded sensitivity script, after that lane's rerun;
  - one new uncertainty-validator run, after the saved recomputes;
  - a repeat of the duplicate-isolated independent review.
- **No rerun:** the all-MIC cohort lanes and their reviews, co-resistance/time (`coresistance-time-all23-20261027T101316Z`), and the genome-source and duplicate audit (`grm-research-provenance-dup-all23-r1`).
- **No rerun while non-lasso outputs are unchanged:** the nested RF/GBT refits, the GBT de-escalation replay and the all-23 MIC baseline (`grm-all23-mic-baseline-20261027-kimi`). If the field-level diff shows any non-lasso bytes or MIC predictions changed, these three join the rerun set.
- **Historical:** the revision-5 lanes, including `grm-research-selective-calling-r1`, `grm-research-deescalation-r1` and `grm-research-clindamycin-imlsb-r1`.

An earlier version of this list predated the map. It named the all-MIC, MIC-baseline, co-resistance/time and clindamycin iMLSb lanes as consumers, and asked only for "revision-8-source reruns" of the two refit lanes. The map supersedes it.

Older revision-5 lanes remain historical. AMRFinder/GWAS work reads the revision-5 snapshot and is not on this migration list; this does not authorize restarting hard-stopped work.

**Completion of these prerequisites does not authorize the large-species scan, the full revision-8 run, downstream reruns, staging, or publication.**
