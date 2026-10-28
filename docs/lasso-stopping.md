# Lasso stopping fix and revision-8 migration prerequisites

Status: the solver change is approved, the code/test/checker prerequisites are implemented and independently reviewed, and the review's findings are closed. The large-species scan harness is prepared and its stop rule is registered below. **The large-species scan and the full 275-report revision-8 evaluation have not run.** No scan, biological evaluation, report rebuild, deployment, or publication is authorized by this document.

The blocking review that set these prerequisites wrote its verdict to `data/research/grm-lasso-stopping-fix-20261027-r8kq3v/MIGRATION-REVIEW.md`. That file was lost in the 2026-10-28 deletion and is not on Vast either, so its findings survive only in the review's thread and in this summary: it approved the solver change and required evaluator enforcement, a revision bump, schema/checker migration, edge guards, and stronger tests before any revision-8 run. The guards were then reviewed independently; that review, its mutation logs and the revision-7 timing tables it measured are in [`data/research/grm-lasso-r8-review-20261028/`](../data/research/grm-lasso-r8-review-20261028/README.md).

## Approved solver change

`src/learn/lasso.rs` formerly stopped coordinate descent after 300 sweeps without preserving a usable status. On ill-conditioned panels this stopped short of the optimum. The approved change retains 25 IRLS steps and raises the coordinate-descent limit to 100,000 sweeps per IRLS step (glmnet's default `maxit`). The objective, unpenalized intercept, standardization, class weights, penalty grid, tolerances (1e-6 sweep and 1e-5 IRLS), 1e-5 weight floor, warm starts, `max_nonzero`, and refit semantics are unchanged.

`Convergence { irls_steps, sweep_cap_hits, converged }` records how each penalty point stopped. A point converges only when:

- the IRLS change is below 1e-5;
- that step's coordinate descent settled rather than reaching the sweep cap; and
- the intercept and every coefficient are finite.

`fit_path_with_status` and `refit_path_with_status` return `LassoFit { model, convergence }`. The compatibility wrappers `fit_path` and `refit_path` retain their model-only signatures and report unconverged points on stderr; the evaluator no longer relies on those warnings.

### One-class and non-finite edges

Both public status-returning entry points return an **empty path by design** when the fitting mask contains only one class. One-class logistic regression has no finite intercept optimum; a long divergent path must not eventually claim convergence because of floating-point behavior. The legacy model-only wrappers inherit the same empty-path behavior.

The evaluator never hands the lasso one-class input: its two-class outer-fold, inner-fold, and tuning-holdout gates (`tuning_holdout` in `src/learn/mod.rs`) refuse such a partition first. Those gates are what keep one-class input from becoming an ordinary no-call: inside the evaluator's lasso code, an empty path from the one-class guard would look the same as a two-class null fit.

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

Any recorded unconverged two-class lasso path point makes `evaluation.valid` false, even when a different converged candidate could be selected. Missing per-fold/final solver evidence, or the wrong number of per-fold records on a split, also makes a newly generated report invalid. Thus an excluded, non-selected failure remains visible and blocks a valid revision-8 report rather than being silently accepted.

An invalid report is written for review but is not a finished dataset. `grm learn` logs no `done` line for it: it stops the species with an error naming the dataset and the audit's reasons, as it does on a lasso error. A resumed run recomputes it rather than reusing it (see below). The decision is made where the report is written and where reuse is decided, not left to the checker.

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
- rejects missing fields, wrong types, negative or out-of-range counts and indices, a wrong number of fold records, a non-empty path with no selected point (or the reverse), inconsistent selected/refit presence, and a final model inconsistent with an empty-path record;
- checks each fold's record against that fold's scores: a lasso fold is called exactly when its `fold_complexity` (the refit's weight count) is positive, and a fold whose record has no selected point has no calls; and
- requires zero unconverged two-class lasso path points and true selected/final convergence before accepting `valid: true`.

The checker tests cover a valid revision-8 report, a valid empty-path final fit with a null model, a valid empty-path record on an uncalled fold, and future-revision rejection.

Negative cases cover missing, malformed and unconverged records; a saved model whose record has no refit; a dropped or extra fold record; a 30-point path with no selected point; an empty-path record on a called fold; a weight count that disagrees with the fold's calls; and a short `fold_complexity`.

Each check's removal fails at least one test (mutations P9, P11 and N4a–c in the close-out lane). Revision-7 reports do not need to be retrofitted with these fields, and the record-against-score check applies only from revision 8.

### Stale and reuse behavior

`input_fingerprint` deliberately hashes inputs and settings, not code. The revision bump is therefore part of the same change as evaluator enforcement. A saved revision-7 report with an otherwise identical fingerprint is stale with reason `evaluation audit revision changed`; it cannot be reused inside a revision-8 report.

A current-revision report whose `evaluation.valid` is not `true` is stale with reason `saved report's evaluation audit is not valid`, so an invalid report is recomputed on resume, never passed off as done. Recomputing it on the same code and inputs gives the same invalid report and stops the species again. Tests cover both rules in `stale_reason` and in the `analyze_species` reuse path, and `done_line` refuses an invalid audit. No forced partial real-data rebuild was run.

## Test coverage and required revision-8 gate

The saved cross-check panels are in `tests/data/lasso-stopping-fixtures.txt`. They preserve the 16 fixtures and candidate paths from the earlier cross-check, the best known objective at each point, and the 300-sweep trace. Bit-exact fixture checks remain pinned to the Linux platform because `powf` and `exp` can differ by one ULP on other platforms.

The lasso suite now has eleven tests: the original seven plus public-wrapper, refit-and-model-only-wrapper, long-one-class, and finite-parameter regressions. In particular:

- `public_wrapper_uses_the_raised_cap_on_a_saved_panel` is a default, non-ignored test through `fit_path_with_status` on `streptococcus_agalactiae__clindamycin_lineage_fold0`. Its old 300-sweep path has unconverged points; reverting the inner cap to 300 fails this default test.
- `public_refit_and_model_only_wrappers_use_the_raised_cap_on_a_saved_panel` is its default counterpart for the refit, which the evaluator scores and saves. Reverting the refit cap to 300, or pointing `fit_path` at the old limits, fails it (mutations L7 and L8).
- `a_cap_that_ends_the_solve_is_flagged` checks wrapper equivalence on a synthetic panel: path length, penalties/models with `same()`, and statuses for both `fit_path` and `refit_path`.
- `one_class_long_paths_are_empty_by_design` exercises both public wrappers with 4,001 penalties for all-resistant and all-susceptible fitting masks.
- `non_finite_parameters_never_report_convergence` covers NaN and both infinities and exercises the solver with a NaN working weight.
- Evaluator tests prove that a lower-loss unconverged candidate cannot be selected, a path with no converged candidate cannot be refit, and an unconverged refit returns an error with dataset/scheme/fold/stage context. An unconverged or empty final-model refit names the final model instead of a fold.
- Analysis tests prove that each of these makes `valid` false: an unconverged path point, selected point or refit; a missing final record; a split with no per-fold records or no lasso row; and a split with too few or too many records. An unavailable split needs none.
- Further analysis tests prove that revision-7 results cannot be reused at revision 8, and that an invalid current report is neither logged as done nor reused.
- `lasso_only_repeats_the_evaluator_lasso_fits` and `the_lasso_scan_fits_the_report_inputs_as_the_evaluator_does` pin the scan harness to the evaluator (see the scan section).

The refit test refits the saved panel's converged penalties through `refit_path_with_status` and requires the converged path back, point for point. It also requires the model-only `fit_path` and `refit_path`, which nothing in production calls, to return the same models. Before it, only the ignored C. jejuni release test caught a refit-cap revert.

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

Under the new limits, all 84 two-class saved points converge, with independently recomputed KKT residual at most 1e-6 and objective within 1e-6 of the best known. Four C. jejuni points still hit the 100,000-sweep cap in an intermediate IRLS step and settle later. The seven raw one-class fixture points remain unconverged in the internal trace used for accounting; public one-class calls now return empty paths instead of running that trace.

At the three prior C. jejuni failures, the old limits stop at 25 IRLS steps with 25 sweep-cap hits, KKT above 1e-4, and excess objective above 1e-4. The new limits converge in 7, 9, and 18 IRLS steps, which the ignored release test now pins. Refitting at the same penalties on the same mask returns the identical model and status.

## What revision 8 can change

The small deployment exercise (four datasets, 99–508 samples) reproduced every old-solver score/model on replay and found zero unconverged points among 2,107 after the fix. Thirty-four of 44 models changed numerically without a selected-penalty or call change, and lasso time increased from 15.0 s to 28.8 s (1.92×). This is evidence for the revision bump, not a prediction of full-run results or runtime. Its evidence directory, `grm-lasso-stopping-fix-20261027-r8kq3v`, is not available locally, so these figures cannot be rechecked here; the blocking review checked them on 2026-10-28.

A full revision-8 field-level diff must treat the following as liable to change:

- lasso scores, metrics, ROC, fold balanced accuracy, fold complexity, no-call counts, and timing;
- `models.lasso` and the new solver-status fields;
- `features` and annotations, because the list includes the top 30 lasso weights (newly included loci need annotation downloads);
- best-learner labels and corresponding log wording;
- `evaluation.revision`, report validity evidence, timing fields, and the report index.

SCM, CART, forest, boosting, MIC, catalog, and probe outputs should remain identical. The earlier cross-check reproduced all 278 non-lasso rows; the full-run diff must enforce that rather than assume it. Boosting is unaffected by this solver fix. XGBoost and the default `colsample = 0.5` remain a separate verification gap and do not block this revision-8 prerequisite work.

## Large-species lasso scan

The new limits permit 25 × 100,000 sweeps per point, or 75 million sweeps for a 30-point path: 333× the old theoretical worst case. The observed 1.92× runtime covers only four small datasets, and intermediate cap hits persist on the motivating C. jejuni panels. No large-species runtime bound has been measured. The scan measures one, under the rule registered here, before any 275-report run. **Running it needs separate user approval.**

### The harness

`scripts/lasso-scan/` runs the shipped revision-8 code, not a copy. The earlier `lasso_deploy.rs` harness, which replayed a copied revision-7 `fit_lasso`, is not used.

- `main.rs`, built as the Cargo example `lasso-scan`, calls `grm::analysis::lasso_scan`. That function loads a species with `load_inputs` and builds each dataset's cohort, projection and lineages with `eval_inputs`, the same function `analyze_dataset` uses. It reports the dataset's input fingerprint, then runs `eval::lasso_only`.
- `eval::lasso_only` runs only the lasso fits of `evaluate`, in its order: each available split's outer folds, then the final model. It records each penalty point's `Convergence`, for the tuning path and the refit.
- `scan.py` is the driver. It checks the list against each dataset's revision-7 report, runs the binary on one dataset at a time, and enforces the rule below. It writes only its output directory: `plan.json`, `events.jsonl`, one stderr log per dataset, `summary.json`, and `stop.json`.
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

At the measured 1.92× the scan should take about 2.1 h of lasso work, plus loading.

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
- at most 5% of the solved penalty points, tuning and refit together, hit the sweep cap in some IRLS step. Above that, a review comes first. The saved fixtures' rate is 4 of 84.

The projection is (25,088 s + (median r − 1) × 1.1 × 14,989 s) ÷ 3,600, where 25,088 s and 14,989 s are the 275 revision-7 reports' total report and lasso seconds. It gives 11.2 h at 1.92×, 16.1 h at 3× and 25.3 h at 5×.

`scan.py` exits 0 for go, 1 for no-go, 2 when the rule stops the scan (`stop.json` says why) and 3 when it refuses to start.

### Disk

The scan copies no input and writes no temporary file. It reads about 4.4 GiB in place under the input root: the five species' matrices, counts and structures, the curated phenotypes and the BV-BRC genome list.

Each fit line it writes holds at most 30 tuning and 30 refit statuses and one weight per screened feature, so it stays under 256 KiB. Add a 1 MiB stderr allowance per dataset and 4 MiB of fixed files, and the 18 datasets need at most **75,249,664 bytes (71.8 MiB)**. For scale, the synthetic smoke run wrote 51,444 bytes for one dataset. Staging the release binary (7,725,344 bytes as built for `znver2` from this lane's source, needing glibc 2.34 or later), `scan.py` and the list adds 7,759,021 bytes (7.4 MiB). The scan therefore needs at most **83,008,685 bytes (79.2 MiB)** of new disk in all; the inputs-only check writes a few kilobytes more. Building on Vast instead would add a release target directory of about 443 MB, which is why the binary is built elsewhere.

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

1. Before the scan: rerun the release lasso tests and the r8kq3v `xcheck.sh` on the frozen commit, with a separate Python interpreter rather than the hard-stopped GWAS lane's. Then run and review the scan under this rule.
2. Decide how to respond if a lasso error or an invalid report stops `grm learn` for a species, since no downstream consumer handles fewer than 275 reports.
3. Only after the scan passes, run the full 275-report evaluation with revision-8 source into a **new data root**, leaving revision-7 candidate, backup, and published roots untouched.
4. Require every report to pass the revision-8 checker with zero unconverged two-class lasso fits, and run `metric-audit --revision 8`. Then do the field-level migration diff and the completion-log audit.
5. Rerun the downstream consumers listed below.
6. Rebuild and review the report/web bundle. Publication is a further explicit decision.

`models.lasso_fit` sits between `models.lasso` and `rf_importance`, so whole revision-8 files will not be byte-identical to revision-7 ones. The migration diff has to compare non-lasso fields field by field.

## Deployment settings and downstream work

The current published release remains revision 7. `scripts/deployment/` retains its revision-7 defaults and 237-dataset published-release expectations; this lane deliberately did not edit it while its independent review is active.

At a future, authorized revision-8 staging/publishing operation, the operator must explicitly set:

```bash
GRM_EXPECT_REVISION=8
GRM_EXPECT_DATASETS=275
```

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
