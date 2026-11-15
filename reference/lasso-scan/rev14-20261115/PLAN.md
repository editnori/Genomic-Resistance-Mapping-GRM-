# Lane lassoscan: plan, fixed before any result is read

Written 2026-11-15 between 14:10Z and 14:20Z. Its sha256 is recorded in `REPORT.md` and in `results/plan.sha256` before the scan or the dump binary runs. Anything decided after a result is read is labelled exploratory in `REPORT.md`.

## Sources and code

- Repository `main` at `852ebc9`, exported with `git archive` into `/home/user/scratch/lassoscan/src` and built there with Rust 1.98.1 and the repository's `.cargo/config.toml` (`x86-64-v3`): `grm` and the example `lasso-scan`. `git diff --stat 469a4f9 852ebc9 -- src/ Cargo.lock Cargo.toml .cargo/ scripts/lasso-scan/` is empty, so this is the source the revision 14 release was built from (`469a4f9`, release notes).
- Inputs: the public data root `/home/user/grm-h1raw`, whose checksum list `MASTER-MANIFEST.sha256` has sha256 `7dd00779...7887`, the list the revision 14 release names. Matrices are built with `grm count`, `grm structure`, `grm matrix` at default settings (k = 31) into `/home/user/scratch/lassoscan/data`, which links the read-only inputs.
- Release compared against: `/home/user/grm-ops/owner2/work/rev14-full/rev14-274-20261112/datasets/`.

## Species and datasets

The scan's set is the registered list `scripts/lasso-scan/scan-list.json`: the five species with the most revision 7 lasso seconds (*M. tuberculosis*, *K. pneumoniae*, *S. enterica*, *E. coli*, *S. aureus*), and in each its three costliest revision 7 datasets plus its datasets with the most samples and the most features: 18 datasets. All 18 are in the revision 14 release. "Large species" means these five, as the list defines them.

**Order, in case time runs out.** Species in the order each first appears in the registered list (most expensive revision 7 dataset first): *M. tuberculosis*, *S. enterica*, *S. aureus*, *K. pneumoniae*, *E. coli*. Within a species, the registered order. The first dataset, `mycobacterium_tuberculosis__delamanid`, is also item 5's. A species not reached is reported as not scanned.

## Settings at revision 14

- `lasso-scan` built from `852ebc9` (`EVALUATION_REVISION` 14, `EvalConfig::default()`: 5 outer folds, 4,000 screened features, 30 penalties from lambda_max down to 1e-3 of it, path stop above 400 non-zero, seed 42, holdout 0.2). Each dataset runs the 10 outer fits (random and lineage splits, 5 folds each, or fewer if the binary reports a split unavailable) and the final fit, in `evaluate`'s order.
- One dataset per binary call, as `scan.py` runs it. `RAYON_NUM_THREADS=16`: this machine has 16 cores and no other lane (revision 7 and the registered rule used 64 on a 64-thread host).
- Pass 1 (timed) runs alone: no matrix build or other job of mine runs during it. `nproc`, CPU model and load average are written beside it.
- `scan.py` itself cannot run here: it refuses unless each dataset's revision 7 report is present with the revision 7 input fingerprint, and the revision 14 cohorts differ. My driver `scripts/scan14.py` runs the same binary and applies the same per-fit test (`unconverged(audit)` of `scan.py`, which is `LassoFitAudit::has_unconverged`) and the same input check, against the revision 14 report's input fingerprint instead of revision 7's.

## Convergence: the scan's own rule

A fit converges when the binary returns no lasso error and its audit has `path_unconverged == 0` and converged selected and refit points. A dataset converges when every one of its fits converges and the binary exits 0 with the expected number of fits. Recorded for every fit: path points, unconverged points, selected point, IRLS steps per point (selected, refit, largest, total), `sweep_cap_hits`, seconds.

The registered rule stops the scan at the first unconverged fit or budget overrun. Here the scan does not stop: it continues to the next dataset and records where the registered rule would have stopped, because the aim is a per-dataset answer and the 274-report run it gated has already run. Budgets: the registered fit and dataset budgets at h = 1 are recorded and compared, but they stop nothing; h for this host cannot be measured (the revision 7 host is gone). Each dataset has an operational cap of 30 minutes of wall-clock, after which it is killed and reported as not finished.

The registered go/no-go (worst ratio at most 5, median at most 3, projection at most 16 h, at most 5% of points with `sweep_cap_hits > 0`) is computed at h = 1 and reported as exploratory, since the host differs.

## The release's fit against the scan's

For every dataset:

1. Input fingerprint from the binary equals the released report's `input_fingerprint`; samples and features equal the report's `counts`.
2. Each fold's audit (path points, unconverged, selected point, selected and refit convergence) equals the report's `lasso_fits` entry, and the refit's weight count equals `fold_complexity`. The final fit's audit equals `models.lasso_fit`, and its model (intercept, every coefficient's feature and weight) equals `models.lasso`, compared as JSON numbers.
3. Pass 2 runs the instrumented binary (scratch copy `src-dump`: one function added to `src/learn/eval.rs` that writes each fit's partition, path, refit and scores when `GRM_LASSO_DUMP` is set; nothing else changed). Its audits, traces and models must equal pass 1's exactly, or its scores are not the scan's. From its scores, each split's out-of-fold scores are assembled (test samples = samples outside the fold's training mask; a fold whose refit keeps no weight is scored at its training prevalence, as revision 14 does) and compared with the released lasso `scores`.

**Released measures compared, per split**: calls (score >= 0.5 is resistant), recall (sensitivity), log loss, balanced accuracy, VME, ME, AUROC. Log loss is the mean binary log loss of the pooled out-of-fold scores, probabilities clipped to [1e-15, 1 - 1e-15] (the definition of `scripts/hyperopt/canaries.py`); the release has no log loss field, so both sides are computed the same way from scores. Recall is shown beside log loss for each dataset, for the lasso (the only learner the scan fits).

**When the released lasso figures need a caveat.** Any one of:

- a listed dataset fails the scan's convergence rule (any unconverged path point, selected point or refit, or a lasso error), since its released report says valid;
- the scan's input fingerprint differs from the released one;
- a fit's audit, weight count or final model differs from the release;
- a held-out call differs, or recall or log loss differs by more than 1e-9, when recomputed from the scan's scores.

If none holds, the released lasso figures of the scanned datasets are reproduced by an independent rerun of the same code on the same inputs, and the difference in every released measure is reported as 0.

## Item 5: the lost diagnosis, measured again at revision 14

Dataset `mycobacterium_tuberculosis__delamanid`, split random, outer fold 0, revision 14 cohort. The October figures are history whose files are lost; they are shown beside, never mixed in.

1. **Near twins.** Among the screened columns of fold 0's tuning path (the `features` the path is fitted on), on the fold's fitting samples (the tuning holdout's fitting mask), count the columns that have another screened column differing from them in at most one fitting sample (Hamming distance 0 or 1). Also reported: the number of fitting samples, screened columns, columns with an exact copy, redundant copies (columns whose lowest-index identical column is another, which `representatives` collapses), and columns present on every fitting sample. October: 2,820 of 4,000 on 798 fitting samples; 410 redundant copies.
2. **Time of the tuning path.** From pass 1's stderr progress lines: the sum of the path points' seconds, the number of points, and the fit's seconds (path plus refit, from the fit line). A second run of this dataset alone (`scripts/time-delamanid.sh`), so there are two timings as in October. Machine: `nproc`, CPU model, `RAYON_NUM_THREADS`, load average before and after. October: 13.4 s and 11.5 s on an i5-10400 with 12 threads at load 13 to 21.
3. **Gap to an independent solve.** For every path point, with grm's lambda: the objective (mean class-weighted log loss on the fitting samples plus lambda times the sd-weighted L1 norm, grm's definition in `src/learn/lasso.rs`) evaluated in numpy at grm's model, P_grm, and at the solution of an independent solver, P_ind. The independent solver is Clarabel 0.9.0 through CVXPY 1.6.0, installed in scratch from PyPI wheels pinned by version and sha256 (checked against PyPI's published digests), with gap and feasibility tolerances 1e-10. Columns with no carrier on the fitting samples, or present on all of them (which grm clears), are left out of the independent problem; exact copies are kept, so the collapse is tested, not assumed. Reported per point: P_grm - P_ind; the KKT residual of the independent solution by grm's definition, after coefficients with |b| x sd below 1e-8 are set to zero (the objective at that thresholded point is reported too); grm's own KKT residual; and a duality-gap bound: from each solution a dual-feasible point is built (the residual vector, centred to sum zero with the class weights, scaled into the box |X_k' u| <= lambda sd_k), giving a certified lower bound D on the optimum; P_grm - max(D) bounds how far grm's point is above the optimum, whatever any solver did. Second solver, if time allows: scikit-learn 1.5.2 `LogisticRegression(penalty="l1", solver="saga")` on sd-scaled columns with class weights as sample weights, time-boxed to 30 minutes in all; its P - P_ind is reported with its status.
4. **Held-out calls.** The independent selection takes the path point with the least class-balanced validation log loss (grm's `lasso::log_loss`, clamp 1e-9) among the independent points, and solves the refit at that lambda on the outer training fold (sd and class weights on the training fold). Its scores on the fold's test samples are compared with grm's refit scores from the dump: calls that differ (score >= 0.5), the largest score change, and each side's test log loss. If the two selections differ, the comparison is also made at grm's selected lambda. October: point 3, 53 refit weights, 0 of 251 calls changed, largest score change 1.7e-5.

Item 5 reading: the page's claims hold at revision 14 if every path point's duality-gap bound is at most 1e-6 (the KKT tolerance and the ignored test's acceptance) and no held-out call differs. Whether the page should quote revision 14's numbers is the owner's decision; I change no document.

## Outputs

`results/scan/` (events and stderr per dataset, `summary.json`), `results/compare/` (per-dataset comparison with the release), `results/item5/` (twins, timing, per-point table, calls), `scripts/`, `MANIFEST.sha256`. Large dumps stay in scratch; their sha256 go in the report.
