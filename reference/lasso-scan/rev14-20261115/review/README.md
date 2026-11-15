# Independent review of the revision 14 lasso scan

On 15 November 2026 a second lane on the same machine reviewed the scan in the folder above, independently of the lane that ran it. It recounted every figure from the scan's files with code of its own, refitted 3 of the datasets with its own build of the same commit, and checked the one-fit measurement of MTB delamanid random fold 0 again. Its verdict: every count reproduces, and the data support a narrower reading than some of the scan's sentences. Its findings are named F1 to F5 (should-fix) and N1 to N9 (notes) on the page. The section "Revision 14 scan" of [docs/lasso-stopping.md](../../../../docs/lasso-stopping.md) carries the narrower reading.

This folder holds the review's result files that the page quotes, the scripts that wrote them, and the review's plan with its sha256. Every file here except this README was copied unchanged from the review's output folder, keeping its paths, and `../SHA256SUMS` gives the sha256 of each. None has one row per sample or genome.

## Files

In the commands, `L` is the scan lane's output folder, `D` the scan's pass 2 dump folder, `REL` the unpacked release `rev14-274-20261112` and `S` the review's scratch folder. As in the folder above, the scripts are the record of what ran, not a kit to run as they stand.

| File | What it is | Written by |
|---|---|---|
| `PLAN.md` | The review's plan for its replication of 3 datasets, chosen by release fields alone: what is compared, the measures and the pass rule | the review, before its runs |
| `results/plan.sha256` | `PLAN.md`'s sha256 and the time it was recorded, 16:35:52 UTC, before the replication ran | `sha256sum PLAN.md` |
| `results/recount.json` | The recount of the whole scan: per dataset and split, fits, penalty points, IRLS steps, scores compared with the release, calls changed, and recall, log loss, balanced accuracy and AUROC from the scan, from the release's scores and from the release's own fields; the totals; and the four folds whose refit keeps no weight, with their training prevalence | `python3 -I scripts/recount.py $L/results/scan $D $REL $S/src/scripts/lasso-scan/scan-list.json results/recount.json` (standard library only) |
| `results/testrecall.json` | Recall beside log loss on the 250 test samples of MTB delamanid random fold 0, from grm's refit scores: the counts of true and false calls, recall, specificity and log loss | `python3 -I scripts/testrecall.py $D/mycobacterium_tuberculosis__delamanid__random__0.json results/testrecall.json` |
| `results/item5check.json` | The one-fit measurement checked again: the near twins recounted with integers as bit patterns, and per path point the certificate's dual point rebuilt and checked against the dual box, the bound after the point is shrunk into the box, and a second certificate that starts from zero rather than from grm's solution | `scripts/item5check.py DUMP.json LANE_RESULTS_ITEM5_DIR LANECODE_DIR OUT_JSON SCRATCH_DIR`, run with the scan's solver environment; `LANECODE_DIR` holds copies of the scan's `item5.py` and `item5-polish.py` (the second named `item5_polish.py` there), equal by sha256 to those in `../scripts/` |
| `results/replication.json` | The replication of `PLAN.md` for `mycobacterium_tuberculosis__clofazimine`, `salmonella_enterica__ampicillin` and `escherichia_coli__ampicillin`: inputs, events, dumps, scores and calls against the scan and the release | `python3 -I scripts/rep_compare.py $S/runs $L/results/scan $D $REL results/replication.json mycobacterium_tuberculosis__clofazimine salmonella_enterica__ampicillin escherichia_coli__ampicillin` |
| `scripts/recount.py`, `scripts/testrecall.py`, `scripts/item5check.py`, `scripts/rep_compare.py` | The four scripts, as the review left them. One line of `item5check.py`'s docstring was reworded after its run; its code is unchanged | the review |

## What the page quotes from these files

| Page sentence | File and field |
|---|---|
| 57,573 fitted scores and 1,045 at training prevalence, and the four folds that keep no weight | `results/recount.json`: `totals.scores_from_binary`, `totals.scores_from_binary_identical`, `totals.null_fold_scores`, `null_folds` |
| Test recall 0.4111 (37 of 90 resistant) beside test log loss 0.670059 on MTB delamanid random fold 0 | `results/testrecall.json` |
| The certificate's dual point outside the box at 13 of the 30 points, by a factor of at most 1 + 3.18e-13; the shrink changes the bound by at most 4.3e-15; the largest bound after it, 1.29906e-9 at point 29; the same largest bound from the certificate that starts from zero | `results/item5check.json`: `certificate[].rebuilt_check.max_ratio`, `rebuilt_D_minus_lane_D2`, `rebuilt_bound`, `zero_start_bound` |
| The replication of 3 datasets, equal to the scan and the release | `results/replication.json` |

## Sentences of the copied scan files that the review narrowed

The scan's `PLAN.md` and `results/tables.md` above stay as they were written. Their sentences below say more than the review found, and the page uses the corrected form. The folder's own `README.md` was written for the repository, so its two sentences are corrected in place.

| File and line | As written | Corrected form |
|---|---|---|
| `PLAN.md:22` | "My driver `scripts/scan14.py` runs the same binary and applies the same per-fit test ... and the same input check, against the revision 14 report's input fingerprint instead of revision 7's." | The driver applies the same per-fit test. Its input check is a different check of the same field: it compares the fingerprint with the revision 14 report's after the run and stops nothing, and the budgets are computed after the run at h = 1 and stop nothing. "18 of 18 datasets pass the registered scan's convergence test, applied by a different driver; `scan.py` itself did not run." (N1) |
| `PLAN.md:38` | "a fold whose refit keeps no weight is scored at its training prevalence, as revision 14 does" | Those scores come from no fit. "All 57,573 out-of-fold scores from fits that keep a weight are bit-identical to the release, and the 1,045 scores of the four folds that keep none equal their training prevalence in both." (F3) |
| `PLAN.md:40` | "**Released measures compared, per split**: calls ..., recall (sensitivity), log loss, ..." | "The release has no log loss field; log loss computed from the identical scores is equal." It is not a separate check. (N2) |
| `PLAN.md:42-49` | "When the released lasso figures need a caveat", and its conclusion that the figures "are reproduced by an independent rerun of the same code on the same inputs" | The rule tests reproduction and solver stopping only. "The released lasso figures of these 18 datasets are reproduced bit for bit by a rerun of the same code on rebuilt inputs, and every fit met the solver's convergence test. This says nothing about the lasso's accuracy or calibration." (F2) |
| `PLAN.md:57` | The planned dual point "scaled into the box", "giving a certified lower bound D on the optimum" | "On the one fold measured, the planned dual construction gave no bound at the 14 smallest penalties and one bound above 1e-6; the exploratory construction bounded every point after a final shrink into the dual box." (F5) |
| `PLAN.md:58` | "calls that differ ..., the largest score change, and each side's test log loss" | Recall beside log loss: "test recall 0.411 (37 of 90 resistant) and test log loss 0.670 for both solves." (F4) |
| `PLAN.md:60` | "Item 5 reading: the page's claims hold at revision 14 if every path point's duality-gap bound is at most 1e-6 ... and no held-out call differs." | This reading failed: 1.10e-6 at point 11 and no bound at points 16 to 29. The certificate that replaced it was chosen after that result. (N6) |
| `results/tables.md:22` | The column "Largest change in a released measure", whose maximum includes log loss | Log loss is not a released field; with identical scores it is equal by construction. (N2) |
| `results/tables.md:43` | "Scores equal to release", for example "1248 of 1248" for `MTB__linezolid` lineage | The count includes the scores of the four folds that keep no weight, which are training prevalence, not a fit's: 250 of that split's 1,248. (F3) |
| `results/tables.md:82` | The columns "P_grm - P_ind" and "Bound for grm, exploratory" | The first rests on one solver, Clarabel, and at its 4 `optimal_inaccurate` points it is not a gap to an optimum. The second is the bound before the shrink into the dual box; after it, the largest is 1.29906e-9 at point 29. (N4, N5) |
| `README.md:5` | "Equal scores show that this build and the release compute the same lasso on these datasets." | Corrected in place: 57,573 of the equal scores come from fits, and 1,045 are training prevalence in folds that keep no weight. (F3) |
| `README.md`, the row of `results/item5/item5-polish.json` | "a certified optimum at each path point, and grm's distance above it" | Corrected in place: its dual point lies outside the dual box at 13 of the 30 points and certifies only after a final shrink into it. (N4) |

## Not kept

- The review's report, its builds, its replication runs and its logs. They stayed on the machine that ran the review.
- The dump that `scripts/testrecall.py` and `scripts/item5check.py` read, as in the folder above.
