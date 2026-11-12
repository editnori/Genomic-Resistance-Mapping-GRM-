# Every learner scored on the same genomes (revision 14)

From evaluation revision 14, each learner row of a dataset report covers every genome of the split. Before it, revisions 7 to 13 left the genomes of a null-fit fold without a call, so a learner with a null fold was scored on fewer genomes than the others, and its balanced accuracy could not be compared with theirs.

## The rule

A null fit is an outer fold whose fitted model learned nothing beyond the class balance: a boosted model that keeps no tree, or a lasso that keeps no weight (`fold_complexity` 0). The other three learners always fit.

- **Score.** Each test genome of a null-fit fold gets that fold's training prevalence as its score: the resistant share of the samples the fold trained on (`learn::train_prevalence`, one definition shared by `learn::eval`, `learn::planned` and `ablate::arm`). It is not the whole dataset's prevalence and not the test fold's own share.
- **Call.** The score is called at the usual 0.5 cutoff with `>=`, so a training prevalence of exactly 0.5 calls every genome of the fold resistant, as a fitted score of 0.5 does.
- **Metrics.** `metrics`, `roc`, `fold_balanced_accuracy` and everything derived from them cover every genome. A null fold's balanced accuracy is that of its constant call, 0.5 when its test fold holds both classes. AUROC ranks the null fold's genomes with the fitted folds' scores, ties at half credit.
- **`no_call`.** Keeps its name and meaning on a learner row: the count, by class, of genomes in null-fit folds. They now have a score and a call.
- **`called_metrics`.** A new optional field, written only on a row that has such genomes. It holds the revision-13 figure: the metrics over the genomes of fitted folds only, or `null` when every fold is null. A row with no null fit has no `called_metrics`, and its numbers are the ones revision 13 wrote.

Nothing else changes: no fit, fold, learner, hyperparameter or tuning (production tunes inside each fit, as before). The report schema stays 2, since the new field is optional and a schema-2 reader that ignores it reads every number right. A same-fingerprint revision-13 report is stale against revision 14 (`stale_reason`).

## What follows from it

- **Best learner.** `report::best_learner` picks the best random-split balanced accuracy among learners scored on every isolate on every split that ran. From revision 14 a row with `called_metrics` (index flag `null_fits_scored`) counts as scored on every isolate, so a learner with a null fold now competes. On older rows the rule is unchanged: a learner with a no-call is left out.
- **Input comparison (H-1).** `grm arm` reads `called_metrics` when present, else `metrics`, and writes a null fold's prediction rows as before: score the training prevalence, `called` false. Its arms, statistics and prediction rows are therefore unchanged by this revision (`ablate::arm::tests::run_metrics_and_prediction_rows_keep_a_null_fits_samples_uncalled`).
- **Conformal and threshold diagnostics.** No report carries them: `pool_folds` writes an empty `conformal` list, `learn::eval::conformal` is called only by its test, and no code sweeps the cutoff over a report's scores. Nothing there changes. A future diagnostic that reads `scores` gets every genome scored, the null folds at a constant.
- **Checkers.** `scripts/check-results.py` and `scripts/metric-audit/audit.py` accept revision 14 and check, on each learner row: no null score; each null fold (by `fold_complexity`) scored at its training prevalence; `no_call`; `metrics` over every genome; `called_metrics` present exactly when `no_call` is non-zero and equal to the fitted folds' figure; each fold's balanced accuracy. `nullimpact.py` and `transitions.py` find null folds by null scores, so they refuse revision 14.
- **Dashboard.** A revision-14 report counts no learner no-calls. The learner table's Null folds column counts the null folds, and its note says each one predicts its training isolates' resistant share and that every rate includes them. The best-learner wording reads "scored on every isolate" instead of "calls every isolate".
