# Null-fold refit experiment (off by default)

Revision 14 scores an outer fold whose boosted model keeps no tree, or whose lasso keeps no weight, at the fold's training share of resistant samples ([equal-scoring.md](equal-scoring.md)). This experiment asks whether such a fold hid a usable model or had no signal for that learner. It is an experiment flag. Without it nothing changes.

## The flag

```
grm learn --species SLUG --only DATASET[,DATASET...] --force --exp-null-refit FILE
```

With `--exp-null-refit FILE`, `grm learn` runs exactly as without it: the same fits, the same scores, the same reports. In addition, for each boosting or lasso outer fold whose fit is null, it refits the fold's training samples under the fallback settings below and appends one JSON line to FILE. The evaluator calls the experiment after the fold's production fit and outside its timer (`learn::eval::evaluate`), and only for a null fit.

`EVALUATION_REVISION` and the report schema do not change, because no report field changes. `tests/null_fold_refit.rs` checks that an evaluation with the flag equals one without it on every field except wall-clock `seconds`.

## What each fallback refits

Both redo the production fit's own steps on the outer fold's training samples, with the fold's seed, the same tuning holdout and the same chi-square-screened features (`learn::null_refit`).

- **Lasso.** The same tuning path (30 penalties from the largest that keeps every weight at zero down to 1e-3 of it), then the refit on the whole training fold down through every path point. Production selects the point of lowest class-balanced holdout log loss; on a null fold that is point 0, where every weight is zero. The record keeps every point's holdout log loss and balanced accuracy and every point's refit scores on the fold's test samples, so a fixed rule can be read off afterwards:
  - `smallest_point`: the last point converged on both the path and the refit, the smallest penalty reached.
  - `next_best_point`: the production rule over the converged points after point 0 that keep a weight.
- **Boosted trees.** The same fit with early stopping off (`gbt::fit_traced` with `early_stop` false): all `rounds` trees on the fitting samples, with the class-weighted holdout loss of the 0.5 prior and after each round. `best_round` is the number of trees early stopping would keep with no patience limit (`gbt::best_round`), 0 when no round beats the prior. The record has the test scores of the model with every tree, and of the first `best_round` trees when that is not 0.

`gbt::fit` is `fit_traced` with early stopping on, so the production model is unchanged.

## A record

One line per null fold:

| Field | Meaning |
|---|---|
| `dataset`, `scheme`, `fold`, `learner` | the outer fold |
| `train_n`, `train_resistant` | its training samples |
| `test` | its test samples, every sample outside the training fold, in order |
| `lasso.holdout_n`, `lasso.holdout_resistant` | the tuning holdout |
| `lasso.production_point` | the point production selects, 0 on a null fold |
| `lasso.smallest_point`, `lasso.next_best_point` | as above, or null |
| `lasso.points[k]` | `lambda`, `path_weights`, `path_converged`, `holdout_log_loss`, `holdout_balanced_accuracy`, `refit_weights`, `refit_converged`, `test_scores` |
| `gbt.holdout_n`, `gbt.holdout_resistant` | the early-stopping holdout |
| `gbt.holdout_loss` | prior, then after each round (`rounds` + 1 values) |
| `gbt.best_round`, `gbt.full_test_scores`, `gbt.best_round_test_scores` | as above |

## Costs

The experiment adds one lasso path and refit, or one boosted fit of `rounds` trees, per null fold. A dataset with no null fold costs nothing more.

## What it is not

It is not a production setting and selects nothing for a report. A fallback that scores better on some folds would still need its own plan and a revision before it could change a report. Every reading is cross-validation within the release's genomes.
