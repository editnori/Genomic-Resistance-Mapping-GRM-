# Null-fold refit experiment (off by default)

Revision 14 scores an outer fold whose boosted model keeps no tree, or whose lasso keeps no weight, at the fold's training share of resistant samples ([equal-scoring.md](equal-scoring.md)). This experiment asks whether such a fold hid a usable model or had no signal for that learner. It is an experiment flag. Without it no output changes; the time can (see Costs).

## The flag

```
grm learn --species SLUG --only DATASET[,DATASET...] --force --exp-null-refit FILE
```

With `--exp-null-refit FILE`, `grm learn` writes the same reports as without it, as long as FILE takes every line: the same fits, the same scores. In addition, for each boosting or lasso outer fold whose fit is null, it refits the fold's training samples under the fallback settings below and appends one JSON line to FILE. The evaluator calls the experiment after the fold's production fit and outside its timer (`learn::eval::fit_outer_folds`, which `evaluate` calls for each split), only for a null fit, and with the fold's own seeds: `--seed` for every partition, and `--learner-seed`, when given, for boosting's column samples.

The flag is refused beside `--fold-file` ([fold-file.md](fold-file.md)): the argument parser stops before FILE is opened, and `learn::eval::evaluate_assigned` refuses to run while the experiment is on. A record names every sample outside a fold's training samples as its test samples, which a fold file's unused genomes would make wrong, and the records are read against the drawn splits of the dataset reports.

`EVALUATION_REVISION` and the report schema do not change, because no report field changes. `tests/null_fold_refit.rs` checks that an evaluation with the flag equals one without it on every field except wall-clock `seconds`, and that FILE keeps the lines it held. A test in `learn::eval` checks that each record is what the fallbacks fit when called directly with that fold's candidates, seeds and tuning groups.

## FILE

- FILE is opened for appending, and created if missing, when the learn starts, before any species is read. So it exists, perhaps empty, even when the learn then fails.
- A line FILE does not take stops the learn, because a skipped record would bias every count taken from FILE. The error names FILE, the dataset, the split, the learner and the fold, for example `appending the null-fold record of DATASET, lineage split, gbt, outer fold 1 to FILE: No space left on device (os error 28)`. That dataset's report is not written. Datasets run one at a time, so the reports finished before it are written, and later datasets and species do not run. A write cut short can leave part of a line at the end of FILE, so an analysis must refuse a line that does not parse.
- A refit that fails stops the learn the same way, with an error that names the dataset, the split, the learner and the fold.
- A record names no run and no seed. A second run of a dataset into the same FILE adds a second record for each of its null folds, so an analysis must refuse a repeated key (`dataset`, `scheme`, `fold`, `learner`). A run at another `--seed` or `--learner-seed` needs a FILE of its own.
- Without `--force`, a dataset whose report is up to date is skipped and writes no record. That is why the usage line has `--force`.
- One process writes one FILE. Inside it, each record is one write under a lock. Nothing locks FILE between processes, and one record can run to hundreds of kilobytes, so two processes must not write the same FILE.

## What each fallback refits

Both redo the production fit's own steps on the outer fold's training samples, with the fold's seeds, the same tuning holdout and the same chi-square-screened features (`learn::null_refit`).

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

The experiment adds one lasso path and refit, or one boosted fit of `rounds` trees, per null fold. It runs outside the fold's timer, so a report's `seconds` leave it out.

Without the flag no output changes, but the time can: the production lasso's speed turns on how the compiler splits the crate into codegen units. On `salmonella_enterica__tetracycline` with the flag off (15 Nov 2026, release builds, 5 threads, three runs of each build in turn), this branch after its merge of `main` at `41cd616` ran the lasso in 30.0 to 31.1 s, `main` at `41cd616` in 29.9 to 31.2 s, and `main` at `7020d4c` in 7.9 to 8.1 s. So `main` itself slowed between `7020d4c` and `41cd616`, and the branch adds no lasso time that these runs separate from `main`'s. Earlier that day, the review of this branch timed it 3.9 times slower than `main` at `852ebc9` (30.2 s and 32.4 s against 7.8 s and 8.0 s), and found that `main` built with one codegen unit was as slow (31.2 s): the cause is the split into units, not the experiment's code.

## What it is not

It is not a production setting and selects nothing for a report. A fallback that scores better on some folds would still need its own plan and a revision before it could change a report. Every reading is cross-validation within the release's genomes.
