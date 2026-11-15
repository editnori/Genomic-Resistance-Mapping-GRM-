# Fold draws

A release scores each dataset on one assignment of genomes to folds. `grm learn --seed` scores the same dataset on another, so the spread of a figure across draws can be measured. `grm pipeline`, the release's command, takes no seed and always draws at 42.

## What the seed moves

`EvalConfig::seed` (42 in every release) is the one source of randomness in an evaluation. It draws:

- both splits' outer folds: the stratified folds of the random split and the grouped folds of the lineage split (`learn::eval::outer_splits`);
- inside outer fold `f`, at seed + 101 `f`: the inner folds SCM and CART choose their settings on, boosting's early-stopping holdout and the lasso's penalty holdout;
- the forest's subsamples and node features and boosting's column samples, at the same per-fold seed;
- the final models' tuning partitions and draws (seed + 7, 11, 13 and 17);
- the MIC regression's folds, holdouts and boosting draws.

It does not move the species matrix, the dataset's genomes and labels, the support floor, the lineage clusters or the duplicate groups. None of them reads it.

The chi-square screen and the per-fit candidate features have no draw of their own. Each is a fixed function of a fit's own fitting genomes, with ties broken by feature index, so they follow the folds and nothing else. The same holds for every setting a learner selects inside a fit.

So another seed is another draw of the folds on the same features and cohort, and it also redraws the forest's and boosting's own randomness. `--learner-seed` separates the two: it seeds the forest's and boosting's own draws and leaves every partition as `--seed` draws it. SCM, CART and the lasso have no draw of their own, so a run that changes only `--learner-seed` reproduces their scores exactly.

## Where the reports go

A run at the release's seeds writes `<data>/species/<slug>/results/`, as before. Any other run writes beside it:

| Run | Folder |
|---|---|
| `grm learn --species S` | `species/S/results/` |
| `grm learn --species S --seed 43` | `species/S/draws/seed-43/` |
| `grm learn --species S --learner-seed 7` | `species/S/draws/seed-42-learner-7/` |

`grm report`, `grm intervals` and `grm provenance` read `results/` alone, so a draw is never indexed or served, and it never replaces a released report. Each draw's report records its seed in `eval_config.seed`, and `eval_config.learner_seed` when one was set. The field is left out when it is not set, so a release report keeps its bytes and its input key.

A draw needs what `grm learn` needs: the species counted and built, and its k-mer sets on disk for a species with a probe assay. A new draw can leave a fold or a tuning partition without a class on a small dataset. The lineage split is then recorded as unavailable in that draw's report. An invalid random split is an error, as it is at seed 42 (docs/duplicate-isolation.md): the species' run stops at that dataset and writes no report for it. The datasets it had not reached can be run one at a time with `--only <id>` and the same seeds.

## Repeat seeds and the input comparison

The input comparison's `grm arm --repeats` takes repeat `r`'s seed from `ablate::plan::repeat_seed(42, r)`. `grm learn --seed` can take the same values, so a set of draws can use the arm's repeat seeds:

| Repeat | Seed |
|---|---|
| 0 | 42 |
| 1 | 10622993013839509822 |
| 2 | 6315510376195533293 |
| 3 | 1291059630886940069 |

The seeds match; the folds do not. A draw of `grm learn` at a repeat's seed does not score the partition the arm scores at that repeat, at any repeat, repeat 0 included:

- **Folds keyed on genome identity.** A reference arm draws each repeat's plan with `ablate::plan::keyed_fold_plan`, which sends a genome's unit to the outer fold given by a hash of its key and the repeat seed. `grm learn` draws its folds over sample positions (`stratified_folds` and `grouped_folds`, which `ablate::plan::fold_plan` reproduces). The two give different partitions from the same seed.
- **The random split at repeat 0 only.** The arm scores the lineage split at every repeat and the random split at repeat 0 alone (`ablate::arm::split_repeats`).
- **Replaced repeats.** A lineage repeat the arm refuses is replaced by a later keyed draw (`ablate::arm::MAX_REPLACEMENTS`), so the arm's draw `r` need not be repeat `r`.

A figure from a `grm learn` draw and a figure from an arm repeat with the same seed are therefore not paired on one partition. The fold-draw experiment of pull request #115 used the repeat seeds only to choose its seeds, and paired none of its draws with an arm.

## Evaluation revision

Neither flag changes a report at its default, so `EVALUATION_REVISION` and the report schema stay as they are.
