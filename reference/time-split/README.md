# Time split record

This directory holds the record of the split by first-public day inside BV-BRC, rerun after the review. The page that reads it is [`docs/time-split.md`](../../docs/time-split.md). `scripts/time-split/record.py` writes every file here from the input set, and its `--check` mode fails on any byte it would not write.

## Files

| File | Rows | Contents |
|---|---|---|
| `fits.tsv` | 360 | Row file. One row per dataset, arm and learner |
| `datasets.tsv` | 24 | Row file. One row per dataset |
| `species.tsv` | 4 | Row file. One row per species |
| `review.tsv` | 5 | Row file. One row per learner: what the review's fix changed |
| `totals.tsv` | 7 | Row file. Counts over all datasets |
| `order.tsv` | 7 | Row file. One row per event of the plans and runs |
| `settings.tsv` | 24 | Row file. The design's minimums, the fits' settings, the settings the script applies, and the sha256 of the two plans and two scripts |
| `checks.tsv` | 10 | The checks the build ran on the input set |
| `inputs.sha256` | 155 | The sha256 of each file of the input set |
| `contrasts.tsv` | 360 | Derived. One row per dataset, learner and contrast |
| `readings.tsv` | 60 | Derived. One row per learner, contrast and measure |
| `species-medians.tsv` | 240 | Derived. One row per learner, contrast, measure and species |
| `rule.tsv` | 5 | Derived. The plan's rule per learner |
| `arms.tsv` | 15 | Derived. One row per learner and arm, pooled over datasets |
| `mix.tsv` | 3 | Derived. One row per arm: where its training genomes come from |
| `summary.tsv` | 40 | Derived. Single figures the page quotes |

`SHA256SUMS` gives the sha256 of every other file, and this README is what the script writes.

## Columns

- `fits.tsv`: `dataset`, `species`, `arm` and `learner` name the fit. `test` and `test_resistant` count its test genomes and the resistant ones, and `train` and `train_resistant` its training genomes and the resistant ones. `tp`, `fn`, `tn` and `fp` count resistant genomes called resistant and called susceptible, and susceptible genomes called susceptible and called resistant, at a score of 0.5. `wrong_at_clip` counts wrong calls at the clip and `scores_0_or_1` scores of exactly 0 or 1. `loss_sum` is the sum of clipped log loss over the test genomes and `clip_loss_sum` the part of it from wrong calls at the clip. `null_fit` is 1 for a gbt or lasso fit of complexity 0.
- `datasets.tsv`: `dataset`, `species` and `antibiotic` name the dataset. `cohort` and `cohort_resistant` count the genomes `grm learn` evaluates for it and the resistant ones. `cut_day` is its cut day and `share_public_by_cut` the share of the cohort public by that day. `drawn_test_studies` counts the test studies drawn and `studies_with_a_tested_genome` those with a tested genome. `test`, `test_resistant`, `train` and `train_resistant` count the test and training genomes and the resistant ones. `past_from_past`, `past_from_later_other_studies` and `past_from_later_test_studies` count the past arm's training genomes that are past, later genomes of other studies and later genomes of the test studies. `random_from_past`, `random_from_later_other_studies`, `random_from_later_test_studies`, `study_out_from_past`, `study_out_from_later_other_studies` and `study_out_from_later_test_studies` count the same for the random and study-out arms.
- `species.tsv`: `species` names the species and `datasets` counts its datasets. `tested_pairs` and `tested_resistant_pairs` count its tested genome and dataset pairs and the resistant ones, `distinct_test_genomes` its distinct tested genomes, `distinct_test_studies` its distinct test studies, and `largest_study_tested_genomes` the tested genomes of its largest test study.
- `review.tsv`: `learner` names the learner. `fits` and `fits_changed` count its fits and those in which the review's fix changed a score, and `scores` and `scores_changed` its scores and those that changed. `readings` and `readings_changed` count its readings and those that changed, and `changed_readings` names each changed one by contrast and measure, with `;` between them.
- `totals.tsv`: `name` names a count over all datasets and `value` gives it.
- `order.tsv`: `event` names the event, `time` gives its UTC time and `source` the file it is read from.
- `settings.tsv`: `name` names a setting and `value` gives it: the design's `cut_share`, its class minimums `min_test_class` and `min_train_class`, its pool minimums `min_reserve` and `min_later_other_studies`, the fits' `seed`, `screened_features` and `inner_folds`, the settings the script applies, and the sha256 of the two plans and two scripts.
- `checks.tsv`: `check` names a check the build ran on the input set and `count` gives how many passed it. `genome_ids_of_the_inputs_in_the_record` gives the count of the input set's genome ids found in the record.
- `inputs.sha256`: Each line gives the sha256 of a file of the input set and its path there.
- `contrasts.tsv`: `dataset`, `species`, `learner` and `contrast` name the row. `d_recall`, `d_log_loss`, `d_specificity` and `d_balanced_accuracy` give the paired difference in each measure.
- `readings.tsv`: `learner`, `contrast` and `measure` name the reading, and `minuend` and `subtrahend` the two arms of the contrast. `median` is the median difference over datasets, `low` and `high` its interval from 10,000 species draws, and `reading` what that interval reads. `exact_low`, `exact_high` and `exact_reading` give the same over every ordered species draw. `datasets` counts the datasets, and `worse`, `better` and `equal` those whose difference is on each side of zero or at zero. `minuend_median` and `subtrahend_median` give each arm's median, and `plan` what the plan fixed.
- `species-medians.tsv`: `learner`, `contrast`, `measure` and `species` name the row. `datasets` counts the species' datasets and `median` gives their median difference.
- `rule.tsv`: `learner` names the learner. `recall` and `log_loss` give its total contrast readings, and `rule` what the plan's rule makes of them.
- `arms.tsv`: `learner` and `arm` name the row and `datasets` counts its datasets. `tested_pairs` and `resistant_pairs` count the pooled test pairs and the resistant ones, and `tp`, `fn`, `tn` and `fp` the pooled confusion counts. `recall` and `specificity` follow from them, and `wrong_calls` is `fn` plus `fp`. `wrong_at_clip`, `scores_0_or_1` and `null_fits` add up the fits' cells. `pooled_log_loss` is the summed `loss_sum` divided by the summed `test`, and `clip_share` the share of it from wrong calls at the clip.
- `mix.tsv`: `arm` names the arm and `datasets` counts the datasets. `train` counts its training genomes over all datasets, `from_past`, `from_later_other_studies` and `from_later_test_studies` how many are past, later genomes of other studies and later genomes of the test studies, and `past_share` the share that is past.
- `summary.tsv`: `name` names a single figure the page quotes and `value` gives it.

## Rules

**Measures.** A score at or above 0.5 calls resistant. Recall and specificity are fractions of the confusion counts and balanced accuracy is their mean. Log loss clips the score given to the true class to [0.0001, 1 - 0.0001], and a fit's log loss is `loss_sum` divided by `test`. A wrong call at the clip is a resistant genome scored at most 0.0001 or a susceptible genome scored at least 1 - 0.0001.

**Contrasts.** total is past minus random, study is study-out minus random, and time is past minus study-out, each a paired difference on one dataset's test genomes.

**Intervals.** `low` and `high` are the 251st and 9,750th of 10,000 medians, each over the datasets of the species drawn with replacement, from one generator seeded 20261114 and drawn in the order of `readings.tsv`, as the plan's analysis script draws them. `exact_low` and `exact_high` take the same two quantiles over every ordered draw of the species, each drawn as often as any other.

**Readings.** A median reads worse when its interval lies wholly on the worse side of zero, better when wholly on the other side, and not shown otherwise. Worse is lower recall, specificity or balanced accuracy, or a higher log loss. In `readings.tsv`, `plan` is `rule` for the two readings the rule reads, `fixed` for recall and log loss on the other two contrasts, and `supporting` for specificity and balanced accuracy, which the plan named as supporting figures.

## Limits

Each row is per dataset, species, arm, learner, contrast, measure or event, or a single count. No row is about one genome, and no cell names a genome or a study. 3 datasets have one drawn and one tested study (`drawn_test_studies` and `studies_with_a_tested_genome` are 1), so their rows in `fits.tsv`, `contrasts.tsv` and `datasets.tsv` are results on one unnamed study. In `species.tsv`, `largest_study_tested_genomes` is the count of the tested genomes of one unnamed study. The result files hold no study, so the study counts rest on the input set's study count file alone. `--check` alone cannot hold a row file to the result files; that needs the input set.

## Guards

`--check` with the input set rebuilds every file and fails on any byte that differs. `--check` alone reads the row files, derives every derived file and this README again, and fails on any byte that differs or on a file that differs from `SHA256SUMS`. The input set is not in the repository: it holds the result files, whose scores are per genome. `--out` writes only to a new path, an empty directory or an earlier record of this script whose files agree with its own `SHA256SUMS`, and refuses a symlink.

## Commands

```sh
python3 -I scripts/time-split/record.py --inputs <dir> --out reference/time-split
python3 -I scripts/time-split/record.py --inputs <dir> --check reference/time-split
python3 -I scripts/time-split/record.py --check reference/time-split
python3 -I -m unittest discover -s tests/time_split -t .
```
