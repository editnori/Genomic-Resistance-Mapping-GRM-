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

Each line of `inputs.sha256` gives the sha256 of a file of the input set and its path there. The other files' columns, in file order:

| File | Column | Gives |
|---|---|---|
| `fits.tsv` | `dataset` | the dataset, as species and antibiotic |
| `fits.tsv` | `species` | the species |
| `fits.tsv` | `arm` | past, random or study-out |
| `fits.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `fits.tsv` | `test` | test genomes |
| `fits.tsv` | `test_resistant` | resistant test genomes |
| `fits.tsv` | `train` | training genomes |
| `fits.tsv` | `train_resistant` | resistant training genomes |
| `fits.tsv` | `tp` | resistant genomes called resistant at a score of 0.5 |
| `fits.tsv` | `fn` | resistant genomes called susceptible |
| `fits.tsv` | `tn` | susceptible genomes called susceptible |
| `fits.tsv` | `fp` | susceptible genomes called resistant |
| `fits.tsv` | `wrong_at_clip` | wrong calls at the clip |
| `fits.tsv` | `scores_0_or_1` | scores of exactly 0 or 1 |
| `fits.tsv` | `loss_sum` | the sum of clipped log loss over the test genomes |
| `fits.tsv` | `clip_loss_sum` | the part of `loss_sum` from wrong calls at the clip |
| `fits.tsv` | `null_fit` | 1 for a gbt or lasso fit of complexity 0 |
| `datasets.tsv` | `dataset` | the dataset, as species and antibiotic |
| `datasets.tsv` | `species` | the species |
| `datasets.tsv` | `antibiotic` | the antibiotic |
| `datasets.tsv` | `cohort` | genomes `grm learn` evaluates for the dataset |
| `datasets.tsv` | `cohort_resistant` | resistant genomes of the cohort |
| `datasets.tsv` | `cut_day` | the cut day |
| `datasets.tsv` | `share_public_by_cut` | the share of the cohort public by the cut day |
| `datasets.tsv` | `drawn_test_studies` | test studies drawn |
| `datasets.tsv` | `studies_with_a_tested_genome` | test studies with a tested genome |
| `datasets.tsv` | `test` | test genomes |
| `datasets.tsv` | `test_resistant` | resistant test genomes |
| `datasets.tsv` | `train` | training genomes |
| `datasets.tsv` | `train_resistant` | resistant training genomes |
| `datasets.tsv` | `past_from_past` | past genomes in the past arm's training |
| `datasets.tsv` | `past_from_later_other_studies` | later genomes of other studies in the past arm's training |
| `datasets.tsv` | `past_from_later_test_studies` | later genomes of the test studies in the past arm's training |
| `datasets.tsv` | `random_from_past` | past genomes in the random arm's training |
| `datasets.tsv` | `random_from_later_other_studies` | later genomes of other studies in the random arm's training |
| `datasets.tsv` | `random_from_later_test_studies` | later genomes of the test studies in the random arm's training |
| `datasets.tsv` | `study_out_from_past` | past genomes in the study-out arm's training |
| `datasets.tsv` | `study_out_from_later_other_studies` | later genomes of other studies in the study-out arm's training |
| `datasets.tsv` | `study_out_from_later_test_studies` | later genomes of the test studies in the study-out arm's training |
| `species.tsv` | `species` | the species |
| `species.tsv` | `datasets` | datasets of the species |
| `species.tsv` | `tested_pairs` | tested genome and dataset pairs |
| `species.tsv` | `tested_resistant_pairs` | resistant tested pairs |
| `species.tsv` | `distinct_test_genomes` | distinct tested genomes |
| `species.tsv` | `distinct_test_studies` | distinct test studies |
| `species.tsv` | `largest_study_tested_genomes` | tested genomes of the largest test study |
| `review.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `review.tsv` | `fits` | fits |
| `review.tsv` | `fits_changed` | fits in which the review's fix changed a score |
| `review.tsv` | `scores` | scores |
| `review.tsv` | `scores_changed` | scores the fix changed |
| `review.tsv` | `readings` | readings |
| `review.tsv` | `readings_changed` | readings the fix changed |
| `review.tsv` | `changed_readings` | each changed reading as contrast and measure, with `;` between them |
| `totals.tsv` | `name` | the name of a count over all datasets |
| `totals.tsv` | `value` | the figure |
| `order.tsv` | `event` | the event |
| `order.tsv` | `time` | its UTC time |
| `order.tsv` | `source` | the file the time is read from |
| `settings.tsv` | `name` | the name of a setting |
| `settings.tsv` | `value` | the setting |
| `checks.tsv` | `check` | a check the build ran on the input set |
| `checks.tsv` | `count` | how many passed it |
| `contrasts.tsv` | `dataset` | the dataset, as species and antibiotic |
| `contrasts.tsv` | `species` | the species |
| `contrasts.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `contrasts.tsv` | `contrast` | total, study or time |
| `contrasts.tsv` | `d_recall` | the paired difference in recall |
| `contrasts.tsv` | `d_log_loss` | the paired difference in log loss |
| `contrasts.tsv` | `d_specificity` | the paired difference in specificity |
| `contrasts.tsv` | `d_balanced_accuracy` | the paired difference in balanced accuracy |
| `readings.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `readings.tsv` | `contrast` | total, study or time |
| `readings.tsv` | `minuend` | the arm the contrast takes from |
| `readings.tsv` | `subtrahend` | the arm the contrast subtracts |
| `readings.tsv` | `measure` | the measure read |
| `readings.tsv` | `median` | the median difference over datasets |
| `readings.tsv` | `low` | the 251st of 10,000 sorted species draws |
| `readings.tsv` | `high` | the 9,750th of 10,000 sorted species draws |
| `readings.tsv` | `reading` | worse, better or not shown, from `low` and `high` |
| `readings.tsv` | `exact_low` | the same quantile as `low` over every ordered species draw |
| `readings.tsv` | `exact_high` | the same quantile as `high` over every ordered species draw |
| `readings.tsv` | `exact_reading` | the reading from `exact_low` and `exact_high` |
| `readings.tsv` | `datasets` | datasets the median is over |
| `readings.tsv` | `worse` | datasets whose difference is on the worse side of zero |
| `readings.tsv` | `better` | datasets whose difference is on the other side |
| `readings.tsv` | `equal` | datasets whose difference is zero |
| `readings.tsv` | `minuend_median` | the median of the measure in the first arm |
| `readings.tsv` | `subtrahend_median` | the median of the measure in the second arm |
| `readings.tsv` | `plan` | rule, fixed or supporting: what the plan fixed for this reading |
| `species-medians.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `species-medians.tsv` | `contrast` | total, study or time |
| `species-medians.tsv` | `measure` | the measure read |
| `species-medians.tsv` | `species` | the species |
| `species-medians.tsv` | `datasets` | datasets of the species |
| `species-medians.tsv` | `median` | the median difference over the species' datasets |
| `rule.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `rule.tsv` | `recall` | the total contrast's recall reading |
| `rule.tsv` | `log_loss` | the total contrast's log loss reading |
| `rule.tsv` | `rule` | cost, benefit, mixed or not shown, by the plan's rule |
| `arms.tsv` | `learner` | scm, cart, rf, gbt or lasso |
| `arms.tsv` | `arm` | past, random or study-out |
| `arms.tsv` | `datasets` | datasets pooled |
| `arms.tsv` | `tested_pairs` | tested genome and dataset pairs |
| `arms.tsv` | `resistant_pairs` | resistant tested pairs |
| `arms.tsv` | `tp` | resistant genomes called resistant at a score of 0.5 |
| `arms.tsv` | `fn` | resistant genomes called susceptible |
| `arms.tsv` | `tn` | susceptible genomes called susceptible |
| `arms.tsv` | `fp` | susceptible genomes called resistant |
| `arms.tsv` | `recall` | `tp` over `tp` plus `fn` |
| `arms.tsv` | `specificity` | `tn` over `tn` plus `fp` |
| `arms.tsv` | `wrong_calls` | `fn` plus `fp` |
| `arms.tsv` | `wrong_at_clip` | wrong calls at the clip |
| `arms.tsv` | `scores_0_or_1` | scores of exactly 0 or 1 |
| `arms.tsv` | `pooled_log_loss` | the summed `loss_sum` over the summed `test` |
| `arms.tsv` | `clip_share` | the share of the pooled log loss from wrong calls at the clip |
| `arms.tsv` | `null_fits` | fits with `null_fit` 1 |
| `mix.tsv` | `arm` | past, random or study-out |
| `mix.tsv` | `datasets` | datasets pooled |
| `mix.tsv` | `train` | training genomes, all datasets |
| `mix.tsv` | `from_past` | past genomes in training |
| `mix.tsv` | `from_later_other_studies` | later genomes of other studies in training |
| `mix.tsv` | `from_later_test_studies` | later genomes of the test studies in training |
| `mix.tsv` | `past_share` | the share of training genomes that are past |
| `summary.tsv` | `name` | the name of a figure the page quotes |
| `summary.tsv` | `value` | the figure |

## Rules

**Measures.** A score at or above 0.5 calls resistant. Recall and specificity are fractions of the confusion counts and balanced accuracy is their mean. Log loss clips the score given to the true class to [0.0001, 1 - 0.0001], and a fit's log loss is `loss_sum` divided by `test`. A wrong call at the clip is a resistant genome scored at most 0.0001 or a susceptible genome scored at least 1 - 0.0001.

**Contrasts.** total is past minus random, study is study-out minus random, and time is past minus study-out, each a paired difference on one dataset's test genomes.

**Intervals.** `low` and `high` are the 251st and 9,750th of 10,000 medians, each over the datasets of the species drawn with replacement, from one generator seeded 20261114 and drawn in the order of `readings.tsv`, as the plan's analysis script draws them. `exact_low` and `exact_high` take the same two quantiles over every ordered draw of the species, each drawn as often as any other.

**Readings.** A median reads worse when its interval lies wholly on the worse side of zero, better when wholly on the other side, and not shown otherwise. Worse is lower recall, specificity or balanced accuracy, or a higher log loss. In `readings.tsv`, `plan` is `rule` for the two readings the rule reads, `fixed` for recall and log loss on the other two contrasts, and `supporting` for specificity and balanced accuracy, which the plan named as supporting figures.

## Limits

Each row is per dataset, species, arm, learner, contrast, measure or event, or a single count. No row is about one genome, and no cell names a genome or a study. 3 datasets have one drawn and one tested study (`drawn_test_studies` and `studies_with_a_tested_genome` are 1), so their rows in `fits.tsv`, `contrasts.tsv` and `datasets.tsv` are results on one unnamed study. In `species.tsv`, `largest_study_tested_genomes` is the count of the tested genomes of one unnamed study.

The result files hold no study, so the study counts rest on the input set's study count file alone. `--check` alone cannot hold a row file to the result files; that needs the input set.

## Guards

`--check` with the input set rebuilds every file and fails on any byte that differs. `--check` alone reads the row files, derives every derived file and this README again, and fails on any byte that differs or on a file that differs from `SHA256SUMS`. The input set is not in the repository: it holds the result files, whose scores are per genome.

`--out` writes only to a new path, an empty directory or an earlier record of this script whose files agree with its own `SHA256SUMS`, and refuses a symlink.

## Commands

```sh
python3 -I scripts/time-split/record.py --inputs <dir> --out reference/time-split
python3 -I scripts/time-split/record.py --inputs <dir> --check reference/time-split
python3 -I scripts/time-split/record.py --check reference/time-split
python3 -I -m unittest discover -s tests/time_split -t .
```
