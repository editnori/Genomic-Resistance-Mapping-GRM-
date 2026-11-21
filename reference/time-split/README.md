# Time split record

This directory holds the record of the split by first-public day inside BV-BRC, rerun after the review. The page that reads it is [`docs/time-split.md`](../../docs/time-split.md). `scripts/time-split/record.py` writes every file here from the input set, and its `--check` mode fails on any byte it would not write.

## Files

| File | Rows | Contents |
|---|---|---|
| `fits.tsv` | 360 | Row file. One row per dataset, arm and learner: test genomes and resistant ones, training genomes and resistant ones, the confusion counts at a score of 0.5, wrong calls at the clip, scores of exactly 0 or 1, the sum of clipped log loss over the test genomes and the part of it from wrong calls at the clip, and a null fit flag |
| `datasets.tsv` | 24 | Row file. One row per dataset: species, antibiotic, cohort, cut day, share public by the cut day, test studies drawn and test studies with a tested genome, test and training counts, and where each arm's training genomes come from |
| `species.tsv` | 4 | Row file. One row per species: datasets, tested pairs, distinct tested genomes, distinct test studies, and the tested genomes of its largest test study |
| `review.tsv` | 5 | Row file. One row per learner: fits and scores, and how many the review's fix changed; readings, how many changed and which |
| `totals.tsv` | 7 | Row file. Counts over all datasets: distinct tested genomes and studies, studies tested in two species, result files and datasets with a changed fit |
| `order.tsv` | 7 | Row file. The times of the plans' hashes and of the first and last result of each run, with their source |
| `settings.tsv` | 22 | Row file. The design's cut share and class minimums, the fits' settings, the settings the script applies, and the sha256 of the two plans and two scripts |
| `checks.tsv` | 10 | The checks the build ran on the input set, with the count that passed each, and the count of the input set's genome ids found in the record |
| `inputs.sha256` | 155 | The sha256 of each file of the input set |
| `contrasts.tsv` | 360 | Derived. One row per dataset, learner and contrast: the paired difference in recall, log loss, specificity and balanced accuracy |
| `readings.tsv` | 60 | Derived. One row per learner, contrast and measure: median over datasets, the interval from 10,000 species draws and its reading, the same interval over every ordered species draw and its reading, datasets worse, better and equal, each arm's median, and what the plan fixed |
| `species-medians.tsv` | 240 | Derived. The median difference per species, for each learner, contrast and measure |
| `rule.tsv` | 5 | Derived. The plan's rule per learner, from the total contrast's recall and log loss readings |
| `arms.tsv` | 15 | Derived. One row per learner and arm, pooled over datasets: confusion counts, recall, specificity, wrong calls, wrong calls at the clip, scores of 0 or 1, pooled log loss and the share of it from wrong calls at the clip, and null fits |
| `mix.tsv` | 3 | Derived. One row per arm: training genomes over all datasets, and how many are past, later genomes of other studies and later genomes of the test studies |
| `summary.tsv` | 40 | Derived. Single figures the page quotes |

`SHA256SUMS` gives the sha256 of every other file, and this README is what the script writes.

## Rules

**Measures.** A score at or above 0.5 calls resistant. Recall and specificity are fractions of the confusion counts and balanced accuracy is their mean. Log loss clips the score given to the true class to [0.0001, 1 - 0.0001] and is `loss_sum` over the test genomes. A wrong call at the clip is a resistant genome scored at most 0.0001 or a susceptible genome scored at least 1 - 0.0001.

**Contrasts.** total is past minus random, study is study-out minus random, and time is past minus study-out, each a paired difference on one dataset's test genomes.

**Intervals.** `low` and `high` are the 251st and 9,750th of 10,000 medians, each over the datasets of the species drawn with replacement, from one generator seeded 20261114 and drawn in the order of `readings.tsv`, as the plan's analysis script draws them. `exact_low` and `exact_high` take the same two quantiles over every ordered draw of the species, each drawn as often as any other.

**Readings.** A median reads worse when its interval lies wholly on the worse side of zero, better when wholly on the other side, and not shown otherwise. Worse is lower recall, specificity or balanced accuracy, or a higher log loss. In `readings.tsv`, `plan` is `rule` for the two readings the rule reads, `fixed` for recall and log loss on the other two contrasts, and `supporting` for specificity and balanced accuracy, which the plan named as supporting figures.

## Guards

`--check` with the input set rebuilds every file and fails on any byte that differs. `--check` alone reads the row files, derives every derived file and this README again, and fails on any byte that differs or on a file that differs from `SHA256SUMS`. The input set is not in the repository: it holds the result files, whose scores are per genome.

## Commands

```sh
python3 -I scripts/time-split/record.py --inputs <dir> --out reference/time-split
python3 -I scripts/time-split/record.py --inputs <dir> --check reference/time-split
python3 -I scripts/time-split/record.py --check reference/time-split
python3 -I -m unittest discover -s tests/time_split -t .
```
