# Kover comparison record

This directory is the record behind [`docs/kover-comparison.md`](../../docs/kover-comparison.md): Kover 2.0.7 against the released lineage calls of GRM's `scm` learner on 12 public datasets of release `rev14-274-20261112`. `scripts/kover-comparison/record.py` writes every file of it, and its `--check` mode fails on any byte it would not write.

It is not `reference/kover-summary.json`. That file holds the 73 historical Kover rows of the legacy application, on other cohorts and random splits, and `grm report` reads it. This record is read by no program.

## Files

| File | Rows | Contents |
|---|---|---|
| `README.md` |  | this list |
| `SHA256SUMS` |  | the sha256 of every other file of this directory |
| `cohort.tsv` | 8 | the 12 datasets against all 268 eligible datasets of the release: sizes, class balance, lineages and GRM's median lineage balanced accuracy |
| `datasets.tsv` | 12 | the 12 datasets as the plan chose them, byte for byte the file its hash list names: genomes, resistant, susceptible, lineages and test genomes per lineage fold |
| `folds.tsv` | 60 | one row per dataset and lineage fold: test genomes, the confusion counts and rule count of each side, Kover's kind and p, and the calls on which the two agree |
| `mapping.tsv` | 12 | the plan's table of how GRM's scm settings map to Kover's, with whether each is the same |
| `order.tsv` | 8 | the plan's hash list with its time, the first Kover fit, the review's plan and the review's tie summary |
| `per_dataset.tsv` | 12 | derived: per dataset, the mean over its folds of Kover minus GRM in each measure, a t interval for balanced accuracy, and the share of calls agreeing |
| `settings.tsv` | 19 | the release, the GRM learner, the Kover version and commands, the clip and the bootstrap |
| `steps.tsv` | 19 | the fits on all genomes with GRM's final kind and p held for Kover: covered and lost genomes at each greedy step, and whether the two rules' patterns over the genomes are equal |
| `summary.tsv` | 10 | derived: per measure, the mean over datasets of Kover minus GRM with its 95% bootstrap interval and seed, and the datasets and folds where Kover is higher, equal and lower |
| `ties.tsv` | 41 | exploratory: Kover's one-rule folds, the rules tied at the first step, and the test balanced accuracy each tied rule would give, as counts of rules per value |
| `totals.tsv` | 26 | derived: counts the page quotes |

Every file is plain text, and each TSV has a header row. All differences are Kover minus GRM. A file holds counts and measures per dataset, per fold, per greedy step or per measure: no genome, no sequence, no rule and no importance value.

## Columns

- `cohort.tsv`: `measure`, `chosen`, `eligible`.
- `datasets.tsv`: `order`, `dataset`, `species`, `genomes`, `resistant`, `susceptible`, `lineages`, `lineage_folds`, `genomes_listed_in_data_root`, `fold_sizes`.
- `folds.tsv`: `dataset`, `fold`, `test_genomes`, `resistant`, `susceptible`, `grm_tp`, `grm_fn`, `grm_tn`, `grm_fp`, `grm_rules`, `kover_tp`, `kover_fn`, `kover_tn`, `kover_fp`, `kover_rules`, `kover_kind`, `kover_p`, `agreeing`.
- `mapping.tsv`: `setting`, `grm`, `kover`, `same`.
- `order.tsv`: `utc`, `event`, `file`, `sha256`.
- `per_dataset.tsv`: `dataset`, `genomes`, `folds`, `balanced_accuracy`, `balanced_accuracy_t_low`, `balanced_accuracy_t_high`, `recall`, `specificity`, `log_loss`, `error_rate`, `rules`, `agreeing_share`.
- `settings.tsv`: `key`, `value`.
- `steps.tsv`: `dataset`, `step`, `kind`, `p`, `grm_covered`, `grm_lost`, `kover_covered`, `kover_lost`, `grm_listed_tied`, `counts_equal`, `pattern_equal`, `kover_equals_a_listed_tied`.
- `summary.tsv`: `measure`, `reading`, `value`, `low`, `high`, `seed`, `datasets_higher`, `datasets_equal`, `datasets_lower`, `folds_higher`, `folds_equal`, `folds_lower`.
- `ties.tsv`: `dataset`, `fold`, `rules_tied`, `rules_kept`, `patterns`, `first_index_equals_run`, `first_index_test_balanced_accuracy`, `tied_test_balanced_accuracy`.
- `totals.tsv`: `key`, `value`.

## Rules

- **Measures.** From each fold's confusion counts: balanced accuracy is the mean of recall and specificity, and log loss takes each 0 or 1 call as a probability clipped to [0.0001, 0.9999], so it is 0.0001 plus 9.2102 times the error rate.
- **Over datasets.** Each dataset's value is the mean of its 5 fold differences. The figure over datasets is their mean, with a 95% percentile interval from 10,000 resamples of whole datasets at the seed `summary.tsv` gives, drawn as the comparison drew them.
- **Readings.** `fixed` in `summary.tsv` is a measure the plan fixed before any Kover fit. `exploratory` is one chosen later: the error rate, and the three tie readings, which the review chose after reading the results per fold.
- **Ties.** In each of Kover's one-rule folds, `ties_equal` counts every tied first rule as equally likely, and `ties_lowest` and `ties_highest` take the lowest and highest test result among them, bounds chosen by test result. GRM's side is the same in all three.
- **Steps.** `pattern_equal` is an equality of hashes of the two rules' patterns over the genomes, as the comparison computed them; the record holds none of those hashes.

## Commands

```sh
python3 -I scripts/kover-comparison/record.py --check reference/kover-comparison
python3 -I scripts/kover-comparison/record.py --source "$KOVER" --review "$REVIEW" --check reference/kover-comparison
```

The first rebuilds the derived files, this README and `SHA256SUMS` from the row files. The second rebuilds every file from the comparison's and the review's output directories, which are not public.
