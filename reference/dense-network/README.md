# Dense network record

This directory is the record behind [`docs/dense-network.md`](../../docs/dense-network.md): one dense network fitted on 24 public datasets of release `rev14-274-20261112` and compared with the saved calls of GRM's five learners. `scripts/dense-network/record.py` writes every file of it, and its `--check` mode fails on any byte it would not write. No program reads it.

## Files

| File | Rows | Contents |
|---|---|---|
| `README.md` |  | this list |
| `SHA256SUMS` |  | the sha256 of every other file of this directory |
| `cost.tsv` | 4 | the network's fits and processes, its CPU seconds (user plus system) and the peak memory of a process |
| `datasets.tsv` | 24 | the 24 datasets in the plan's order: species, group, genomes, resistant genomes, lineage gap and whether the export passed gates A and B |
| `measures.tsv` | 520 | one row per dataset, split and learner: balanced accuracy, recall and log loss over the dataset's genomes; best5 and site_best name the learner they pick |
| `order.tsv` | 5 | the plan's file time and sha256, the fitting script's last change, the first network result, the analysis run on the first datasets, and the analysis script's last change |
| `rule.tsv` | 7 | derived: the four parts of the plan's rule, the value each takes and whether it holds |
| `settings.tsv` | 21 | the network, its inputs and fitting, the datasets, the comparator, the measures, the interval and the rule, as the plan fixed them |
| `summary.tsv` | 56 | derived: per split, comparator and measure, the medians over datasets and the median of the network minus the comparator with its interval over species |
| `totals.tsv` | 19 | derived: counts the page quotes |

Every file is plain text, and each TSV has a header row. A file holds measures per dataset, split and learner, or medians and counts over datasets: no genome, no sequence, no feature and no fitted parameter.

## Columns

- `cost.tsv`: `key`, `value`.
- `datasets.tsv`: `order`, `dataset`, `species`, `group`, `genomes`, `resistant`, `lineage_gap`, `gate_a`, `gate_b`.
- `measures.tsv`: `dataset`, `split`, `learner`, `picks`, `balanced_accuracy`, `recall`, `log_loss`.
- `order.tsv`: `utc`, `event`, `sha256`.
- `rule.tsv`: `rule`, `statistic`, `value`, `needs`, `holds`.
- `settings.tsv`: `key`, `value`.
- `summary.tsv`: `split`, `comparator`, `measure`, `candidate_median`, `comparator_median`, `difference_median`, `low`, `high`, `datasets`, `species`, `positive`, `negative`.
- `totals.tsv`: `key`, `value`.

## Rules

- **Network and comparator.** `candidate` is the network, the mean of its 3 repeats. `best5` is the best of the five in hindsight on that split, the plan's comparator on the lineage split.
- **Other comparators.** `site_best` is the learner the release names best, chosen on the random split; `mean4` averages the scores of `cart`, `rf`, `gbt` and `lasso`; `baseline` calls by lineage identity and has balanced accuracy only; `catalog` is the catalog rule, where the release scores it.
- **Over datasets.** `summary.tsv` gives the median over datasets of the network minus the comparator, with the 2.5% and 97.5% percentiles of that median over 2,000 resamples of species at seed 20261115, drawn as the experiment drew them.
- **Draw.** Each resample takes every dataset of each drawn species.
- **Log loss.** It is clipped at 0.0001; a positive log loss difference means the network's is higher.
- **Readings.** The plan fixed the rule on the lineage split against `best5`. Every other comparator, the random split, recall and log loss are reported beside the rule and do not enter it.
- **Exploratory.** `below_best5_lineage_*` in `totals.tsv` and the counts of `site_best` were read by the review after the results, so they are exploratory.

## Commands

```sh
python3 -I scripts/dense-network/record.py --check reference/dense-network
python3 -I scripts/dense-network/record.py --source "$RUN" --review "$REVIEW" --check reference/dense-network
```

The first rebuilds the derived files, this README and `SHA256SUMS` from the row files. The second rebuilds every file from the experiment's and the review's output directories, which are not public, and fails if a derived figure parts from the experiment's own by more than 1e-12.
