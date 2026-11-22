# The D0 and C8 reading of 22 November 2026: the record

This directory holds the tables behind [docs/refit-result.md](../../docs/refit-result.md): the assumption plan's diagnostic D0 and its block C8, read over a refit of the release's rule on 238 datasets in 23 species. Every row is an aggregate over a species, a dataset, a learner or a split. No file holds a genome, a sequence or a feature.

`scripts/assumption-d0-c8/read_refit.py` wrote every file here but this README and `SHA256SUMS`. It ran once, at commit `0f5cc5c`, at 21:40Z on 22 Nov 2026, with `--equal-tables datasets.tsv,learners.tsv --c8-one-fold zero --c8-over-five pairs --c8-bands none --order-rest default-after`. The rules of the reading were fixed and public before it ran (`docs/assumption-plan-corrections.md`, entries "D0 causes" and "Refit reading").

| File | Holds |
|---|---|
| `gates.txt` | the five gates the program passed before it wrote a table |
| `summary.txt` | the gates, the rules chosen, D0's ordering line with its coverage and stability words |
| `page-d0-order.tsv` | each cause's miss pairs of the 177, and the cause that sets the order |
| `page-d0-learners.tsv` | per learner: miss pairs, datasets in each cause, each cause's share of the gap, ties |
| `page-d0-datasets.tsv` | each of the 177 miss pairs by species and antibiotic: both balanced accuracies, the gap, the missed folds and the cause |
| `page-c8-summary.tsv` | per learner and split: the median overlap with its quartiles, the units counted, the features summed |
| `page-c8-datasets.tsv` | each dataset, learner and split: folds that chose, pairs, both overlaps and three feature counts |
| `page-species.tsv` | the 23 species with their antibiotics in full, D0's datasets, and the refit's check counts |
| `d0/d0-summary.tsv`, `d0/d0-ordering.txt` | D0's own sum over the species, as `d0.py combine` wrote it |

One change was made to two files. `summary.txt` and `gates.txt` named a local path, the label count's table given on the command line. Here that path reads "the label count's datasets.tsv". No other byte of a file differs from the program's output.

The program also wrote `c8/c8-combined-summary.tsv` and `c8/c8-combined-datasets.tsv`. They are equal in bytes to `page-c8-summary.tsv` and `page-c8-datasets.tsv` and are not kept twice.

`d0/d0-summary.tsv` holds the gap share pooled over the five learners. By the reading's last ruling it is not on the page.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/assumption-d0-c8/page_tables.py` prints the page's tables from these files, and `tests/assumption_d0_c8/test_page.py` checks the page against them.
