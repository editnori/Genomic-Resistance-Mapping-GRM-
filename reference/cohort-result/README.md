# The small-cohort reading of 22 November 2026: the record

This directory holds the tables behind [docs/cohort-result.md](../../docs/cohort-result.md): plan block C1, part B, read on 89 datasets in 10 species. Every row is an aggregate over a dataset, a learner or a line. No file holds a genome, a sequence or a feature.

`scripts/c1-cutoff/read.py` wrote every file here but this README, `part-a-lines.tsv` and `SHA256SUMS`. It ran once, at commit `b030d4e`, from 20:49:31Z to 20:50:42Z on 22 Nov 2026, with `--measures-over shared --line-50 beside --twin-rule beside --dataset-differences yes`. The rules of the reading were fixed and public before it ran (`docs/assumption-plan-corrections.md`, entry "C1 reading").

| File | Holds |
|---|---|
| `gates.txt` | the eight gates the program passed before it wrote a table |
| `summary.txt` | the reading of each lower line and of the twin, in the fixed words |
| `page-deciding.tsv` | the rows that decide: lines 10, 15 and 20 against 25, per learner and per line |
| `page-twin.tsv` | the twin: 25 at fold seed 45 against 25 at seed 42 |
| `page-other.tsv` | every row that decides nothing, line 50 among them |
| `page-measures.tsv` | recall, specificity, both error rates and log loss, for the estimate and the truth |
| `page-entry.tsv` | each dataset's mean difference per line and learner, by species and antibiotic |
| `page-species.tsv`, `page-datasets.tsv` | the 10 species and the 89 datasets by name |
| `count/count.tsv` | the scorer's own count, as `score.py --count` wrote it |
| `part-a-lines.tsv` | part A's line counts of 21 Nov 2026, taken with no fit |

The scorer writes a `candidate` cell on the twin's row and on line 50's row as on any line. `count/count.tsv` keeps those cells. Neither the twin nor line 50 is ever a candidate, and the page's tables say so.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c1-cutoff/page_tables.py` prints the page's tables from these files, and `tests/c1_cutoff/test_page.py` checks the page against them.
