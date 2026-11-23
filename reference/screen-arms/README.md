# The screen and floor reading: the record

This directory holds the tables behind [docs/screen-arms-result.md](../../docs/screen-arms-result.md): plan block C4, the Cochran-Mantel-Haenszel screen over lineage strata at budgets 4,000 and 16,000 (`cmh-4000`, `cmh-16000`) and the candidate floor fixed at 3 (`floor-3`), each against the release's chi-square screen at 4,000 with the floor max(3, n/200), beside a twin, read once on the first run's 233 datasets in 23 species, 55 of them with a study split. Every row is an aggregate over a dataset, a species, an arm or a learner. No file holds a genome, a sequence or a feature.

`scripts/c4-screen/score.py --count` wrote the four tables under `count/`, and `withdrawn.txt` beside them names the arms it took as withdrawn, one a line. `checks.txt` is the printout of `scripts/c4-screen/run_checks.py`, the checks of the run's copies made before any score was read. `scripts/c4-screen/page_tables.py write` checked the reading against the 16 gates of `gates.txt`, and wrote every other file here once they passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed at 05:52Z on 23 Nov 2026, before any C4 score was read (`docs/assumption-plan-corrections.md`, entry "C4 design").

| File | Holds |
|---|---|
| `checks.txt` | the checks of the run's copies before the reading, counts only |
| `gates.txt` | the gates the table script passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/count.tsv`, `count/conditions.tsv` | the count that decides, and the plan's conditions per arm as `score.py` read them |
| `count/count-no-constant.tsv`, `count/count-floor-all.tsv` | the counts beside it, each deciding nothing |
| `page-conditions.tsv` | the conditions per arm, `decides` 1 for the reading that decides and 0 for the stricter reading beside it; condition 6 is the maintainer's |
| `page-count.tsv` | every row of the count that decides, with the 21 rows of the Holm family marked and conditions 1 and 2 read on each |
| `page-beside.tsv` | every row of the two counts beside it, each marked as deciding nothing |
| `page-species.tsv` | per species, named in full, split, pair, set and learner: better, worse and the same on the primary measure; it decides nothing |
| `page-datasets.tsv` | per dataset, species and antibiotic in full, split, pair and learner: the default's value and the arm's on balanced accuracy and recall, the differences, the draws scored, the genome count, whether the floor arm moves the dataset and whether the cell is counted, and a study run refused |
| `page-checks.tsv` | `score.py`'s checks of every result as checks passed of checks made, and the table script's own check that a screen arm's `scm` and `cart`, and the floor arm below 800 genomes, equal the default |
| `page-fit-times.tsv` | per arm and learner, the outer fits' wall and CPU seconds and peak memory, as sums and medians |
| `page-binaries.tsv` | each species' commit, binary and code sha256, and its restarts |

All 23 species ran one binary by its code sha256, `5a39538947ab78961ac5c07e79adfb46851489f1782a466941df5c50aace3c63`, built from 1 commit: `2ee9e6b`.

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. A screen arm is counted on `rf`, `gbt` and `lasso`, the learners its screen reaches; its `scm` and `cart` read no screen and are checked equal to the default's, not counted. The floor arm is counted on all five learners over the datasets of 800 genomes or more, the ones its floor can move. A split by lineage or study is a retrospective split inside one public database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c4-screen/page_tables.py print` prints the page's tables from these files.
