# The served-learner reading of 23 November 2026: the record

This directory holds the tables behind [docs/served-pick-result.md](../../docs/served-pick-result.md): plan block C9, the learner picked inside each outer training fold against the release's pick, read once on the first run's 233 datasets in 23 species, 55 of them with a study split. Every row is an aggregate over a dataset, a species, a split or a learner. No file holds a genome, a sequence or a feature.

`scripts/c9-pick/score.py` wrote the files under `count/`; `checks.txt` is the printout of the run's six checks, made before any score was read, and of the genome-id count over the reading and these tables, added beside the program's files. `scripts/c9-pick/page_tables.py write` wrote every other file here from the run's own `score/` folders and the reading's three scorings beside them, after the eight gates of `gates.txt` passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed at 20:36Z on 22 Nov 2026, before any C9 score was read (`docs/assumption-plan-corrections.md`, entry "C9 design").

| File | Holds |
|---|---|
| `checks.txt` | the run's six checks before the reading, and the genome-id count after it, counts only |
| `gates.txt` | the eight gates the program passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/read/` | the count that decides: `--tie later --no-pick default --default per-draw` |
| `count/tie-earlier/`, `count/tie-default/`, `count/default-seed42/` | the counts beside it, each deciding nothing |
| `page-conditions.tsv` | the plan's change conditions 1 to 5 under each answer, as the scorer read them; 6 is the maintainer's |
| `page-count.tsv` | every count row of the four answers, with conditions 1 and 2 read on the lineage primary rows |
| `page-species.tsv` | per species, named in full: datasets, and better, worse and same on the primary measure per split and pair |
| `page-datasets.tsv` | per dataset and split, species and antibiotic in full: the default's learner per draw and figure, the inner pick's figure, the difference, the twin's, the draws scored and the folds |
| `page-picks.tsv` | per answer and split: the folds each learner serves under the default and under the inner pick, folds with no pick, tied and differing from the default's learner |
| `page-pick-pairs.tsv` | the read answer, per split and default's learner: the folds the inner pick serves to each learner |
| `page-pooled.tsv` | the read answer, per split and arm: the calls summed over datasets and draws, with Rules' measures and Wilson intervals |

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. The pooled calls count each genome once a draw. The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism and the inner pick's cannot. A same or a worse is never read as showing that the release's pick is right.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c9-pick/page_tables.py print` prints the page's tables from these files, and `tests/c9_pick/test_page.py` checks the page against them.
