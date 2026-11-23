# The input-arm reading of plan block C3: the record

This directory holds the tables behind [docs/input-arms-result.md](../../docs/input-arms-result.md): plan block C3, five other inputs beside the presence of 31-mers (K31), each against K31 refitted in the same run, read once on the first run's 233 datasets in 23 species, 55 of them with a study split. The two program arms run on A's first run, 204 datasets in 21 species; the three matrix arms, the default and the twin on L's. Every row is an aggregate over a dataset, a species, a split, an arm or a learner. No file holds a genome, a k-mer, a unitig, a sequence, a position or a column name of a matrix.

`scripts/c3-inputs/score.py --count` wrote the three files under `count/`. They are score.py's own tables: `scripts/c3-inputs/page_tables.py write` counted them again from the species' `datasets.tsv` files (gate 13), copied them here as they are, and sums them in `SHA256SUMS`. `checks.txt` is the printout of `scripts/c3-inputs/run_checks.py` over the run's copies, made before any score was read; gate 1 holds that it says every check passed. It wrote every other file here from the run's own `score/` tables and each species' run files, after the 13 gates of `gates.txt` passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed before any C3 score was read (`docs/assumption-plan-corrections.md`, entries "C3 design" and "C3 arms").

One binary: every species' `binary.txt` gives the code sha256 `d3de72fcc10f6ef9503eb43f8cb5fd1449bca137d5f0a1594caf6fda0c614363`, 23 of 23 species.

| File | Holds |
|---|---|
| `gates.txt` | the gates the program passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the printout of `scripts/c3-inputs/run_checks.py`, the checks of the run's copies made before any score was read, counts only |
| `count/count.tsv` | the count that decides, per split, arm row, learner and measure |
| `count/count-no-constant.tsv` | the count with the constant-model datasets left out, which decides nothing |
| `count/conditions.tsv` | the plan's conditions per arm, read by ruling 15, which decides, and by the stricter reading, which decides nothing |
| `page-conditions.tsv` | `conditions.tsv` with a column `decides`: 1 for ruling 15's rows, 0 for the stricter reading's |
| `page-count.tsv` | every row of both counts, with `decides`, and conditions 1 and 2 read on the lineage primary rows |
| `page-arms.tsv` | per arm row: its set, its first-run datasets and species, its study datasets, whether it was withdrawn, its build and its study runs refused for no candidate |
| `page-species.tsv` | per species, named in full: datasets in L, A and S, datasets of M waiting, constant-model datasets and the binary's code sha256 |
| `page-species-count.tsv` | per species, split, arm row and learner: better, worse and same on the primary measure, and the median difference |
| `page-datasets.tsv` | per dataset, split, arm row and learner, species and antibiotic in full: the arm's balanced accuracy and recall, the default's, the differences, the draws scored, constant-model draws and a study run refused |
| `page-builds.tsv` | per species and matrix build, the counts its `build.json` gives, and nothing else |
| `page-binaries.tsv` | per species: the commit, binary and code sha256 of its run, restarts, and the outside programs' versions and sha256 |

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. A difference is always the arm's value less the default's. A dataset with a draw not scored is not scored. An arm's study run refused for no candidate is scored as the constant model on the default's folds (ruling 16), enters the count, and is marked in `constant`. "No difference shown" is the rule's reading and no statement that two inputs are equal.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c3-inputs/page_tables.py print` prints the page's tables from these files.
