# The outside-learner reading of plan block C7: the record

This directory holds the tables behind [docs/outside-learners-result.md](../../docs/outside-learners-result.md): plan block C7, XGBoost 3.0.5 and LightGBM 4.6.0 on the 4,000 rows `grm export` screens, each read against GRM's boosted trees (`gbt`), and a scikit-learn 1.7.2 elastic net on unitigs read against GRM's lasso (`lasso`), beside a twin, read once on the first run's 233 datasets in 23 species, 55 of them with a study split. Every row is an aggregate over a dataset, a species, a split, an arm, a pair or a fold count. No file holds a genome, a k-mer, a unitig, a sequence, a position, a feature name, an elastic net weight or a booster's feature importance.

`scripts/c7-outside/score.py --count` wrote the four files under `count/`. They are score.py's own tables: `scripts/c7-outside/page_tables.py write` counted them again from the species' `datasets.tsv` files (gate 17), copied them here as they are, and sums them in `SHA256SUMS`. It wrote every other file here from the run's own `score/` tables and each species' run files, after the 17 gates of `gates.txt` passed. `checks.txt` is the printout of `scripts/c7-outside/run_checks.py check`. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed before any C7 score was read (`docs/assumption-plan-corrections.md`, entry "C7 design").

No arm was withdrawn: all three arms ran.

One binary: every species' `DONE` gives the code sha256 `d3de72fcc10f6ef9503eb43f8cb5fd1449bca137d5f0a1594caf6fda0c614363` of the C3 run's `binary.txt`, 23 of 23 species.

One version of each library in every species: XGBoost 3.0.5, LightGBM 4.6.0, scikit-learn 1.7.2, numpy 1.26.4 and scipy 1.15.3, on python 3.12.3.

| File | Holds |
|---|---|
| `gates.txt` | the gates the program passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the checks of the run, made before any score was read |
| `count/count.tsv` | the count that decides, per split, pair, learner and measure |
| `count/count-no-constant.tsv` | the count with the datasets that hold a constant fold left out, which decides nothing |
| `count/count-no-unconverged.tsv` | the count with the datasets whose elastic net draws hold a refit that reached 2,000 epochs left out, which decides nothing |
| `count/conditions.tsv` | the plan's conditions per arm; the arm against its own learner's twin pair alone stands beside with `decides` 0 |
| `page-conditions.tsv` | `conditions.tsv` per arm, with its pair and whether the arm was withdrawn |
| `page-count.tsv` | every row of the three counts, with `decides` 1 for `count.tsv` alone, the pair as the page names it and, in `score_pair`, as `score.py` names it, the row it is read against, and conditions 1 and 2 on the deciding lineage rows of the primary measure |
| `page-withdrawn.tsv` | each arm withdrawn, with the reason and the time the lead gave, and how many species had run a second pass with it |
| `page-species.tsv` | per species, named in full, split and pair: datasets better, worse and the same of those entered, which decides nothing |
| `page-datasets.tsv` | per dataset, split and pair, species and antibiotic in full: the arm's balanced accuracy and recall, the default learner's at the arm's seeds, the differences, the draws scored and the dataset's genomes |
| `page-folds.tsv` | per arm and split: the folds of its reports, and of them the boosters' one-tree folds and the elastic net's unconverged refits, null fits and constant folds |
| `page-libraries.tsv` | the library versions each species' run read back: one row when they agree, a row per species when they do not |
| `page-binaries.tsv` | per species: the binary, code and export binary sha256 its `DONE` names, the elastic net's pass, and its result count |
| `page-checks.tsv` | `score.py`'s checks per name and `plans.py`'s, passed of made |
| `page-times.tsv` | per arm: its reports, the fits timed, and the sums and medians of their wall and CPU seconds |

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, and one study run at seed 42. Each arm is read against the default's own learner at the arm's seeds (`k31-at-arm-seeds`); the twin, the default at seeds 45, 46 and 47, against the default. A difference is always the arm's value less the default's. A dataset with a draw not scored is not scored. "No difference shown" is the rule's reading and no statement that two learners are equal.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c7-outside/page_tables.py print` prints the page's tables from these files.
