# The cap-arms reading: the record

This directory holds the tables behind [docs/cap-arms-result.md](../../docs/cap-arms-result.md): plan block C2, `cap-all` (the top 4,096 genomes a species by the release's selection rule) against `cap-2000` (the release's 2,000), scored on the default's held-out genomes, beside a twin, with `cap-4000` withdrawn before any fit. The first run holds 143 datasets in 8 species, 39 of them with a study split. Every row is an aggregate over a dataset, a species, a split, an arm or a learner. No file holds a genome, a BioProject, a sequence or a feature.

`scripts/c2-cap/score.py --count` wrote the files under `count/`, copied here as they are; `scripts/c2-cap/page_tables.py write` checked them against the scorer's own count of the species' `datasets.tsv` files and wrote every other file here from the run's `score/` and `results/` tables, after the 14 gates of `gates.txt` passed. The rules of the reading were fixed before any C2 score was read (`docs/assumption-plan-corrections.md`, entries "C2 design", "C2 rulings" and "C2 review").

| File | Holds |
|---|---|
| `gates.txt` | the gates the program passed before it wrote a table |
| `checks.txt` | the printout of `scripts/c2-cap/run_checks.py`, the checks of the run's copies made before any score was read, counts only; gate 1 holds that each check passed |
| `count/count.tsv`, `count/conditions.tsv` | the count that decides, and the conditions as the scorer reads them |
| `page-reading.tsv` | `cap-all`'s deciding conditions under ruling 15, 1 and 2, 3 and 4, and whether all are met; 6 is the maintainer's |
| `page-conditions.tsv` | the conditions per arm and reading, with which of them decide: ruling 15's 1 and 2, 3 and 4 for `cap-all`; condition 5 reads the study split and decides nothing; the stricter reading decides nothing; `cap-4000`'s rows are withdrawn and not read; 6 is the maintainer's |
| `page-count.tsv` | every row of `count.tsv`, with conditions 1 and 2 read on the lineage primary rows; study rows and `cap-4000`'s rows marked |
| `page-count-without-few.tsv` | beside it, deciding nothing: the same count with `cap-all`'s rows of the datasets that add fewer than 20 labelled genomes left out, Holm over its own rows |
| `page-few-added.tsv` | the datasets that add fewer than 20 labelled genomes under `cap-all`, named in full, with the added genomes by class; they stay in the count |
| `page-species.tsv` | per species, named in full: the cohorts' sizes, the added genomes, and the dataset and learner pairs better, worse and the same; deciding nothing |
| `page-datasets.tsv` | per dataset, split and learner, species and antibiotic in full: `cap-2000`'s and `cap-all`'s balanced accuracy and recall, the differences, the twin's, the draws scored, the genomes and the added genomes by class, and the role |
| `page-not-fitted.tsv` | the datasets an arm does not score, named in full with the reason |
| `page-study.tsv` | per study dataset: the added genomes that may train under `cap-all` and those that trained; deciding nothing |
| `page-checks.tsv` | the scorer's checks per name, passed over made |
| `page-binaries.tsv` | each binary by commit and code sha256, with its species named, and under it each binary sha256 with its species |
| `page-times.tsv` | the fits' wall and CPU seconds and peak memory per arm and split, as sums and medians |

A dataset's figure is the mean of its draws: `grm arm` slots 0 to 2 for each arm and 3 to 5 for the twin, at evaluation seed 42; the study split has one fit per dataset and arm. An added genome only trains and is never scored, so every arm is read on the default's held-out genomes. A split by lineage or study is a retrospective split inside one database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c2-cap/page_tables.py print` prints the page's tables from these files.
