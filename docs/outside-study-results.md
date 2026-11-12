# Outside-study test: results

This is the result of a retrospective test on newly deposited studies from the same database. The plan is [outside-study-test.md](outside-study-test.md). It was frozen at commit `e3e1c34` before any of these genomes was fetched. The test is not a lockbox, and it is not external validation. The plan says why.

## The question

Each dataset report gives a learner's sensitivity and specificity from cross-validation on the cohort. The lineage split is the stricter of the two splits, because it holds out whole lineages. The test asks whether the lineage-split figures hold when the report's final model is applied, unchanged, to genomes from studies the cohort does not hold.

## Result

On genomes from studies the cohort does not hold, the figures hold on fewer than half of the datasets. On new genomes from studies the cohort already holds, they hold on most.

The primary population is 7,258 genomes on 82 datasets in 7 species. Every genome was scored and none was excluded. One dataset has no lineage-split figure in its report, so 81 datasets can be read.

A reading is "below" when the report's lineage-split figure is above the whole 95% interval of the new estimate. The plan reads sensitivity and specificity one at a time. This document also joins them: a learner holds on a dataset when neither its sensitivity nor its specificity reads below.

- The set covering machine holds on 32 of the 81 datasets and reads below on the other 49.
- The decision tree holds on 23 of the 81 datasets.
- The lasso can be scored and read on 25 of the 81 datasets, and holds on 7 of them.

Three facts belong beside those counts.

**The figures move both ways.** For the set covering machine, sensitivity reads above the report's figure on 36 of the 81 datasets and below it on 27. Its specificity reads below on 35 datasets and above on 25. A model that calls more genomes resistant moves both figures this way.

**The known-gene rule moves too.** It is a fixed rule with no training, so it cannot have learned the cohort. It can be read on 74 datasets, and it still reads below on at least one measure on 35 of them. Where a fixed rule's figure moves, the new genomes and their labels differ from the cohort's.

**One study holds most of a batch.** On 69 of the 82 datasets, a single study holds at least half of the resistant genomes. An estimate on such a dataset is mostly that study's.

So a report's figure describes the cohort it was measured on. It did not predict the figure on the next deposited study, in either direction. The test does not say why a figure moved.

The same-study population is 3,155 genomes on 49 datasets in 6 species. Its genomes are new, and they come from studies the cohort already holds. 3,144 of the 3,155 genomes were scored. The other 11 are *K. pneumoniae* genomes that failed the species check, and the plan excludes a genome that fails it.

- The set covering machine holds on 36 of the 49 datasets.
- The decision tree holds on 31 of the 49 datasets.
- The lasso can be scored and read on 15 of the 49 datasets, and holds on 9 of them.

34 datasets are in both populations. On each of them one final model and one report figure are read on two sets of new genomes. The plan scores the two populations apart and does not compare them, so this comparison was added after both results were read. The set covering machine holds on 29 of the 34 in the same-study population and on 10 of the 34 in the primary population. No dataset holds in the primary population alone.

The same-study batches are smaller, and a smaller batch has a wider interval, so it reads below less often. A count with no interval gives the same direction. The set covering machine's sensitivity is higher in the same-study population on 27 of the 34 datasets, and its specificity on 29.

The lineage split holds out lineages. It does not hold out studies. Its figure fits new genomes from a study the cohort already holds more often than it fits a new study.

## Runs

Each species ran once, on 12 Nov 2026, on one machine of its own, from a fresh clone of this repository at the kit commit. The run fetched the frozen genomes from BV-BRC, checked their species, counted their k-mers, and applied each report's final model with `grm predict`. The reports are the ones the site serves, at evaluation revision 12, each pinned by SHA-256.

- **Learners.** The set covering machine and the decision tree are scored on all 82 primary datasets. The lasso is scored on 26 of them, and the known-gene rule on 75. On the 49 same-study datasets the lasso is scored on 15 and the other three on all 49.
- **Not scored.** The lasso needs a complete published model, and the known-gene rule needs a published rule. Random forests and boosted trees are never scored, because a report does not hold those models. AMRFinderPlus was not run.
- **Replay first.** Before a new genome was called, every final model was run on its own cohort genomes, and its calls had to equal the report's. The known-gene rule was replayed on the cohort the same way.
- **Scored twice.** Each run's output files were made a second time from another fresh clone of the kit commit, on the saved work directory. For all 13 runs the second set is byte-identical to the first.
- **Repeated once.** *S. aureus* was run a second time from the fetch on, on another machine. The run of record is the first run.

In the repeat, the 46 files compared by bytes are identical, and they include every prediction and every k-mer table. The stage receipts differ in the hashes of 79 of the 803 fetched FASTA files. BV-BRC returned those genomes' contigs in another order, and the k-mer sets are the same.

The two populations are scored apart and are never added together. The primary population holds genomes from studies the cohort does not hold. The same-study population holds newly added genomes from studies the cohort already holds.

## Tables

`scripts/outside-study-results/tables.py` prints these tables from the output files in `reference/outside-study/`. It estimates nothing: each figure is a cell of those files, a count of their rows, or a cell shown as a percentage.

<!-- tables:begin -->
**Primary population**

| Species | Frozen genomes | Scored | Borderline in the species check | Excluded | Datasets | With estimates | Studies on a dataset |
|---|---:|---:|---:|---:|---:|---:|---|
| *E. coli* | 601 | 601 | 328 | 0 | 18 | 18 | 2 to 4 |
| *K. pneumoniae* | 2,014 | 2,014 | 0 | 0 | 21 | 21 | 1 to 12 |
| *M. tuberculosis* | 1,152 | 1,152 | 1 | 0 | 11 | 11 | 4 to 7 |
| *S. enterica* | 2,107 | 2,107 | 0 | 0 | 17 | 17 | 2 to 39 |
| *S. flexneri* | 78 | 78 | 0 | 0 | 1 | 1 | 1 |
| *S. aureus* | 803 | 803 | 0 | 0 | 7 | 7 | 1 to 3 |
| *S. pneumoniae* | 503 | 503 | 3 | 0 | 7 | 7 | 2 to 3 |
| All 7 | 7,258 | 7,258 | 332 | 0 | 82 | 82 | |

Scored genomes whose k-mer set equals another genome's (listed, never excluded), and the attempts each stage of the run took.

| Species | Same k-mer set as a cohort genome | Same k-mer set as another scored genome | Attempts per stage |
|---|---:|---:|---|
| *E. coli* | 0 | 0 | 1 |
| *K. pneumoniae* | 0 | 0 | 1 |
| *M. tuberculosis* | 0 | 0 | 1 |
| *S. enterica* | 0 | 10 | 1 |
| *S. flexneri* | 0 | 0 | 1 |
| *S. aureus* | 0 | 0 | 1 |
| *S. pneumoniae* | 0 | 2 | 1 |

Where the scored genomes of a dataset with estimates come from. A study is a BioProject.

| Species | Datasets with estimates | One study holds at least half of the resistant genomes | One study holds at least half of the susceptible genomes | Most resistant and most susceptible genomes come from different studies | No study can be left out |
|---|---:|---:|---:|---:|---:|
| *E. coli* | 18 | 18 | 13 | 2 | 0 |
| *K. pneumoniae* | 21 | 11 | 12 | 17 | 1 |
| *M. tuberculosis* | 11 | 11 | 11 | 1 | 0 |
| *S. enterica* | 17 | 14 | 4 | 16 | 2 |
| *S. flexneri* | 1 | 1 | 1 | 0 | 1 |
| *S. aureus* | 7 | 7 | 7 | 1 | 1 |
| *S. pneumoniae* | 7 | 7 | 7 | 7 | 0 |
| All 7 | 82 | 69 | 55 | 44 | 5 |

The study that most often holds the most resistant genomes of a dataset with estimates, and its share of them there.

| Species | Study | Datasets where it holds the most resistant genomes | Its share of them, lowest to highest |
|---|---|---:|---|
| *E. coli* | PRJNA841811 | 18 of 18 | 71.6% to 90.6% |
| *K. pneumoniae* | PRJEB29742 | 11 of 21 | 24.0% to 90.6% |
| *M. tuberculosis* | PRJEB48275 | 9 of 11 | 51.0% to 95.1% |
| *S. enterica* | PRJEB67737 | 7 of 17 | 49.7% to 86.7% |
| *S. flexneri* | PRJNA1069659 | 1 of 1 | 100.0% |
| *S. aureus* | PRJNA1050907 | 7 of 7 | 74.5% to 100.0% |
| *S. pneumoniae* | PRJEB82317 | 5 of 7 | 53.8% to 100.0% |

Sensitivity against the report's lineage-split figure, in datasets with estimates. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. coli* | 10 / 3 / 5 | 10 / 4 / 4 | 3 / 0 / 0 | 2 / 9 / 7 |
| *K. pneumoniae* | 11 / 2 / 8 | 7 / 1 / 13 | 2 / 1 / 1 | 11 / 5 / 5 |
| *M. tuberculosis* | 1 / 4 / 6 | 0 / 2 / 9 | 0 / 3 / 1 | 0 / 5 / 0 |
| *S. enterica* | 4 / 5 / 8 | 9 / 4 / 4 | 0 / 1 / 3 | 4 / 2 / 10 |
| *S. flexneri* | 0 / 1 / 0 | 0 / 1 / 0 | 0 / 1 / 0 | 1 / 0 / 0 |
| *S. aureus* | 1 / 1 / 4 / 1 | 1 / 1 / 4 / 1 | 0 / 1 / 1 / 1 | 0 / 4 / 2 / 1 |
| *S. pneumoniae* | 0 / 2 / 5 | 0 / 3 / 4 | 0 / 4 / 3 | 1 / 3 / 3 |
| All 7 | 27 / 18 / 36 / 1 | 27 / 16 / 38 / 1 | 5 / 11 / 9 / 1 | 19 / 28 / 27 / 1 |

Specificity against the report's lineage-split figure, in datasets with estimates. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. coli* | 13 / 3 / 2 | 13 / 4 / 1 | 2 / 1 / 0 | 11 / 7 / 0 |
| *K. pneumoniae* | 9 / 8 / 4 | 15 / 4 / 2 | 2 / 1 / 1 | 6 / 8 / 7 |
| *M. tuberculosis* | 3 / 2 / 6 | 4 / 2 / 5 | 3 / 0 / 1 | 0 / 5 / 0 |
| *S. enterica* | 3 / 5 / 9 | 2 / 4 / 11 | 2 / 0 / 2 | 5 / 4 / 7 |
| *S. flexneri* | 1 / 0 / 0 | 0 / 1 / 0 | 1 / 0 / 0 | 0 / 0 / 1 |
| *S. aureus* | 0 / 2 / 4 / 1 | 2 / 1 / 3 / 1 | 0 / 1 / 1 / 1 | 0 / 5 / 1 / 1 |
| *S. pneumoniae* | 6 / 1 / 0 | 6 / 0 / 1 | 6 / 0 / 1 | 4 / 3 / 0 |
| All 7 | 35 / 21 / 25 / 1 | 42 / 16 / 23 / 1 | 16 / 3 / 6 / 1 | 26 / 32 / 16 / 1 |

The two readings of a row together, over all the species above. A row is one dataset with estimates and one learner.

| Learner | Rows read | Neither reading below | Sensitivity below only | Specificity below only | Both below |
|---|---:|---:|---:|---:|---:|
| Set covering machine | 81 | 32 | 14 | 22 | 13 |
| Decision tree | 81 | 23 | 16 | 31 | 11 |
| Lasso | 25 | 7 | 2 | 13 | 3 |
| Known-gene rule | 74 | 39 | 9 | 16 | 10 |

Datasets where a learner holds, of the datasets read for it. A learner holds when neither its sensitivity nor its specificity reads below.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. coli* | 1 of 18 | 0 of 18 | 0 of 3 | 6 of 18 |
| *K. pneumoniae* | 8 of 21 | 4 of 21 | 1 of 4 | 9 of 21 |
| *M. tuberculosis* | 7 of 11 | 7 of 11 | 1 of 4 | 5 of 5 |
| *S. enterica* | 10 of 17 | 7 of 17 | 2 of 4 | 10 of 16 |
| *S. flexneri* | 0 of 1 | 1 of 1 | 0 of 1 | 0 of 1 |
| *S. aureus* | 5 of 6 | 3 of 6 | 2 of 2 | 6 of 6 |
| *S. pneumoniae* | 1 of 7 | 1 of 7 | 1 of 7 | 3 of 7 |
| All 7 | 32 of 81 | 23 of 81 | 7 of 25 | 39 of 74 |

Datasets with estimates where the known-gene rule calls every scored genome the same way.

| Species | Rule rows with estimates | Calls no genome resistant | Calls every genome resistant |
|---|---:|---:|---:|
| *E. coli* | 18 | 0 | 0 |
| *K. pneumoniae* | 21 | 0 | 0 |
| *M. tuberculosis* | 5 | 3 | 1 |
| *S. enterica* | 16 | 0 | 0 |
| *S. flexneri* | 1 | 1 | 0 |
| *S. aureus* | 7 | 3 | 1 |
| *S. pneumoniae* | 7 | 3 | 0 |
| All 7 | 75 | 10 | 2 |

**Same-study population**

| Species | Frozen genomes | Scored | Borderline in the species check | Excluded | Datasets | With estimates | Studies on a dataset |
|---|---:|---:|---:|---:|---:|---:|---|
| *E. cloacae* | 71 | 71 | 25 | 0 | 2 | 2 | 1 |
| *E. faecium* | 368 | 368 | 10 | 0 | 5 | 5 | 1 |
| *E. coli* | 87 | 87 | 8 | 0 | 13 | 13 | 1 |
| *K. pneumoniae* | 470 | 459 | 0 | 11 | 17 | 17 | 1 |
| *S. enterica* | 1,741 | 1,741 | 0 | 0 | 8 | 8 | 1 to 4 |
| *S. pneumoniae* | 418 | 418 | 0 | 0 | 4 | 4 | 1 |
| All 6 | 3,155 | 3,144 | 43 | 11 | 49 | 49 | |

Scored genomes whose k-mer set equals another genome's (listed, never excluded), and the attempts each stage of the run took.

| Species | Same k-mer set as a cohort genome | Same k-mer set as another scored genome | Attempts per stage |
|---|---:|---:|---|
| *E. cloacae* | 0 | 0 | 1 |
| *E. faecium* | 0 | 0 | 1 |
| *E. coli* | 0 | 0 | 1 |
| *K. pneumoniae* | 0 | 0 | 1 |
| *S. enterica* | 0 | 0 | 1 |
| *S. pneumoniae* | 0 | 0 | 1 |

Where the scored genomes of a dataset with estimates come from. A study is a BioProject.

| Species | Datasets with estimates | One study holds at least half of the resistant genomes | One study holds at least half of the susceptible genomes | Most resistant and most susceptible genomes come from different studies | No study can be left out |
|---|---:|---:|---:|---:|---:|
| *E. cloacae* | 2 | 2 | 2 | 0 | 2 |
| *E. faecium* | 5 | 5 | 5 | 0 | 5 |
| *E. coli* | 13 | 13 | 13 | 0 | 13 |
| *K. pneumoniae* | 17 | 17 | 17 | 0 | 17 |
| *S. enterica* | 8 | 8 | 8 | 0 | 2 |
| *S. pneumoniae* | 4 | 4 | 4 | 0 | 4 |
| All 6 | 49 | 49 | 49 | 0 | 43 |

The study that most often holds the most resistant genomes of a dataset with estimates, and its share of them there.

| Species | Study | Datasets where it holds the most resistant genomes | Its share of them, lowest to highest |
|---|---|---:|---|
| *E. cloacae* | PRJNA288601 | 2 of 2 | 100.0% |
| *E. faecium* | PRJNA771404 | 5 of 5 | 100.0% |
| *E. coli* | PRJNA288601 | 13 of 13 | 100.0% |
| *K. pneumoniae* | PRJNA288601 | 17 of 17 | 100.0% |
| *S. enterica* | PRJNA248792 | 8 of 8 | 99.9% to 100.0% |
| *S. pneumoniae* | PRJEB3084 | 4 of 4 | 100.0% |

Sensitivity against the report's lineage-split figure, in datasets with estimates. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. cloacae* | 2 / 0 / 0 | 1 / 1 / 0 | 1 / 1 / 0 | 2 / 0 / 0 |
| *E. faecium* | 2 / 2 / 1 | 2 / 3 / 0 | 1 / 0 / 1 | 0 / 4 / 1 |
| *E. coli* | 2 / 11 / 0 | 3 / 10 / 0 | 0 / 2 / 0 | 1 / 12 / 0 |
| *K. pneumoniae* | 2 / 8 / 7 | 4 / 2 / 11 | 0 / 3 / 0 | 7 / 10 / 0 |
| *S. enterica* | 2 / 2 / 4 | 2 / 2 / 4 | 1 / 0 / 1 | 2 / 2 / 4 |
| *S. pneumoniae* | 0 / 3 / 1 | 0 / 2 / 2 | 1 / 1 / 2 | 0 / 4 / 0 |
| All 6 | 10 / 26 / 13 | 12 / 20 / 17 | 4 / 7 / 4 | 12 / 32 / 5 |

Specificity against the report's lineage-split figure, in datasets with estimates. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. cloacae* | 0 / 1 / 1 | 0 / 2 / 0 | 0 / 1 / 1 | 1 / 1 / 0 |
| *E. faecium* | 0 / 2 / 3 | 0 / 4 / 1 | 0 / 1 / 1 | 2 / 2 / 1 |
| *E. coli* | 1 / 9 / 3 | 1 / 9 / 3 | 0 / 1 / 1 | 2 / 10 / 1 |
| *K. pneumoniae* | 1 / 7 / 9 | 3 / 6 / 8 | 0 / 0 / 3 | 2 / 3 / 12 |
| *S. enterica* | 1 / 1 / 6 | 0 / 1 / 7 | 0 / 1 / 1 | 0 / 0 / 8 |
| *S. pneumoniae* | 1 / 2 / 1 | 2 / 1 / 1 | 2 / 1 / 1 | 1 / 3 / 0 |
| All 6 | 4 / 22 / 23 | 6 / 23 / 20 | 2 / 5 / 8 | 8 / 19 / 22 |

The two readings of a row together, over all the species above. A row is one dataset with estimates and one learner.

| Learner | Rows read | Neither reading below | Sensitivity below only | Specificity below only | Both below |
|---|---:|---:|---:|---:|---:|
| Set covering machine | 49 | 36 | 9 | 3 | 1 |
| Decision tree | 49 | 31 | 12 | 6 | 0 |
| Lasso | 15 | 9 | 4 | 2 | 0 |
| Known-gene rule | 49 | 30 | 11 | 7 | 1 |

Datasets where a learner holds, of the datasets read for it. A learner holds when neither its sensitivity nor its specificity reads below.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. cloacae* | 0 of 2 | 1 of 2 | 1 of 2 | 0 of 2 |
| *E. faecium* | 3 of 5 | 3 of 5 | 1 of 2 | 3 of 5 |
| *E. coli* | 11 of 13 | 9 of 13 | 2 of 2 | 10 of 13 |
| *K. pneumoniae* | 14 of 17 | 10 of 17 | 3 of 3 | 8 of 17 |
| *S. enterica* | 5 of 8 | 6 of 8 | 1 of 2 | 6 of 8 |
| *S. pneumoniae* | 3 of 4 | 2 of 4 | 1 of 4 | 3 of 4 |
| All 6 | 36 of 49 | 31 of 49 | 9 of 15 | 30 of 49 |

Datasets with estimates where the known-gene rule calls every scored genome the same way.

| Species | Rule rows with estimates | Calls no genome resistant | Calls every genome resistant |
|---|---:|---:|---:|
| *E. cloacae* | 2 | 0 | 0 |
| *E. faecium* | 5 | 1 | 1 |
| *E. coli* | 13 | 0 | 0 |
| *K. pneumoniae* | 17 | 0 | 0 |
| *S. enterica* | 8 | 0 | 0 |
| *S. pneumoniae* | 4 | 2 | 0 |
| All 6 | 49 | 3 | 1 |

**Datasets scored in both populations**

A dataset in both populations has one final model and one report figure, read on two sets of new genomes. A learner holds when neither its sensitivity nor its specificity reads below. No genome counts are added. The plan does not make this comparison; it was added after both results were read.

| Learner | Rows read in both | Holds in both | Holds in the same-study population only | Holds in the primary population only | Holds in neither |
|---|---:|---:|---:|---:|---:|
| Set covering machine | 34 | 10 | 19 | 0 | 5 |
| Decision tree | 34 | 4 | 19 | 1 | 10 |
| Lasso | 8 | 1 | 4 | 0 | 3 |
| Known-gene rule | 34 | 12 | 13 | 3 | 6 |

The same rows by the estimate alone, with no interval. A small batch has a wide interval, so it reads below less often; this count does not depend on the batch's size. Each cell counts rows as higher in the same-study population / equal / higher in the primary population.

| Learner | Rows read in both | Sensitivity | Specificity |
|---|---:|---|---|
| Set covering machine | 34 | 27 / 0 / 7 | 29 / 0 / 5 |
| Decision tree | 34 | 25 / 0 / 9 | 30 / 0 / 4 |
| Lasso | 8 | 5 / 0 / 3 | 4 / 0 / 4 |
| Known-gene rule | 34 | 13 / 2 / 19 | 23 / 2 / 9 |

**Rows cited in the text**

Each estimate is followed by its 95% interval and its counts. The report's figure is the lineage-split figure for a learner, and the whole-cohort figure for the known-gene rule.

| Population | Species | Drug | Learner | Sensitivity | Report | Reading | Specificity | Report | Reading |
|---|---|---|---|---|---|---|---|---|---|
| primary | *S. aureus* | oxacillin | Set covering machine | 98.9% (97.8% to 99.6%), 652 of 659 | 90.3% | above | 94.6% (85.1% to 98.9%), 53 of 56 | 98.3% | consistent |
| primary | *S. aureus* | tetracycline | Known-gene rule | 100.0% (96.5% to 100.0%), 104 of 104 | 100.0% | consistent | 0.0% (0.0% to 0.6%), 0 of 627 | 0.0% | consistent |
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | 96.7% (90.6% to 99.3%), 87 of 90 | 68.7% | above | 69.7% (65.0% to 74.1%), 287 of 412 | 97.7% | below |
| primary | *S. pneumoniae* | penicillin | Set covering machine | 97.7% (88.0% to 99.9%), 43 of 44 | 95.9% | consistent | 76.0% (70.8% to 80.6%), 234 of 308 | 92.4% | below |
| primary | *E. coli* | amikacin | Set covering machine | 9.4% (3.5% to 19.3%), 6 of 64 | 46.7% | below | 93.9% (91.6% to 95.8%), 496 of 528 | 98.4% | below |
| primary | *E. coli* | ceftazidime | Set covering machine | 22.9% (17.0% to 29.8%), 41 of 179 | 87.2% | below | 94.8% (92.0% to 96.9%), 331 of 349 | 92.1% | consistent |
| primary | *E. coli* | sulfamethoxazole/trimethoprim | Set covering machine | 93.4% (90.8% to 95.5%), 428 of 458 | 89.8% | above | 39.3% (30.4% to 48.8%), 46 of 117 | 83.5% | below |
| primary | *M. tuberculosis* | ethambutol | Set covering machine | 86.0% (81.9% to 89.5%), 296 of 344 | 79.6% | above | 31.0% (26.7% to 35.5%), 138 of 445 | 60.9% | below |
| primary | *M. tuberculosis* | pyrazinamide | Set covering machine | 21.7% (18.0% to 25.9%), 96 of 442 | 58.0% | below | 68.4% (63.6% to 72.9%), 279 of 408 | 55.9% | above |
| primary | *M. tuberculosis* | levofloxacin | Set covering machine | 94.3% (88.5% to 97.7%), 115 of 122 | 59.8% | above | 94.0% (91.8% to 95.7%), 606 of 645 | 79.2% | above |
| primary | *S. enterica* | sulfamethoxazole | Set covering machine | 30.6% (25.2% to 36.5%), 83 of 271 | 87.0% | below | 97.5% (92.9% to 99.5%), 117 of 120 | 93.4% | consistent |
| primary | *S. enterica* | sulfamethoxazole | Lasso | 93.0% (89.3% to 95.7%), 252 of 271 | 75.8% | above | 98.3% (94.1% to 99.8%), 118 of 120 | 88.2% | above |
| primary | *S. enterica* | sulfamethoxazole | Known-gene rule | 98.9% (96.8% to 99.8%), 268 of 271 | 93.5% | above | 100.0% (97.0% to 100.0%), 120 of 120 | 99.5% | consistent |
| primary | *S. enterica* | tetracycline | Set covering machine | 98.9% (97.4% to 99.6%), 435 of 440 | 77.0% | above | 43.0% (38.0% to 48.1%), 169 of 393 | 96.1% | below |
| primary | *S. enterica* | tetracycline | Known-gene rule | 98.9% (97.4% to 99.6%), 435 of 440 | 92.8% | above | 98.5% (96.7% to 99.4%), 387 of 393 | 91.7% | above |
| primary | *K. pneumoniae* | meropenem | Set covering machine | 61.3% (58.3% to 64.2%), 658 of 1,074 | 85.6% | below | 91.5% (88.9% to 93.7%), 517 of 565 | 95.7% | below |
| primary | *K. pneumoniae* | meropenem | Known-gene rule | 90.5% (88.6% to 92.2%), 972 of 1,074 | 95.0% | below | 88.8% (86.0% to 91.3%), 502 of 565 | 91.5% | below |
| primary | *K. pneumoniae* | ertapenem | Decision tree | 93.5% (91.2% to 95.4%), 533 of 570 | 82.9% | above | 12.6% (9.6% to 16.1%), 54 of 428 | 89.0% | below |
| primary | *K. pneumoniae* | tigecycline | Set covering machine | 97.7% (94.7% to 99.3%), 213 of 218 | 23.7% | above | 3.7% (2.3% to 5.8%), 19 of 509 | 66.4% | below |
| primary | *S. flexneri* | ciprofloxacin | Set covering machine | 74.2% (55.4% to 88.1%), 23 of 31 | 69.7% | consistent | 83.0% (69.2% to 92.4%), 39 of 47 | 95.7% | below |
| same-study | *K. pneumoniae* | meropenem | Set covering machine | 92.9% (88.9% to 95.8%), 223 of 240 | 85.6% | above | 97.6% (94.4% to 99.2%), 201 of 206 | 95.7% | consistent |
| same-study | *E. coli* | ceftazidime | Set covering machine | 80.6% (64.0% to 91.8%), 29 of 36 | 87.2% | consistent | 95.8% (85.7% to 99.5%), 46 of 48 | 92.1% | consistent |
| same-study | *S. pneumoniae* | penicillin | Set covering machine | 98.8% (95.8% to 99.9%), 169 of 171 | 95.9% | consistent | 97.9% (95.2% to 99.3%), 234 of 239 | 92.4% | above |
| same-study | *S. pneumoniae* | penicillin | Lasso | 29.8% (23.1% to 37.3%), 51 of 171 | 95.4% | below | 99.2% (97.0% to 99.9%), 237 of 239 | 98.3% | consistent |
| same-study | *E. faecium* | gentamicin | Set covering machine | 13.3% (3.8% to 30.7%), 4 of 30 | 84.6% | below | 96.6% (94.0% to 98.3%), 313 of 324 | 89.1% | above |
| same-study | *E. cloacae* | sulfamethoxazole/trimethoprim | Set covering machine | 59.3% (38.8% to 77.6%), 16 of 27 | 100.0% | below | 90.9% (78.3% to 97.5%), 40 of 44 | 94.9% | consistent |
| same-study | *S. enterica* | amoxicillin/clavulanate | Set covering machine | 0.0% (0.0% to 1.2%), 0 of 301 | 91.7% | below | 100.0% (99.5% to 100.0%), 716 of 716 | 99.7% | consistent |
| same-study | *S. enterica* | amoxicillin/clavulanate | Known-gene rule | 15.9% (12.0% to 20.6%), 48 of 301 | 93.6% | below | 100.0% (99.5% to 100.0%), 716 of 716 | 91.2% | above |
| same-study | *S. enterica* | trimethoprim | Set covering machine | 4.0% (2.2% to 6.8%), 13 of 322 | 82.0% | below | 100.0% (99.7% to 100.0%), 1,382 of 1,382 | 78.0% | above |
| same-study | *S. enterica* | trimethoprim | Decision tree | 99.7% (98.3% to 100.0%), 321 of 322 | 87.9% | above | 99.9% (99.6% to 100.0%), 1,381 of 1,382 | 66.7% | above |

**Cited rows, study by study**

The same counts as the row above them in `metrics.tsv`, split by the genome's BioProject.

| Population | Species | Drug | Learner | Study | Resistant genomes | Of them called resistant | Susceptible genomes | Of them called resistant |
|---|---|---|---|---|---:|---:|---:|---:|
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | PRJEB82317 | 42 | 42 | 95 | 79 |
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | PRJNA1357085 | 48 | 45 | 55 | 45 |
| primary | *S. pneumoniae* | chloramphenicol | Set covering machine | PRJNA320861 | 0 | 0 | 262 | 1 |
| primary | *S. pneumoniae* | penicillin | Set covering machine | PRJEB82317 | 44 | 43 | 36 | 3 |
| primary | *S. pneumoniae* | penicillin | Set covering machine | PRJNA1357085 | 0 | 0 | 87 | 71 |
| primary | *S. pneumoniae* | penicillin | Set covering machine | PRJNA320861 | 0 | 0 | 185 | 0 |
<!-- tables:end -->

## Primary population, species by species

The rows named here are in the table "Rows cited in the text". Every dataset's row is in the species' `metrics.tsv`, and `reading.txt` gives the same rows in plain words.

***S. aureus*, 803 genomes, 7 datasets.** One study, PRJNA1050907, holds 74.5% to 100% of the resistant genomes on every dataset. Sensitivity mostly holds or rises: for the set covering machine, oxacillin has 98.9% against 90.3%. Levofloxacin has no lineage-split figure in its report, so its rows are not read. The known-gene rule is one-sided on four of the seven datasets. On tetracycline it calls all 731 genomes resistant.

***S. pneumoniae*, 503 genomes, 7 datasets.** Sensitivity never reads below for a learner. Specificity reads below on 18 of the 21 learner rows. The false positives sit in two of the three studies. On chloramphenicol the set covering machine calls 1 of 262 susceptible genomes resistant in PRJNA320861, 79 of 95 in PRJEB82317 and 45 of 55 in PRJNA1357085.

On every *S. pneumoniae* dataset, most resistant genomes come from one study and most susceptible genomes from another. So these rows compare studies as much as they measure a model.

***E. coli*, 601 genomes, 18 datasets.** This is the weakest batch for the learners. The set covering machine reads below on sensitivity on 10 datasets and on specificity on 13. Its amikacin sensitivity is 9.4% against 46.7%, and its ceftazidime sensitivity is 22.9% against 87.2%. Two facts limit this batch. One study, PRJNA841811, holds 71.6% to 90.6% of the resistant genomes on every dataset. The species check marks 328 of the 601 genomes as borderline, and the plan keeps and flags a borderline genome, so they are scored.

***M. tuberculosis*, 1,152 genomes, 11 datasets.** Sensitivity mostly holds or rises. Specificity falls on some drugs: for the set covering machine, ethambutol has 31.0% against 60.9%. Pyrazinamide goes the other way, with 21.7% sensitivity against 58.0%. One study, PRJEB48275, holds at least half of the resistant genomes on 9 of the 11 datasets.

***Salmonella enterica*, 2,107 genomes, 17 datasets.** This batch has the most studies, up to 39 on a dataset. On two drugs the learned rules fail where the known-gene rule holds. On sulfamethoxazole the set covering machine has 30.6% sensitivity against 87.0%, and the known-gene rule has 98.9%. On tetracycline the set covering machine's specificity is 43.0% against 96.1%, and the known-gene rule's is 98.5%. The lasso holds on sulfamethoxazole, with 93.0%.

***K. pneumoniae*, 2,014 genomes, 21 datasets.** No single study holds this batch: PRJEB29742 holds the most resistant genomes on 11 of the 21 datasets, with 24.0% to 90.6% of them. The set covering machine reads below on sensitivity on 11 datasets: meropenem has 61.3% against 85.6%, while the known-gene rule has 90.5%.

The decision tree reads below on specificity on 15 datasets, and on ertapenem its specificity is 12.6% against 89.0%. On tigecycline the set covering machine's specificity is 3.7% against 66.4%, with sensitivity at 97.7% against 23.7%.

***S. flexneri*, 78 genomes, 1 dataset.** Every genome is from one study. Sensitivity is consistent with the report's figure and specificity reads below it, on 31 resistant and 47 susceptible genomes.

## Same-study population, species by species

On every dataset here, one study holds at least 99.9% of the resistant genomes. *E. cloacae*, *E. coli* and *K. pneumoniae* share one study, PRJNA288601, so their 32 datasets are one study's.

***E. cloacae*, 71 genomes, 2 datasets.** This is the smallest batch, and the species check marks 25 of its 71 genomes as borderline. The set covering machine reads below on sensitivity on both datasets. On sulfamethoxazole/trimethoprim it has 59.3% against 100.0%.

***E. faecium*, 368 genomes, 5 datasets.** The study is PRJNA771404. The set covering machine holds on 3 of the 5 datasets. Its gentamicin sensitivity is 13.3% against 84.6%, on 30 resistant genomes.

***E. coli*, 87 genomes, 13 datasets.** The set covering machine holds on 11 of the 13 datasets. Its ceftazidime sensitivity is 80.6% against 87.2%, where the primary population gave 22.9%. A dataset here has 28 to 56 resistant genomes, so its intervals are wide.

***K. pneumoniae*, 459 genomes scored of 470, 17 datasets.** The 11 genomes that failed the species check are each nearer the *K. quasipneumoniae* reference. The set covering machine holds on 14 of the 17 datasets. Its meropenem sensitivity is 92.9% against 85.6%, where the primary population gave 61.3%.

***Salmonella enterica*, 1,741 genomes, 8 datasets.** One study, PRJNA248792, holds 99.9% to 100% of the resistant genomes on every dataset. The set covering machine holds on 5 of the 8 datasets, and two of its rows fail outright. On amoxicillin/clavulanate it calls none of the 301 resistant genomes resistant, and the known-gene rule finds 15.9% of them against 93.6%. On trimethoprim its sensitivity is 4.0% against 82.0%, where the decision tree has 99.7%.

***S. pneumoniae*, 418 genomes, 4 datasets.** The study is PRJEB3084. The set covering machine holds on 3 of the 4 datasets. Its penicillin specificity is 97.9% against 92.4%, where the primary population gave 76.0%. The lasso's penicillin sensitivity is 29.8% against 95.4%.

## Limits

- **Same database.** The genomes come from BV-BRC, like the cohort, with labels made under the same practices. The batch's metadata and label totals were read before the plan was frozen.
- **One model, no refit.** Each estimate comes from one final model for each dataset. Its interval covers the sampling of these genomes only. It does not cover the differences between studies.
- **No test.** The reading compares one figure with one interval. There is no hypothesis test and no correction for the number of rows. The datasets of a species share genomes and studies, so their rows are not independent.
- **Few studies.** A study that holds most of a batch decides its estimate. On 44 of the 82 primary datasets, most resistant genomes come from one study and most susceptible genomes from another, so the two measures come from different studies.
- **Batch size.** A row reads below only when the report's figure is above the whole interval. A small batch has a wide interval, so it reads below less often. The same-study batches of *E. cloacae* and *E. coli* have 71 and 87 genomes.
- **Two populations.** The populations differ in species and in datasets. Only the 34 datasets in both compare one model on two sets of genomes, and there the batches still differ in size and in study.
- **Not in the plan.** The plan fixes the reading of each measure on each dataset. The count of datasets where a learner holds, and the comparison of the two populations, were added after the results were read. No genome is pooled.
- **Two kinds of figure.** The report's figure is for models refitted in each fold. The new figure is for the one final model. For the known-gene rule the report's figure is the fixed rule on the whole cohort.
- **Precision is not read.** It depends on the share of resistant genomes, which differs between the cohort and the new genomes.
- **Three learners.** Random forests and boosted trees are not scored, and the lasso is scored on part of the datasets.
- **No cause.** A k-mer that marks resistance in the cohort can mark a lineage or a plasmid there. The test shows that a figure moved. It does not show why.

## Supported statements

It supports four statements.

1. A report's cross-validation figure is an estimate for the cohort's own population. On newly deposited studies it was often outside the new estimate's interval, with the lineage split too.
2. The direction is not fixed. Sensitivity read above the report's figure more often than below it, and specificity read below more often than above.
3. The known-gene rule is not a stable yardstick either, and it is one-sided on 12 of its 75 datasets.
4. The figures held more often on new genomes from a study the cohort already holds. On the 34 datasets in both populations, the set covering machine holds on 29 with same-study genomes and on 10 with genomes from other studies.

It does not support a claim that a model works, or fails, in a clinical lab. It names no gene as a cause.

## Files

`reference/outside-study/<population>/<species>/` holds the kit's 11 output files for each run, as the run wrote them. `results.json` is stored gzipped. `reference/outside-study/SHA256SUMS` lists every file.

```sh
python3 -I scripts/outside-study-results/tables.py reference/outside-study
python3 -I scripts/outside-study-results/tables.py reference/outside-study --check docs/outside-study-results.md
python3 -I -m unittest discover -s tests/outside_study_results -t . -v
```

The second command exits 1 unless the tables above are exactly what the files give. The suite runs it in CI. It also checks each run's record against the frozen population: the kit commit, the population, the species, and the genome and dataset lists.
