# Outside-study test: results

This is the result of a retrospective test on newly deposited studies from the same database. The plan is [outside-study-test.md](outside-study-test.md). It was frozen at commit `e3e1c34` before any of these genomes was fetched. The test is not a lockbox, and it is not external validation. The plan says why.

## The question

Each dataset report gives a learner's sensitivity and specificity from cross-validation on the cohort. The lineage split is the stricter of the two splits, because it holds out whole lineages. The test asks whether the lineage-split figures hold when the report's final model is applied, unchanged, to genomes from studies the cohort does not hold.

## Result

On genomes from studies the cohort does not hold, each learner's figures hold on fewer than half of the datasets read for it. On new genomes from studies the cohort already holds, each learner's figures hold on more than half.

Those two sentences are about the three learners: the set covering machine, the decision tree and the lasso. The two populations cover different species and datasets. The plan scores them apart and does not compare them, so setting them side by side was added after both results were read.

Both sentences read each estimate against the report its run pinned, at evaluation revision 12. The site has served revision 14 since 13 Nov 2026. A revision 14 report can hold two figures for the lasso. The plan chooses neither, because it reads the revision 12 reports only.

Against either figure the first sentence stands for each learner. The second stands for the set covering machine and the decision tree. For the lasso it depends on the figure read. The lasso holds on 9 of the 15 same-study datasets read for it against the figure over fitted folds, which the runs read, and on 7 of the 15 against the figure over every genome. [Revision 14](#revision-14) gives the rows.

The primary population is 7,258 genomes on 82 datasets in 7 species. Every genome was scored and none was excluded. One dataset has no lineage-split figure in its report, so 81 datasets can be read.

A reading is "below" when the report's figure is above the whole 95% interval of the new estimate. For a learner the report's figure is the lineage-split figure. For the known-gene rule it is the fixed rule's figure on the whole cohort. The plan reads sensitivity and specificity one at a time. This document also joins them: a learner holds on a dataset when neither its sensitivity nor its specificity reads below.

- The set covering machine holds on 32 of the 81 datasets and reads below on the other 49.
- The decision tree holds on 23 of the 81 datasets.
- The lasso can be scored and read on 25 of the 81 datasets, and holds on 7 of them.

Three facts belong beside those counts.

**The figures move both ways.** For the set covering machine, sensitivity reads above the report's figure on 36 of the 81 datasets and below it on 27. Its specificity reads below on 35 datasets and above on 25. Only its specificity reads below on 22 datasets, only its sensitivity on 14, and both on 13.

**The known-gene rule moves too.** It is a fixed rule with no training, and its report figure is the rule's figure on the whole cohort. It can be read on 74 datasets. It holds on 39 of them and reads below on at least one measure on the other 35. On 12 of the 39 its report figure is exactly 0% or 100% on a measure, and such a figure can be read one way only.

**One study holds most of a batch.** On 69 of the 82 datasets, a single study holds at least half of the resistant genomes. A sensitivity estimate on such a dataset is mostly that study's.

For the set covering machine, both readings are consistent on 3 of the 81 datasets. On the other 78 at least one figure moved, up or down. The test does not say why a figure moved.

The same-study population is 3,155 genomes on 49 datasets in 6 species. Its genomes are new, and they come from studies the cohort already holds. 3,144 of the 3,155 genomes were scored. The other 11 are *K. pneumoniae* genomes that failed the species check, and the plan excludes a genome that fails it.

- The set covering machine holds on 36 of the 49 datasets.
- The decision tree holds on 31 of the 49 datasets.
- The lasso can be scored and read on 15 of the 49 datasets, and holds on 9 of them.
- The known-gene rule holds on 30 of the 49 datasets.

34 datasets are in both populations. On each of them one final model and one report figure are read on two sets of new genomes. The plan scores the two populations apart and does not compare them, so this comparison was added after both results were read.

The set covering machine holds on 29 of the 34 in the same-study population and on 10 of the 34 in the primary population. For the set covering machine, no dataset holds in the primary population alone. The decision tree has 1 such dataset and the known-gene rule has 3.

The same-study batches are smaller, and a smaller batch has a wider interval, so it reads below less often. A count with no interval gives the same direction for the set covering machine and the decision tree. The set covering machine's sensitivity is higher in the same-study population on 27 of the 34 datasets, and its specificity on 29.

The lasso is read on 8 of the 34, and its specificity is higher in each population on 4 of them. The known-gene rule's sensitivity goes the other way: it is higher in the same-study population on 13 of the 34 and higher in the primary population on 19.

The lineage split holds out lineages. It does not hold out studies. On these 34 datasets a learner's lineage-split figure held more often on new genomes from a study the cohort already holds than on genomes from other studies.

In the same-study population, one study has at least 99.9% of the resistant genomes on each of these 34 datasets, and there are three such studies: PRJNA288601 on 23 datasets, PRJNA248792 on 7 and PRJEB3084 on 4. So the same-study side of this comparison is three studies, and the comparison was added after the results were read.

## Runs

Each run of record ran on 12 Nov 2026, on one machine of its own, from a fresh clone of this repository at the kit commit. There are 13 of them, one for each species of a population. The run fetched the frozen genomes from BV-BRC, checked their species, counted their k-mers, and applied each report's final model with `grm predict`. The reports are the ones the site served on that day, at evaluation revision 12, each pinned by SHA-256.

- **Learners.** The set covering machine and the decision tree are scored on every dataset. The lasso is scored on 26 of the 82 primary datasets and on 15 of the 49 same-study datasets. The known-gene rule is scored on 75 of the 82 and on all 49.
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

Where the scored genomes of a dataset with estimates come from. A study is a BioProject. The last column counts the datasets with no reportable leave-one-study-out fold, and a dataset with one study only is always among them.

| Species | Datasets with estimates | One study holds at least half of the resistant genomes | One study holds at least half of the susceptible genomes | More than half of the resistant genomes are from one study and more than half of the susceptible genomes from another | One study only | No study can be left out |
|---|---:|---:|---:|---:|---:|---:|
| *E. coli* | 18 | 18 | 13 | 1 | 0 | 0 |
| *K. pneumoniae* | 21 | 11 | 12 | 6 | 1 | 1 |
| *M. tuberculosis* | 11 | 11 | 11 | 1 | 0 | 0 |
| *S. enterica* | 17 | 14 | 4 | 4 | 0 | 2 |
| *S. flexneri* | 1 | 1 | 1 | 0 | 1 | 1 |
| *S. aureus* | 7 | 7 | 7 | 1 | 1 | 1 |
| *S. pneumoniae* | 7 | 7 | 7 | 7 | 0 | 0 |
| All 7 | 82 | 69 | 55 | 20 | 3 | 5 |

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

Sensitivity against the report's figure, in datasets with estimates. The report's figure is the lineage-split figure for a learner, and the fixed rule's whole-cohort figure for the known-gene rule. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

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

Specificity against the report's figure, in datasets with estimates. The report's figure is the lineage-split figure for a learner, and the fixed rule's whole-cohort figure for the known-gene rule. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

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

The two readings of a row together, over all the species above. A row is one dataset with estimates and one learner. The rows with both readings consistent are part of the rows with neither reading below.

| Learner | Rows read | Neither reading below | Of them, both consistent | Sensitivity below only | Specificity below only | Both below |
|---|---:|---:|---:|---:|---:|---:|
| Set covering machine | 81 | 32 | 3 | 14 | 22 | 13 |
| Decision tree | 81 | 23 | 2 | 16 | 31 | 11 |
| Lasso | 25 | 7 | 1 | 2 | 13 | 3 |
| Known-gene rule | 74 | 39 | 21 | 9 | 16 | 10 |

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

Rows read that sit at an edge, over all the species above. A report figure of exactly 0% or 100% can be read one way only: nothing reads above 100% or below 0%. A row that gives every scored genome the same call has a new estimate of 0% on one measure and 100% on the other. Each "Of them, holds" column counts the rows of the column to its left where the learner holds.

| Learner | Rows read | Report figure at 0% or 100% on a measure | Of them, holds | Same call for every scored genome | Of them, holds |
|---|---:|---:|---:|---:|---:|
| Set covering machine | 81 | 1 | 0 | 0 | 0 |
| Decision tree | 81 | 0 | 0 | 0 | 0 |
| Lasso | 25 | 1 | 0 | 0 | 0 |
| Known-gene rule | 74 | 12 | 12 | 11 | 10 |

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

Where the scored genomes of a dataset with estimates come from. A study is a BioProject. The last column counts the datasets with no reportable leave-one-study-out fold, and a dataset with one study only is always among them.

| Species | Datasets with estimates | One study holds at least half of the resistant genomes | One study holds at least half of the susceptible genomes | More than half of the resistant genomes are from one study and more than half of the susceptible genomes from another | One study only | No study can be left out |
|---|---:|---:|---:|---:|---:|---:|
| *E. cloacae* | 2 | 2 | 2 | 0 | 2 | 2 |
| *E. faecium* | 5 | 5 | 5 | 0 | 5 | 5 |
| *E. coli* | 13 | 13 | 13 | 0 | 13 | 13 |
| *K. pneumoniae* | 17 | 17 | 17 | 0 | 17 | 17 |
| *S. enterica* | 8 | 8 | 8 | 0 | 2 | 2 |
| *S. pneumoniae* | 4 | 4 | 4 | 0 | 4 | 4 |
| All 6 | 49 | 49 | 49 | 0 | 43 | 43 |

The study that most often holds the most resistant genomes of a dataset with estimates, and its share of them there.

| Species | Study | Datasets where it holds the most resistant genomes | Its share of them, lowest to highest |
|---|---|---:|---|
| *E. cloacae* | PRJNA288601 | 2 of 2 | 100.0% |
| *E. faecium* | PRJNA771404 | 5 of 5 | 100.0% |
| *E. coli* | PRJNA288601 | 13 of 13 | 100.0% |
| *K. pneumoniae* | PRJNA288601 | 17 of 17 | 100.0% |
| *S. enterica* | PRJNA248792 | 8 of 8 | 99.9% to 100.0% |
| *S. pneumoniae* | PRJEB3084 | 4 of 4 | 100.0% |

Sensitivity against the report's figure, in datasets with estimates. The report's figure is the lineage-split figure for a learner, and the fixed rule's whole-cohort figure for the known-gene rule. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. cloacae* | 2 / 0 / 0 | 1 / 1 / 0 | 1 / 1 / 0 | 2 / 0 / 0 |
| *E. faecium* | 2 / 2 / 1 | 2 / 3 / 0 | 1 / 0 / 1 | 0 / 4 / 1 |
| *E. coli* | 2 / 11 / 0 | 3 / 10 / 0 | 0 / 2 / 0 | 1 / 12 / 0 |
| *K. pneumoniae* | 2 / 8 / 7 | 4 / 2 / 11 | 0 / 3 / 0 | 7 / 10 / 0 |
| *S. enterica* | 2 / 2 / 4 | 2 / 2 / 4 | 1 / 0 / 1 | 2 / 2 / 4 |
| *S. pneumoniae* | 0 / 3 / 1 | 0 / 2 / 2 | 1 / 1 / 2 | 0 / 4 / 0 |
| All 6 | 10 / 26 / 13 | 12 / 20 / 17 | 4 / 7 / 4 | 12 / 32 / 5 |

Specificity against the report's figure, in datasets with estimates. The report's figure is the lineage-split figure for a learner, and the fixed rule's whole-cohort figure for the known-gene rule. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.

| Species | Set covering machine | Decision tree | Lasso | Known-gene rule |
|---|---|---|---|---|
| *E. cloacae* | 0 / 1 / 1 | 0 / 2 / 0 | 0 / 1 / 1 | 1 / 1 / 0 |
| *E. faecium* | 0 / 2 / 3 | 0 / 4 / 1 | 0 / 1 / 1 | 2 / 2 / 1 |
| *E. coli* | 1 / 9 / 3 | 1 / 9 / 3 | 0 / 1 / 1 | 2 / 10 / 1 |
| *K. pneumoniae* | 1 / 7 / 9 | 3 / 6 / 8 | 0 / 0 / 3 | 2 / 3 / 12 |
| *S. enterica* | 1 / 1 / 6 | 0 / 1 / 7 | 0 / 1 / 1 | 0 / 0 / 8 |
| *S. pneumoniae* | 1 / 2 / 1 | 2 / 1 / 1 | 2 / 1 / 1 | 1 / 3 / 0 |
| All 6 | 4 / 22 / 23 | 6 / 23 / 20 | 2 / 5 / 8 | 8 / 19 / 22 |

The two readings of a row together, over all the species above. A row is one dataset with estimates and one learner. The rows with both readings consistent are part of the rows with neither reading below.

| Learner | Rows read | Neither reading below | Of them, both consistent | Sensitivity below only | Specificity below only | Both below |
|---|---:|---:|---:|---:|---:|---:|
| Set covering machine | 49 | 36 | 17 | 9 | 3 | 1 |
| Decision tree | 49 | 31 | 14 | 12 | 6 | 0 |
| Lasso | 15 | 9 | 1 | 4 | 2 | 0 |
| Known-gene rule | 49 | 30 | 16 | 11 | 7 | 1 |

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

Rows read that sit at an edge, over all the species above. A report figure of exactly 0% or 100% can be read one way only: nothing reads above 100% or below 0%. A row that gives every scored genome the same call has a new estimate of 0% on one measure and 100% on the other. Each "Of them, holds" column counts the rows of the column to its left where the learner holds.

| Learner | Rows read | Report figure at 0% or 100% on a measure | Of them, holds | Same call for every scored genome | Of them, holds |
|---|---:|---:|---:|---:|---:|
| Set covering machine | 49 | 2 | 0 | 1 | 0 |
| Decision tree | 49 | 1 | 0 | 1 | 0 |
| Lasso | 15 | 0 | 0 | 1 | 0 |
| Known-gene rule | 49 | 6 | 4 | 4 | 4 |

**Datasets scored in both populations**

A dataset in both populations has one final model and one report figure, read on two sets of new genomes. A learner holds when neither its sensitivity nor its specificity reads below. No genome counts are added. The plan does not make this comparison; it was added after both results were read.

| Learner | Rows read in both | Holds in both | Holds in the same-study population only | Holds in the primary population only | Holds in neither |
|---|---:|---:|---:|---:|---:|
| Set covering machine | 34 | 10 | 19 | 0 | 5 |
| Decision tree | 34 | 4 | 19 | 1 | 10 |
| Lasso | 8 | 1 | 4 | 0 | 3 |
| Known-gene rule | 34 | 12 | 13 | 3 | 6 |

The same rows by the estimate alone. A small batch has a wide interval, so it reads below less often; this count uses no interval. Each cell counts rows as higher in the same-study population / equal / higher in the primary population.

| Learner | Rows read in both | Sensitivity | Specificity |
|---|---:|---|---|
| Set covering machine | 34 | 27 / 0 / 7 | 29 / 0 / 5 |
| Decision tree | 34 | 25 / 0 / 9 | 30 / 0 / 4 |
| Lasso | 8 | 5 / 0 / 3 | 4 / 0 / 4 |
| Known-gene rule | 34 | 13 / 2 / 19 | 23 / 2 / 9 |

The same-study genomes of the datasets with estimates in both populations, by the study that holds the most resistant genomes there. Each share is the lowest over the study's datasets.

| Study | Datasets | Lowest share of the resistant genomes | Lowest share of the susceptible genomes, for the study that holds the most of them |
|---|---:|---:|---:|
| PRJNA288601 | 23 | 100.0% | 100.0% |
| PRJNA248792 | 7 | 99.9% | 98.4% |
| PRJEB3084 | 4 | 100.0% | 100.0% |
| All 3 | 34 | | |
<!-- tables:end -->

## Revision 14

The site has served evaluation revision 14 since 13 Nov 2026. The runs of record pinned the revision 12 reports. Here the same new estimates are read against the revision 14 figures. Nothing was run again, and each new estimate and interval is the run's.

**Same models.** The release notes of revision 14 say that no fold, final model or known-gene baseline differs from revision 12 in any report. For the final models and the known-gene rule, the repository checks this on the 97 datasets read here. `reference/outside-study-models/` holds one digest of each final model in each of the two releases, for the five learners. The digest is the same in both revisions on 485 of the 485 pairs of a dataset and a learner. The known-gene rule's digest is the same on the 90 datasets whose reports publish a rule, and the other 7 publish none in either revision. The suite fails if one differs. A digest covers every field a report publishes for the model. For the forest and the boosted trees a report publishes their largest importances and no tree, so their digests cover the importances only. The runs scored neither of those two. Each new estimate is therefore an estimate for a model the site serves now. No check in the repository covers the folds. `reference/outside-study-figures/rev14-274-20261112.tsv` holds the lineage-split row of each revision 14 report read here, 97 reports in all. For the set covering machine, the decision tree and the known-gene rule, every row has the counts the runs read.

**The lasso's figure.** Revision 14 scores a null fit, a cross-validation fold where the lasso keeps no weight ([equal-scoring.md](equal-scoring.md)). From revision 7 to revision 13 that fold's genomes had no call and stayed out of the figure. Revision 14 scores each of them at the fold's training share of resistant genomes and calls it at the usual cutoff. The runs scored the lasso on 33 of the 97 datasets, and its lineage-split counts differ on 8 of them.

`scripts/outside-study-results/revision.py` prints these three tables. A row whose counts are the runs' keeps the runs' readings. A row whose counts differ is read with the kit's own rule on the revision 14 figure and the run's interval.

<!-- revision:begin -->
**Datasets where a learner holds, against each report figure**

| Population | Learner | Datasets read | Report figure differs | Holds, runs of record | Holds, revision 14 |
|---|---|---:|---:|---:|---:|
| primary | Set covering machine | 81 | 0 | 32 | 32 |
| primary | Decision tree | 81 | 0 | 23 | 23 |
| primary | Lasso | 25 | 7 | 7 | 8 |
| primary | Known-gene rule | 74 | 0 | 39 | 39 |
| same-study | Set covering machine | 49 | 0 | 36 | 36 |
| same-study | Decision tree | 49 | 0 | 31 | 31 |
| same-study | Lasso | 15 | 3 | 9 | 7 |
| same-study | Known-gene rule | 49 | 0 | 30 | 30 |

**Rows read whose report figure differs**

Each new estimate is followed by its 95% interval and its counts, and each report figure by its counts.

| Population | Species | Drug | Learner | Measure | New estimate | Report, runs of record | Reading | Report, revision 14 | Reading |
|---|---|---|---|---|---|---|---|---|---|
| primary | *K. pneumoniae* | ceftazidime/avibactam | Lasso | sensitivity | 87.1% (70.2% to 96.4%), 27 of 31 | 0.0%, 0 of 7 | above | 0.0%, 0 of 79 | above |
| primary | *K. pneumoniae* | ceftazidime/avibactam | Lasso | specificity | 75.0% (55.1% to 89.3%), 21 of 28 | 89.6%, 69 of 77 | below | 95.6%, 172 of 180 | below |
| primary | *K. pneumoniae* | colistin | Lasso | sensitivity | 33.3% (24.6% to 43.1%), 36 of 108 | 60.5%, 23 of 38 | below | 20.5%, 23 of 112 | above |
| primary | *K. pneumoniae* | colistin | Lasso | specificity | 74.6% (61.6% to 85.0%), 44 of 59 | 94.8%, 110 of 116 | below | 97.3%, 218 of 224 | below |
| primary | *K. pneumoniae* | doripenem | Lasso | sensitivity | 76.6% (68.8% to 83.2%), 111 of 145 | 94.1%, 16 of 17 | below | 99.5%, 216 of 217 | below |
| primary | *K. pneumoniae* | doripenem | Lasso | specificity | 57.7% (50.0% to 65.1%), 101 of 175 | 25.0%, 1 of 4 | above | 3.8%, 1 of 26 | above |
| primary | *S. enterica* | cefazolin | Lasso | sensitivity | 94.3% (89.7% to 97.2%), 164 of 174 | 40.8%, 29 of 71 | above | 64.7%, 77 of 119 | above |
| primary | *S. enterica* | cefazolin | Lasso | specificity | 41.4% (23.5% to 61.1%), 12 of 29 | 61.8%, 21 of 34 | below | 31.8%, 21 of 66 | consistent |
| primary | *S. enterica* | ceftazidime | Lasso | sensitivity | 91.6% (87.3% to 94.8%), 218 of 238 | 25.0%, 3 of 12 | above | 10.0%, 3 of 30 | above |
| primary | *S. enterica* | ceftazidime | Lasso | specificity | 19.1% (14.7% to 24.3%), 53 of 277 | 96.0%, 311 of 324 | below | 98.3%, 740 of 753 | below |
| primary | *S. enterica* | sulfamethoxazole | Lasso | sensitivity | 93.0% (89.3% to 95.7%), 252 of 271 | 75.8%, 47 of 62 | above | 54.3%, 50 of 92 | above |
| primary | *S. enterica* | sulfamethoxazole | Lasso | specificity | 98.3% (94.1% to 99.8%), 118 of 120 | 88.2%, 30 of 34 | above | 37.7%, 69 of 183 | above |
| primary | *S. pneumoniae* | meropenem | Lasso | sensitivity | 97.9% (88.7% to 99.9%), 46 of 47 | 95.2%, 59 of 62 | consistent | 81.9%, 59 of 72 | above |
| primary | *S. pneumoniae* | meropenem | Lasso | specificity | 88.1% (83.9% to 91.5%), 266 of 302 | 98.3%, 1,365 of 1,389 | below | 98.6%, 1,712 of 1,736 | below |
| same-study | *E. cloacae* | ciprofloxacin | Lasso | sensitivity | 64.0% (42.5% to 82.0%), 16 of 25 | 58.1%, 36 of 62 | consistent | 61.2%, 41 of 67 | consistent |
| same-study | *E. cloacae* | ciprofloxacin | Lasso | specificity | 100.0% (92.3% to 100.0%), 46 of 46 | 88.6%, 31 of 35 | above | 72.1%, 31 of 43 | above |
| same-study | *K. pneumoniae* | doripenem | Lasso | sensitivity | 93.1% (89.0% to 96.0%), 216 of 232 | 94.1%, 16 of 17 | consistent | 99.5%, 216 of 217 | below |
| same-study | *K. pneumoniae* | doripenem | Lasso | specificity | 50.0% (43.0% to 57.0%), 104 of 208 | 25.0%, 1 of 4 | above | 3.8%, 1 of 26 | above |
| same-study | *S. enterica* | ceftazidime | Lasso | sensitivity | 94.1% (83.8% to 98.8%), 48 of 51 | 25.0%, 3 of 12 | above | 10.0%, 3 of 30 | above |
| same-study | *S. enterica* | ceftazidime | Lasso | specificity | 95.6% (94.5% to 96.6%), 1,581 of 1,653 | 96.0%, 311 of 324 | consistent | 98.3%, 740 of 753 | below |

**Genomes of null-fit folds in the revision 14 figure**

Each dataset read whose report keeps a second figure, over the fitted folds only. The calls are the ones the genomes of null-fit folds get.

| Species | Drug | Learner | Genomes | In null-fit folds | Share | Called resistant | Called susceptible |
|---|---|---|---:|---:|---:|---:|---:|
| *E. cloacae* | ciprofloxacin | Lasso | 110 | 13 | 11.8% | 13 | 0 |
| *S. pneumoniae* | meropenem | Lasso | 1,808 | 357 | 19.7% | 0 | 357 |
| *S. enterica* | cefazolin | Lasso | 185 | 80 | 43.2% | 80 | 0 |
| *K. pneumoniae* | colistin | Lasso | 336 | 182 | 54.2% | 0 | 182 |
| *S. enterica* | ceftazidime | Lasso | 783 | 447 | 57.1% | 0 | 447 |
| *S. enterica* | sulfamethoxazole | Lasso | 275 | 179 | 65.1% | 113 | 66 |
| *K. pneumoniae* | ceftazidime/avibactam | Lasso | 259 | 175 | 67.6% | 0 | 175 |
| *K. pneumoniae* | doripenem | Lasso | 243 | 222 | 91.4% | 222 | 0 |
<!-- revision:end -->

**Three datasets move the count.** Against revision 14 the lasso holds on 8 of the 25 primary datasets read for it, where the runs read 7. It holds on 7 of the 15 same-study datasets, where the runs read 9.

- *S. enterica* cefazolin, primary: specificity read below and now reads consistent, so the lasso holds there.
- *K. pneumoniae* doripenem, same-study: sensitivity read consistent and now reads below.
- *S. enterica* ceftazidime, same-study: specificity read consistent and now reads below.

On 2 more primary datasets a reading changes and the lasso holds against neither figure: *K. pneumoniae* colistin and *S. pneumoniae* meropenem.

**Null-fit folds.** On these 8 datasets the revision 14 figure also counts the genomes of null-fit folds, and the third table gives their number. It runs from 13 of 110 genomes to 222 of 243, and it is more than half of the cohort on 5 of the 8. Each of those genomes is called from its fold's training share of resistant genomes, with no fitted lasso. On 7 of the 8 datasets they all get one call. On *S. enterica* sulfamethoxazole 113 are called resistant and 66 susceptible.

**One dataset.** For *K. pneumoniae* doripenem the runs' figure counted 21 of the cohort's 243 genomes, the ones a fitted fold called, and 19 of the 21 were called resistant. The revision 14 figure counts all 243. It has 216 of the 217 resistant genomes called resistant and 1 of the 26 susceptible genomes called susceptible. The 222 genomes it adds are all called resistant. So 222 of the 243 calls behind that figure come from a training share, and 21 from a fitted lasso.

**Two figures.** The plan pins the revision 12 reports and reads no other, so it does not say which figure of a revision 14 report to read. A revision 12 report has one lineage-split figure for the lasso, over the genomes of fitted folds. A revision 14 report puts the figure over every genome in its place and keeps the earlier one as `called_metrics`. On all 8 datasets the counts of `called_metrics` are the ones the runs read.

**Choice.** The tables above read the figure over every genome, the field the kit's code takes (`score.reference_row`). This page made that choice after the results were read. So the lasso's count depends on the figure read: 7 or 8 of 25, and 9 or 7 of 15. The first sentence of the result stands against both figures. The second stands for the lasso against the figure over fitted folds only.

**Direction.** The second supported statement stands for the lasso against revision 14 too. In the primary population its sensitivity reads above on 11 datasets and below on 4, and its specificity reads below on 15 and above on 6.

## Species by species

[outside-study-results-by-species.md](outside-study-results-by-species.md) reads each species of each population, with the rows it cites and their figures.

## Limits

- **Same database.** The genomes come from BV-BRC, like the cohort, with labels made under the same practices. The batch's metadata and label totals were read before the plan was frozen.
- **One model, no refit.** Each estimate comes from one final model for each dataset. Its interval covers the sampling of these genomes only. It does not cover the differences between studies.
- **No test.** The reading compares one figure with one interval. There is no hypothesis test and no correction for the number of rows. The datasets of a species share genomes and studies, so their rows are not independent.
- **Few studies.** A study that holds most of a batch decides its estimate. On 20 of the 82 primary datasets, more than half of the resistant genomes come from one study and more than half of the susceptible genomes from another.
- **Batch size.** A row reads below only when the report's figure is above the whole interval. A small batch has a wide interval, so it reads below less often. The same-study batches of *E. cloacae* and *E. coli* have 71 and 87 genomes.
- **Two populations.** The populations differ in species and in datasets. Only the 34 datasets in both compare one model on two sets of genomes, and there the batches still differ in size and in study.
- **Not in the plan.** The plan fixes the reading of each measure on each dataset. The count of datasets where a learner holds, and the comparison of the two populations, were added after the results were read. So was the reading against revision 14, and with it the choice between the two figures of a revision 14 report. No genome is pooled.
- **Two kinds of figure.** The report's figure is for models refitted in each fold. The new figure is for the one final model. For the known-gene rule the report's figure is the fixed rule on the whole cohort.
- **Precision is not read.** It depends on the share of resistant genomes, which differs between the cohort and the new genomes.
- **Three learners.** Random forests and boosted trees are not scored, and the lasso is scored on part of the datasets.
- **No cause.** A k-mer that marks resistance in the cohort can mark a lineage or a plasmid there. The test shows that a figure moved. It does not show why.

The other limits come from the genomes and from the figures themselves.

**The same-study population is four studies.** On each of its 49 datasets one study holds at least 99.9% of the resistant genomes, and 43 of the 49 have one study only. Four studies fill that place: PRJNA288601, PRJNA771404, PRJNA248792 and PRJEB3084. Its result is almost wholly the result on those four.

**Figures at an edge.** A report figure of exactly 0% or 100% can be read one way only. In the primary population 12 of the 74 known-gene rows read have such a figure, and the rule holds on all 12. On 11 of the 74 it gives every scored genome the same call, and it holds on 10 of those.

**Borderline genomes.** The species check marks 332 of the 7,258 primary genomes as borderline, 328 of them *E. coli*. It marks 43 of the 3,144 scored same-study genomes, 25 of them *E. cloacae*. The plan keeps and flags a borderline genome, so they are scored.

**Repeated k-mer sets.** 12 scored primary genomes have the k-mer set of another scored genome: 10 *S. enterica* and 2 *S. pneumoniae*. The plan lists them and does not exclude them, so each counts as a genome.

**Older genome records.** The plan records that 260 of the 7,258 primary genomes have a genome record from before July 2025. There the sequence record is older and the lab record is new.

## Supported statements

It supports four statements.

1. A report's cross-validation figure is an estimate for the cohort's own population. On newly deposited studies it was often outside the new estimate's interval, with the lineage split too.
2. The direction is not fixed. In the primary population, for each learner and for the known-gene rule, sensitivity read above the report's figure more often than below it, and specificity read below more often than above.
3. The known-gene rule's figure moves too. In the primary population it reads below on at least one measure on 35 of its 74 datasets read.
4. On the 34 datasets in both populations, the set covering machine holds on 29 with same-study genomes and on 10 with genomes from other studies.

Three facts qualify them. In the same-study population specificity read above more often than below. The known-gene rule gives every scored genome the same call on 12 of its 75 primary datasets with estimates. The same-study genomes of the 34 datasets come from three studies, and that comparison was added after the results were read.

It does not support a claim that a model works, or fails, in a clinical lab. It names no gene as a cause.

## Files

`reference/outside-study/<population>/<species>/` holds the kit's 11 output files for each run, as the run wrote them. `results.json` is stored gzipped. `reference/outside-study/SHA256SUMS` lists every file.

```sh
python3 -I scripts/outside-study-results/tables.py reference/outside-study
python3 -I scripts/outside-study-results/tables.py reference/outside-study --check docs/outside-study-results.md
python3 -I scripts/outside-study-results/tables.py reference/outside-study --species --check docs/outside-study-results-by-species.md
python3 -I -m unittest discover -s tests/outside_study_results -t . -v
```

The second and third commands exit 1 unless a page's tables are exactly what the files give. The suite runs both in CI. It also rebuilds from the files each sentence of the two pages that carries a count written in digits, and it fails when the sentence is not in the text as written. Any other digit in a page's text fails the suite too, unless the suite lists it as a fixed term. The fixed terms are the level of the interval, the two ends a figure can sit at, the day of the runs, the day the site changed revision, the names of the revisions and a command's exit status.

`reference/outside-study-figures/` holds the figure file of the section Revision 14. The file gives the hash of each report it was read from, and its README names the release and the command that wrote it. A second file there holds the release's own checksum line for each of the 97 reports, and the suite compares the two. The suite checks that section's tables too:

```sh
python3 -I scripts/outside-study-results/revision.py reference/outside-study reference/outside-study-figures/rev14-274-20261112.tsv
python3 -I scripts/outside-study-results/revision.py reference/outside-study reference/outside-study-figures/rev14-274-20261112.tsv --check docs/outside-study-results.md
```

`reference/outside-study-models/` holds the two digest files of the section Revision 14, one for each release. Its README names the releases, the command that wrote the files and what a digest covers. The suite compares the two files. It also compares each file's report hashes with hashes the repository already holds: the kit's pins for revision 12, and the checksum lines above for revision 14. The repository does not hold the reports, so the suite cannot take a digest again from its report. The second command does that where a release is unpacked, with `R` its directory:

```sh
python3 -I scripts/outside-study-results/models.py reference/outside-study-models/rev12-274-20261105.tsv reference/outside-study-models/rev14-274-20261112.tsv
python3 -I scripts/outside-study-results/models.py reference/outside-study-models/rev14-274-20261112.tsv --reports "$R/datasets"
```

Two counts are outside that rule: the 46 files and the 79 FASTA files of the *S. aureus* repeat. They come from the run logs, which this repository does not hold, and the suite lists them as such.

The suite checks each run's record against the frozen population too: the kit commit, the population, the species, and the genome and dataset lists.
