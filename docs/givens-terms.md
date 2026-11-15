# Givens terms

GRM gives [Givens](https://github.com/editnori/givens) the terms below, the assumptions a report figure rests on and their counts. `scripts/givens-export/export.py` writes this page and `reference/givens-export/grm.json` from the same sources, so the two hold the same words and numbers. [givens-export.md](givens-export.md) describes the file.

The page holds 22 terms, 24 links between them and 3 assumptions. Their 36 rows of counts split into 28 read from a file and 8 quoted from a document. In 2 entries GRM says what it has not counted, after the words Not counted.

A change made to this page by hand fails the check. Change `scripts/givens-export/spec.json` or a source, then run the script with `--write`.

## Key

A row gives a count, what it counts, the evaluation revision of the count and its source.

- **File.** The script read the number from a committed summary file. The row names the file, and Sources gives the path for each revision.
- **Document.** The script found the number in a sentence of that document, and the export holds the sentence. No committed file rebuilds the number unless the row names a suite.
- **Median.** A median has no count. Its row gives the value under each revision, in the order of the Revision cell, and the number of splits the median covers.
- **Note.** A letter in brackets points to the note under the table.

A term names its kind in brackets. A state may add a truth state, one of the values Givens keeps apart, such as missing or inferred. A term's links name other terms of this page, and an assumption sits under the terms it bears on.

## Terms

**Dataset** (concept). One species with one antibiotic, and the genomes that have a reported laboratory call for it. Stated in `README.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 274 of 277 | curated datasets that have a report in the release the site serves | 14 | `README.md` |

**Report figure** (measure). A learner's sensitivity, specificity, balanced accuracy, very major error rate or major error rate from cross-validation, as a dataset report states it. It is measured by Split, is part of Dataset, is an input of Reading and is part of Learner row. Stated in `docs/outside-study-results.md` and `README.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 0 of 6,424 | error-rate counts in the release's interval files that differ from the count rebuilt from the released index | 14 | `index-counts.json` |

**Split** (concept). How a dataset's genomes are divided into folds for cross-validation. Each learner is cross-validated on two splits. Stated in `README.md`.

**Random split** (state). Stratified 5-fold cross-validation, with each fold keeping the dataset's mix of resistant and susceptible genomes. It is a value of Split. Stated in `README.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 274 of 274 | datasets whose report has a random split | 14 | `guidance.json` and `summary.json` |

**Lineage split** (state). Folds that hold out whole lineage clusters, so a test genome has no near-identical relative in training. It is the stricter of the two splits. It is a value of Split. Stated in `README.md` and `docs/outside-study-results.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 268 of 274 | datasets whose report has a lineage split | 14 | `summary.json` |

**Duplicate group** (concept). The genomes of one dataset that have the same k-mer set. From revision 10 the random split keeps each group whole in one test fold. It is an input of Random split. Stated in `docs/duplicate-isolation.md`.

**Evaluation revision** (concept). The numbered version of the scoring rules a report was computed under. A figure is read with its split and its revision. It guards Report figure and guards Fold scored at prevalence. Stated in `README.md`.

**Outside-study population** (concept). Newly deposited genomes from the same database, scored by a report's final model, unchanged. It is a retrospective test on newly deposited studies from the same database, and it is no external validation. It is an input of Reading. Stated in `docs/outside-study-results.md`.

**Studies the cohort does not hold** (state). New genomes from studies the cohort does not hold. The plan's first population. It is a value of Outside-study population. Stated in `docs/outside-study-results.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 81 of 82 | datasets of the population that can be read against a report figure | 12 | `outside-study.json` |
| 7,258 of 7,258 | genomes scored, of the genomes the population froze | 12 | `outside-study.json` |
| 7 | species in the population (a) | 12 | `outside-study.json` |

**Note a.** A total with no denominator. The export names no species.

**Studies the cohort already holds** (state). New genomes from studies the cohort already holds. The plan scores this population apart from the first. It is a value of Outside-study population. Stated in `docs/outside-study-results.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 49 of 49 | datasets of the population that can be read against a report figure | 12 | `outside-study.json` |
| 3,144 of 3,155 | genomes scored, of the genomes the population froze | 12 | `outside-study.json` |
| 6 | species in the population (a) | 12 | `outside-study.json` |

**Note a.** A total with no denominator. The export names no species.

**Reading** (concept). How a report figure compares with the estimate on new genomes: below, consistent or above. Sensitivity and specificity are read one at a time. Stated in `docs/outside-study-results.md`.

**Holds** (state, asserted). Neither the learner's sensitivity nor its specificity reads below on the dataset. It is a value of Reading and contrasts with Below. Stated in `docs/outside-study-results.md`.

**Below** (state, asserted). The report's figure is above the whole 95% interval of the new estimate. It is a measured comparison, and figures move both ways. It is a value of Reading. Stated in `docs/outside-study-results.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 36 of 81 | datasets on which the set covering machine's sensitivity reads above the report's figure, on genomes from studies the cohort does not hold | 12 | `outside-study.json` |
| 27 of 81 | datasets on which the set covering machine's sensitivity reads below the report's figure, same population | 12 | `outside-study.json` |

**Not readable** (state, missing). The report has no lineage-split figure to read against. The random split is not put in its place. It is a value of Reading. Stated in `docs/outside-study-test.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 1 of 82 | datasets that cannot be read, on genomes from studies the cohort does not hold | 12 | `outside-study.json` |

**Antibiotic vocabulary** (artifact). The names GRM resolves antibiotics through. GRM's bundled vocabulary works alone, a SigDrug vocabulary is an optional added file, and the report records which one was used. Stated in `README.md`.

**Learner row** (concept). One learner's figures on one split of one dataset. From revision 14 each learner row covers every genome of the split. It is part of Dataset and is an input of Best learner. Stated in `docs/equal-scoring.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 2,710 of 2,710 | learner rows in the release, of the rows that five learners on each split that ran would give | 14 | `index-counts.json`, `summary.json` and `guidance.json` |

**Best learner** (concept). The learner with the best random-split balanced accuracy among learners scored on every genome on every split that ran. A dataset report names it. It is part of Dataset. Stated in `docs/equal-scoring.md`.

**Called metrics** (measure). A learner's figures over the genomes of fitted folds only. From revision 14 a row with a null-fit fold keeps this earlier figure beside the figure over every genome. It is part of Learner row. Stated in `docs/equal-scoring.md`.

**Fold scored at prevalence** (state, inferred). A fold whose fit learned nothing beyond the class balance. Revision 14 scores each of its test genomes at the fold's training prevalence. It is part of Learner row. Stated in `docs/equal-scoring.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 122 of 2,710 | learner rows that carry the flag for such a fold | 14 | `index-counts.json` |

**No report** (state, missing). A curated dataset with no report, because one class is under the class minimum. The release lists each one with its reason. It is part of Dataset. Stated in `README.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 3 of 277 | curated datasets with no report in the release the site serves | 14 | `README.md` |

**Split unavailable** (state, missing). A dataset whose lineage split could not be formed. The report records why. It is a value of Lineage split. Stated in `README.md`.

| Count | Counted | Revision | Source |
|---|---|---|---|
| 6 of 274 | datasets with no lineage split | 14 | `summary.json` |

**Label inferred from a MIC and left out** (state, inferred). A label inferred from a measured concentration. It stays in the curated data, marked, and no learner trains or is scored on it. It is part of Dataset. Stated in `README.md`.

**Not counted.** No committed file counts the datasets that hold such labels. The released index has the field, and the repository keeps no count of it.

## Assumptions

**Assumption.** A genome in a test fold has no copy in the training folds.

| Part | Text |
|---|---|
| Hidden in | Stratified random folds drawn without looking for identical genomes. |
| If false | A model is scored on a genome it trained on, and the figure reads better than the model is. |
| Under | Random split, Duplicate group and Report figure |
| Stated in | `docs/duplicate-isolation.md` |

| Count | Counted | Revision | Source |
|---|---|---|---|
| 58 of 277 | datasets with at least one duplicate group, on the cohorts of revision 9 | 9 | `docs/duplicate-isolation.md` |
| 55 of 277 | datasets with a group that crossed the saved random folds, on the cohorts of revision 9, before groups were kept whole | 9 | `docs/duplicate-isolation.md` |

**Not counted.** No committed file counts duplicate groups on the cohorts of revision 10 or later. The two counts are of the revision 9 cohorts, and the document says to count again after the revision 10 curation. How far a figure moved is given for named species only, so it is left out.

**Assumption.** A cross-validated figure holds for genomes from a newly deposited study the cohort never saw.

| Part | Text |
|---|---|
| Hidden in | Reading a report's sensitivity as what the model will do next. |
| If false | The figure is trusted on new genomes where it reads below more often than it holds. |
| Under | Report figure, Lineage split, Studies the cohort does not hold and Reading |
| Stated in | `docs/outside-study-results.md` |

| Count | Counted | Revision | Source |
|---|---|---|---|
| 32 of 81 | datasets on which the set covering machine holds, on genomes from studies the cohort does not hold | 12 | `outside-study.json` |
| 23 of 81 | datasets on which the decision tree holds, same population | 12 | `outside-study.json` |
| 7 of 25 | datasets on which the lasso holds, of the datasets it can be scored and read on, same population | 12 | `outside-study.json` |
| 8 of 25 | datasets on which the lasso holds against the revision 14 figure, same population and same new estimates | 14 | `docs/outside-study-results.md`, rebuilt by `tests/outside_study_results` |
| 36 of 49 | datasets on which the set covering machine holds, on new genomes from studies the cohort already holds (a) | 12 | `outside-study.json` |
| 31 of 49 | datasets on which the decision tree holds, on new genomes from studies the cohort already holds (a) | 12 | `outside-study.json` |
| 9 of 15 | datasets on which the lasso holds, of the datasets it can be scored and read on, on new genomes from studies the cohort already holds (a) | 12 | `outside-study.json` |
| 7 of 15 | datasets on which the lasso holds against the revision 14 figure, on new genomes from studies the cohort already holds (b) | 14 | `docs/outside-study-results.md`, rebuilt by `tests/outside_study_results` |

**Note a.** A separate population. The plan scores the two apart, and the document set them side by side after both results were read.

**Note b.** A separate population, and the same new estimates as the revision 12 line.

**Assumption.** Every learner row of a dataset is scored on the same genomes.

| Part | Text |
|---|---|
| Hidden in | A table with one balanced-accuracy column for each learner, and the best-learner label. |
| If false | A row compares figures over different genomes, and the best learner is chosen among fewer learners. |
| Under | Learner row, Best learner, Report figure, Fold scored at prevalence and Called metrics |
| Stated in | `docs/equal-scoring.md` |

| Count | Counted | Revision | Source |
|---|---|---|---|
| 122 of 2,710 | learner rows that left genomes without a call, in the revision 12 release | 12 | `index-counts.json` |
| 65 of 122 | of those rows, the boosted trees' rows on either split | 12 | `guidance.json` and `index-counts.json` |
| 57 of 122 | of those rows, the lasso's rows on either split | 12 | `guidance.json` and `index-counts.json` |
| 122 of 2,710 | learner rows scored on every genome under the flag, in the revision 14 release | 14 | `index-counts.json` |
| 0 of 2,710 | learner rows that still count fewer genomes than their class, in the revision 14 release | 14 | `index-counts.json` |
| 73 of 274 | datasets whose revision 14 report differs from revision 12 in the rows of a learner with a null fit: 48 there only, and 25 there and in panel rows | 14 | `README.md` |
| 21 of 274 | datasets whose best learner changed from revision 12 to revision 14 | 14 | `README.md` |
| 0.8459 and 0.8409, over 268 | lineage splits over which the boosted trees' median balanced accuracy is taken, under each revision (a) | 12 and 14 | `summary.json` |
| 0.8567 and 0.8565, over 268 | lineage splits over which the lasso's median balanced accuracy is taken, under each revision (a) | 12 and 14 | `summary.json` |

**Note a.** A median, so the line has no count.

## Rules

The export carries these rules for whoever shows its entries. In them, the file is the export.

- GRM writes this file with scripts/givens-export/export.py. Givens stores it unchanged and names the GRM commit it copied it from.
- A measured line is an aggregate, and its text says what it counts. The file names no dataset, species or antibiotic.
- A line with basis file is read from a committed summary file, and the export lists that file's sha256.
- A line with basis document is quoted from the repository's text. No committed file rebuilds it unless the line names a suite.
- An entry with unmeasured says what GRM has not counted. Givens shows that as missing and fills in no number.
- A figure is read with its split, its learner and its evaluation revision. A line of revision 12 is never shown as revision 14's.
- GRM reads nothing from Givens.

## Left out

The export and this page never hold:

- Dataset ids, species names and antibiotic names.
- Rows for one species or one drug class.
- GRM's features, models, panels and samples.
- The known-gene rule's rows and its catalog.
- MIC values and the MIC regression's figures.
- Sequences, k-mers, variants and genes.
- The error-limit screen's counts.
- How far a figure moved where the document gives it for named species only.

## Sources

A row names a summary file by its file name. It means the file of the row's revision.

| File | Revision | Path |
|---|---|---|
| `guidance.json` | 12 | `reference/release-summary/rev12-274-20261105/guidance.json` |
| `index-counts.json` | 12 | `reference/release-summary/rev12-274-20261105/index-counts.json` |
| `summary.json` | 12 | `reference/release-summary/rev12-274-20261105/summary.json` |
| `guidance.json` | 14 | `reference/release-summary/rev14-274-20261112/guidance.json` |
| `index-counts.json` | 14 | `reference/release-summary/rev14-274-20261112/index-counts.json` |
| `summary.json` | 14 | `reference/release-summary/rev14-274-20261112/summary.json` |
| `outside-study.json` | 12 | `web/src/data/outside-study.json` |

The export lists each file's hash. It also holds the sentence each term and each quoted number rests on. The quoted documents are `README.md`, `docs/duplicate-isolation.md`, `docs/equal-scoring.md`, `docs/outside-study-results.md` and `docs/outside-study-test.md`.
