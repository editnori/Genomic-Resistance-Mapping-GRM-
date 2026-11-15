# Other database

The released reports of `rev14-274-20261112` (evaluation revision 14) were scored on 13,377 NCBI Pathogen Detection genomes with antibiogram records: 3,897 *E. coli* genomes in a pilot and 9,480 in 8 more species. SCM, CART and the lasso are scored; RF and GBT are not, as the release contains their importances only.

It is a test on data from another public database that shares isolates with the training database. The genomes are outside BV-BRC's labelled set by BioSample and assembly accession, yet in *S. enterica* and *C. coli* most scored genomes come from a BioProject that also contains training genomes. It is not external validation.

## Plan

The test ran in two parts, `outside` for the *E. coli* pilot and `outside2` for eight more species, and a review recomputed each part from its saved files. Neither the parts' output directories nor the reviews' reports are public; the record contains every figure this page takes from them, and `--check` needs those directories.

The pilot's plan, `PLAN.md`, was fixed at 2026-11-15T12:30:24Z, before the pilot ran. The eight species' addenda were fixed after it: `ADDENDUM.md` at 2026-11-15T16:04:00Z and `ADDENDUM2.md`, which adds readings L, C and D, at 2026-11-15T17:02:14Z. Each time is the one written in the sha256 file beside the plan file, and `plans.tsv` in the record gives each file's sha256.

The review of the eight species assigned each of the 199 label requests (193 answered with status 200) to a species by the BioSamples it names. In all 8 species the BioSamples asked equal the scored BioSamples, and the first label request came at or after the time the species' predictions were written, and after `ADDENDUM2.md` was fixed.

The pilot's 78 label requests ran from 14:30:42Z to 14:32:37Z, at or after its predictions were written (14:30:42Z). Before `PLAN.md` was fixed, 6 requests of the pilot went to NCBI Pathogen Detection's search service between 12:04:27Z and 12:05:23Z with no parameters saved, and answered 67, 304, 67, 271, 285 and 64 bytes. Whether any carried a phenotype cannot be shown.

In 6 of the 9 species the first label request came in the same second as the predictions were written, so the clock alone cannot order the two. Each run's log shows the label step starting after the prediction step ended.

The pilot's readings other than the plan's rule are exploratory: the pilot's review chose them after its results were seen. For the eight species, `ADDENDUM2.md` fixed readings L, C and D before any of their labels was read.

## Genomes

Each species runs through the same funnel, from NCBI AST isolates (BioSamples) to scored genomes. The first table gives the isolates each step removes and the genomes scored; steps D1, D2, D4, F1 and Q2 removed none in any species.

C2 removes isolates with no release drug, C3 those with no GenBank assembly, C4 those in BV-BRC's labelled set by BioSample or assembly accession, O1 the outside-study candidates, P2 those past the per-species cap, Q1 those failing QC, and D3 those at 0.0005 or less from a training genome on the nearest-training distance.

| Species | Isolates | C2 | C3 | C4 | O1 | P2 | Q1 | D3 | Genomes |
|---|---|---|---|---|---|---|---|---|---|
| *E. coli* | 10,986 | 0 | 215 | 393 |  | 6,378 | 0 | 103 | 3,897 |
| *S. enterica* | 9,976 | 0 | 102 | 4,635 | 0 | 1,239 | 0 | 1,623 | 2,377 |
| *K. pneumoniae* | 2,412 | 9 | 228 | 581 | 0 | 0 | 0 | 59 | 1,535 |
| *S. aureus* | 1,522 | 3 | 103 | 41 | 0 | 0 | 0 | 276 | 1,099 |
| *C. jejuni* | 2,995 | 0 | 11 | 473 | 0 | 0 | 0 | 80 | 2,431 |
| *C. coli* | 1,426 | 0 | 7 | 292 | 0 | 0 | 0 | 125 | 1,002 |
| *N. gonorrhoeae* | 794 | 0 | 33 | 0 | 441 | 0 | 0 | 77 | 243 |
| *P. aeruginosa* | 979 | 1 | 132 | 157 | 0 | 0 | 1 | 162 | 526 |
| *A. baumannii* | 1,142 | 0 | 44 | 702 | 0 | 0 | 0 | 129 | 267 |

Share is the largest BioProject's part of the scored genomes and Training the training genomes it contains. Overlap counts scored genomes in a BioProject that also contains a training genome, Neighbours those within 0.001 of a training genome, and Accession gap the share of BV-BRC's labelled genomes with neither a BioSample nor an assembly accession. Each overlap count is a lower bound.

| Species | Genomes | BioProject | Share | Training | Overlap | Neighbours | Accession gap |
|---|---|---|---|---|---|---|---|
| *E. coli* | 3,897 | PRJNA292663 | 0.266 | 0 | 179 | 316 | 0.347 |
| *S. enterica* | 2,377 | PRJNA292661 | 0.744 | 191 | 1,768 | 1,017 | 0.109 |
| *K. pneumoniae* | 1,535 | PRJNA1322038 | 0.154 | 0 | 127 | 148 | 0.061 |
| *S. aureus* | 1,099 | PRJNA1470711 | 0.264 | 0 | 260 | 607 | 0.239 |
| *C. jejuni* | 2,431 | PRJNA292664 | 0.789 | 0 | 0 | 172 | 0.000 |
| *C. coli* | 1,002 | PRJNA292664 | 0.946 | 1 | 948 | 152 | 0.000 |
| *N. gonorrhoeae* | 243 | PRJNA1217803 | 0.543 | 0 | 0 | 170 | 0.557 |
| *P. aeruginosa* | 526 | PRJNA940483 | 0.253 | 0 | 9 | 116 | 0.026 |
| *A. baumannii* | 267 | PRJNA1322038 | 0.296 | 0 | 5 | 55 | 0.341 |

In *S. enterica*, 1,768 of 2,377 scored genomes are in PRJNA292661, which also contains 191 of the 2,000 training genomes. In *C. coli*, 948 of 1,002 are in PRJNA292664, which contains 1 of its 391 training genomes. So in these two species the scored genomes are different isolates, most of them from BioProjects that also contain training genomes.

In *C. jejuni*, 1,918 of 2,431 scored genomes are in the same PRJNA292664, which contains 0 of its 508 training genomes. So the difference between its overlap count and that of *C. coli* turns on a single training genome.

## Rule

A dataset is estimable for a learner when it has at least 25 resistant and 25 susceptible labelled genomes with a call, and the release has a lineage figure for it. It holds when recall is not below and log loss is not worse than the release's lineage split: recall by a 95% Clopper-Pearson interval, log loss by a 95% bootstrap of 2,000 genome resamples with seed 1.

The plan's rule asks for 25 of each class; the condition that the release has a lineage figure was added after the *S. aureus* labels were read, when 4 of its reports were found to have none, and it moves the lasso in *S. aureus* from meeting the rule to not assessable.

A learner meets the plan's rule in a species when it holds on more than half of at least 3 estimable datasets; with fewer it is not assessable. The resistant and susceptible calls are an adapted research screen: SCM's p is its call, 1 or 0, and CART and the lasso call a genome resistant at p 0.5 or more. Log loss clips p at 0.0001, so SCM's log loss is its error rate times 9.21.

The rule reads recall and log loss only. Specificity is not part of it, so a dataset can hold while its learner calls every susceptible genome resistant.

## Results

The review of the eight species gives its reading in four sentences, kept whole here.

Read them with three limits given in full below: the rule does not read specificity, intervals over genomes treat genomes as independent although genomes cluster by BioProject, and log loss depends on the class mix.

This is a test on data from another public database: NCBI Pathogen Detection genomes outside BV-BRC's labelled set by BioSample and assembly accession, scored with the released reports of `rev14-274-20261112`. In *S. enterica* (0.744) and *C. coli* (0.946), most scored genomes come from a BioProject that also holds training genomes. It is not external validation.

SCM meets the plan's rule in 6 of 8 species. It fails in *C. coli* (2 of 7 datasets hold) and *N. gonorrhoeae* (0 of 3). Over its 84 estimable datasets, recall is below the lineage figure in 22 and log loss is worse in 23. In *K. pneumoniae* the result rests on one dataset (14 of 26 hold), polymyxin B, where SCM calls every susceptible genome resistant.

CART meets the rule in 5 of 8 species. It fails in *S. enterica* (6 of 12), *K. pneumoniae* (12 of 26) and *N. gonorrhoeae* (0 of 3). Over 84 datasets, recall is below in 19 and log loss is worse in 26.

The lasso meets the rule in 3 of the 4 species where it can be read: *K. pneumoniae* 4 of 6 (one of them polymyxin B, specificity 0.0183), *C. jejuni* 3 of 3 and *P. aeruginosa* 3 of 3. It fails in *C. coli* (2 of 6). It is not assessable in 4 species. Over 24 datasets, recall is below in 5 and log loss is worse in 3. RF and GBT are not scored.

Each table below gives one learner alone, with recall, log loss and specificity side by side. Datasets counts the estimable datasets and Hold count those that hold. Recall counts the datasets with recall below the lineage figure and Log loss those with log loss worse. Specificity gives the lowest specificity with its antibiotic, Rule the result under the plan's rule, Margin the datasets that must change for the other result, and Readings the results under L, C and D.

**SCM.** Over its 84 estimable datasets in the eight species it holds on 52, recall is below the lineage figure in 22, log loss is worse in 23, and specificity is under 0.5 in 5. It meets the plan's rule in 6 of 8 species.

| Species | Learner | Datasets | Hold count | Recall | Log loss | Specificity | Rule | Margin | Readings |
|---|---|---|---|---|---|---|---|---|---|
| *S. enterica* | SCM | 12 | 10 | 1 | 1 | 0.3783 (nalidixic acid) | meets | 4 | L meets, C meets, D meets |
| *K. pneumoniae* | SCM | 26 | 14 | 10 | 10 | 0.0000 (cefpodoxime) | meets | 1 | L meets, C meets, D meets |
| *S. aureus* | SCM | 9 | 7 | 2 | 0 | 0.7714 (penicillin) | meets | 3 | L meets, C meets, D meets |
| *C. jejuni* | SCM | 3 | 3 | 0 | 0 | 0.9970 (tetracycline) | meets | 2 | L meets, C meets, D meets |
| *C. coli* | SCM | 7 | 2 | 3 | 3 | 0.9666 (clindamycin) | fails | 2 | L meets, C meets, D meets |
| *N. gonorrhoeae* | SCM | 3 | 0 | 0 | 3 | 0.0000 (tetracycline) | fails | 2 | L fails, C fails, D meets |
| *P. aeruginosa* | SCM | 10 | 8 | 1 | 1 | 0.3333 (cefepime) | meets | 3 | L meets, C meets, D meets |
| *A. baumannii* | SCM | 14 | 8 | 5 | 5 | 0.5714 (sulfamethoxazole/trimethoprim) | meets | 1 | L meets, C meets, D meets |

**CART.** Over its 84 estimable datasets in the eight species it holds on 50, recall is below the lineage figure in 19, log loss is worse in 26, and specificity is under 0.5 in 5. It meets the plan's rule in 5 of 8 species.

| Species | Learner | Datasets | Hold count | Recall | Log loss | Specificity | Rule | Margin | Readings |
|---|---|---|---|---|---|---|---|---|---|
| *S. enterica* | CART | 12 | 6 | 3 | 3 | 0.7237 (nalidixic acid) | fails | 1 | L meets, C fails, D meets |
| *K. pneumoniae* | CART | 26 | 12 | 7 | 13 | 0.0000 (cefpodoxime) | fails | 2 | L meets, C fails, D meets |
| *S. aureus* | CART | 9 | 9 | 0 | 0 | 0.8000 (penicillin) | meets | 5 | L meets, C meets, D meets |
| *C. jejuni* | CART | 3 | 3 | 0 | 0 | 0.9970 (tetracycline) | meets | 2 | L meets, C meets, D meets |
| *C. coli* | CART | 7 | 4 | 3 | 0 | 0.9720 (clindamycin) | meets | 1 | L meets, C meets, D meets |
| *N. gonorrhoeae* | CART | 3 | 0 | 0 | 3 | 0.4906 (tetracycline) | fails | 2 | L fails, C fails, D fails |
| *P. aeruginosa* | CART | 10 | 7 | 2 | 3 | 0.2963 (imipenem) | meets | 2 | L meets, C meets, D meets |
| *A. baumannii* | CART | 14 | 9 | 4 | 4 | 0.6790 (cefepime) | meets | 2 | L meets, C meets, D meets |

**Lasso.** Over its 24 estimable datasets in the eight species it holds on 17, recall is below the lineage figure in 5, log loss is worse in 3, and specificity is under 0.5 in 1. It meets the plan's rule in 3 of the 4 species where it can be read.

| Species | Learner | Datasets | Hold count | Recall | Log loss | Specificity | Rule | Margin | Readings |
|---|---|---|---|---|---|---|---|---|---|
| *S. enterica* | Lasso | 2 | 2 | 0 | 0 | 0.9995 (amoxicillin/clavulanate) | not assessable |  | L not assessable, C not assessable, D not assessable |
| *K. pneumoniae* | Lasso | 6 | 4 | 1 | 2 | 0.0183 (polymyxin b) | meets | 1 | L meets, C meets, D meets |
| *S. aureus* | Lasso | 2 | 2 | 0 | 0 | 0.9797 (gentamicin) | not assessable |  | L not assessable, C not assessable, D not assessable |
| *C. jejuni* | Lasso | 3 | 3 | 0 | 0 | 0.9970 (tetracycline) | meets | 2 | L meets, C meets, D meets |
| *C. coli* | Lasso | 6 | 2 | 3 | 1 | 0.9718 (erythromycin) | fails | 2 | L meets, C fails, D meets |
| *N. gonorrhoeae* | Lasso | 0 | 0 | 0 | 0 |  | not assessable |  | L not assessable, C not assessable, D not assessable |
| *P. aeruginosa* | Lasso | 3 | 3 | 0 | 0 | 0.5312 (piperacillin/tazobactam) | meets | 2 | L meets, C meets, D meets |
| *A. baumannii* | Lasso | 2 | 1 | 1 | 0 | 0.9759 (ceftazidime) | not assessable |  | L not assessable, C not assessable, D not assessable |

In *K. pneumoniae*, polymyxin B holds for all three learners with specificity 0.0000 for SCM, 0.0459 for CART and 0.0183 for the lasso, against 0.0192, 0.0577 and 0.4615 in the lineage split. SCM calls all 109 susceptible genomes resistant, and the dataset holds because its log loss, 3.4148, is under the lineage figure, 5.0657. With SCM on 14 of 26 and the lasso on 4 of 6, polymyxin B decides both results.

*P. aeruginosa* cefepime also holds with specificity 0.3333 for SCM and CART, against 0.7317 and 0.6016 in the lineage split, and decides neither result (margins 3 and 2).

5 of the 20 results under the rule are one dataset from the other result: SCM in *K. pneumoniae* (14 of 26) and *A. baumannii* (8 of 14), CART in *S. enterica* (6 of 12, fails) and *C. coli* (4 of 7), and the lasso in *K. pneumoniae* (4 of 6).

At seeds 2 to 6 in place of the plan's seed 1, one SCM dataset, *P. aeruginosa* ceftazidime, holds at 4 of the 5, and no species result changes.

## Clear failures

Some datasets fail clearly, as the review of the eight species lists them. Each cell gives the figure on these genomes, with the lineage figure in brackets.

| Dataset | Learner | Recall | Specificity | Log loss |
|---|---|---|---|---|
| *S. enterica* ciprofloxacin | SCM | 0.0962 (0.8272) | 0.9995 (0.8930) | 0.3782 (1.1276) |
| *S. enterica* ciprofloxacin | CART | 0.0288 (0.8215) | 0.9787 (0.9417) | 0.3192 (0.3311) |
| *K. pneumoniae* ceftazidime/avibactam | SCM | 0.2951 (0.8987) | 0.7109 (1.0000) | 3.8986 (0.2846) |
| *K. pneumoniae* ceftazidime/avibactam | CART | 0.2787 (0.9241) | 0.7266 (0.8778) | 3.1455 (0.5049) |
| *N. gonorrhoeae* tetracycline | SCM | 1.0000 (0.9248) | 0.0000 (0.9269) | 4.0680 (0.6786) |
| *K. pneumoniae* cefpodoxime | SCM | 0.9636 (0.8944) | 0.0000 (0.9268) | 3.7256 (0.9170) |
| *S. enterica* nalidixic acid | SCM | 0.8810 (0.8560) | 0.3783 (0.9454) | 5.6397 (0.7604) |
| *A. baumannii* ampicillin/sulbactam | SCM | 0.4255 (0.7423) | 0.9877 (0.8664) | 2.0148 (1.8719) |
| *A. baumannii* ampicillin/sulbactam | CART | 0.4468 (0.8090) | 0.9630 (0.6843) | 0.9238 (1.7110) |

## Readings

Reading L keeps the lineage figure with its own interval, C uses a class-balanced log loss with a stratified bootstrap, and D resamples whole BioProjects. Each keeps the rest of the rule.

Under the plan's rule and readings L, C and D, SCM meets the rule in 6, 7, 7 and 8 of 8 species, CART in 5, 7, 5 and 7 of 8, and the lasso in 3, 4, 3 and 4 of the 4 species where it can be read.

No difference between two learners holds in all four readings. Of the 27 comparisons, 24 of a pair in one species and 3 of a pair counted over the species, 0 put the same learner ahead under the plan's rule and all three readings.

Reading D cannot show a difference where a dataset contains few BioProjects. The median count of BioProjects over SCM's estimable datasets is 16.5 in *S. enterica*, 16.5 in *K. pneumoniae*, 10 in *S. aureus*, 5 in *C. jejuni*, 2 in *C. coli*, 7 in *N. gonorrhoeae*, 15 in *P. aeruginosa* and 11.5 in *A. baumannii*.

In *C. coli*, D turns SCM's 2 of 7 into 6 of 7, and in 6 of its 7 datasets the cluster interval of log loss reaches 0.0001, the value of a resample with no error.

## Class mix

The resistant share of the scored genomes differs from the training cohort's in most datasets. Per learner, the table counts the estimable datasets whose share is below or above the cohort's.

| Species | SCM | CART | Lasso |
|---|---|---|---|
| *S. enterica* | 12 below, 0 above | 12 below, 0 above | 2 below, 0 above |
| *K. pneumoniae* | 13 below, 13 above | 13 below, 13 above | 3 below, 3 above |
| *S. aureus* | 8 below, 1 above | 8 below, 1 above | 2 below, 0 above |
| *C. jejuni* | 1 below, 2 above | 1 below, 2 above | 1 below, 2 above |
| *C. coli* | 7 below, 0 above | 7 below, 0 above | 6 below, 0 above |
| *N. gonorrhoeae* | 1 below, 2 above | 1 below, 2 above |  |
| *P. aeruginosa* | 1 below, 9 above | 1 below, 9 above | 0 below, 3 above |
| *A. baumannii* | 12 below, 2 above | 12 below, 2 above | 2 below, 0 above |

SCM's log loss counts its errors, which favours a species with few resistant genomes. In *S. enterica* the resistant share is below the cohort's in 12 of 12 SCM datasets; under reading C, which balances the classes, SCM still holds on 10 of 12.

In *N. gonorrhoeae* SCM fails. In tetracycline and ciprofloxacin its errors are false resistant calls, which the mix does not explain. In tetracycline SCM calls all 53 susceptible genomes resistant: specificity 0.0000, against 0.9269 in the lineage split, and log loss 4.0680. In ciprofloxacin 9 of its 12 errors are false resistant calls.

In azithromycin 44 of its 54 errors are false susceptible calls, with recall 0.5368 against 0.4240 in the lineage split. The resistant share is 0.531 against 0.215 in azithromycin, 0.542 against 0.622 in ciprofloxacin and 0.558 against 0.262 in tetracycline.

In azithromycin the mix does explain the worse log loss. At the lineage split's recall, 0.4240, and specificity, 0.9732, a resistant share of 0.531 gives an error rate of 0.318 and a log loss of 2.93, above the 2.7786 measured. Under reading C, which balances the classes, SCM holds on azithromycin.

## Pilot

The pilot's review gives its reading in one sentence, here in two paragraphs. Each reading of the pilot other than the plan's rule is exploratory.

On 3,897 *E. coli* genomes from NCBI Pathogen Detection whose BioSample and assembly accession match no genome of BV-BRC's labelled set, SCM recall was not below its lineage-split figure on 15 of 25 estimable datasets and CART recall on 16 of 25; under the plan's fixed rule SCM holds on 14 and CART on 12, one dataset either side of the line.

This margin changes with how the lineage figure's own uncertainty and the clustering of genomes by BioProject are handled, so the pilot does not show that the two learners differ. The pilot has fewer resistant genomes than the training cohorts (median share 0.25 against 0.50), which lowers log loss, and SCM's log loss is an error rate times 9.21. The lasso was scored on 4 estimable datasets, all sparse fits.

The table gives the datasets that hold for each learner under the plan's rule and under the four other readings of the pilot's review.

| Reading | SCM | CART | Lasso |
|---|---|---|---|
| Plan's rule | 14 of 25 | 12 of 25 | 3 of 4 |
| L: lineage interval | 17 of 25 | 15 of 25 | 4 of 4 |
| T: two-sample interval | 15 of 25 | 13 of 25 | 3 of 4 |
| C: balanced log loss | 14 of 25 | 11 of 25 | 1 of 4 |
| D: BioProject resamples | 21 of 25 | 18 of 25 | 4 of 4 |

SCM's estimable pilot datasets contain a median of 15 BioProjects, at least 2, and the pilot's resistant share is below the cohort's in 21 of the 25.

## Limits

These limits stand beside every result:

- It is a test on data from another public database that shares isolates with the training database. It is not external validation, and it is not the outside-study test.
- The rule does not read specificity: in the eight species 5 dataset results hold with specificity under 0.5, and in the pilot 4 more hold with specificity under 0.5 where the lineage split is above it.
- Every overlap count is a lower bound: 0.000 to 0.557 of labelled BV-BRC genomes have no accession, and isolates BV-BRC took in after the data root's copy are not excluded.
- In *S. enterica* and *C. coli* most scored genomes share a BioProject with training genomes.
- The class mix differs from the training cohorts, and log loss depends on it.
- The calls are an adapted research screen.
- RF and GBT are not scored, and the lasso cannot be read in 4 species.
- 4 *S. aureus* reports have no lineage figure in the release, so no learner is read on them.
- Reading D shows little where a dataset contains few BioProjects.
- Intervals over genomes treat genomes as independent, and genomes cluster by BioProject.

The two reviews could not check:

- The bytes of the requests rebuilt from the eight species' scripts as they stand, and what their answers contained.
- Whether the saved label answers are complete against NCBI as it stands now: no label was fetched again.
- Isolates BV-BRC contains under no accession or another BioSample, or took in after the data root's copy.
- The parameters and answers of the pilot's requests before `PLAN.md` that the log records without parameters.
- Steps Q2 and D3 for every genome: the reviews recounted them from the two parts' per-genome columns.
- The times the plan files were fixed, beyond the file times of one machine.
- Breakpoint versions: NCBI records no standard year in the fields the two parts fetched.

## Files

- `reference/other-database/` contains the record: the rows per species, per dataset and learner, and per species and learner, the results at other bootstrap seeds, the scored genomes as accessions only, the plan files with their sha256, the settings, and the sha256 of every input.
- `scripts/other-database/record.py` writes the record from the two parts' saved predictions, labels, funnel files and request logs, the release's reports and the data root's drug names. It fits nothing and downloads nothing.
- `tests/other_database` checks this page against the record: it rebuilds each line of the page that carries a number, and keeps the opening's second paragraph, the limits and the definitions of readings L, C and D word for word.

The record contains no label and no prediction per genome. A reader can fetch the scored genomes by the accessions in `genomes.tsv` and redo the overlap, and cannot recompute a cell.

`LANES` is the directory that contains the two parts' output directories, `outside/` and `outside2/`, `R` the unpacked release and `ROOT` the public data root. `--check` succeeds only when the record is what the script writes from them. Neither output directory is public, so only a machine that has them can run the first two commands.

```sh
python3 -I scripts/other-database/record.py --lanes "$LANES" --release "$R" --data-root "$ROOT" --out reference/other-database
python3 -I scripts/other-database/record.py --lanes "$LANES" --release "$R" --data-root "$ROOT" --check reference/other-database
python3 -I -m unittest discover -s tests/other_database -t .
```
