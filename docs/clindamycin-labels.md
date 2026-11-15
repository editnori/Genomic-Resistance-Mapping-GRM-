# Clindamycin figures and the study behind their labels

The release `rev14-274-20261112` (evaluation revision 14) has six clindamycin datasets. In each of them, one deposited study, named here by its PMID, holds most or a large share of the samples on which GRM's known-gene rule and the reported phenotype disagree. What that means differs by dataset, so this page gives each dataset its kind:

- **Measured dependence:** `staphylococcus_aureus__clindamycin` and `campylobacter_coli__clindamycin`. Their figures move when the study's samples leave the cohort.
- **One study:** `enterococcus_faecium__clindamycin`. Every sample comes from one study.
- **Class balance, dependence not measured:** `campylobacter_jejuni__clindamycin` and `streptococcus_agalactiae__clindamycin`. One study holds so much of one class that the dataset cannot be scored without it.
- **No dependence under the rule:** `streptococcus_pneumoniae__clindamycin`.

Under the plan fixed before any count was read, the caveat rule returns a caveat for five of the six datasets. The rule does not tell the three kinds of caveat apart.

The files this page quotes are in [`reference/clindamycin-labels/rev14-20261115/`](../reference/clindamycin-labels/rev14-20261115/README.md), and a path such as `results/chosen.tsv` or `PLAN.md` below is a file of that directory. A path under `review/` is a file of the independent review of the check, made on the same day: every count, refit score and rule verdict of the check reproduced there, and the review narrowed parts of the reading. Every figure is a retrospective reading of public BV-BRC records under internal cross-validation.

## The question

Does one deposited study hold a large share of the disagreements between the reported clindamycin phenotype and GRM's known-gene rule, and do the released figures depend on that study?

The known-gene rule is GRM's own reader of the AMRFinderPlus catalog. It calls a genome resistant from the catalog genes GRM finds in it, and it is fitted on nothing. A disagreement is a sample on which the rule's call and the reported phenotype differ. A disagreement is a difference between two readings, and neither is called the error.

This is a check of labels and scores. It names no gene and ranks nothing.

## The plan and the release read

| | |
|---|---|
| Plan | `PLAN.md`, sha256 `b14ec01aaef5471aa1b8f0e3cf53fd7a3a9485db0dca8a6f5c7d9d8cffd7915b`, recorded 2026-11-14T18:11:52Z |
| Addendum | `PLAN-ADDENDUM-20261115.md`, sha256 `85deaefb48bf88c466f8348363d5bf0a0ef6b0ce8b02256488ed96c6377f2afb`, recorded 2026-11-15T11:36:42Z |
| Release read | `rev14-274-20261112`, evaluation revision 14: the samples of its six clindamycin reports |
| Data | the public data root's `curated/phenotypes.jsonl` and the raw BV-BRC records of `bvbrc/genome_amr.jsonl` |
| Program | `grm` built at `852ebc9` with Rust 1.98.1, binary sha256 `b20b5f1415a011293a63f9fc8b90ef54b3f0935952cbcfc67aa0a36c5bcc310e` |

Both plan files were hashed before any count or score of a clindamycin dataset was read, and both are in the reference directory as they were hashed. The addendum changes no rule of the plan. It fixes what a restart on another machine needed: the commit, the thread count, a run with no network, and the scratch data roots.

**The study of a sample.** A sample's study is the `pmid` of the raw BV-BRC records that carry its reported clindamycin label. The curated tables keep no study field. A sample whose label records carry no PMID has no study and is never picked. Every sample of the six reports has a curated clindamycin row whose reported call equals the report's label: 0 mismatches over 5376 samples (`results/count-checks.txt`).

**The study the plan picks.** In each dataset, the PMID with the most disagreements, on counts alone and before any fit. Ties go to more samples called resistant by the rule and reported susceptible, then to more samples, then to the smaller PMID.

**The caveat rule.** A released dataset needs a caveat when either holds:

- **Existence.** Removing the study leaves a class under 25 samples, so GRM writes no report for the dataset without it.
- **Dependence.** For the known-gene rule or any of the five learners, on either split, the change in balanced accuracy, recall or specificity when the study leaves is at least 0.05 in size, and its 95% interval excludes 0. Each refit dataset has 36 tested rows: six rows, two splits, three measures. They are 33 distinct tests, because the known-gene rule's three rows are the same on both splits.

The rule returns a caveat or none. It does not say which kind of caveat, so the second table below adds the kind. The 0.05 threshold is an adapted research screen chosen for this check, not a validation criterion of any regulator. So are the plan's words for a share of the disagreements: 50% or more is "most", and 25% to under 50% is "a large share". The words are read beside the study's share of the samples.

## The six datasets

From `results/chosen.tsv`, `results/caveat.tsv` and the review's `review/results/counts/chosen.tsv` and `review/results/counts/by-key.tsv`. Each share of the disagreements stands beside the study's share of the samples twice: over all samples, and over the samples that carry a PMID. Shares over all samples carry an exact (Clopper-Pearson) 95% interval. Genomes of one study are not independent, so the intervals are descriptive.

| Dataset | PMID the plan picks | Its disagreements | Its samples | Samples with no PMID | Its disagreements, samples with a PMID | Its samples, samples with a PMID | Disagreement rate, study / rest |
|---|---|---|---|---|---|---|---|
| `campylobacter_coli__clindamycin` | 26519386 | 32 of 89, 36.0% (26.1% to 46.8%) | 81 of 389, 20.8% | 294 | 32 of 45, 71.1% | 81 of 95, 85.3% | 39.5% / 18.5% |
| `campylobacter_jejuni__clindamycin` | 26519386 | 17 of 24, 70.8% (48.9% to 87.4%) | 31 of 508, 6.1% | 475 | 17 of 18, 94.4% | 31 of 33, 93.9% | 54.8% / 1.5% |
| `enterococcus_faecium__clindamycin` | 38219758 | 277 of 277, 100.0% (98.7% to 100.0%) | 901 of 901, 100.0% | 0 | 277 of 277, 100.0% | 901 of 901, 100.0% | 30.7% / no rest |
| `staphylococcus_aureus__clindamycin` | 34100643 | 182 of 219, 83.1% (77.5% to 87.8%) | 659 of 1543, 42.7% | 72 | 182 of 218, 83.5% | 659 of 1471, 44.8% | 27.6% / 4.2% |
| `streptococcus_agalactiae__clindamycin` | 26135871 | 1 of 1, 100.0% (2.5% to 100.0%) | 85 of 141, 60.3% | 1 | 1 of 1, 100.0% | 85 of 140, 60.7% | 1.2% / 0.0% |
| `streptococcus_pneumoniae__clindamycin` | 37823632 | 18 of 22, 81.8% (59.7% to 94.8%) | 1461 of 1894, 77.1% | 301 | 18 of 18, 100.0% | 1461 of 1593, 91.7% | 1.2% / 0.9% |

"Left without it" is the count of reported resistant and reported susceptible samples outside the study.

| Dataset | Left without it, R / S | Scored without it | Plan's verdict | Kind |
|---|---|---|---|---|
| `campylobacter_coli__clindamycin` | 60 / 248 | yes | Caveat, dependence: 9 of 36 rows pass (9 of 33 distinct) | Measured dependence |
| `campylobacter_jejuni__clindamycin` | 8 / 469 | no | Caveat, existence: 8 reported resistant left, under 25 | Class balance, dependence not measured |
| `enterococcus_faecium__clindamycin` | 0 / 0 | no | Caveat, existence: the study is the whole dataset | One study |
| `staphylococcus_aureus__clindamycin` | 325 / 559 | yes | Caveat, dependence: 22 of 36 rows pass (19 of 33 distinct) | Measured dependence |
| `streptococcus_agalactiae__clindamycin` | 38 / 18 | no | Caveat, existence: 18 reported susceptible left, under 25 | Class balance, dependence not measured |
| `streptococcus_pneumoniae__clindamycin` | 41 / 392 | yes | No caveat: 0 of 36 rows pass (0 of 33 distinct) | No dependence under the rule |

In `staphylococcus_aureus__clindamycin`, PMID 34100643 holds 182 of the 219 disagreements (83.1%) on 42.7% of the samples, and 83.5% of the disagreements on 44.8% of the samples that carry a PMID. In the two Campylobacter datasets the study's share of the disagreements is well above its share of the samples only because most samples carry no PMID (294 of 389 and 475 of 508), and those disagree less often (44 of 294, 15.0%, and 6 of 475, 1.3%). Among samples with a PMID the two shares are close or reversed: PMID 26519386 holds 32 of 45 disagreements on 81 of 95 samples in `campylobacter_coli__clindamycin`, and 17 of 18 on 31 of 33 in `campylobacter_jejuni__clindamycin`. In `campylobacter_coli__clindamycin` the only other PMID, 31665255, disagrees on 13 of its 14 samples (92.9%), against 39.5% for PMID 26519386. So the Campylobacter contrast is between one study and samples with no PMID, not between studies. In `streptococcus_pneumoniae__clindamycin` and `streptococcus_agalactiae__clindamycin` the share mostly reflects the study's size, and the second has a single disagreement. `results/by-study.tsv` and `results/by-pmid.tsv` give every study of every dataset, with its counts.

## The refits

The three datasets that keep 25 samples of each class without their study were refit twice at revision 14's settings, with GRM's own commands: `grm count`, `grm structure` and `grm matrix`, then `grm learn --only <dataset>`. Arm "with" is the released cohort. Arm "without" uses a data root whose `curated/phenotypes.jsonl` lacks only the study's clindamycin rows (`results/logs/root-checks.txt`). The study's genomes stay in the species matrix and the lineage clustering, so every remaining sample keeps its matrix columns and its lineage id.

For all three datasets, arm "with" reproduces the release: the same sample ids, labels, known-gene calls, lineage ids, input fingerprint and `eval_config`, and every out-of-fold score of all six rows on both splits (`results/scores/<dataset>/acceptance.txt`). What changes by construction is the sample set, the fold plans drawn over it, the per-fold feature screen and setting selection, and the minimum support, which GRM derives from the sample count.

In the tables below a value reads "released cohort to cohort without the study", and a change is the second minus the first. Its 95% interval is a percentile bootstrap of 2000 replicates that resamples whole lineage clusters. Values are given to four decimals, as `results/refit-table.tsv` holds them.

### staphylococcus_aureus__clindamycin without PMID 34100643

1543 samples to 884, 43 lineage clusters.

| Split | Row | Recall | Log loss | Balanced accuracy | Change in recall (95% interval) | Change in balanced accuracy (95% interval) |
|---|---|---|---|---|---|---|
| random | scm | 0.7965 to 0.9692 | 1.3849 to 0.3752 * | 0.8447 to 0.9614 | +0.1727 (+0.0798 to +0.3933) | +0.1166 (+0.0725 to +0.2181) |
| random | cart | 0.8470 to 0.9631 | 0.6313 to 0.2419 | 0.8412 to 0.9529 | +0.1160 (+0.0419 to +0.2741) | +0.1117 (+0.0713 to +0.1797) |
| random | rf | 0.8283 to 0.9046 | 0.3881 to 0.2327 | 0.8500 to 0.9237 | +0.0763 (-0.1879 to +0.1831) | +0.0737 (-0.0437 to +0.1154) |
| random | gbt | 0.8326 to 0.9108 | 0.3549 to 0.2006 | 0.8610 to 0.9268 | +0.0782 (-0.1785 to +0.1966) | +0.0658 (-0.0524 to +0.1093) |
| random | lasso | 0.8268 to 0.9138 | 0.3689 to 0.2122 | 0.8569 to 0.9292 | +0.0870 (-0.1746 to +0.2189) | +0.0722 (-0.0485 to +0.1188) |
| both | known-gene rule | 0.9163 to 0.9877 | 1.3073 to 0.3856 * | 0.8634 to 0.9643 | +0.0714 (+0.0266 to +0.1782) | +0.1009 (+0.0598 to +0.1565) |
| lineage | scm | 0.4040 to 0.8062 | 3.3070 to 1.2608 * | 0.6191 to 0.8512 | +0.4021 (+0.1672 to +0.6339) | +0.2321 (+0.1071 to +0.3691) |
| lineage | cart | 0.2525 to 0.8154 | 1.5569 to 0.5074 | 0.5451 to 0.8442 | +0.5629 (+0.0427 to +0.8245) | +0.2991 (+0.0457 to +0.4587) |
| lineage | rf | 0.2121 to 0.3631 | 0.6499 to 0.4552 | 0.5637 to 0.6368 | +0.1510 (-0.1397 to +0.6577) | +0.0731 (-0.0685 to +0.3245) |
| lineage | gbt | 0.0620 to 0.7477 | 0.7461 to 0.5091 | 0.4640 to 0.7978 | +0.6856 (+0.0401 to +0.8885) | +0.3339 (+0.0326 to +0.4362) |
| lineage | lasso | 0.1472 to 0.8185 | 0.7306 to 0.5891 | 0.5159 to 0.8153 | +0.6713 (+0.2180 to +0.8125) | +0.2994 (+0.0829 to +0.4001) |

\* The set covering machine and the known-gene rule give scores of 0 or 1, so their log loss is a function of their error count, with scores clipped at 1e-4.

Balanced accuracy is higher without the study on every row, by +0.0658 to +0.3339, and 22 of the 36 tested rows pass the plan's rule (19 of the 33 distinct rows), on both splits and for the known-gene rule. Part of each change comes from no longer scoring the study's samples (composition) and part from no longer fitting on them (training). On the random split, the known-gene rule's whole change is composition, because it is fitted on nothing. For the learners there, composition is +0.0577 to +0.0669 of each change in balanced accuracy, and training is +0.0536 for scm, +0.0494 for cart and +0.0068 to +0.0081 for rf, gbt and lasso. On the lineage split most of the change for gbt (+0.3307 of +0.3339) and lasso (+0.2559 of +0.2994) is training. The training part is more than the study's labels: it also holds a new fold plan, a smaller training set (884 samples against 1543) and a lower minimum support (7 to 4). No placebo arm could be drawn: outside the study the cohort holds 325 reported resistant samples, and the study holds 368.

### staphylococcus_aureus__clindamycin over 11 fold plans (exploratory)

The intervals above hold each arm's fold plan fixed. The review measured how far the fold plan alone moves these figures. This part is exploratory: the check's plan did not foresee it. The review's plan for this part is `review/ADDENDUM.md`, sha256 `dccb2cccb6b8313c2432d0b4b16330e88ff6a60973e44b2c39d62da06fc37b3e`, recorded 2026-11-15T16:38:12Z, before any of its runs. It used `grm` built at `c1f36fc`, whose default seed reproduces the check's reports of both arms (18516 and 10608 out-of-fold scores, 0 different), and fitted each arm at that seed and at seeds 1 to 10: 11 distinct fold plans per arm and split (`review/results/foldplans/check.txt`). A change stands outside the spread when the two arms' ranges over the 11 fold plans do not overlap and lie the way the change points.

Of the 22 passing rows (19 distinct), 20 stand outside the spread (17 distinct). The addendum fixed the test for balanced accuracy and recall, and 14 of those 16 passing rows (12 of 14 distinct) meet it; the 6 passing specificity rows (5 distinct), shown beside them, all meet it. The two that do not are the set covering machine's rows on the lineage split, so its change there (+0.2321 in balanced accuracy, +0.4021 in recall) lies inside the spread. Over the 11 fold plans its balanced accuracy runs from 0.6191 to 0.6591 with the study and from 0.6104 to 0.8782 without it, and its change at the same seed runs from -0.0165 to +0.2475 in balanced accuracy and from -0.0348 to +0.4727 in recall (`review/results/foldplans/tests.tsv` and `review/results/foldplans/spread.tsv`).

### campylobacter_coli__clindamycin without PMID 26519386, in short

389 samples to 308, 226 lineage clusters. The study holds 32 of the 83 disagreements where the rule calls susceptible and the reported phenotype is resistant, and none of the 6 the other way. Without it, recall and balanced accuracy fall for every learner on both splits, and log loss rises for every learner on both splits. The known-gene rule's balanced accuracy rises by +0.0146, which does not pass.

Ten placebo arms each removed 34 reported resistant and 47 reported susceptible samples at random from outside the study, and were refit the same way. 9 of the 36 rows pass (9 of the 33 distinct rows), all of them decreases, and each lies below the lowest of its 10 placebo changes. Ten placebo arms bound how strong that can be: a change below all 10 ranks 1 of 11, the lowest rank they allow (`review/results/refit/placebo-check.txt`).

| Split | Row | Measure | Released to without | Change (95% interval) | 10 placebo changes, lowest to highest | Log loss |
|---|---|---|---|---|---|---|
| random | rf | recall | 0.5426 to 0.2833 | -0.2592 (-0.4096 to -0.1183) | -0.2259 to +0.0241 | 0.3667 to 0.4016 |
| random | rf | balanced accuracy | 0.7492 to 0.6316 | -0.1177 (-0.1915 to -0.0474) | -0.0990 to +0.0260 | 0.3667 to 0.4016 |
| lineage | scm | recall | 0.8511 to 0.7500 | -0.1011 (-0.2039 to -0.0131) | -0.0511 to +0.0489 | 0.4973 to 0.5384 * |
| lineage | cart | recall | 0.7979 to 0.7167 | -0.0812 (-0.1582 to -0.0114) | +0.0021 to +0.1021 | 0.2554 to 0.2891 |
| lineage | rf | recall | 0.5213 to 0.1500 | -0.3713 (-0.4879 to -0.2454) | -0.1379 to -0.0046 | 0.3715 to 0.4101 |
| lineage | rf | balanced accuracy | 0.7420 to 0.5589 | -0.1831 (-0.2439 to -0.1176) | -0.0604 to +0.0083 | 0.3715 to 0.4101 |
| lineage | gbt | recall | 0.8404 to 0.7333 | -0.1071 (-0.1969 to -0.0307) | -0.0238 to +0.0762 | 0.2549 to 0.3383 |
| lineage | gbt | balanced accuracy | 0.9016 to 0.8445 | -0.0571 (-0.1017 to -0.0163) | -0.0194 to +0.0507 | 0.2549 to 0.3383 |
| lineage | lasso | recall | 0.8191 to 0.7000 | -0.1191 (-0.2120 to -0.0400) | -0.0025 to +0.0975 | 0.2755 to 0.3450 |

Every row of the dataset, with recall beside log loss and balanced accuracy, is in `results/refit-table.tsv`, and every tested row with its placebo range is in `results/scores/campylobacter_coli__clindamycin/caveat.tsv`.

### streptococcus_pneumoniae__clindamycin without PMID 37823632

1894 samples to 433, 103 lineage clusters: no row passes, every change in balanced accuracy lies between -0.0660 and +0.0120 with an interval that includes 0, released balanced accuracy is 0.9405 to 0.9637 on every row, and recall is 0.9036 to 0.9337 with the study and 0.8293 to 0.9512 without it.

## What a released clindamycin figure can be read as

- **`staphylococcus_aureus__clindamycin`: measured dependence.** The figures measure agreement with reported labels of which 659 of 1543 come from PMID 34100643, a study that holds 182 of the dataset's 219 disagreements. Without that study, 19 of 33 distinct rows (22 of 36 counting the known-gene rule on both splits) rise by 0.05 or more with lineage-cluster intervals that exclude 0, and all but the set covering machine on the lineage split stand outside the spread of 11 fold plans. Those 19 rows are associated with that one study's labels or genomes. The other 14 distinct rows do not pass the rule, and neither result shows which reading is wrong.
- **`campylobacter_coli__clindamycin`: measured dependence.** The released learner figures include the 81 samples of PMID 26519386, and they are higher with them than without them. Without it, 9 of 36 rows (9 of 33 distinct) fall by 0.05 or more: recall for all five learners on the lineage split and for rf on the random split, and balanced accuracy for rf on both splits and for gbt on the lineage split. The known-gene rule's figure does not move by the plan's rule. Each of the 9 changes that pass is lower than the lowest change of 10 random removals of the same size. The study's share of the disagreements stands out only against samples with no PMID.
- **`enterococcus_faecium__clindamycin`: one study.** Every sample comes from PMID 38219758 (901 of 901), so the figures describe that study's genomes and labels alone. No figure without that study can exist. The rule and the reported phenotype disagree on 277 of the 901 samples.
- **`campylobacter_jejuni__clindamycin`: class balance, dependence not measured.** 17 of the dataset's 25 reported resistant samples come from PMID 26519386. Without it 8 remain, under the class minimum of 25 this release uses, so whether the figures depend on that study could not be measured, and the resistant class of every released figure rests mostly on that one study.
- **`streptococcus_agalactiae__clindamycin`: class balance, dependence not measured.** 8 of the dataset's 26 reported susceptible samples come from PMID 26135871. Without it 18 remain, under 25, so the dataset cannot be scored. The dataset has a single disagreement, and the study is picked on it.
- **`streptococcus_pneumoniae__clindamycin`: no dependence under the rule.** PMID 37823632 holds 18 of 22 disagreements, close to its share of the samples (1461 of 1894). Without it 0 of 36 rows pass (0 of 33 distinct). That does not show that the figures do not depend on that study: some intervals are wide (the change in rf's recall on the lineage split is -0.0924, with an interval of -0.4318 to +0.0516), and the fold plans of this dataset were not varied.

## Limits

- **Association only.** The change when a study leaves is an association with that study's labels or genomes. It does not show which label is wrong. It does not separate label quality from the study's lineages, testing method or sampling. It is no statement about any gene.
- **Composition and training.** Part of each change comes from no longer scoring the study's samples and part from no longer fitting on them. The training part also includes a new fold plan, a smaller training set and, for `staphylococcus_aureus__clindamycin`, a lower minimum support (7 to 4). The bootstrap intervals do not cover the fold plan. Over 11 fold plans (exploratory), every passing row of `staphylococcus_aureus__clindamycin` but the set covering machine on the lineage split stands outside that spread. The fold plans of `streptococcus_pneumoniae__clindamycin` were not varied. In `campylobacter_coli__clindamycin` each of the 10 placebo arms draws its own fold plan, so the placebo range also reflects the fold plan.
- **Placebo arms for one dataset only.** For `staphylococcus_aureus__clindamycin` and `streptococcus_pneumoniae__clindamycin`, the rest of the cohort holds fewer reported resistant samples than the study (325 against 368, and 41 against 125), so the plan's matched placebo could not be drawn, and their training effects have no placebo reading.
- **Wide lineage intervals.** `staphylococcus_aureus__clindamycin` has 43 lineage clusters, so its lineage-split intervals are wide, for example +0.0326 to +0.4362 for the change in gbt's balanced accuracy. No refit dataset falls under the plan's flag of 20 clusters.
- **Many rows.** 36 rows are tested per refit dataset (33 distinct), and 22 and 9 pass (19 and 9 distinct). In `streptococcus_pneumoniae__clindamycin` the scm and cart rows also give the same values on both splits, though each split is its own fit, and none of them passes. Rows within a dataset share samples and are not independent.
- **The study field.** Samples whose label records carry no PMID are not a study: 294 in `campylobacter_coli__clindamycin`, 475 in `campylobacter_jejuni__clindamycin`, 72 in `staphylococcus_aureus__clindamycin`, 1 in `streptococcus_agalactiae__clindamycin` and 301 in `streptococcus_pneumoniae__clindamycin` (`results/by-study.tsv`). In the two Campylobacter datasets they are most of the rest that the study is compared with.
- **Existence rests on class counts.** For `campylobacter_jejuni__clindamycin` and `streptococcus_agalactiae__clindamycin` the caveat says the dataset depends on one study for its class balance. Whether their figures depend on that study was not measured. The study picked in `streptococcus_agalactiae__clindamycin` has 1 disagreement.
- **Not refit, not rechecked.** `enterococcus_faecium__clindamycin`, `campylobacter_jejuni__clindamycin` and `streptococcus_agalactiae__clindamycin` were not refit in either arm, so their released figures were not rechecked here.
- **The labels are the database's.** A reported phenotype is the label BV-BRC holds, not the result of a reference method run for this check.
- **The October record is not reproduced.** Its scripts and files are lost, and it states no measure. It comes from revision 7, with a cohort of the same size whose members are not known.

**History (exploratory, checked after the counts were read).** The October record of a revision 7 cohort of 1543 samples (quoted in `PLAN.md`, section 0, and [docs/work-program-20261027.md](work-program-20261027.md)) gives 161 known-gene resistant calls among 850 reported susceptible samples, 128 of them from one study, and a set covering machine comparison of 0.772809 against 0.870555. Revision 14's `staphylococcus_aureus__clindamycin` has the same four counts, the 128 sit in PMID 34100643, and the set covering machine's lineage-split balanced accuracy on the 884 samples outside that study, scored from the released fit before any refit, is 0.772809 (`review/results/exploratory/october-match.txt`). The October run very likely studied this dataset and this study, and compared the samples outside the study fitted with and without it. The record itself names neither the species nor the PMID. Its second figure, 0.870555, matches none of the figures searched, and lies inside the range of the arm without the study over 11 fold plans (0.6104 to 0.8782).

## Files and commands

[`reference/clindamycin-labels/rev14-20261115/README.md`](../reference/clindamycin-labels/rev14-20261115/README.md) lists each file, the command that wrote it and its sha256, and [`review/README.md`](../reference/clindamycin-labels/rev14-20261115/review/README.md) does the same for the review's files. The check ran on 15 November 2026 on a machine that holds the public data root and the release, with the scripts there. The per-sample tables and the lists of removed genome ids are not in the repository.

`scripts/clindamycin-labels/site_data.py` writes `web/src/data/clindamycin-labels.json` from the reference directory's files and the release's report checksums in `reference/lineage-control/rev14-274-20261112.reports.sha256`, the web app shows each of the six datasets' kind with the counts behind it on that dataset's report from that file, and `tests/clindamycin_labels` fails if the file differs from what the script writes or if one of its counts differs from this page.
