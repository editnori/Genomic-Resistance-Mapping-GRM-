# Small-cohort result

The release fits a dataset only when each class holds at least 25 genomes. A dataset that passes gets an estimate from its own lineage split. This page asks how far such an estimate lies from the truth it stands for when the cohort is smaller than 25 a class. It reads block C1, part B, of the assumption plan (`docs/assumption-plan.md`).

The reading's rules were fixed at 20:03Z on 22 Nov 2026, with no score read, and were public before the reading program ran (`docs/assumption-plan-corrections.md`, entries "C1 design", "C1 scorer" and "C1 reading"). The program ran once. The code is in `docs/cohort-draw.md`.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margins of 0.02 and 0.05 are an adapted research screen. They are no applicable FDA validation criterion.

## Reading

**Result.** No lower line is a candidate on the 89 datasets. A line needs 3 of the 5 learners to meet the plan's release row. Line 10 has 0, line 15 has 0 and line 20 has 1 (`rf`). The release keeps its 25-per-class line, and nothing changes in the release or on the site.

**Lines 10 and 15.** At 10 a class the estimate lies further from the truth than at 25 for `scm`, `cart` and `rf`: a median of 0.031 to 0.039 further in balanced accuracy, each interval wholly above 0. At 15 the same three are 0.015 to 0.022 further.

**Line 20.** Each learner's median difference is under 0.01, and 4 of the 5 intervals end below 0.02. Three of those four fail the row on coverage, which is more than 0.05 above 25's. The section Seen after says what that leaves open.

**Twin.** The r = 25 cohorts fitted again at fold seed 45 meet the row on 3 of 5 learners (`cart`, `rf` and `gbt`). So the row holds against a change of fold seed alone.

**Line 50.** At 50 a class the estimate lies closer to the truth than at 25 for every learner, by a median of 0.017 to 0.029. Line 50 decides nothing.

The reading is the first run's: the 89 datasets of BIG outside M. The plan's set for part B is BIG, 108 datasets. The count over the 108 is taken when the maintainer's label rule lets the other 19 run, with these rules, this scorer and these arguments unchanged.

## Data

**Datasets.** 89 datasets in 10 species. They are the release's datasets with a lineage split whose smaller class holds 200 or more genomes (set BIG), less the 19 that wait for the maintainer's label rule (set M). The section Datasets names all 89.

**Cohorts.** From each dataset, cohorts of r resistant and r susceptible genomes are drawn at r = 10, 15, 20, 25 and 50, five draws each. A draw takes genomes round robin over the dataset's lineages. Its cohort at a smaller r lies inside its cohort at a larger one.

**Fits.** Each cohort is fitted with the release's rule at fold seed 42. The r = 25 cohorts are fitted again at fold seed 45: the twin, the same cohorts with only the fold seed changed. That makes 30 results a dataset and 2,670 in all.

**Estimate.** The cohort's own lineage split gives the estimate: 5 outer folds, 3 inner folds and the five learners (`scm`, `cart`, `rf`, `gbt` and `lasso`). The cohort also runs the random split, and the study split on the 27 datasets with a reviewed study plan.

**Truth.** The cohort's final models score the dataset's other genomes: all of them, and those in lineages the cohort lacks. Genomes k-mer-identical to a cohort genome are left out. The scored genomes come from the release's dataset, inside its 2,000-genome cap.

**Runs.** Each species ran on a machine of its own. Each was then scored again, from the same fitted results, by the scorer with a review's fixes (commit `461a552`). Only that second scoring is read.

**Part A.** Part A counts what each line would admit, with no fit. It was taken on 21 Nov 2026, before any cohort fit, from the rebuilt revision 14 table:

```
grm --data <data root> curate --cap 2000 --min-class <line>
```

A line admits a species and antibiotic pair with at least that many genomes in each class. At 25 the count is the 277 curated datasets, of which 274 have a report in the release. The six rows are separate lines over one table and do not add.

| Per class | Datasets | Added | Removed | Species | New species | Their datasets |
|---|---|---|---|---|---|---|
| 5 | 457 | 180 | 0 | 35 | 16 | 68 |
| 10 | 373 | 96 | 0 | 29 | 10 | 26 |
| 15 | 332 | 55 | 0 | 25 | 8 | 14 |
| 20 | 301 | 24 | 0 | 18 | 4 | 6 |
| 25 | 277 | 0 | 0 | 0 | 0 | 0 |
| 50 | 212 | 0 | 65 | 0 | 0 | 0 |

Added and Removed are against the 277. Species counts the species with an added dataset. New species counts those among them with no dataset at 25, and the last column their added datasets. So line 20 would add 24 datasets in 18 species, 6 of them in 4 species the release does not cover. Line 10 would add 96 in 29 species.

## Count

The measure is the estimate's error: its absolute difference from the truth in balanced accuracy, per learner. Each lower line is read against 25 on its own, and the page names no line as the one to take.

**Pairing.** A dataset's difference is the mean, over the draws available at both lines, of the error at the line less the error at 25. A dataset enters with 3 or more such draws. All 89 entered at every line, with 5 draws each, so each row counts 445 results.

**Three counts.** Among the datasets that enter, the error at the line is smaller by more than 0.02 (better), larger by more than 0.02 (worse), or within 0.02 (same).

**Median.** The median of the differences has a 95% interval from 4,000 species resamples, read at the 101st and 3,900th values. It resamples the species with an entered dataset, 10 here.

**Coverage.** How often the levelled interval on the estimate covers the truth, at the line and at 25, over the same results. The levelled interval is built from the split's pooled recall and specificity and their Wilson 95% intervals, so it is built to the same nominal level at every line.

**Release row.** A learner meets the row when its interval lies wholly below 0.02 and its coverage is within 0.05 of 25's. A line is a candidate on the 89 when 3 or more of the 5 learners meet. For a learner that does not meet, the last column gives the part that failed.

One row per lower line decides: the lineage split, the truth on all scored genomes, the levelled interval. The gap is the coverage at the line less the coverage at 25.

| Line | Learner | Better | Worse | Same | Median | Interval | At line | At 25 | Gap | Row |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | `scm` | 14 | 50 | 25 | 0.031 | 0.017 to 0.039 | 331 of 445 | 258 of 445 | 0.164 | the interval's upper end and the coverage |
| 10 | `cart` | 14 | 58 | 17 | 0.038 | 0.024 to 0.056 | 328 of 445 | 270 of 445 | 0.130 | the interval's upper end and the coverage |
| 10 | `rf` | 13 | 53 | 23 | 0.039 | 0.024 to 0.046 | 364 of 445 | 313 of 445 | 0.115 | the interval's upper end and the coverage |
| 10 | `gbt` | 33 | 35 | 21 | 0.004 | -0.015 to 0.039 | 348 of 445 | 282 of 445 | 0.148 | the interval's upper end and the coverage |
| 10 | `lasso` | 34 | 29 | 26 | -0.002 | -0.024 to 0.015 | 357 of 445 | 257 of 445 | 0.225 | the coverage |
| 15 | `scm` | 20 | 42 | 27 | 0.018 | 0.008 to 0.025 | 303 of 445 | 258 of 445 | 0.101 | the interval's upper end and the coverage |
| 15 | `cart` | 18 | 41 | 30 | 0.015 | 0.004 to 0.029 | 301 of 445 | 270 of 445 | 0.070 | the interval's upper end and the coverage |
| 15 | `rf` | 14 | 47 | 28 | 0.022 | 0.013 to 0.038 | 333 of 445 | 313 of 445 | 0.045 | the interval's upper end |
| 15 | `gbt` | 25 | 34 | 30 | 0.005 | -0.004 to 0.022 | 315 of 445 | 282 of 445 | 0.074 | the interval's upper end and the coverage |
| 15 | `lasso` | 27 | 35 | 27 | 0.009 | -0.005 to 0.021 | 311 of 445 | 257 of 445 | 0.121 | the interval's upper end and the coverage |
| 20 | `scm` | 19 | 30 | 40 | 0.006 | 0.001 to 0.017 | 283 of 445 | 258 of 445 | 0.056 | the coverage |
| 20 | `cart` | 19 | 37 | 33 | 0.009 | -0.006 to 0.023 | 279 of 445 | 270 of 445 | 0.020 | the interval's upper end |
| 20 | `rf` | 13 | 33 | 43 | 0.009 | 0.002 to 0.018 | 307 of 445 | 313 of 445 | -0.013 | meets |
| 20 | `gbt` | 26 | 30 | 33 | 0.005 | -0.004 to 0.014 | 307 of 445 | 282 of 445 | 0.056 | the coverage |
| 20 | `lasso` | 24 | 32 | 33 | 0.005 | -0.004 to 0.013 | 280 of 445 | 257 of 445 | 0.052 | the coverage |

A candidate would be a line the maintainer may then decide to fetch and fit what part A counts for. It would change nothing in the release by itself: a change of the release's line needs those fits and the maintainer's ruling.

**Twin.** The twin is the r = 25 cohorts fitted again at fold seed 45, read against r = 25 at seed 42 as if it were a line. Only the fold seed changed, so its count is the count with no change: how far a seed alone moves the error. It is never a candidate.

| Learner | Seed better | Seed worse | Seed same | Seed median | Seed interval | At seed 45 | At seed 42 | Seed gap | Twin row |
|---|---|---|---|---|---|---|---|---|---|
| `scm` | 28 | 21 | 40 | -0.005 | -0.012 to 0.000 | 286 of 445 | 258 of 445 | 0.063 | the coverage |
| `cart` | 26 | 17 | 46 | -0.004 | -0.008 to 0.000 | 289 of 445 | 270 of 445 | 0.043 | meets |
| `rf` | 9 | 19 | 61 | 0.004 | -0.001 to 0.006 | 319 of 445 | 313 of 445 | 0.013 | meets |
| `gbt` | 28 | 18 | 43 | -0.004 | -0.012 to 0.010 | 295 of 445 | 282 of 445 | 0.029 | meets |
| `lasso` | 29 | 20 | 40 | -0.006 | -0.015 to 0.000 | 287 of 445 | 257 of 445 | 0.067 | the coverage |

The fixed reading: when 3 or more learners meet on the twin's deciding row, the row holds against a change of fold seed alone. Here 3 meet. Had fewer met, a line that is no candidate would be written "not shown to be within the margin". Rules' condition 3 is written for a gain and is not applied here.

The scorer writes a `candidate` cell on the twin's row, and `decides` and `candidate` cells on line 50's row, as on any line. `count/count.tsv` keeps them. The page's tables carry the twin's cell as `three_or_more_meet`.

**Line 50.** Line 50 is the legacy line and lies above 25. Its row is written, decides nothing and is never a candidate. It is read only as the other end of the trend over r.

| Learner | Better at 50 | Worse at 50 | Same at 50 | Median at 50 | Interval at 50 | Covered at 50 | Covered at 25 |
|---|---|---|---|---|---|---|---|
| `scm` | 49 | 9 | 31 | -0.023 | -0.034 to -0.015 | 248 of 445 | 258 of 445 |
| `cart` | 45 | 8 | 36 | -0.021 | -0.026 to -0.014 | 242 of 445 | 270 of 445 |
| `rf` | 43 | 15 | 31 | -0.017 | -0.026 to -0.006 | 251 of 445 | 313 of 445 |
| `gbt` | 50 | 10 | 29 | -0.025 | -0.036 to -0.018 | 253 of 445 | 282 of 445 |
| `lasso` | 49 | 10 | 30 | -0.029 | -0.038 to -0.014 | 242 of 445 | 257 of 445 |

**Other rows.** These rows are written beside the count and decide nothing. Each cell gives the learners that meet, of 5, and in brackets the datasets that entered, of the datasets that row has. The first row is the one that decides, repeated for comparison.

| Split | Truth | Interval | Line 10 | Line 15 | Line 20 | Twin |
|---|---|---|---|---|---|---|
| the lineage split | all scored genomes | levelled | 0 of 5 (89 of 89) | 0 of 5 (89 of 89) | 1 of 5 (89 of 89) | 3 of 5 (89 of 89) |
| the lineage split | all scored genomes | fold range | 0 of 5 (89 of 89) | 0 of 5 (89 of 89) | 4 of 5 (89 of 89) | 4 of 5 (89 of 89) |
| the lineage split | new lineages | levelled | 0 of 5 (78 of 89) | 0 of 5 (78 of 89) | 4 of 5 (78 of 89) | 3 of 5 (78 of 89) |
| the lineage split | new lineages | fold range | 0 of 5 (78 of 89) | 0 of 5 (78 of 89) | 4 of 5 (78 of 89) | 5 of 5 (78 of 89) |
| the random split | all scored genomes | levelled | 0 of 5 (89 of 89) | 0 of 5 (89 of 89) | 3 of 5 (89 of 89) | 5 of 5 (89 of 89) |
| the random split | all scored genomes | fold range | 0 of 5 (89 of 89) | 0 of 5 (89 of 89) | 4 of 5 (89 of 89) | 5 of 5 (89 of 89) |
| the random split | new lineages | levelled | 0 of 5 (78 of 89) | 0 of 5 (78 of 89) | 2 of 5 (78 of 89) | 5 of 5 (78 of 89) |
| the random split | new lineages | fold range | 0 of 5 (78 of 89) | 0 of 5 (78 of 89) | 4 of 5 (78 of 89) | 5 of 5 (78 of 89) |
| the study split | all scored genomes | levelled | 0 of 5 (0 of 27) | 0 of 5 (6 of 27) | 0 of 5 (10 of 27) | 1 of 5 (16 of 27) |
| the study split | all scored genomes | fold range | 0 of 5 (0 of 27) | 0 of 5 (6 of 27) | 0 of 5 (10 of 27) | 3 of 5 (16 of 27) |
| the study split | new lineages | levelled | 0 of 5 (0 of 27) | 0 of 5 (6 of 27) | 0 of 5 (10 of 27) | 2 of 5 (14 of 27) |
| the study split | new lineages | fold range | 0 of 5 (0 of 27) | 0 of 5 (6 of 27) | 0 of 5 (10 of 27) | 1 of 5 (14 of 27) |

**Fold range.** The fold range is the lowest to the highest of the five fold balanced accuracies, which is what a report shows. It has no nominal level and widens as r falls, so its coverage rises with its width. A result with a fold that has no balanced accuracy leaves its count.

**Study split.** A small cohort seldom has a study split that runs: no dataset of the 27 entered at line 10, 6 at line 15 and 10 at line 20.

**Each learner.** Every learner's own row is in the count file as well, and decides nothing alone.

**Width.** The levelled interval is wider at a smaller r. The table gives its median width over the 445 results, in balanced accuracy.

| Learner | Width at 10 | Width at 15 | Width at 20 | Width at 25 | Width at 50 |
|---|---|---|---|---|---|
| `scm` | 0.353 | 0.298 | 0.257 | 0.225 | 0.147 |
| `cart` | 0.354 | 0.298 | 0.256 | 0.226 | 0.151 |
| `rf` | 0.359 | 0.303 | 0.264 | 0.233 | 0.154 |
| `gbt` | 0.341 | 0.291 | 0.257 | 0.226 | 0.153 |
| `lasso` | 0.339 | 0.289 | 0.250 | 0.221 | 0.148 |

## Beside it

Recall, specificity and log loss are given for the estimate and for the truth, on the lineage split with the truth on all scored genomes. They run over the same 445 results as the count, with 25 beside on the same draws. They stand beside the error and decide nothing.

The pooled figure comes first, as a count over its genomes. The median beside it is over the 89 datasets, where a dataset's value is the mean over its results of each result's rate. The last column gives in how many results the truth's rate lies inside the estimate's Wilson 95% interval.

| Line | Learner | Estimate recall | Estimate median | Truth recall | Truth median | Inside |
|---|---|---|---|---|---|---|
| 10 | `scm` | 0.590 (2,625 of 4,450) | 0.580 | 0.638 (161,020 of 252,418) | 0.629 | 297 of 445 |
| 10 | `cart` | 0.608 (2,706 of 4,450) | 0.600 | 0.670 (168,998 of 252,418) | 0.651 | 278 of 445 |
| 10 | `rf` | 0.569 (2,532 of 4,450) | 0.560 | 0.711 (179,519 of 252,418) | 0.692 | 263 of 445 |
| 10 | `gbt` | 0.778 (3,464 of 4,450) | 0.780 | 0.790 (199,425 of 252,418) | 0.782 | 225 of 445 |
| 10 | `lasso` | 0.810 (3,606 of 4,450) | 0.820 | 0.844 (213,108 of 252,418) | 0.842 | 209 of 445 |
| 15 | `scm` | 0.630 (4,207 of 6,675) | 0.613 | 0.702 (175,550 of 250,168) | 0.717 | 257 of 445 |
| 15 | `cart` | 0.642 (4,287 of 6,675) | 0.627 | 0.709 (177,269 of 250,168) | 0.696 | 255 of 445 |
| 15 | `rf` | 0.646 (4,310 of 6,675) | 0.640 | 0.749 (187,368 of 250,168) | 0.764 | 245 of 445 |
| 15 | `gbt` | 0.739 (4,931 of 6,675) | 0.733 | 0.791 (197,926 of 250,168) | 0.818 | 238 of 445 |
| 15 | `lasso` | 0.763 (5,091 of 6,675) | 0.760 | 0.833 (208,302 of 250,168) | 0.837 | 206 of 445 |
| 20 | `scm` | 0.665 (5,917 of 8,900) | 0.660 | 0.723 (179,328 of 247,920) | 0.747 | 256 of 445 |
| 20 | `cart` | 0.670 (5,960 of 8,900) | 0.660 | 0.729 (180,641 of 247,920) | 0.745 | 247 of 445 |
| 20 | `rf` | 0.666 (5,926 of 8,900) | 0.670 | 0.786 (194,966 of 247,920) | 0.818 | 221 of 445 |
| 20 | `gbt` | 0.742 (6,605 of 8,900) | 0.750 | 0.823 (204,072 of 247,920) | 0.819 | 223 of 445 |
| 20 | `lasso` | 0.763 (6,793 of 8,900) | 0.750 | 0.820 (203,353 of 247,920) | 0.845 | 223 of 445 |
| 25 | `scm` | 0.686 (7,629 of 11,125) | 0.672 | 0.753 (185,000 of 245,670) | 0.765 | 212 of 445 |
| 25 | `cart` | 0.685 (7,626 of 11,125) | 0.696 | 0.757 (185,972 of 245,670) | 0.744 | 248 of 445 |
| 25 | `rf` | 0.689 (7,663 of 11,125) | 0.704 | 0.811 (199,222 of 245,670) | 0.849 | 208 of 445 |
| 25 | `gbt` | 0.744 (8,275 of 11,125) | 0.744 | 0.825 (202,779 of 245,670) | 0.838 | 235 of 445 |
| 25 | `lasso` | 0.759 (8,442 of 11,125) | 0.752 | 0.845 (207,495 of 245,670) | 0.860 | 217 of 445 |

The very major error rate (resistant called susceptible) is 1 less the recall, over the same resistant genomes. The major error rate (susceptible called resistant) is 1 less the specificity, over the same susceptible genomes.

| Line | Learner | Estimate specificity | Median estimate | Truth specificity | Median truth | Truth inside |
|---|---|---|---|---|---|---|
| 10 | `scm` | 0.687 (3,057 of 4,450) | 0.700 | 0.704 (209,485 of 297,613) | 0.715 | 298 of 445 |
| 10 | `cart` | 0.656 (2,917 of 4,450) | 0.660 | 0.664 (197,479 of 297,613) | 0.669 | 293 of 445 |
| 10 | `rf` | 0.654 (2,911 of 4,450) | 0.660 | 0.588 (175,020 of 297,613) | 0.598 | 349 of 445 |
| 10 | `gbt` | 0.401 (1,784 of 4,450) | 0.380 | 0.461 (137,201 of 297,613) | 0.458 | 177 of 445 |
| 10 | `lasso` | 0.380 (1,690 of 4,450) | 0.340 | 0.413 (122,946 of 297,613) | 0.414 | 158 of 445 |
| 15 | `scm` | 0.727 (4,851 of 6,675) | 0.733 | 0.728 (215,024 of 295,347) | 0.761 | 262 of 445 |
| 15 | `cart` | 0.717 (4,789 of 6,675) | 0.720 | 0.741 (218,752 of 295,347) | 0.778 | 278 of 445 |
| 15 | `rf` | 0.665 (4,438 of 6,675) | 0.653 | 0.628 (185,353 of 295,347) | 0.636 | 333 of 445 |
| 15 | `gbt` | 0.570 (3,808 of 6,675) | 0.547 | 0.584 (172,530 of 295,347) | 0.623 | 201 of 445 |
| 15 | `lasso` | 0.552 (3,684 of 6,675) | 0.520 | 0.539 (159,061 of 295,347) | 0.554 | 178 of 445 |
| 20 | `scm` | 0.752 (6,697 of 8,900) | 0.770 | 0.775 (227,181 of 293,088) | 0.813 | 261 of 445 |
| 20 | `cart` | 0.759 (6,754 of 8,900) | 0.780 | 0.769 (225,325 of 293,088) | 0.815 | 299 of 445 |
| 20 | `rf` | 0.692 (6,156 of 8,900) | 0.680 | 0.651 (190,788 of 293,088) | 0.676 | 317 of 445 |
| 20 | `gbt` | 0.632 (5,624 of 8,900) | 0.660 | 0.621 (181,884 of 293,088) | 0.635 | 237 of 445 |
| 20 | `lasso` | 0.635 (5,655 of 8,900) | 0.660 | 0.643 (188,369 of 293,088) | 0.703 | 208 of 445 |
| 25 | `scm` | 0.772 (8,592 of 11,125) | 0.800 | 0.785 (228,405 of 290,810) | 0.840 | 253 of 445 |
| 25 | `cart` | 0.771 (8,577 of 11,125) | 0.792 | 0.784 (228,073 of 290,810) | 0.824 | 260 of 445 |
| 25 | `rf` | 0.729 (8,106 of 11,125) | 0.720 | 0.664 (193,153 of 290,810) | 0.689 | 303 of 445 |
| 25 | `gbt` | 0.695 (7,731 of 11,125) | 0.728 | 0.656 (190,765 of 290,810) | 0.760 | 237 of 445 |
| 25 | `lasso` | 0.695 (7,732 of 11,125) | 0.720 | 0.684 (199,028 of 290,810) | 0.760 | 214 of 445 |

The log loss is a mean over genomes, given with its genome count.

| Line | Learner | Estimate loss | Estimate genomes | Estimate's median | Truth loss | Truth genomes | Truth's median |
|---|---|---|---|---|---|---|---|
| 10 | `scm` | 3.330 | 8,900 | 3.500 | 3.006 | 550,031 | 2.990 |
| 10 | `cart` | 2.260 | 8,900 | 2.220 | 1.811 | 550,031 | 1.828 |
| 10 | `rf` | 0.653 | 8,900 | 0.669 | 0.644 | 550,031 | 0.647 |
| 10 | `gbt` | 0.656 | 8,900 | 0.680 | 0.627 | 550,031 | 0.662 |
| 10 | `lasso` | 0.807 | 8,900 | 0.782 | 0.807 | 550,031 | 0.711 |
| 15 | `scm` | 2.961 | 13,350 | 3.009 | 2.616 | 545,515 | 2.632 |
| 15 | `cart` | 1.703 | 13,350 | 1.782 | 1.229 | 545,515 | 1.256 |
| 15 | `rf` | 0.627 | 13,350 | 0.644 | 0.608 | 545,515 | 0.611 |
| 15 | `gbt` | 0.610 | 13,350 | 0.637 | 0.571 | 545,515 | 0.612 |
| 15 | `lasso` | 0.773 | 13,350 | 0.777 | 0.652 | 545,515 | 0.659 |
| 20 | `scm` | 2.683 | 17,800 | 2.717 | 2.290 | 541,008 | 2.063 |
| 20 | `cart` | 1.367 | 17,800 | 1.320 | 0.999 | 541,008 | 0.870 |
| 20 | `rf` | 0.600 | 17,800 | 0.609 | 0.571 | 541,008 | 0.562 |
| 20 | `gbt` | 0.582 | 17,800 | 0.592 | 0.536 | 541,008 | 0.560 |
| 20 | `lasso` | 0.670 | 17,800 | 0.675 | 0.589 | 541,008 | 0.599 |
| 25 | `scm` | 2.496 | 22,250 | 2.432 | 2.113 | 536,480 | 1.830 |
| 25 | `cart` | 1.212 | 22,250 | 1.116 | 0.840 | 536,480 | 0.699 |
| 25 | `rf` | 0.572 | 22,250 | 0.583 | 0.541 | 536,480 | 0.521 |
| 25 | `gbt` | 0.553 | 22,250 | 0.573 | 0.512 | 536,480 | 0.524 |
| 25 | `lasso` | 0.616 | 22,250 | 0.624 | 0.528 | 536,480 | 0.532 |

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above.

**Coverage side.** The row asks that coverage be within 0.05 of 25's, and the scorer reads that on both sides. Every learner that fails on coverage at lines 10, 15 and 20 has coverage higher than 25's, and none has lower. The interval is wider at a lower line (Width, above), and a wider interval covers the truth more often.

**Line 20.** At line 20 the interval's upper end is below 0.02 for `scm`, `rf`, `gbt` and `lasso`. Of these, `scm`, `gbt` and `lasso` fail on coverage alone, at 0.056, 0.056 and 0.052 above 25's. If higher coverage passed, line 20 would count 4 of 5, line 15 would count 0 of 5 and line 10 would count 1 of 5.

**Open question.** That one-sided reading was not fixed before the scores were read, so it decides nothing here. Whether the row should pass higher coverage is the maintainer's to rule.

**Seed alone.** On the twin, `scm` and `lasso` fail on coverage at 0.063 and 0.067 above, with the fold seed alone changed. So a coverage gap near 0.05 is inside what a fold seed moves for those two learners.

**Nominal level.** The levelled interval is built to a nominal 95%. At 25 it covers the truth in 257 to 313 of 445 results, 0.578 to 0.703. At 10, its widest, it covers the truth in 328 to 364 of 445 results. No line reaches 95%.

**Recall.** At every line, each learner's pooled recall on the cohort's own split is below its recall on the other genomes, and its pooled log loss is above. At 10 a class `gbt` and `lasso` hold recall near 0.8 with a truth specificity of 0.461 and 0.413.

## Limits

**Labels.** The 89 are outside M, so label rules (a), (b) and (d) give them the same labels. Under rule (c), the most resistant call, a dataset holding a pair with an R and an S record is checked again. The 19 datasets in M are left out:

| Species | Left out (M) |
|---|---|
| *Acinetobacter baumannii* | amikacin, ampicillin/sulbactam, imipenem, meropenem, tetracycline, tobramycin |
| *Escherichia coli* | amoxicillin/clavulanate, ampicillin, ciprofloxacin, gentamicin, sulfamethoxazole/trimethoprim |
| *Klebsiella pneumoniae* | amikacin, gentamicin, meropenem, piperacillin/tazobactam |
| *Mycobacterium tuberculosis* | isoniazid |
| *Neisseria gonorrhoeae* | azithromycin |
| *Pseudomonas aeruginosa* | ciprofloxacin |
| *Staphylococcus aureus* | oxacillin |

**Breakpoints.** The run gives each dataset the release's breakpoints. The rebuilt table moves them on 36 of the 89 datasets, and no fit reads the moved ones. No count is split by them, and the per-dataset file has no column for them.

| Species | Moved breakpoint |
|---|---|
| *Escherichia coli* | aztreonam, cefepime, cefotaxime, cefoxitin, ceftazidime, ceftriaxone, cefuroxime, ertapenem, piperacillin/tazobactam, tetracycline, tobramycin, trimethoprim |
| *Klebsiella pneumoniae* | ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefoxitin, ceftazidime, ceftriaxone, ciprofloxacin, ertapenem, imipenem, levofloxacin, sulfamethoxazole/trimethoprim, tetracycline, tobramycin |
| *Mycobacterium tuberculosis* | amikacin, delamanid, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, para-aminosalicylic acid, rifabutin, rifampin |

**Truth.** The truth is the release dataset's other genomes inside its 2,000-genome cap. It is less than every labelled genome of the species and antibiotic.

**Error size.** The tables give each line's error against 25's. No table gives a line's error by itself, so the page does not say how far an estimate at 25 lies from its truth.

**Seeds.** Each line has one fold seed (42) and five draws, in place of Rules' fold seeds 42, 43 and 44. The twin has seed 45 alone.

**Species.** The species interval resamples 10 species, and the datasets are unevenly spread over them, from 3 to 15 a species.

**One database.** Every genome comes from one public database, so a small cohort's estimate is tested against genomes from the same collections.

**Conditions.** Rules' six conditions are for a change to the release and are not applied here. Conditions 2 (the sign test after Holm), 4 (the recall interval) and 5 (the study-split sign) are not computed. No marker is found, named or ranked.

## Method

The reading program, `scripts/c1-cutoff/read.py`, read each species' two directories: the run's copy and the second scoring. It ran once, at commit `b030d4e`, with `--measures-over shared --line-50 beside --twin-rule beside --dataset-differences yes`. It stops and writes no table unless every gate passes. All passed, and no gate stopped the run.

| Gate | Check |
|---|---|
| 1 | the scorer file has the sha256 given |
| 2 | each species holds `DONE` with the expected code sha256, and its binary and file list agree |
| 3 | every dataset holds its 30 results: the lines, the draws and the twin |
| 4 | each second scoring's files equal its `SHA256SUMS` |
| 5 | each `checks:` line says every result passes every check |
| 6 | the score and line files hold exactly the expected rows, and antibiotics equal the names file |
| 7 | 10 species, 89 datasets and 2,670 results, none twice |
| 8 | each species' `count.tsv` is what the scorer's `--count` writes |

It also checks that each species' line file holds the draws and errors its score file gives. The count is the scorer's: `scripts/c1-cutoff/score.py --count` over the ten species. The reading program copies those cells, lists each dataset's entry and checks it against the scorer's counts, and pools Rules' measures from the score files. A difference, a median and an interval's end are compared at 12 decimals, and coverage in exact fractions.

**Files.** The reading's tables are in `reference/cohort-result/`, with part A's line counts. Each dataset's mean difference is listed there by name, in `page-entry.tsv`. The differences decide nothing.

**Tables.** `scripts/c1-cutoff/page_tables.py` prints this page's tables from those files, rounded to 3 decimals, and a test checks that each printed row stands here.

## Datasets

| Species | Antibiotics | Datasets | Study plans |
|---|---|---:|---:|
| *Enterococcus faecium* | ampicillin, clindamycin, gentamicin, kanamycin, streptomycin, teicoplanin, tetracycline, vancomycin | 8 | 0 |
| *Escherichia coli* | aztreonam, cefepime, cefotaxime, cefoxitin, ceftazidime, ceftriaxone, cefuroxime, chloramphenicol, ertapenem, piperacillin/tazobactam, tetracycline, tobramycin, trimethoprim | 13 | 8 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefoxitin, ceftazidime, ceftriaxone, ciprofloxacin, ertapenem, imipenem, levofloxacin, sulfamethoxazole/trimethoprim, tetracycline, tigecycline, tobramycin | 15 | 14 |
| *Mycobacterium tuberculosis* | amikacin, capreomycin, delamanid, ethambutol, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, ofloxacin, para-aminosalicylic acid, rifabutin, rifampin | 13 | 0 |
| *Neisseria gonorrhoeae* | ciprofloxacin, penicillin, tetracycline | 3 | 0 |
| *Pseudomonas aeruginosa* | ceftazidime, imipenem, levofloxacin, meropenem | 4 | 1 |
| *Salmonella enterica* | amoxicillin/clavulanate, ampicillin, cefotaxime, cefoxitin, ceftiofur, ceftriaxone, chloramphenicol, ciprofloxacin, gentamicin, nalidixic acid, streptomycin, sulfamethoxazole/trimethoprim, tetracycline, trimethoprim | 14 | 0 |
| *Shigella sonnei* | amikacin, ampicillin, azithromycin, ceftriaxone, ciprofloxacin, gentamicin, nalidixic acid | 7 | 0 |
| *Staphylococcus aureus* | cefoxitin, ciprofloxacin, clindamycin, erythromycin, fusidic acid, gentamicin, penicillin, sulfamethoxazole/trimethoprim, tetracycline | 9 | 4 |
| *Streptococcus pneumoniae* | erythromycin, sulfamethoxazole/trimethoprim, tetracycline | 3 | 0 |
| All | | 89 | 27 |

The 27 study plans are those of:

| Species | With a plan |
|---|---|
| *Escherichia coli* | aztreonam, cefepime, cefoxitin, ceftriaxone, ertapenem, piperacillin/tazobactam, tetracycline, tobramycin |
| *Klebsiella pneumoniae* | every antibiotic above but tetracycline |
| *Pseudomonas aeruginosa* | imipenem |
| *Staphylococcus aureus* | clindamycin, erythromycin, penicillin, sulfamethoxazole/trimethoprim |
