# Served-learner result

The release serves, per dataset, the learner with the best random-split balanced accuracy, a tie going to the later learner. This page asks whether a learner picked inside each outer training fold, from inner folds grouped by lineage, serves better on the outer test folds. It reads block C9 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/served-pick.md`.

The reading's rules were fixed at 20:36Z on 22 Nov 2026, before any C9 score was read (`docs/assumption-plan-corrections.md`, entry "C9 design"). The count was read once, from 06:43:51Z to 06:44:30Z on 23 Nov 2026, with no other option, subset or measure.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism, and the inner pick's cannot. Wherever the two are compared on this page, a same or a worse for the inner pick is never read as showing that the release's pick is right (ruling 9).

## Reading

**Result.** The inner pick meets neither of the first two change conditions on the 233 datasets, so the release's pick does not change. Of the plan's six conditions, 3 and 4 are met and 1, 2 and 5 are not. Condition 6 is the maintainer's.

**Condition 1.** On the lineage split the inner pick's median difference in balanced accuracy is -0.001, with a species interval of -0.003 to 0.000. The interval does not lie wholly above 0.

**Condition 2.** The inner pick is better on 37 of 232 datasets, worse on 39 and the same on 156. The one-sided sign test gives p 0.634, and 1.000 after Holm.

**Condition 3.** The twin, the default fitted again at fold seeds 45, 46 and 47, is better on 34 of 232, worse on 36 and the same on 162, with a median of 0.000, an interval of 0.000 to 0.001 and a Holm p of 1.000. It meets neither 1 nor 2, so a change of fold seed alone does not meet the change rule here.

**Condition 4.** The inner pick's recall interval is -0.002 to 0.000. It does not lie wholly below -0.02.

**Condition 5.** The study split's median is 0.001 on 55 datasets, against the lineage split's -0.001. The signs differ.

| Condition | Reads | Read | Ties to the earlier | Ties to the default's | Seed-42 default |
|---|---|---|---|---|---|
| 1 | the lineage median difference's interval lies wholly above 0 | not met | not met | not met | not met |
| 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm | not met | not met | not met | not met |
| 3 | the twin meets neither 1 nor 2 | met | met | met | met |
| 4 | the recall interval does not lie wholly below -0.02 | met | met | met | met |
| 5 | the study median has the lineage median's sign | not met | not met | not met | not met |
| 6 | the maintainer rules the change | the maintainer's | the maintainer's | the maintainer's | the maintainer's |

**Ruling 9.** The default's served figure can carry the pick's optimism, and the inner pick's cannot. So the count says only that the inner pick did not beat the default by the plan's rule. It does not show that the release's pick is right.

**Beside it.** The three answers scored beside the read ones give the same reading: conditions 3 and 4 met, 1, 2 and 5 not. They decide nothing.

**Folds.** The inner pick served another learner than the default's on 2,247 of 3,490 lineage folds and on 190 of 275 study folds. On 170 lineage folds it had no pick and served the default's learner.

The release keeps its rule whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run: 233 datasets in 23 species, the 268 of L, the reports whose lineage split ran, less the 35 in M, which wait for the maintainer's label rule. Of the 233, 55 in 8 species are in S and have a study split. The section Datasets names them all.

**Arms.** The default serves the release's pick on every outer lineage fold, read from the same run's random split at that run's seed. The inner pick serves, per outer fold, the learner with the best balanced accuracy over 3 inner folds of that fold's training genomes, each lineage whole. Both arms are scored with the release's own outer fits.

**Draws.** The default and the inner pick come from the same runs at fold seeds 42, 43 and 44. A dataset's value is the mean of its draws. The twin is the default at seeds 45, 46 and 47.

**Study split.** The 55 datasets in S get one study run each, on their reviewed study plans, at the release's seeds, with no twin. Their default is the seed-42 lineage run's pick. One study run and no twin is a deviation from the plan's Rules, Draws.

**Results.** 1,453 results: 6 lineage results for each of the 233 datasets and 1 study result for each of the 55.

**Not scored.** *Klebsiella pneumoniae* with polymyxin b leaves the lineage count, so 232 of 233 enter it. Its lineage split did not run at seed 43, nor at one of the twin's seeds, and a dataset with a draw unscored is unscored. The run recorded at seed 43: `invalid lineage outer evaluation: fold 4 validation lacks one or both classes`. Its study split ran and is counted.

**Runs.** Each species ran on a machine of its own and wrote its scorer's tables there, under the answers read. The read count is the scorer's count over those 23 tables, copied as the runs wrote them.

## Count

The primary measure is the served learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals.

The median difference is the inner pick's value less the default's, with a 95% interval from 4,000 species resamples read at the 101st and 3,900th values. The one-sided sign test of better against worse runs with Holm over two pairs on the lineage primary measure: the inner pick and the twin, each against the default.

Wherever the inner pick is compared with the default below, the default's figure can carry its pick's optimism (ruling 9).

| Pair | Condition 1 | Condition 2 | Better | Worse | Same | Median | Interval | Sign p | Holm p |
|---|---|---|---|---|---|---|---|---|---|
| inner pick less default | not met | not met | 37 of 232 | 39 of 232 | 156 of 232 | -0.001 | -0.003 to 0.000 | 0.634 | 1.000 |
| twin less default | not met | not met | 34 of 232 | 36 of 232 | 162 of 232 | 0.000 | 0.000 to 0.001 | 0.640 | 1.000 |

Every measure of the count follows. For balanced accuracy, recall and specificity higher is better. For log loss and both error rates lower is better, so better counts a lower value. The median and its interval stay the inner pick's value less the default's.

| Split | Pair | Measure | Better is | Entered | Better | Worse | Same | Median | Interval | Species |
|---|---|---|---|---|---|---|---|---|---|---|
| lineage | inner pick less default | balanced accuracy | higher | 232 of 233 | 37 | 39 | 156 | -0.001 | -0.003 to 0.000 | 23 |
| lineage | inner pick less default | recall | higher | 232 of 233 | 58 | 50 | 124 | 0.000 | -0.002 to 0.000 | 23 |
| lineage | inner pick less default | specificity | higher | 232 of 233 | 30 | 59 | 143 | -0.001 | -0.004 to 0.000 | 23 |
| lineage | inner pick less default | very major error rate | lower | 232 of 233 | 58 | 50 | 124 | 0.000 | 0.000 to 0.002 | 23 |
| lineage | inner pick less default | major error rate | lower | 232 of 233 | 30 | 59 | 143 | 0.001 | 0.000 to 0.004 | 23 |
| lineage | inner pick less default | log loss | lower | 232 of 233 | 81 | 110 | 41 | 0.014 | 0.000 to 0.058 | 23 |
| lineage | twin less default | balanced accuracy | higher | 232 of 233 | 34 | 36 | 162 | 0.000 | 0.000 to 0.001 | 23 |
| lineage | twin less default | recall | higher | 232 of 233 | 51 | 45 | 136 | 0.000 | 0.000 to 0.003 | 23 |
| lineage | twin less default | specificity | higher | 232 of 233 | 38 | 55 | 139 | 0.000 | -0.001 to 0.000 | 23 |
| lineage | twin less default | very major error rate | lower | 232 of 233 | 51 | 45 | 136 | 0.000 | -0.003 to 0.000 | 23 |
| lineage | twin less default | major error rate | lower | 232 of 233 | 38 | 55 | 139 | 0.000 | 0.000 to 0.001 | 23 |
| lineage | twin less default | log loss | lower | 232 of 233 | 78 | 71 | 83 | 0.000 | -0.005 to 0.003 | 23 |
| study | inner pick less default | balanced accuracy | higher | 55 of 55 | 18 | 8 | 29 | 0.001 | -0.003 to 0.012 | 8 |
| study | inner pick less default | recall | higher | 55 of 55 | 21 | 13 | 21 | 0.000 | -0.014 to 0.015 | 8 |
| study | inner pick less default | specificity | higher | 55 of 55 | 17 | 17 | 21 | 0.000 | -0.013 to 0.003 | 8 |
| study | inner pick less default | very major error rate | lower | 55 of 55 | 21 | 13 | 21 | 0.000 | -0.015 to 0.014 | 8 |
| study | inner pick less default | major error rate | lower | 55 of 55 | 17 | 17 | 21 | 0.000 | -0.003 to 0.013 | 8 |
| study | inner pick less default | log loss | lower | 55 of 55 | 17 | 27 | 11 | 0.016 | 0.004 to 0.080 | 8 |

**Recall and log loss.** On the lineage split the inner pick's recall is better on 58 of 232 datasets and worse on 50, with a median of 0.000. Its log loss is better on 81 of 232 and worse on 110, with a median of 0.014, the inner pick's higher.

## Beside it

Three answers of the design's questions 3 and 5 were scored beside the read ones, from the same results, each with its own count. They decide nothing.

| Answer | Question | Changes |
|---|---|---|
| `--tie earlier` | 3. Ties | an inner tie goes to the earlier learner |
| `--tie default` | 3. Ties | an inner tie goes to the default's learner when it is tied |
| `--default seed42` | 5. Default | the default serves the seed-42 pick at every draw |

Each answer meets conditions 3 and 4 and none of 1, 2 and 5, as the conditions table under Reading gives. Ruling 9 holds under each of them: the default's figure can carry its pick's optimism.

| Answer | Split | Pair | Measure | Better | Worse | Same | Median | Interval | Conditions met |
|---|---|---|---|---|---|---|---|---|---|
| read: ties to the later learner, each draw's own default | lineage | inner pick less default | balanced accuracy | 37 | 39 | 156 | -0.001 | -0.003 to 0.000 | 3, 4 |
| read: ties to the later learner, each draw's own default | lineage | inner pick less default | recall | 58 | 50 | 124 | 0.000 | -0.002 to 0.000 | 3, 4 |
| read: ties to the later learner, each draw's own default | lineage | twin less default | balanced accuracy | 34 | 36 | 162 | 0.000 | 0.000 to 0.001 | 3, 4 |
| read: ties to the later learner, each draw's own default | lineage | twin less default | recall | 51 | 45 | 136 | 0.000 | 0.000 to 0.003 | 3, 4 |
| read: ties to the later learner, each draw's own default | study | inner pick less default | balanced accuracy | 18 | 8 | 29 | 0.001 | -0.003 to 0.012 | 3, 4 |
| read: ties to the later learner, each draw's own default | study | inner pick less default | recall | 21 | 13 | 21 | 0.000 | -0.014 to 0.015 | 3, 4 |
| beside: ties to the earlier learner | lineage | inner pick less default | balanced accuracy | 38 | 37 | 157 | 0.000 | -0.002 to 0.000 | 3, 4 |
| beside: ties to the earlier learner | lineage | inner pick less default | recall | 59 | 49 | 124 | 0.000 | -0.001 to 0.000 | 3, 4 |
| beside: ties to the earlier learner | lineage | twin less default | balanced accuracy | 34 | 36 | 162 | 0.000 | 0.000 to 0.001 | 3, 4 |
| beside: ties to the earlier learner | lineage | twin less default | recall | 51 | 45 | 136 | 0.000 | 0.000 to 0.003 | 3, 4 |
| beside: ties to the earlier learner | study | inner pick less default | balanced accuracy | 18 | 8 | 29 | 0.001 | -0.003 to 0.012 | 3, 4 |
| beside: ties to the earlier learner | study | inner pick less default | recall | 20 | 13 | 22 | 0.000 | -0.014 to 0.013 | 3, 4 |
| beside: ties to the default's learner | lineage | inner pick less default | balanced accuracy | 37 | 37 | 158 | 0.000 | -0.002 to 0.000 | 3, 4 |
| beside: ties to the default's learner | lineage | inner pick less default | recall | 58 | 49 | 125 | 0.000 | -0.001 to 0.000 | 3, 4 |
| beside: ties to the default's learner | lineage | twin less default | balanced accuracy | 34 | 36 | 162 | 0.000 | 0.000 to 0.001 | 3, 4 |
| beside: ties to the default's learner | lineage | twin less default | recall | 51 | 45 | 136 | 0.000 | 0.000 to 0.003 | 3, 4 |
| beside: ties to the default's learner | study | inner pick less default | balanced accuracy | 18 | 8 | 29 | 0.001 | -0.003 to 0.012 | 3, 4 |
| beside: ties to the default's learner | study | inner pick less default | recall | 20 | 13 | 22 | 0.000 | -0.014 to 0.013 | 3, 4 |
| beside: the seed-42 default at every draw | lineage | inner pick less default | balanced accuracy | 35 | 50 | 147 | -0.002 | -0.003 to 0.000 | 3, 4 |
| beside: the seed-42 default at every draw | lineage | inner pick less default | recall | 59 | 58 | 115 | 0.000 | -0.002 to 0.000 | 3, 4 |
| beside: the seed-42 default at every draw | lineage | twin less default | balanced accuracy | 31 | 26 | 175 | 0.000 | 0.000 to 0.000 | 3, 4 |
| beside: the seed-42 default at every draw | lineage | twin less default | recall | 42 | 42 | 148 | 0.000 | 0.000 to 0.000 | 3, 4 |
| beside: the seed-42 default at every draw | study | inner pick less default | balanced accuracy | 18 | 8 | 29 | 0.001 | -0.003 to 0.012 | 3, 4 |
| beside: the seed-42 default at every draw | study | inner pick less default | recall | 21 | 13 | 21 | 0.000 | -0.014 to 0.015 | 3, 4 |

**Pooled calls.** The calls summed over the scored datasets and their draws give Rules' measures with their Wilson intervals. A genome counts once a draw. They decide nothing.

| Split | Arm | Datasets scored | Balanced accuracy | Recall | Recall interval | Specificity | Very major error rate | Major error rate |
|---|---|---|---|---|---|---|---|---|
| lineage | default | 232 of 233 | 0.874 | 0.836 (184,408 of 220,656) | 0.834 to 0.837 | 0.913 (264,703 of 290,004) | 0.164 (36,248 of 220,656) | 0.087 (25,301 of 290,004) |
| lineage | inner | 232 of 233 | 0.874 | 0.840 (185,400 of 220,656) | 0.839 to 0.842 | 0.908 (263,291 of 290,004) | 0.160 (35,256 of 220,656) | 0.092 (26,713 of 290,004) |
| lineage | twin | 232 of 233 | 0.873 | 0.838 (184,991 of 220,656) | 0.837 to 0.840 | 0.908 (263,319 of 290,004) | 0.162 (35,665 of 220,656) | 0.092 (26,685 of 290,004) |
| study | default | 55 of 55 | 0.844 | 0.828 (16,265 of 19,650) | 0.822 to 0.833 | 0.861 (16,353 of 19,002) | 0.172 (3,385 of 19,650) | 0.139 (2,649 of 19,002) |
| study | inner | 55 of 55 | 0.851 | 0.843 (16,556 of 19,650) | 0.837 to 0.848 | 0.860 (16,349 of 19,002) | 0.157 (3,094 of 19,650) | 0.140 (2,653 of 19,002) |

**Picks.** The folds each arm serves to each learner, over the scored draws. A fold with no pick serves the default's learner and is counted under its learner and under No pick.

| Answer | Split | Serves | Draws | Folds | `scm` | `cart` | `rf` | `gbt` | `lasso` | No pick | Tied | Same as default | Differs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| read: ties to the later learner, each draw's own default | lineage | default | 698 | 3,490 | 820 | 530 | 455 | 790 | 895 |  |  |  |  |
| read: ties to the later learner, each draw's own default | lineage | inner | 698 | 3,490 | 903 | 755 | 463 | 563 | 806 | 170 | 486 | 1,073 | 2,247 of 3,490 |
| read: ties to the later learner, each draw's own default | study | default | 55 | 275 | 65 | 35 | 40 | 65 | 70 |  |  |  |  |
| read: ties to the later learner, each draw's own default | study | inner | 55 | 275 | 77 | 46 | 32 | 58 | 62 | 0 | 15 | 85 | 190 of 275 |
| beside: ties to the earlier learner | lineage | default | 698 | 3,490 | 820 | 530 | 455 | 790 | 895 |  |  |  |  |
| beside: ties to the earlier learner | lineage | inner | 698 | 3,490 | 1,325 | 590 | 480 | 508 | 587 | 170 | 486 | 978 | 2,342 of 3,490 |
| beside: ties to the earlier learner | study | default | 55 | 275 | 65 | 35 | 40 | 65 | 70 |  |  |  |  |
| beside: ties to the earlier learner | study | inner | 55 | 275 | 88 | 44 | 31 | 58 | 54 | 0 | 15 | 85 | 190 of 275 |
| beside: ties to the default's learner | lineage | default | 698 | 3,490 | 820 | 530 | 455 | 790 | 895 |  |  |  |  |
| beside: ties to the default's learner | lineage | inner | 698 | 3,490 | 948 | 783 | 475 | 556 | 728 | 170 | 486 | 1,195 | 2,125 of 3,490 |
| beside: ties to the default's learner | study | default | 55 | 275 | 65 | 35 | 40 | 65 | 70 |  |  |  |  |
| beside: ties to the default's learner | study | inner | 55 | 275 | 79 | 50 | 31 | 59 | 56 | 0 | 15 | 92 | 183 of 275 |
| beside: the seed-42 default at every draw | lineage | default | 698 | 3,490 | 685 | 585 | 405 | 840 | 975 |  |  |  |  |
| beside: the seed-42 default at every draw | lineage | inner | 698 | 3,490 | 908 | 748 | 445 | 569 | 820 | 170 | 486 | 1,086 | 2,234 of 3,490 |
| beside: the seed-42 default at every draw | study | default | 55 | 275 | 65 | 35 | 40 | 65 | 70 |  |  |  |  |
| beside: the seed-42 default at every draw | study | inner | 55 | 275 | 77 | 46 | 32 | 58 | 62 | 0 | 15 | 85 | 190 of 275 |

For the answers read, the folds the inner pick serves to each learner, by the default's learner of the draw:

| Split | Default's learner | Draws | Folds | No pick | Inner `scm` | Inner `cart` | Inner `rf` | Inner `gbt` | Inner `lasso` |
|---|---|---|---|---|---|---|---|---|---|
| lineage | `scm` | 164 | 820 | 39 | 375 | 178 | 60 | 71 | 136 |
| lineage | `cart` | 106 | 530 | 16 | 111 | 191 | 42 | 68 | 118 |
| lineage | `rf` | 91 | 455 | 27 | 77 | 89 | 151 | 62 | 76 |
| lineage | `gbt` | 158 | 790 | 24 | 166 | 139 | 123 | 206 | 156 |
| lineage | `lasso` | 179 | 895 | 64 | 174 | 158 | 87 | 156 | 320 |
| study | `scm` | 13 | 65 | 0 | 33 | 10 | 5 | 6 | 11 |
| study | `cart` | 7 | 35 | 0 | 7 | 10 | 2 | 6 | 10 |
| study | `rf` | 8 | 40 | 0 | 7 | 5 | 8 | 12 | 8 |
| study | `gbt` | 13 | 65 | 0 | 8 | 12 | 14 | 16 | 15 |
| study | `lasso` | 14 | 70 | 0 | 22 | 9 | 3 | 18 | 18 |

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above, and decide nothing. Where they compare the inner pick with the default, the default's figure can carry its pick's optimism (ruling 9).

**Seed alone.** The twin moves 27 of 232 lineage datasets by more than 0.05 in balanced accuracy, from -0.192 (*Klebsiella oxytoca* with cefpodoxime) to 0.100 (*Enterobacter cloacae* with imipenem). The inner pick moves 28 of 232 by more than 0.05, from -0.132 (*Salmonella enterica* with spectinomycin) to 0.134 (*Klebsiella pneumoniae* with ofloxacin). A single dataset's difference is of the size a fold seed moves.

**Specificity.** On the lineage split the inner pick's specificity is better on 30 of 232 datasets and worse on 59, with a median of -0.001 and an interval of -0.004 to 0.000. Its recall is better on 58 of 232 and worse on 50.

**Pooled.** Pooled over datasets and draws, the inner pick's specificity is 0.908 (263,291 of 290,004 susceptible genome calls, a genome once a draw) against the default's 0.913 (264,703 of 290,004), and the twin's is 0.908 as well. Its pooled recall is 0.840 (185,400 of 220,656 resistant genome calls) against the default's 0.836 (184,408 of 220,656), the twin's 0.838. These do not show either pick right.

**Study split.** The study split's sign test alone, outside the Holm family, gives p 0.038 for 18 better against 8 worse of 55. The plan reads only the study median's sign, and the study split has one run per dataset and no twin.

**Species.** The species medians on the lineage split run from -0.036 (*Staphylococcus haemolyticus*, 3 datasets) to 0.022 (*Klebsiella aerogenes*, 3 datasets). The release does not change on one species.

## Limits

**Ruling 9.** The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism, and the inner pick's cannot. A same or a worse for the inner pick is never read as showing that the release's pick is right.

**Study split.** One study run per dataset at the release's seeds, and no study twin. That departs from the plan's Rules, Draws, which give each arm seeds 42, 43 and 44.

**Binaries.** 22 species ran on the binary built from main at `0fd05e9`. *Klebsiella pneumoniae*'s first run stopped, and it ran again on the binary of main at `dd17f98`, which writes the result of a draw whose lineage split cannot be made, with its reason.

Under `src/` and `scripts/c9-pick/` the two commits differ only in `src/learn/inner_pick.rs`, and `score.py` is the same file in both. The first run's 32 seed-42 results equal the second run's in 32 of 32, with the timing field `seconds_inner` left out. Each species' results come from one run on one binary.

**Labels.** The 35 datasets of M wait for the maintainer's label rule. The count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged. On the 233, the label rules `r-or-s-wins` (the release's), `intermediate-aside` and `intermediate-class` give the same labels. 32 of the 233 hold a genome pair of S and R calls, and their measures depend on the most resistant call's rule against the other three (`docs/served-pick.md`).

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits and learners only.

## Method

**Checks of the run.** Before any score was read, a check program read each species' copy and printed counts only. All six passed. `checks.txt` in `reference/served-pick/` holds its printout.

| Check | Result |
|---|---|
| species | 23 copies, the 23 species of the first run |
| `SHA256SUMS` | holds in 23 of 23 species, 1,683 of 1,683 lines |
| `DONE` and `binary.txt` | 23 of 23 name their species, `SHA256SUMS`, result count, binary and code, and none is `FAILED`; 22 species on `0fd05e9`, 1 on `dd17f98` |
| `expect.py` | 23 of 23 `expect.tsv` files are `expect.py`'s from `scripts/c9-pick/l.tsv`; 1,453 of 1,453 results found and none more |
| genome ids | 0 in 2,671 files, against the 23,757 genome ids of the release's 23 species lists |
| *Klebsiella pneumoniae* | 32 of 32 first-run seed-42 results equal the second run's, `seconds_inner` left out |

**Reading.** The answers read are the run's own: each species' `score/` tables, copied as its run wrote them. The count ran once over them, and then each answer beside was scored per species from the same results and counted:

```
python3 -I scripts/c9-pick/score.py --count read/count read/score/<species>/datasets.tsv      # 23 files
python3 -I scripts/c9-pick/score.py <copy>/expect.tsv <copy>/results/picks <release> tie-earlier/score/<species> --tie earlier
python3 -I scripts/c9-pick/score.py <copy>/expect.tsv <copy>/results/picks <release> tie-default/score/<species> --tie default
python3 -I scripts/c9-pick/score.py <copy>/expect.tsv <copy>/results/picks <release> default-seed42/score/<species> --default seed42
python3 -I scripts/c9-pick/score.py --count <answer>/count <answer>/score/<species>/datasets.tsv   # each answer beside
```

Each count printed `count: 809 dataset rows from 23 files`: 3 lineage rows for each of the 233 datasets and 2 study rows for each of the 55. Each of the 69 scorings beside passed every check of every result. `reading.txt` lists every command with the sha256 of each file it wrote.

**Tables.** `scripts/c9-pick/page_tables.py write` checked the reading against eight gates and wrote `reference/served-pick/`. All eight passed:

| Gate | What passed |
|---|---|
| 1. answers | the reading holds reading.txt and the four answers, each with count/count.tsv, count/conditions.tsv and score/ |
| 2. species | each answer's score/ holds exactly the species of l.tsv's first run, each with its four files |
| 3. options | each species' options.txt holds its answer's options |
| 4. datasets | each datasets.tsv holds exactly its species' first-run datasets, none in M: the default, the inner pick and the twin on the lineage split, and the default and the inner pick on the study split for those in S; species and antibiotic equal l.tsv's |
| 5. draws | each draws.tsv holds each dataset's draws at seeds 42, 43 and 44, its twin at 45, 46 and 47 and its study run at seed 42, and nothing else; datasets.tsv is their mean; picks.tsv is their inner rows |
| 6. beside | the tie answers' default and twin rows equal the read answer's |
| 7. count | each answer's count.tsv and conditions.tsv are what score.py --count writes from its datasets.tsv files |
| 8. conditions | the inner pick's conditions 1 and 2 as read here off count.tsv, and condition 3 off the twin's row, equal conditions.tsv's |

A count of the release's 23,757 genome ids over the reading and those tables found 0 in 488 files (`checks.txt`).

**Tests.** `tests/c9_pick/test_page_tables.py` holds 20 SYNTHETIC tests of the table script, on three made-up species, and `tests/c9_pick/test_page.py` holds 14 tests of this page against `reference/served-pick/`. `scripts/c9-pick/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it. All 7 were caught.

| Fault | Caught by |
|---|---|
| 7. A twin row read as the inner pick's | `Tables.test_the_inner_figure_is_the_inner_picks_and_the_twin_stands_apart` |
| 8. A beside answer read as the read one | `Tables.test_the_read_answer_is_the_reading_and_the_beside_answers_stand_apart` |
| 8b. A read folder scored under a beside answer's options passed | `Gates.test_an_answer_scored_under_other_options_is_refused` |
| 9. A species dropped | `Gates.test_a_species_missing_from_the_reading_is_refused` |
| 10. A dataset of M counted | `Gates.test_a_dataset_of_m_is_refused` |
| 11. A study row counted on the lineage split | `Tables.test_a_study_row_stays_on_the_study_split` |
| 12. A missing draw read as zero | `Tables.test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero` |

**Files.** The reading's tables are in `reference/served-pick/`. `scripts/c9-pick/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals, and `tests/c9_pick/test_page.py` checks that each printed row stands here.

## Datasets

Per species, the datasets of the first run, those in S, and better, worse and same on the primary measure by split and pair. A species' median decides nothing.

| Species | Datasets | In S | Lineage inner better / worse / same | Lineage inner median | Twin better / worse / same | Twin median | Study inner better / worse / same | Study inner median |
|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | 7 | 2 | 0 / 3 / 4 | -0.019 | 1 / 3 / 3 | 0.000 | 1 / 0 / 1 | 0.025 |
| *Campylobacter coli* | 8 | 0 | 0 / 1 / 7 | 0.000 | 0 / 1 / 7 | 0.000 | 0 / 0 / 0 | none |
| *Campylobacter jejuni* | 4 | 0 | 0 / 1 / 3 | -0.001 | 1 / 0 / 3 | 0.005 | 0 / 0 / 0 | none |
| *Clostridioides difficile* | 4 | 0 | 0 / 0 / 4 | -0.009 | 0 / 0 / 4 | -0.004 | 0 / 0 / 0 | none |
| *Corynebacterium diphtheriae* | 3 | 0 | 0 / 1 / 2 | -0.003 | 0 / 0 / 3 | 0.000 | 0 / 0 / 0 | none |
| *Enterobacter cloacae* | 9 | 8 | 4 / 2 / 3 | 0.010 | 2 / 4 / 3 | -0.011 | 4 / 3 / 1 | 0.018 |
| *Enterococcus faecium* | 11 | 0 | 3 / 1 / 7 | 0.002 | 2 / 0 / 9 | 0.000 | 0 / 0 / 0 | none |
| *Escherichia coli* | 29 | 15 | 5 / 5 / 19 | -0.001 | 2 / 4 / 23 | 0.001 | 3 / 3 / 9 | 0.000 |
| *Klebsiella aerogenes* | 3 | 0 | 2 / 1 / 0 | 0.022 | 1 / 1 / 1 | 0.000 | 0 / 0 / 0 | none |
| *Klebsiella oxytoca* | 4 | 0 | 1 / 1 / 2 | 0.010 | 2 / 1 / 1 | 0.002 | 0 / 0 / 0 | none |
| *Klebsiella pneumoniae* | 32 | 18 | 6 / 3 / 22 | 0.000 | 5 / 7 / 19 | -0.001 | 8 / 1 / 9 | 0.013 |
| *Klebsiella quasipneumoniae* | 3 | 0 | 0 / 0 / 3 | -0.003 | 0 / 0 / 3 | 0.011 | 0 / 0 / 0 | none |
| *Mycobacterium tuberculosis* | 20 | 1 | 0 / 4 / 16 | -0.005 | 3 / 3 / 14 | 0.000 | 0 / 0 / 1 | -0.012 |
| *Neisseria gonorrhoeae* | 6 | 0 | 2 / 1 / 3 | 0.000 | 1 / 0 / 5 | -0.002 | 0 / 0 / 0 | none |
| *Pseudomonas aeruginosa* | 12 | 6 | 2 / 2 / 8 | -0.002 | 3 / 1 / 8 | 0.000 | 2 / 1 / 3 | -0.003 |
| *Salmonella enterica* | 24 | 1 | 6 / 4 / 14 | 0.000 | 7 / 5 / 12 | -0.001 | 0 / 0 / 1 | -0.006 |
| *Shigella flexneri* | 4 | 0 | 1 / 1 / 2 | 0.003 | 1 / 0 / 3 | 0.003 | 0 / 0 / 0 | none |
| *Shigella sonnei* | 20 | 0 | 0 / 4 / 16 | -0.006 | 1 / 1 / 18 | 0.002 | 0 / 0 / 0 | none |
| *Staphylococcus aureus* | 13 | 4 | 4 / 0 / 9 | 0.000 | 1 / 1 / 11 | -0.003 | 0 / 0 / 4 | -0.002 |
| *Staphylococcus haemolyticus* | 3 | 0 | 0 / 2 / 1 | -0.036 | 0 / 0 / 3 | -0.001 | 0 / 0 / 0 | none |
| *Streptococcus agalactiae* | 3 | 0 | 0 / 0 / 3 | -0.002 | 0 / 1 / 2 | -0.003 | 0 / 0 / 0 | none |
| *Streptococcus pneumoniae* | 8 | 0 | 1 / 2 / 5 | -0.009 | 0 / 2 / 6 | 0.001 | 0 / 0 / 0 | none |
| *Vibrio cholerae* | 3 | 0 | 0 / 0 / 3 | 0.007 | 1 / 1 / 1 | 0.010 | 0 / 0 / 0 | none |

Per dataset and split, with the default's learner at each scored draw in seed order. The default's figure can carry its pick's optimism (ruling 9). Each difference is the inner pick's value less the default's, and the twin's difference is the twin's value less the default's. A dataset with a draw unscored has no figure.

| Species | Antibiotic | Split | Default's learner | Default | Inner pick | Difference | Twin difference | Draws scored | No pick | Folds differing |
|---|---|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | cefepime | lineage | `cart`, `scm`, `scm` | 0.574 | 0.555 | -0.019 | -0.049 | 3 of 3 | 0 | 9 of 15 |
| *Acinetobacter baumannii* | cefepime | study | `cart` | 0.696 | 0.696 | 0.000 |  | 1 of 1 | 0 | 1 of 5 |
| *Acinetobacter baumannii* | cefotaxime | lineage | `scm`, `scm`, `scm` | 0.658 | 0.665 | 0.007 | 0.011 | 3 of 3 | 0 | 11 of 15 |
| *Acinetobacter baumannii* | doripenem | lineage | `gbt`, `gbt`, `lasso` | 0.711 | 0.722 | 0.011 | -0.078 | 3 of 3 | 0 | 14 of 15 |
| *Acinetobacter baumannii* | doxycycline | lineage | `rf`, `rf`, `rf` | 0.594 | 0.531 | -0.063 | 0.000 | 3 of 3 | 0 | 9 of 15 |
| *Acinetobacter baumannii* | minocycline | lineage | `lasso`, `gbt`, `gbt` | 0.505 | 0.482 | -0.023 | 0.038 | 3 of 3 | 0 | 11 of 15 |
| *Acinetobacter baumannii* | minocycline | study | `lasso` | 0.503 | 0.553 | 0.050 |  | 1 of 1 | 0 | 4 of 5 |
| *Acinetobacter baumannii* | moxifloxacin | lineage | `cart`, `lasso`, `cart` | 0.959 | 0.901 | -0.058 | 0.000 | 3 of 3 | 0 | 9 of 15 |
| *Acinetobacter baumannii* | sulfamethoxazole/trimethoprim | lineage | `cart`, `cart`, `cart` | 0.773 | 0.775 | 0.002 | -0.063 | 3 of 3 | 0 | 12 of 15 |
| *Campylobacter coli* | azithromycin | lineage | `lasso`, `gbt`, `cart` | 0.981 | 0.983 | 0.002 | 0.000 | 3 of 3 | 0 | 10 of 15 |
| *Campylobacter coli* | ciprofloxacin | lineage | `lasso`, `lasso`, `lasso` | 1.000 | 1.000 | 0.000 | 0.000 | 3 of 3 | 0 | 0 of 15 |
| *Campylobacter coli* | clindamycin | lineage | `scm`, `scm`, `scm` | 0.905 | 0.894 | -0.011 | 0.001 | 3 of 3 | 0 | 11 of 15 |
| *Campylobacter coli* | erythromycin | lineage | `cart`, `lasso`, `gbt` | 0.989 | 0.989 | 0.000 | -0.002 | 3 of 3 | 0 | 10 of 15 |
| *Campylobacter coli* | gentamicin | lineage | `scm`, `scm`, `scm` | 0.956 | 0.885 | -0.071 | -0.052 | 3 of 3 | 0 | 9 of 15 |
| *Campylobacter coli* | nalidixic acid | lineage | `lasso`, `lasso`, `cart` | 0.984 | 0.984 | 0.000 | 0.003 | 3 of 3 | 0 | 11 of 15 |
| *Campylobacter coli* | telithromycin | lineage | `lasso`, `scm`, `lasso` | 0.977 | 0.976 | -0.001 | -0.001 | 3 of 3 | 0 | 8 of 15 |
| *Campylobacter coli* | tetracycline | lineage | `gbt`, `rf`, `rf` | 0.981 | 0.981 | 0.000 | -0.001 | 3 of 3 | 0 | 7 of 15 |
| *Campylobacter jejuni* | ciprofloxacin | lineage | `cart`, `cart`, `lasso` | 0.985 | 0.981 | -0.004 | 0.000 | 3 of 3 | 0 | 8 of 15 |
| *Campylobacter jejuni* | clindamycin | lineage | `scm`, `lasso`, `lasso` | 0.773 | 0.728 | -0.044 | 0.053 | 3 of 3 | 0 | 8 of 15 |
| *Campylobacter jejuni* | nalidixic acid | lineage | `cart`, `cart`, `lasso` | 0.978 | 0.984 | 0.006 | 0.006 | 3 of 3 | 0 | 8 of 15 |
| *Campylobacter jejuni* | tetracycline | lineage | `lasso`, `lasso`, `lasso` | 0.980 | 0.982 | 0.003 | 0.004 | 3 of 3 | 0 | 4 of 15 |
| *Clostridioides difficile* | azithromycin | lineage | `gbt`, `lasso`, `scm` | 0.914 | 0.912 | -0.001 | 0.000 | 3 of 3 | 0 | 11 of 15 |
| *Clostridioides difficile* | ceftriaxone | lineage | `cart`, `scm`, `cart` | 0.919 | 0.901 | -0.018 | -0.016 | 3 of 3 | 0 | 7 of 15 |
| *Clostridioides difficile* | clarithromycin | lineage | `lasso`, `gbt`, `cart` | 0.911 | 0.894 | -0.017 | -0.003 | 3 of 3 | 0 | 12 of 15 |
| *Clostridioides difficile* | moxifloxacin | lineage | `gbt`, `lasso`, `cart` | 0.982 | 0.984 | 0.002 | -0.004 | 3 of 3 | 0 | 13 of 15 |
| *Corynebacterium diphtheriae* | sulfamethoxazole/trimethoprim | lineage | `cart`, `cart`, `cart` | 0.924 | 0.924 | 0.000 | 0.000 | 3 of 3 | 0 | 4 of 15 |
| *Corynebacterium diphtheriae* | tetracycline | lineage | `lasso`, `lasso`, `lasso` | 0.907 | 0.904 | -0.003 | -0.017 | 3 of 3 | 0 | 3 of 15 |
| *Corynebacterium diphtheriae* | trimethoprim | lineage | `lasso`, `lasso`, `scm` | 0.749 | 0.678 | -0.071 | 0.013 | 3 of 3 | 0 | 8 of 15 |
| *Enterobacter cloacae* | cefepime | lineage | `cart`, `lasso`, `cart` | 0.820 | 0.858 | 0.038 | 0.012 | 3 of 3 | 9 | 5 of 15 |
| *Enterobacter cloacae* | cefepime | study | `cart` | 0.869 | 0.875 | 0.006 |  | 1 of 1 | 0 | 3 of 5 |
| *Enterobacter cloacae* | ciprofloxacin | lineage | `gbt`, `lasso`, `lasso` | 0.669 | 0.719 | 0.051 | -0.011 | 3 of 3 | 9 | 5 of 15 |
| *Enterobacter cloacae* | ciprofloxacin | study | `gbt` | 0.852 | 0.746 | -0.106 |  | 1 of 1 | 0 | 5 of 5 |
| *Enterobacter cloacae* | gentamicin | lineage | `gbt`, `scm`, `lasso` | 0.671 | 0.647 | -0.024 | -0.060 | 3 of 3 | 0 | 10 of 15 |
| *Enterobacter cloacae* | gentamicin | study | `gbt` | 0.634 | 0.663 | 0.029 |  | 1 of 1 | 0 | 4 of 5 |
| *Enterobacter cloacae* | imipenem | lineage | `cart`, `rf`, `cart` | 0.472 | 0.550 | 0.078 | 0.100 | 3 of 3 | 0 | 11 of 15 |
| *Enterobacter cloacae* | imipenem | study | `cart` | 0.551 | 0.622 | 0.071 |  | 1 of 1 | 0 | 4 of 5 |
| *Enterobacter cloacae* | levofloxacin | lineage | `lasso`, `lasso`, `scm` | 0.542 | 0.552 | 0.010 | -0.045 | 3 of 3 | 9 | 3 of 15 |
| *Enterobacter cloacae* | levofloxacin | study | `lasso` | 0.642 | 0.764 | 0.122 |  | 1 of 1 | 0 | 3 of 5 |
| *Enterobacter cloacae* | meropenem | lineage | `lasso`, `cart`, `lasso` | 0.498 | 0.507 | 0.010 | 0.011 | 3 of 3 | 0 | 11 of 15 |
| *Enterobacter cloacae* | meropenem | study | `lasso` | 0.623 | 0.669 | 0.046 |  | 1 of 1 | 0 | 2 of 5 |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | lineage | `cart`, `cart`, `cart` | 0.974 | 0.974 | 0.000 | -0.020 | 3 of 3 | 0 | 4 of 15 |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | study | `cart` | 0.974 | 0.899 | -0.076 |  | 1 of 1 | 0 | 4 of 5 |
| *Enterobacter cloacae* | tobramycin | lineage | `gbt`, `gbt`, `gbt` | 0.801 | 0.900 | 0.099 | 0.060 | 3 of 3 | 0 | 12 of 15 |
| *Enterobacter cloacae* | tobramycin | study | `gbt` | 0.959 | 0.876 | -0.083 |  | 1 of 1 | 0 | 5 of 5 |
| *Enterobacter cloacae* | trimethoprim | lineage | `cart`, `cart`, `cart` | 0.904 | 0.857 | -0.047 | -0.034 | 3 of 3 | 0 | 8 of 15 |
| *Enterococcus faecium* | ampicillin | lineage | `gbt`, `lasso`, `scm` | 0.988 | 0.989 | 0.001 | 0.000 | 3 of 3 | 0 | 13 of 15 |
| *Enterococcus faecium* | clindamycin | lineage | `lasso`, `lasso`, `gbt` | 0.619 | 0.727 | 0.108 | 0.060 | 3 of 3 | 0 | 10 of 15 |
| *Enterococcus faecium* | dalfopristin/quinupristin | lineage | `cart`, `cart`, `scm` | 0.632 | 0.535 | -0.098 | -0.004 | 3 of 3 | 0 | 8 of 15 |
| *Enterococcus faecium* | doxycycline | lineage | `gbt`, `lasso`, `gbt` | 0.890 | 0.892 | 0.002 | 0.000 | 3 of 3 | 0 | 12 of 15 |
| *Enterococcus faecium* | erythromycin | lineage | `scm`, `cart`, `lasso` | 0.791 | 0.824 | 0.033 | 0.046 | 3 of 3 | 0 | 11 of 15 |
| *Enterococcus faecium* | gentamicin | lineage | `gbt`, `rf`, `rf` | 0.884 | 0.920 | 0.036 | 0.015 | 3 of 3 | 0 | 15 of 15 |
| *Enterococcus faecium* | kanamycin | lineage | `cart`, `cart`, `cart` | 0.888 | 0.892 | 0.004 | 0.003 | 3 of 3 | 0 | 7 of 15 |
| *Enterococcus faecium* | streptomycin | lineage | `rf`, `rf`, `rf` | 0.925 | 0.908 | -0.017 | -0.012 | 3 of 3 | 0 | 15 of 15 |
| *Enterococcus faecium* | teicoplanin | lineage | `gbt`, `lasso`, `gbt` | 0.959 | 0.950 | -0.010 | -0.003 | 3 of 3 | 0 | 11 of 15 |
| *Enterococcus faecium* | tetracycline | lineage | `gbt`, `gbt`, `rf` | 0.921 | 0.924 | 0.003 | 0.002 | 3 of 3 | 0 | 12 of 15 |
| *Enterococcus faecium* | vancomycin | lineage | `scm`, `scm`, `scm` | 0.987 | 0.979 | -0.008 | -0.002 | 3 of 3 | 0 | 8 of 15 |
| *Escherichia coli* | amikacin | lineage | `lasso`, `scm`, `lasso` | 0.752 | 0.746 | -0.006 | 0.012 | 3 of 3 | 12 | 1 of 15 |
| *Escherichia coli* | amikacin | study | `lasso` | 0.737 | 0.734 | -0.002 |  | 1 of 1 | 0 | 2 of 5 |
| *Escherichia coli* | ampicillin/sulbactam | lineage | `lasso`, `gbt`, `rf` | 0.643 | 0.671 | 0.028 | 0.006 | 3 of 3 | 0 | 11 of 15 |
| *Escherichia coli* | ampicillin/sulbactam | study | `lasso` | 0.671 | 0.621 | -0.050 |  | 1 of 1 | 0 | 5 of 5 |
| *Escherichia coli* | aztreonam | lineage | `lasso`, `gbt`, `gbt` | 0.856 | 0.881 | 0.025 | -0.015 | 3 of 3 | 0 | 13 of 15 |
| *Escherichia coli* | aztreonam | study | `lasso` | 0.868 | 0.851 | -0.017 |  | 1 of 1 | 0 | 2 of 5 |
| *Escherichia coli* | cefazolin | lineage | `lasso`, `cart`, `lasso` | 0.772 | 0.754 | -0.018 | -0.076 | 3 of 3 | 0 | 13 of 15 |
| *Escherichia coli* | cefazolin | study | `lasso` | 0.735 | 0.768 | 0.033 |  | 1 of 1 | 0 | 3 of 5 |
| *Escherichia coli* | cefepime | lineage | `gbt`, `scm`, `scm` | 0.825 | 0.733 | -0.092 | 0.001 | 3 of 3 | 0 | 10 of 15 |
| *Escherichia coli* | cefepime | study | `gbt` | 0.850 | 0.824 | -0.025 |  | 1 of 1 | 0 | 5 of 5 |
| *Escherichia coli* | cefotaxime | lineage | `scm`, `scm`, `scm` | 0.918 | 0.915 | -0.002 | 0.005 | 3 of 3 | 0 | 6 of 15 |
| *Escherichia coli* | cefoxitin | lineage | `rf`, `lasso`, `rf` | 0.855 | 0.845 | -0.009 | -0.010 | 3 of 3 | 0 | 9 of 15 |
| *Escherichia coli* | cefoxitin | study | `rf` | 0.834 | 0.835 | 0.001 |  | 1 of 1 | 0 | 4 of 5 |
| *Escherichia coli* | cefpodoxime | lineage | `scm`, `cart`, `scm` | 0.736 | 0.743 | 0.006 | -0.052 | 3 of 3 | 0 | 11 of 15 |
| *Escherichia coli* | cefpodoxime/clavulanate | lineage | `scm`, `scm`, `scm` | 0.613 | 0.507 | -0.106 | -0.068 | 3 of 3 | 0 | 14 of 15 |
| *Escherichia coli* | ceftazidime | lineage | `gbt`, `gbt`, `rf` | 0.897 | 0.903 | 0.006 | -0.005 | 3 of 3 | 0 | 11 of 15 |
| *Escherichia coli* | ceftiofur | lineage | `scm`, `lasso`, `lasso` | 0.918 | 0.921 | 0.004 | 0.024 | 3 of 3 | 0 | 11 of 15 |
| *Escherichia coli* | ceftriaxone | lineage | `scm`, `scm`, `scm` | 0.921 | 0.927 | 0.006 | 0.018 | 3 of 3 | 0 | 8 of 15 |
| *Escherichia coli* | ceftriaxone | study | `scm` | 0.906 | 0.912 | 0.005 |  | 1 of 1 | 0 | 3 of 5 |
| *Escherichia coli* | cefuroxime | lineage | `rf`, `gbt`, `lasso` | 0.876 | 0.900 | 0.024 | -0.001 | 3 of 3 | 0 | 12 of 15 |
| *Escherichia coli* | cephalothin | lineage | `lasso`, `lasso`, `lasso` | 0.908 | 0.899 | -0.009 | 0.002 | 3 of 3 | 0 | 12 of 15 |
| *Escherichia coli* | chloramphenicol | lineage | `gbt`, `lasso`, `scm` | 0.890 | 0.925 | 0.035 | 0.005 | 3 of 3 | 0 | 9 of 15 |
| *Escherichia coli* | colistin | lineage | `cart`, `cart`, `cart` | 0.939 | 0.916 | -0.023 | 0.000 | 3 of 3 | 0 | 6 of 15 |
| *Escherichia coli* | doxycycline | lineage | `cart`, `lasso`, `lasso` | 0.782 | 0.788 | 0.006 | 0.039 | 3 of 3 | 0 | 9 of 15 |
| *Escherichia coli* | ertapenem | lineage | `scm`, `scm`, `scm` | 0.864 | 0.862 | -0.002 | -0.016 | 3 of 3 | 0 | 5 of 15 |
| *Escherichia coli* | ertapenem | study | `scm` | 0.842 | 0.842 | 0.000 |  | 1 of 1 | 0 | 0 of 5 |
| *Escherichia coli* | erythromycin | lineage | `lasso`, `lasso`, `lasso` | 1.000 | 1.000 | 0.000 | -0.010 | 3 of 3 | 9 | 5 of 15 |
| *Escherichia coli* | imipenem | lineage | `scm`, `scm`, `scm` | 0.935 | 0.921 | -0.014 | 0.000 | 3 of 3 | 0 | 3 of 15 |
| *Escherichia coli* | levofloxacin | lineage | `lasso`, `gbt`, `lasso` | 0.886 | 0.903 | 0.017 | 0.002 | 3 of 3 | 0 | 11 of 15 |
| *Escherichia coli* | levofloxacin | study | `lasso` | 0.941 | 0.936 | -0.005 |  | 1 of 1 | 0 | 4 of 5 |
| *Escherichia coli* | nalidixic acid | lineage | `lasso`, `gbt`, `gbt` | 0.956 | 0.955 | -0.001 | -0.014 | 3 of 3 | 0 | 12 of 15 |
| *Escherichia coli* | nitrofurantoin | lineage | `lasso`, `scm`, `scm` | 0.767 | 0.735 | -0.032 | -0.059 | 3 of 3 | 0 | 10 of 15 |
| *Escherichia coli* | nitrofurantoin | study | `lasso` | 0.518 | 0.570 | 0.052 |  | 1 of 1 | 0 | 4 of 5 |
| *Escherichia coli* | piperacillin/tazobactam | lineage | `scm`, `scm`, `lasso` | 0.850 | 0.864 | 0.014 | 0.009 | 3 of 3 | 0 | 8 of 15 |
| *Escherichia coli* | piperacillin/tazobactam | study | `scm` | 0.876 | 0.844 | -0.032 |  | 1 of 1 | 0 | 1 of 5 |
| *Escherichia coli* | streptomycin | lineage | `scm`, `scm`, `scm` | 0.867 | 0.893 | 0.026 | 0.014 | 3 of 3 | 0 | 13 of 15 |
| *Escherichia coli* | streptomycin | study | `scm` | 0.822 | 0.844 | 0.022 |  | 1 of 1 | 0 | 4 of 5 |
| *Escherichia coli* | sulfamethoxazole | lineage | `cart`, `lasso`, `gbt` | 0.889 | 0.904 | 0.015 | 0.008 | 3 of 3 | 0 | 14 of 15 |
| *Escherichia coli* | sulfamethoxazole | study | `cart` | 0.910 | 0.921 | 0.011 |  | 1 of 1 | 0 | 5 of 5 |
| *Escherichia coli* | tetracycline | lineage | `lasso`, `lasso`, `lasso` | 0.938 | 0.922 | -0.016 | 0.006 | 3 of 3 | 0 | 10 of 15 |
| *Escherichia coli* | tetracycline | study | `lasso` | 0.938 | 0.928 | -0.010 |  | 1 of 1 | 0 | 4 of 5 |
| *Escherichia coli* | tobramycin | lineage | `scm`, `gbt`, `scm` | 0.925 | 0.917 | -0.008 | -0.002 | 3 of 3 | 0 | 9 of 15 |
| *Escherichia coli* | tobramycin | study | `scm` | 0.916 | 0.925 | 0.009 |  | 1 of 1 | 0 | 5 of 5 |
| *Escherichia coli* | trimethoprim | lineage | `scm`, `scm`, `scm` | 0.884 | 0.845 | -0.039 | 0.009 | 3 of 3 | 0 | 12 of 15 |
| *Klebsiella aerogenes* | ertapenem | lineage | `rf`, `lasso`, `gbt` | 0.575 | 0.597 | 0.022 | 0.041 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella aerogenes* | imipenem | lineage | `rf`, `rf`, `lasso` | 0.532 | 0.481 | -0.051 | -0.058 | 3 of 3 | 0 | 14 of 15 |
| *Klebsiella aerogenes* | meropenem | lineage | `lasso`, `rf`, `cart` | 0.500 | 0.583 | 0.083 | 0.000 | 3 of 3 | 9 | 4 of 15 |
| *Klebsiella oxytoca* | cefotaxime | lineage | `lasso`, `rf`, `gbt` | 0.801 | 0.821 | 0.020 | 0.024 | 3 of 3 | 0 | 12 of 15 |
| *Klebsiella oxytoca* | cefoxitin | lineage | `rf`, `rf`, `rf` | 0.852 | 0.815 | -0.037 | -0.019 | 3 of 3 | 0 | 7 of 15 |
| *Klebsiella oxytoca* | cefpodoxime | lineage | `scm`, `cart`, `gbt` | 0.713 | 0.714 | 0.001 | -0.192 | 3 of 3 | 0 | 12 of 15 |
| *Klebsiella oxytoca* | ertapenem | lineage | `rf`, `gbt`, `cart` | 0.606 | 0.666 | 0.060 | 0.073 | 3 of 3 | 0 | 9 of 15 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | lineage | `lasso`, `gbt`, `gbt` | 0.892 | 0.909 | 0.018 | -0.002 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | study | `lasso` | 0.520 | 0.533 | 0.013 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | aztreonam | lineage | `rf`, `rf`, `gbt` | 0.742 | 0.769 | 0.027 | 0.038 | 3 of 3 | 0 | 10 of 15 |
| *Klebsiella pneumoniae* | aztreonam | study | `rf` | 0.672 | 0.730 | 0.058 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | cefazolin | lineage | `lasso`, `gbt`, `scm` | 0.908 | 0.910 | 0.002 | 0.013 | 3 of 3 | 0 | 12 of 15 |
| *Klebsiella pneumoniae* | cefazolin | study | `lasso` | 0.862 | 0.907 | 0.045 |  | 1 of 1 | 0 | 5 of 5 |
| *Klebsiella pneumoniae* | cefepime | lineage | `gbt`, `rf`, `rf` | 0.720 | 0.807 | 0.087 | -0.038 | 3 of 3 | 0 | 10 of 15 |
| *Klebsiella pneumoniae* | cefepime | study | `gbt` | 0.811 | 0.800 | -0.011 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | cefotetan | lineage | `rf`, `rf`, `rf` | 0.527 | 0.508 | -0.019 | -0.048 | 3 of 3 | 9 | 6 of 15 |
| *Klebsiella pneumoniae* | cefoxitin | lineage | `scm`, `scm`, `scm` | 0.790 | 0.763 | -0.027 | -0.046 | 3 of 3 | 0 | 5 of 15 |
| *Klebsiella pneumoniae* | cefoxitin | study | `scm` | 0.770 | 0.781 | 0.011 |  | 1 of 1 | 0 | 2 of 5 |
| *Klebsiella pneumoniae* | cefpodoxime | lineage | `scm`, `gbt`, `scm` | 0.845 | 0.831 | -0.014 | -0.036 | 3 of 3 | 0 | 12 of 15 |
| *Klebsiella pneumoniae* | cefpodoxime/ceftriaxone | lineage | `cart`, `gbt`, `scm` | 0.899 | 0.845 | -0.054 | -0.012 | 3 of 3 | 0 | 14 of 15 |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | lineage | `cart`, `scm`, `scm` | 0.802 | 0.777 | -0.025 | 0.012 | 3 of 3 | 0 | 12 of 15 |
| *Klebsiella pneumoniae* | ceftazidime | lineage | `gbt`, `gbt`, `gbt` | 0.856 | 0.877 | 0.020 | 0.003 | 3 of 3 | 0 | 15 of 15 |
| *Klebsiella pneumoniae* | ceftazidime | study | `gbt` | 0.825 | 0.863 | 0.038 |  | 1 of 1 | 0 | 5 of 5 |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | lineage | `lasso`, `lasso`, `lasso` | 0.756 | 0.774 | 0.019 | 0.047 | 3 of 3 | 3 | 11 of 15 |
| *Klebsiella pneumoniae* | ceftriaxone | lineage | `lasso`, `scm`, `scm` | 0.953 | 0.958 | 0.005 | 0.003 | 3 of 3 | 0 | 6 of 15 |
| *Klebsiella pneumoniae* | ceftriaxone | study | `lasso` | 0.864 | 0.875 | 0.011 |  | 1 of 1 | 0 | 5 of 5 |
| *Klebsiella pneumoniae* | cefuroxime | lineage | `lasso`, `gbt`, `scm` | 0.891 | 0.884 | -0.007 | -0.003 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | chloramphenicol | lineage | `lasso`, `gbt`, `lasso` | 0.671 | 0.684 | 0.014 | -0.024 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | ciprofloxacin | lineage | `gbt`, `gbt`, `gbt` | 0.914 | 0.907 | -0.007 | 0.005 | 3 of 3 | 0 | 8 of 15 |
| *Klebsiella pneumoniae* | ciprofloxacin | study | `gbt` | 0.932 | 0.930 | -0.002 |  | 1 of 1 | 0 | 2 of 5 |
| *Klebsiella pneumoniae* | colistin | lineage | `scm`, `scm`, `scm` | 0.748 | 0.748 | 0.000 | 0.002 | 3 of 3 | 6 | 3 of 15 |
| *Klebsiella pneumoniae* | colistin | study | `scm` | 0.786 | 0.641 | -0.145 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | doripenem | lineage | `cart`, `lasso`, `gbt` | 0.519 | 0.508 | -0.010 | -0.029 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | doripenem | study | `cart` | 0.405 | 0.499 | 0.094 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | doxycycline | lineage | `lasso`, `rf`, `gbt` | 0.556 | 0.556 | 0.000 | -0.006 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | ertapenem | lineage | `scm`, `scm`, `scm` | 0.951 | 0.942 | -0.009 | -0.004 | 3 of 3 | 0 | 2 of 15 |
| *Klebsiella pneumoniae* | ertapenem | study | `scm` | 0.898 | 0.896 | -0.002 |  | 1 of 1 | 0 | 1 of 5 |
| *Klebsiella pneumoniae* | fosfomycin | lineage | `lasso`, `lasso`, `gbt` | 0.522 | 0.518 | -0.004 | 0.036 | 3 of 3 | 0 | 11 of 15 |
| *Klebsiella pneumoniae* | imipenem | lineage | `scm`, `scm`, `scm` | 0.945 | 0.943 | -0.002 | 0.001 | 3 of 3 | 0 | 2 of 15 |
| *Klebsiella pneumoniae* | imipenem | study | `scm` | 0.960 | 0.958 | -0.002 |  | 1 of 1 | 0 | 1 of 5 |
| *Klebsiella pneumoniae* | levofloxacin | lineage | `rf`, `gbt`, `lasso` | 0.946 | 0.960 | 0.014 | 0.015 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | levofloxacin | study | `rf` | 0.924 | 0.963 | 0.039 |  | 1 of 1 | 0 | 5 of 5 |
| *Klebsiella pneumoniae* | nitrofurantoin | lineage | `rf`, `lasso`, `gbt` | 0.576 | 0.677 | 0.102 | 0.060 | 3 of 3 | 0 | 10 of 15 |
| *Klebsiella pneumoniae* | nitrofurantoin | study | `rf` | 0.555 | 0.692 | 0.136 |  | 1 of 1 | 0 | 5 of 5 |
| *Klebsiella pneumoniae* | ofloxacin | lineage | `lasso`, `rf`, `gbt` | 0.450 | 0.584 | 0.134 | -0.001 | 3 of 3 | 9 | 5 of 15 |
| *Klebsiella pneumoniae* | polymyxin b | lineage | `scm`, `cart` | none | none | none | none | 2 of 3 | 10 | 0 of 10 |
| *Klebsiella pneumoniae* | polymyxin b | study | `scm` | 0.694 | 0.831 | 0.137 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | streptomycin | lineage | `cart`, `lasso`, `rf` | 0.466 | 0.515 | 0.049 | 0.028 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | lineage | `gbt`, `lasso`, `gbt` | 0.912 | 0.912 | 0.000 | 0.002 | 3 of 3 | 0 | 11 of 15 |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | study | `gbt` | 0.908 | 0.904 | -0.004 |  | 1 of 1 | 0 | 2 of 5 |
| *Klebsiella pneumoniae* | tetracycline | lineage | `scm`, `scm`, `scm` | 0.836 | 0.831 | -0.005 | -0.019 | 3 of 3 | 0 | 4 of 15 |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | lineage | `lasso`, `lasso`, `rf` | 0.637 | 0.652 | 0.014 | -0.048 | 3 of 3 | 9 | 5 of 15 |
| *Klebsiella pneumoniae* | tigecycline | lineage | `gbt`, `lasso`, `gbt` | 0.616 | 0.612 | -0.004 | 0.005 | 3 of 3 | 0 | 11 of 15 |
| *Klebsiella pneumoniae* | tigecycline | study | `gbt` | 0.567 | 0.593 | 0.027 |  | 1 of 1 | 0 | 4 of 5 |
| *Klebsiella pneumoniae* | tobramycin | lineage | `gbt`, `gbt`, `gbt` | 0.899 | 0.892 | -0.008 | -0.012 | 3 of 3 | 0 | 11 of 15 |
| *Klebsiella pneumoniae* | tobramycin | study | `gbt` | 0.865 | 0.877 | 0.012 |  | 1 of 1 | 0 | 2 of 5 |
| *Klebsiella pneumoniae* | trimethoprim | lineage | `rf`, `lasso`, `rf` | 0.920 | 0.927 | 0.006 | -0.011 | 3 of 3 | 0 | 13 of 15 |
| *Klebsiella quasipneumoniae* | cefotaxime | lineage | `rf`, `lasso`, `rf` | 0.946 | 0.946 | 0.000 | 0.013 | 3 of 3 | 0 | 8 of 15 |
| *Klebsiella quasipneumoniae* | ceftazidime | lineage | `rf`, `rf`, `rf` | 0.960 | 0.957 | -0.003 | 0.011 | 3 of 3 | 0 | 4 of 15 |
| *Klebsiella quasipneumoniae* | sulfamethoxazole/trimethoprim | lineage | `gbt`, `cart`, `gbt` | 0.896 | 0.890 | -0.007 | 0.009 | 3 of 3 | 0 | 11 of 15 |
| *Mycobacterium tuberculosis* | amikacin | lineage | `gbt`, `gbt`, `gbt` | 0.880 | 0.880 | 0.000 | 0.002 | 3 of 3 | 0 | 9 of 15 |
| *Mycobacterium tuberculosis* | amoxicillin | lineage | `gbt`, `scm`, `gbt` | 0.761 | 0.732 | -0.029 | -0.069 | 3 of 3 | 0 | 15 of 15 |
| *Mycobacterium tuberculosis* | capreomycin | lineage | `gbt`, `gbt`, `gbt` | 0.840 | 0.831 | -0.008 | -0.010 | 3 of 3 | 0 | 8 of 15 |
| *Mycobacterium tuberculosis* | ciprofloxacin | lineage | `lasso`, `scm`, `rf` | 0.758 | 0.722 | -0.037 | 0.043 | 3 of 3 | 0 | 13 of 15 |
| *Mycobacterium tuberculosis* | cycloserine | lineage | `lasso`, `gbt`, `lasso` | 0.606 | 0.566 | -0.040 | -0.052 | 3 of 3 | 0 | 13 of 15 |
| *Mycobacterium tuberculosis* | delamanid | lineage | `lasso`, `cart`, `cart` | 0.542 | 0.550 | 0.007 | 0.008 | 3 of 3 | 0 | 13 of 15 |
| *Mycobacterium tuberculosis* | ethambutol | lineage | `gbt`, `gbt`, `gbt` | 0.712 | 0.711 | -0.001 | -0.002 | 3 of 3 | 0 | 13 of 15 |
| *Mycobacterium tuberculosis* | ethionamide | lineage | `gbt`, `rf`, `gbt` | 0.734 | 0.731 | -0.004 | -0.008 | 3 of 3 | 0 | 8 of 15 |
| *Mycobacterium tuberculosis* | kanamycin | lineage | `lasso`, `lasso`, `lasso` | 0.894 | 0.886 | -0.007 | 0.000 | 3 of 3 | 0 | 7 of 15 |
| *Mycobacterium tuberculosis* | levofloxacin | lineage | `lasso`, `scm`, `cart` | 0.698 | 0.695 | -0.003 | 0.002 | 3 of 3 | 0 | 10 of 15 |
| *Mycobacterium tuberculosis* | linezolid | lineage | `gbt`, `cart`, `lasso` | 0.562 | 0.547 | -0.016 | 0.025 | 3 of 3 | 0 | 12 of 15 |
| *Mycobacterium tuberculosis* | moxifloxacin | lineage | `gbt`, `cart`, `gbt` | 0.821 | 0.822 | 0.001 | 0.002 | 3 of 3 | 0 | 14 of 15 |
| *Mycobacterium tuberculosis* | nicotinamide | lineage | `cart`, `lasso`, `lasso` | 0.823 | 0.813 | -0.010 | -0.018 | 3 of 3 | 0 | 11 of 15 |
| *Mycobacterium tuberculosis* | nicotinamide | study | `cart` | 0.813 | 0.801 | -0.012 |  | 1 of 1 | 0 | 4 of 5 |
| *Mycobacterium tuberculosis* | ofloxacin | lineage | `lasso`, `cart`, `cart` | 0.913 | 0.902 | -0.011 | 0.002 | 3 of 3 | 0 | 10 of 15 |
| *Mycobacterium tuberculosis* | para-aminosalicylic acid | lineage | `lasso`, `gbt`, `lasso` | 0.734 | 0.732 | -0.003 | 0.005 | 3 of 3 | 0 | 5 of 15 |
| *Mycobacterium tuberculosis* | prothionamide | lineage | `lasso`, `rf`, `scm` | 0.694 | 0.667 | -0.027 | 0.031 | 3 of 3 | 0 | 11 of 15 |
| *Mycobacterium tuberculosis* | pyrazinamide | lineage | `gbt`, `gbt`, `gbt` | 0.643 | 0.638 | -0.005 | -0.057 | 3 of 3 | 0 | 9 of 15 |
| *Mycobacterium tuberculosis* | rifabutin | lineage | `gbt`, `rf`, `gbt` | 0.885 | 0.880 | -0.005 | -0.004 | 3 of 3 | 0 | 11 of 15 |
| *Mycobacterium tuberculosis* | rifampin | lineage | `scm`, `scm`, `cart` | 0.905 | 0.904 | -0.001 | 0.000 | 3 of 3 | 0 | 5 of 15 |
| *Mycobacterium tuberculosis* | streptomycin | lineage | `cart`, `scm`, `cart` | 0.800 | 0.808 | 0.008 | -0.005 | 3 of 3 | 0 | 11 of 15 |
| *Neisseria gonorrhoeae* | cefixime | lineage | `lasso`, `rf`, `lasso` | 0.851 | 0.851 | 0.000 | -0.007 | 3 of 3 | 15 | 0 of 15 |
| *Neisseria gonorrhoeae* | ceftriaxone | lineage | `gbt`, `gbt`, `gbt` | 0.501 | 0.555 | 0.054 | -0.001 | 3 of 3 | 3 | 11 of 15 |
| *Neisseria gonorrhoeae* | ciprofloxacin | lineage | `scm`, `scm`, `scm` | 0.987 | 0.986 | -0.001 | 0.002 | 3 of 3 | 0 | 15 of 15 |
| *Neisseria gonorrhoeae* | erythromycin | lineage | `cart`, `scm`, `cart` | 0.850 | 0.824 | -0.025 | -0.003 | 3 of 3 | 0 | 13 of 15 |
| *Neisseria gonorrhoeae* | penicillin | lineage | `gbt`, `gbt`, `gbt` | 0.589 | 0.684 | 0.095 | 0.048 | 3 of 3 | 0 | 13 of 15 |
| *Neisseria gonorrhoeae* | tetracycline | lineage | `lasso`, `gbt`, `gbt` | 0.930 | 0.926 | -0.004 | -0.004 | 3 of 3 | 0 | 12 of 15 |
| *Pseudomonas aeruginosa* | amikacin | lineage | `lasso`, `gbt`, `lasso` | 0.688 | 0.765 | 0.077 | 0.014 | 3 of 3 | 0 | 15 of 15 |
| *Pseudomonas aeruginosa* | aztreonam | lineage | `rf`, `gbt`, `scm` | 0.575 | 0.586 | 0.011 | 0.017 | 3 of 3 | 0 | 14 of 15 |
| *Pseudomonas aeruginosa* | aztreonam | study | `rf` | 0.506 | 0.542 | 0.036 |  | 1 of 1 | 0 | 3 of 5 |
| *Pseudomonas aeruginosa* | cefepime | lineage | `rf`, `rf`, `rf` | 0.600 | 0.556 | -0.044 | -0.048 | 3 of 3 | 0 | 9 of 15 |
| *Pseudomonas aeruginosa* | cefepime | study | `rf` | 0.527 | 0.523 | -0.004 |  | 1 of 1 | 0 | 1 of 5 |
| *Pseudomonas aeruginosa* | ceftazidime | lineage | `gbt`, `lasso`, `lasso` | 0.705 | 0.714 | 0.008 | -0.002 | 3 of 3 | 0 | 14 of 15 |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | lineage | `scm`, `cart`, `scm` | 0.630 | 0.611 | -0.019 | 0.024 | 3 of 3 | 0 | 8 of 15 |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | study | `scm` | 0.769 | 0.748 | -0.022 |  | 1 of 1 | 0 | 4 of 5 |
| *Pseudomonas aeruginosa* | ceftolozane/tazobactam | lineage | `gbt`, `rf`, `gbt` | 0.554 | 0.651 | 0.097 | 0.035 | 3 of 3 | 0 | 12 of 15 |
| *Pseudomonas aeruginosa* | gentamicin | lineage | `gbt`, `scm`, `lasso` | 0.877 | 0.868 | -0.009 | 0.002 | 3 of 3 | 0 | 12 of 15 |
| *Pseudomonas aeruginosa* | gentamicin | study | `gbt` | 0.859 | 0.857 | -0.003 |  | 1 of 1 | 0 | 2 of 5 |
| *Pseudomonas aeruginosa* | imipenem | lineage | `rf`, `gbt`, `lasso` | 0.712 | 0.699 | -0.014 | -0.016 | 3 of 3 | 0 | 10 of 15 |
| *Pseudomonas aeruginosa* | imipenem | study | `rf` | 0.740 | 0.734 | -0.006 |  | 1 of 1 | 0 | 5 of 5 |
| *Pseudomonas aeruginosa* | levofloxacin | lineage | `cart`, `cart`, `cart` | 0.908 | 0.915 | 0.007 | -0.001 | 3 of 3 | 0 | 15 of 15 |
| *Pseudomonas aeruginosa* | meropenem | lineage | `cart`, `scm`, `scm` | 0.734 | 0.723 | -0.011 | -0.004 | 3 of 3 | 0 | 12 of 15 |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | lineage | `gbt`, `lasso`, `gbt` | 0.647 | 0.606 | -0.041 | 0.026 | 3 of 3 | 0 | 11 of 15 |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | study | `gbt` | 0.638 | 0.720 | 0.083 |  | 1 of 1 | 0 | 5 of 5 |
| *Pseudomonas aeruginosa* | tobramycin | lineage | `rf`, `rf`, `rf` | 0.896 | 0.902 | 0.006 | -0.006 | 3 of 3 | 0 | 13 of 15 |
| *Salmonella enterica* | amoxicillin/clavulanate | lineage | `scm`, `rf`, `rf` | 0.956 | 0.952 | -0.004 | -0.001 | 3 of 3 | 0 | 12 of 15 |
| *Salmonella enterica* | ampicillin | lineage | `gbt`, `scm`, `scm` | 0.879 | 0.868 | -0.011 | -0.002 | 3 of 3 | 0 | 8 of 15 |
| *Salmonella enterica* | ampicillin/sulbactam | lineage | `scm`, `gbt`, `scm` | 0.574 | 0.574 | -0.001 | 0.001 | 3 of 3 | 9 | 6 of 15 |
| *Salmonella enterica* | azithromycin | lineage | `rf`, `lasso`, `scm` | 0.895 | 0.894 | -0.001 | -0.002 | 3 of 3 | 0 | 12 of 15 |
| *Salmonella enterica* | cefazolin | lineage | `gbt`, `gbt`, `scm` | 0.484 | 0.519 | 0.034 | 0.058 | 3 of 3 | 9 | 3 of 15 |
| *Salmonella enterica* | cefotaxime | lineage | `rf`, `rf`, `rf` | 0.766 | 0.771 | 0.005 | 0.028 | 3 of 3 | 0 | 14 of 15 |
| *Salmonella enterica* | cefoxitin | lineage | `rf`, `rf`, `rf` | 0.965 | 0.962 | -0.003 | -0.006 | 3 of 3 | 0 | 7 of 15 |
| *Salmonella enterica* | ceftazidime | lineage | `scm`, `scm`, `scm` | 0.621 | 0.622 | 0.001 | -0.072 | 3 of 3 | 0 | 7 of 15 |
| *Salmonella enterica* | ceftiofur | lineage | `cart`, `rf`, `scm` | 0.936 | 0.957 | 0.021 | 0.016 | 3 of 3 | 0 | 9 of 15 |
| *Salmonella enterica* | ceftriaxone | lineage | `cart`, `scm`, `cart` | 0.923 | 0.926 | 0.003 | 0.033 | 3 of 3 | 0 | 10 of 15 |
| *Salmonella enterica* | chloramphenicol | lineage | `scm`, `scm`, `scm` | 0.905 | 0.876 | -0.029 | -0.020 | 3 of 3 | 0 | 4 of 15 |
| *Salmonella enterica* | ciprofloxacin | lineage | `cart`, `cart`, `cart` | 0.886 | 0.860 | -0.026 | -0.013 | 3 of 3 | 0 | 13 of 15 |
| *Salmonella enterica* | gentamicin | lineage | `gbt`, `gbt`, `rf` | 0.753 | 0.756 | 0.003 | -0.047 | 3 of 3 | 0 | 10 of 15 |
| *Salmonella enterica* | kanamycin | lineage | `scm`, `scm`, `scm` | 0.833 | 0.827 | -0.006 | -0.023 | 3 of 3 | 0 | 5 of 15 |
| *Salmonella enterica* | nalidixic acid | lineage | `gbt`, `gbt`, `scm` | 0.897 | 0.897 | 0.000 | 0.001 | 3 of 3 | 0 | 14 of 15 |
| *Salmonella enterica* | pefloxacin | lineage | `lasso`, `lasso`, `cart` | 0.987 | 0.987 | 0.000 | 0.000 | 3 of 3 | 0 | 9 of 15 |
| *Salmonella enterica* | spectinomycin | lineage | `scm`, `cart`, `lasso` | 0.837 | 0.706 | -0.132 | -0.006 | 3 of 3 | 0 | 14 of 15 |
| *Salmonella enterica* | streptomycin | lineage | `lasso`, `lasso`, `lasso` | 0.842 | 0.862 | 0.021 | 0.031 | 3 of 3 | 0 | 12 of 15 |
| *Salmonella enterica* | sulfamethoxazole | lineage | `scm`, `lasso`, `rf` | 0.672 | 0.718 | 0.045 | 0.097 | 3 of 3 | 12 | 2 of 15 |
| *Salmonella enterica* | sulfamethoxazole/trimethoprim | lineage | `gbt`, `gbt`, `lasso` | 0.689 | 0.720 | 0.032 | 0.014 | 3 of 3 | 0 | 12 of 15 |
| *Salmonella enterica* | sulfathiazole | lineage | `scm`, `scm`, `scm` | 0.716 | 0.662 | -0.054 | -0.053 | 3 of 3 | 0 | 12 of 15 |
| *Salmonella enterica* | sulfisoxazole | lineage | `scm`, `scm`, `scm` | 0.949 | 0.938 | -0.011 | -0.083 | 3 of 3 | 0 | 7 of 15 |
| *Salmonella enterica* | sulfisoxazole | study | `scm` | 0.955 | 0.949 | -0.006 |  | 1 of 1 | 0 | 3 of 5 |
| *Salmonella enterica* | tetracycline | lineage | `scm`, `scm`, `scm` | 0.877 | 0.877 | 0.000 | 0.022 | 3 of 3 | 0 | 1 of 15 |
| *Salmonella enterica* | trimethoprim | lineage | `cart`, `cart`, `scm` | 0.779 | 0.805 | 0.026 | 0.059 | 3 of 3 | 0 | 12 of 15 |
| *Shigella flexneri* | ampicillin | lineage | `lasso`, `scm`, `scm` | 0.901 | 0.910 | 0.009 | 0.003 | 3 of 3 | 0 | 12 of 15 |
| *Shigella flexneri* | azithromycin | lineage | `gbt`, `lasso`, `scm` | 0.935 | 0.931 | -0.004 | 0.002 | 3 of 3 | 0 | 11 of 15 |
| *Shigella flexneri* | ciprofloxacin | lineage | `cart`, `cart`, `cart` | 0.904 | 0.832 | -0.072 | 0.021 | 3 of 3 | 0 | 6 of 15 |
| *Shigella flexneri* | trimethoprim | lineage | `cart`, `scm`, `gbt` | 0.922 | 0.948 | 0.026 | 0.002 | 3 of 3 | 0 | 8 of 15 |
| *Shigella sonnei* | amikacin | lineage | `cart`, `cart`, `cart` | 0.972 | 0.955 | -0.017 | -0.005 | 3 of 3 | 0 | 9 of 15 |
| *Shigella sonnei* | ampicillin | lineage | `rf`, `rf`, `gbt` | 0.938 | 0.937 | -0.001 | 0.004 | 3 of 3 | 0 | 10 of 15 |
| *Shigella sonnei* | azithromycin | lineage | `rf`, `rf`, `rf` | 0.971 | 0.976 | 0.005 | 0.000 | 3 of 3 | 0 | 15 of 15 |
| *Shigella sonnei* | cefoperazone | lineage | `scm`, `scm`, `cart` | 0.898 | 0.791 | -0.108 | 0.052 | 3 of 3 | 6 | 8 of 15 |
| *Shigella sonnei* | cefotaxime | lineage | `scm`, `lasso`, `scm` | 0.944 | 0.951 | 0.007 | 0.003 | 3 of 3 | 0 | 9 of 15 |
| *Shigella sonnei* | cefoxitin | lineage | `gbt`, `lasso`, `lasso` | 0.905 | 0.888 | -0.017 | 0.004 | 3 of 3 | 0 | 8 of 15 |
| *Shigella sonnei* | cefpirome | lineage | `gbt`, `rf`, `gbt` | 0.955 | 0.934 | -0.021 | 0.002 | 3 of 3 | 0 | 13 of 15 |
| *Shigella sonnei* | ceftazidime | lineage | `gbt`, `gbt`, `gbt` | 0.956 | 0.953 | -0.003 | -0.007 | 3 of 3 | 0 | 9 of 15 |
| *Shigella sonnei* | ceftriaxone | lineage | `cart`, `lasso`, `lasso` | 0.991 | 0.988 | -0.003 | -0.010 | 3 of 3 | 0 | 9 of 15 |
| *Shigella sonnei* | ciprofloxacin | lineage | `gbt`, `gbt`, `gbt` | 0.965 | 0.967 | 0.002 | 0.008 | 3 of 3 | 0 | 12 of 15 |
| *Shigella sonnei* | colistin | lineage | `lasso`, `lasso`, `lasso` | 0.986 | 0.986 | 0.000 | 0.000 | 3 of 3 | 3 | 8 of 15 |
| *Shigella sonnei* | ertapenem | lineage | `gbt`, `rf`, `gbt` | 0.950 | 0.930 | -0.019 | -0.002 | 3 of 3 | 0 | 12 of 15 |
| *Shigella sonnei* | gentamicin | lineage | `gbt`, `lasso`, `cart` | 0.966 | 0.955 | -0.011 | 0.004 | 3 of 3 | 0 | 12 of 15 |
| *Shigella sonnei* | nalidixic acid | lineage | `scm`, `scm`, `scm` | 0.980 | 0.979 | 0.000 | 0.000 | 3 of 3 | 0 | 2 of 15 |
| *Shigella sonnei* | streptomycin | lineage | `gbt`, `scm`, `lasso` | 0.897 | 0.903 | 0.006 | -0.009 | 3 of 3 | 0 | 10 of 15 |
| *Shigella sonnei* | sulfamethoxazole/trimethoprim | lineage | `cart`, `scm`, `scm` | 0.893 | 0.880 | -0.013 | 0.003 | 3 of 3 | 0 | 9 of 15 |
| *Shigella sonnei* | temocillin | lineage | `gbt`, `rf`, `gbt` | 0.955 | 0.934 | -0.021 | 0.002 | 3 of 3 | 0 | 13 of 15 |
| *Shigella sonnei* | tetracycline | lineage | `scm`, `scm`, `scm` | 0.957 | 0.948 | -0.009 | 0.005 | 3 of 3 | 0 | 9 of 15 |
| *Shigella sonnei* | tobramycin | lineage | `gbt`, `gbt`, `gbt` | 0.863 | 0.836 | -0.027 | -0.060 | 3 of 3 | 0 | 5 of 15 |
| *Shigella sonnei* | trimethoprim | lineage | `gbt`, `gbt`, `gbt` | 0.900 | 0.901 | 0.001 | 0.008 | 3 of 3 | 0 | 11 of 15 |
| *Staphylococcus aureus* | cefoxitin | lineage | `cart`, `gbt`, `cart` | 0.886 | 0.879 | -0.006 | -0.010 | 3 of 3 | 0 | 10 of 15 |
| *Staphylococcus aureus* | chloramphenicol | lineage | `lasso`, `lasso`, `lasso` | 0.795 | 0.788 | -0.008 | -0.010 | 3 of 3 | 0 | 10 of 15 |
| *Staphylococcus aureus* | ciprofloxacin | lineage | `lasso`, `lasso`, `lasso` | 0.935 | 0.940 | 0.005 | 0.002 | 3 of 3 | 0 | 11 of 15 |
| *Staphylococcus aureus* | clindamycin | lineage | `gbt`, `gbt`, `gbt` | 0.488 | 0.560 | 0.072 | 0.019 | 3 of 3 | 0 | 15 of 15 |
| *Staphylococcus aureus* | clindamycin | study | `gbt` | 0.713 | 0.710 | -0.003 |  | 1 of 1 | 0 | 4 of 5 |
| *Staphylococcus aureus* | daptomycin | lineage | `lasso`, `gbt`, `lasso` | 0.577 | 0.581 | 0.004 | -0.008 | 3 of 3 | 0 | 6 of 15 |
| *Staphylococcus aureus* | erythromycin | lineage | `scm`, `scm`, `scm` | 0.877 | 0.869 | -0.008 | -0.032 | 3 of 3 | 0 | 5 of 15 |
| *Staphylococcus aureus* | erythromycin | study | `scm` | 0.896 | 0.896 | 0.000 |  | 1 of 1 | 0 | 0 of 5 |
| *Staphylococcus aureus* | fusidic acid | lineage | `scm`, `scm`, `scm` | 0.943 | 0.943 | 0.000 | -0.019 | 3 of 3 | 0 | 0 of 15 |
| *Staphylococcus aureus* | gentamicin | lineage | `gbt`, `gbt`, `rf` | 0.882 | 0.926 | 0.043 | 0.008 | 3 of 3 | 0 | 12 of 15 |
| *Staphylococcus aureus* | mupirocin | lineage | `gbt`, `scm`, `gbt` | 0.906 | 0.906 | 0.000 | 0.000 | 3 of 3 | 9 | 3 of 15 |
| *Staphylococcus aureus* | penicillin | lineage | `lasso`, `gbt`, `cart` | 0.912 | 0.905 | -0.008 | -0.003 | 3 of 3 | 0 | 11 of 15 |
| *Staphylococcus aureus* | penicillin | study | `lasso` | 0.905 | 0.901 | -0.004 |  | 1 of 1 | 0 | 5 of 5 |
| *Staphylococcus aureus* | rifampin | lineage | `lasso`, `scm`, `lasso` | 0.977 | 0.977 | 0.000 | 0.000 | 3 of 3 | 0 | 9 of 15 |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | lineage | `rf`, `rf`, `rf` | 0.886 | 0.914 | 0.028 | -0.014 | 3 of 3 | 0 | 15 of 15 |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | study | `rf` | 0.884 | 0.888 | 0.004 |  | 1 of 1 | 0 | 5 of 5 |
| *Staphylococcus aureus* | tetracycline | lineage | `lasso`, `lasso`, `lasso` | 0.742 | 0.789 | 0.047 | 0.023 | 3 of 3 | 0 | 15 of 15 |
| *Staphylococcus haemolyticus* | ciprofloxacin | lineage | `lasso`, `scm`, `gbt` | 0.918 | 0.926 | 0.008 | -0.006 | 3 of 3 | 0 | 12 of 15 |
| *Staphylococcus haemolyticus* | fusidic acid | lineage | `lasso`, `lasso`, `lasso` | 0.690 | 0.618 | -0.071 | 0.012 | 3 of 3 | 0 | 13 of 15 |
| *Staphylococcus haemolyticus* | tetracycline | lineage | `rf`, `rf`, `rf` | 0.743 | 0.707 | -0.036 | -0.001 | 3 of 3 | 0 | 7 of 15 |
| *Streptococcus agalactiae* | clindamycin | lineage | `cart`, `cart`, `cart` | 0.933 | 0.932 | -0.001 | -0.001 | 3 of 3 | 0 | 8 of 15 |
| *Streptococcus agalactiae* | erythromycin | lineage | `cart`, `cart`, `cart` | 0.923 | 0.921 | -0.002 | -0.022 | 3 of 3 | 0 | 12 of 15 |
| *Streptococcus agalactiae* | tetracycline | lineage | `lasso`, `scm`, `cart` | 0.983 | 0.980 | -0.003 | -0.003 | 3 of 3 | 0 | 9 of 15 |
| *Streptococcus pneumoniae* | azithromycin | lineage | `lasso`, `lasso`, `scm` | 0.928 | 0.888 | -0.040 | 0.015 | 3 of 3 | 0 | 11 of 15 |
| *Streptococcus pneumoniae* | cefuroxime | lineage | `lasso`, `gbt`, `lasso` | 0.987 | 0.986 | -0.002 | 0.000 | 3 of 3 | 0 | 9 of 15 |
| *Streptococcus pneumoniae* | chloramphenicol | lineage | `lasso`, `lasso`, `lasso` | 0.852 | 0.807 | -0.045 | -0.030 | 3 of 3 | 0 | 9 of 15 |
| *Streptococcus pneumoniae* | clindamycin | lineage | `cart`, `lasso`, `cart` | 0.961 | 0.945 | -0.016 | 0.002 | 3 of 3 | 0 | 11 of 15 |
| *Streptococcus pneumoniae* | erythromycin | lineage | `lasso`, `lasso`, `lasso` | 0.947 | 0.957 | 0.010 | 0.005 | 3 of 3 | 0 | 6 of 15 |
| *Streptococcus pneumoniae* | meropenem | lineage | `scm`, `scm`, `scm` | 0.917 | 0.947 | 0.030 | -0.032 | 3 of 3 | 1 | 11 of 15 |
| *Streptococcus pneumoniae* | sulfamethoxazole/trimethoprim | lineage | `lasso`, `gbt`, `lasso` | 0.967 | 0.951 | -0.016 | -0.011 | 3 of 3 | 0 | 11 of 15 |
| *Streptococcus pneumoniae* | tetracycline | lineage | `scm`, `cart`, `cart` | 0.946 | 0.946 | 0.001 | 0.002 | 3 of 3 | 0 | 12 of 15 |
| *Vibrio cholerae* | nalidixic acid | lineage | `lasso`, `lasso`, `rf` | 0.937 | 0.944 | 0.007 | 0.010 | 3 of 3 | 0 | 10 of 15 |
| *Vibrio cholerae* | sulfamethoxazole | lineage | `rf`, `rf`, `rf` | 0.986 | 0.968 | -0.018 | -0.023 | 3 of 3 | 0 | 13 of 15 |
| *Vibrio cholerae* | trimethoprim | lineage | `lasso`, `lasso`, `lasso` | 0.943 | 0.953 | 0.010 | 0.034 | 3 of 3 | 0 | 9 of 15 |
