# Screen and floor result

The release ranks the candidates of the forest (`rf`), boosted trees (`gbt`) and the lasso (`lasso`) by the chi-square statistic and keeps the top 4,000, after a candidate floor of max(3, n/200) for a dataset of n genomes. This page asks whether the Cochran-Mantel-Haenszel screen over lineage strata, at budgets of 4,000 and 16,000 (`cmh-4000`, `cmh-16000`), or the floor fixed at 3 (`floor-3`), serves better on the lineage split than the release's rule. It reads block C4 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/screen-arms.md`.

The reading's rules were fixed at 04:03Z on 23 Nov 2026 and ruled at 05:52Z, before any C4 score was read (`docs/assumption-plan-corrections.md`, entry "C4 design").

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

The Cochran-Mantel-Haenszel screen adjusts the ranking of the candidates for lineage strata and names no marker. It changes which candidates a learner reads, and this page counts only how the learners then score.

"No difference shown" is the rule's reading: an arm that does not meet the change conditions leaves the release as it is. It is no statement that two screens are equal.

## Reading

**Result.** No arm meets conditions 1 and 2, so the release's screen, budget and floor do not change. `cmh-4000` meets 3 of the 5 conditions: 3, 4 and 5. `cmh-16000` meets 3 of 5, the same three: `gbt` and `lasso` meet 1 and 2, but `rf` does not, and both conditions ask all three. `floor-3`, read by ruling 15 over the 92 datasets of 800 genomes or more, meets 2 of 5, conditions 3 and 4: no learner carries 1 and 2, so none can carry 5. The twin meets 1 and 2 on no learner. Condition 6 is the maintainer's.

**Conditions.** A screen arm is read on `rf`, `gbt` and `lasso`, the learners its screen reaches, over every dataset of the first run. It meets condition 1 when, for each of the three, the lineage median difference's interval lies wholly above 0, and condition 2 when, for each of the three, better beats worse by the one-sided sign test at p < 0.05 after Holm over the 21 lineage pairs. Condition 3 holds when the twin, over every dataset, meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. Condition 4 holds when no recall interval of the three lies wholly below -0.02, and condition 5 when each of the three has a study median of its lineage median's sign.

The floor arm is read on all five learners over the datasets of 800 genomes or more, by C3's ruling 15: a learner carries the change when it meets 1 and 2, and at least 3 of 5 must carry it. Condition 3 reads the twin over the same datasets, condition 4 holds for all five, and condition 5 for each learner that carries. The stricter reading of C3 stands beside it and decides nothing. Condition 6, for every arm, is the maintainer's.

The conditions per arm, as `reference/screen-arms/page-conditions.tsv` holds them. Learners counts the learners that meet the condition.

| Arm | Reading | Decides | Condition | Reads | Met | Learners |
|---|---|---|---|---|---|---|
| `cmh-4000` | the learners it touches, every dataset | yes | 1 | the lineage median difference's interval lies wholly above 0, for each of rf, gbt and lasso | not met | 0 |
| `cmh-4000` | the learners it touches, every dataset | yes | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 21, for each of the three | not met | 0 |
| `cmh-4000` | the learners it touches, every dataset | yes | 3 | the twin, over every dataset, meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| `cmh-4000` | the learners it touches, every dataset | yes | 4 | no recall interval of the three lies wholly below -0.02 | met | 3 |
| `cmh-4000` | the learners it touches, every dataset | yes | 5 | each of the three has a study median of its lineage median's sign | met | 3 |
| `cmh-4000` | the learners it touches, every dataset | yes | 6 | the maintainer rules the change | the maintainer's | none |
| `cmh-16000` | the learners it touches, every dataset | yes | 1 | the lineage median difference's interval lies wholly above 0, for each of rf, gbt and lasso | not met | 2 |
| `cmh-16000` | the learners it touches, every dataset | yes | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 21, for each of the three | not met | 2 |
| `cmh-16000` | the learners it touches, every dataset | yes | 3 | the twin, over every dataset, meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| `cmh-16000` | the learners it touches, every dataset | yes | 4 | no recall interval of the three lies wholly below -0.02 | met | 3 |
| `cmh-16000` | the learners it touches, every dataset | yes | 5 | each of the three has a study median of its lineage median's sign | met | 3 |
| `cmh-16000` | the learners it touches, every dataset | yes | 6 | the maintainer rules the change | the maintainer's | none |
| `floor-3` | ruling 15, 800 genomes or more | yes | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 21; at least 3 of 5 carry | not met | 0 |
| `floor-3` | ruling 15, 800 genomes or more | yes | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| `floor-3` | ruling 15, 800 genomes or more | yes | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| `floor-3` | ruling 15, 800 genomes or more | yes | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 |
| `floor-3` | ruling 15, 800 genomes or more | yes | 6 | the maintainer rules the change | the maintainer's | none |
| `floor-3` | stricter, 800 genomes or more | no, beside | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 0 |
| `floor-3` | stricter, 800 genomes or more | no, beside | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 21, for all five | not met | 0 |
| `floor-3` | stricter, 800 genomes or more | no, beside | 3 | the twin meets neither 1 nor 2 as this reading reads them | met | none |
| `floor-3` | stricter, 800 genomes or more | no, beside | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| `floor-3` | stricter, 800 genomes or more | no, beside | 5 | each of the five learners has a study median of its lineage median's sign | not met | 0 |
| `floor-3` | stricter, 800 genomes or more | no, beside | 6 | the maintainer rules the change | the maintainer's | none |

The release keeps its screen, its budget and its floor whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run: 233 datasets in 23 species, the 268 of L less the 35 in M. The datasets of M wait for the label rule. Of the 233, 55 in 8 species are in S and have a study split. 92 of the 233, in 11 species, have 800 genomes or more; only there is max(3, n/200) above 3, so only there can the floor arm change the candidates. The section Datasets names them all.

**Arms.** All four are fitted in the same run with the same binary.

| Arm | Rule | Budget | Floor |
|---|---|---|---|
| default (`k31`) | chi-square | 4,000 | max(3, n/200) |
| `cmh-4000` | Cochran-Mantel-Haenszel over lineage strata | 4,000 | max(3, n/200) |
| `cmh-16000` | Cochran-Mantel-Haenszel over lineage strata | 16,000 | max(3, n/200) |
| `floor-3` | chi-square | 4,000 | 3 |

**Draws.** The default and each arm come from runs at fold seeds 42, 43 and 44. A dataset's value is the mean of its draws, and a dataset with a draw not scored is not scored: a missing draw is never read as zero. The twin is the default at seeds 45, 46 and 47.

**Study split.** The 55 datasets in S get one study run per arm and the default, on their reviewed `study-whole` plans, at the release's seeds, with no twin. One study run and no twin is a deviation from the plan's Rules, Draws.

**Results.** `scripts/c4-screen/expect.py` lists 3,715 results: 15 lineage results for each of the 233 datasets (the default and three arms at 3 seeds, and the twin at 3) and 4 study results for each of the 55.

**Run.** All 23 species ended `DONE` and none `FAILED`. Every species ran one binary, code sha256 `5a39538947ab78961ac5c07e79adfb46851489f1782a466941df5c50aace3c63`, built from commit `2ee9e6b`, with no restart. No arm was withdrawn, so `withdrawn.txt` is empty.

**Not scored.** *Klebsiella pneumoniae* with polymyxin b, 100 genomes, leaves the lineage count, so 232 of 233 datasets enter it. Its lineage split was unavailable at seed 43, for the default and each arm, and at the twin's seed 46, and a dataset with a draw not scored is not scored. The section Datasets shows its 2 of 3 draws. Its study split ran and is counted.

## Count

The primary measure is each learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. The median difference is the arm's value less the default's, with a 95% interval from 4,000 species resamples, as C3's count takes it.

The Holm family holds 21 pairs on the lineage primary measure: `cmh-4000` and `cmh-16000` each on `rf`, `gbt` and `lasso` over every dataset; `floor-3` on all five learners over the datasets of 800 genomes or more; the twin on all five over every dataset, and again over those of 800 genomes or more. An arm withdrawn before it ran keeps its pairs in the family at p = 1.

A screen arm's `scm` and `cart` read no screen. They are checked equal to the default's in every fold and are not counted.

Recall stands beside log loss wherever log loss is given.

The 21 pairs of the Holm family, on the lineage primary measure. Better is the datasets where the arm is better by more than 0.02, of those entered.

| Pair | Datasets | Learner | Withdrawn | Better | Worse | Same | Median | Interval | Sign p | Holm p | Condition 1 | Condition 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `cmh-4000` less default | every dataset | `rf` | no | 57 of 232 | 50 | 125 | 0.001 | 0.000 to 0.004 | 0.281 | 1.000 | not met | not met |
| `cmh-4000` less default | every dataset | `gbt` | no | 62 of 232 | 45 | 125 | 0.002 | 0.000 to 0.005 | 0.061 | 1.000 | not met | not met |
| `cmh-4000` less default | every dataset | `lasso` | no | 58 of 232 | 42 | 132 | 0.002 | 0.000 to 0.007 | 0.067 | 1.000 | not met | not met |
| `cmh-16000` less default | every dataset | `rf` | no | 36 of 232 | 71 | 125 | -0.003 | -0.008 to 0.000 | 1.000 | 1.000 | not met | not met |
| `cmh-16000` less default | every dataset | `gbt` | no | 59 of 232 | 29 | 144 | 0.003 | 0.000 to 0.008 | < 0.001 | 0.018 | met | met |
| `cmh-16000` less default | every dataset | `lasso` | no | 59 of 232 | 28 | 145 | 0.004 | 0.001 to 0.007 | < 0.001 | 0.012 | met | met |
| twin less default | every dataset | `scm` | no | 29 of 232 | 27 | 176 | 0.000 | 0.000 to 0.001 | 0.447 | 1.000 | not met | not met |
| twin less default | every dataset | `cart` | no | 27 of 232 | 30 | 175 | 0.000 | 0.000 to 0.000 | 0.702 | 1.000 | not met | not met |
| twin less default | every dataset | `rf` | no | 18 of 232 | 22 | 192 | -0.001 | -0.003 to 0.000 | 0.785 | 1.000 | not met | not met |
| twin less default | every dataset | `gbt` | no | 40 of 232 | 33 | 159 | 0.000 | 0.000 to 0.001 | 0.241 | 1.000 | not met | not met |
| twin less default | every dataset | `lasso` | no | 28 of 232 | 19 | 185 | 0.000 | -0.001 to 0.000 | 0.121 | 1.000 | not met | not met |
| `floor-3` less default | 800 genomes or more | `scm` | no | 2 of 92 | 3 | 87 | 0.000 | 0.000 to 0.000 | 0.812 | 1.000 | not met | not met |
| `floor-3` less default | 800 genomes or more | `cart` | no | 0 of 92 | 1 | 91 | 0.000 | 0.000 to 0.000 | 1.000 | 1.000 | not met | not met |
| `floor-3` less default | 800 genomes or more | `rf` | no | 1 of 92 | 0 | 91 | 0.000 | 0.000 to 0.000 | 0.500 | 1.000 | not met | not met |
| `floor-3` less default | 800 genomes or more | `gbt` | no | 2 of 92 | 1 | 89 | 0.000 | 0.000 to 0.000 | 0.500 | 1.000 | not met | not met |
| `floor-3` less default | 800 genomes or more | `lasso` | no | 0 of 92 | 1 | 91 | 0.000 | 0.000 to 0.000 | 1.000 | 1.000 | not met | not met |
| twin less default | 800 genomes or more | `scm` | no | 9 of 92 | 8 | 75 | 0.000 | -0.001 to 0.002 | 0.500 | 1.000 | not met | not met |
| twin less default | 800 genomes or more | `cart` | no | 9 of 92 | 9 | 74 | 0.000 | 0.000 to 0.000 | 0.593 | 1.000 | not met | not met |
| twin less default | 800 genomes or more | `rf` | no | 3 of 92 | 4 | 85 | 0.000 | -0.002 to 0.000 | 0.773 | 1.000 | not met | not met |
| twin less default | 800 genomes or more | `gbt` | no | 10 of 92 | 9 | 73 | 0.000 | 0.000 to 0.001 | 0.500 | 1.000 | not met | not met |
| twin less default | 800 genomes or more | `lasso` | no | 6 of 92 | 4 | 82 | 0.000 | -0.002 to 0.000 | 0.377 | 1.000 | not met | not met |

The count that decides, on balanced accuracy and recall. Entered counts the datasets scored of those the set holds, and Constant model the datasets with a study run scored as the constant model.

| Split | Pair | Datasets | Learner | Measure | Better is | Entered | Better | Worse | Same | Constant model | Median | Interval |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lineage | `cmh-4000` less default | every dataset | `rf` | balanced accuracy | higher | 232 of 233 | 57 | 50 | 125 | 0 | 0.001 | 0.000 to 0.004 |
| lineage | `cmh-4000` less default | every dataset | `rf` | recall | higher | 232 of 233 | 63 | 58 | 111 | 0 | 0.000 | 0.000 to 0.005 |
| lineage | `cmh-4000` less default | every dataset | `gbt` | balanced accuracy | higher | 232 of 233 | 62 | 45 | 125 | 0 | 0.002 | 0.000 to 0.005 |
| lineage | `cmh-4000` less default | every dataset | `gbt` | recall | higher | 232 of 233 | 66 | 62 | 104 | 0 | 0.000 | 0.000 to 0.004 |
| lineage | `cmh-4000` less default | every dataset | `lasso` | balanced accuracy | higher | 232 of 233 | 58 | 42 | 132 | 0 | 0.002 | 0.000 to 0.007 |
| lineage | `cmh-4000` less default | every dataset | `lasso` | recall | higher | 232 of 233 | 71 | 55 | 106 | 0 | 0.000 | 0.000 to 0.008 |
| lineage | `cmh-16000` less default | every dataset | `rf` | balanced accuracy | higher | 232 of 233 | 36 | 71 | 125 | 0 | -0.003 | -0.008 to 0.000 |
| lineage | `cmh-16000` less default | every dataset | `rf` | recall | higher | 232 of 233 | 44 | 73 | 115 | 0 | 0.000 | -0.008 to 0.000 |
| lineage | `cmh-16000` less default | every dataset | `gbt` | balanced accuracy | higher | 232 of 233 | 59 | 29 | 144 | 0 | 0.003 | 0.000 to 0.008 |
| lineage | `cmh-16000` less default | every dataset | `gbt` | recall | higher | 232 of 233 | 70 | 38 | 124 | 0 | 0.003 | 0.000 to 0.007 |
| lineage | `cmh-16000` less default | every dataset | `lasso` | balanced accuracy | higher | 232 of 233 | 59 | 28 | 145 | 0 | 0.004 | 0.001 to 0.007 |
| lineage | `cmh-16000` less default | every dataset | `lasso` | recall | higher | 232 of 233 | 72 | 39 | 121 | 0 | 0.004 | 0.000 to 0.010 |
| lineage | twin less default | every dataset | `scm` | balanced accuracy | higher | 232 of 233 | 29 | 27 | 176 | 0 | 0.000 | 0.000 to 0.001 |
| lineage | twin less default | every dataset | `scm` | recall | higher | 232 of 233 | 44 | 43 | 145 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | twin less default | every dataset | `cart` | balanced accuracy | higher | 232 of 233 | 27 | 30 | 175 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | twin less default | every dataset | `cart` | recall | higher | 232 of 233 | 40 | 44 | 148 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | twin less default | every dataset | `rf` | balanced accuracy | higher | 232 of 233 | 18 | 22 | 192 | 0 | -0.001 | -0.003 to 0.000 |
| lineage | twin less default | every dataset | `rf` | recall | higher | 232 of 233 | 25 | 39 | 168 | 0 | 0.000 | -0.001 to 0.000 |
| lineage | twin less default | every dataset | `gbt` | balanced accuracy | higher | 232 of 233 | 40 | 33 | 159 | 0 | 0.000 | 0.000 to 0.001 |
| lineage | twin less default | every dataset | `gbt` | recall | higher | 232 of 233 | 53 | 47 | 132 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | twin less default | every dataset | `lasso` | balanced accuracy | higher | 232 of 233 | 28 | 19 | 185 | 0 | 0.000 | -0.001 to 0.000 |
| lineage | twin less default | every dataset | `lasso` | recall | higher | 232 of 233 | 41 | 32 | 159 | 0 | 0.000 | 0.000 to 0.001 |
| study | `cmh-4000` less default | every dataset | `rf` | balanced accuracy | higher | 55 of 55 | 11 | 8 | 36 | 0 | 0.001 | -0.009 to 0.007 |
| study | `cmh-4000` less default | every dataset | `rf` | recall | higher | 55 of 55 | 13 | 14 | 28 | 0 | 0.000 | -0.003 to 0.015 |
| study | `cmh-4000` less default | every dataset | `gbt` | balanced accuracy | higher | 55 of 55 | 15 | 9 | 31 | 0 | 0.005 | 0.000 to 0.009 |
| study | `cmh-4000` less default | every dataset | `gbt` | recall | higher | 55 of 55 | 18 | 9 | 28 | 0 | 0.000 | -0.002 to 0.005 |
| study | `cmh-4000` less default | every dataset | `lasso` | balanced accuracy | higher | 55 of 55 | 19 | 11 | 25 | 0 | 0.007 | 0.000 to 0.012 |
| study | `cmh-4000` less default | every dataset | `lasso` | recall | higher | 55 of 55 | 14 | 12 | 29 | 0 | 0.000 | -0.005 to 0.017 |
| study | `cmh-16000` less default | every dataset | `rf` | balanced accuracy | higher | 55 of 55 | 8 | 13 | 34 | 0 | -0.001 | -0.015 to 0.007 |
| study | `cmh-16000` less default | every dataset | `rf` | recall | higher | 55 of 55 | 11 | 21 | 23 | 0 | 0.000 | -0.027 to 0.010 |
| study | `cmh-16000` less default | every dataset | `gbt` | balanced accuracy | higher | 55 of 55 | 15 | 10 | 30 | 0 | 0.001 | -0.013 to 0.009 |
| study | `cmh-16000` less default | every dataset | `gbt` | recall | higher | 55 of 55 | 16 | 15 | 24 | 0 | 0.000 | -0.023 to 0.014 |
| study | `cmh-16000` less default | every dataset | `lasso` | balanced accuracy | higher | 55 of 55 | 18 | 7 | 30 | 0 | 0.012 | 0.005 to 0.015 |
| study | `cmh-16000` less default | every dataset | `lasso` | recall | higher | 55 of 55 | 19 | 5 | 31 | 0 | 0.010 | 0.002 to 0.016 |
| lineage | `floor-3` less default | 800 genomes or more | `scm` | balanced accuracy | higher | 92 of 92 | 2 | 3 | 87 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `scm` | recall | higher | 92 of 92 | 1 | 5 | 86 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `cart` | balanced accuracy | higher | 92 of 92 | 0 | 1 | 91 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `cart` | recall | higher | 92 of 92 | 1 | 6 | 85 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `rf` | balanced accuracy | higher | 92 of 92 | 1 | 0 | 91 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `rf` | recall | higher | 92 of 92 | 2 | 1 | 89 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `gbt` | balanced accuracy | higher | 92 of 92 | 2 | 1 | 89 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `gbt` | recall | higher | 92 of 92 | 3 | 2 | 87 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `lasso` | balanced accuracy | higher | 92 of 92 | 0 | 1 | 91 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | `floor-3` less default | 800 genomes or more | `lasso` | recall | higher | 92 of 92 | 0 | 1 | 91 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | twin less default | 800 genomes or more | `scm` | balanced accuracy | higher | 92 of 92 | 9 | 8 | 75 | 0 | 0.000 | -0.001 to 0.002 |
| lineage | twin less default | 800 genomes or more | `scm` | recall | higher | 92 of 92 | 12 | 14 | 66 | 0 | 0.000 | -0.001 to 0.003 |
| lineage | twin less default | 800 genomes or more | `cart` | balanced accuracy | higher | 92 of 92 | 9 | 9 | 74 | 0 | 0.000 | 0.000 to 0.000 |
| lineage | twin less default | 800 genomes or more | `cart` | recall | higher | 92 of 92 | 10 | 15 | 67 | 0 | 0.000 | 0.000 to 0.001 |
| lineage | twin less default | 800 genomes or more | `rf` | balanced accuracy | higher | 92 of 92 | 3 | 4 | 85 | 0 | 0.000 | -0.002 to 0.000 |
| lineage | twin less default | 800 genomes or more | `rf` | recall | higher | 92 of 92 | 8 | 11 | 73 | 0 | 0.000 | -0.001 to 0.002 |
| lineage | twin less default | 800 genomes or more | `gbt` | balanced accuracy | higher | 92 of 92 | 10 | 9 | 73 | 0 | 0.000 | 0.000 to 0.001 |
| lineage | twin less default | 800 genomes or more | `gbt` | recall | higher | 92 of 92 | 15 | 15 | 62 | 0 | 0.000 | -0.001 to 0.000 |
| lineage | twin less default | 800 genomes or more | `lasso` | balanced accuracy | higher | 92 of 92 | 6 | 4 | 82 | 0 | 0.000 | -0.002 to 0.000 |
| lineage | twin less default | 800 genomes or more | `lasso` | recall | higher | 92 of 92 | 11 | 9 | 72 | 0 | 0.000 | -0.002 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `scm` | balanced accuracy | higher | 24 of 24 | 1 | 1 | 22 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `scm` | recall | higher | 24 of 24 | 1 | 3 | 20 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `cart` | balanced accuracy | higher | 24 of 24 | 1 | 0 | 23 | 0 | 0.000 | 0.000 to 0.006 |
| study | `floor-3` less default | 800 genomes or more | `cart` | recall | higher | 24 of 24 | 1 | 0 | 23 | 0 | 0.000 | 0.000 to 0.007 |
| study | `floor-3` less default | 800 genomes or more | `rf` | balanced accuracy | higher | 24 of 24 | 0 | 0 | 24 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `rf` | recall | higher | 24 of 24 | 0 | 0 | 24 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `gbt` | balanced accuracy | higher | 24 of 24 | 0 | 0 | 24 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `gbt` | recall | higher | 24 of 24 | 0 | 1 | 23 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `lasso` | balanced accuracy | higher | 24 of 24 | 0 | 1 | 23 | 0 | 0.000 | 0.000 to 0.000 |
| study | `floor-3` less default | 800 genomes or more | `lasso` | recall | higher | 24 of 24 | 0 | 1 | 23 | 0 | 0.000 | 0.000 to 0.000 |

Every measure, as the median difference with better / worse in brackets. For balanced accuracy, recall and specificity higher is better. For log loss and both error rates lower is better, so better counts a lower value.

| Split | Pair | Datasets | Learner | Balanced accuracy | Recall | Specificity | Very major error rate | Major error rate | Log loss |
|---|---|---|---|---|---|---|---|---|---|
| lineage | `cmh-4000` less default | every dataset | `rf` | 0.001 (57 / 50) | 0.000 (63 / 58) | 0.000 (46 / 46) | 0.000 (63 / 58) | 0.000 (46 / 46) | -0.009 (69 / 43) |
| lineage | `cmh-4000` less default | every dataset | `gbt` | 0.002 (62 / 45) | 0.000 (66 / 62) | 0.000 (58 / 49) | 0.000 (66 / 62) | 0.000 (58 / 49) | -0.004 (75 / 54) |
| lineage | `cmh-4000` less default | every dataset | `lasso` | 0.002 (58 / 42) | 0.000 (71 / 55) | 0.000 (47 / 55) | 0.000 (71 / 55) | 0.000 (47 / 55) | -0.004 (86 / 56) |
| lineage | `cmh-16000` less default | every dataset | `rf` | -0.003 (36 / 71) | 0.000 (44 / 73) | 0.000 (48 / 48) | 0.000 (44 / 73) | 0.000 (48 / 48) | 0.003 (49 / 68) |
| lineage | `cmh-16000` less default | every dataset | `gbt` | 0.003 (59 / 29) | 0.003 (70 / 38) | 0.004 (59 / 42) | -0.003 (70 / 38) | -0.004 (59 / 42) | -0.012 (84 / 35) |
| lineage | `cmh-16000` less default | every dataset | `lasso` | 0.004 (59 / 28) | 0.004 (72 / 39) | 0.000 (59 / 42) | -0.004 (72 / 39) | 0.000 (59 / 42) | -0.006 (84 / 45) |
| lineage | twin less default | every dataset | `scm` | 0.000 (29 / 27) | 0.000 (44 / 43) | 0.000 (45 / 36) | 0.000 (44 / 43) | 0.000 (45 / 36) | 0.000 (85 / 79) |
| lineage | twin less default | every dataset | `cart` | 0.000 (27 / 30) | 0.000 (40 / 44) | 0.000 (37 / 38) | 0.000 (40 / 44) | 0.000 (37 / 38) | 0.000 (81 / 56) |
| lineage | twin less default | every dataset | `rf` | -0.001 (18 / 22) | 0.000 (25 / 39) | 0.000 (21 / 28) | 0.000 (25 / 39) | 0.000 (21 / 28) | 0.001 (10 / 18) |
| lineage | twin less default | every dataset | `gbt` | 0.000 (40 / 33) | 0.000 (53 / 47) | 0.000 (37 / 38) | 0.000 (53 / 47) | 0.000 (37 / 38) | 0.001 (36 / 48) |
| lineage | twin less default | every dataset | `lasso` | 0.000 (28 / 19) | 0.000 (41 / 32) | 0.000 (33 / 31) | 0.000 (41 / 32) | 0.000 (33 / 31) | 0.000 (52 / 54) |
| study | `cmh-4000` less default | every dataset | `rf` | 0.001 (11 / 8) | 0.000 (13 / 14) | 0.006 (14 / 12) | 0.000 (13 / 14) | -0.006 (14 / 12) | -0.008 (19 / 11) |
| study | `cmh-4000` less default | every dataset | `gbt` | 0.005 (15 / 9) | 0.000 (18 / 9) | 0.007 (17 / 10) | 0.000 (18 / 9) | -0.007 (17 / 10) | -0.009 (22 / 13) |
| study | `cmh-4000` less default | every dataset | `lasso` | 0.007 (19 / 11) | 0.000 (14 / 12) | 0.008 (22 / 12) | 0.000 (14 / 12) | -0.008 (22 / 12) | -0.005 (22 / 15) |
| study | `cmh-16000` less default | every dataset | `rf` | -0.001 (8 / 13) | 0.000 (11 / 21) | 0.000 (12 / 8) | 0.000 (11 / 21) | 0.000 (12 / 8) | -0.010 (19 / 11) |
| study | `cmh-16000` less default | every dataset | `gbt` | 0.001 (15 / 10) | 0.000 (16 / 15) | 0.006 (18 / 10) | 0.000 (16 / 15) | -0.006 (18 / 10) | -0.012 (24 / 11) |
| study | `cmh-16000` less default | every dataset | `lasso` | 0.012 (18 / 7) | 0.010 (19 / 5) | 0.011 (25 / 10) | -0.010 (19 / 5) | -0.011 (25 / 10) | -0.014 (25 / 13) |
| lineage | `floor-3` less default | 800 genomes or more | `scm` | 0.000 (2 / 3) | 0.000 (1 / 5) | 0.000 (2 / 1) | 0.000 (1 / 5) | 0.000 (2 / 1) | 0.000 (7 / 11) |
| lineage | `floor-3` less default | 800 genomes or more | `cart` | 0.000 (0 / 1) | 0.000 (1 / 6) | 0.000 (3 / 2) | 0.000 (1 / 6) | 0.000 (3 / 2) | 0.000 (10 / 14) |
| lineage | `floor-3` less default | 800 genomes or more | `rf` | 0.000 (1 / 0) | 0.000 (2 / 1) | 0.000 (0 / 0) | 0.000 (2 / 1) | 0.000 (0 / 0) | 0.000 (0 / 0) |
| lineage | `floor-3` less default | 800 genomes or more | `gbt` | 0.000 (2 / 1) | 0.000 (3 / 2) | 0.000 (1 / 1) | 0.000 (3 / 2) | 0.000 (1 / 1) | 0.000 (0 / 1) |
| lineage | `floor-3` less default | 800 genomes or more | `lasso` | 0.000 (0 / 1) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (0 / 0) |
| lineage | twin less default | 800 genomes or more | `scm` | 0.000 (9 / 8) | 0.000 (12 / 14) | 0.000 (9 / 9) | 0.000 (12 / 14) | 0.000 (9 / 9) | 0.000 (30 / 27) |
| lineage | twin less default | 800 genomes or more | `cart` | 0.000 (9 / 9) | 0.000 (10 / 15) | 0.000 (17 / 6) | 0.000 (10 / 15) | 0.000 (17 / 6) | -0.002 (28 / 13) |
| lineage | twin less default | 800 genomes or more | `rf` | 0.000 (3 / 4) | 0.000 (8 / 11) | 0.000 (4 / 7) | 0.000 (8 / 11) | 0.000 (4 / 7) | 0.001 (0 / 3) |
| lineage | twin less default | 800 genomes or more | `gbt` | 0.000 (10 / 9) | 0.000 (15 / 15) | 0.001 (12 / 9) | 0.000 (15 / 15) | -0.001 (12 / 9) | 0.000 (10 / 13) |
| lineage | twin less default | 800 genomes or more | `lasso` | 0.000 (6 / 4) | 0.000 (11 / 9) | -0.001 (2 / 9) | 0.000 (11 / 9) | 0.001 (2 / 9) | -0.001 (15 / 16) |
| study | `floor-3` less default | 800 genomes or more | `scm` | 0.000 (1 / 1) | 0.000 (1 / 3) | 0.000 (0 / 0) | 0.000 (1 / 3) | 0.000 (0 / 0) | 0.000 (2 / 2) |
| study | `floor-3` less default | 800 genomes or more | `cart` | 0.000 (1 / 0) | 0.000 (1 / 0) | 0.000 (0 / 1) | 0.000 (1 / 0) | 0.000 (0 / 1) | 0.000 (3 / 1) |
| study | `floor-3` less default | 800 genomes or more | `rf` | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) |
| study | `floor-3` less default | 800 genomes or more | `gbt` | 0.000 (0 / 0) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (1 / 0) |
| study | `floor-3` less default | 800 genomes or more | `lasso` | 0.000 (0 / 1) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (0 / 0) |

**Recall and log loss.** On the lineage split no arm's recall median lies below 0. The interval reaching lowest is `cmh-16000`'s on `rf`, -0.008 to 0.000, with a median of 0.000. Log loss is lower for `cmh-4000` on all three learners (medians -0.009, -0.004 and -0.004) and for `cmh-16000` on `gbt` and `lasso` (-0.012 and -0.006), and higher for `cmh-16000` on `rf` (0.003).

## Per arm

Each arm is read on its own conditions, against the same default and the same twin.

**`cmh-4000`.** The Cochran-Mantel-Haenszel screen at the release's budget of 4,000. It meets 3 of 5 conditions.

- Condition 1, not met. On 232 of 233 datasets `rf` is better on 57, worse on 50 and the same on 125, with a median of 0.001 and an interval of 0.000 to 0.004. `gbt` is better on 62 and worse on 45, median 0.002 (0.000 to 0.005). `lasso` is better on 58 and worse on 42, median 0.002 (0.000 to 0.007). Each interval starts at 0, so none lies wholly above it.
- Condition 2, not met. The sign tests give p 0.281, 0.061 and 0.067, and 1.000 each after Holm.
- Condition 3, met. The twin meets 1 and 2 on none of its five learners.
- Condition 4, met. The recall intervals are 0.000 to 0.005, 0.000 to 0.004 and 0.000 to 0.008.
- Condition 5, met. The study medians over 55 datasets, 0.001, 0.005 and 0.007, have the lineage medians' sign.

**`cmh-16000`.** The Cochran-Mantel-Haenszel screen at a budget of 16,000. It meets 3 of 5 conditions.

- Condition 1, not met. On 232 of 233 datasets `gbt` is better on 59 and worse on 29, median 0.003, with an interval from 0.00047 to 0.008 (0.000 to 0.008 at 3 decimals). `lasso` is better on 59 and worse on 28, median 0.004 (0.001 to 0.007). Both intervals lie wholly above 0. `rf` is better on 36 and worse on 71, median -0.003 (-0.008 to 0.000). The condition asks all three.
- Condition 2, not met. `gbt` and `lasso` pass the sign test at 0.018 and 0.012 after Holm. `rf` gives p 1.000, and 1.000 after Holm.
- Condition 3, met. The twin meets 1 and 2 on none of its five learners.
- Condition 4, met. The recall intervals are -0.008 to 0.000 for `rf`, 0.000 to 0.007 for `gbt` and 0.000 to 0.010 for `lasso`.
- Condition 5, met. The study medians over 55 datasets, -0.001, 0.001 and 0.012, have the lineage medians' sign.

**`floor-3`.** The release's chi-square screen at 4,000, with the floor fixed at 3. It is read over the datasets of 800 genomes or more. It meets 2 of 5 conditions.

- Conditions 1 and 2, not met. On 92 of 92 datasets `scm` is better on 2 and worse on 3, `cart` better on 0 and worse on 1, `rf` 1 and 0, `gbt` 2 and 1, `lasso` 0 and 1. Every median is 0.000 with an interval of 0.000 to 0.000. The sign tests give p 0.812, 1.000, 0.500, 0.500 and 1.000, and 1.000 each after Holm. No learner carries, against the 3 of 5 that ruling 15 asks.
- Condition 3, met. The twin over the same 92 datasets meets 1 and 2 on none of its five learners.
- Condition 4, met. Every recall interval is 0.000 to 0.000.
- Condition 5, not met. No learner carries the change.

## Beside it

These counts stand beside the count that decides, and decide nothing.

**Floor arm over every dataset.** `count-floor-all.tsv` counts the floor arm over the whole first run, the 141 datasets of under 800 genomes entering at a difference of 0. With 141 of 233 differences fixed at 0 its median is 0 by construction. Under 800 genomes the floor arm's 860 figures all equal the default's, so over 232 of 233 datasets it is better and worse on the same datasets as over the 92: `scm` 2 and 3, `cart` 0 and 1, `rf` 1 and 0, `gbt` 2 and 1, `lasso` 0 and 1. Every median is 0.000.

| Beside | Split | Pair | Datasets | Learner | Measure | Entered | Better | Worse | Same | Median | Interval | Sign p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| floor-all | lineage | `floor-3` less default | every dataset | `scm` | balanced accuracy | 232 of 233 | 2 | 3 | 227 | 0.000 | 0.000 to 0.000 | 0.812 |
| floor-all | lineage | `floor-3` less default | every dataset | `scm` | recall | 232 of 233 | 1 | 5 | 226 | 0.000 | 0.000 to 0.000 | 0.984 |
| floor-all | lineage | `floor-3` less default | every dataset | `cart` | balanced accuracy | 232 of 233 | 0 | 1 | 231 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | lineage | `floor-3` less default | every dataset | `cart` | recall | 232 of 233 | 1 | 6 | 225 | 0.000 | 0.000 to 0.000 | 0.992 |
| floor-all | lineage | `floor-3` less default | every dataset | `rf` | balanced accuracy | 232 of 233 | 1 | 0 | 231 | 0.000 | 0.000 to 0.000 | 0.500 |
| floor-all | lineage | `floor-3` less default | every dataset | `rf` | recall | 232 of 233 | 2 | 1 | 229 | 0.000 | 0.000 to 0.000 | 0.500 |
| floor-all | lineage | `floor-3` less default | every dataset | `gbt` | balanced accuracy | 232 of 233 | 2 | 1 | 229 | 0.000 | 0.000 to 0.000 | 0.500 |
| floor-all | lineage | `floor-3` less default | every dataset | `gbt` | recall | 232 of 233 | 3 | 2 | 227 | 0.000 | 0.000 to 0.000 | 0.500 |
| floor-all | lineage | `floor-3` less default | every dataset | `lasso` | balanced accuracy | 232 of 233 | 0 | 1 | 231 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | lineage | `floor-3` less default | every dataset | `lasso` | recall | 232 of 233 | 0 | 1 | 231 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | study | `floor-3` less default | every dataset | `scm` | balanced accuracy | 55 of 55 | 1 | 1 | 53 | 0.000 | 0.000 to 0.000 | 0.750 |
| floor-all | study | `floor-3` less default | every dataset | `scm` | recall | 55 of 55 | 1 | 3 | 51 | 0.000 | 0.000 to 0.000 | 0.938 |
| floor-all | study | `floor-3` less default | every dataset | `cart` | balanced accuracy | 55 of 55 | 1 | 0 | 54 | 0.000 | 0.000 to 0.000 | 0.500 |
| floor-all | study | `floor-3` less default | every dataset | `cart` | recall | 55 of 55 | 1 | 0 | 54 | 0.000 | 0.000 to 0.000 | 0.500 |
| floor-all | study | `floor-3` less default | every dataset | `rf` | balanced accuracy | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | study | `floor-3` less default | every dataset | `rf` | recall | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | study | `floor-3` less default | every dataset | `gbt` | balanced accuracy | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | study | `floor-3` less default | every dataset | `gbt` | recall | 55 of 55 | 0 | 1 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | study | `floor-3` less default | every dataset | `lasso` | balanced accuracy | 55 of 55 | 0 | 1 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| floor-all | study | `floor-3` less default | every dataset | `lasso` | recall | 55 of 55 | 0 | 1 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |

**No constant model.** `count-no-constant.tsv` leaves out every dataset with a study run refused for no candidate and scored as the constant model (ruling 16). Ruling 16 cannot arise in C4: every arm screens the default's candidates or more. No study run was refused, so `count-no-constant.tsv` equals `count.tsv` row for row.

**Stricter reading.** C3's stricter reading of the floor arm's conditions is written beside ruling 15. It meets 3 and 4 and not 1, 2 or 5: no interval lies wholly above 0, no sign test passes Holm, and every lineage median is 0, which meets no condition 5.

**Species.** The counts per species are dataset and learner pairs over each species' own datasets, better, worse and the same on the primary measure.

| Species | Split | Pair | Datasets | Learners | Better / worse / same, dataset and learner pairs | Not scored |
|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | lineage | `cmh-4000` less default | every dataset | 3 | 11 / 4 / 6 of 21 | 0 |
| *Acinetobacter baumannii* | lineage | `cmh-16000` less default | every dataset | 3 | 7 / 2 / 12 of 21 | 0 |
| *Acinetobacter baumannii* | lineage | `floor-3` less default | 800 genomes or more | 5 | 1 / 0 / 9 of 10 | 0 |
| *Acinetobacter baumannii* | lineage | twin less default | every dataset | 5 | 12 / 6 / 17 of 35 | 0 |
| *Acinetobacter baumannii* | lineage | twin less default | 800 genomes or more | 5 | 3 / 1 / 6 of 10 | 0 |
| *Acinetobacter baumannii* | study | `cmh-4000` less default | every dataset | 3 | 2 / 2 / 2 of 6 | 0 |
| *Acinetobacter baumannii* | study | `cmh-16000` less default | every dataset | 3 | 1 / 4 / 1 of 6 | 0 |
| *Acinetobacter baumannii* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter coli* | lineage | `cmh-4000` less default | every dataset | 3 | 3 / 4 / 17 of 24 | 0 |
| *Campylobacter coli* | lineage | `cmh-16000` less default | every dataset | 3 | 2 / 7 / 15 of 24 | 0 |
| *Campylobacter coli* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter coli* | lineage | twin less default | every dataset | 5 | 2 / 4 / 34 of 40 | 0 |
| *Campylobacter coli* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter coli* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter coli* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter coli* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter jejuni* | lineage | `cmh-4000` less default | every dataset | 3 | 0 / 4 / 8 of 12 | 0 |
| *Campylobacter jejuni* | lineage | `cmh-16000` less default | every dataset | 3 | 1 / 2 / 9 of 12 | 0 |
| *Campylobacter jejuni* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter jejuni* | lineage | twin less default | every dataset | 5 | 3 / 1 / 16 of 20 | 0 |
| *Campylobacter jejuni* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter jejuni* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter jejuni* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Campylobacter jejuni* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Clostridioides difficile* | lineage | `cmh-4000` less default | every dataset | 3 | 0 / 10 / 2 of 12 | 0 |
| *Clostridioides difficile* | lineage | `cmh-16000` less default | every dataset | 3 | 0 / 6 / 6 of 12 | 0 |
| *Clostridioides difficile* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Clostridioides difficile* | lineage | twin less default | every dataset | 5 | 0 / 2 / 18 of 20 | 0 |
| *Clostridioides difficile* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Clostridioides difficile* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Clostridioides difficile* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Clostridioides difficile* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Corynebacterium diphtheriae* | lineage | `cmh-4000` less default | every dataset | 3 | 3 / 2 / 4 of 9 | 0 |
| *Corynebacterium diphtheriae* | lineage | `cmh-16000` less default | every dataset | 3 | 1 / 3 / 5 of 9 | 0 |
| *Corynebacterium diphtheriae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Corynebacterium diphtheriae* | lineage | twin less default | every dataset | 5 | 3 / 3 / 9 of 15 | 0 |
| *Corynebacterium diphtheriae* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Corynebacterium diphtheriae* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Corynebacterium diphtheriae* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Corynebacterium diphtheriae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Enterobacter cloacae* | lineage | `cmh-4000` less default | every dataset | 3 | 10 / 5 / 12 of 27 | 0 |
| *Enterobacter cloacae* | lineage | `cmh-16000` less default | every dataset | 3 | 6 / 10 / 11 of 27 | 0 |
| *Enterobacter cloacae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Enterobacter cloacae* | lineage | twin less default | every dataset | 5 | 11 / 5 / 29 of 45 | 0 |
| *Enterobacter cloacae* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Enterobacter cloacae* | study | `cmh-4000` less default | every dataset | 3 | 5 / 7 / 12 of 24 | 0 |
| *Enterobacter cloacae* | study | `cmh-16000` less default | every dataset | 3 | 3 / 8 / 13 of 24 | 0 |
| *Enterobacter cloacae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Enterococcus faecium* | lineage | `cmh-4000` less default | every dataset | 3 | 6 / 5 / 22 of 33 | 0 |
| *Enterococcus faecium* | lineage | `cmh-16000` less default | every dataset | 3 | 6 / 2 / 25 of 33 | 0 |
| *Enterococcus faecium* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 1 / 49 of 50 | 0 |
| *Enterococcus faecium* | lineage | twin less default | every dataset | 5 | 11 / 2 / 42 of 55 | 0 |
| *Enterococcus faecium* | lineage | twin less default | 800 genomes or more | 5 | 10 / 2 / 38 of 50 | 0 |
| *Enterococcus faecium* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Enterococcus faecium* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Enterococcus faecium* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Escherichia coli* | lineage | `cmh-4000` less default | every dataset | 3 | 22 / 10 / 55 of 87 | 0 |
| *Escherichia coli* | lineage | `cmh-16000` less default | every dataset | 3 | 17 / 21 / 49 of 87 | 0 |
| *Escherichia coli* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 2 / 48 of 50 | 0 |
| *Escherichia coli* | lineage | twin less default | every dataset | 5 | 9 / 10 / 126 of 145 | 0 |
| *Escherichia coli* | lineage | twin less default | 800 genomes or more | 5 | 2 / 1 / 47 of 50 | 0 |
| *Escherichia coli* | study | `cmh-4000` less default | every dataset | 3 | 9 / 2 / 34 of 45 | 0 |
| *Escherichia coli* | study | `cmh-16000` less default | every dataset | 3 | 8 / 5 / 32 of 45 | 0 |
| *Escherichia coli* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 2 / 33 of 35 | 0 |
| *Klebsiella aerogenes* | lineage | `cmh-4000` less default | every dataset | 3 | 2 / 2 / 5 of 9 | 0 |
| *Klebsiella aerogenes* | lineage | `cmh-16000` less default | every dataset | 3 | 1 / 2 / 6 of 9 | 0 |
| *Klebsiella aerogenes* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella aerogenes* | lineage | twin less default | every dataset | 5 | 4 / 5 / 6 of 15 | 0 |
| *Klebsiella aerogenes* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella aerogenes* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella aerogenes* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella aerogenes* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella oxytoca* | lineage | `cmh-4000` less default | every dataset | 3 | 5 / 4 / 3 of 12 | 0 |
| *Klebsiella oxytoca* | lineage | `cmh-16000` less default | every dataset | 3 | 4 / 4 / 4 of 12 | 0 |
| *Klebsiella oxytoca* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella oxytoca* | lineage | twin less default | every dataset | 5 | 8 / 7 / 5 of 20 | 0 |
| *Klebsiella oxytoca* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella oxytoca* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella oxytoca* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella oxytoca* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella pneumoniae* | lineage | `cmh-4000` less default | every dataset | 3 | 38 / 20 / 35 of 96 | 3 |
| *Klebsiella pneumoniae* | lineage | `cmh-16000` less default | every dataset | 3 | 36 / 12 / 45 of 96 | 3 |
| *Klebsiella pneumoniae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 65 of 65 | 0 |
| *Klebsiella pneumoniae* | lineage | twin less default | every dataset | 5 | 18 / 26 / 111 of 160 | 5 |
| *Klebsiella pneumoniae* | lineage | twin less default | 800 genomes or more | 5 | 3 / 7 / 55 of 65 | 0 |
| *Klebsiella pneumoniae* | study | `cmh-4000` less default | every dataset | 3 | 19 / 9 / 26 of 54 | 0 |
| *Klebsiella pneumoniae* | study | `cmh-16000` less default | every dataset | 3 | 17 / 7 / 30 of 54 | 0 |
| *Klebsiella pneumoniae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 65 of 65 | 0 |
| *Klebsiella quasipneumoniae* | lineage | `cmh-4000` less default | every dataset | 3 | 1 / 1 / 7 of 9 | 0 |
| *Klebsiella quasipneumoniae* | lineage | `cmh-16000` less default | every dataset | 3 | 1 / 1 / 7 of 9 | 0 |
| *Klebsiella quasipneumoniae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella quasipneumoniae* | lineage | twin less default | every dataset | 5 | 6 / 0 / 9 of 15 | 0 |
| *Klebsiella quasipneumoniae* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella quasipneumoniae* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella quasipneumoniae* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Klebsiella quasipneumoniae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Mycobacterium tuberculosis* | lineage | `cmh-4000` less default | every dataset | 3 | 5 / 9 / 46 of 60 | 0 |
| *Mycobacterium tuberculosis* | lineage | `cmh-16000` less default | every dataset | 3 | 2 / 7 / 51 of 60 | 0 |
| *Mycobacterium tuberculosis* | lineage | `floor-3` less default | 800 genomes or more | 5 | 1 / 0 / 54 of 55 | 0 |
| *Mycobacterium tuberculosis* | lineage | twin less default | every dataset | 5 | 10 / 11 / 79 of 100 | 0 |
| *Mycobacterium tuberculosis* | lineage | twin less default | 800 genomes or more | 5 | 1 / 1 / 53 of 55 | 0 |
| *Mycobacterium tuberculosis* | study | `cmh-4000` less default | every dataset | 3 | 0 / 1 / 2 of 3 | 0 |
| *Mycobacterium tuberculosis* | study | `cmh-16000` less default | every dataset | 3 | 0 / 1 / 2 of 3 | 0 |
| *Mycobacterium tuberculosis* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Neisseria gonorrhoeae* | lineage | `cmh-4000` less default | every dataset | 3 | 5 / 7 / 6 of 18 | 0 |
| *Neisseria gonorrhoeae* | lineage | `cmh-16000` less default | every dataset | 3 | 3 / 6 / 9 of 18 | 0 |
| *Neisseria gonorrhoeae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 1 / 1 / 18 of 20 | 0 |
| *Neisseria gonorrhoeae* | lineage | twin less default | every dataset | 5 | 4 / 4 / 22 of 30 | 0 |
| *Neisseria gonorrhoeae* | lineage | twin less default | 800 genomes or more | 5 | 2 / 3 / 15 of 20 | 0 |
| *Neisseria gonorrhoeae* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Neisseria gonorrhoeae* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Neisseria gonorrhoeae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Pseudomonas aeruginosa* | lineage | `cmh-4000` less default | every dataset | 3 | 7 / 9 / 20 of 36 | 0 |
| *Pseudomonas aeruginosa* | lineage | `cmh-16000` less default | every dataset | 3 | 4 / 10 / 22 of 36 | 0 |
| *Pseudomonas aeruginosa* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 10 of 10 | 0 |
| *Pseudomonas aeruginosa* | lineage | twin less default | every dataset | 5 | 8 / 8 / 44 of 60 | 0 |
| *Pseudomonas aeruginosa* | lineage | twin less default | 800 genomes or more | 5 | 0 / 1 / 9 of 10 | 0 |
| *Pseudomonas aeruginosa* | study | `cmh-4000` less default | every dataset | 3 | 5 / 6 / 7 of 18 | 0 |
| *Pseudomonas aeruginosa* | study | `cmh-16000` less default | every dataset | 3 | 6 / 4 / 8 of 18 | 0 |
| *Pseudomonas aeruginosa* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Salmonella enterica* | lineage | `cmh-4000` less default | every dataset | 3 | 30 / 8 / 34 of 72 | 0 |
| *Salmonella enterica* | lineage | `cmh-16000` less default | every dataset | 3 | 33 / 7 / 32 of 72 | 0 |
| *Salmonella enterica* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 60 of 60 | 0 |
| *Salmonella enterica* | lineage | twin less default | every dataset | 5 | 14 / 18 / 88 of 120 | 0 |
| *Salmonella enterica* | lineage | twin less default | 800 genomes or more | 5 | 8 / 9 / 43 of 60 | 0 |
| *Salmonella enterica* | study | `cmh-4000` less default | every dataset | 3 | 3 / 0 / 0 of 3 | 0 |
| *Salmonella enterica* | study | `cmh-16000` less default | every dataset | 3 | 3 / 0 / 0 of 3 | 0 |
| *Salmonella enterica* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Shigella flexneri* | lineage | `cmh-4000` less default | every dataset | 3 | 2 / 3 / 7 of 12 | 0 |
| *Shigella flexneri* | lineage | `cmh-16000` less default | every dataset | 3 | 3 / 1 / 8 of 12 | 0 |
| *Shigella flexneri* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Shigella flexneri* | lineage | twin less default | every dataset | 5 | 2 / 0 / 18 of 20 | 0 |
| *Shigella flexneri* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Shigella flexneri* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Shigella flexneri* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Shigella flexneri* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Shigella sonnei* | lineage | `cmh-4000` less default | every dataset | 3 | 5 / 15 / 40 of 60 | 0 |
| *Shigella sonnei* | lineage | `cmh-16000` less default | every dataset | 3 | 8 / 12 / 40 of 60 | 0 |
| *Shigella sonnei* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 40 of 40 | 0 |
| *Shigella sonnei* | lineage | twin less default | every dataset | 5 | 3 / 4 / 93 of 100 | 0 |
| *Shigella sonnei* | lineage | twin less default | 800 genomes or more | 5 | 1 / 0 / 39 of 40 | 0 |
| *Shigella sonnei* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Shigella sonnei* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Shigella sonnei* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Staphylococcus aureus* | lineage | `cmh-4000` less default | every dataset | 3 | 14 / 5 / 20 of 39 | 0 |
| *Staphylococcus aureus* | lineage | `cmh-16000` less default | every dataset | 3 | 9 / 6 / 24 of 39 | 0 |
| *Staphylococcus aureus* | lineage | `floor-3` less default | 800 genomes or more | 5 | 2 / 1 / 57 of 60 | 0 |
| *Staphylococcus aureus* | lineage | twin less default | every dataset | 5 | 3 / 6 / 56 of 65 | 0 |
| *Staphylococcus aureus* | lineage | twin less default | 800 genomes or more | 5 | 2 / 6 / 52 of 60 | 0 |
| *Staphylococcus aureus* | study | `cmh-4000` less default | every dataset | 3 | 2 / 1 / 9 of 12 | 0 |
| *Staphylococcus aureus* | study | `cmh-16000` less default | every dataset | 3 | 3 / 1 / 8 of 12 | 0 |
| *Staphylococcus aureus* | study | `floor-3` less default | 800 genomes or more | 5 | 2 / 0 / 18 of 20 | 0 |
| *Staphylococcus haemolyticus* | lineage | `cmh-4000` less default | every dataset | 3 | 1 / 4 / 4 of 9 | 0 |
| *Staphylococcus haemolyticus* | lineage | `cmh-16000` less default | every dataset | 3 | 2 / 2 / 5 of 9 | 0 |
| *Staphylococcus haemolyticus* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Staphylococcus haemolyticus* | lineage | twin less default | every dataset | 5 | 4 / 3 / 8 of 15 | 0 |
| *Staphylococcus haemolyticus* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Staphylococcus haemolyticus* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Staphylococcus haemolyticus* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Staphylococcus haemolyticus* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus agalactiae* | lineage | `cmh-4000` less default | every dataset | 3 | 2 / 0 / 7 of 9 | 0 |
| *Streptococcus agalactiae* | lineage | `cmh-16000` less default | every dataset | 3 | 1 / 0 / 8 of 9 | 0 |
| *Streptococcus agalactiae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus agalactiae* | lineage | twin less default | every dataset | 5 | 0 / 1 / 14 of 15 | 0 |
| *Streptococcus agalactiae* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus agalactiae* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus agalactiae* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus agalactiae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus pneumoniae* | lineage | `cmh-4000` less default | every dataset | 3 | 2 / 5 / 17 of 24 | 0 |
| *Streptococcus pneumoniae* | lineage | `cmh-16000` less default | every dataset | 3 | 5 / 4 / 15 of 24 | 0 |
| *Streptococcus pneumoniae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 1 / 39 of 40 | 0 |
| *Streptococcus pneumoniae* | lineage | twin less default | every dataset | 5 | 5 / 3 / 32 of 40 | 0 |
| *Streptococcus pneumoniae* | lineage | twin less default | 800 genomes or more | 5 | 5 / 3 / 32 of 40 | 0 |
| *Streptococcus pneumoniae* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus pneumoniae* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Streptococcus pneumoniae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Vibrio cholerae* | lineage | `cmh-4000` less default | every dataset | 3 | 3 / 1 / 5 of 9 | 0 |
| *Vibrio cholerae* | lineage | `cmh-16000` less default | every dataset | 3 | 2 / 1 / 6 of 9 | 0 |
| *Vibrio cholerae* | lineage | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Vibrio cholerae* | lineage | twin less default | every dataset | 5 | 2 / 2 / 11 of 15 | 0 |
| *Vibrio cholerae* | lineage | twin less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |
| *Vibrio cholerae* | study | `cmh-4000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Vibrio cholerae* | study | `cmh-16000` less default | every dataset | 3 | 0 / 0 / 0 of 0 | 0 |
| *Vibrio cholerae* | study | `floor-3` less default | 800 genomes or more | 5 | 0 / 0 / 0 of 0 | 0 |

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**Single datasets.** The screen arms move single datasets much more than a fold seed does, in both directions. On `rf`, `gbt` and `lasso`, `cmh-4000` moves 136 of 696 lineage dataset and learner pairs by more than 0.05 in balanced accuracy, `cmh-16000` 124 of 696, and the twin 42 of 696. The largest moves are on *Staphylococcus aureus*: -0.367 with sulfamethoxazole/trimethoprim on `rf` and 0.259 with clindamycin on `lasso` for `cmh-4000`, and -0.386 with sulfamethoxazole/trimethoprim on `rf` and 0.308 with erythromycin on `lasso` for `cmh-16000`. The twin's moves run from -0.141 (*Klebsiella pneumoniae* with ticarcillin/clavulanate, `gbt`) to 0.169 (*Mycobacterium tuberculosis* with nicotinamide, `gbt`). The medians stay within 0.004 of 0.

**`rf` at 16,000.** At the budget of 16,000, `rf` is worse on 71 of 232 datasets and better on 36, while `gbt` and `lasso` are better on 59 and worse on 29 and 28. At 4,000 `rf` is better on 57 and worse on 50.

**Floor.** On the 92 datasets it can move, `floor-3` changes 123 of 460 dataset and learner figures by any amount, and moves from 1 (`cart`, `rf`, `lasso`) to 5 (`scm`) of 92 by more than 0.02.

**Study split.** On the study split, `cmh-16000`'s `lasso` has a median of 0.012 with an interval of 0.005 to 0.015 over 55 datasets, and its recall a median of 0.010 (0.002 to 0.016). The plan reads only the study median's sign, and the study split has one run per dataset and no twin.

**Fit time.** At 16,000 the outer fits of `lasso` take 9,471.663 seconds in all, against 2,187.286 at 4,000, and those of `gbt` 878.197 against 580.667. The default's and the twin's fit times are not in the table. The median peak memory runs from 1,659,816 to 1,692,756 KB across the three arms.

| Arm | Learner | Fits | Seconds, sum | Seconds, median | CPU seconds, sum | CPU seconds, median | Peak KB, median |
|---|---|---|---|---|---|---|---|
| `cmh-4000` | `scm` | 7,260 | 3273.690 | 0.187 | 45637.934 | 2.428 | 1,679,212 |
| `cmh-4000` | `cart` | 7,260 | 497.057 | 0.019 | 7281.892 | 0.271 | 1,679,212 |
| `cmh-4000` | `rf` | 7,260 | 155.504 | 0.010 | 2109.704 | 0.115 | 1,679,212 |
| `cmh-4000` | `gbt` | 7,260 | 580.667 | 0.050 | 6634.385 | 0.664 | 1,679,212 |
| `cmh-4000` | `lasso` | 7,260 | 2187.286 | 0.152 | 3640.277 | 0.271 | 1,679,212 |
| `cmh-16000` | `scm` | 7,260 | 3095.167 | 0.176 | 44661.596 | 2.401 | 1,659,816 |
| `cmh-16000` | `cart` | 7,260 | 463.814 | 0.018 | 7048.155 | 0.261 | 1,659,816 |
| `cmh-16000` | `rf` | 7,260 | 165.200 | 0.012 | 2202.136 | 0.128 | 1,659,816 |
| `cmh-16000` | `gbt` | 7,260 | 878.197 | 0.086 | 11185.376 | 1.137 | 1,659,816 |
| `cmh-16000` | `lasso` | 7,260 | 9471.663 | 0.567 | 10887.641 | 0.686 | 1,659,816 |
| `floor-3` | `scm` | 7,260 | 3358.757 | 0.191 | 48019.334 | 2.695 | 1,692,756 |
| `floor-3` | `cart` | 7,260 | 505.324 | 0.020 | 7631.238 | 0.277 | 1,692,756 |
| `floor-3` | `rf` | 7,260 | 44.285 | 0.004 | 397.542 | 0.032 | 1,692,756 |
| `floor-3` | `gbt` | 7,260 | 440.608 | 0.040 | 5094.471 | 0.541 | 1,692,756 |
| `floor-3` | `lasso` | 7,260 | 1912.862 | 0.143 | 2072.771 | 0.159 | 1,692,756 |

## Limits

**Study split.** One study run per dataset at the release's seeds, and no study twin. That departs from the plan's Rules, Draws, which give each arm seeds 42, 43 and 44.

**Binaries.** All 23 species ran one binary, code sha256 `5a39538947ab78961ac5c07e79adfb46851489f1782a466941df5c50aace3c63`, built from commit `2ee9e6bbe20096ef0601a6b1d36ee44c04d1e6ea`, with no restart. `reference/screen-arms/page-binaries.tsv` names each species' commit, binary and code sha256, and the restarts its run made.

**Labels.** The 35 datasets of M wait for the maintainer's label rule. The count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged.

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Markers.** No marker is found, named or ranked. The Cochran-Mantel-Haenszel screen adjusts the ranking for lineage strata and names no marker. The page holds aggregates over datasets, species, arms and learners only.

**Programs.** No AMRFinderPlus program is installed or run in this block.

## Method

**Checks of the run.** Before any score is read, `scripts/c4-screen/run_checks.py check` reads each species' copy and prints counts only. `checks.txt` in `reference/screen-arms/` holds its printout. Its checks:

| Check | What must hold |
|---|---|
| species | the copies hold the first run's species, each once |
| `DONE` | each species ended `DONE`; a `FAILED` older than its `DONE` is named and set aside |
| `SHA256SUMS` | each species' list holds every file of the run's rule, each at its sha256 |
| `score.log` | `score/score.log` stays out of `SHA256SUMS` and is never opened |
| `DONE` and `binary.txt` | each `DONE` names its species, its `binary.txt`'s binary and code, its results and its `SHA256SUMS` |
| one binary | every species ran one commit and one code sha256 |
| `lists.txt` | each names its species' first-run datasets, seeds 42, 43 and 44, twin 45, 46 and 47 and the three arms |
| `expect.py` | each `expect.tsv` is `expect.py`'s list from `scripts/c4-screen/lists.tsv` |
| results | each species holds every listed result or refusal record and inputs file, and no other |
| refusals | each refusal record is a study run scored as the constant model on five learners, and each constant-model row a refusal record's |
| genome ids | no genome id of the release's species lists in any file that will be read or published |

All 11 passed. The copies came to the reading without any `score/score.log` and without the run's `back.tar`. The printout:

```
PASS the first run's species, each once: 23 folders for 23 species; missing [], more []
PASS each species ended DONE, any FAILED older than its DONE: 23 of 23
PASS each species' SHA256SUMS holds every file of the run's rule, each at its sha256: 6207 files
PASS score/score.log stays out of SHA256SUMS
PASS each DONE names its species, its binary.txt's binary and code, its results and its SHA256SUMS
PASS every species ran one commit and one code sha256: 1 commits, 1 code sha256, 1 binaries
PASS each lists.txt names its species' first-run datasets, seeds 42, 43 and 44, twin 45, 46 and 47 and the three arms
PASS each expect.tsv is expect.py's list from lists.tsv: 3715 results
PASS each species holds every listed result or refusal record and inputs file, and no other
PASS each refusal record is a study run scored as the constant model on five learners, and each constant-model row a refusal record's: 0 refusal records, 0 constant-model rows
PASS no genome id in any file that will be read or published: 0 genome ids in 6253 files, against 23757 ids
CHECKS 11 of 11
```

**Reading.** The answers read are the run's own: each species' `score/` tables, copied as its run wrote them.

```
python3 -I scripts/c4-screen/run_checks.py check copies scripts/c4-screen/lists.tsv rev14-274-20261112 read/checks.txt
python3 -I scripts/c4-screen/run_checks.py gather copies scripts/c4-screen/lists.tsv read/checks.txt read
python3 -I scripts/c4-screen/score.py --count read/count read/score/*/datasets.tsv
python3 -I scripts/c4-screen/page_tables.py write read scripts/c4-screen/lists.tsv read/checks.txt reference/screen-arms
python3 -I scripts/c4-screen/page_tables.py print reference/screen-arms
```

`copies` holds the run's 23 species copies, each as `run_species.sh` left it. `rev14-274-20261112` is the rebuilt evaluation revision 14 release the run read, read only. `read/score/*/datasets.tsv` is the 23 species' files. Each command ran once, in this order, on 23 Nov 2026. The check and `gather` ran back to back, so the copies the checks read are the copies gathered.

| Command | Start | End | Last line |
|---|---|---|---|
| `run_checks.py check` | 14:24:09Z | 14:24:11Z | `CHECKS 11 of 11` |
| `run_checks.py gather` | 14:24:11Z | 14:24:11Z | `gathered 115 files` |
| `score.py --count` | 14:24:11Z | 14:24:53Z | `count: 6925 dataset rows from 23 files` |
| `page_tables.py write` | 14:24:53Z | 14:25:37Z | `gates: 16 of 16 pass; wrote 18 files` |
| `page_tables.py print` | 14:25:37Z | 14:25:37Z | the page's tables |

`read/count/withdrawn.txt` names each arm passed to `score.py --count` with `--withdrawn`, one a line, and is empty when none is. `reading.txt` lists every command with the sha256 of each file it wrote.

**Tables.** `scripts/c4-screen/page_tables.py write` checks the reading against 16 gates and writes `reference/screen-arms/` only when all pass. Where `score.py`'s tables already are the page's, it checks and copies them.

| Gate | What must hold |
|---|---|
| 1. files | the reading holds reading.txt, count/ with score.py --count's four tables and withdrawn.txt, and score/; no score.log is in it; the run's checks.txt ends CHECKS 11 of 11 with no FAIL line |
| 2. species | score/ holds each species of lists.tsv's first run once and no other, each with its five files |
| 3. datasets | each datasets.tsv holds exactly its species' first-run datasets, none in M: the default, each arm not withdrawn and the twin on the lineage split, and the default and each arm not withdrawn on the study split for those in S, each for the five learners; species, antibiotic and genome count equal lists.tsv's, and the floor arm moves the datasets of 800 genomes or more |
| 4. study | every study row is a study run at seed 42, and no study row stands on the lineage split |
| 5. twin | the twin stands on the lineage split at seeds 45, 46 and 47 alone, counted as the pair twin-k31 over every dataset and over those of 800 genomes or more, and conditions.tsv reads it as no arm |
| 6. missing | a draw not scored has no figure, and a dataset with a draw not scored is not scored and has no figure: a missing draw is never read as zero |
| 7. constant | each constant-model draw is an arm's study run refused for no candidate, scored and marked, and count-no-constant.tsv counts none |
| 8. draws | each draws.tsv holds each dataset's draws at seeds 42, 43 and 44, the twin's at 45, 46 and 47 and the study run at seed 42, each for the five learners, and nothing else; datasets.tsv is their mean |
| 9. checks | every check of each checks.tsv passed, and the figures agree: a screen arm's scm and cart, and the floor arm below 800 genomes, equal the default's in datasets.tsv |
| 10. screen | count.tsv and count-no-constant.tsv count a screen arm on rf, gbt and lasso only |
| 11. floor | count.tsv counts the floor arm over the datasets of 800 genomes or more only, the ones lists.tsv says it moves; count-floor-all.tsv holds the floor arm over every dataset with no Holm p, and decides nothing |
| 12. withdrawn | withdrawn.txt names arms of the block only, and the Holm family is the 21 lineage pairs score.py holds, a withdrawn arm's at p = 1 before Holm |
| 13. binary | each binary.txt names a commit, a binary and a code sha256; the page names each binary with its species unless every species ran one code sha256 |
| 14. conditions | conditions 1, 2, 3 and 4 of each screen arm, and 1 and 2, 3 and 4 of the floor arm's ruling 15, as read here off count.tsv, equal conditions.tsv's |
| 15. count | count.tsv, count-no-constant.tsv, conditions.tsv and count-floor-all.tsv are what score.py --count writes from the datasets.tsv files with withdrawn.txt's arms |
| 16. fits | each fits.tsv is fit_times.py's table of the arms not withdrawn and names neither the default nor the twin |

All 16 passed. `score.py`'s checks of every result, and the table script's own check of the figures in `datasets.tsv`, as checks passed of checks made:

| Source | Check | Passed |
|---|---|---|
| score.py checks.tsv | a lineage result names the default's pick | 3,490 of 3,490 |
| score.py checks.tsv | a study result names no default of its own | 220 of 220 |
| score.py checks.tsv | a twin seed is the twin's and no draw's | 699 of 699 |
| score.py checks.tsv | an unavailable split has no fold | 5 of 5 |
| score.py checks.tsv | below 800 genomes the floor arm's result is the default's, timing apart | 454 of 454 |
| score.py checks.tsv | cart: the outer calls give the split's balanced accuracy | 3,710 of 3,710 |
| score.py checks.tsv | dataset, scheme, seed and inner folds | 3,715 of 3,715 |
| score.py checks.tsv | each fold's outer calls are its test genomes' | 3,710 of 3,710 |
| score.py checks.tsv | exists | 3,715 of 3,715 |
| score.py checks.tsv | fold 0: no inner fit | 3,710 of 3,710 |
| score.py checks.tsv | fold 1: no inner fit | 3,710 of 3,710 |
| score.py checks.tsv | fold 2: no inner fit | 3,710 of 3,710 |
| score.py checks.tsv | fold 3: no inner fit | 3,710 of 3,710 |
| score.py checks.tsv | fold 4: no inner fit | 3,710 of 3,710 |
| score.py checks.tsv | fold tests sum to the tested genomes | 3,710 of 3,710 |
| score.py checks.tsv | folds in order | 3,710 of 3,710 |
| score.py checks.tsv | format and revision | 3,715 of 3,715 |
| score.py checks.tsv | gbt: the outer calls give the split's balanced accuracy | 3,710 of 3,710 |
| score.py checks.tsv | inputs: a fit time for each learner of each outer fold | 2,262 of 2,262 |
| score.py checks.tsv | inputs: a fold file's run fits no final model | 165 of 165 |
| score.py checks.tsv | inputs: each adjusted screen names its strata, no fewer than the informative ones | 1,508 of 1,508 |
| score.py checks.tsv | inputs: each outer screen reads its fold's candidates, rf its training genomes, gbt and lasso no more | 2,262 of 2,262 |
| score.py checks.tsv | inputs: each screen keeps the budget, or every candidate if fewer | 2,262 of 2,262 |
| score.py checks.tsv | inputs: every fit is of the run's splits | 2,262 of 2,262 |
| score.py checks.tsv | inputs: every screen by the arm's rule | 2,262 of 2,262 |
| score.py checks.tsv | inputs: format, revision, dataset, tag, seed and arm | 2,262 of 2,262 |
| score.py checks.tsv | inputs: rf, gbt and lasso screen once in each outer fold | 2,262 of 2,262 |
| score.py checks.tsv | inputs: the arm's floor, and the result's | 2,262 of 2,262 |
| score.py checks.tsv | inputs: the arm's rule and budget | 2,262 of 2,262 |
| score.py checks.tsv | inputs: the chi-square screens name no strata | 754 of 754 |
| score.py checks.tsv | inputs: the final models screen once each on the dataset's rows, rf on every genome | 2,097 of 2,097 |
| score.py checks.tsv | inputs: the result's outer folds, with its training genomes | 2,262 of 2,262 |
| score.py checks.tsv | lasso: the outer calls give the split's balanced accuracy | 3,710 of 3,710 |
| score.py checks.tsv | release report | 3,715 of 3,715 |
| score.py checks.tsv | rf: the outer calls give the split's balanced accuracy | 3,710 of 3,710 |
| score.py checks.tsv | scm and cart, which read no screen, make the default's calls in every fold | 1,508 of 1,508 |
| score.py checks.tsv | scm: the outer calls give the split's balanced accuracy | 3,710 of 3,710 |
| score.py checks.tsv | seed 42: each learner's lineage figure is the release report's | 233 of 233 |
| score.py checks.tsv | seed 42: the default is held to the release report | 233 of 233 |
| score.py checks.tsv | seed 42: the default's learner is the one the release serves | 233 of 233 |
| score.py checks.tsv | seed 42: the default's random-split figure is the release report's | 233 of 233 |
| score.py checks.tsv | the default and the twin write no inputs file | 1,453 of 1,453 |
| score.py checks.tsv | the default's learner is the later of its tied set, as the release breaks a tie | 3,490 of 3,490 |
| score.py checks.tsv | the five learners in order | 3,715 of 3,715 |
| score.py checks.tsv | the folds are the default's at the same dataset, split and seed, by class | 2,262 of 2,262 |
| score.py checks.tsv | the result's scheme is its split's | 3,715 of 3,715 |
| score.py checks.tsv | tie rule | 3,710 of 3,710 |
| datasets.tsv | a screen arm's scm and cart figures equal the default's | 1,152 of 1,152 |
| datasets.tsv | below 800 genomes the floor arm's figures equal the default's | 860 of 860 |

**Tests.** `tests/c4_screen/test_page_tables.py` holds SYNTHETIC tests of the table script and `tests/c4_screen/test_run_checks.py` of the check script, on three made-up species. `scripts/c4-screen/page_faults.py` plants each fault below in memory and runs the test meant to catch it, which must pass as the script stands and fail with the fault; `tests/c4_screen/test_page_faults.py` holds every fault to that, and each fault's text to one place in its file.

| Fault | Planted | Gate | Caught by |
|---|---|---|---|
| P1 | a score.log in the reading passes | 1 files | `Gates.test_a_score_log_in_the_reading_is_refused` |
| P24 | a printout of fewer checks than 11 passes | 1 files | `Gates.test_a_printout_of_fewer_checks_than_eleven_is_refused` |
| P2 | the species are those found, not those of the first run | 2 species | `Gates.test_a_species_missing_from_the_reading_is_refused` |
| P3 | a dataset more than the first run's, as one of M, passes | 3 datasets | `Gates.test_a_dataset_of_m_in_the_reading_is_refused` |
| P4 | a study row on the lineage split passes | 4 study | `Gates.test_a_study_row_on_the_lineage_split_is_refused` |
| P5 | the twin read as an arm passes | 5 twin | `Gates.test_the_twin_read_as_an_arm_is_refused` |
| P6 | a draw not scored keeps a figure | 6 missing | `Gates.test_a_draw_not_scored_with_a_figure_is_refused` |
| P7 | a constant-model draw not marked passes | 7 constant | `Gates.test_a_constant_model_draw_not_marked_is_refused` |
| P8 | datasets.tsv held to its draws on the names alone | 8 draws | `Gates.test_a_dataset_figure_that_is_not_the_mean_of_its_draws_is_refused` |
| P9 | a check of the run that failed passes | 9 checks | `Gates.test_a_check_of_the_run_that_failed_is_refused` |
| P10 | a screen arm's cart held to the default's on whether it is scored alone | 9 checks | `CartOff.test_a_screen_arm_whose_cart_differs_from_the_default_is_refused` |
| P11 | a screen arm counted on any learner | 10 screen | `Gates.test_a_count_row_of_a_screen_arm_on_scm_is_refused` |
| P12 | the floor arm counted over every dataset passes | 11 floor | `Gates.test_the_floor_arm_counted_over_every_dataset_is_refused` |
| P13 | a withdrawn arm's pairs dropped from the family pass | 12 withdrawn | `Withdrawn.test_a_withdrawn_arm_dropped_from_the_family_is_refused` |
| P14 | a binary.txt with no code sha256 passes | 13 binary | `Gates.test_a_binary_with_no_code_sha256_is_refused` |
| P15 | two binaries written as one | 13 binary | `Binaries.test_two_binaries_are_named_each_with_its_species` |
| P16 | a condition the count does not bear passes | 14 conditions | `Gates.test_a_condition_that_the_count_does_not_bear_is_refused` |
| P17 | the count held to score.py's by its length alone | 15 count | `Gates.test_a_count_figure_that_score_py_does_not_write_is_refused` |
| P18 | fit times of the twin pass | 16 fits | `Gates.test_fit_times_of_the_twin_are_refused` |
| P19 | a figure not scored read as zero in the dataset table | a table | `Tables.test_a_missing_draw_is_never_read_as_zero` |
| P20 | a refused study run not marked in the dataset table | a table | `Tables.test_a_refused_study_run_is_kept_and_marked` |
| P21 | the counts beside it marked as deciding | a table | `Tables.test_the_floor_arm_decides_over_800_genomes_or_more_and_floor_all_stands_beside` |
| P22 | the stricter reading marked as deciding | a table | `Tables.test_the_stricter_reading_decides_nothing_and_no_twin_arm_is_read` |
| P23 | a species named by its slug | a table | `Tables.test_species_are_named_in_full_and_their_counts_decide_nothing` |
| R1 | a species folder more passes | check script | `Check.test_a_species_folder_more_is_refused` |
| R2 | a FAILED newer than its DONE set aside | check script | `Check.test_a_failed_newer_than_done_is_refused` |
| R3 | SHA256SUMS held to the files' names alone | check script | `Check.test_a_file_that_differs_from_its_sha256_is_refused` |
| R4 | score.log in SHA256SUMS passes | check script | `Check.test_score_log_in_sha256sums_is_refused` |
| R5 | DONE's SHA256SUMS sha256 not read | check script | `Check.test_a_done_naming_other_sha256sums_is_refused` |
| R6 | two code sha256 pass | check script | `TwoBinaries.test_two_code_sha256_are_refused` |
| R7 | expect.tsv held to its header alone | check script | `Check.test_an_expect_tsv_that_is_not_expect_py_s_is_refused` |
| R8 | a result more than listed passes | check script | `Check.test_a_result_more_than_listed_is_refused` |
| R9 | a refusal record with no constant-model rows passes | check script | `Check.test_a_refusal_record_with_no_constant_model_rows_is_refused` |
| R10 | the genome id count skips the tables | check script | `Check.test_a_genome_id_in_a_table_is_refused` |
| R11 | gather copies score.log | check script | `Gather.test_gather_copies_five_files_a_species_and_never_score_log` |

**Files.** The reading's tables are in `reference/screen-arms/`. `scripts/c4-screen/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets

Each dataset of the first run is named with its species and antibiotic in full, its genome count and whether the floor arm moves it, with each arm's difference from the default on each learner it is counted on, and a study run refused, if any. `reference/screen-arms/page-datasets.tsv` holds each figure with the draws scored.

Each cell is the arm's balanced accuracy less the default's, the mean over the draws scored. The floor arm is counted only where it moves the dataset; elsewhere its figures equal the default's. `none` marks a split with a draw not scored, which has no figure. The twin's cell is empty on the study split, which has no twin.

| Species | Antibiotic | Genomes | Floor arm moves it | Split | `cmh-4000`: `rf` / `gbt` / `lasso` | `cmh-16000`: `rf` / `gbt` / `lasso` | `floor-3`: five learners | Twin: five learners | Study run refused |
|---|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | cefepime | 437 | no | lineage | 0.002 / 0.067 / 0.053 | 0.002 / 0.004 / -0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.005 / -0.007 / 0.000 / 0.029 / 0.053 | none |
| *Acinetobacter baumannii* | cefepime | 437 | no | study | -0.043 / -0.062 / 0.012 | -0.022 / -0.064 / 0.025 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Acinetobacter baumannii* | cefotaxime | 808 | yes | lineage | 0.052 / 0.043 / -0.061 | 0.055 / 0.103 / -0.009 | 0.020 / 0.000 / 0.000 / 0.000 / 0.000 | 0.011 / 0.026 / 0.001 / 0.031 / -0.009 | none |
| *Acinetobacter baumannii* | doripenem | 213 | no | lineage | 0.086 / 0.024 / -0.003 | 0.097 / 0.007 / 0.015 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.082 / -0.006 / -0.013 / -0.062 / -0.055 | none |
| *Acinetobacter baumannii* | doxycycline | 167 | no | lineage | 0.018 / 0.014 / 0.032 | 0.017 / 0.044 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.015 / 0.024 / 0.042 / 0.042 / -0.030 | none |
| *Acinetobacter baumannii* | minocycline | 139 | no | lineage | -0.024 / 0.053 / -0.007 | -0.039 / 0.004 / -0.030 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.072 / 0.027 / 0.000 / 0.000 / 0.014 | none |
| *Acinetobacter baumannii* | minocycline | 139 | no | study | 0.040 / 0.027 / -0.018 | -0.054 / -0.036 / -0.018 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Acinetobacter baumannii* | moxifloxacin | 91 | no | lineage | -0.055 / -0.027 / 0.000 | 0.000 / -0.016 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.058 / 0.000 / -0.011 / -0.045 / -0.085 | none |
| *Acinetobacter baumannii* | sulfamethoxazole/trimethoprim | 1,114 | yes | lineage | 0.117 / 0.110 / 0.134 | 0.082 / 0.126 / 0.140 | 0.000 / 0.004 / 0.000 / 0.000 / 0.000 | 0.027 / 0.009 / 0.000 / -0.018 / -0.024 | none |
| *Campylobacter coli* | azithromycin | 389 | no | lineage | 0.019 / 0.005 / 0.001 | -0.090 / 0.001 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.026 / 0.002 / -0.006 | none |
| *Campylobacter coli* | ciprofloxacin | 391 | no | lineage | -0.029 / 0.000 / 0.000 | -0.150 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.012 / 0.000 / 0.000 | none |
| *Campylobacter coli* | clindamycin | 389 | no | lineage | 0.048 / 0.003 / 0.001 | -0.057 / 0.001 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.001 / 0.006 / 0.002 / -0.010 / -0.006 | none |
| *Campylobacter coli* | erythromycin | 391 | no | lineage | 0.073 / 0.003 / 0.000 | -0.079 / 0.003 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.046 / 0.003 / -0.006 | none |
| *Campylobacter coli* | gentamicin | 391 | no | lineage | -0.069 / 0.057 / 0.010 | -0.049 / 0.065 / 0.022 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.015 / -0.039 / -0.004 / 0.088 / 0.000 | none |
| *Campylobacter coli* | nalidixic acid | 375 | no | lineage | -0.067 / 0.000 / 0.002 | -0.141 / 0.000 / -0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.024 / 0.000 / 0.003 | none |
| *Campylobacter coli* | telithromycin | 374 | no | lineage | -0.030 / 0.003 / -0.006 | -0.079 / -0.006 / -0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / -0.001 / 0.009 / -0.021 / 0.001 | none |
| *Campylobacter coli* | tetracycline | 391 | no | lineage | 0.000 / -0.003 / 0.002 | 0.000 / -0.002 / 0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.001 / -0.006 / 0.000 / -0.001 / 0.002 | none |
| *Campylobacter jejuni* | ciprofloxacin | 508 | no | lineage | -0.028 / 0.000 / 0.003 | -0.186 / -0.001 / 0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.016 / 0.000 / 0.004 | none |
| *Campylobacter jejuni* | clindamycin | 508 | no | lineage | -0.001 / -0.077 / -0.143 | 0.000 / 0.002 / 0.063 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.018 / 0.040 / 0.000 / 0.036 / 0.060 | none |
| *Campylobacter jejuni* | nalidixic acid | 506 | no | lineage | -0.045 / 0.000 / 0.005 | -0.167 / -0.004 / 0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.024 / 0.000 / 0.005 | none |
| *Campylobacter jejuni* | tetracycline | 508 | no | lineage | 0.000 / 0.000 / -0.002 | 0.000 / 0.000 / 0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.002 | none |
| *Clostridioides difficile* | azithromycin | 203 | no | lineage | -0.048 / -0.071 / -0.055 | -0.027 / -0.013 / -0.011 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.007 / -0.010 / 0.010 / 0.001 / 0.004 | none |
| *Clostridioides difficile* | ceftriaxone | 136 | no | lineage | -0.057 / -0.077 / -0.102 | -0.026 / -0.031 / -0.016 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.004 / -0.025 / 0.005 / -0.007 / 0.007 | none |
| *Clostridioides difficile* | clarithromycin | 203 | no | lineage | -0.042 / -0.079 / -0.057 | -0.027 / -0.027 / -0.011 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.006 / 0.015 / 0.002 / -0.020 / -0.006 | none |
| *Clostridioides difficile* | moxifloxacin | 204 | no | lineage | -0.030 / -0.018 / -0.016 | -0.027 / -0.017 / -0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / -0.004 | none |
| *Corynebacterium diphtheriae* | sulfamethoxazole/trimethoprim | 236 | no | lineage | 0.036 / 0.146 / 0.001 | 0.000 / 0.093 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.040 / 0.000 / 0.000 / 0.048 / -0.056 | none |
| *Corynebacterium diphtheriae* | tetracycline | 243 | no | lineage | -0.007 / 0.001 / -0.055 | -0.032 / -0.021 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.016 / -0.055 / 0.008 / 0.019 / -0.017 | none |
| *Corynebacterium diphtheriae* | trimethoprim | 229 | no | lineage | 0.000 / 0.033 / -0.084 | 0.000 / 0.000 / -0.040 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.006 / 0.013 / 0.000 / 0.029 / 0.022 | none |
| *Enterobacter cloacae* | cefepime | 103 | no | lineage | 0.000 / 0.001 / -0.092 | -0.078 / -0.002 / -0.035 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.009 / -0.005 / -0.007 / 0.001 / 0.047 | none |
| *Enterobacter cloacae* | cefepime | 103 | no | study | -0.037 / -0.021 / -0.022 | -0.029 / 0.001 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | ciprofloxacin | 110 | no | lineage | -0.032 / 0.015 / 0.044 | -0.069 / 0.003 / -0.117 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.019 / -0.030 / 0.037 / 0.048 / -0.001 | none |
| *Enterobacter cloacae* | ciprofloxacin | 110 | no | study | -0.016 / 0.019 / -0.041 | -0.015 / -0.034 / -0.039 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | gentamicin | 104 | no | lineage | -0.012 / 0.021 / 0.014 | -0.029 / -0.002 / 0.024 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.007 / -0.011 / -0.016 / 0.030 / 0.012 | none |
| *Enterobacter cloacae* | gentamicin | 104 | no | study | -0.014 / 0.006 / -0.083 | -0.025 / -0.019 / -0.064 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | imipenem | 86 | no | lineage | -0.079 / -0.042 / 0.068 | -0.099 / -0.073 / -0.023 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.013 / 0.062 / 0.092 / 0.014 / -0.015 | none |
| *Enterobacter cloacae* | imipenem | 86 | no | study | -0.061 / 0.071 / 0.158 | -0.017 / 0.097 / 0.079 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | levofloxacin | 104 | no | lineage | -0.013 / 0.005 / 0.101 | 0.045 / 0.070 / -0.077 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.036 / -0.022 / -0.046 / 0.009 / -0.092 | none |
| *Enterobacter cloacae* | levofloxacin | 104 | no | study | -0.006 / -0.017 / 0.019 | -0.006 / -0.098 / 0.017 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | meropenem | 99 | no | lineage | 0.004 / -0.016 / -0.017 | 0.004 / -0.015 / 0.008 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / 0.020 / -0.001 / -0.014 / 0.057 | none |
| *Enterobacter cloacae* | meropenem | 99 | no | study | -0.017 / -0.013 / 0.054 | 0.000 / -0.042 / 0.038 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | 80 | no | lineage | 0.053 / 0.069 / 0.041 | 0.008 / 0.001 / 0.041 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.008 / 0.033 / -0.008 | none |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | 80 | no | study | 0.037 / 0.000 / 0.026 | -0.024 / 0.000 / 0.013 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | tobramycin | 99 | no | lineage | -0.102 / 0.137 / 0.015 | -0.121 / 0.026 / 0.009 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / -0.006 / 0.099 / 0.060 / 0.010 | none |
| *Enterobacter cloacae* | tobramycin | 99 | no | study | -0.009 / -0.013 / -0.057 | -0.009 / -0.013 / -0.011 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Enterobacter cloacae* | trimethoprim | 124 | no | lineage | 0.035 / 0.005 / 0.024 | -0.005 / 0.018 / 0.020 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.030 / -0.001 / -0.017 | none |
| *Enterococcus faecium* | ampicillin | 1,680 | yes | lineage | 0.015 / 0.003 / 0.000 | 0.008 / 0.003 / 0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.003 / 0.000 / -0.002 | none |
| *Enterococcus faecium* | clindamycin | 901 | yes | lineage | 0.052 / -0.063 / 0.037 | 0.047 / -0.061 / 0.042 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.039 / 0.061 / 0.018 / -0.022 / 0.036 | none |
| *Enterococcus faecium* | dalfopristin/quinupristin | 1,118 | yes | lineage | -0.002 / -0.018 / -0.008 | 0.000 / 0.088 / 0.005 | 0.015 / 0.000 / 0.000 / 0.000 / -0.005 | 0.021 / 0.031 / -0.004 / 0.062 / 0.000 | none |
| *Enterococcus faecium* | doxycycline | 345 | no | lineage | 0.004 / -0.010 / 0.014 | -0.001 / -0.002 / 0.023 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.037 / 0.005 / -0.006 / -0.004 / -0.002 | none |
| *Enterococcus faecium* | erythromycin | 1,151 | yes | lineage | 0.036 / 0.036 / 0.055 | 0.011 / 0.034 / 0.021 | 0.000 / 0.001 / 0.000 / 0.000 / 0.006 | 0.043 / 0.014 / -0.015 / 0.038 / -0.023 | none |
| *Enterococcus faecium* | gentamicin | 1,345 | yes | lineage | 0.022 / 0.002 / 0.002 | 0.008 / 0.005 / -0.015 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.031 / 0.000 / -0.010 / -0.011 / -0.015 | none |
| *Enterococcus faecium* | kanamycin | 944 | yes | lineage | -0.055 / -0.016 / 0.013 | -0.029 / 0.014 / 0.015 | 0.000 / 0.003 / 0.000 / 0.000 / 0.000 | 0.002 / 0.003 / -0.009 / -0.010 / 0.018 | none |
| *Enterococcus faecium* | streptomycin | 1,239 | yes | lineage | 0.001 / 0.011 / -0.045 | 0.007 / -0.008 / -0.016 | 0.000 / -0.027 / 0.000 / 0.000 / 0.000 | 0.034 / -0.020 / -0.012 / 0.013 / -0.010 | none |
| *Enterococcus faecium* | teicoplanin | 1,192 | yes | lineage | 0.006 / 0.003 / 0.007 | 0.004 / 0.000 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.001 / 0.000 / 0.002 / 0.000 / 0.003 | none |
| *Enterococcus faecium* | tetracycline | 958 | yes | lineage | -0.022 / -0.022 / 0.007 | -0.012 / 0.005 / 0.005 | -0.002 / 0.000 / 0.000 / 0.000 / 0.000 | -0.007 / -0.011 / -0.006 / 0.005 / 0.019 | none |
| *Enterococcus faecium* | vancomycin | 1,741 | yes | lineage | 0.000 / -0.011 / 0.000 | 0.001 / 0.003 / -0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / -0.001 / 0.000 | none |
| *Escherichia coli* | amikacin | 1,125 | yes | lineage | -0.033 / 0.061 / -0.004 | -0.038 / -0.033 / 0.006 | -0.022 / 0.001 / 0.011 / -0.018 / -0.049 | 0.000 / 0.020 / 0.006 / -0.007 / 0.000 | none |
| *Escherichia coli* | amikacin | 1,125 | yes | study | -0.069 / 0.010 / 0.005 | -0.067 / 0.013 / 0.005 | -0.050 / 0.000 / 0.000 / -0.011 / -0.030 |  | none |
| *Escherichia coli* | ampicillin/sulbactam | 283 | no | lineage | 0.009 / 0.003 / -0.097 | -0.032 / 0.056 / -0.061 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.018 / -0.012 / -0.027 / 0.059 / 0.012 | none |
| *Escherichia coli* | ampicillin/sulbactam | 283 | no | study | -0.014 / -0.014 / 0.025 | 0.014 / 0.007 / 0.056 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | aztreonam | 699 | no | lineage | 0.172 / 0.033 / 0.036 | 0.103 / 0.021 / 0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.044 / 0.007 / -0.006 / 0.002 / 0.000 | none |
| *Escherichia coli* | aztreonam | 699 | no | study | 0.007 / 0.037 / -0.003 | 0.016 / 0.029 / 0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | cefazolin | 414 | no | lineage | 0.114 / 0.053 / -0.026 | 0.084 / 0.060 / -0.099 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.038 / -0.017 / 0.007 / -0.007 / -0.005 | none |
| *Escherichia coli* | cefazolin | 414 | no | study | 0.080 / 0.004 / 0.036 | 0.019 / 0.002 / 0.017 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | cefepime | 1,124 | yes | lineage | 0.005 / -0.003 / 0.048 | -0.022 / -0.081 / 0.035 | -0.005 / 0.000 / 0.000 / 0.000 / 0.000 | 0.027 / -0.060 / 0.004 / 0.006 / 0.089 | none |
| *Escherichia coli* | cefepime | 1,124 | yes | study | -0.001 / 0.001 / 0.008 | -0.011 / 0.001 / 0.010 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | cefotaxime | 1,149 | yes | lineage | 0.006 / 0.016 / 0.003 | 0.008 / 0.011 / 0.012 | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 | 0.005 / -0.003 / -0.001 / 0.003 / -0.001 | none |
| *Escherichia coli* | cefoxitin | 1,066 | yes | lineage | 0.007 / 0.018 / 0.007 | 0.018 / 0.031 / 0.010 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.002 / 0.006 / 0.000 / 0.001 / 0.001 | none |
| *Escherichia coli* | cefoxitin | 1,066 | yes | study | -0.001 / 0.020 / -0.001 | 0.004 / 0.034 / 0.008 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | cefpodoxime | 175 | no | lineage | 0.001 / 0.002 / 0.004 | 0.014 / -0.006 / -0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.065 / 0.007 / 0.017 / 0.001 / 0.006 | none |
| *Escherichia coli* | cefpodoxime/clavulanate | 136 | no | lineage | 0.043 / -0.018 / 0.004 | 0.007 / -0.060 / -0.052 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.068 / 0.017 / 0.011 / -0.024 / 0.031 | none |
| *Escherichia coli* | ceftazidime | 1,466 | yes | lineage | 0.019 / -0.006 / -0.003 | 0.013 / 0.009 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.004 / -0.003 / 0.000 / -0.001 / -0.006 | none |
| *Escherichia coli* | ceftiofur | 203 | no | lineage | 0.029 / -0.011 / -0.155 | -0.131 / 0.004 / -0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.010 / 0.029 / -0.003 / -0.032 / -0.002 | none |
| *Escherichia coli* | ceftriaxone | 923 | yes | lineage | 0.010 / 0.000 / 0.031 | 0.006 / -0.023 / 0.032 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.018 / 0.002 / -0.015 / -0.011 / -0.013 | none |
| *Escherichia coli* | ceftriaxone | 923 | yes | study | -0.001 / 0.009 / 0.002 | -0.001 / 0.015 / 0.020 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | cefuroxime | 548 | no | lineage | 0.008 / 0.025 / 0.035 | -0.003 / 0.025 / 0.030 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.002 / 0.000 / 0.000 / 0.003 / -0.006 | none |
| *Escherichia coli* | cephalothin | 472 | no | lineage | 0.001 / 0.009 / -0.003 | -0.013 / -0.001 / -0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / -0.023 / 0.003 / 0.013 / -0.001 | none |
| *Escherichia coli* | chloramphenicol | 709 | no | lineage | -0.019 / 0.005 / 0.021 | -0.036 / 0.011 / 0.022 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.001 / 0.000 / -0.003 / 0.007 / 0.008 | none |
| *Escherichia coli* | colistin | 239 | no | lineage | 0.000 / 0.031 / 0.020 | -0.032 / 0.007 / 0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / -0.008 | none |
| *Escherichia coli* | doxycycline | 227 | no | lineage | 0.004 / -0.017 / -0.032 | -0.033 / 0.019 / -0.044 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.014 / -0.007 / 0.020 / 0.002 / 0.038 | none |
| *Escherichia coli* | ertapenem | 820 | yes | lineage | 0.032 / 0.041 / 0.036 | 0.005 / 0.049 / 0.054 | 0.000 / -0.004 / 0.000 / 0.000 / 0.000 | 0.006 / 0.001 / 0.010 / 0.000 / 0.011 | none |
| *Escherichia coli* | ertapenem | 820 | yes | study | 0.005 / 0.077 / 0.072 | -0.033 / 0.062 / 0.097 | 0.007 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | erythromycin | 205 | no | lineage | 0.000 / 0.000 / 0.000 | 0.000 / 0.004 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.001 / -0.010 | none |
| *Escherichia coli* | imipenem | 614 | no | lineage | -0.011 / 0.002 / 0.001 | -0.076 / -0.003 / 0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.003 / 0.000 / -0.005 / -0.007 | none |
| *Escherichia coli* | levofloxacin | 403 | no | lineage | 0.079 / 0.008 / 0.025 | -0.008 / -0.002 / 0.021 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.011 / 0.009 / 0.001 / 0.010 / 0.011 | none |
| *Escherichia coli* | levofloxacin | 403 | no | study | 0.027 / 0.009 / -0.002 | -0.018 / 0.005 / -0.009 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | nalidixic acid | 406 | no | lineage | 0.017 / 0.003 / 0.007 | 0.005 / -0.004 / 0.010 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.002 / 0.001 / 0.008 / -0.014 / -0.002 | none |
| *Escherichia coli* | nitrofurantoin | 482 | no | lineage | 0.014 / 0.041 / -0.038 | -0.025 / 0.019 / -0.034 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.021 / -0.017 / -0.007 / -0.022 / 0.013 | none |
| *Escherichia coli* | nitrofurantoin | 482 | no | study | -0.001 / -0.008 / -0.019 | -0.001 / -0.005 / -0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | piperacillin/tazobactam | 979 | yes | lineage | 0.000 / -0.024 / -0.027 | -0.004 / -0.023 / 0.048 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.002 / -0.004 / 0.001 / -0.007 / -0.016 | none |
| *Escherichia coli* | piperacillin/tazobactam | 979 | yes | study | 0.001 / 0.013 / 0.020 | -0.010 / -0.008 / 0.070 | 0.000 / -0.002 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | streptomycin | 386 | no | lineage | -0.020 / 0.005 / -0.018 | -0.012 / 0.000 / -0.030 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.014 / 0.000 / 0.006 / -0.009 / -0.011 | none |
| *Escherichia coli* | streptomycin | 386 | no | study | -0.020 / 0.005 / -0.002 | -0.050 / -0.012 / -0.053 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | sulfamethoxazole | 341 | no | lineage | 0.002 / -0.004 / -0.022 | -0.004 / -0.001 / -0.021 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / 0.026 / -0.002 / 0.009 / -0.008 | none |
| *Escherichia coli* | sulfamethoxazole | 341 | no | study | 0.027 / 0.018 / -0.027 | 0.007 / 0.023 / -0.019 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | tetracycline | 782 | no | lineage | 0.007 / 0.014 / -0.003 | 0.003 / 0.018 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.004 / -0.007 / -0.003 / 0.012 / -0.004 | none |
| *Escherichia coli* | tetracycline | 782 | no | study | 0.000 / -0.010 / -0.003 | -0.020 / -0.007 / -0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | tobramycin | 1,070 | yes | lineage | 0.053 / -0.007 / 0.009 | 0.032 / 0.013 / 0.007 | 0.000 / 0.004 / 0.000 / 0.000 / 0.000 | -0.005 / 0.000 / -0.014 / -0.003 / -0.006 | none |
| *Escherichia coli* | tobramycin | 1,070 | yes | study | 0.006 / 0.019 / 0.017 | -0.005 / 0.017 / 0.017 | 0.000 / 0.007 / 0.000 / 0.000 / 0.000 |  | none |
| *Escherichia coli* | trimethoprim | 861 | yes | lineage | 0.051 / 0.019 / 0.004 | 0.046 / 0.009 / 0.018 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.009 / 0.005 / 0.011 / -0.001 / 0.003 | none |
| *Klebsiella aerogenes* | ertapenem | 83 | no | lineage | 0.093 / 0.011 / 0.028 | 0.052 / -0.130 / -0.037 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.048 / 0.027 / 0.087 / -0.013 / 0.081 | none |
| *Klebsiella aerogenes* | imipenem | 68 | no | lineage | 0.004 / 0.002 / -0.023 | 0.008 / -0.005 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.026 / -0.034 / -0.008 / -0.004 / -0.032 | none |
| *Klebsiella aerogenes* | meropenem | 78 | no | lineage | -0.003 / -0.022 / -0.013 | -0.003 / 0.000 / -0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.071 / -0.019 / 0.006 / -0.010 / 0.160 | none |
| *Klebsiella oxytoca* | cefotaxime | 98 | no | lineage | -0.063 / -0.084 / -0.108 | 0.022 / -0.030 / -0.028 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.045 / 0.036 / 0.026 / 0.053 / 0.000 | none |
| *Klebsiella oxytoca* | cefoxitin | 99 | no | lineage | 0.004 / -0.024 / 0.016 | -0.003 / 0.003 / 0.050 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.015 / -0.042 / -0.006 / -0.022 / 0.032 | none |
| *Klebsiella oxytoca* | cefpodoxime | 74 | no | lineage | 0.255 / 0.113 / 0.224 | -0.047 / -0.050 / -0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.068 / -0.021 / -0.048 / -0.019 / -0.050 | none |
| *Klebsiella oxytoca* | ertapenem | 112 | no | lineage | 0.066 / 0.085 / -0.011 | 0.084 / 0.133 / 0.020 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.082 / -0.052 / -0.046 / 0.140 / 0.018 | none |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 744 | no | lineage | -0.005 / -0.034 / -0.048 | 0.004 / 0.008 / 0.014 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.031 / -0.008 / -0.001 / 0.004 / -0.012 | none |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 744 | no | study | -0.044 / 0.031 / 0.097 | -0.035 / -0.039 / 0.080 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | aztreonam | 1,050 | yes | lineage | 0.022 / -0.019 / -0.003 | 0.026 / -0.021 / 0.009 | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 | 0.021 / 0.013 / -0.005 / -0.010 / -0.048 | none |
| *Klebsiella pneumoniae* | aztreonam | 1,050 | yes | study | 0.034 / 0.022 / 0.045 | 0.024 / 0.056 / 0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | cefazolin | 1,141 | yes | lineage | 0.006 / 0.010 / 0.035 | 0.011 / 0.013 / 0.059 | 0.000 / 0.010 / 0.000 / 0.000 / 0.000 | -0.004 / -0.009 / 0.001 / 0.014 / 0.032 | none |
| *Klebsiella pneumoniae* | cefazolin | 1,141 | yes | study | 0.008 / 0.001 / 0.012 | 0.008 / 0.009 / 0.016 | 0.000 / -0.009 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | cefepime | 1,233 | yes | lineage | 0.006 / 0.051 / -0.006 | -0.005 / -0.005 / -0.009 | 0.000 / -0.013 / 0.000 / 0.000 / 0.000 | -0.040 / -0.037 / 0.007 / -0.018 / -0.002 | none |
| *Klebsiella pneumoniae* | cefepime | 1,233 | yes | study | 0.006 / 0.004 / -0.024 | 0.007 / -0.006 / -0.004 | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | cefotetan | 61 | no | lineage | 0.000 / 0.049 / 0.071 | -0.007 / 0.074 / 0.157 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.014 / 0.098 / -0.020 / -0.024 / 0.012 | none |
| *Klebsiella pneumoniae* | cefoxitin | 1,181 | yes | lineage | 0.063 / 0.048 / 0.062 | 0.051 / 0.054 / 0.063 | 0.001 / 0.000 / 0.000 / 0.000 / 0.000 | -0.046 / -0.094 / 0.000 / -0.007 / -0.006 | none |
| *Klebsiella pneumoniae* | cefoxitin | 1,181 | yes | study | 0.024 / 0.023 / 0.010 | 0.011 / 0.014 / 0.005 | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | cefpodoxime | 221 | no | lineage | 0.120 / 0.077 / 0.093 | 0.079 / 0.066 / 0.019 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.017 / -0.047 / 0.014 / -0.004 / 0.076 | none |
| *Klebsiella pneumoniae* | cefpodoxime/ceftriaxone | 56 | no | lineage | 0.042 / -0.048 / 0.036 | 0.089 / -0.018 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.018 / 0.024 / 0.012 / 0.030 / 0.036 | none |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | 207 | no | lineage | 0.148 / 0.080 / 0.049 | 0.070 / 0.130 / -0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.099 / -0.036 / -0.006 / 0.015 / -0.019 | none |
| *Klebsiella pneumoniae* | ceftazidime | 1,310 | yes | lineage | 0.030 / 0.040 / 0.050 | 0.033 / 0.061 / 0.047 | 0.000 / -0.016 / 0.000 / 0.000 / 0.000 | 0.008 / -0.005 / 0.005 / 0.003 / 0.009 | none |
| *Klebsiella pneumoniae* | ceftazidime | 1,310 | yes | study | 0.002 / 0.033 / 0.015 | 0.021 / 0.028 / 0.028 | 0.000 / 0.002 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | 259 | no | lineage | -0.006 / -0.150 / 0.152 | -0.001 / 0.125 / 0.017 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.014 / 0.013 / 0.000 / 0.003 / -0.122 | none |
| *Klebsiella pneumoniae* | ceftriaxone | 1,208 | yes | lineage | -0.003 / -0.002 / -0.003 | 0.012 / 0.028 / 0.022 | 0.000 / 0.010 / 0.000 / 0.000 / 0.000 | -0.001 / 0.000 / -0.005 / 0.002 / -0.002 | none |
| *Klebsiella pneumoniae* | ceftriaxone | 1,208 | yes | study | 0.009 / -0.003 / 0.024 | 0.011 / 0.027 / 0.074 | 0.000 / 0.001 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | cefuroxime | 565 | no | lineage | 0.016 / 0.021 / 0.040 | -0.010 / -0.005 / -0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.012 / -0.034 / -0.003 / 0.004 / 0.012 | none |
| *Klebsiella pneumoniae* | chloramphenicol | 409 | no | lineage | -0.051 / -0.015 / -0.009 | -0.060 / 0.004 / -0.017 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.019 / -0.001 / 0.001 / 0.021 / -0.034 | none |
| *Klebsiella pneumoniae* | ciprofloxacin | 1,706 | yes | lineage | -0.010 / 0.014 / -0.093 | 0.003 / 0.011 / -0.005 | 0.006 / -0.003 / 0.000 / 0.000 / 0.000 | 0.004 / -0.001 / -0.003 / 0.005 / -0.009 | none |
| *Klebsiella pneumoniae* | ciprofloxacin | 1,706 | yes | study | -0.015 / -0.040 / -0.014 | -0.001 / -0.004 / -0.003 | 0.000 / 0.002 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | colistin | 336 | no | lineage | -0.047 / -0.022 / 0.002 | -0.048 / -0.048 / 0.073 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.002 / -0.002 / -0.038 / -0.004 / -0.001 | none |
| *Klebsiella pneumoniae* | colistin | 336 | no | study | -0.016 / 0.025 / 0.038 | -0.049 / 0.013 / 0.018 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | doripenem | 243 | no | lineage | -0.002 / -0.104 / -0.024 | 0.016 / -0.019 / -0.009 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.002 / -0.017 / 0.000 / -0.004 / 0.016 | none |
| *Klebsiella pneumoniae* | doripenem | 243 | no | study | 0.007 / -0.098 / 0.075 | 0.005 / -0.104 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | doxycycline | 141 | no | lineage | 0.038 / 0.091 / -0.015 | 0.030 / 0.033 / -0.014 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.040 / -0.036 / -0.014 / 0.082 / 0.002 | none |
| *Klebsiella pneumoniae* | ertapenem | 952 | yes | lineage | 0.057 / 0.030 / 0.010 | 0.061 / 0.044 / 0.044 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.004 / 0.080 / 0.017 / -0.003 / -0.013 | none |
| *Klebsiella pneumoniae* | ertapenem | 952 | yes | study | 0.004 / 0.038 / -0.031 | 0.044 / 0.053 / 0.077 | 0.000 / 0.001 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | fosfomycin | 294 | no | lineage | -0.097 / -0.030 / -0.033 | -0.049 / 0.054 / 0.008 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.033 / 0.001 / 0.019 / 0.054 / 0.042 | none |
| *Klebsiella pneumoniae* | imipenem | 1,158 | yes | lineage | 0.063 / 0.067 / 0.055 | 0.076 / 0.065 / 0.056 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.001 / -0.013 / 0.005 / 0.000 / -0.005 | none |
| *Klebsiella pneumoniae* | imipenem | 1,158 | yes | study | 0.039 / 0.029 / 0.021 | 0.042 / 0.028 / 0.030 | -0.013 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | levofloxacin | 871 | yes | lineage | 0.020 / -0.003 / 0.020 | -0.025 / -0.002 / 0.017 | 0.001 / 0.000 / 0.000 / 0.000 / 0.000 | 0.013 / -0.001 / -0.002 / 0.003 / -0.002 | none |
| *Klebsiella pneumoniae* | levofloxacin | 871 | yes | study | -0.016 / -0.042 / 0.024 | -0.003 / 0.005 / 0.013 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | nitrofurantoin | 562 | no | lineage | -0.043 / -0.031 / 0.040 | -0.059 / -0.039 / 0.122 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.014 / 0.030 / -0.005 / -0.048 / 0.066 | none |
| *Klebsiella pneumoniae* | nitrofurantoin | 562 | no | study | 0.001 / 0.035 / 0.024 | 0.007 / 0.075 / 0.055 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | ofloxacin | 71 | no | lineage | 0.043 / 0.018 / 0.208 | 0.008 / -0.051 / -0.015 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.025 / 0.014 / 0.000 / -0.041 / -0.012 | none |
| *Klebsiella pneumoniae* | polymyxin b | 100 | no | lineage | none / none / none | none / none / none | none / none / none / none / none | none / none / none / none / none | none |
| *Klebsiella pneumoniae* | polymyxin b | 100 | no | study | 0.010 / -0.106 / -0.116 | 0.000 / -0.115 / -0.096 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | streptomycin | 80 | no | lineage | -0.045 / -0.049 / -0.019 | -0.045 / 0.018 / -0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.011 / 0.016 / -0.045 / -0.044 / 0.000 | none |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 1,546 | yes | lineage | 0.002 / -0.010 / -0.006 | -0.003 / -0.015 / -0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.005 / 0.000 / -0.002 / -0.001 / 0.002 | none |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 1,546 | yes | study | 0.008 / 0.001 / -0.003 | -0.005 / -0.006 / -0.023 | -0.006 / 0.008 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | tetracycline | 741 | no | lineage | 0.012 / 0.049 / 0.049 | -0.017 / 0.042 / 0.048 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.003 / -0.006 / 0.003 / -0.009 / -0.006 | none |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | 250 | no | lineage | 0.032 / -0.023 / 0.088 | 0.034 / 0.108 / 0.227 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.029 / -0.039 / -0.027 / -0.141 / 0.062 | none |
| *Klebsiella pneumoniae* | tigecycline | 924 | yes | lineage | -0.033 / -0.040 / 0.012 | -0.041 / -0.011 / 0.027 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.011 / -0.020 / 0.018 / -0.025 / 0.002 | none |
| *Klebsiella pneumoniae* | tigecycline | 924 | yes | study | 0.007 / -0.004 / -0.030 | 0.002 / 0.010 / -0.009 | 0.000 / -0.003 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | tobramycin | 1,230 | yes | lineage | 0.000 / -0.015 / 0.051 | -0.014 / -0.031 / 0.052 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.003 / 0.013 / -0.001 / -0.009 / -0.006 | none |
| *Klebsiella pneumoniae* | tobramycin | 1,230 | yes | study | 0.014 / -0.018 / 0.003 | 0.007 / 0.001 / 0.036 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Klebsiella pneumoniae* | trimethoprim | 258 | no | lineage | 0.017 / 0.004 / 0.007 | 0.016 / 0.004 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.004 / -0.001 / -0.014 / -0.003 / -0.026 | none |
| *Klebsiella quasipneumoniae* | cefotaxime | 101 | no | lineage | -0.010 / 0.007 / 0.004 | -0.010 / -0.003 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.019 / 0.032 / 0.013 / 0.011 / -0.002 | none |
| *Klebsiella quasipneumoniae* | ceftazidime | 102 | no | lineage | 0.000 / 0.011 / 0.002 | 0.000 / 0.016 / 0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.020 / 0.024 / 0.011 / 0.010 / 0.010 | none |
| *Klebsiella quasipneumoniae* | sulfamethoxazole/trimethoprim | 101 | no | lineage | 0.004 / 0.026 / -0.072 | 0.014 / 0.048 / -0.055 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.020 / -0.010 / 0.028 / 0.021 / -0.009 | none |
| *Mycobacterium tuberculosis* | amikacin | 1,813 | yes | lineage | 0.016 / -0.001 / 0.002 | -0.034 / 0.006 / 0.002 | 0.002 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / -0.003 / -0.002 / 0.005 / 0.000 | none |
| *Mycobacterium tuberculosis* | amoxicillin | 57 | no | lineage | -0.066 / -0.154 / -0.082 | -0.012 / -0.038 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.032 / -0.070 / -0.064 / 0.032 / 0.045 | none |
| *Mycobacterium tuberculosis* | capreomycin | 670 | no | lineage | 0.025 / 0.007 / -0.002 | -0.006 / 0.008 / 0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / 0.011 / -0.009 / -0.010 / 0.003 | none |
| *Mycobacterium tuberculosis* | ciprofloxacin | 101 | no | lineage | -0.032 / 0.030 / -0.022 | -0.013 / 0.000 / -0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.107 / 0.028 / -0.042 / 0.035 / -0.003 | none |
| *Mycobacterium tuberculosis* | cycloserine | 401 | no | lineage | -0.015 / -0.051 / -0.070 | -0.005 / -0.056 / -0.022 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.024 / -0.011 / -0.033 / -0.037 / -0.045 | none |
| *Mycobacterium tuberculosis* | delamanid | 1,248 | yes | lineage | -0.020 / 0.007 / -0.004 | -0.022 / -0.002 / -0.004 | 0.000 / 0.001 / -0.011 / 0.006 / -0.008 | -0.009 / 0.009 / -0.010 / 0.000 / -0.007 | none |
| *Mycobacterium tuberculosis* | ethambutol | 1,947 | yes | lineage | -0.012 / -0.003 / 0.003 | -0.019 / -0.003 / 0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.001 / -0.004 / -0.003 / -0.002 / -0.001 | none |
| *Mycobacterium tuberculosis* | ethionamide | 1,774 | yes | lineage | 0.016 / 0.002 / 0.027 | -0.003 / 0.000 / 0.029 | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 | -0.015 / -0.008 / -0.009 / -0.001 / 0.002 | none |
| *Mycobacterium tuberculosis* | kanamycin | 1,800 | yes | lineage | 0.013 / 0.013 / -0.004 | -0.044 / 0.018 / -0.003 | 0.000 / -0.004 / 0.000 / 0.000 / 0.000 | 0.004 / -0.001 / 0.006 / 0.001 / 0.000 | none |
| *Mycobacterium tuberculosis* | levofloxacin | 1,248 | yes | lineage | 0.001 / -0.006 / -0.006 | -0.018 / -0.004 / -0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.002 / -0.008 / -0.013 / 0.000 | none |
| *Mycobacterium tuberculosis* | linezolid | 1,248 | yes | lineage | -0.008 / 0.036 / 0.007 | -0.008 / 0.031 / 0.009 | -0.002 / -0.001 / -0.002 / 0.033 / -0.006 | -0.017 / -0.003 / 0.003 / 0.038 / -0.009 | none |
| *Mycobacterium tuberculosis* | moxifloxacin | 1,524 | yes | lineage | 0.012 / 0.001 / -0.002 | -0.015 / 0.002 / -0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.002 / 0.001 / 0.003 / -0.003 | none |
| *Mycobacterium tuberculosis* | nicotinamide | 166 | no | lineage | 0.004 / -0.004 / 0.000 | 0.004 / -0.004 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.014 / -0.012 / 0.002 / 0.169 / -0.010 | none |
| *Mycobacterium tuberculosis* | nicotinamide | 166 | no | study | -0.006 / -0.036 / 0.000 | -0.006 / -0.036 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Mycobacterium tuberculosis* | ofloxacin | 643 | no | lineage | -0.005 / 0.005 / 0.004 | -0.096 / 0.010 / 0.010 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.004 / 0.000 / -0.006 / -0.005 / 0.008 | none |
| *Mycobacterium tuberculosis* | para-aminosalicylic acid | 1,650 | yes | lineage | -0.009 / -0.004 / -0.016 | -0.015 / -0.002 / -0.015 | 0.000 / -0.004 / 0.000 / 0.000 / 0.000 | 0.002 / -0.025 / -0.001 / 0.004 / 0.003 | none |
| *Mycobacterium tuberculosis* | prothionamide | 170 | no | lineage | -0.019 / -0.005 / -0.018 | 0.008 / 0.008 / -0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.021 / 0.006 / 0.007 / -0.011 / 0.026 | none |
| *Mycobacterium tuberculosis* | pyrazinamide | 492 | no | lineage | -0.013 / -0.028 / -0.015 | 0.002 / -0.002 / 0.011 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.007 / -0.008 / -0.006 / -0.057 / 0.011 | none |
| *Mycobacterium tuberculosis* | rifabutin | 1,317 | yes | lineage | -0.005 / 0.004 / 0.001 | -0.011 / 0.003 / 0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.002 / 0.002 / -0.006 / 0.001 / 0.000 | none |
| *Mycobacterium tuberculosis* | rifampin | 1,967 | yes | lineage | 0.015 / 0.012 / 0.016 | -0.003 / 0.014 / 0.017 | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 | -0.001 / 0.001 / -0.007 / 0.002 / -0.003 | none |
| *Mycobacterium tuberculosis* | streptomycin | 698 | no | lineage | 0.010 / -0.023 / 0.025 | 0.017 / 0.008 / 0.010 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.007 / -0.027 / -0.015 / 0.030 / -0.017 | none |
| *Neisseria gonorrhoeae* | cefixime | 600 | no | lineage | -0.027 / 0.008 / 0.040 | -0.044 / 0.002 / 0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.004 / 0.004 / -0.006 / -0.002 / -0.003 | none |
| *Neisseria gonorrhoeae* | ceftriaxone | 1,926 | yes | lineage | 0.000 / 0.058 / 0.015 | 0.000 / 0.008 / 0.010 | -0.015 / 0.000 / 0.000 / 0.026 / 0.004 | -0.001 / -0.082 / 0.000 / -0.002 / -0.001 | none |
| *Neisseria gonorrhoeae* | ciprofloxacin | 1,746 | yes | lineage | -0.036 / 0.002 / -0.160 | -0.116 / 0.002 / -0.006 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.003 / 0.001 / -0.087 / 0.002 / 0.001 | none |
| *Neisseria gonorrhoeae* | erythromycin | 179 | no | lineage | -0.037 / -0.072 / -0.146 | -0.027 / -0.022 / -0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.012 / -0.004 / 0.021 / 0.037 / -0.032 | none |
| *Neisseria gonorrhoeae* | penicillin | 1,354 | yes | lineage | 0.168 / 0.169 / 0.066 | 0.129 / 0.254 / 0.105 | -0.022 / 0.012 / 0.000 / 0.000 / 0.000 | -0.037 / 0.037 / 0.003 / 0.022 / -0.014 | none |
| *Neisseria gonorrhoeae* | tetracycline | 1,371 | yes | lineage | 0.004 / -0.016 / -0.041 | -0.009 / -0.023 / -0.036 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.003 / -0.018 / -0.018 / -0.006 / -0.012 | none |
| *Pseudomonas aeruginosa* | amikacin | 708 | no | lineage | 0.006 / 0.007 / 0.034 | -0.043 / -0.019 / 0.048 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.011 / 0.008 / 0.000 / 0.019 / 0.004 | none |
| *Pseudomonas aeruginosa* | aztreonam | 235 | no | lineage | -0.086 / -0.029 / -0.063 | -0.070 / -0.011 / -0.020 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.041 / -0.038 / -0.029 / -0.041 / -0.019 | none |
| *Pseudomonas aeruginosa* | aztreonam | 235 | no | study | 0.023 / -0.029 / -0.077 | 0.030 / 0.000 / -0.037 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Pseudomonas aeruginosa* | cefepime | 245 | no | lineage | -0.067 / -0.044 / 0.075 | -0.003 / -0.039 / 0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.041 / -0.035 / 0.026 / 0.022 / 0.014 | none |
| *Pseudomonas aeruginosa* | cefepime | 245 | no | study | -0.004 / -0.041 / 0.118 | -0.008 / -0.020 / 0.041 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Pseudomonas aeruginosa* | ceftazidime | 833 | yes | lineage | 0.011 / 0.036 / 0.011 | 0.012 / 0.032 / 0.016 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.005 / 0.011 / -0.006 / 0.017 / -0.019 | none |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 309 | no | lineage | 0.022 / -0.033 / 0.012 | -0.034 / 0.012 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.023 / 0.076 / 0.002 / -0.008 / -0.018 | none |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 309 | no | study | -0.022 / 0.026 / 0.068 | -0.001 / 0.000 / 0.083 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Pseudomonas aeruginosa* | ceftolozane/tazobactam | 259 | no | lineage | 0.071 / -0.036 / -0.020 | 0.000 / -0.056 / -0.030 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.092 / 0.072 / 0.034 / -0.120 / 0.001 | none |
| *Pseudomonas aeruginosa* | gentamicin | 346 | no | lineage | -0.004 / -0.004 / 0.019 | -0.013 / 0.005 / 0.010 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.005 / 0.000 / -0.008 / 0.003 / 0.016 | none |
| *Pseudomonas aeruginosa* | gentamicin | 346 | no | study | 0.000 / 0.002 / 0.008 | 0.007 / -0.003 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Pseudomonas aeruginosa* | imipenem | 447 | no | lineage | 0.027 / 0.005 / 0.018 | 0.027 / 0.006 / 0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.005 / 0.003 / -0.005 / 0.019 / -0.008 | none |
| *Pseudomonas aeruginosa* | imipenem | 447 | no | study | 0.032 / 0.014 / -0.016 | 0.026 / 0.040 / 0.029 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Pseudomonas aeruginosa* | levofloxacin | 702 | no | lineage | 0.004 / 0.004 / 0.009 | -0.009 / 0.015 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / -0.001 / -0.004 / 0.002 / -0.001 | none |
| *Pseudomonas aeruginosa* | meropenem | 883 | yes | lineage | 0.011 / -0.006 / 0.009 | 0.003 / -0.011 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.004 / -0.010 / -0.008 / -0.024 / -0.002 | none |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 210 | no | lineage | -0.095 / -0.003 / -0.079 | -0.057 / -0.102 / -0.038 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.010 / -0.003 / 0.010 / -0.076 / 0.011 | none |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 210 | no | study | -0.050 / 0.008 / -0.046 | -0.036 / -0.016 / -0.077 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Pseudomonas aeruginosa* | tobramycin | 645 | no | lineage | -0.010 / -0.002 / 0.031 | -0.007 / -0.017 / 0.030 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.001 / -0.002 / -0.012 / -0.003 / -0.018 | none |
| *Salmonella enterica* | amoxicillin/clavulanate | 799 | no | lineage | 0.000 / -0.004 / -0.003 | -0.002 / -0.001 / 0.008 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.003 / 0.000 / 0.009 | none |
| *Salmonella enterica* | ampicillin | 1,997 | yes | lineage | -0.013 / -0.046 / 0.005 | 0.012 / -0.005 / 0.034 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.006 / 0.036 / -0.027 / -0.001 / -0.027 | none |
| *Salmonella enterica* | ampicillin/sulbactam | 183 | no | lineage | -0.005 / 0.025 / 0.009 | 0.003 / 0.009 / 0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.020 / -0.010 / 0.001 / 0.016 / 0.016 | none |
| *Salmonella enterica* | azithromycin | 1,355 | yes | lineage | -0.020 / 0.014 / 0.000 | -0.132 / 0.022 / 0.007 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.010 / -0.019 / -0.001 / 0.014 / 0.001 | none |
| *Salmonella enterica* | cefazolin | 185 | no | lineage | 0.028 / 0.044 / 0.017 | 0.006 / 0.077 / 0.026 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.011 / -0.027 / 0.005 / 0.025 / -0.019 | none |
| *Salmonella enterica* | cefotaxime | 941 | yes | lineage | 0.052 / 0.002 / -0.096 | 0.084 / 0.088 / 0.038 | 0.000 / 0.001 / 0.000 / 0.000 / 0.000 | -0.067 / -0.061 / 0.013 / 0.055 / -0.010 | none |
| *Salmonella enterica* | cefoxitin | 1,080 | yes | lineage | 0.000 / -0.002 / 0.003 | -0.001 / 0.000 / -0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / -0.003 | none |
| *Salmonella enterica* | ceftazidime | 783 | no | lineage | -0.072 / 0.054 / 0.001 | -0.079 / 0.037 / 0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.004 / 0.053 / -0.031 / 0.008 / 0.001 | none |
| *Salmonella enterica* | ceftiofur | 689 | no | lineage | 0.026 / -0.003 / 0.022 | 0.018 / 0.014 / 0.026 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.005 / 0.000 / 0.000 / -0.005 / -0.008 | none |
| *Salmonella enterica* | ceftriaxone | 1,155 | yes | lineage | 0.032 / 0.020 / 0.041 | 0.034 / 0.033 / 0.040 | 0.003 / 0.001 / 0.000 / 0.000 / 0.000 | -0.010 / 0.035 / 0.002 / -0.020 / -0.005 | none |
| *Salmonella enterica* | chloramphenicol | 1,989 | yes | lineage | -0.028 / 0.079 / 0.032 | 0.008 / 0.084 / 0.054 | 0.000 / 0.005 / 0.000 / 0.000 / 0.000 | -0.020 / 0.011 / -0.020 / 0.009 / -0.010 | none |
| *Salmonella enterica* | ciprofloxacin | 1,503 | yes | lineage | -0.034 / 0.017 / 0.025 | -0.016 / 0.012 / 0.032 | -0.006 / 0.000 / 0.000 / 0.000 / 0.000 | -0.015 / -0.003 / -0.014 / 0.005 / 0.014 | none |
| *Salmonella enterica* | gentamicin | 1,675 | yes | lineage | 0.036 / 0.011 / 0.033 | 0.048 / 0.048 / 0.079 | 0.003 / 0.006 / 0.000 / 0.000 / 0.000 | 0.011 / -0.002 / -0.044 / -0.025 / 0.000 | none |
| *Salmonella enterica* | kanamycin | 584 | no | lineage | 0.043 / 0.028 / 0.054 | -0.011 / 0.018 / 0.051 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.002 / -0.004 / -0.013 / -0.031 / 0.033 | none |
| *Salmonella enterica* | nalidixic acid | 1,599 | yes | lineage | -0.009 / -0.011 / 0.017 | -0.027 / 0.015 / 0.055 | -0.002 / 0.000 / 0.000 / 0.000 / 0.000 | -0.015 / 0.000 / 0.003 / 0.004 / 0.011 | none |
| *Salmonella enterica* | pefloxacin | 348 | no | lineage | -0.028 / 0.000 / 0.000 | -0.072 / -0.003 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.002 / 0.000 / 0.051 / 0.000 / 0.000 | none |
| *Salmonella enterica* | spectinomycin | 396 | no | lineage | 0.182 / 0.009 / -0.019 | 0.090 / 0.031 / -0.037 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.013 / 0.004 / -0.011 / -0.006 / 0.011 | none |
| *Salmonella enterica* | streptomycin | 1,269 | yes | lineage | 0.002 / 0.029 / 0.012 | 0.006 / 0.037 / 0.015 | -0.005 / -0.005 / 0.000 / 0.000 / 0.000 | 0.013 / -0.035 / -0.002 / 0.033 / 0.028 | none |
| *Salmonella enterica* | sulfamethoxazole | 275 | no | lineage | -0.007 / 0.135 / 0.005 | -0.071 / 0.053 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / -0.045 / 0.000 | none |
| *Salmonella enterica* | sulfamethoxazole/trimethoprim | 1,145 | yes | lineage | -0.016 / -0.054 / 0.009 | -0.006 / -0.028 / 0.028 | 0.001 / -0.010 / 0.000 / 0.000 / 0.000 | -0.033 / 0.025 / 0.005 / -0.008 / 0.038 | none |
| *Salmonella enterica* | sulfathiazole | 225 | no | lineage | 0.032 / 0.111 / 0.072 | 0.019 / 0.156 / 0.065 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.015 / -0.019 / 0.009 / -0.095 / -0.074 | none |
| *Salmonella enterica* | sulfisoxazole | 544 | no | lineage | 0.106 / 0.106 / 0.148 | 0.112 / 0.120 / 0.148 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.031 / -0.021 / 0.005 / -0.018 / 0.012 | none |
| *Salmonella enterica* | sulfisoxazole | 544 | no | study | 0.073 / 0.076 / 0.074 | 0.066 / 0.078 / 0.073 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Salmonella enterica* | tetracycline | 1,869 | yes | lineage | 0.122 / 0.080 / 0.050 | 0.127 / 0.092 / 0.072 | 0.000 / -0.003 / 0.000 / 0.000 / 0.000 | 0.022 / 0.001 / -0.054 / -0.008 / -0.009 | none |
| *Salmonella enterica* | trimethoprim | 634 | no | lineage | 0.004 / -0.026 / -0.020 | -0.005 / 0.018 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.016 / 0.056 / 0.003 / 0.046 / -0.011 | none |
| *Shigella flexneri* | ampicillin | 358 | no | lineage | -0.098 / 0.016 / -0.008 | -0.017 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.010 / -0.001 / -0.004 / 0.015 / -0.001 | none |
| *Shigella flexneri* | azithromycin | 361 | no | lineage | 0.009 / 0.021 / 0.023 | 0.013 / 0.021 / 0.019 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.002 / 0.013 / 0.006 / -0.003 / 0.007 | none |
| *Shigella flexneri* | ciprofloxacin | 361 | no | lineage | -0.007 / -0.016 / -0.017 | -0.008 / -0.001 / 0.005 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / 0.010 / -0.009 / -0.003 / 0.040 | none |
| *Shigella flexneri* | trimethoprim | 362 | no | lineage | -0.089 / 0.010 / -0.056 | -0.115 / 0.082 / 0.043 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.009 / -0.006 / -0.004 / 0.033 / -0.009 | none |
| *Shigella sonnei* | amikacin | 580 | no | lineage | 0.041 / 0.025 / 0.023 | 0.040 / 0.019 / 0.016 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.011 / -0.005 / -0.009 / 0.002 / -0.009 | none |
| *Shigella sonnei* | ampicillin | 1,674 | yes | lineage | 0.006 / -0.003 / 0.001 | 0.010 / 0.005 / 0.012 | 0.001 / 0.009 / 0.000 / 0.000 / 0.000 | -0.001 / 0.010 / 0.005 / 0.005 / 0.007 | none |
| *Shigella sonnei* | azithromycin | 1,307 | yes | lineage | -0.003 / 0.004 / 0.000 | 0.001 / 0.003 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.001 / -0.003 | none |
| *Shigella sonnei* | cefoperazone | 257 | no | lineage | 0.006 / -0.169 / 0.000 | -0.016 / -0.143 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / 0.003 / -0.005 / -0.093 / -0.006 | none |
| *Shigella sonnei* | cefotaxime | 593 | no | lineage | -0.032 / -0.022 / -0.023 | -0.043 / -0.017 / -0.035 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.014 / -0.002 / -0.004 / -0.003 / -0.016 | none |
| *Shigella sonnei* | cefoxitin | 593 | no | lineage | -0.038 / -0.035 / -0.007 | -0.047 / -0.013 / -0.009 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.040 / -0.018 / -0.020 / 0.005 / -0.001 | none |
| *Shigella sonnei* | cefpirome | 333 | no | lineage | -0.006 / -0.027 / -0.010 | -0.021 / -0.001 / 0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.019 / 0.018 / -0.005 / -0.012 / 0.006 | none |
| *Shigella sonnei* | ceftazidime | 582 | no | lineage | -0.033 / -0.008 / -0.007 | -0.022 / -0.003 / -0.021 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.017 / -0.014 / -0.012 / -0.011 / -0.010 | none |
| *Shigella sonnei* | ceftriaxone | 1,163 | yes | lineage | 0.003 / 0.006 / 0.000 | 0.000 / 0.005 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.001 / -0.005 / 0.000 | none |
| *Shigella sonnei* | ciprofloxacin | 1,397 | yes | lineage | 0.001 / 0.000 / 0.002 | 0.002 / 0.006 / 0.009 | 0.000 / 0.007 / 0.000 / 0.000 / 0.000 | 0.005 / 0.006 / 0.002 / 0.003 / 0.003 | none |
| *Shigella sonnei* | colistin | 268 | no | lineage | -0.025 / 0.000 / 0.000 | -0.001 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.009 / -0.010 / 0.000 | none |
| *Shigella sonnei* | ertapenem | 338 | no | lineage | -0.022 / -0.006 / -0.033 | -0.028 / -0.011 / -0.017 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.006 / 0.008 / 0.000 / -0.006 / -0.006 | none |
| *Shigella sonnei* | gentamicin | 1,288 | yes | lineage | 0.000 / -0.019 / -0.034 | 0.034 / 0.016 / 0.012 | 0.000 / 0.002 / 0.000 / 0.000 / 0.000 | 0.010 / -0.004 / 0.001 / 0.004 / -0.003 | none |
| *Shigella sonnei* | nalidixic acid | 854 | yes | lineage | 0.010 / 0.011 / 0.006 | 0.011 / 0.020 / 0.023 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.001 / -0.002 / 0.002 / 0.002 | none |
| *Shigella sonnei* | streptomycin | 595 | no | lineage | 0.002 / 0.011 / 0.011 | -0.008 / 0.012 / 0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.003 / 0.015 / -0.014 / 0.003 / -0.002 | none |
| *Shigella sonnei* | sulfamethoxazole/trimethoprim | 617 | no | lineage | -0.010 / -0.019 / 0.013 | -0.036 / -0.020 / 0.011 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.003 / 0.000 / -0.020 / -0.006 / 0.021 | none |
| *Shigella sonnei* | temocillin | 333 | no | lineage | -0.006 / -0.027 / -0.010 | -0.021 / -0.001 / 0.012 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.019 / 0.018 / -0.005 / -0.012 / 0.006 | none |
| *Shigella sonnei* | tetracycline | 878 | yes | lineage | -0.015 / 0.003 / -0.015 | -0.021 / 0.000 / 0.005 | -0.001 / 0.000 / 0.000 / 0.000 / 0.000 | 0.005 / 0.006 / 0.001 / 0.008 / 0.006 | none |
| *Shigella sonnei* | tobramycin | 541 | no | lineage | 0.023 / -0.035 / -0.033 | 0.028 / 0.028 / 0.051 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.052 / -0.015 / -0.034 / -0.060 / -0.014 | none |
| *Shigella sonnei* | trimethoprim | 1,104 | yes | lineage | 0.098 / 0.008 / -0.004 | 0.107 / 0.010 / -0.025 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.004 / -0.003 / 0.030 / -0.005 / -0.019 | none |
| *Staphylococcus aureus* | cefoxitin | 1,407 | yes | lineage | 0.001 / 0.022 / -0.010 | 0.002 / 0.029 / 0.018 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.005 / -0.002 / 0.003 / 0.005 | none |
| *Staphylococcus aureus* | chloramphenicol | 1,038 | yes | lineage | -0.010 / -0.054 / -0.042 | -0.010 / -0.038 / 0.014 | -0.017 / 0.000 / 0.000 / 0.000 / 0.000 | -0.016 / 0.003 / 0.018 / -0.037 / 0.014 | none |
| *Staphylococcus aureus* | ciprofloxacin | 1,624 | yes | lineage | -0.001 / 0.006 / 0.000 | -0.009 / -0.013 / 0.001 | 0.000 / 0.002 / 0.000 / 0.000 / 0.000 | -0.006 / 0.000 / 0.000 / 0.009 / 0.002 | none |
| *Staphylococcus aureus* | clindamycin | 1,543 | yes | lineage | 0.034 / 0.161 / 0.259 | 0.008 / 0.180 / 0.206 | -0.001 / -0.010 / 0.000 / 0.000 / 0.000 | 0.010 / 0.005 / 0.003 / -0.035 / 0.000 | none |
| *Staphylococcus aureus* | clindamycin | 1,543 | yes | study | 0.018 / 0.023 / 0.020 | 0.021 / 0.034 / 0.015 | 0.000 / -0.006 / 0.000 / 0.000 / 0.000 |  | none |
| *Staphylococcus aureus* | daptomycin | 904 | yes | lineage | 0.000 / 0.042 / -0.001 | 0.000 / -0.004 / -0.001 | -0.021 / 0.002 / 0.000 / -0.011 / 0.007 | 0.001 / 0.000 / 0.000 / 0.043 / 0.000 | none |
| *Staphylococcus aureus* | erythromycin | 1,957 | yes | lineage | 0.191 / 0.056 / 0.197 | 0.184 / 0.135 / 0.308 | -0.001 / 0.003 / 0.000 / 0.000 / 0.000 | -0.032 / -0.006 / -0.004 / -0.025 / 0.003 | none |
| *Staphylococcus aureus* | erythromycin | 1,957 | yes | study | -0.031 / -0.006 / -0.009 | -0.024 / -0.001 / 0.013 | 0.024 / 0.031 / 0.000 / 0.000 / 0.000 |  | none |
| *Staphylococcus aureus* | fusidic acid | 1,101 | yes | lineage | -0.020 / 0.006 / -0.001 | -0.083 / -0.002 / -0.010 | 0.004 / -0.011 / 0.000 / 0.000 / 0.000 | -0.005 / -0.022 / 0.019 / -0.023 / -0.004 | none |
| *Staphylococcus aureus* | gentamicin | 1,913 | yes | lineage | 0.002 / 0.074 / 0.010 | -0.002 / 0.019 / 0.027 | 0.021 / 0.009 / 0.000 / 0.000 / 0.000 | 0.001 / 0.000 / 0.000 / 0.016 / -0.009 | none |
| *Staphylococcus aureus* | mupirocin | 736 | no | lineage | 0.049 / 0.000 / 0.000 | -0.093 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.033 / 0.000 / 0.000 | none |
| *Staphylococcus aureus* | penicillin | 1,508 | yes | lineage | 0.000 / 0.003 / -0.002 | 0.002 / 0.002 / -0.008 | -0.001 / 0.009 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.003 / 0.001 / 0.002 | none |
| *Staphylococcus aureus* | penicillin | 1,508 | yes | study | 0.006 / -0.002 / 0.007 | 0.004 / 0.006 / -0.003 | 0.000 / 0.011 / -0.004 / 0.008 / 0.000 |  | none |
| *Staphylococcus aureus* | rifampin | 931 | yes | lineage | 0.020 / -0.013 / 0.000 | -0.040 / -0.003 / 0.000 | 0.000 / 0.000 / 0.031 / -0.007 / 0.000 | 0.000 / 0.000 / 0.024 / 0.000 / 0.000 | none |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | 1,389 | yes | lineage | -0.367 / -0.113 / 0.004 | -0.386 / -0.032 / 0.004 | 0.000 / 0.000 / 0.000 / -0.001 / 0.000 | 0.000 / -0.002 / -0.014 / -0.017 / 0.001 | none |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | 1,389 | yes | study | -0.017 / 0.009 / 0.001 | -0.020 / 0.022 / 0.004 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 |  | none |
| *Staphylococcus aureus* | tetracycline | 1,647 | yes | lineage | 0.066 / 0.021 / 0.078 | 0.067 / -0.002 / 0.079 | -0.006 / 0.005 / 0.000 / 0.000 / 0.000 | 0.005 / 0.013 / -0.001 / 0.004 / 0.000 | none |
| *Staphylococcus haemolyticus* | ciprofloxacin | 124 | no | lineage | -0.030 / -0.080 / -0.079 | -0.061 / -0.001 / -0.013 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.002 / -0.002 / -0.018 / -0.029 / -0.033 | none |
| *Staphylococcus haemolyticus* | fusidic acid | 114 | no | lineage | 0.006 / -0.058 / -0.010 | -0.019 / 0.039 / 0.043 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.009 / 0.019 / -0.026 / 0.049 / 0.004 | none |
| *Staphylococcus haemolyticus* | tetracycline | 104 | no | lineage | -0.012 / 0.025 / 0.000 | -0.036 / 0.004 / 0.009 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.042 / -0.009 / -0.001 / 0.040 / 0.041 | none |
| *Streptococcus agalactiae* | clindamycin | 141 | no | lineage | 0.000 / 0.024 / 0.017 | -0.001 / 0.022 / 0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.006 / -0.001 / 0.000 / 0.004 / 0.008 | none |
| *Streptococcus agalactiae* | erythromycin | 141 | no | lineage | 0.001 / 0.025 / 0.004 | 0.000 / 0.011 / 0.003 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.048 / 0.005 / -0.015 / -0.014 / -0.008 | none |
| *Streptococcus agalactiae* | tetracycline | 141 | no | lineage | -0.014 / -0.016 / -0.018 | -0.005 / 0.002 / 0.002 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | -0.006 / -0.011 / -0.003 / -0.003 / 0.005 | none |
| *Streptococcus pneumoniae* | azithromycin | 1,297 | yes | lineage | 0.013 / -0.006 / 0.017 | 0.012 / 0.024 / 0.009 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.001 / 0.033 / 0.017 | none |
| *Streptococcus pneumoniae* | cefuroxime | 1,732 | yes | lineage | -0.001 / -0.003 / 0.003 | -0.002 / -0.007 / -0.001 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.014 / -0.007 / 0.000 | none |
| *Streptococcus pneumoniae* | chloramphenicol | 1,859 | yes | lineage | -0.021 / 0.024 / 0.041 | -0.134 / 0.035 / 0.061 | -0.019 / -0.017 / 0.000 / 0.000 / 0.000 | -0.025 / 0.013 / -0.002 / 0.018 / 0.001 | none |
| *Streptococcus pneumoniae* | clindamycin | 1,894 | yes | lineage | 0.000 / 0.001 / 0.016 | -0.014 / 0.001 / 0.024 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.001 / 0.000 / 0.016 | none |
| *Streptococcus pneumoniae* | erythromycin | 1,991 | yes | lineage | -0.003 / -0.001 / -0.025 | -0.002 / -0.003 / -0.049 | 0.000 / -0.011 / 0.000 / 0.000 / 0.000 | -0.004 / 0.025 / 0.000 / 0.001 / 0.005 | none |
| *Streptococcus pneumoniae* | meropenem | 1,808 | yes | lineage | -0.044 / -0.140 / -0.038 | -0.166 / 0.074 / 0.000 | -0.002 / 0.000 / 0.000 / -0.040 / 0.000 | -0.032 / -0.001 / 0.030 / 0.046 / 0.021 | none |
| *Streptococcus pneumoniae* | sulfamethoxazole/trimethoprim | 1,777 | yes | lineage | 0.008 / -0.008 / -0.012 | 0.002 / -0.007 / -0.010 | 0.000 / -0.004 / 0.000 / 0.000 / 0.000 | -0.001 / -0.001 / 0.001 / -0.022 / 0.007 | none |
| *Streptococcus pneumoniae* | tetracycline | 1,971 | yes | lineage | 0.001 / -0.011 / -0.003 | 0.001 / 0.011 / -0.020 | 0.001 / 0.000 / 0.000 / 0.000 / 0.000 | -0.007 / 0.001 / 0.001 / 0.008 / -0.004 | none |
| *Vibrio cholerae* | nalidixic acid | 91 | no | lineage | 0.003 / -0.006 / -0.039 | 0.007 / -0.011 / -0.031 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.003 / 0.003 / 0.000 / 0.003 / 0.007 | none |
| *Vibrio cholerae* | sulfamethoxazole | 74 | no | lineage | 0.000 / 0.032 / 0.000 | 0.000 / 0.032 / 0.000 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.023 / 0.032 / -0.063 | none |
| *Vibrio cholerae* | trimethoprim | 74 | no | lineage | 0.034 / 0.010 / 0.025 | 0.018 / 0.015 / 0.025 | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.020 / 0.020 / -0.005 / -0.015 / 0.044 | none |
