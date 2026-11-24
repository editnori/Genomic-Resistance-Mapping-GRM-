# Settings-arm result

The release fits each learner with its own settings, some of them chosen inside each fit, and calls a genome resistant when its score is at least 0.5. This page asks whether a search of each learner's settings, made inside each outer training fold, or a call cutoff chosen on inner scores, serves better on the outer test folds. It reads block C5 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/settings-arms.md`.

The design was fixed at 04:08Z on 23 Nov 2026 and the lead's rulings on its questions at 05:52Z, before any C5 score was read (`docs/assumption-plan-corrections.md`, entry "C5 design"). The table script, its gates and this page's form were committed before any score of the run was opened.

The scripts that read the run were fixed at 00:58:43Z on 24 Nov 2026, before any score of the run was opened. The count was read once, from 01:05:26Z to 01:06:07Z on 24 Nov 2026, with no other option, subset or measure.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

"No difference shown" is the rule's reading and no statement that two rules are equal. Recall stands beside log loss wherever log loss is given. The cutoff arm moves no score, so its log loss equals the default's by construction.

## Reading

**Result.** By ruling 15, neither arm meets conditions 1 and 2, so the rule reads no difference shown for the search and for the cutoff arm. On the lineage split 54 of the 55 datasets in 8 species enter the count, and 0 of 5 learners carry either arm: no learner's median difference in balanced accuracy has an interval wholly above 0, and better beats worse after Holm on none of the 15 rows of the family. Each of the 15 Holm p values is 1.000. The twin meets condition 1 in 0 of 5 learners and condition 2 in 0 of 5, so condition 3 holds. No learner's recall interval lies wholly below -0.02 for either arm, so condition 4 holds for all five. Condition 5 is not met for either arm, as no learner carries. Condition 6 is the maintainer's.

**Search.** On lineage balanced accuracy, better and worse of 54 datasets: `scm` 4 and 9, `cart` 7 and 4, `rf` 8 and 5, `gbt` 5 and 7, `lasso` 0 and 0. Each learner's median difference is 0.000. Recall: each learner's median difference is 0.000. Log loss: higher than the default's for `cart` on 31 datasets and lower on 5 (median 0.049, 0.021 to 0.074), and for `gbt` higher on 11 and lower on 2 (median 0.004, 0.003 to 0.008).

**Cutoff arm.** On lineage balanced accuracy, better and worse of 54 datasets: `scm` 0 and 0, `cart` 4 and 12, `rf` 15 and 9, `gbt` 8 and 12, `lasso` 7 and 11. The median differences run from -0.006 (`gbt`) to 0.000; for `lasso` the median is -0.004 and its interval, -0.013 to -0.001, lies wholly below 0. Recall: the median difference is 0.024 for `rf` (0.006 to 0.104), 0.018 for `gbt` (0.011 to 0.109) and 0.009 for `lasso` (0.000 to 0.068), and 0.000 for `scm` and `cart`. Log loss: the cutoff arm's equals the default's by construction.

**The stricter reading.** Neither arm meets its conditions 1 and 2. It stands beside and decides nothing.

| Arm | Reading | Decides | Condition | Reads | Met | Learners |
|---|---|---|---|---|---|---|
| search | ruling 15 | yes | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 15; at least 3 of 5 carry | not met | 0 |
| search | ruling 15 | yes | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| search | ruling 15 | yes | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| search | ruling 15 | yes | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 |
| search | ruling 15 | yes | 6 | the maintainer rules the change | the maintainer's | none |
| search | stricter | no, beside | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 0 |
| search | stricter | no, beside | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 15, for all five | not met | 0 |
| search | stricter | no, beside | 3 | the twin meets neither 1 nor 2 as this reading reads them | met | none |
| search | stricter | no, beside | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| search | stricter | no, beside | 5 | each of the five learners has a study median of its lineage median's sign | not met | 0 |
| search | stricter | no, beside | 6 | the maintainer rules the change | the maintainer's | none |
| cutoff | ruling 15 | yes | 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 15; at least 3 of 5 carry | not met | 0 |
| cutoff | ruling 15 | yes | 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart (learners: the larger of the twin's two counts) | met | 0 |
| cutoff | ruling 15 | yes | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| cutoff | ruling 15 | yes | 5 | each learner that carries has a study median of its lineage median's sign | not met | 0 |
| cutoff | ruling 15 | yes | 6 | the maintainer rules the change | the maintainer's | none |
| cutoff | stricter | no, beside | 1 | the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners | not met | 0 |
| cutoff | stricter | no, beside | 2 | better beats worse by the one-sided sign test at p < 0.05 after Holm over 15, for all five | not met | 0 |
| cutoff | stricter | no, beside | 3 | the twin meets neither 1 nor 2 as this reading reads them | met | none |
| cutoff | stricter | no, beside | 4 | no learner's recall interval lies wholly below -0.02, for all five | met | 5 |
| cutoff | stricter | no, beside | 5 | each of the five learners has a study median of its lineage median's sign | not met | 3 |
| cutoff | stricter | no, beside | 6 | the maintainer rules the change | the maintainer's | none |

**Arms.** Two arms are read against the default, the release's rule refitted in the same run with the same binary. The search picks, per learner and outer fold, one of 31 trials (trial 0 is the release's own setting) by the default-anchored rule at 0.04, from the mean balanced accuracy of 3 inner folds of the outer fold's training genomes. The cutoff arm keeps the release's fit and calls at the point of a 19-point grid, 0.05 to 0.95, with the best inner balanced accuracy.

**Conditions.** The conditions are read as ruling 15 reads them, with condition 3 as corrected: a learner carries a change when it meets conditions 1 and 2, and at least 3 of 5 learners must carry it; condition 3 holds when the twin meets 1 in fewer than 3 learners and 2 in fewer than 3, each read apart; condition 4 holds when no learner's recall interval lies wholly below -0.02; condition 5 holds when each learner that carries has a study median of its lineage median's sign. Condition 6 is the maintainer's. The stricter reading stands beside ruling 15 and decides nothing.

**The cutoff arm's log loss.** The cutoff arm moves no score, so its log loss equals the default's by construction on every dataset. Its log loss rows count every dataset as the same and decide nothing. `scm` scores 0 or 1, so its cutoff arm is the default on every fold.

The release keeps its settings and its 0.5 cutoff whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The search covers 55 datasets in 8 species, those in both S and L and not in M. The 14 datasets of M wait for the label rule. Each of the 55 has a reviewed study plan and a study split. The section Datasets names them all.

**Draws.** The default and both arms come from the same runs at fold seeds 42, 43 and 44 on the lineage split. A dataset's value is the mean of its draws. The twin is the default at seeds 45, 46 and 47, on the lineage split alone.

**Study split.** Each of the 55 datasets gets one study run per arm at the release's seeds, with no twin, as blocks C3 and C9 did.

**Plain best.** The plain best, no anchor, is the trial with the highest inner mean, a tie going to the lowest trial id, refitted on the outer training genomes and scored beside the two arms. It stands in no Holm family and no condition, and decides nothing.

**Not scored.** A dataset with a draw not scored is not scored, and a missing draw is never read as zero. *Klebsiella pneumoniae* with polymyxin B leaves the lineage count, as it left block C3's (`docs/input-arms-result.md`, Not scored). Its default draw on the lineage split at seed 43 came back unavailable, with no outer fold, and the result gives the reason as `invalid lineage outer evaluation: fold <k> validation lacks one or both classes`. The twin's draw at seed 46 is unavailable for the same reason. Each arm's row at seed 43 is not scored, for the default's reason, as C3's scorer marks the default's own row, and no choices file is owed there. Its draws at seeds 42 and 44, the twin's at 45 and 47 and its study run each hold 5 outer folds.

The rescored `checks.tsv` gives the reason with `<k>` as 4: `invalid lineage outer evaluation: fold 4 validation lacks one or both classes`, as block C3's results gave it.

**Entered.** On the lineage split, 54 of the 55 datasets enter the count of the default, both arms, the twin and the plain best; the one left out is *Klebsiella pneumoniae* with polymyxin B, with 2 of its 3 draws scored for the default, each arm and the twin. On the study split, 55 of 55 enter.

**Results.** `expect.py` owes 1,045 results, and the copies hold them: 385 default and twin results and 1,095 choices files, one per outer fold of each draw and study run of the default whose arms are scored, 1,480 files with no other. The scorer checked every result, and every check passed (Method). No draw is missing. No study run was refused, and no draw is a constant model.

## Count

The primary measure is each learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. Every difference is the arm's value less the default's; for the very major error rate, the major error rate and log loss, lower is better, so better and worse turn round. The interval is a 95% interval of the median over 4,000 resamples of species. The sign test is one-sided, and Holm runs over the 15 lineage rows of the primary measure: the search, the cutoff arm and the twin, each by five learners.

The family: lineage balanced accuracy, the 15 rows Holm runs over, each of the 54 datasets entered.

| Pair | Learner | Better | Worse | Same | Median | Interval | Sign p | Holm p | Condition 1 | Condition 2 |
|---|---|---|---|---|---|---|---|---|---|---|
| search less default | `scm` | 4 of 54 | 9 | 41 | 0.000 | -0.001 to 0.000 | 0.954 | 1.000 | not met | not met |
| search less default | `cart` | 7 of 54 | 4 | 43 | 0.000 | -0.001 to 0.000 | 0.274 | 1.000 | not met | not met |
| search less default | `rf` | 8 of 54 | 5 | 41 | 0.000 | -0.002 to 0.002 | 0.291 | 1.000 | not met | not met |
| search less default | `gbt` | 5 of 54 | 7 | 42 | 0.000 | -0.001 to 0.002 | 0.806 | 1.000 | not met | not met |
| search less default | `lasso` | 0 of 54 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 | 1.000 | not met | not met |
| cutoff less default | `scm` | 0 of 54 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 | 1.000 | not met | not met |
| cutoff less default | `cart` | 4 of 54 | 12 | 38 | 0.000 | -0.006 to 0.000 | 0.989 | 1.000 | not met | not met |
| cutoff less default | `rf` | 15 of 54 | 9 | 30 | 0.000 | -0.018 to 0.013 | 0.154 | 1.000 | not met | not met |
| cutoff less default | `gbt` | 8 of 54 | 12 | 34 | -0.006 | -0.012 to 0.000 | 0.868 | 1.000 | not met | not met |
| cutoff less default | `lasso` | 7 of 54 | 11 | 36 | -0.004 | -0.013 to -0.001 | 0.881 | 1.000 | not met | not met |
| twin less default | `scm` | 7 of 54 | 8 | 39 | 0.002 | 0.000 to 0.005 | 0.696 | 1.000 | not met | not met |
| twin less default | `cart` | 6 of 54 | 9 | 39 | -0.002 | -0.006 to 0.000 | 0.849 | 1.000 | not met | not met |
| twin less default | `rf` | 4 of 54 | 4 | 46 | 0.000 | -0.002 to 0.000 | 0.637 | 1.000 | not met | not met |
| twin less default | `gbt` | 8 of 54 | 7 | 39 | 0.000 | -0.007 to 0.009 | 0.500 | 1.000 | not met | not met |
| twin less default | `lasso` | 6 of 54 | 2 | 46 | 0.000 | -0.002 to 0.010 | 0.145 | 1.000 | not met | not met |

## Per arm

Each measure's median difference, with better and worse in brackets, per split, arm and learner:

| Split | Pair | Learner | Balanced accuracy | Recall | Specificity | Very major error rate | Major error rate | Log loss |
|---|---|---|---|---|---|---|---|---|
| lineage | search less default | `scm` | 0.000 (4 / 9) | 0.000 (10 / 4) | -0.001 (6 / 12) | 0.000 (10 / 4) | 0.001 (6 / 12) | 0.000 (12 / 22) |
| lineage | search less default | `cart` | 0.000 (7 / 4) | 0.000 (15 / 7) | -0.001 (5 / 13) | 0.000 (15 / 7) | 0.001 (5 / 13) | 0.049 (5 / 31) |
| lineage | search less default | `rf` | 0.000 (8 / 5) | 0.000 (10 / 9) | 0.000 (7 / 8) | 0.000 (10 / 9) | 0.000 (7 / 8) | 0.000 (5 / 10) |
| lineage | search less default | `gbt` | 0.000 (5 / 7) | 0.000 (10 / 12) | 0.000 (10 / 10) | 0.000 (10 / 12) | 0.000 (10 / 10) | 0.004 (2 / 11) |
| lineage | search less default | `lasso` | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 1) | 0.000 (0 / 0) | 0.000 (0 / 1) | 0.000 (0 / 0) |
| lineage | cutoff less default | `scm` | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) |
| lineage | cutoff less default | `cart` | 0.000 (4 / 12) | 0.000 (12 / 5) | 0.000 (4 / 16) | 0.000 (12 / 5) | 0.000 (4 / 16) | 0.000 (0 / 0) |
| lineage | cutoff less default | `rf` | 0.000 (15 / 9) | 0.024 (28 / 14) | -0.028 (10 / 30) | -0.024 (28 / 14) | 0.028 (10 / 30) | 0.000 (0 / 0) |
| lineage | cutoff less default | `gbt` | -0.006 (8 / 12) | 0.018 (27 / 8) | -0.037 (5 / 33) | -0.018 (27 / 8) | 0.037 (5 / 33) | 0.000 (0 / 0) |
| lineage | cutoff less default | `lasso` | -0.004 (7 / 11) | 0.009 (23 / 9) | -0.025 (5 / 29) | -0.009 (23 / 9) | 0.025 (5 / 29) | 0.000 (0 / 0) |
| lineage | twin less default | `scm` | 0.002 (7 / 8) | 0.003 (14 / 8) | 0.000 (10 / 12) | -0.003 (14 / 8) | 0.000 (10 / 12) | -0.005 (22 / 19) |
| lineage | twin less default | `cart` | -0.002 (6 / 9) | 0.000 (10 / 13) | -0.003 (9 / 12) | 0.000 (10 / 13) | 0.003 (9 / 12) | -0.002 (20 / 16) |
| lineage | twin less default | `rf` | 0.000 (4 / 4) | -0.002 (6 / 9) | 0.000 (6 / 4) | 0.002 (6 / 9) | 0.000 (6 / 4) | -0.001 (2 / 3) |
| lineage | twin less default | `gbt` | 0.000 (8 / 7) | -0.002 (12 / 14) | -0.001 (10 / 12) | 0.002 (12 / 14) | 0.001 (10 / 12) | -0.003 (12 / 8) |
| lineage | twin less default | `lasso` | 0.000 (6 / 2) | 0.000 (14 / 5) | 0.000 (11 / 9) | 0.000 (14 / 5) | 0.000 (11 / 9) | 0.001 (15 / 13) |
| study | search less default | `scm` | 0.000 (8 / 10) | 0.000 (10 / 10) | 0.000 (10 / 16) | 0.000 (10 / 10) | 0.000 (10 / 16) | 0.000 (14 / 18) |
| study | search less default | `cart` | -0.002 (5 / 12) | 0.000 (9 / 8) | 0.000 (7 / 14) | 0.000 (9 / 8) | 0.000 (7 / 14) | 0.000 (17 / 24) |
| study | search less default | `rf` | 0.000 (12 / 4) | 0.000 (10 / 10) | 0.000 (13 / 12) | 0.000 (10 / 10) | 0.000 (13 / 12) | 0.000 (12 / 11) |
| study | search less default | `gbt` | 0.000 (6 / 10) | 0.000 (9 / 10) | 0.000 (5 / 14) | 0.000 (9 / 10) | 0.000 (5 / 14) | 0.000 (4 / 11) |
| study | search less default | `lasso` | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (1 / 0) | 0.000 (0 / 0) | 0.000 (1 / 0) | 0.000 (0 / 0) |
| study | cutoff less default | `scm` | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) | 0.000 (0 / 0) |
| study | cutoff less default | `cart` | 0.000 (1 / 11) | 0.000 (7 / 6) | 0.000 (7 / 11) | 0.000 (7 / 6) | 0.000 (7 / 11) | 0.000 (0 / 0) |
| study | cutoff less default | `rf` | 0.007 (20 / 6) | 0.023 (29 / 9) | -0.008 (16 / 22) | -0.023 (29 / 9) | 0.008 (16 / 22) | 0.000 (0 / 0) |
| study | cutoff less default | `gbt` | -0.007 (5 / 14) | 0.017 (25 / 9) | -0.024 (11 / 29) | -0.017 (25 / 9) | 0.024 (11 / 29) | 0.000 (0 / 0) |
| study | cutoff less default | `lasso` | -0.002 (8 / 12) | 0.013 (25 / 12) | -0.031 (10 / 32) | -0.013 (25 / 12) | 0.031 (10 / 32) | 0.000 (0 / 0) |

Every row of the count, with its interval and whether it decides:

| Split | Pair | Learner | Measure | Better is | Entered | Better | Worse | Same | Median | Interval | Decides |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lineage | search less default | `scm` | balanced accuracy | higher | 54 of 55 | 4 | 9 | 41 | 0.000 | -0.001 to 0.000 | yes |
| lineage | search less default | `scm` | recall | higher | 54 of 55 | 10 | 4 | 40 | 0.000 | 0.000 to 0.004 | yes |
| lineage | search less default | `scm` | log loss | lower | 54 of 55 | 12 | 22 | 20 | 0.000 | 0.000 to 0.026 | yes |
| lineage | search less default | `cart` | balanced accuracy | higher | 54 of 55 | 7 | 4 | 43 | 0.000 | -0.001 to 0.000 | yes |
| lineage | search less default | `cart` | recall | higher | 54 of 55 | 15 | 7 | 32 | 0.000 | 0.000 to 0.016 | yes |
| lineage | search less default | `cart` | log loss | lower | 54 of 55 | 5 | 31 | 18 | 0.049 | 0.021 to 0.074 | yes |
| lineage | search less default | `rf` | balanced accuracy | higher | 54 of 55 | 8 | 5 | 41 | 0.000 | -0.002 to 0.002 | yes |
| lineage | search less default | `rf` | recall | higher | 54 of 55 | 10 | 9 | 35 | 0.000 | -0.006 to 0.000 | yes |
| lineage | search less default | `rf` | log loss | lower | 54 of 55 | 5 | 10 | 39 | 0.000 | 0.000 to 0.003 | yes |
| lineage | search less default | `gbt` | balanced accuracy | higher | 54 of 55 | 5 | 7 | 42 | 0.000 | -0.001 to 0.002 | yes |
| lineage | search less default | `gbt` | recall | higher | 54 of 55 | 10 | 12 | 32 | 0.000 | -0.004 to 0.010 | yes |
| lineage | search less default | `gbt` | log loss | lower | 54 of 55 | 2 | 11 | 41 | 0.004 | 0.003 to 0.008 | yes |
| lineage | search less default | `lasso` | balanced accuracy | higher | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | yes |
| lineage | search less default | `lasso` | recall | higher | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | yes |
| lineage | search less default | `lasso` | log loss | lower | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | yes |
| lineage | cutoff less default | `scm` | balanced accuracy | higher | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | yes |
| lineage | cutoff less default | `scm` | recall | higher | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | yes |
| lineage | cutoff less default | `scm` | log loss | lower | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| lineage | cutoff less default | `cart` | balanced accuracy | higher | 54 of 55 | 4 | 12 | 38 | 0.000 | -0.006 to 0.000 | yes |
| lineage | cutoff less default | `cart` | recall | higher | 54 of 55 | 12 | 5 | 37 | 0.000 | 0.000 to 0.000 | yes |
| lineage | cutoff less default | `cart` | log loss | lower | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| lineage | cutoff less default | `rf` | balanced accuracy | higher | 54 of 55 | 15 | 9 | 30 | 0.000 | -0.018 to 0.013 | yes |
| lineage | cutoff less default | `rf` | recall | higher | 54 of 55 | 28 | 14 | 12 | 0.024 | 0.006 to 0.104 | yes |
| lineage | cutoff less default | `rf` | log loss | lower | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| lineage | cutoff less default | `gbt` | balanced accuracy | higher | 54 of 55 | 8 | 12 | 34 | -0.006 | -0.012 to 0.000 | yes |
| lineage | cutoff less default | `gbt` | recall | higher | 54 of 55 | 27 | 8 | 19 | 0.018 | 0.011 to 0.109 | yes |
| lineage | cutoff less default | `gbt` | log loss | lower | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| lineage | cutoff less default | `lasso` | balanced accuracy | higher | 54 of 55 | 7 | 11 | 36 | -0.004 | -0.013 to -0.001 | yes |
| lineage | cutoff less default | `lasso` | recall | higher | 54 of 55 | 23 | 9 | 22 | 0.009 | 0.000 to 0.068 | yes |
| lineage | cutoff less default | `lasso` | log loss | lower | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| lineage | twin less default | `scm` | balanced accuracy | higher | 54 of 55 | 7 | 8 | 39 | 0.002 | 0.000 to 0.005 | yes |
| lineage | twin less default | `scm` | recall | higher | 54 of 55 | 14 | 8 | 32 | 0.003 | 0.000 to 0.005 | yes |
| lineage | twin less default | `scm` | log loss | lower | 54 of 55 | 22 | 19 | 13 | -0.005 | -0.021 to 0.028 | yes |
| lineage | twin less default | `cart` | balanced accuracy | higher | 54 of 55 | 6 | 9 | 39 | -0.002 | -0.006 to 0.000 | yes |
| lineage | twin less default | `cart` | recall | higher | 54 of 55 | 10 | 13 | 31 | 0.000 | -0.004 to 0.000 | yes |
| lineage | twin less default | `cart` | log loss | lower | 54 of 55 | 20 | 16 | 18 | -0.002 | -0.019 to 0.005 | yes |
| lineage | twin less default | `rf` | balanced accuracy | higher | 54 of 55 | 4 | 4 | 46 | 0.000 | -0.002 to 0.000 | yes |
| lineage | twin less default | `rf` | recall | higher | 54 of 55 | 6 | 9 | 39 | -0.002 | -0.008 to -0.001 | yes |
| lineage | twin less default | `rf` | log loss | lower | 54 of 55 | 2 | 3 | 49 | -0.001 | -0.004 to 0.001 | yes |
| lineage | twin less default | `gbt` | balanced accuracy | higher | 54 of 55 | 8 | 7 | 39 | 0.000 | -0.007 to 0.009 | yes |
| lineage | twin less default | `gbt` | recall | higher | 54 of 55 | 12 | 14 | 28 | -0.002 | -0.008 to 0.014 | yes |
| lineage | twin less default | `gbt` | log loss | lower | 54 of 55 | 12 | 8 | 34 | -0.003 | -0.011 to 0.002 | yes |
| lineage | twin less default | `lasso` | balanced accuracy | higher | 54 of 55 | 6 | 2 | 46 | 0.000 | -0.002 to 0.010 | yes |
| lineage | twin less default | `lasso` | recall | higher | 54 of 55 | 14 | 5 | 35 | 0.000 | -0.004 to 0.016 | yes |
| lineage | twin less default | `lasso` | log loss | lower | 54 of 55 | 15 | 13 | 26 | 0.001 | -0.021 to 0.007 | yes |
| study | search less default | `scm` | balanced accuracy | higher | 55 of 55 | 8 | 10 | 37 | 0.000 | 0.000 to 0.000 | yes |
| study | search less default | `scm` | recall | higher | 55 of 55 | 10 | 10 | 35 | 0.000 | 0.000 to 0.005 | yes |
| study | search less default | `scm` | log loss | lower | 55 of 55 | 14 | 18 | 23 | 0.000 | 0.000 to 0.000 | yes |
| study | search less default | `cart` | balanced accuracy | higher | 55 of 55 | 5 | 12 | 38 | -0.002 | -0.004 to 0.000 | yes |
| study | search less default | `cart` | recall | higher | 55 of 55 | 9 | 8 | 38 | 0.000 | 0.000 to 0.000 | yes |
| study | search less default | `cart` | log loss | lower | 55 of 55 | 17 | 24 | 14 | 0.000 | 0.000 to 0.187 | yes |
| study | search less default | `rf` | balanced accuracy | higher | 55 of 55 | 12 | 4 | 39 | 0.000 | -0.004 to 0.014 | yes |
| study | search less default | `rf` | recall | higher | 55 of 55 | 10 | 10 | 35 | 0.000 | 0.000 to 0.000 | yes |
| study | search less default | `rf` | log loss | lower | 55 of 55 | 12 | 11 | 32 | 0.000 | -0.001 to 0.007 | yes |
| study | search less default | `gbt` | balanced accuracy | higher | 55 of 55 | 6 | 10 | 39 | 0.000 | -0.018 to 0.000 | yes |
| study | search less default | `gbt` | recall | higher | 55 of 55 | 9 | 10 | 36 | 0.000 | -0.012 to 0.000 | yes |
| study | search less default | `gbt` | log loss | lower | 55 of 55 | 4 | 11 | 40 | 0.000 | 0.000 to 0.008 | yes |
| study | search less default | `lasso` | balanced accuracy | higher | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | yes |
| study | search less default | `lasso` | recall | higher | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | yes |
| study | search less default | `lasso` | log loss | lower | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | yes |
| study | cutoff less default | `scm` | balanced accuracy | higher | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | yes |
| study | cutoff less default | `scm` | recall | higher | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | yes |
| study | cutoff less default | `scm` | log loss | lower | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| study | cutoff less default | `cart` | balanced accuracy | higher | 55 of 55 | 1 | 11 | 43 | 0.000 | 0.000 to 0.000 | yes |
| study | cutoff less default | `cart` | recall | higher | 55 of 55 | 7 | 6 | 42 | 0.000 | 0.000 to 0.001 | yes |
| study | cutoff less default | `cart` | log loss | lower | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| study | cutoff less default | `rf` | balanced accuracy | higher | 55 of 55 | 20 | 6 | 29 | 0.007 | -0.003 to 0.011 | yes |
| study | cutoff less default | `rf` | recall | higher | 55 of 55 | 29 | 9 | 17 | 0.023 | 0.005 to 0.052 | yes |
| study | cutoff less default | `rf` | log loss | lower | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| study | cutoff less default | `gbt` | balanced accuracy | higher | 55 of 55 | 5 | 14 | 36 | -0.007 | -0.013 to -0.004 | yes |
| study | cutoff less default | `gbt` | recall | higher | 55 of 55 | 25 | 9 | 21 | 0.017 | 0.000 to 0.020 | yes |
| study | cutoff less default | `gbt` | log loss | lower | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | no, equal by construction |
| study | cutoff less default | `lasso` | balanced accuracy | higher | 55 of 55 | 8 | 12 | 35 | -0.002 | -0.004 to 0.001 | yes |
| study | cutoff less default | `lasso` | recall | higher | 55 of 55 | 25 | 12 | 18 | 0.013 | 0.000 to 0.029 | yes |
| study | cutoff less default | `lasso` | log loss | lower | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | no, equal by construction |

**Search.** The search's rows decide with the cutoff arm's, as the conditions read them.

**Cutoff.** The cutoff arm's log loss rows are equal by construction and decide nothing. Its other rows decide.

**What the search chose.** Counts only: per learner, the folds that kept trial 0 of the folds searched, the folds where the plain best differs from the pick, the cutoffs chosen as a count per grid point, and the folds where trial 0 or every trial failed inside. A trial is named by its id alone, and no trial's settings are given.

Per split and learner, the outer folds the search ran on and what it chose:

| Split | Learner | Folds | Kept trial 0 | Picked another | Plain best trial 0 | Plain best is the pick | Plain best differs | Trial 0 failed | Every trial failed | Cutoff moved from 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| lineage | `scm` | 820 | 499 of 820 | 321 | 208 | 529 | 236 of 820 | 55 | 55 | 0 of 820 |
| lineage | `cart` | 820 | 470 of 820 | 350 | 134 | 484 | 281 of 820 | 55 | 55 | 371 of 820 |
| lineage | `rf` | 820 | 515 of 820 | 305 | 47 | 352 | 413 of 820 | 55 | 55 | 668 of 820 |
| lineage | `gbt` | 820 | 436 of 820 | 384 | 23 | 407 | 358 of 820 | 55 | 55 | 633 of 820 |
| lineage | `lasso` | 820 | 771 of 820 | 49 | 240 | 289 | 476 of 820 | 55 | 55 | 634 of 820 |
| study | `scm` | 275 | 160 of 275 | 115 | 73 | 188 | 87 of 275 | 0 | 0 | 0 of 275 |
| study | `cart` | 275 | 155 of 275 | 120 | 42 | 162 | 113 of 275 | 0 | 0 | 143 of 275 |
| study | `rf` | 275 | 162 of 275 | 113 | 12 | 125 | 150 of 275 | 0 | 0 | 237 of 275 |
| study | `gbt` | 275 | 154 of 275 | 121 | 11 | 132 | 143 of 275 | 0 | 0 | 235 of 275 |
| study | `lasso` | 275 | 260 of 275 | 15 | 67 | 82 | 193 of 275 | 0 | 0 | 241 of 275 |

The cutoffs chosen, as folds per grid point:

| Split | Learner | 0.05 | 0.1 | 0.15 | 0.2 | 0.25 | 0.3 | 0.35 | 0.4 | 0.45 | 0.5 | 0.55 | 0.6 | 0.65 | 0.7 | 0.75 | 0.8 | 0.85 | 0.9 | 0.95 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lineage | `scm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 820 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| lineage | `cart` | 16 | 31 | 32 | 22 | 29 | 25 | 18 | 18 | 6 | 449 | 9 | 8 | 10 | 16 | 21 | 40 | 21 | 26 | 23 |
| lineage | `rf` | 4 | 11 | 30 | 37 | 48 | 63 | 75 | 82 | 76 | 152 | 77 | 68 | 33 | 14 | 14 | 17 | 3 | 6 | 10 |
| lineage | `gbt` | 3 | 2 | 19 | 39 | 51 | 70 | 67 | 93 | 86 | 187 | 69 | 50 | 27 | 17 | 17 | 11 | 8 | 2 | 2 |
| lineage | `lasso` | 6 | 7 | 7 | 26 | 41 | 82 | 76 | 87 | 102 | 186 | 71 | 45 | 29 | 26 | 7 | 7 | 9 | 1 | 5 |
| study | `scm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 275 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| study | `cart` | 9 | 11 | 11 | 19 | 6 | 6 | 1 | 5 | 2 | 132 | 1 | 4 | 7 | 4 | 11 | 10 | 8 | 17 | 11 |
| study | `rf` | 5 | 9 | 11 | 10 | 11 | 25 | 25 | 28 | 25 | 38 | 27 | 16 | 9 | 11 | 9 | 4 | 8 | 1 | 3 |
| study | `gbt` | 4 | 4 | 9 | 14 | 23 | 21 | 26 | 29 | 29 | 40 | 21 | 20 | 17 | 6 | 7 | 2 | 2 | 0 | 1 |
| study | `lasso` | 2 | 2 | 3 | 10 | 16 | 32 | 36 | 35 | 31 | 34 | 22 | 14 | 14 | 13 | 2 | 3 | 2 | 2 | 2 |

## Beside it

These counts stand beside the count. They are in no Holm family and decide nothing.

**Plain best, no anchor.** The plain best against the default, counted as an arm is counted, with no Holm figure.

**Without constant models.** The count with every dataset that holds a constant-model draw left out.

**SCM.** On synthetic cohorts the default-anchored rule lost 0.0126 (0.0088 to 0.0163) for SCM on the lineage split (`docs/hyperopt-selection.md`), so SCM's search result is read beside that loss. The synthetic figure decides nothing.

**The SCM reading.** `scm`'s search on the lineage split is better on 4 of 54 datasets and worse on 9, median 0.000 (-0.001 to 0.000), beside the synthetic loss of 0.0126.

| Beside | Split | Pair | Learner | Measure | Entered | Better | Worse | Same | Median | Interval | Sign p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| plain-best | lineage | plain best, no anchor less default | `scm` | balanced accuracy | 54 of 55 | 8 | 11 | 35 | 0.000 | -0.011 to 0.001 | 0.820 |
| plain-best | lineage | plain best, no anchor less default | `scm` | recall | 54 of 55 | 16 | 8 | 30 | 0.004 | -0.002 to 0.007 | 0.076 |
| plain-best | lineage | plain best, no anchor less default | `scm` | log loss | 54 of 55 | 20 | 22 | 12 | 0.000 | -0.010 to 0.142 | 0.678 |
| plain-best | lineage | plain best, no anchor less default | `cart` | balanced accuracy | 54 of 55 | 8 | 6 | 40 | -0.001 | -0.003 to 0.000 | 0.395 |
| plain-best | lineage | plain best, no anchor less default | `cart` | recall | 54 of 55 | 20 | 7 | 27 | 0.006 | 0.000 to 0.025 | 0.010 |
| plain-best | lineage | plain best, no anchor less default | `cart` | log loss | 54 of 55 | 4 | 39 | 11 | 0.094 | 0.061 to 0.211 | 1.000 |
| plain-best | lineage | plain best, no anchor less default | `rf` | balanced accuracy | 54 of 55 | 11 | 3 | 40 | 0.001 | -0.002 to 0.005 | 0.029 |
| plain-best | lineage | plain best, no anchor less default | `rf` | recall | 54 of 55 | 15 | 14 | 25 | -0.006 | -0.011 to 0.019 | 0.500 |
| plain-best | lineage | plain best, no anchor less default | `rf` | log loss | 54 of 55 | 4 | 17 | 33 | 0.007 | 0.005 to 0.010 | 0.999 |
| plain-best | lineage | plain best, no anchor less default | `gbt` | balanced accuracy | 54 of 55 | 5 | 10 | 39 | -0.003 | -0.006 to 0.001 | 0.941 |
| plain-best | lineage | plain best, no anchor less default | `gbt` | recall | 54 of 55 | 11 | 16 | 27 | -0.001 | -0.006 to 0.006 | 0.876 |
| plain-best | lineage | plain best, no anchor less default | `gbt` | log loss | 54 of 55 | 5 | 20 | 29 | 0.013 | 0.007 to 0.024 | 1.000 |
| plain-best | lineage | plain best, no anchor less default | `lasso` | balanced accuracy | 54 of 55 | 1 | 1 | 52 | 0.000 | -0.001 to 0.000 | 0.750 |
| plain-best | lineage | plain best, no anchor less default | `lasso` | recall | 54 of 55 | 1 | 2 | 51 | 0.000 | 0.000 to 0.000 | 0.875 |
| plain-best | lineage | plain best, no anchor less default | `lasso` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.001 | 1.000 |
| plain-best | study | plain best, no anchor less default | `scm` | balanced accuracy | 55 of 55 | 9 | 16 | 30 | 0.000 | -0.009 to 0.008 | 0.946 |
| plain-best | study | plain best, no anchor less default | `scm` | recall | 55 of 55 | 13 | 16 | 26 | 0.002 | 0.000 to 0.011 | 0.771 |
| plain-best | study | plain best, no anchor less default | `scm` | log loss | 55 of 55 | 21 | 28 | 6 | 0.053 | -0.027 to 0.080 | 0.874 |
| plain-best | study | plain best, no anchor less default | `cart` | balanced accuracy | 55 of 55 | 11 | 19 | 25 | -0.009 | -0.014 to -0.001 | 0.951 |
| plain-best | study | plain best, no anchor less default | `cart` | recall | 55 of 55 | 19 | 14 | 22 | 0.002 | 0.000 to 0.010 | 0.243 |
| plain-best | study | plain best, no anchor less default | `cart` | log loss | 55 of 55 | 16 | 33 | 6 | 0.066 | 0.017 to 0.266 | 0.995 |
| plain-best | study | plain best, no anchor less default | `rf` | balanced accuracy | 55 of 55 | 16 | 8 | 31 | 0.000 | -0.006 to 0.015 | 0.076 |
| plain-best | study | plain best, no anchor less default | `rf` | recall | 55 of 55 | 13 | 17 | 25 | 0.000 | -0.003 to 0.007 | 0.819 |
| plain-best | study | plain best, no anchor less default | `rf` | log loss | 55 of 55 | 18 | 12 | 25 | -0.006 | -0.011 to 0.007 | 0.181 |
| plain-best | study | plain best, no anchor less default | `gbt` | balanced accuracy | 55 of 55 | 13 | 11 | 31 | 0.003 | -0.018 to 0.009 | 0.419 |
| plain-best | study | plain best, no anchor less default | `gbt` | recall | 55 of 55 | 14 | 13 | 28 | 0.000 | -0.013 to 0.002 | 0.500 |
| plain-best | study | plain best, no anchor less default | `gbt` | log loss | 55 of 55 | 12 | 14 | 29 | 0.004 | -0.004 to 0.015 | 0.721 |
| plain-best | study | plain best, no anchor less default | `lasso` | balanced accuracy | 55 of 55 | 1 | 0 | 54 | 0.000 | 0.000 to 0.000 | 0.500 |
| plain-best | study | plain best, no anchor less default | `lasso` | recall | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| plain-best | study | plain best, no anchor less default | `lasso` | log loss | 55 of 55 | 0 | 1 | 54 | 0.000 | 0.000 to 0.001 | 1.000 |
| no-constant | lineage | search less default | `scm` | balanced accuracy | 54 of 55 | 4 | 9 | 41 | 0.000 | -0.001 to 0.000 | 0.954 |
| no-constant | lineage | search less default | `scm` | recall | 54 of 55 | 10 | 4 | 40 | 0.000 | 0.000 to 0.004 | 0.090 |
| no-constant | lineage | search less default | `scm` | log loss | 54 of 55 | 12 | 22 | 20 | 0.000 | 0.000 to 0.026 | 0.971 |
| no-constant | lineage | search less default | `cart` | balanced accuracy | 54 of 55 | 7 | 4 | 43 | 0.000 | -0.001 to 0.000 | 0.274 |
| no-constant | lineage | search less default | `cart` | recall | 54 of 55 | 15 | 7 | 32 | 0.000 | 0.000 to 0.016 | 0.067 |
| no-constant | lineage | search less default | `cart` | log loss | 54 of 55 | 5 | 31 | 18 | 0.049 | 0.021 to 0.074 | 1.000 |
| no-constant | lineage | search less default | `rf` | balanced accuracy | 54 of 55 | 8 | 5 | 41 | 0.000 | -0.002 to 0.002 | 0.291 |
| no-constant | lineage | search less default | `rf` | recall | 54 of 55 | 10 | 9 | 35 | 0.000 | -0.006 to 0.000 | 0.500 |
| no-constant | lineage | search less default | `rf` | log loss | 54 of 55 | 5 | 10 | 39 | 0.000 | 0.000 to 0.003 | 0.941 |
| no-constant | lineage | search less default | `gbt` | balanced accuracy | 54 of 55 | 5 | 7 | 42 | 0.000 | -0.001 to 0.002 | 0.806 |
| no-constant | lineage | search less default | `gbt` | recall | 54 of 55 | 10 | 12 | 32 | 0.000 | -0.004 to 0.010 | 0.738 |
| no-constant | lineage | search less default | `gbt` | log loss | 54 of 55 | 2 | 11 | 41 | 0.004 | 0.003 to 0.008 | 0.998 |
| no-constant | lineage | search less default | `lasso` | balanced accuracy | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | search less default | `lasso` | recall | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | search less default | `lasso` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `scm` | balanced accuracy | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `scm` | recall | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `scm` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `cart` | balanced accuracy | 54 of 55 | 4 | 12 | 38 | 0.000 | -0.006 to 0.000 | 0.989 |
| no-constant | lineage | cutoff less default | `cart` | recall | 54 of 55 | 12 | 5 | 37 | 0.000 | 0.000 to 0.000 | 0.072 |
| no-constant | lineage | cutoff less default | `cart` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `rf` | balanced accuracy | 54 of 55 | 15 | 9 | 30 | 0.000 | -0.018 to 0.013 | 0.154 |
| no-constant | lineage | cutoff less default | `rf` | recall | 54 of 55 | 28 | 14 | 12 | 0.024 | 0.006 to 0.104 | 0.022 |
| no-constant | lineage | cutoff less default | `rf` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `gbt` | balanced accuracy | 54 of 55 | 8 | 12 | 34 | -0.006 | -0.012 to 0.000 | 0.868 |
| no-constant | lineage | cutoff less default | `gbt` | recall | 54 of 55 | 27 | 8 | 19 | 0.018 | 0.011 to 0.109 | < 0.001 |
| no-constant | lineage | cutoff less default | `gbt` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | cutoff less default | `lasso` | balanced accuracy | 54 of 55 | 7 | 11 | 36 | -0.004 | -0.013 to -0.001 | 0.881 |
| no-constant | lineage | cutoff less default | `lasso` | recall | 54 of 55 | 23 | 9 | 22 | 0.009 | 0.000 to 0.068 | 0.010 |
| no-constant | lineage | cutoff less default | `lasso` | log loss | 54 of 55 | 0 | 0 | 54 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | lineage | twin less default | `scm` | balanced accuracy | 54 of 55 | 7 | 8 | 39 | 0.002 | 0.000 to 0.005 | 0.696 |
| no-constant | lineage | twin less default | `scm` | recall | 54 of 55 | 14 | 8 | 32 | 0.003 | 0.000 to 0.005 | 0.143 |
| no-constant | lineage | twin less default | `scm` | log loss | 54 of 55 | 22 | 19 | 13 | -0.005 | -0.021 to 0.028 | 0.378 |
| no-constant | lineage | twin less default | `cart` | balanced accuracy | 54 of 55 | 6 | 9 | 39 | -0.002 | -0.006 to 0.000 | 0.849 |
| no-constant | lineage | twin less default | `cart` | recall | 54 of 55 | 10 | 13 | 31 | 0.000 | -0.004 to 0.000 | 0.798 |
| no-constant | lineage | twin less default | `cart` | log loss | 54 of 55 | 20 | 16 | 18 | -0.002 | -0.019 to 0.005 | 0.309 |
| no-constant | lineage | twin less default | `rf` | balanced accuracy | 54 of 55 | 4 | 4 | 46 | 0.000 | -0.002 to 0.000 | 0.637 |
| no-constant | lineage | twin less default | `rf` | recall | 54 of 55 | 6 | 9 | 39 | -0.002 | -0.008 to -0.001 | 0.849 |
| no-constant | lineage | twin less default | `rf` | log loss | 54 of 55 | 2 | 3 | 49 | -0.001 | -0.004 to 0.001 | 0.812 |
| no-constant | lineage | twin less default | `gbt` | balanced accuracy | 54 of 55 | 8 | 7 | 39 | 0.000 | -0.007 to 0.009 | 0.500 |
| no-constant | lineage | twin less default | `gbt` | recall | 54 of 55 | 12 | 14 | 28 | -0.002 | -0.008 to 0.014 | 0.721 |
| no-constant | lineage | twin less default | `gbt` | log loss | 54 of 55 | 12 | 8 | 34 | -0.003 | -0.011 to 0.002 | 0.252 |
| no-constant | lineage | twin less default | `lasso` | balanced accuracy | 54 of 55 | 6 | 2 | 46 | 0.000 | -0.002 to 0.010 | 0.145 |
| no-constant | lineage | twin less default | `lasso` | recall | 54 of 55 | 14 | 5 | 35 | 0.000 | -0.004 to 0.016 | 0.032 |
| no-constant | lineage | twin less default | `lasso` | log loss | 54 of 55 | 15 | 13 | 26 | 0.001 | -0.021 to 0.007 | 0.425 |
| no-constant | study | search less default | `scm` | balanced accuracy | 55 of 55 | 8 | 10 | 37 | 0.000 | 0.000 to 0.000 | 0.760 |
| no-constant | study | search less default | `scm` | recall | 55 of 55 | 10 | 10 | 35 | 0.000 | 0.000 to 0.005 | 0.588 |
| no-constant | study | search less default | `scm` | log loss | 55 of 55 | 14 | 18 | 23 | 0.000 | 0.000 to 0.000 | 0.811 |
| no-constant | study | search less default | `cart` | balanced accuracy | 55 of 55 | 5 | 12 | 38 | -0.002 | -0.004 to 0.000 | 0.975 |
| no-constant | study | search less default | `cart` | recall | 55 of 55 | 9 | 8 | 38 | 0.000 | 0.000 to 0.000 | 0.500 |
| no-constant | study | search less default | `cart` | log loss | 55 of 55 | 17 | 24 | 14 | 0.000 | 0.000 to 0.187 | 0.894 |
| no-constant | study | search less default | `rf` | balanced accuracy | 55 of 55 | 12 | 4 | 39 | 0.000 | -0.004 to 0.014 | 0.038 |
| no-constant | study | search less default | `rf` | recall | 55 of 55 | 10 | 10 | 35 | 0.000 | 0.000 to 0.000 | 0.588 |
| no-constant | study | search less default | `rf` | log loss | 55 of 55 | 12 | 11 | 32 | 0.000 | -0.001 to 0.007 | 0.500 |
| no-constant | study | search less default | `gbt` | balanced accuracy | 55 of 55 | 6 | 10 | 39 | 0.000 | -0.018 to 0.000 | 0.895 |
| no-constant | study | search less default | `gbt` | recall | 55 of 55 | 9 | 10 | 36 | 0.000 | -0.012 to 0.000 | 0.676 |
| no-constant | study | search less default | `gbt` | log loss | 55 of 55 | 4 | 11 | 40 | 0.000 | 0.000 to 0.008 | 0.982 |
| no-constant | study | search less default | `lasso` | balanced accuracy | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | search less default | `lasso` | recall | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | search less default | `lasso` | log loss | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `scm` | balanced accuracy | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `scm` | recall | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `scm` | log loss | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `cart` | balanced accuracy | 55 of 55 | 1 | 11 | 43 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `cart` | recall | 55 of 55 | 7 | 6 | 42 | 0.000 | 0.000 to 0.001 | 0.500 |
| no-constant | study | cutoff less default | `cart` | log loss | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `rf` | balanced accuracy | 55 of 55 | 20 | 6 | 29 | 0.007 | -0.003 to 0.011 | 0.005 |
| no-constant | study | cutoff less default | `rf` | recall | 55 of 55 | 29 | 9 | 17 | 0.023 | 0.005 to 0.052 | < 0.001 |
| no-constant | study | cutoff less default | `rf` | log loss | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `gbt` | balanced accuracy | 55 of 55 | 5 | 14 | 36 | -0.007 | -0.013 to -0.004 | 0.990 |
| no-constant | study | cutoff less default | `gbt` | recall | 55 of 55 | 25 | 9 | 21 | 0.017 | 0.000 to 0.020 | 0.005 |
| no-constant | study | cutoff less default | `gbt` | log loss | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |
| no-constant | study | cutoff less default | `lasso` | balanced accuracy | 55 of 55 | 8 | 12 | 35 | -0.002 | -0.004 to 0.001 | 0.868 |
| no-constant | study | cutoff less default | `lasso` | recall | 55 of 55 | 25 | 12 | 18 | 0.013 | 0.000 to 0.029 | 0.024 |
| no-constant | study | cutoff less default | `lasso` | log loss | 55 of 55 | 0 | 0 | 55 | 0.000 | 0.000 to 0.000 | 1.000 |

| Split | Pair | Entered | Better | Worse | Same | Median | Interval | Holm p | Synthetic loss of the rule, default less rule |
|---|---|---|---|---|---|---|---|---|---|
| lineage | search less default | 54 of 55 | 4 | 9 | 41 | 0.000 | -0.001 to 0.000 | 1.000 | 0.0126 (0.0088 to 0.0163) |
| study | search less default | 55 of 55 | 8 | 10 | 37 | 0.000 | 0.000 to 0.000 | none | none |
| lineage | plain best, no anchor less default | 54 of 55 | 8 | 11 | 35 | 0.000 | -0.011 to 0.001 | none | 0.0126 (0.0088 to 0.0163) |
| study | plain best, no anchor less default | 55 of 55 | 9 | 16 | 30 | 0.000 | -0.009 to 0.008 | none | none |

## Seen after

These points are seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**Species.** Per species, named in full, the rows are dataset and learner pairs, with their denominator. They decide nothing.

| Species | Split | Pair | Datasets | Better / worse / same, dataset and learner pairs | Not scored |
|---|---|---|---|---|---|
| *Acinetobacter baumannii* | lineage | search less default | 2 | 2 / 3 / 5 of 10 | 0 |
| *Acinetobacter baumannii* | lineage | cutoff less default | 2 | 3 / 2 / 5 of 10 | 0 |
| *Acinetobacter baumannii* | lineage | twin less default | 2 | 4 / 0 / 6 of 10 | 0 |
| *Acinetobacter baumannii* | lineage | plain best, no anchor less default | 2 | 3 / 2 / 5 of 10 | 0 |
| *Acinetobacter baumannii* | study | search less default | 2 | 4 / 1 / 5 of 10 | 0 |
| *Acinetobacter baumannii* | study | cutoff less default | 2 | 1 / 1 / 8 of 10 | 0 |
| *Acinetobacter baumannii* | study | plain best, no anchor less default | 2 | 5 / 1 / 4 of 10 | 0 |
| *Enterobacter cloacae* | lineage | search less default | 8 | 3 / 6 / 31 of 40 | 0 |
| *Enterobacter cloacae* | lineage | cutoff less default | 8 | 9 / 8 / 23 of 40 | 0 |
| *Enterobacter cloacae* | lineage | twin less default | 8 | 10 / 5 / 25 of 40 | 0 |
| *Enterobacter cloacae* | lineage | plain best, no anchor less default | 8 | 5 / 7 / 28 of 40 | 0 |
| *Enterobacter cloacae* | study | search less default | 8 | 4 / 13 / 23 of 40 | 0 |
| *Enterobacter cloacae* | study | cutoff less default | 8 | 6 / 9 / 25 of 40 | 0 |
| *Enterobacter cloacae* | study | plain best, no anchor less default | 8 | 5 / 15 / 20 of 40 | 0 |
| *Escherichia coli* | lineage | search less default | 15 | 8 / 3 / 64 of 75 | 0 |
| *Escherichia coli* | lineage | cutoff less default | 15 | 11 / 10 / 54 of 75 | 0 |
| *Escherichia coli* | lineage | twin less default | 15 | 5 / 5 / 65 of 75 | 0 |
| *Escherichia coli* | lineage | plain best, no anchor less default | 15 | 13 / 6 / 56 of 75 | 0 |
| *Escherichia coli* | study | search less default | 15 | 10 / 8 / 57 of 75 | 0 |
| *Escherichia coli* | study | cutoff less default | 15 | 8 / 8 / 59 of 75 | 0 |
| *Escherichia coli* | study | plain best, no anchor less default | 15 | 10 / 15 / 50 of 75 | 0 |
| *Klebsiella pneumoniae* | lineage | search less default | 18 | 3 / 3 / 79 of 90 | 5 |
| *Klebsiella pneumoniae* | lineage | cutoff less default | 18 | 3 / 9 / 73 of 90 | 5 |
| *Klebsiella pneumoniae* | lineage | twin less default | 18 | 6 / 9 / 70 of 90 | 5 |
| *Klebsiella pneumoniae* | lineage | plain best, no anchor less default | 18 | 5 / 3 / 77 of 90 | 5 |
| *Klebsiella pneumoniae* | study | search less default | 18 | 7 / 5 / 78 of 90 | 0 |
| *Klebsiella pneumoniae* | study | cutoff less default | 18 | 12 / 15 / 63 of 90 | 0 |
| *Klebsiella pneumoniae* | study | plain best, no anchor less default | 18 | 13 / 12 / 65 of 90 | 0 |
| *Mycobacterium tuberculosis* | lineage | search less default | 1 | 1 / 1 / 3 of 5 | 0 |
| *Mycobacterium tuberculosis* | lineage | cutoff less default | 1 | 1 / 2 / 2 of 5 | 0 |
| *Mycobacterium tuberculosis* | lineage | twin less default | 1 | 1 / 0 / 4 of 5 | 0 |
| *Mycobacterium tuberculosis* | lineage | plain best, no anchor less default | 1 | 1 / 2 / 2 of 5 | 0 |
| *Mycobacterium tuberculosis* | study | search less default | 1 | 0 / 1 / 4 of 5 | 0 |
| *Mycobacterium tuberculosis* | study | cutoff less default | 1 | 0 / 1 / 4 of 5 | 0 |
| *Mycobacterium tuberculosis* | study | plain best, no anchor less default | 1 | 1 / 2 / 2 of 5 | 0 |
| *Pseudomonas aeruginosa* | lineage | search less default | 6 | 5 / 4 / 21 of 30 | 0 |
| *Pseudomonas aeruginosa* | lineage | cutoff less default | 6 | 1 / 10 / 19 of 30 | 0 |
| *Pseudomonas aeruginosa* | lineage | twin less default | 6 | 5 / 6 / 19 of 30 | 0 |
| *Pseudomonas aeruginosa* | lineage | plain best, no anchor less default | 6 | 3 / 5 / 22 of 30 | 0 |
| *Pseudomonas aeruginosa* | study | search less default | 6 | 4 / 4 / 22 of 30 | 0 |
| *Pseudomonas aeruginosa* | study | cutoff less default | 6 | 2 / 6 / 22 of 30 | 0 |
| *Pseudomonas aeruginosa* | study | plain best, no anchor less default | 6 | 10 / 5 / 15 of 30 | 0 |
| *Salmonella enterica* | lineage | search less default | 1 | 0 / 2 / 3 of 5 | 0 |
| *Salmonella enterica* | lineage | cutoff less default | 1 | 0 / 1 / 4 of 5 | 0 |
| *Salmonella enterica* | lineage | twin less default | 1 | 0 / 2 / 3 of 5 | 0 |
| *Salmonella enterica* | lineage | plain best, no anchor less default | 1 | 0 / 2 / 3 of 5 | 0 |
| *Salmonella enterica* | study | search less default | 1 | 0 / 2 / 3 of 5 | 0 |
| *Salmonella enterica* | study | cutoff less default | 1 | 0 / 2 / 3 of 5 | 0 |
| *Salmonella enterica* | study | plain best, no anchor less default | 1 | 0 / 1 / 4 of 5 | 0 |
| *Staphylococcus aureus* | lineage | search less default | 4 | 2 / 3 / 15 of 20 | 0 |
| *Staphylococcus aureus* | lineage | cutoff less default | 4 | 6 / 2 / 12 of 20 | 0 |
| *Staphylococcus aureus* | lineage | twin less default | 4 | 0 / 3 / 17 of 20 | 0 |
| *Staphylococcus aureus* | lineage | plain best, no anchor less default | 4 | 3 / 4 / 13 of 20 | 0 |
| *Staphylococcus aureus* | study | search less default | 4 | 2 / 2 / 16 of 20 | 0 |
| *Staphylococcus aureus* | study | cutoff less default | 4 | 5 / 1 / 14 of 20 | 0 |
| *Staphylococcus aureus* | study | plain best, no anchor less default | 4 | 6 / 3 / 11 of 20 | 0 |

**Recall and specificity under the cutoff arm.** On the lineage split the cutoff arm raises recall for `rf`, `gbt` and `lasso` (medians 0.024, 0.018 and 0.009, better on 28, 27 and 23 of 54 datasets) and lowers specificity (medians -0.028, -0.037 and -0.025, worse on 30, 33 and 29). The study split turns the same way. The cutoff moved from 0.5 on 668 of 820 lineage folds for `rf`, 633 for `gbt`, 634 for `lasso` and 371 for `cart`, and on 0 of 820 for `scm`.

**Trial 0.** The search kept trial 0, the release's own setting, on 771 of 820 lineage folds for `lasso` and on 436 to 515 of 820 for the other four. Trial 0 and every trial failed inside on 55 of 820 lineage folds for each learner, and on 0 of 275 study folds. These tables do not say why.

**Largest moves.** On lineage balanced accuracy the search moves 13 of 270 dataset and learner pairs by more than 0.05, from -0.094 (*Staphylococcus aureus* with clindamycin, `scm`, 0.633 to 0.538) to 0.303 (*Escherichia coli* with cefepime, `rf`, 0.527 to 0.830). The cutoff arm moves 38 of 270, from -0.171 (*Klebsiella pneumoniae* with aztreonam, `lasso`, 0.734 to 0.563, recall 0.693 to 0.964) to 0.311 (*Escherichia coli* with cefepime, `rf`, 0.527 to 0.838). The twin, the seed alone, moves 17 of 270, from -0.094 (*Klebsiella pneumoniae* with cefoxitin, `cart`) to 0.169 (*Mycobacterium tuberculosis* with nicotinamide, `gbt`).

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Study split.** One study run per dataset at the release's seeds, and no study twin.

**Labels.** The 14 datasets of M wait for the label rule. The count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged.

**The cutoff arm.** The cutoff arm moves no score, so its log loss equals the default's by construction. A same on log loss for the cutoff arm says nothing about its calls.

**Times.** A time can equal a genome id as text, so the id check counts apart, by file kind and key or label, a token that is the whole value of `seconds_inner` or `seconds_outer` in a choices file or the value on a `User time (seconds)`, `System time (seconds)` or `Elapsed (wall clock) time` line of a drive's `.time` file, and every other hit fails it.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits, arms, folds and learners only.

**Binaries.** All 8 species ran one binary by its code sha256, `343172462952f291`. Every species' `binary.txt` names commit `0d987ba` and code sha256 `343172462952f291381d765a54530c7f913bf101b04b4f9ac5a6ad27a16f9ee8`.

## Method

**Checks of the run.** Before any score was read, `scripts/c5-settings/run_checks.py check` read each species' copy and printed counts only: the first run's species, each once; each species ended `DONE`, a `FAILED` older than its `DONE` named and set aside; each species' `SHA256SUMS`; `score/score.log` kept out of it; each `DONE` against its `binary.txt`, its results and its choices files; one commit and one code sha256; each `lists.txt`; each `expect.tsv` against `expect.py`; every default and twin result and one choices file per outer fold, and no other; the refusal records and the constant-model rows, counted; and the genome ids over every file that will be read or published and every log. `score/score.log` stays out of the record.

*Klebsiella pneumoniae*'s run stopped at its score step and wrote no `DONE` and no `SHA256SUMS`. `run_checks.py check --failed-at-score klebsiella_pneumoniae` holds its record to four conditions instead, each printed with its counts: its `FAILED` is the at score line in the run's time form, with no `DONE` and no `SHA256SUMS`; `logs/run.log` names as many drives ended ok as it had to run, and none ended otherwise; every line of `SHA256SUMS.lead` holds, its own sha256 is `bd0d9b38edcc6b2cfd0a9c1e80827699d015308ddca1177dac094fadeae2cf7c`, and it names no `score/score.log`; every file of the copy, wherever it lies, is named in it, other than `score/score.log`, `back.tar` and the list itself, and every name in it is a file of the copy; and its `binary.txt` names the run's commit and code sha256. Without the allowance, a `FAILED` with no newer `DONE` fails.

**The scorer after the run.** `scripts/c5-settings/score.py` changed after the run and before any score was opened. The run's scorer, at commit `0d987ba`, stopped *Klebsiella pneumoniae* at its score step on one failed check, "the default's result at the same dataset, split and seed is there and available", for polymyxin B on the lineage split at seed 43, whose default result is unavailable: it read as a failed check a draw that C3's scorer reads as not scored. The changed scorer makes each arm's row there not scored, for the default's reason, by handing C3's scorer the same unavailable result it reads for the default, and a choices file beside such a draw fails a check. For a species with no unavailable default draw it writes `checks.tsv`, `draws.tsv` and `datasets.tsv` byte for byte as `0d987ba`'s scorer does, as a SYNTHETIC test holds on made-up species.

The 7 species that ended `DONE`, scored again by the changed scorer from their copies, give 21 of 21 tables the same byte for byte as their own `SHA256SUMS` lines: `rescore.txt` ends `SAME 21 of 21`. *Klebsiella pneumoniae*'s three tables are the changed scorer's alone, as its run wrote none.

**One commit.** The fits are the run's, at commit `0d987ba`, one code sha256 for all 8 species. The scoring is this branch's `scripts/c5-settings/score.py`.

The checks' printout, `reference/settings-arms/checks.txt`:

```
PASS the first run's species, each once: 8 folders for 8 species; missing [], more []
PASS each species ended DONE, any FAILED older than its DONE: 7 of 8; klebsiella_pneumoniae stopped at its score step, under the allowance
PASS each species' SHA256SUMS holds every file of the run's rule, each at its sha256: 1069 files
PASS score/score.log stays out of SHA256SUMS
PASS each DONE names its species, its binary.txt's binary and code, its results, its choices files and its SHA256SUMS
PASS every species ran one commit and one code sha256: 1 commits, 1 code sha256, 1 binaries
PASS klebsiella_pneumoniae: FAILED is its at score line in the run's time form, with no DONE and no SHA256SUMS: FAILED at score 2026-11-23T22:15:43Z; DONE 0, SHA256SUMS 0
PASS klebsiella_pneumoniae: logs/run.log names as many drives ended ok as it had to run, and none ended otherwise: 355 to run over 1 drives lines, 355 ended ok (355 named once), 0 ended otherwise
PASS klebsiella_pneumoniae: every line of SHA256SUMS.lead holds, its own sha256 is the one given, and it names no score/score.log: 1269 lines, 0 differ or absent, score/score.log named 0; its sha256 is the one given
PASS klebsiella_pneumoniae: every file of its copy is named in SHA256SUMS.lead, other than score/score.log, back.tar and the list itself: 1269 files, 0 not named
PASS klebsiella_pneumoniae: every name in SHA256SUMS.lead is a file of its copy: 1269 names other than score/score.log, 0 with no file
PASS klebsiella_pneumoniae: binary.txt names the run's commit and code sha256: commit 0d987ba8572430588e738aa873b6690a24058d97, code 343172462952f291381d765a54530c7f913bf101b04b4f9ac5a6ad27a16f9ee8; the other species name 1 commits and 1 code sha256
PASS each lists.txt names its species' first-run datasets, seeds 42, 43 and 44 and twin 45, 46 and 47
PASS each expect.tsv is expect.py's list from lists.tsv: 1045 results
PASS each species holds every default and twin result and one choices file per outer fold of each draw and study run of the default, and no other: 1480 files
PASS the refusal records and the constant-model rows, counted: five constant-model rows per refusal record: 0 refusal records, 0 constant-model rows
PASS no genome id in any file that will be read or published: 0 genome ids in 4138 files, against 13393 ids; times equal to an id as text, counted apart: time User time (seconds) 1
CHECKS 17 of 17
```

**Reading.** The answers read are each species' three tables as the changed scorer writes them from its copy, the 7 finished species' tables held first to their own `SHA256SUMS` lines, and each species' choices files as its run wrote them. `check --tables` and `gather --tables` read no rescored table unless `rescore.txt` ends `SAME <n> of <n>` with no `FAIL` line and n is 3 tables for each species that ended `DONE`. The count runs once over them:

```
python3 -I scripts/c5-settings/rescore.py <copies> scripts/c5-settings/lists.tsv <release> rescored
python3 -I scripts/c5-settings/run_checks.py check <copies> scripts/c5-settings/lists.tsv <release> read/checks.txt --failed-at-score klebsiella_pneumoniae --lead-sums bd0d9b38edcc6b2cfd0a9c1e80827699d015308ddca1177dac094fadeae2cf7c --tables rescored
python3 -I scripts/c5-settings/run_checks.py gather <copies> scripts/c5-settings/lists.tsv read/checks.txt read --tables rescored
python3 -I scripts/c5-settings/score.py --count read/count read/score/<species>/datasets.tsv   # 8 files
python3 -I scripts/c5-settings/page_tables.py write read scripts/c5-settings/lists.tsv read/checks.txt reference/settings-arms
python3 -I scripts/c5-settings/page_tables.py print reference/settings-arms
```

**Tables.** `scripts/c5-settings/page_tables.py write` checks the reading against sixteen gates and writes `reference/settings-arms/` only when every gate passes:

| Gate | What must hold |
|---|---|
| 1. files | the reading holds reading.txt, count/ with score.py --count's four tables, and score/; no score.log is in it; the run's checks.txt says every check passed |
| 2. species | score/ holds each species of lists.tsv's first run once and no other, each with its five files and choices/, and their datasets sum to the first run's |
| 3. datasets | each datasets.tsv holds exactly its species' first-run datasets, none in M, with each arm on each split for the five learners; species and antibiotic equal lists.tsv's |
| 4. study | every study row is a study run at seed 42, and no study row stands on the lineage split |
| 5. twin | the twin stands on the lineage split at seeds 45, 46 and 47 alone, no choices file is at a twin seed, and the twin is never read as an arm |
| 6. plain | the plain best stands in count-plain-best.tsv alone, with no Holm p, and in no condition and no other count |
| 7. missing | a draw not scored has no figure, and a dataset with a draw not scored is not scored: a missing draw is never read as zero |
| 8. constant | each constant-model draw is an arm's study run refused for no candidate, kept, scored and marked |
| 9. draws | each draws.tsv holds exactly the draws of the run, and datasets.tsv is their mean |
| 10. checks | every check of each checks.tsv passed |
| 11. cutoff | the cutoff arm's log loss is the default's on every dataset, split and learner; scm's cutoff arm is the default on every fold |
| 12. choices | each draw and study run of the default whose arms are scored has one choices file per outer fold, as many as the default's result has (folds.tsv, which `gather` writes), each pick by the default-anchored rule and each cutoff a grid point; a draw with no choices file is a missing draw |
| 13. settings | the counts of what the search chose hold trial ids, grid points and fold counts alone: no trial's settings and no row per genome |
| 14. binary | each binary.txt names a commit, a binary and a code sha256; the page names each binary with its species unless every species ran one |
| 15. conditions | conditions 1 and 2, 3 and 4, as read here off count.tsv, equal conditions.tsv's |
| 16. count | the four count tables are what score.py --count writes from the datasets.tsv files |

**Tests.** `tests/c5_settings/test_page_tables.py` and `tests/c5_settings/test_run_checks.py` hold SYNTHETIC tests of the table script and the check script, on three made-up species. `scripts/c5-settings/page_faults.py` plants each fault in memory, in turn, and runs the test that must catch it: every gate of the table script has a fault, and `tests/c5_settings/test_page_faults.py` checks that each fault's text occurs once in its file and that each fault is caught. `tests/c5_settings/test_unavailable.py` holds the scorer at an unavailable default draw, and the changed scorer against `0d987ba`'s own (`tests/c5_settings/score_at_0d987ba.txt`) on made-up species.

`tests/c5_settings` holds 156 SYNTHETIC tests, run at `PYTHONHASHSEED` 0, 1 and 4242 and with `-I`. `scripts/c5-settings/page_faults.py` caught 49 of 49 faults: 20 of the table script's gates (P), 9 of its tables (T) and 20 of the run's checks and `rescore.py` (R). The block's own `scripts/c5-settings/planted_faults.py` caught 18 of 18 faults of the driver, the scorer and the binary (1 to 11, with 2a and 2b, 6a and 6b, 7a and 7b, and 12a to 12d).

| Fault | Planted | Caught by | Caught |
|---|---|---|---|
| P1 | a score.log in the reading passes | `Gates.test_a_score_log_in_the_reading_is_refused` | yes |
| P2 | the species are those found, not those of the first run | `Gates.test_a_species_missing_from_the_reading_is_refused` | yes |
| P3 | a dataset of M on the first run's list is counted | `Gates.test_a_dataset_of_m_on_the_list_is_never_counted` | yes |
| P4 | a study row on the lineage split passes | `Gates.test_a_study_row_on_the_lineage_split_is_refused` | yes |
| P5 | the twin read as an arm passes | `Gates.test_the_twin_read_as_an_arm_is_refused` | yes |
| P6 | the plain best in a Holm family passes | `Gates.test_the_plain_best_in_a_holm_family_is_refused` | yes |
| P7 | a draw not scored keeps a figure | `Gates.test_a_missing_draw_read_as_zero_is_refused` | yes |
| P8 | a constant-model row left unmarked passes | `Gates.test_an_unmarked_constant_model_row_is_refused` | yes |
| P9 | a draw listed twice passes | `Gates.test_a_draw_twice_is_refused` | yes |
| P10 | a failed check passes | `Gates.test_a_failed_check_is_refused` | yes |
| P11 | the cutoff arm's log loss apart from the default's passes | `Gates.test_the_cutoff_arms_log_loss_apart_from_the_defaults_is_refused` | yes |
| P11b | scm's cutoff moved from 0.5 passes | `Gates.test_scms_cutoff_moved_from_half_is_refused` | yes |
| P12 | a choices file missing inside a scored draw passes | `Gates.test_a_choices_file_missing_inside_a_scored_draw_is_refused` | yes |
| P12c | a scored draw short of its default's outer folds passes | `Gates.test_a_scored_draw_short_of_its_defaults_folds_is_refused` | yes |
| P12b | a pick against the default-anchored rule passes | `Gates.test_a_pick_against_the_default_anchored_rule_is_refused` | yes |
| P13 | a setting in the counts of what the search chose passes | `Gates.test_a_setting_in_the_counts_of_what_the_search_chose_is_refused` | yes |
| P14 | a binary.txt without its code passes | `Gates.test_a_binary_txt_without_its_code_is_refused` | yes |
| P15 | conditions.tsv apart from the count passes | `Gates.test_conditions_tsv_apart_from_the_count_is_refused` | yes |
| P16 | a count that is not score.py's passes | `Gates.test_a_count_that_is_not_score_pys_is_refused` | yes |
| T1 | the stricter reading decides | `Tables.test_ruling_15_decides_and_the_stricter_reading_stands_beside` | yes |
| T2 | the cutoff arm's log loss rows decide | `Tables.test_the_cutoff_arms_log_loss_rows_decide_nothing` | yes |
| T3 | the plain best decides | `Tables.test_the_plain_best_and_the_count_without_constant_models_decide_nothing` | yes |
| T4 | a missing value read as zero in the datasets table | `Tables.test_a_missing_draw_is_never_read_as_zero_in_the_tables` | yes |
| T5 | a missing figure printed as zero | `Tables.test_a_missing_draw_is_never_read_as_zero_in_the_tables` | yes |
| T6 | a species named by its slug | `Tables.test_each_species_is_named_in_full_and_decides_nothing` | yes |
| T7 | a trial's settings ride with the counts | `Tables.test_the_counts_of_what_the_search_chose_hold_no_setting` | yes |
| T8 | two binaries named as one | `Tables.test_one_binary_is_named_once_and_two_with_their_species` | yes |
| R1 | a species more than the first run's passes | `Check.test_a_species_more_than_the_first_runs_is_refused` | yes |
| R2 | a FAILED newer than its DONE is set aside | `Check.test_a_failed_file_newer_than_done_is_refused` | yes |
| R3 | a FAILED older than its DONE is not named | `Check.test_a_failed_file_older_than_done_is_named_and_set_aside` | yes |
| R4 | a file changed after its SHA256SUMS passes | `Check.test_a_file_changed_after_its_sha256sums_is_refused` | yes |
| R5 | score.log in SHA256SUMS passes | `Check.test_score_log_in_sha256sums_is_refused` | yes |
| R6 | two code sha256 pass | `Check.test_two_code_sha256_are_refused` | yes |
| R7 | an expect.tsv that is not expect.py's passes | `Check.test_an_expect_tsv_that_is_not_expect_pys_is_refused` | yes |
| R8 | a missing choices file passes | `Check.test_a_missing_choices_file_is_refused` | yes |
| R9 | the genome id count skips the choices files | `Check.test_a_genome_id_in_a_choices_file_is_refused` | yes |
| R10 | a refusal record with no constant-model rows passes | `Check.test_a_refusal_record_with_no_constant_model_rows_is_refused` | yes |
| R11 | gather copies score.log | `Gather.test_gather_copies_four_files_a_species_and_the_choices_and_never_score_log` | yes |
| P12d | an unavailable default that owes choices files passes | `UnavailableDraw.test_an_unavailable_default_that_owes_choices_files_is_refused` | yes |
| T9 | a column of 0 in every row printed | `UnavailableDraw.test_a_column_of_0_in_every_row_is_not_printed` | yes |
| R12 | a species with no DONE passes as stopped at score without the allowance | `Allowance.test_a_species_stopped_at_score_fails_without_the_allowance` | yes |
| R13 | SHA256SUMS.lead of another sha256 passes | `Allowance.test_the_allowance_refuses_a_lead_sums_of_another_sha256` | yes |
| R14 | a drive ended otherwise passes the allowance | `Allowance.test_the_allowance_refuses_a_drive_ended_otherwise` | yes |
| R15 | an id on any line of a drive's .time file counted as a time | `Times.test_an_id_on_another_line_of_a_time_file_still_fails` | yes |
| R16 | a seconds_inner value equal to an id fails as an id | `Times.test_a_seconds_inner_value_equal_to_an_id_is_counted_apart` | yes |
| R17 | gather with tables copies the copy's tables | `Tables.test_gather_with_tables_takes_the_rescored_tables` | yes |
| R18 | a rescored table unlike its SHA256SUMS line passes | `Tables.test_rescore_refuses_a_table_unlike_its_sha256sums_line` | yes |
| R19 | a file of the stopped species' copy its SHA256SUMS.lead does not name passes | `Allowance.test_the_allowance_refuses_a_results_file_the_lead_sums_do_not_name` | yes |
| R20 | a held count other than 3 tables for each species that ended DONE passes | `Tables.test_check_and_gather_refuse_a_held_count_other_than_3_tables_for_each_species_done` | yes |

**Files.** The reading's tables are in `reference/settings-arms/`. `scripts/c5-settings/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

**Columns in the files only.** A column that reads 0 in every row is not printed and stays in its file only: where no draw is a constant model, the constant-model column of `page-count.tsv` and the constant-model draws and study run refused columns of `page-datasets.tsv`.

## Datasets

The first run: 55 datasets in 8 species, with each dataset's genomes as `scripts/c5-settings/lists.tsv` gives them.

| Species | Datasets | Antibiotics, with genomes |
|---|---:|---|
| *Acinetobacter baumannii* | 2 | cefepime (437), minocycline (139) |
| *Enterobacter cloacae* | 8 | cefepime (103), ciprofloxacin (110), gentamicin (104), imipenem (86), levofloxacin (104), meropenem (99), sulfamethoxazole/trimethoprim (80), tobramycin (99) |
| *Escherichia coli* | 15 | amikacin (1,125), ampicillin/sulbactam (283), aztreonam (699), cefazolin (414), cefepime (1,124), cefoxitin (1,066), ceftriaxone (923), ertapenem (820), levofloxacin (403), nitrofurantoin (482), piperacillin/tazobactam (979), streptomycin (386), sulfamethoxazole (341), tetracycline (782), tobramycin (1,070) |
| *Klebsiella pneumoniae* | 18 | ampicillin/sulbactam (744), aztreonam (1,050), cefazolin (1,141), cefepime (1,233), cefoxitin (1,181), ceftazidime (1,310), ceftriaxone (1,208), ciprofloxacin (1,706), colistin (336), doripenem (243), ertapenem (952), imipenem (1,158), levofloxacin (871), nitrofurantoin (562), polymyxin b (100), sulfamethoxazole/trimethoprim (1,546), tigecycline (924), tobramycin (1,230) |
| *Mycobacterium tuberculosis* | 1 | nicotinamide (166) |
| *Pseudomonas aeruginosa* | 6 | aztreonam (235), cefepime (245), ceftazidime/avibactam (309), gentamicin (346), imipenem (447), piperacillin/tazobactam (210) |
| *Salmonella enterica* | 1 | sulfisoxazole (544) |
| *Staphylococcus aureus* | 4 | clindamycin (1,543), erythromycin (1,957), penicillin (1,508), sulfamethoxazole/trimethoprim (1,389) |

Per dataset, each arm's difference from the default on balanced accuracy, for `scm`, `cart`, `rf`, `gbt` and `lasso` in that order; `none` where a draw is not scored. The values, recall, log loss and the draws scored are in `reference/settings-arms/page-datasets.tsv`.

| Species | Antibiotic | Genomes | Split | Search: five learners | Cutoff: five learners | Twin: five learners | Plain best: five learners |
|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | cefepime | 437 | lineage | -0.039 / 0.038 / 0.029 / -0.005 / 0.000 | 0.000 / -0.007 / 0.095 / 0.067 / 0.028 | 0.005 / -0.007 / 0.000 / 0.029 / 0.053 | -0.039 / 0.038 / 0.032 / -0.003 / -0.001 |
| *Acinetobacter baumannii* | cefepime | 437 | study | -0.012 / -0.004 / 0.123 / -0.011 / 0.000 | 0.000 / 0.000 / -0.018 / -0.016 / 0.000 |  | -0.012 / -0.001 / 0.123 / 0.045 / 0.000 |
| *Acinetobacter baumannii* | minocycline | 139 | lineage | 0.017 / -0.083 / -0.027 / 0.001 / 0.000 | 0.000 / -0.027 / -0.053 / 0.006 / 0.004 | 0.072 / 0.027 / 0.000 / 0.000 / 0.014 | 0.034 / -0.090 / -0.014 / 0.001 / 0.000 |
| *Acinetobacter baumannii* | minocycline | 139 | study | 0.147 / -0.111 / 0.106 / 0.049 / 0.000 | 0.000 / -0.152 / 0.021 / -0.001 / 0.014 |  | 0.133 / -0.111 / 0.106 / 0.049 / 0.000 |
| *Enterobacter cloacae* | cefepime | 103 | lineage | -0.021 / 0.005 / -0.002 / 0.010 / 0.000 | 0.000 / 0.000 / -0.019 / -0.012 / -0.020 | 0.009 / -0.005 / -0.007 / 0.001 / 0.047 | -0.054 / 0.000 / -0.002 / 0.005 / 0.000 |
| *Enterobacter cloacae* | cefepime | 103 | study | 0.000 / 0.000 / 0.037 / -0.037 / 0.000 | 0.000 / 0.000 / 0.114 / 0.015 / 0.023 |  | 0.000 / -0.029 / 0.045 / -0.052 / 0.000 |
| *Enterobacter cloacae* | ciprofloxacin | 110 | lineage | -0.031 / -0.006 / -0.027 / 0.002 / 0.000 | 0.000 / 0.064 / -0.027 / -0.038 / -0.008 | 0.019 / -0.030 / 0.037 / 0.048 / -0.001 | -0.044 / -0.006 / -0.027 / 0.010 / -0.010 |
| *Enterobacter cloacae* | ciprofloxacin | 110 | study | -0.118 / -0.012 / -0.016 / -0.019 / 0.000 | 0.000 / -0.026 / -0.037 / -0.027 / 0.000 |  | -0.118 / -0.012 / -0.012 / -0.041 / 0.000 |
| *Enterobacter cloacae* | gentamicin | 104 | lineage | 0.013 / 0.006 / -0.022 / 0.052 / 0.000 | 0.000 / -0.033 / -0.020 / 0.026 / -0.025 | 0.007 / -0.011 / -0.016 / 0.030 / 0.012 | -0.009 / -0.003 / -0.029 / 0.054 / 0.000 |
| *Enterobacter cloacae* | gentamicin | 104 | study | 0.022 / -0.049 / -0.061 / -0.018 / 0.000 | 0.000 / 0.000 / -0.020 / 0.056 / -0.026 |  | 0.022 / -0.070 / -0.043 / -0.018 / 0.000 |
| *Enterobacter cloacae* | imipenem | 86 | lineage | 0.009 / -0.016 / 0.031 / -0.072 / -0.006 | 0.000 / 0.048 / -0.058 / 0.001 / 0.081 | -0.013 / 0.062 / 0.092 / 0.014 / -0.015 | 0.010 / 0.038 / 0.028 / -0.072 / -0.003 |
| *Enterobacter cloacae* | imipenem | 86 | study | -0.167 / -0.009 / -0.035 / 0.016 / 0.018 | 0.000 / 0.000 / 0.026 / -0.027 / 0.121 |  | -0.167 / -0.009 / -0.035 / 0.033 / 0.018 |
| *Enterobacter cloacae* | levofloxacin | 104 | lineage | -0.001 / 0.034 / 0.003 / 0.009 / 0.000 | 0.000 / 0.000 / 0.119 / -0.001 / 0.102 | -0.036 / -0.022 / -0.046 / 0.009 / -0.092 | 0.032 / 0.007 / 0.003 / 0.009 / 0.000 |
| *Enterobacter cloacae* | levofloxacin | 104 | study | -0.004 / -0.084 / 0.011 / -0.089 / 0.000 | 0.000 / 0.000 / 0.006 / -0.124 / -0.013 |  | 0.002 / -0.084 / -0.006 / -0.089 / 0.000 |
| *Enterobacter cloacae* | meropenem | 99 | lineage | -0.013 / -0.017 / -0.004 / 0.011 / 0.000 | 0.000 / -0.006 / 0.011 / -0.096 / 0.003 | 0.003 / 0.020 / -0.001 / -0.014 / 0.057 | -0.013 / -0.013 / -0.007 / 0.001 / -0.003 |
| *Enterobacter cloacae* | meropenem | 99 | study | -0.029 / -0.059 / -0.017 / -0.029 / 0.000 | 0.000 / -0.039 / -0.079 / -0.006 / 0.011 |  | -0.029 / -0.064 / -0.017 / -0.029 / 0.008 |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | 80 | lineage | 0.000 / 0.000 / 0.004 / 0.008 / 0.000 | 0.000 / 0.000 / -0.038 / 0.108 / -0.004 | 0.000 / 0.000 / -0.008 / 0.033 / -0.008 | 0.000 / 0.000 / -0.013 / -0.073 / 0.000 |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | 80 | study | 0.000 / 0.000 / -0.023 / 0.000 / 0.000 | 0.000 / 0.000 / -0.013 / -0.024 / 0.000 |  | 0.000 / 0.000 / -0.023 / 0.000 / 0.000 |
| *Enterobacter cloacae* | tobramycin | 99 | lineage | 0.007 / -0.062 / -0.010 / -0.003 / 0.000 | 0.000 / -0.092 / 0.064 / 0.096 / 0.004 | 0.000 / -0.006 / 0.099 / 0.060 / 0.010 | 0.007 / -0.100 / 0.025 / 0.000 / 0.000 |
| *Enterobacter cloacae* | tobramycin | 99 | study | 0.044 / 0.028 / 0.009 / -0.022 / 0.000 | 0.000 / 0.000 / 0.046 / 0.000 / -0.019 |  | 0.044 / 0.028 / 0.009 / -0.022 / 0.000 |
| *Escherichia coli* | amikacin | 1,125 | lineage | 0.005 / 0.000 / 0.000 / -0.008 / 0.000 | 0.000 / 0.000 / 0.038 / -0.020 / -0.010 | 0.000 / 0.020 / 0.006 / -0.007 / 0.000 | 0.000 / 0.007 / 0.000 / -0.005 / 0.000 |
| *Escherichia coli* | amikacin | 1,125 | study | -0.023 / 0.000 / 0.034 / 0.034 / 0.000 | 0.000 / 0.000 / 0.150 / -0.039 / -0.075 |  | -0.056 / 0.020 / 0.050 / 0.034 / 0.000 |
| *Escherichia coli* | ampicillin/sulbactam | 283 | lineage | -0.037 / -0.033 / 0.075 / -0.019 / 0.000 | 0.000 / -0.039 / 0.126 / 0.024 / -0.101 | 0.018 / -0.012 / -0.027 / 0.059 / 0.012 | -0.040 / -0.030 / 0.066 / -0.026 / -0.015 |
| *Escherichia coli* | ampicillin/sulbactam | 283 | study | -0.032 / -0.007 / -0.003 / -0.007 / 0.000 | 0.000 / -0.011 / 0.099 / -0.124 / -0.050 |  | -0.032 / -0.007 / -0.003 / 0.004 / -0.007 |
| *Escherichia coli* | aztreonam | 699 | lineage | -0.013 / -0.006 / 0.163 / -0.002 / 0.000 | 0.000 / -0.023 / 0.176 / 0.025 / -0.010 | 0.044 / 0.007 / -0.006 / 0.002 / 0.000 | 0.004 / -0.006 / 0.164 / -0.006 / 0.000 |
| *Escherichia coli* | aztreonam | 699 | study | 0.019 / -0.014 / 0.019 / 0.002 / 0.000 | 0.000 / 0.000 / -0.003 / -0.002 / 0.003 |  | 0.019 / -0.024 / 0.019 / 0.010 / 0.000 |
| *Escherichia coli* | cefazolin | 414 | lineage | 0.037 / 0.016 / 0.019 / 0.025 / 0.000 | 0.000 / 0.000 / 0.071 / 0.006 / -0.099 | -0.038 / -0.017 / 0.007 / -0.007 / -0.005 | 0.048 / -0.052 / 0.034 / 0.040 / 0.002 |
| *Escherichia coli* | cefazolin | 414 | study | -0.046 / 0.043 / 0.106 / -0.017 / 0.000 | 0.000 / -0.002 / 0.147 / -0.007 / -0.005 |  | -0.049 / 0.043 / 0.109 / 0.004 / -0.011 |
| *Escherichia coli* | cefepime | 1,124 | lineage | 0.009 / -0.002 / 0.303 / 0.023 / 0.000 | 0.000 / -0.013 / 0.311 / -0.011 / 0.029 | 0.027 / -0.060 / 0.004 / 0.006 / 0.089 | 0.018 / -0.002 / 0.303 / 0.029 / -0.001 |
| *Escherichia coli* | cefepime | 1,124 | study | -0.022 / -0.026 / 0.056 / -0.020 / 0.000 | 0.000 / 0.002 / 0.017 / -0.010 / -0.002 |  | -0.024 / -0.031 / 0.056 / -0.020 / 0.000 |
| *Escherichia coli* | cefoxitin | 1,066 | lineage | 0.005 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / -0.003 / 0.000 / 0.005 / -0.001 | -0.002 / 0.006 / 0.000 / 0.001 / 0.001 | 0.009 / 0.011 / 0.004 / -0.004 / -0.002 |
| *Escherichia coli* | cefoxitin | 1,066 | study | 0.000 / -0.003 / 0.000 / 0.000 / 0.000 | 0.000 / 0.015 / 0.012 / 0.014 / 0.018 |  | 0.000 / -0.003 / 0.001 / -0.004 / 0.001 |
| *Escherichia coli* | ceftriaxone | 923 | lineage | 0.000 / 0.000 / -0.002 / -0.004 / 0.000 | 0.000 / 0.000 / -0.016 / -0.005 / -0.014 | 0.018 / 0.002 / -0.015 / -0.011 / -0.013 | 0.011 / -0.001 / -0.007 / -0.013 / 0.003 |
| *Escherichia coli* | ceftriaxone | 923 | study | -0.014 / 0.000 / -0.005 / -0.002 / 0.000 | 0.000 / -0.001 / -0.011 / -0.015 / -0.005 |  | -0.010 / 0.005 / -0.022 / 0.004 / 0.002 |
| *Escherichia coli* | ertapenem | 820 | lineage | -0.006 / -0.011 / -0.001 / 0.000 / 0.000 | 0.000 / 0.000 / -0.014 / -0.005 / 0.004 | 0.006 / 0.001 / 0.010 / 0.000 / 0.011 | -0.002 / -0.018 / -0.007 / -0.002 / 0.000 |
| *Escherichia coli* | ertapenem | 820 | study | 0.009 / 0.000 / 0.054 / 0.020 / 0.000 | 0.000 / -0.001 / -0.002 / -0.001 / -0.015 |  | 0.010 / -0.028 / -0.020 / -0.007 / 0.000 |
| *Escherichia coli* | levofloxacin | 403 | lineage | 0.000 / 0.000 / -0.051 / 0.013 / 0.000 | 0.000 / 0.000 / -0.093 / -0.009 / 0.010 | 0.011 / 0.009 / 0.001 / 0.010 / 0.011 | 0.000 / 0.000 / -0.079 / 0.013 / -0.001 |
| *Escherichia coli* | levofloxacin | 403 | study | 0.000 / 0.000 / 0.038 / 0.007 / 0.000 | 0.000 / 0.000 / 0.049 / 0.011 / 0.004 |  | 0.000 / -0.015 / 0.023 / 0.014 / 0.005 |
| *Escherichia coli* | nitrofurantoin | 482 | lineage | -0.019 / -0.012 / 0.158 / 0.027 / 0.000 | 0.000 / -0.053 / 0.210 / -0.091 / -0.001 | -0.021 / -0.017 / -0.007 / -0.022 / 0.013 | -0.020 / -0.007 / 0.158 / 0.030 / 0.024 |
| *Escherichia coli* | nitrofurantoin | 482 | study | -0.026 / -0.038 / 0.017 / -0.009 / 0.000 | 0.000 / 0.000 / -0.023 / -0.017 / -0.028 |  | -0.026 / -0.043 / 0.017 / 0.003 / -0.003 |
| *Escherichia coli* | piperacillin/tazobactam | 979 | lineage | -0.001 / -0.013 / 0.000 / 0.003 / 0.000 | 0.000 / 0.000 / -0.005 / -0.036 / 0.017 | 0.002 / -0.004 / 0.001 / -0.007 / -0.016 | -0.014 / -0.020 / -0.004 / 0.004 / -0.002 |
| *Escherichia coli* | piperacillin/tazobactam | 979 | study | 0.000 / -0.004 / 0.001 / 0.000 / 0.000 | 0.000 / 0.002 / -0.013 / -0.002 / 0.030 |  | -0.026 / -0.032 / -0.050 / 0.008 / 0.003 |
| *Escherichia coli* | streptomycin | 386 | lineage | 0.000 / 0.000 / -0.001 / -0.005 / 0.000 | 0.000 / 0.000 / 0.014 / -0.010 / -0.025 | 0.014 / 0.000 / 0.006 / -0.009 / -0.011 | -0.019 / -0.002 / 0.012 / -0.006 / 0.001 |
| *Escherichia coli* | streptomycin | 386 | study | 0.009 / 0.002 / 0.096 / 0.089 / 0.000 | 0.000 / 0.019 / 0.091 / 0.051 / -0.047 |  | 0.009 / -0.012 / 0.097 / 0.098 / -0.004 |
| *Escherichia coli* | sulfamethoxazole | 341 | lineage | 0.001 / 0.004 / -0.003 / 0.006 / 0.000 | 0.000 / 0.000 / 0.002 / 0.003 / -0.011 | 0.003 / 0.026 / -0.002 / 0.009 / -0.008 | -0.014 / 0.027 / -0.005 / 0.005 / 0.001 |
| *Escherichia coli* | sulfamethoxazole | 341 | study | 0.000 / -0.006 / 0.014 / 0.000 / 0.000 | 0.000 / 0.000 / 0.061 / -0.010 / 0.011 |  | -0.026 / -0.011 / 0.015 / 0.018 / 0.006 |
| *Escherichia coli* | tetracycline | 782 | lineage | 0.000 / 0.013 / 0.000 / 0.000 / 0.000 | 0.000 / 0.023 / 0.001 / 0.009 / 0.002 | -0.004 / -0.007 / -0.003 / 0.012 / -0.004 | 0.000 / 0.042 / 0.007 / 0.009 / 0.000 |
| *Escherichia coli* | tetracycline | 782 | study | 0.000 / -0.001 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.007 / -0.005 / -0.002 |  | 0.000 / -0.002 / 0.010 / 0.000 / 0.000 |
| *Escherichia coli* | tobramycin | 1,070 | lineage | 0.000 / 0.000 / -0.002 / 0.000 / 0.000 | 0.000 / -0.001 / 0.046 / -0.008 / -0.002 | -0.005 / 0.000 / -0.014 / -0.003 / -0.006 | 0.002 / 0.007 / 0.031 / -0.006 / 0.004 |
| *Escherichia coli* | tobramycin | 1,070 | study | 0.000 / 0.000 / -0.002 / -0.002 / 0.000 | 0.000 / 0.000 / -0.009 / -0.004 / -0.056 |  | 0.000 / 0.033 / 0.000 / -0.005 / 0.001 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 744 | lineage | -0.006 / 0.007 / 0.000 / 0.000 / 0.000 | 0.000 / -0.019 / 0.000 / 0.007 / 0.004 | 0.031 / -0.008 / -0.001 / 0.004 / -0.012 | 0.028 / -0.001 / -0.018 / 0.001 / -0.002 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 744 | study | 0.001 / 0.003 / 0.000 / 0.007 / 0.000 | 0.000 / 0.000 / 0.038 / -0.008 / 0.047 |  | 0.004 / -0.015 / 0.040 / 0.003 / -0.001 |
| *Klebsiella pneumoniae* | aztreonam | 1,050 | lineage | -0.006 / 0.002 / 0.000 / 0.019 / 0.000 | 0.000 / -0.017 / -0.002 / -0.015 / -0.171 | 0.021 / 0.013 / -0.005 / -0.010 / -0.048 | 0.006 / 0.017 / -0.002 / 0.001 / -0.024 |
| *Klebsiella pneumoniae* | aztreonam | 1,050 | study | -0.006 / -0.036 / 0.100 / 0.007 / 0.000 | 0.000 / 0.009 / 0.068 / 0.059 / 0.038 |  | -0.032 / -0.044 / 0.093 / 0.004 / 0.001 |
| *Klebsiella pneumoniae* | cefazolin | 1,141 | lineage | 0.000 / -0.002 / 0.000 / 0.000 / 0.000 | 0.000 / -0.026 / -0.004 / -0.001 / -0.011 | -0.004 / -0.009 / 0.001 / 0.014 / 0.032 | 0.001 / -0.001 / -0.003 / -0.008 / 0.001 |
| *Klebsiella pneumoniae* | cefazolin | 1,141 | study | 0.000 / -0.003 / 0.000 / 0.009 / 0.000 | 0.000 / 0.000 / 0.012 / -0.022 / -0.012 |  | 0.001 / 0.039 / 0.007 / 0.012 / 0.000 |
| *Klebsiella pneumoniae* | cefepime | 1,233 | lineage | 0.001 / -0.001 / 0.000 / 0.000 / 0.000 | 0.000 / -0.029 / 0.008 / 0.064 / -0.014 | -0.040 / -0.037 / 0.007 / -0.018 / -0.002 | 0.001 / 0.003 / 0.008 / -0.001 / -0.004 |
| *Klebsiella pneumoniae* | cefepime | 1,233 | study | -0.009 / -0.003 / 0.000 / 0.000 / 0.000 | 0.000 / 0.008 / 0.008 / -0.022 / -0.021 |  | -0.009 / -0.006 / -0.002 / -0.011 / -0.002 |
| *Klebsiella pneumoniae* | cefoxitin | 1,181 | lineage | 0.002 / -0.016 / 0.000 / -0.010 / 0.000 | 0.000 / -0.017 / -0.008 / -0.012 / -0.001 | -0.046 / -0.094 / 0.000 / -0.007 / -0.006 | -0.003 / -0.008 / 0.000 / -0.005 / -0.002 |
| *Klebsiella pneumoniae* | cefoxitin | 1,181 | study | 0.000 / -0.002 / 0.000 / 0.003 / 0.000 | 0.000 / 0.000 / 0.006 / 0.001 / -0.028 |  | 0.013 / -0.072 / 0.000 / 0.003 / -0.002 |
| *Klebsiella pneumoniae* | ceftazidime | 1,310 | lineage | 0.001 / -0.001 / 0.000 / -0.001 / 0.000 | 0.000 / 0.002 / -0.006 / -0.002 / 0.002 | 0.008 / -0.005 / 0.005 / 0.003 / 0.009 | 0.028 / 0.011 / -0.004 / -0.006 / -0.003 |
| *Klebsiella pneumoniae* | ceftazidime | 1,310 | study | -0.004 / 0.004 / 0.000 / 0.000 / 0.000 | 0.000 / -0.017 / -0.015 / -0.006 / -0.009 |  | -0.015 / 0.009 / -0.005 / 0.024 / 0.004 |
| *Klebsiella pneumoniae* | ceftriaxone | 1,208 | lineage | -0.002 / 0.001 / 0.000 / -0.001 / 0.000 | 0.000 / 0.001 / -0.012 / -0.013 / -0.001 | -0.001 / 0.000 / -0.005 / 0.002 / -0.002 | 0.001 / 0.003 / -0.005 / 0.003 / -0.002 |
| *Klebsiella pneumoniae* | ceftriaxone | 1,208 | study | 0.000 / -0.005 / 0.000 / 0.000 / 0.000 | 0.000 / -0.001 / 0.003 / 0.003 / -0.049 |  | -0.010 / 0.001 / 0.023 / 0.001 / 0.001 |
| *Klebsiella pneumoniae* | ciprofloxacin | 1,706 | lineage | -0.004 / -0.005 / 0.000 / 0.000 / 0.001 | 0.000 / 0.001 / 0.001 / -0.007 / -0.004 | 0.004 / -0.001 / -0.003 / 0.005 / -0.009 | -0.009 / -0.004 / 0.001 / -0.005 / 0.002 |
| *Klebsiella pneumoniae* | ciprofloxacin | 1,706 | study | 0.000 / -0.002 / 0.000 / 0.000 / 0.000 | 0.000 / 0.006 / -0.001 / -0.008 / 0.005 |  | -0.014 / 0.002 / 0.005 / -0.023 / -0.007 |
| *Klebsiella pneumoniae* | colistin | 336 | lineage | 0.000 / 0.000 / 0.000 / 0.001 / 0.000 | 0.000 / 0.000 / -0.066 / -0.001 / 0.000 | 0.002 / -0.002 / -0.038 / -0.004 / -0.001 | -0.004 / 0.002 / 0.000 / 0.001 / 0.000 |
| *Klebsiella pneumoniae* | colistin | 336 | study | 0.004 / 0.078 / 0.004 / 0.016 / 0.000 | 0.000 / -0.038 / 0.036 / 0.020 / -0.002 |  | 0.000 / 0.078 / -0.009 / 0.018 / -0.007 |
| *Klebsiella pneumoniae* | doripenem | 243 | lineage | -0.010 / 0.001 / -0.026 / -0.035 / 0.000 | 0.000 / 0.000 / -0.007 / -0.047 / -0.020 | 0.002 / -0.017 / 0.000 / -0.004 / 0.016 | -0.036 / 0.003 / -0.015 / -0.040 / 0.000 |
| *Klebsiella pneumoniae* | doripenem | 243 | study | 0.082 / 0.051 / -0.003 / -0.127 / 0.000 | 0.000 / 0.074 / 0.021 / -0.092 / 0.033 |  | 0.099 / 0.074 / -0.018 / -0.191 / 0.000 |
| *Klebsiella pneumoniae* | ertapenem | 952 | lineage | 0.000 / 0.028 / 0.003 / 0.002 / 0.000 | 0.000 / 0.001 / 0.001 / -0.007 / -0.006 | -0.004 / 0.080 / 0.017 / -0.003 / -0.013 | -0.004 / 0.031 / 0.010 / -0.009 / -0.001 |
| *Klebsiella pneumoniae* | ertapenem | 952 | study | 0.000 / -0.021 / -0.002 / 0.000 / 0.000 | 0.000 / 0.000 / 0.000 / -0.020 / -0.004 |  | -0.025 / -0.024 / -0.017 / 0.004 / 0.003 |
| *Klebsiella pneumoniae* | imipenem | 1,158 | lineage | -0.003 / -0.007 / 0.000 / -0.005 / 0.000 | 0.000 / 0.000 / -0.009 / -0.003 / -0.007 | 0.001 / -0.013 / 0.005 / 0.000 / -0.005 | -0.003 / -0.009 / -0.010 / -0.013 / -0.002 |
| *Klebsiella pneumoniae* | imipenem | 1,158 | study | 0.000 / -0.033 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.014 / -0.015 / -0.002 |  | -0.047 / -0.031 / -0.003 / -0.003 / -0.001 |
| *Klebsiella pneumoniae* | levofloxacin | 871 | lineage | -0.001 / 0.000 / -0.005 / -0.004 / 0.000 | 0.000 / 0.000 / 0.014 / -0.004 / -0.002 | 0.013 / -0.001 / -0.002 / 0.003 / -0.002 | -0.001 / 0.000 / -0.009 / -0.002 / 0.000 |
| *Klebsiella pneumoniae* | levofloxacin | 871 | study | -0.010 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / 0.004 / 0.001 / -0.014 |  | -0.010 / 0.000 / -0.001 / -0.002 / 0.000 |
| *Klebsiella pneumoniae* | nitrofurantoin | 562 | lineage | -0.029 / 0.000 / 0.014 / -0.008 / -0.001 | 0.000 / -0.072 / 0.017 / -0.031 / 0.037 | -0.014 / 0.030 / -0.005 / -0.048 / 0.066 | 0.012 / -0.001 / 0.006 / 0.017 / -0.002 |
| *Klebsiella pneumoniae* | nitrofurantoin | 562 | study | 0.208 / 0.025 / 0.083 / 0.015 / -0.001 | 0.000 / -0.029 / 0.072 / 0.017 / 0.012 |  | 0.205 / -0.027 / 0.088 / 0.121 / -0.012 |
| *Klebsiella pneumoniae* | polymyxin b | 100 | lineage | none / none / none / none / none | none / none / none / none / none | none / none / none / none / none | none / none / none / none / none |
| *Klebsiella pneumoniae* | polymyxin b | 100 | study | -0.008 / 0.001 / 0.000 / -0.097 / 0.000 | 0.000 / -0.045 / -0.032 / -0.117 / 0.000 |  | -0.008 / 0.022 / 0.000 / -0.097 / 0.000 |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 1,546 | lineage | -0.003 / -0.004 / 0.000 / 0.000 / 0.000 | 0.000 / -0.002 / -0.006 / -0.008 / -0.004 | 0.005 / 0.000 / -0.002 / -0.001 / 0.002 | -0.012 / -0.004 / -0.001 / -0.005 / 0.002 |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 1,546 | study | 0.000 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / -0.076 / 0.010 / -0.006 / -0.005 |  | 0.002 / 0.005 / 0.006 / 0.002 / -0.002 |
| *Klebsiella pneumoniae* | tigecycline | 924 | lineage | 0.028 / -0.008 / 0.039 / 0.004 / 0.000 | 0.000 / -0.004 / 0.067 / -0.007 / 0.011 | 0.011 / -0.020 / 0.018 / -0.025 / 0.002 | 0.041 / -0.018 / 0.056 / 0.001 / 0.004 |
| *Klebsiella pneumoniae* | tigecycline | 924 | study | -0.004 / -0.002 / -0.018 / -0.004 / -0.001 | 0.000 / -0.060 / 0.022 / 0.001 / -0.001 |  | 0.010 / 0.018 / 0.005 / 0.034 / -0.003 |
| *Klebsiella pneumoniae* | tobramycin | 1,230 | lineage | 0.006 / 0.000 / 0.000 / 0.000 / 0.000 | 0.000 / 0.000 / -0.025 / -0.036 / -0.013 | -0.003 / 0.013 / -0.001 / -0.009 / -0.006 | 0.007 / 0.000 / 0.003 / -0.003 / 0.000 |
| *Klebsiella pneumoniae* | tobramycin | 1,230 | study | 0.000 / -0.015 / 0.000 / 0.000 / 0.000 | 0.000 / -0.039 / -0.008 / -0.010 / -0.007 |  | 0.002 / -0.023 / 0.009 / -0.004 / -0.005 |
| *Mycobacterium tuberculosis* | nicotinamide | 166 | lineage | -0.048 / -0.010 / 0.002 / 0.082 / 0.000 | 0.000 / 0.000 / 0.171 / -0.038 / -0.032 | -0.014 / -0.012 / 0.002 / 0.169 / -0.010 | -0.048 / -0.020 / 0.002 / 0.124 / -0.004 |
| *Mycobacterium tuberculosis* | nicotinamide | 166 | study | 0.000 / 0.000 / -0.018 / -0.024 / 0.000 | 0.000 / 0.000 / -0.006 / -0.036 / 0.012 |  | 0.024 / -0.018 / -0.030 / -0.030 / 0.006 |
| *Pseudomonas aeruginosa* | aztreonam | 235 | lineage | 0.031 / 0.023 / -0.006 / -0.038 / -0.003 | 0.000 / -0.026 / -0.017 / -0.013 / -0.039 | 0.041 / -0.038 / -0.029 / -0.041 / -0.019 | 0.024 / -0.012 / 0.004 / -0.038 / -0.004 |
| *Pseudomonas aeruginosa* | aztreonam | 235 | study | 0.036 / -0.003 / 0.000 / -0.029 / 0.000 | 0.000 / 0.000 / 0.107 / -0.016 / -0.003 |  | 0.036 / 0.011 / 0.116 / -0.029 / 0.004 |
| *Pseudomonas aeruginosa* | cefepime | 245 | lineage | 0.027 / 0.042 / 0.005 / 0.003 / 0.000 | 0.000 / 0.009 / -0.019 / -0.045 / 0.008 | -0.041 / -0.035 / 0.026 / 0.022 / 0.014 | 0.006 / 0.028 / 0.019 / -0.024 / -0.012 |
| *Pseudomonas aeruginosa* | cefepime | 245 | study | 0.061 / 0.004 / 0.000 / -0.012 / 0.000 | 0.000 / 0.000 / 0.020 / -0.053 / 0.004 |  | 0.097 / 0.029 / 0.041 / -0.008 / 0.065 |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 309 | lineage | -0.039 / -0.023 / -0.002 / -0.026 / 0.000 | 0.000 / -0.010 / 0.025 / -0.044 / -0.086 | 0.023 / 0.076 / 0.002 / -0.008 / -0.018 | -0.042 / 0.001 / 0.006 / -0.028 / 0.007 |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 309 | study | 0.004 / -0.004 / -0.001 / 0.005 / 0.000 | 0.000 / -0.004 / -0.075 / -0.013 / -0.022 |  | 0.009 / -0.014 / -0.005 / 0.008 / 0.000 |
| *Pseudomonas aeruginosa* | gentamicin | 346 | lineage | -0.001 / -0.008 / 0.008 / -0.001 / 0.000 | 0.000 / 0.000 / 0.001 / 0.001 / -0.013 | 0.005 / 0.000 / -0.008 / 0.003 / 0.016 | -0.015 / -0.007 / 0.008 / -0.002 / -0.001 |
| *Pseudomonas aeruginosa* | gentamicin | 346 | study | 0.000 / 0.000 / 0.005 / 0.000 / 0.000 | 0.000 / 0.000 / 0.011 / 0.008 / -0.023 |  | 0.008 / -0.014 / 0.027 / 0.025 / 0.000 |
| *Pseudomonas aeruginosa* | imipenem | 447 | lineage | -0.015 / -0.008 / -0.010 / -0.007 / 0.004 | 0.000 / -0.034 / -0.018 / 0.000 / -0.037 | 0.005 / 0.003 / -0.005 / 0.019 / -0.008 | -0.023 / -0.011 / -0.009 / -0.003 / 0.006 |
| *Pseudomonas aeruginosa* | imipenem | 447 | study | -0.032 / -0.041 / -0.018 / -0.004 / 0.000 | 0.000 / 0.000 / 0.001 / 0.014 / -0.012 |  | -0.061 / -0.054 / -0.023 / -0.009 / 0.000 |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 210 | lineage | 0.001 / 0.028 / 0.002 / 0.001 / -0.004 | 0.000 / -0.023 / 0.012 / -0.051 / -0.105 | -0.010 / -0.003 / 0.010 / -0.076 / 0.011 | 0.001 / 0.029 / 0.004 / -0.013 / -0.008 |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 210 | study | 0.091 / -0.028 / -0.016 / 0.039 / 0.000 | 0.000 / -0.114 / -0.092 / -0.012 / 0.000 |  | 0.091 / -0.032 / -0.016 / 0.049 / -0.006 |
| *Salmonella enterica* | sulfisoxazole | 544 | lineage | -0.064 / 0.000 / -0.008 / -0.040 / 0.000 | 0.000 / 0.000 / -0.028 / -0.018 / -0.013 | -0.031 / -0.021 / 0.005 / -0.018 / 0.012 | -0.087 / -0.001 / -0.015 / -0.034 / -0.001 |
| *Salmonella enterica* | sulfisoxazole | 544 | study | -0.022 / 0.000 / -0.003 / -0.026 / 0.000 | 0.000 / 0.000 / -0.004 / -0.025 / -0.021 |  | -0.024 / -0.004 / -0.003 / -0.019 / -0.002 |
| *Staphylococcus aureus* | clindamycin | 1,543 | lineage | -0.094 / 0.023 / 0.004 / -0.024 / 0.000 | 0.000 / 0.088 / 0.053 / 0.078 / 0.067 | 0.010 / 0.005 / 0.003 / -0.035 / 0.000 | -0.061 / 0.048 / 0.005 / -0.024 / 0.011 |
| *Staphylococcus aureus* | clindamycin | 1,543 | study | -0.020 / -0.123 / -0.001 / 0.015 / 0.000 | 0.000 / 0.000 / 0.055 / -0.003 / 0.061 |  | -0.020 / -0.123 / -0.003 / 0.030 / 0.000 |
| *Staphylococcus aureus* | erythromycin | 1,957 | lineage | 0.014 / -0.003 / -0.012 / -0.086 / 0.000 | 0.000 / -0.019 / 0.073 / 0.006 / 0.025 | -0.032 / -0.006 / -0.004 / -0.025 / 0.003 | 0.024 / -0.037 / -0.012 / -0.087 / -0.011 |
| *Staphylococcus aureus* | erythromycin | 1,957 | study | 0.000 / 0.019 / -0.029 / 0.049 / 0.000 | 0.000 / 0.007 / -0.012 / 0.029 / 0.021 |  | -0.006 / 0.023 / -0.029 / 0.044 / 0.001 |
| *Staphylococcus aureus* | penicillin | 1,508 | lineage | 0.000 / 0.000 / 0.005 / 0.000 / 0.000 | 0.000 / 0.002 / 0.004 / -0.002 / -0.018 | 0.000 / 0.000 / -0.003 / 0.001 / 0.002 | 0.002 / -0.005 / 0.007 / 0.000 / -0.001 |
| *Staphylococcus aureus* | penicillin | 1,508 | study | 0.000 / 0.000 / 0.027 / 0.000 / 0.000 | 0.000 / -0.041 / 0.037 / -0.009 / 0.003 |  | -0.009 / -0.002 / 0.048 / 0.001 / 0.000 |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | 1,389 | lineage | 0.000 / 0.000 / 0.021 / 0.007 / 0.000 | 0.000 / 0.000 / -0.039 / -0.019 / -0.169 | 0.000 / -0.002 / -0.014 / -0.017 / 0.001 | 0.000 / 0.001 / 0.021 / -0.002 / 0.000 |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | 1,389 | study | 0.000 / 0.000 / -0.004 / 0.000 / 0.000 | 0.000 / 0.020 / -0.002 / 0.002 / 0.000 |  | 0.008 / 0.022 / -0.006 / 0.035 / -0.002 |
