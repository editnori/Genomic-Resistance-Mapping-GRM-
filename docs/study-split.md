# Study split

At evaluation revision 14, 70 of the release's 274 datasets, in 9 species, pass the first plan's rule for entering a split by study (the Folds section gives it). They hold 57,914 genome and dataset pairs, 28,277 of them resistant, in cohorts of 80 to 1,957 genomes (median 845.5).

Two experiments scored each pair with folds that keep each BioProject whole and with random folds on the same genomes. A BioProject stands in for a study, and an approximate one: the files cannot show whether two BioProjects hold one study, and 3,780 pairs have no BioProject.

The first experiment compares one assignment of BioProjects to five study folds with the release's own random folds and its lineage split. The second draws five study assignments and five random draws per dataset and adds a third arm: matched random folds, which copy each study plan's fold sizes and resistant counts.

Five learners fit each fold plan: scm (set covering machine), cart (one decision tree), rf (random forest), gbt (boosted trees) and lasso. scm scores only 0 or 1, so its log loss is 9.21 times its error rate plus a constant, and outside the first experiment's two tables this page reads scm's loss as its error rate.

Each experiment is a retrospective split by study inside one database, BV-BRC, and never external validation. A test genome's study mates are the other genomes of its BioProject in the dataset, and its lineage-cluster mates those of its cluster in the release's lineage split. The design does not separate study from lineage: study folds keep every study mate of a test genome out of training, and about half of its lineage-cluster mates as well.

## Plans

Each experiment hashed its plan before its first fit. `order.tsv` in the record lists every time on this page with its source, all on 2026-11-15 (UTC). The Rule section defines the contrasts, the species interval and a cost.

| Plan | Hashed | First result | Fixed |
|---|---|---|---|
| First | 11:44:32Z | 12:08:36Z | The eligibility rule, one study assignment, recall and log loss for every learner, the species interval and the rule for a cost |
| Second | 17:33:00Z | 17:45:31Z | Five draws of each arm, the three contrasts and their thresholds, scm's error rate, and each answer with one species left out |

The second experiment's fits started at 17:44:59Z. The last result came at 13:03:17Z in the first experiment and at 19:49:03Z in the second.

**Exploratory.** Every reading the plans did not fix is exploratory: the first experiment's fold cuts, the lineage mates, the distances between assignments, the draw split, the two largest species left out together, the other outcome rules and thresholds, and the shares of the drop. The first plan gave its lineage contrasts no rule, and the second lists the genomes without a BioProject as a description with no rule.

## Folds

The fold builder places groups of genomes. A copy group is a set of genomes the release marks as identical copies. A group joins a BioProject with every copy group that holds one of its genomes, and so with any other BioProject such a copy group reaches; a genome with no BioProject is a group of its own unless a copy group joins it.

A dataset enters when it holds at least 5 BioProjects, at least 0.8 of its genomes carry one, and no group holds over 0.5 of its genomes. Of the release's 274 datasets, 71 pass; 77 fail on BioProjects, 103 on the share and 104 on the largest group, some on more than one.

**Variant count.** The first plan gives 36 variants of its eligibility rule, counted, in the plan's words, before it fixed its numbers. The file times place the count at 11:42:13Z, after the fold builder (11:38:56Z), which already holds 5, 0.8 and 0.5, and after the fold files (11:40:30Z to 11:40:33Z), so the numbers and the folds came first.

The variants admit 29 to 153 datasets in 6 to 13 species. The count reads no score: the script rebuilds all 36 from each dataset's three eligibility counts, and all of it comes before the hash.

One dataset passes and did not enter, mycobacterium_tuberculosis__amoxicillin (57 genomes): its fold check gives "fold 4 test genomes lack a class".

That leaves 70 datasets in 9 species: 23 of *K. pneumoniae*, 19 of *E. coli*, 8 of *E. cloacae*, 6 of *A. baumannii*, 6 of *P. aeruginosa*, 5 of *S. aureus*, and 1 each of *E. faecium*, *M. tuberculosis* and *S. enterica*.

The second experiment holds 1,050 fold plans, five of each kind per dataset. Random draw 0 is the release's own random folds, and each matched plan copies its study draw's five fold sizes and resistant counts; the script checks the first in all 70 datasets and the second in all 350 matched plans.

Two placement rules assign the groups of a study plan to its five folds. Rule A, the first experiment's placement, takes groups largest first, with a keyed order only among groups of equal size; rule B takes every group in a keyed random order. Both then put each group in the fold that least raises the class imbalance summed over folds. Rule A places study draw 0, the first experiment's assignment, and rule B draws 1 to 4.

The movement figure is the share of a cohort that a rule moves from draw 0, from the plan builder's table of draws 1 to 4: a median 0.047 under rule A (7 datasets not at all) and 0.194 under rule B (0.005 to 0.732). The table's file time is 17:23:55.01Z.

**Rule time.** The second plan gives 17:23:34Z as the time it fixed, before measuring rule A, that rule A moving a median under 0.10 of a cohort from draw 0 is too little. No file time or log line records that time; only the plan's text gives it. The plan builder's file time, 17:23:24.89Z, comes before it, and the second experiment's review found rule B already in the builder.

Neither placement rule nor the movement figure reads a score, but the first experiment's results, complete at 13:03:17Z, were on disk when the rule was chosen.

`checks.tsv` lists the script's checks: every plan keeps copy groups whole (1,050 plans), every study plan keeps BioProjects whole (350), and every fold holds both classes in its test and training genomes (1,050 fold files). Its rows that begin `inputs_` give the checks run when the reduced input set was written.

Every arm selects settings inside each training set. Tuning keeps copy groups whole and keeps no study whole in any arm; the release's lineage split, the first experiment's third arm, tunes on lineage groups.

Study folds differ from random folds in size and in class balance. The table gives medians over fold plans with their range. Largest is the largest test fold's share of the cohort and Training the smallest training share; Spread is the range of the resistant share over a plan's five test folds, and Departure the largest gap between a test fold's resistant share and the cohort's. Matched plans equal the study draws on every measure, as they copy them.

| Measure | Release random | Random draws | First study | Study draws |
|---|---|---|---|---|
| Largest | 0.201 (0.200 to 0.212) | 0.201 (0.200 to 0.212) | 0.334 (0.205 to 0.503) | 0.350 (0.205 to 0.552) |
| Training | 0.799 (0.787 to 0.800) | 0.799 (0.787 to 0.800) | 0.666 (0.497 to 0.795) | 0.650 (0.448 to 0.795) |
| Spread | 0.006 (0 to 0.033) | 0.006 (0 to 0.033) | 0.326 (0 to 0.699) | 0.375 (0 to 0.787) |
| Departure | 0.004 (0 to 0.021) | 0.004 (0 to 0.021) | 0.195 (0 to 0.491) | 0.233 (0 to 0.637) |

No test genome in a study plan has a study mate in training. In random and matched plans 0.930 and 0.928 of test genomes have one (median over datasets of the mean over draws), with a median 98.5 to 100.5 and 90 to 95.25 study mates each over the draws.

Exploratory: study plans also keep more lineage-cluster mates out of training. The clusters are the ones the release's lineage split holds whole; the script checks this in 69 datasets. The median count of a test genome's lineage mates in training is 15.5 to 16.5 over the study draws, against 30.5 to 33 in matched and 33 to 34 in random plans.

The share of test genomes with any lineage mate in training is 0.861 in study plans, against 0.910 in matched and 0.914 in random plans. Matched minus study is positive in 66 of 70 datasets (median +0.042), and random minus matched in 54 (median +0.004).

Exploratory: the five study assignments share most of their placement. Between two study draws a median 0.206 of a cohort sits in another fold (0.006 to 0.728, after pairing fold labels for the largest agreement); the share is under 0.10 in 14 of 70 datasets, and in 1 dataset two study draws are identical. Between two random draws the median is 0.767, and between two matched draws 0.709.

So the second experiment's 25 study-random draw pairs, each one study draw beside one random draw, cross five close study assignments with five random draws, all on the same 70 datasets in 9 species.

## Rule

Each contrast is a paired difference per dataset on the same genomes, read as its median over the 70 datasets with a 95% species interval. The script computes each interval exactly, over all 24,310 multisets of the 9 species drawn with replacement, each at the weight of its multinomial count. The first plan fixed 10,000 resamples, and the script uses the exact count, as the first experiment's review did.

A median reads worse when its interval lies wholly on the worse side of zero, better when it lies wholly on the other side, and not shown otherwise; a bound at exactly 0 reads not shown. Worse is lower recall, specificity or balanced accuracy, or a higher loss.

For one learner a contrast gives a cost when recall or the loss reads worse and neither reads better. The first plan's answer for the experiment is cost shown when at least 3 of the 5 learners show a cost and none reads better. The second plan asks, per learner, for a cost in at least 20 of 25 study-random draw pairs, 4 of 5 study-matched draw pairs and 20 of 25 matched-random draw pairs.

The interval holds species resampling only. It holds nothing of the study assignment or the random seed, and nothing of other species or datasets. With 42 of the 70 datasets in two species and 3 species with one dataset each, a reading beyond the interval shows only that the median keeps its sign over resampled sets of these 9 species.

## First experiment

On the 57,914 pairs, five folds that keep each BioProject whole, under one assignment, score lower than the release's random folds for all 5 learners. Log loss and specificity read worse for all 5; recall reads worse for rf and gbt and is lower, but not shown, for scm, cart and lasso. The plan's answer is cost shown, with a cost for 5 of 5 learners, and it stands with any one of the 9 species left out.

| Learner | Recall | Specificity | Log loss |
|---|---|---|---|
| scm | -0.0181 [-0.0230, 0] not shown | -0.0392 [-0.0481, -0.0195] worse | +0.2841 [0.2353, 0.3446] worse |
| cart | -0.0437 [-0.0679, 0.0034] not shown | -0.0475 [-0.0678, -0.0217] worse | +0.2495 [0.1319, 0.3765] worse |
| rf | -0.0381 [-0.0477, -0.0158] worse | -0.0471 [-0.0622, -0.0244] worse | +0.1274 [0.1033, 0.1336] worse |
| gbt | -0.0351 [-0.0364, -0.0227] worse | -0.0500 [-0.0647, -0.0183] worse | +0.1251 [0.0756, 0.1449] worse |
| lasso | -0.0332 [-0.0477, 0] not shown | -0.0586 [-0.0715, -0.0204] worse | +0.1120 [0.0941, 0.1344] worse |

scm's log loss rise is a rise of +0.0308 [0.0255, 0.0374] in its error rate. With one species left out, lasso's recall reading changes without *E. cloacae*, *E. faecium* or *P. aeruginosa*, and no other reading of the table changes.

Study folds also differ from random folds in size and in class balance between training and test. The drop remains where both are small: on the 240 of 350 study folds that train on at least as many genomes as the random folds (33,105 pairs), and on the 88 of those whose resistant share lies within 0.054 of their training genomes' (42 datasets, 14,282 pairs).

On the 88 folds the specificity drop and the log loss rise are smaller than on the 240, and the recall drop is larger, for every learner.

A larger class shift is in line with a larger drop in recall and specificity and a larger log loss rise, for every learner. The lowest and highest thirds of study folds by shift train on similar counts (median ratio to the random folds 1.017 and 1.028). The first experiment's review chose these cuts after reading the results; they carry no rule and no interval and rest on one assignment.

Each row pools one group of study folds and gives the study folds' figure minus the random folds' on the same genomes.

| Learner | Folds | Recall | Specificity | Log loss |
|---|---|---|---|---|
| scm | 240 with larger training | -0.0254 | -0.0365 | +0.2854 |
| scm | 88 with both | -0.0320 | -0.0217 | +0.2360 |
| scm | 115 in the lowest third | -0.0267 | -0.0284 | +0.2556 |
| scm | 117 in the highest third | -0.0594 | -0.0542 | +0.5253 |
| cart | 240 with larger training | -0.0477 | -0.0469 | +0.1833 |
| cart | 88 with both | -0.0481 | -0.0351 | +0.1533 |
| cart | 115 in the lowest third | -0.0446 | -0.0449 | +0.2023 |
| cart | 117 in the highest third | -0.0810 | -0.0549 | +0.3432 |
| rf | 240 with larger training | -0.0466 | -0.0338 | +0.0929 |
| rf | 88 with both | -0.0599 | -0.0196 | +0.0746 |
| rf | 115 in the lowest third | -0.0534 | -0.0307 | +0.0867 |
| rf | 117 in the highest third | -0.0867 | -0.0913 | +0.2242 |
| gbt | 240 with larger training | -0.0363 | -0.0480 | +0.0952 |
| gbt | 88 with both | -0.0454 | -0.0284 | +0.0817 |
| gbt | 115 in the lowest third | -0.0383 | -0.0400 | +0.0924 |
| gbt | 117 in the highest third | -0.0932 | -0.0821 | +0.1949 |
| lasso | 240 with larger training | -0.0402 | -0.0431 | +0.0846 |
| lasso | 88 with both | -0.0531 | -0.0240 | +0.0631 |
| lasso | 115 in the lowest third | -0.0429 | -0.0305 | +0.0835 |
| lasso | 117 in the highest third | -0.0880 | -0.0945 | +0.1994 |

scm's four log loss rises are rises of +0.0310, +0.0256, +0.0278 and +0.0570 in its error rate.

The design does not measure how the drop divides between study mates, lineage mates and the class shift.

The first plan gave its lineage contrasts no rule. On the 69 datasets with a lineage split, study minus lineage is not shown on log loss for 4 of 5 learners.

Study folds keep more recall than the lineage split for gbt and lasso and lose more specificity for scm, rf and lasso, with every species left out. scm's recall, cart's recall and log loss, rf's recall, gbt's specificity and lasso's log loss readings each change with one species left out.

Each cell gives the median of study minus lineage, its reading and the species leave-outs, of 9, that keep the reading.

| Learner | Recall | Specificity | Log loss |
|---|---|---|---|
| scm | +0.0113 better, 3 of 9 | -0.0189 worse, 9 of 9 | +0.0282 not shown, 9 of 9 |
| cart | +0.0143 not shown, 8 of 9 | -0.0197 not shown, 9 of 9 | +0.0710 worse, 5 of 9 |
| rf | +0.0381 not shown, 7 of 9 | -0.0536 worse, 9 of 9 | +0.0051 not shown, 9 of 9 |
| gbt | +0.0517 better, 9 of 9 | -0.0357 worse, 7 of 9 | +0.0146 not shown, 9 of 9 |
| lasso | +0.0379 better, 9 of 9 | -0.0373 worse, 9 of 9 | -0.0148 not shown, 8 of 9 |

## Second experiment

Study folds score lower than random folds on the same genomes for all 5 learners, in each of the 25 draw pairs. The loss reads worse in 25 of 25 draw pairs for every learner. Draw pair medians of scm's error rate rise by +0.0274 to +0.0394, and of log loss by +0.1891 to +0.3201 for cart, +0.1148 to +0.1363 for rf, +0.0975 to +0.1263 for gbt and +0.0825 to +0.1352 for lasso.

Recall is lower in every draw pair as a median, by 0.0034 to 0.0337 for scm, 0.0274 to 0.0570 for cart, 0.0254 to 0.0460 for rf, 0.0220 to 0.0440 for gbt and 0.0189 to 0.0355 for lasso.

Recall reads worse beyond the species interval in 25 of 25 draw pairs for rf and gbt, 19 for cart, 15 for lasso and 6 for scm. Specificity reads worse in 25 of 25 draw pairs for scm and rf, 23 for gbt and lasso and 22 for cart.

The table gives each contrast as the mean over draws: its median over datasets with the species interval and reading, beside the draw pairs in which recall and the loss read worse. Loss is the error rate for scm and log loss for the others.

| Learner | Contrast | Recall | Specificity | Loss | Draw pairs |
|---|---|---|---|---|---|
| scm | study-random | -0.0158 [-0.0362, 0.0003] not shown | -0.0303 [-0.0744, -0.0245] worse | +0.0366 [0.0280, 0.0428] worse | recall 6 of 25, loss 25 of 25 |
| scm | study-matched | -0.0141 [-0.0299, -0.0032] worse | -0.0236 [-0.0497, -0.0164] worse | +0.0337 [0.0208, 0.0462] worse | recall 0 of 5, loss 5 of 5 |
| scm | matched-random | -0.0024 [-0.0050, 0.0019] not shown | -0.0022 [-0.0044, 0] not shown | +0.0024 [-0.0007, 0.0056] not shown | recall 0 of 25, loss 3 of 25 |
| cart | study-random | -0.0397 [-0.0541, -0.0171] worse | -0.0489 [-0.0561, -0.0240] worse | +0.2894 [0.0981, 0.4018] worse | recall 19 of 25, loss 25 of 25 |
| cart | study-matched | -0.0387 [-0.0462, -0.0230] worse | -0.0383 [-0.0682, -0.0161] worse | +0.2482 [0.1188, 0.3841] worse | recall 5 of 5, loss 5 of 5 |
| cart | matched-random | -0.0026 [-0.0067, 0.0006] not shown | -0.0020 [-0.0080, 0.0037] not shown | +0.0280 [0.0115, 0.0430] worse | recall 1 of 25, loss 9 of 25 |
| rf | study-random | -0.0373 [-0.0556, -0.0175] worse | -0.0425 [-0.0574, -0.0188] worse | +0.1207 [0.1029, 0.1425] worse | recall 25 of 25, loss 25 of 25 |
| rf | study-matched | -0.0335 [-0.0467, -0.0022] worse | -0.0298 [-0.0493, -0.0029] worse | +0.0948 [0.0638, 0.1035] worse | recall 3 of 5, loss 5 of 5 |
| rf | matched-random | -0.0052 [-0.0115, -0.0018] worse | -0.0063 [-0.0087, -0.0016] worse | +0.0200 [0.0157, 0.0364] worse | recall 7 of 25, loss 25 of 25 |
| gbt | study-random | -0.0348 [-0.0475, -0.0211] worse | -0.0432 [-0.0802, -0.0219] worse | +0.1087 [0.0664, 0.1292] worse | recall 25 of 25, loss 25 of 25 |
| gbt | study-matched | -0.0293 [-0.0433, -0.0174] worse | -0.0403 [-0.0667, -0.0040] worse | +0.0917 [0.0558, 0.1271] worse | recall 3 of 5, loss 5 of 5 |
| gbt | matched-random | -0.0022 [-0.0089, -0.0014] worse | -0.0050 [-0.0130, -0.0029] worse | +0.0090 [0.0049, 0.0228] worse | recall 7 of 25, loss 22 of 25 |
| lasso | study-random | -0.0257 [-0.0525, -0.0073] worse | -0.0505 [-0.0726, -0.0264] worse | +0.1078 [0.0707, 0.1382] worse | recall 15 of 25, loss 25 of 25 |
| lasso | study-matched | -0.0226 [-0.0491, -0.0049] worse | -0.0448 [-0.0688, -0.0205] worse | +0.0966 [0.0715, 0.1047] worse | recall 3 of 5, loss 5 of 5 |
| lasso | matched-random | -0.0005 [-0.0023, 0.0007] not shown | -0.0066 [-0.0167, -0.0026] worse | +0.0095 [0.0038, 0.0137] worse | recall 1 of 25, loss 17 of 25 |

Under the plan's rule the answers are as below. The first experiment's rule gives cost shown in 25 of 25 study-random draw pairs.

| Learner | Study-random | Study-matched | Matched-random |
|---|---|---|---|
| scm | 25 of 25, shown | 5 of 5, shown | 3 of 25, not shown |
| cart | 25 of 25, shown | 5 of 5, shown | 10 of 25, not shown |
| rf | 25 of 25, shown | 5 of 5, shown | 25 of 25, shown |
| gbt | 25 of 25, shown | 5 of 5, shown | 22 of 25, shown |
| lasso | 25 of 25, shown | 5 of 5, shown | 17 of 25, not shown |

Matched folds score a little lower than random folds. Log loss rises by +0.0280 for cart, +0.0200 for rf, +0.0090 for gbt and +0.0095 for lasso. Recall falls by 0.0052 for rf and 0.0022 for gbt, with no difference shown for scm, cart and lasso. Specificity reads worse for rf, gbt and lasso. scm's error rate shows no difference beyond the species interval.

Study folds score lower than matched folds for every learner on the loss: log loss +0.2482 for cart, +0.0948 for rf, +0.0917 for gbt and +0.0966 for lasso, and scm's error rate +0.0337. Recall reads worse in 5 of 5 draw pairs for cart, 3 for rf, gbt and lasso and 0 for scm. Specificity, over the mean of draws, reads worse for all 5 learners.

Exploratory: per dataset, study minus random equals study minus matched plus matched minus random, exactly. On average over datasets the part not matched by fold shape is 0.752 to 0.943 of the mean drop, per learner and on recall and the loss; on specificity it is 0.693 to 0.871. The table gives each share and the datasets, of 70, where that part is the larger.

| Learner | Recall | Specificity | Loss |
|---|---|---|---|
| scm | 0.912, 51 datasets | 0.829, 51 datasets | 0.907, 62 datasets |
| cart | 0.916, 60 datasets | 0.871, 52 datasets | 0.882, 60 datasets |
| rf | 0.841, 53 datasets | 0.693, 45 datasets | 0.752, 59 datasets |
| gbt | 0.855, 53 datasets | 0.780, 55 datasets | 0.845, 60 datasets |
| lasso | 0.943, 55 datasets | 0.793, 54 datasets | 0.883, 62 datasets |

The part not matched by fold shape comes with test genomes losing genomes of their own study, and about half of the genomes of their own lineage cluster, from training. The design does not separate study from lineage.

The part matched by fold shape, the smaller one, shows a cost for rf in 25 of 25 matched-random draw pairs. gbt shows one in 22 of 25, which holds only with the four more uneven study assignments and not without the two largest species; both are exploratory readings. lasso's 17, cart's 10 and scm's 3 stay under 20.

## Choices

Every reading under Choices is exploratory apart from the one-species leave-outs, which the second plan fixed.

**Thresholds.** The study-random answer holds for every learner at any threshold from 13 to 25 of 25 draw pairs, and the study-matched answer at 3 to 5 of 5.

The matched-random answer holds for rf, gbt and lasso at 13 to 17, for rf and gbt at 18 to 22 and for rf alone at 23 to 25. The plan fixed 20 of 25 before any fit, and gbt's 22 and lasso's 17 are the counts nearest to it.

**Draw split.** Study draw 0 is the first experiment's assignment, and matched draw 0 copies it; draws 1 to 4 use rule B, which gives more uneven folds. Each cell gives the draw pairs reading worse with draw 0 of the study or matched arm, then with draws 1 to 4.

| Learner | Contrast | Recall | Specificity | Loss | Cost |
|---|---|---|---|---|---|
| scm | study-random | 0 of 5, 6 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 |
| scm | matched-random | 0 of 5, 0 of 20 | 0 of 5, 0 of 20 | 1 of 5, 2 of 20 | 1 of 5, 2 of 20 |
| cart | study-random | 1 of 5, 18 of 20 | 5 of 5, 17 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 |
| cart | matched-random | 0 of 5, 1 of 20 | 0 of 5, 0 of 20 | 1 of 5, 8 of 20 | 1 of 5, 9 of 20 |
| rf | study-random | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 |
| rf | matched-random | 0 of 5, 7 of 20 | 4 of 5, 14 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 |
| gbt | study-random | 5 of 5, 20 of 20 | 4 of 5, 19 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 |
| gbt | matched-random | 0 of 5, 7 of 20 | 0 of 5, 4 of 20 | 3 of 5, 19 of 20 | 3 of 5, 19 of 20 |
| lasso | study-random | 4 of 5, 11 of 20 | 5 of 5, 18 of 20 | 5 of 5, 20 of 20 | 5 of 5, 20 of 20 |
| lasso | matched-random | 0 of 5, 1 of 20 | 1 of 5, 13 of 20 | 3 of 5, 14 of 20 | 3 of 5, 14 of 20 |

With rule B draws the study-random loss medians per draw pair overlap those of draw 0. With draw 0 against draws 1 to 4, scm's error rate is +0.0302 to +0.0310 against +0.0274 to +0.0394.

Log loss with draw 0 against draws 1 to 4 is +0.2386 to +0.2949 against +0.1891 to +0.3201 for cart, +0.1269 to +0.1359 against +0.1148 to +0.1363 for rf, +0.1163 to +0.1263 against +0.0975 to +0.1251 for gbt and +0.1022 to +0.1127 against +0.0825 to +0.1352 for lasso.

With rule B draws the matched-random cost shows more often for cart, gbt and lasso and less often for scm, and study-random recall reads worse more often for scm and cart and less often for lasso.

**Two species.** *E. coli* and *K. pneumoniae* hold 42 of the 70 datasets. With both left out (28 datasets in 7 species), gbt's matched-random cost falls to 14 of 25 and is not shown, and rf's to 21, still shown. Every other answer stands, study-random at 22 to 25 of 25 and study-matched at 5 of 5.

With one species left out, 2 answers change: gbt's matched-random answer without *E. coli* (19 of 25, not shown) and lasso's matched-random answer without *E. cloacae* (20 of 25, shown).

**Both measures.** Under a stricter rule, a cost only when recall and the loss both read worse, study-matched holds for cart alone (5 of 5; rf, gbt and lasso 3, scm 0) and study-random for rf and gbt (25 of 25; cart 19, lasso 15, scm 6). Recall alone gives the same counts for both contrasts.

With balanced accuracy in place of recall, lasso's matched-random count is 21 of 25, shown, against 17 under the plan's rule.

**No BioProject.** Exploratory, with no rule in the plan. Genomes without a BioProject have no study mate to lose in any arm. On the 39 datasets with at least 10 of them (median 77 genomes, 27 resistant; 6 with under 5 resistant), study minus matched is closer to zero on those genomes than on genomes with a BioProject, for every learner and measure.

These genomes are a different and smaller set, and the comparison does not separate study from lineage.

Each cell gives the mean over datasets of study minus matched, over the mean of draws, on genomes without a BioProject against genomes with one, in the same datasets: 36 for recall and 38 for specificity, where both classes stand among the genomes without a BioProject, and 39 for the loss.

| Learner | Recall | Specificity | Loss |
|---|---|---|---|
| scm | -0.0012 against -0.0413 | -0.0138 against -0.0432 | +0.0071 against +0.0368 |
| cart | -0.0023 against -0.0709 | -0.0109 against -0.0581 | -0.0144 against +0.2706 |
| rf | +0.0040 against -0.0525 | +0.0257 against -0.0691 | -0.0195 against +0.1136 |
| gbt | +0.0012 against -0.0541 | +0.0056 against -0.0616 | -0.0026 against +0.1150 |
| lasso | -0.0013 against -0.0574 | +0.0183 against -0.0558 | +0.0049 against +0.1054 |

## Limits

- Both experiments use one database, BV-BRC, with one curation of its labels and BioProjects.
- A BioProject is an approximate study, and no file gives the laboratory, place or year behind one.
- The matched arm copies fold sizes and class counts only, and study folds also keep about half of the lineage-cluster mates out of training.
- The five study assignments share most of their placement, and four of them use a more uneven rule than the first.
- The species interval rests on 9 species, with 42 of the 70 datasets in two of them and 3 species with one dataset each.
- Recall and specificity rest on one cutoff: a score of 0.5 or more calls resistant.
- scm's figures are its error rate, apart from the log loss in the first experiment's two tables, which the first plan fixed; both give the same readings.
- Study folds differ from random folds in fold size and in class balance between training and test, and the first experiment has no matched arm to hold both fixed.
- The page gives no study-split figure for a single dataset.
- The reviews chose the first experiment's fold cuts, the lineage mates, the draw split and the two largest species left out together after reading the results.
- The first experiment rests on one study assignment and one random seed, and the second on five close assignments and five random draws.
- Nothing on disk shows what the first experiment opened before its plan, so its statement that it read no score of the release reports rests on its word.
- No file time or log line records when the second plan fixed its line for rule A.
- The files cannot show whether genomes without a BioProject have study mates.

## Files

- `reference/study-split/` holds the record: 30 files that `scripts/study-split/record.py` writes from the reduced inputs. Its README lists each file and what guards it, and `SHA256SUMS` gives each file's sha256.
- The record holds no genome identifier, label, prediction per genome, model, feature or importance, and `checks.tsv` gives 0 as the count of the inputs' genome identifiers found in its files. Each fold plan enters by its sha256 and fold counts alone.
- The two plans stay in the reduced inputs. `planfiles.tsv` gives each plan's sha256 and the 151 and 1,063 files the two hash lists hold, and `--check` hashes the plans again.
- `tests/study_split` rebuilds from the record every paragraph, list item and table of this page that carries a number, in figures or words, and holds its other 11 paragraphs and list items word for word; [fold-file.md](fold-file.md) gives the fold file format.

```sh
python3 -I scripts/study-split/record.py --inputs <dir> --check reference/study-split
python3 -I -m unittest discover -s tests/study_split -t .
```

`<dir>` is the reduced input set, the directory that `study-split-inputs.tar.zst` unpacks to; the repository does not hold that archive. The record's `inputs.sha256` lists the sha256 of each of its files.
