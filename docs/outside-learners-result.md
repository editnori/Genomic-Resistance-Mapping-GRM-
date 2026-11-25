# Outside-learner result

The release fits GRM's own learners. This page asks whether learners from outside GRM, fitted on the same folds, call the same held-out genomes better: XGBoost 3.0.5 and LightGBM 4.6.0 on the 4,000 rows `grm export` screens, each read against GRM's boosted trees (`gbt`), and a scikit-learn 1.7.2 elastic net on unitigs, read against GRM's lasso (`lasso`), beside a twin. It reads block C7 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/outside-learners.md`.

The design was fixed at 07:37Z on 23 Nov 2026, the lead's rulings on its questions at 09:05Z and on the elastic net's stop at 11:26Z, before any C7 score was read (`docs/assumption-plan-corrections.md`, entry "C7 design"). The table script, its gates and this page's form were committed before any score of the run was opened.

The run was read once, from 14:24:49Z to 14:25:18Z on 25 Nov 2026, with no other option, subset or measure: its checks, the gathering, the count, the tables and their printout, each once. The section Method gives each command with its times.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

"No difference shown" is the rule's reading and no statement that two learners are equal. Recall stands beside log loss wherever log loss is given.

## Reading

**Result.** `lightgbm` and `elasticnet` meet conditions 1 to 5. `xgboost` meets 1 to 4 and not 5: its study median has the other sign. Condition 6 is the maintainer's for all three, so nothing changes before the maintainer rules.

**`xgboost`.** `xgboost-gbt` is better on 58 of 232 lineage datasets, worse on 17 and the same on 157. Its median difference in balanced accuracy is 0.003, with a species interval of 0.001 to 0.005 (condition 1 met). The one-sided sign test gives p < 0.001, and < 0.001 after Holm (condition 2 met). Its recall interval is 0.000 to 0.005 (condition 4 met). Its study median is -0.002 on 55 of 55 study datasets, against the lineage median's 0.003: the signs differ (condition 5 not met).

**`lightgbm`.** `lightgbm-gbt` is better on 62 of 232 lineage datasets, worse on 14 and the same on 156. Its median difference in balanced accuracy is 0.003, with a species interval of 0.001 to 0.004 (condition 1 met). The one-sided sign test gives p < 0.001, and < 0.001 after Holm (condition 2 met). Its recall interval is 0.000 to 0.005 (condition 4 met). Its study median is 0.001 on 55 of 55 study datasets, against the lineage median's 0.003: the signs agree (condition 5 met).

**`elasticnet`.** `elasticnet-lasso` is better on 71 of 232 lineage datasets, worse on 19 and the same on 142. Its median difference in balanced accuracy is 0.006, with a species interval of 0.003 to 0.009 (condition 1 met). The one-sided sign test gives p < 0.001, and < 0.001 after Holm (condition 2 met). Its recall interval is 0.000 to 0.007 (condition 4 met). Its study median is 0.004 on 55 of 55 study datasets, against the lineage median's 0.006: the signs agree (condition 5 met).

**Twin.** `twin-gbt` is better on 40 of 232, worse on 33 and the same on 159, with a median of 0.000, an interval of 0.000 to 0.001 and a Holm p of 0.243. `twin-lasso` is better on 28 of 232, worse on 19 and the same on 185, with a median of 0.000, an interval of -0.001 to 0.000 and a Holm p of 0.243. The twin meets condition 1 in neither pair and 2 in neither, so condition 3 is met for every arm.

Condition by condition, as `page-conditions.tsv` gives it:

| Pair | Condition | Role | Reads | Read |
|---|---|---|---|---|
| `xgboost-gbt` | 1 | decides | the lineage median difference's interval lies wholly above 0 | met |
| `xgboost-gbt` | 2 | decides | better beats worse by the one-sided sign test at p < 0.05 after Holm over the 5 lineage pairs | met |
| `xgboost-gbt` | 3 | decides | the twin meets 1 in neither of its pairs (gbt, lasso) and 2 in neither, each read apart | met |
| `xgboost-gbt` | 3 | beside, decides nothing | beside, deciding nothing: the twin meets 1 and 2 in neither on gbt, this arm's learner | met |
| `xgboost-gbt` | 4 | decides | the recall interval does not lie wholly below -0.02 | met |
| `xgboost-gbt` | 5 | decides | the study median has the lineage median's sign | not met |
| `xgboost-gbt` | 6 | decides | the maintainer rules the change | the maintainer's |
| `lightgbm-gbt` | 1 | decides | the lineage median difference's interval lies wholly above 0 | met |
| `lightgbm-gbt` | 2 | decides | better beats worse by the one-sided sign test at p < 0.05 after Holm over the 5 lineage pairs | met |
| `lightgbm-gbt` | 3 | decides | the twin meets 1 in neither of its pairs (gbt, lasso) and 2 in neither, each read apart | met |
| `lightgbm-gbt` | 3 | beside, decides nothing | beside, deciding nothing: the twin meets 1 and 2 in neither on gbt, this arm's learner | met |
| `lightgbm-gbt` | 4 | decides | the recall interval does not lie wholly below -0.02 | met |
| `lightgbm-gbt` | 5 | decides | the study median has the lineage median's sign | met |
| `lightgbm-gbt` | 6 | decides | the maintainer rules the change | the maintainer's |
| `elasticnet-lasso` | 1 | decides | the lineage median difference's interval lies wholly above 0 | met |
| `elasticnet-lasso` | 2 | decides | better beats worse by the one-sided sign test at p < 0.05 after Holm over the 5 lineage pairs | met |
| `elasticnet-lasso` | 3 | decides | the twin meets 1 in neither of its pairs (gbt, lasso) and 2 in neither, each read apart | met |
| `elasticnet-lasso` | 3 | beside, decides nothing | beside, deciding nothing: the twin meets 1 and 2 in neither on lasso, this arm's learner | met |
| `elasticnet-lasso` | 4 | decides | the recall interval does not lie wholly below -0.02 | met |
| `elasticnet-lasso` | 5 | decides | the study median has the lineage median's sign | met |
| `elasticnet-lasso` | 6 | decides | the maintainer rules the change | the maintainer's |

**Arms.** All three arms ran: `xgboost`, `lightgbm` and `elasticnet`. No arm was withdrawn.

**Pairs.** Each arm touches one learner. `xgboost` and `lightgbm` are read against `gbt`, and `elasticnet` against `lasso`, each against the default's own learner at the arm's seeds, 42, 43 and 44 (`k31-at-arm-seeds`). The twin, the default at seeds 45, 46 and 47, is read against the default on each learner, `twin-gbt` and `twin-lasso`.

**Conditions.** An arm carries the change when its pair meets condition 1 (the lineage median difference's interval lies wholly above 0) and condition 2 (better beats worse by the one-sided sign test at p < 0.05 after Holm over the 5 lineage pairs). Condition 3 reads both twin pairs for the whole block: the twin meets 1 in neither pair and 2 in neither, each read apart. Condition 4: the arm's recall interval does not lie wholly below -0.02. Condition 5: its study median has its lineage median's sign. Condition 6 is the maintainer's.

**Settings.** Each booster runs GRM's settings where GRM has them, on the rows GRM's screen kept, so the pair reads the library and leaves the screen as it is. The elastic net stops at tol 1e-4 or 2,000 epochs, and every elastic net result carries that limit as `docs/outside-learners.md` words it: saga's tol of 1e-4 stops on the change in the weights and bounds no distance to the optimum, and a difference between the elastic net and `lasso` reads saga's stop as well as the learner.

The release keeps its five learners and their settings whatever the count says, until the maintainer rules. Nothing changes in the release or on the site. An arm that meets conditions 1 to 5 would also add a library to the release, which waits as well on the maintainer's ruling on adding a dependency (plan block C7, row Release).

## Data

**Datasets.** The first run covers the 233 datasets of L less M in 23 species, each species and antibiotic named in full in the section Datasets. The 35 datasets of M wait for the label rule. A dataset's genomes are its genome count in `scripts/c7-outside/lists.tsv`, which counts each genome once.

**Draws.** The default and the twin come from plan block C3's run, the release's rule refitted with the same binary. Each arm runs three lineage draws a dataset, at fold seeds 42, 43 and 44, on the folds `grm export --seed` writes for that seed, checked against the default's folds by class and by genome before any fit. A dataset's value is the mean of its draws, and a dataset with a draw not scored is not scored.

**Study split.** The 55 first-run datasets in S get one study run per arm at seed 42 on the reviewed `study-whole` fold file, with no twin.

**Passes.** A species whose run started with the elastic net left out restarted with it on the same run folder, keeping its finished steps and booster fits. Its first copy, moved aside when it restarted, is named and set aside by its time and never read. The elastic net was not withdrawn, so every species is read from its current copy, whose `DONE` says the elastic net ran, in 23 of 23 species. All 23 species ran a second pass. Their 23 first copies were set aside by name: none was opened, read or sent to the reading, so the run's checks found 0 first copies beside the current ones. `page-withdrawn.tsv` holds its header alone.

**Checks.** Before any score was read, the run's checks passed 6 of 6, as the section Method gives them: 23 copies for the 23 species of the first run, 2,446 of 2,446 `SHA256SUMS` lines held, 2,262 of 2,262 library arm results found and none more, and 0 genome ids in 9,462 files against the 23,757 genome ids of the release's 23 species lists. Each check `plans.py` and `score.py` made in the run passed, as `page-checks.tsv` gives them:

| Script | Check | Passed |
|---|---|---|
| `plans.py` | the export at seed 42 is read | 233 of 233 |
| `plans.py` | the export at seed 43 is read | 233 of 233 |
| `plans.py` | the export at seed 44 is read | 233 of 233 |
| `plans.py` | the input fingerprint is the release report's | 233 of 233 |
| `plans.py` | the plan's folds by class are the default's | 699 of 699 |
| `plans.py` | the plan's folds by genome are the default's | 699 of 699 |
| `score.py` | a lineage result names the default's pick | 1,396 of 1,396 |
| `score.py` | a study result names no default of its own | 55 of 55 |
| `score.py` | a twin seed is the twin's and no draw's | 699 of 699 |
| `score.py` | an unavailable split has no fold | 2 of 2 |
| `score.py` | arm and library | 2,262 of 2,262 |
| `score.py` | calls are the test genomes' | 2,262 of 2,262 |
| `score.py` | cart: the outer calls give the split's balanced accuracy | 1,451 of 1,451 |
| `score.py` | dataset, scheme, seed and inner folds | 1,453 of 1,453 |
| `score.py` | dataset, split and seed | 2,262 of 2,262 |
| `score.py` | each fold's outer calls are its test genomes' | 1,451 of 1,451 |
| `score.py` | exists | 3,715 of 3,715 |
| `score.py` | fold 0: no inner fit | 1,451 of 1,451 |
| `score.py` | fold 1: no inner fit | 1,451 of 1,451 |
| `score.py` | fold 2: no inner fit | 1,451 of 1,451 |
| `score.py` | fold 3: no inner fit | 1,451 of 1,451 |
| `score.py` | fold 4: no inner fit | 1,451 of 1,451 |
| `score.py` | fold tests sum to the tested genomes | 1,451 of 1,451 |
| `score.py` | folds are the default's | 2,262 of 2,262 |
| `score.py` | folds in order | 1,451 of 1,451 |
| `score.py` | format and revision | 3,715 of 3,715 |
| `score.py` | gbt: the outer calls give the split's balanced accuracy | 1,451 of 1,451 |
| `score.py` | lasso: the outer calls give the split's balanced accuracy | 1,451 of 1,451 |
| `score.py` | listed keys only | 2,262 of 2,262 |
| `score.py` | pairs with its GRM learner | 2,262 of 2,262 |
| `score.py` | release report | 3,715 of 3,715 |
| `score.py` | rf: the outer calls give the split's balanced accuracy | 1,451 of 1,451 |
| `score.py` | scm: the outer calls give the split's balanced accuracy | 1,451 of 1,451 |
| `score.py` | seed 42: each learner's lineage figure is the release report's | 233 of 233 |
| `score.py` | seed 42: the default is held to the release report | 233 of 233 |
| `score.py` | seed 42: the default's learner is the one the release serves | 233 of 233 |
| `score.py` | seed 42: the default's random-split figure is the release report's | 233 of 233 |
| `score.py` | the default and the twin write no inputs file | 1,453 of 1,453 |
| `score.py` | the default's learner is the later of its tied set, as the release breaks a tie | 1,396 of 1,396 |
| `score.py` | the five learners in order | 1,453 of 1,453 |
| `score.py` | tie rule | 1,451 of 1,451 |

**Not scored.** One lineage dataset of the 233 is not scored: *Klebsiella pneumoniae* with polymyxin b (100 genomes). Of its 3 lineage draws, 2 are scored for each arm, for the twin and for the default each pair reads against, so every pair enters 232 of 233 lineage datasets and no missing draw is read as zero. `score.py`'s check that an unavailable split has no fold passed 2 of 2. Every pair enters 55 of 55 study datasets.

## Count

The primary measure is each pair's lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. Every difference is the arm's value less the default's; for the very major error rate, the major error rate and log loss, lower is better, so better and worse turn round. The interval is a 95% interval of the median over 4,000 resamples of species. The sign test is one-sided, and Holm runs over the 5 lineage pairs on balanced accuracy, withdrawn or not.

The five pairs on the lineage split's balanced accuracy, as `count.tsv` gives them:

| Pair | Against | Better | Worse | Same | Median | Interval | Sign p | Holm p | Condition 1 | Condition 2 |
|---|---|---|---|---|---|---|---|---|---|---|
| `xgboost-gbt` | `k31-at-arm-seeds` | 58 of 232 | 17 of 232 | 157 of 232 | 0.003 | 0.001 to 0.005 | < 0.001 | < 0.001 | met | met |
| `lightgbm-gbt` | `k31-at-arm-seeds` | 62 of 232 | 14 of 232 | 156 of 232 | 0.003 | 0.001 to 0.004 | < 0.001 | < 0.001 | met | met |
| `elasticnet-lasso` | `k31-at-arm-seeds` | 71 of 232 | 19 of 232 | 142 of 232 | 0.006 | 0.003 to 0.009 | < 0.001 | < 0.001 | met | met |
| `twin-gbt` | `k31` | 40 of 232 | 33 of 232 | 159 of 232 | 0.000 | 0.000 to 0.001 | 0.241 | 0.243 | not met | not met |
| `twin-lasso` | `k31` | 28 of 232 | 19 of 232 | 185 of 232 | 0.000 | -0.001 to 0.000 | 0.121 | 0.243 | not met | not met |

`score.py` names the pairs `xgboost-k31`, `lightgbm-k31`, `elasticnet-k31`, `twin-k31` on `gbt` and `twin-k31` on `lasso`; `page-count.tsv` keeps its name in the column `score_pair` beside the name this page gives.

## Per arm

**Measures.** Each pair on each measure, over the datasets it enters: 232 of 233 on the lineage split and 55 of 55 on the study split. Recall stands beside log loss in each pair's rows.

| Split | Pair | Measure | Better is | Entered | Better | Worse | Same | Median | Interval |
|---|---|---|---|---|---|---|---|---|---|
| lineage | `xgboost-gbt` | recall | higher | 232 of 233 | 63 | 32 | 137 | 0.002 | 0.000 to 0.005 |
| lineage | `xgboost-gbt` | specificity | higher | 232 of 233 | 47 | 37 | 148 | 0.000 | 0.000 to 0.002 |
| lineage | `xgboost-gbt` | very major error rate | lower | 232 of 233 | 63 | 32 | 137 | -0.002 | -0.005 to 0.000 |
| lineage | `xgboost-gbt` | major error rate | lower | 232 of 233 | 47 | 37 | 148 | 0.000 | -0.002 to 0.000 |
| lineage | `xgboost-gbt` | log loss | lower | 232 of 233 | 53 | 32 | 147 | -0.003 | -0.005 to -0.001 |
| lineage | `lightgbm-gbt` | recall | higher | 232 of 233 | 66 | 31 | 135 | 0.001 | 0.000 to 0.005 |
| lineage | `lightgbm-gbt` | specificity | higher | 232 of 233 | 46 | 37 | 149 | 0.001 | 0.000 to 0.003 |
| lineage | `lightgbm-gbt` | very major error rate | lower | 232 of 233 | 66 | 31 | 135 | -0.001 | -0.005 to 0.000 |
| lineage | `lightgbm-gbt` | major error rate | lower | 232 of 233 | 46 | 37 | 149 | -0.001 | -0.003 to 0.000 |
| lineage | `lightgbm-gbt` | log loss | lower | 232 of 233 | 52 | 39 | 141 | -0.002 | -0.004 to 0.001 |
| lineage | `elasticnet-lasso` | recall | higher | 232 of 233 | 68 | 48 | 116 | 0.002 | 0.000 to 0.007 |
| lineage | `elasticnet-lasso` | specificity | higher | 232 of 233 | 71 | 27 | 134 | 0.003 | 0.000 to 0.010 |
| lineage | `elasticnet-lasso` | very major error rate | lower | 232 of 233 | 68 | 48 | 116 | -0.002 | -0.007 to 0.000 |
| lineage | `elasticnet-lasso` | major error rate | lower | 232 of 233 | 71 | 27 | 134 | -0.003 | -0.010 to 0.000 |
| lineage | `elasticnet-lasso` | log loss | lower | 232 of 233 | 121 | 21 | 90 | -0.023 | -0.037 to -0.013 |
| lineage | `twin-gbt` | recall | higher | 232 of 233 | 53 | 47 | 132 | 0.000 | 0.000 to 0.000 |
| lineage | `twin-gbt` | specificity | higher | 232 of 233 | 37 | 38 | 157 | 0.000 | 0.000 to 0.001 |
| lineage | `twin-gbt` | very major error rate | lower | 232 of 233 | 53 | 47 | 132 | 0.000 | 0.000 to 0.000 |
| lineage | `twin-gbt` | major error rate | lower | 232 of 233 | 37 | 38 | 157 | 0.000 | -0.001 to 0.000 |
| lineage | `twin-gbt` | log loss | lower | 232 of 233 | 36 | 48 | 148 | 0.001 | -0.001 to 0.003 |
| lineage | `twin-lasso` | recall | higher | 232 of 233 | 41 | 32 | 159 | 0.000 | 0.000 to 0.001 |
| lineage | `twin-lasso` | specificity | higher | 232 of 233 | 33 | 31 | 168 | 0.000 | -0.003 to 0.000 |
| lineage | `twin-lasso` | very major error rate | lower | 232 of 233 | 41 | 32 | 159 | 0.000 | -0.001 to 0.000 |
| lineage | `twin-lasso` | major error rate | lower | 232 of 233 | 33 | 31 | 168 | 0.000 | 0.000 to 0.003 |
| lineage | `twin-lasso` | log loss | lower | 232 of 233 | 52 | 54 | 126 | 0.000 | -0.002 to 0.006 |
| study | `xgboost-gbt` | balanced accuracy | higher | 55 of 55 | 11 | 15 | 29 | -0.002 | -0.020 to 0.002 |
| study | `xgboost-gbt` | recall | higher | 55 of 55 | 15 | 16 | 24 | 0.000 | -0.023 to 0.005 |
| study | `xgboost-gbt` | specificity | higher | 55 of 55 | 12 | 17 | 26 | -0.004 | -0.019 to 0.000 |
| study | `xgboost-gbt` | very major error rate | lower | 55 of 55 | 15 | 16 | 24 | 0.000 | -0.005 to 0.023 |
| study | `xgboost-gbt` | major error rate | lower | 55 of 55 | 12 | 17 | 26 | 0.004 | 0.000 to 0.019 |
| study | `xgboost-gbt` | log loss | lower | 55 of 55 | 14 | 23 | 18 | 0.010 | 0.003 to 0.034 |
| study | `lightgbm-gbt` | balanced accuracy | higher | 55 of 55 | 7 | 13 | 35 | 0.001 | -0.012 to 0.007 |
| study | `lightgbm-gbt` | recall | higher | 55 of 55 | 14 | 16 | 25 | 0.000 | -0.004 to 0.009 |
| study | `lightgbm-gbt` | specificity | higher | 55 of 55 | 12 | 15 | 28 | 0.000 | -0.017 to 0.000 |
| study | `lightgbm-gbt` | very major error rate | lower | 55 of 55 | 14 | 16 | 25 | 0.000 | -0.009 to 0.002 |
| study | `lightgbm-gbt` | major error rate | lower | 55 of 55 | 12 | 15 | 28 | 0.000 | 0.000 to 0.017 |
| study | `lightgbm-gbt` | log loss | lower | 55 of 55 | 16 | 19 | 20 | 0.010 | 0.003 to 0.029 |
| study | `elasticnet-lasso` | balanced accuracy | higher | 55 of 55 | 14 | 7 | 34 | 0.004 | 0.000 to 0.007 |
| study | `elasticnet-lasso` | recall | higher | 55 of 55 | 10 | 14 | 31 | 0.000 | -0.007 to 0.002 |
| study | `elasticnet-lasso` | specificity | higher | 55 of 55 | 17 | 8 | 30 | 0.007 | 0.004 to 0.011 |
| study | `elasticnet-lasso` | very major error rate | lower | 55 of 55 | 10 | 14 | 31 | 0.000 | 0.000 to 0.007 |
| study | `elasticnet-lasso` | major error rate | lower | 55 of 55 | 17 | 8 | 30 | -0.007 | -0.011 to -0.004 |
| study | `elasticnet-lasso` | log loss | lower | 55 of 55 | 19 | 13 | 23 | -0.009 | -0.013 to 0.004 |

**Folds.** Per arm and split, counts of folds over the folds of its reports: for each booster the folds that kept one tree; for the elastic net the refits that reached 2,000 epochs, the folds whose fit kept no weight, and the folds scored at their training share.

| Arm | Split | Folds | Of the folds |
|---|---|---|---|
| XGBoost 3.0.5 | lineage | 3,490 | one tree folds 99 of 3,490; constant folds 0 of 3,490 |
| XGBoost 3.0.5 | study | 275 | one tree folds 5 of 275; constant folds 0 of 275 |
| LightGBM 4.6.0 | lineage | 3,490 | one tree folds 96 of 3,490; constant folds 0 of 3,490 |
| LightGBM 4.6.0 | study | 275 | one tree folds 2 of 275; constant folds 0 of 275 |
| scikit-learn 1.7.2 elastic net | lineage | 3,490 | unconverged refits 59 of 3,490; null fit folds 221 of 3,490; constant folds 221 of 3,490 |
| scikit-learn 1.7.2 elastic net | study | 275 | unconverged refits 4 of 275; null fit folds 1 of 275; constant folds 1 of 275 |

**Libraries.** The versions each species' run read back, as one row when every species agrees, or a row per species.

| Species | python | numpy | scipy | scikit-learn | XGBoost | LightGBM |
|---|---|---|---|---|---|---|
| every species | 3.12.3 | 1.26.4 | 1.15.3 | 1.7.2 | 3.0.5 | 4.6.0 |

One code sha256 in 23 of 23 species: `d3de72fcc10f`.

The grm sha256 and the code sha256 that each species' C7 `DONE` names stand in that species' C3 `binary.txt`, in 23 of 23 species, compared at 14:13Z on 25 Nov 2026, before the reading.

**Time.** Per arm, the fits timed, with the sums and medians of their wall and CPU seconds.

| Arm | Fits | Wall s, sum | Wall s, median | CPU s, sum | CPU s, median |
|---|---|---|---|---|---|
| XGBoost 3.0.5 | 754 | 8799.714 | 8.009 | 8715.400 | 7.990 |
| LightGBM 4.6.0 | 754 | 5047.740 | 4.878 | 5005.770 | 4.815 |
| scikit-learn 1.7.2 elastic net | 754 | 5419925.100 | 5210.252 | 5416010.210 | 5208.695 |

## Beside it

These counts stand beside the count. They are in no Holm family and decide nothing.

**Own twin pair.** Each arm against its own learner's twin pair alone, condition 3 read on that learner.

**Without constant folds.** The count with every dataset that holds a draw scored at its training share left out (`count-no-constant.tsv`).

**Without unconverged refits.** The count with every dataset whose elastic net draws hold a refit that reached 2,000 epochs left out (`count-no-unconverged.tsv`).

**Own twin pair.** The rows of the conditions table in the section Reading marked "beside, decides nothing" read it: met for `xgboost` and `lightgbm` on `gbt` and for `elasticnet` on `lasso`.

**Without constant folds.** Of the 232 lineage datasets `elasticnet-lasso` enters, 41 hold a draw scored at its training share. Without them it enters 191 of 233: better on 55, worse on 8 and the same on 128, with a median of 0.006 and an interval of 0.003 to 0.009. The boosters and the twin hold no constant fold, so their rows are the count's.

| Pair | Against | Better | Worse | Same | Median | Interval | Sign p | Holm p |
|---|---|---|---|---|---|---|---|---|
| `xgboost-gbt` | `k31-at-arm-seeds` | 58 of 232 | 17 of 232 | 157 of 232 | 0.003 | 0.001 to 0.005 | < 0.001 | < 0.001 |
| `lightgbm-gbt` | `k31-at-arm-seeds` | 62 of 232 | 14 of 232 | 156 of 232 | 0.003 | 0.001 to 0.004 | < 0.001 | < 0.001 |
| `elasticnet-lasso` | `k31-at-arm-seeds` | 55 of 191 | 8 of 191 | 128 of 191 | 0.006 | 0.003 to 0.009 | < 0.001 | < 0.001 |
| `twin-gbt` | `k31` | 40 of 232 | 33 of 232 | 159 of 232 | 0.000 | 0.000 to 0.001 | 0.241 | 0.243 |
| `twin-lasso` | `k31` | 28 of 232 | 19 of 232 | 185 of 232 | 0.000 | -0.001 to 0.000 | 0.121 | 0.243 |

**Without unconverged refits.** Of the 232, 46 hold an elastic net refit that reached 2,000 epochs. Without them `elasticnet-lasso` enters 186 of 233: better on 54, worse on 14 and the same on 118, with a median of 0.005 and an interval of 0.003 to 0.008.

| Pair | Against | Better | Worse | Same | Median | Interval | Sign p | Holm p |
|---|---|---|---|---|---|---|---|---|
| `xgboost-gbt` | `k31-at-arm-seeds` | 58 of 232 | 17 of 232 | 157 of 232 | 0.003 | 0.001 to 0.005 | < 0.001 | < 0.001 |
| `lightgbm-gbt` | `k31-at-arm-seeds` | 62 of 232 | 14 of 232 | 156 of 232 | 0.003 | 0.001 to 0.004 | < 0.001 | < 0.001 |
| `elasticnet-lasso` | `k31-at-arm-seeds` | 54 of 186 | 14 of 186 | 118 of 186 | 0.005 | 0.003 to 0.008 | < 0.001 | < 0.001 |
| `twin-gbt` | `k31` | 40 of 232 | 33 of 232 | 159 of 232 | 0.000 | 0.000 to 0.001 | 0.241 | 0.243 |
| `twin-lasso` | `k31` | 28 of 232 | 19 of 232 | 185 of 232 | 0.000 | -0.001 to 0.000 | 0.121 | 0.243 |

## Seen after

These points are seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**Species.** Per species, named in full, the rows are dataset and pair counts, with their denominator. They decide nothing.

| Species | `xgboost-gbt`: better / worse / same | `lightgbm-gbt`: better / worse / same | `elasticnet-lasso`: better / worse / same | `twin-gbt`: better / worse / same | `twin-lasso`: better / worse / same |
|---|---|---|---|---|---|
| *Acinetobacter baumannii* | 5 / 0 / 2 of 7 | 5 / 1 / 1 of 7 | 3 / 1 / 3 of 7 | 3 / 2 / 2 of 7 | 1 / 4 / 2 of 7 |
| *Campylobacter coli* | 1 / 0 / 7 of 8 | 1 / 0 / 7 of 8 | 1 / 0 / 7 of 8 | 1 / 1 / 6 of 8 | 0 / 0 / 8 of 8 |
| *Campylobacter jejuni* | 0 / 1 / 3 of 4 | 1 / 0 / 3 of 4 | 0 / 1 / 3 of 4 | 1 / 0 / 3 of 4 | 1 / 0 / 3 of 4 |
| *Clostridioides difficile* | 0 / 0 / 4 of 4 | 0 / 0 / 4 of 4 | 0 / 0 / 4 of 4 | 0 / 1 / 3 of 4 | 0 / 0 / 4 of 4 |
| *Corynebacterium diphtheriae* | 2 / 0 / 1 of 3 | 1 / 0 / 2 of 3 | 0 / 2 / 1 of 3 | 2 / 0 / 1 of 3 | 1 / 1 / 1 of 3 |
| *Enterobacter cloacae* | 5 / 0 / 4 of 9 | 5 / 0 / 4 of 9 | 4 / 0 / 5 of 9 | 4 / 0 / 5 of 9 | 2 / 1 / 6 of 9 |
| *Enterococcus faecium* | 1 / 0 / 10 of 11 | 1 / 0 / 10 of 11 | 2 / 0 / 9 of 11 | 2 / 1 / 8 of 11 | 1 / 1 / 9 of 11 |
| *Escherichia coli* | 7 / 2 / 20 of 29 | 9 / 2 / 18 of 29 | 10 / 1 / 18 of 29 | 1 / 3 / 25 of 29 | 3 / 0 / 26 of 29 |
| *Klebsiella aerogenes* | 1 / 0 / 2 of 3 | 1 / 0 / 2 of 3 | 1 / 0 / 2 of 3 | 0 / 0 / 3 of 3 | 2 / 1 / 0 of 3 |
| *Klebsiella oxytoca* | 2 / 0 / 2 of 4 | 2 / 0 / 2 of 4 | 4 / 0 / 0 of 4 | 2 / 1 / 1 of 4 | 1 / 1 / 2 of 4 |
| *Klebsiella pneumoniae* | 10 / 3 / 18 of 31 | 11 / 1 / 19 of 31 | 16 / 2 / 13 of 31 | 4 / 6 / 21 of 31 | 6 / 4 / 21 of 31 |
| *Klebsiella quasipneumoniae* | 0 / 0 / 3 of 3 | 0 / 0 / 3 of 3 | 0 / 1 / 2 of 3 | 1 / 0 / 2 of 3 | 0 / 0 / 3 of 3 |
| *Mycobacterium tuberculosis* | 5 / 0 / 15 of 20 | 4 / 0 / 16 of 20 | 3 / 2 / 15 of 20 | 5 / 2 / 13 of 20 | 2 / 1 / 17 of 20 |
| *Neisseria gonorrhoeae* | 2 / 0 / 4 of 6 | 2 / 0 / 4 of 6 | 2 / 0 / 4 of 6 | 2 / 0 / 4 of 6 | 0 / 1 / 5 of 6 |
| *Pseudomonas aeruginosa* | 3 / 3 / 6 of 12 | 3 / 2 / 7 of 12 | 5 / 1 / 6 of 12 | 1 / 4 / 7 of 12 | 0 / 0 / 12 of 12 |
| *Salmonella enterica* | 7 / 6 / 11 of 24 | 7 / 5 / 12 of 24 | 9 / 5 / 10 of 24 | 4 / 4 / 16 of 24 | 3 / 2 / 19 of 24 |
| *Shigella flexneri* | 0 / 0 / 4 of 4 | 0 / 0 / 4 of 4 | 1 / 0 / 3 of 4 | 1 / 0 / 3 of 4 | 1 / 0 / 3 of 4 |
| *Shigella sonnei* | 2 / 1 / 17 of 20 | 2 / 1 / 17 of 20 | 4 / 0 / 16 of 20 | 0 / 2 / 18 of 20 | 1 / 0 / 19 of 20 |
| *Staphylococcus aureus* | 2 / 0 / 11 of 13 | 3 / 2 / 8 of 13 | 1 / 3 / 9 of 13 | 1 / 4 / 8 of 13 | 0 / 0 / 13 of 13 |
| *Staphylococcus haemolyticus* | 1 / 0 / 2 of 3 | 1 / 0 / 2 of 3 | 2 / 0 / 1 of 3 | 2 / 1 / 0 of 3 | 1 / 1 / 1 of 3 |
| *Streptococcus agalactiae* | 0 / 0 / 3 of 3 | 0 / 0 / 3 of 3 | 1 / 0 / 2 of 3 | 0 / 0 / 3 of 3 | 0 / 0 / 3 of 3 |
| *Streptococcus pneumoniae* | 2 / 1 / 5 of 8 | 2 / 0 / 6 of 8 | 2 / 0 / 6 of 8 | 2 / 1 / 5 of 8 | 1 / 0 / 7 of 8 |
| *Vibrio cholerae* | 0 / 0 / 3 of 3 | 1 / 0 / 2 of 3 | 0 / 0 / 3 of 3 | 1 / 0 / 2 of 3 | 1 / 1 / 1 of 3 |

**Time.** Over 754 fits an arm, the elastic net's median fit took 5210.252 wall seconds, against 8.009 for XGBoost and 4.878 for LightGBM.

**Log loss.** On the lineage split the elastic net's log loss is better on 121 of 232 datasets and worse on 21, with a median of -0.023 and an interval of -0.037 to -0.013, while its recall is better on 68 and worse on 48, with a median of 0.002. On the study split both boosters' log loss is worse: XGBoost's median is 0.010 (0.003 to 0.034) and LightGBM's 0.010 (0.003 to 0.029), on 55 of 55, with a recall median of 0.000 beside each.

**Twin movement.** By its seeds alone, the twin moves 73 of 232 lineage datasets by more than 0.02 on `gbt` and 47 on `lasso`. The arms move 75 (`xgboost-gbt`), 76 (`lightgbm-gbt`) and 90 (`elasticnet-lasso`) of 232.

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Elastic net stop.** The elastic net stops at tol 1e-4 or 2,000 epochs. saga's tol of 1e-4 stops on the change in the weights and bounds no distance to the optimum, so a difference between the elastic net and `lasso` reads saga's stop as well as the learner.

**Screen.** The boosters fit the rows GRM's screen kept, 4,000 per outer fold, so a pair reads the library on GRM's screen and says nothing of another screen.

**Study split.** One study run per dataset at seed 42, and no study twin.

**Labels.** The 35 datasets of M wait for the label rule.

**Markers.** No feature, unitig or k-mer is named or ranked. No elastic net weight and no booster's feature importance leaves the box. The page holds aggregates over datasets, species, splits, arms and folds only.

## Method

**Checks of the run.** Before any score was read, `scripts/c7-outside/run_checks.py check` read each species' copy and printed counts only: the species, each `SHA256SUMS`, each `DONE` with its elastic net pass, `expect.py` against each species' results, the genome-id count and the count of unitigs and matrix column names over every file that is read or published. `score/score.log` stays out of the record.

| Check | Printed |
|---|---|
| species | 23 copies, 23 species in lists.tsv's first run; 0 first copies set aside by time, 0 beside no newer DONE |
| SHA256SUMS | 23 of 23 species; 2446 of 2446 lines hold; results/ and score.tar's tables are the listed files in 23 of 23; stage 1's c3/score/ in the score.tar of 0 species, skipped unread |
| DONE | 23 of 23 DONE lines name their species, export binary, result count and SHA256SUMS; FAILED not set aside in 0 species; the elastic net run in 23 of 23 |
| expect.py | 23 of 23 species: lists.txt names the first run, the seeds, twin seeds, arm seeds and the arms of its pass, expect.tsv is expect.py's, and results/outside/ is exactly its library arm rows; 2262 of 2262 results found, 0 more |
| genome ids | 0 genome ids in 9462 files, against 23757 release ids of 23 species lists |
| features | 0 runs of 31 or more nucleotides, 0 LightGBM and 0 XGBoost column names in 9462 files |

`CHECKS 6 of 6`: all six passed. The `DONE` lines name export commit `da70f4b`, grm sha256 `d3b53997a02c`, code sha256 `d3de72fcc10f` and export grm sha256 `637b2da2b59b` in 23 of 23 species each, and the elastic net run in 23 of 23. `expect.py` found 233 datasets, 55 study datasets and 23 species.

**Reading.** `scripts/c7-outside/run_checks.py gather` copied each species' `score/` tables, as its run wrote them, and its run files into the reading, and `scripts/c7-outside/score.py --count` counted them once.

Each ran once, in this order, under `env -i` at hash seed 0, and each ended with 0. `<back dir>` holds the 23 current copies and `<release dir>` the release's 23 species lists.

| Step | Command | Started | Ended | Last line |
|---|---|---|---|---|
| check | `python3 -I scripts/c7-outside/run_checks.py check <back dir> scripts/c7-outside/lists.tsv <release dir>` | 14:24:49Z | 14:24:50Z | `CHECKS 6 of 6` |
| gather | `python3 -I scripts/c7-outside/run_checks.py gather <back dir> scripts/c7-outside/lists.tsv <release dir> reading` | 14:24:50Z | 14:24:52Z | `gather: wrote reading; next: python3 -I scripts/c7-outside/score.py --count reading/count reading/score/*/datasets.tsv` |
| count | `python3 -I scripts/c7-outside/score.py --count reading/count reading/score/*/datasets.tsv` | 14:24:52Z | 14:25:05Z | `count: 2482 dataset rows from 23 files, 233 of 233 first-run datasets` |
| write | `python3 -I scripts/c7-outside/page_tables.py write reading scripts/c7-outside/lists.tsv reference/outside-learners` | 14:25:05Z | 14:25:18Z | `gates: 17 of 17 pass; wrote 19 files` |
| print | `python3 -I scripts/c7-outside/page_tables.py print reference/outside-learners` | 14:25:18Z | 14:25:18Z | the tables of this page |

`reading.txt` gives the sha256 of the 2,425 files `gather` copied and of the four files the count wrote:

| File | sha256 |
|---|---|
| `count/conditions.tsv` | `4ccaaba9d9afb960228bd1675ecfd3a2f5db2eca57240bc853c7cbfe1fbd2b27` |
| `count/count-no-constant.tsv` | `2dfed4f1a62a527bc8b2ff2310a727edb323cf132a2d57dd7feb737746a3d265` |
| `count/count-no-unconverged.tsv` | `367bd927bc5c05ed4b92fb18aea40d6906ae6cd050be2c09b83dc9645a1cffde` |
| `count/count.tsv` | `0286758d3d7460964fa872a6a4bcebc546dd281159c8a9a0d289a018b8a5234b` |

**Tables.** `scripts/c7-outside/page_tables.py write` checked the reading against 17 gates and wrote `reference/outside-learners/`. The four count files are `score.py`'s own tables: the script counted them again, copied them and summed them.

| Gate | What it holds |
|---|---|
| 1. reading | the reading holds `reading.txt`, `checks.txt` with every check passed, `withdrawn.tsv`, `count/` with its four files, `score/` and `run/`, and nothing else; each withdrawn arm is named with a reason and a time |
| 2. species | `score/` and `run/` each hold every species of the first run once; a folder named `<species slug>.<UTC time>` is a first copy set aside and never read as a species |
| 3. M | no dataset of M, or outside the first run, stands in any table |
| 4. pairs | `xgboost` and `lightgbm` stand on `gbt` only and `elasticnet` on `lasso` only, each read against `k31-at-arm-seeds` |
| 5. twin | the twin's rows are lineage draws at seeds 45, 46 and 47 only, and the conditions never read the twin as an arm |
| 6. study | a study row is a study run at seed 42 on the study split, and is never counted on the lineage split |
| 7. family | each count holds the five lineage pairs, withdrawn or not, each with its Holm p |
| 8. withdrawn | a withdrawn arm holds no row, no report and no fold, and its pair stays at p = 1 |
| 9. passes | a species whose `DONE` gives its first pass alone holds no elastic net result |
| 10. rows | each species holds exactly its first-run datasets' rows |
| 11. missing | a draw not scored holds no figure and is never read as zero |
| 12. draws | `datasets.tsv` is the mean of `draws.tsv`'s draws |
| 13. reports | one report per library arm draw and no other, each of its listed keys only, its fold counts the sums over its folds |
| 14. checks | every `score.py` check and every `plans.py` check passed |
| 15. binary | one code sha256 and one version of each library in every species, or each named with its species |
| 16. beside | `count-no-constant.tsv`, `count-no-unconverged.tsv` and each arm's own twin pair are written as deciding nothing |
| 17. count | the four count files are what `score.py --count` writes from the `datasets.tsv` files |

The fold counts that read 0 by construction, a booster's unconverged refits and null fits and the elastic net's one-tree folds, stay in `page-folds.tsv` and off this page.

**Tests.** `tests/c7_outside/test_page_tables.py` and `tests/c7_outside/test_run_checks.py` hold SYNTHETIC tests of the table script and the run's checks, on three made-up species. `scripts/c7-outside/page_faults.py` plants each fault of its table in a copy of the tree and runs the test meant to catch it. `tests/c7_outside/test_result_record.py` holds 15 tests of this page against `reference/outside-learners/`: each table stands as `page_tables.py print` prints it from those files, and each figure in this page's sentences is rebuilt from them.

| Fault | Gate | Plants | Test | Read |
|---|---|---|---|---|
| 1 | reading | a folder outside the reading passed | `test_page_tables.Gates.test_a_folder_outside_the_reading_is_refused` | caught |
| 2a | species | a first copy set aside read as its species | `test_page_tables.Gates.test_a_first_copy_set_aside_is_never_read_as_a_species` | caught |
| 2b | species | a species outside the first run passed | `test_page_tables.Gates.test_a_species_outside_the_first_run_is_refused` | caught |
| 3 | M | a dataset of M counted | `test_page_tables.Gates.test_a_dataset_of_m_is_never_counted` | caught |
| 4a | pairs | the elastic net read on gbt | `test_page_tables.Gates.test_an_elastic_net_read_against_gbt_is_refused` | caught |
| 4b | pairs | an arm read against the default's every draw, not at its own seeds | `test_page_tables.Tables.test_an_arm_is_read_against_the_default_at_its_own_seeds` | caught |
| 5 | twin | the twin read as an arm in conditions.tsv | `test_page_tables.Gates.test_conditions_that_read_the_twin_as_an_arm_are_refused` | caught |
| 6a | study | a study row counted on the lineage split by the gate | `test_page_tables.Gates.test_a_study_run_moved_onto_the_lineage_split_is_refused` | caught |
| 6b | study | a study row counted on the lineage split in the species table | `test_page_tables.Tables.test_a_study_row_stays_on_the_study_split` | caught |
| 7 | family | a family short of a pair passed | `test_page_tables.Gates.test_a_family_short_of_a_pair_is_refused` | caught |
| 8a | withdrawn | a withdrawn arm's reports let through | `test_page_tables.Gates.test_a_withdrawn_arms_report_is_refused` | caught |
| 8b | withdrawn | a withdrawn arm's pair read off p = 1 | `test_page_tables.Gates.test_a_withdrawn_pair_off_p_one_is_refused` | caught |
| 9 | passes | a first pass alone read as a whole run | `test_page_tables.Gates.test_a_first_pass_alone_is_refused_unless_the_elastic_net_is_withdrawn` | caught |
| 10 | rows | a first-run dataset left out of datasets.tsv | `test_page_tables.Gates.test_a_dataset_left_out_of_datasets_tsv_is_refused` | caught |
| 11a | missing | a missing draw's figures let through | `test_page_tables.Gates.test_a_missing_draw_with_figures_is_refused` | caught |
| 11b | missing | a dataset not scored read as zero in the dataset table | `test_page_tables.Tables.test_a_dataset_not_scored_holds_no_difference` | caught |
| 12 | draws | a dataset figure that is not its draws' mean | `test_page_tables.Gates.test_a_dataset_figure_that_is_not_its_draws_mean_is_refused` | caught |
| 13a | reports | a report's fold counts not held to its folds | `test_page_tables.Gates.test_a_report_whose_fold_count_is_not_its_folds_sum_is_refused` | caught |
| 13b | reports | a report with no draw let through | `test_page_tables.Gates.test_a_report_with_no_draw_is_refused` | caught |
| 14 | checks | a failed score check let through | `test_page_tables.Gates.test_a_failed_score_check_is_refused` | caught |
| 15a | binary | library versions that differ written as one row | `test_page_tables.Tables.test_versions_that_differ_are_named_with_each_species` | caught |
| 15b | binary | a DONE's result count not held to the reports | `test_page_tables.Gates.test_a_done_with_another_result_count_is_refused` | caught |
| 16a | beside | the no-constant and no-unconverged counts written as deciding | `test_page_tables.Tables.test_the_count_rows_decide_for_count_tsv_alone` | caught |
| 16b | beside | the arm against its own learner's twin pair alone written as deciding | `test_page_tables.Tables.test_the_beside_condition_decides_nothing` | caught |
| 17 | count | count files score.py did not write let through | `test_page_tables.Gates.test_a_count_figure_score_py_did_not_write_is_refused` | caught |
| R1 | run: species | a first copy newer than its DONE set aside | `test_run_checks.RunChecks.test_a_first_copy_newer_than_its_done_fails` | caught |
| R2 | run: SHA256SUMS | a file that differs from its SHA256SUMS line let through | `test_run_checks.RunChecks.test_a_changed_result_fails_its_sha256sums` | caught |
| R3 | run: DONE | a first pass alone beside second passes let through | `test_run_checks.RunChecks.test_a_first_pass_beside_second_passes_fails_done` | caught |
| R4 | run: DONE | a FAILED newer than its DONE set aside | `test_run_checks.RunChecks.test_a_failed_newer_than_its_done_fails` | caught |
| R5 | run: expect.py | an expect.tsv that is not expect.py's let through | `test_run_checks.RunChecks.test_an_expect_tsv_that_is_not_expect_pys_fails` | caught |
| R6 | run: genome ids | a genome id in a log let through | `test_run_checks.RunChecks.test_a_genome_id_in_a_log_fails` | caught |
| R7 | run: features | a matrix column name in a log let through | `test_run_checks.RunChecks.test_a_column_name_in_a_log_fails` | caught |
| R9 | run: species | a second-pass copy read when the elastic net is withdrawn after some second passes ran | `test_run_checks.RunChecks.test_late_withdrawal_reads_each_first_copy_and_never_opens_a_second_pass` | caught |
| R10 | reading | the count of species that ran a second pass left off the page | `test_page_tables.Gates.test_the_second_pass_count_goes_to_the_page` | caught |
| R8 | run: score.log | score/score.log read with score.py's tables | `test_run_checks.ScoreLog.test_score_log_is_neither_read_nor_copied` | caught |

All 35 were caught, before any score was opened.

**Files.** The reading's tables are in `reference/outside-learners/`. `scripts/c7-outside/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets

The first run's 233 datasets, per species: the datasets, those in S, and each antibiotic with the dataset's genomes. Per-dataset figures are in `reference/outside-learners/page-datasets.tsv`, with species and antibiotic in full.

| Species | Datasets | In S | Antibiotics, with genomes |
|---|---:|---:|---|
| *Acinetobacter baumannii* | 7 | 2 | cefepime (437), cefotaxime (808), doripenem (213), doxycycline (167), minocycline (139), moxifloxacin (91), sulfamethoxazole/trimethoprim (1,114) |
| *Campylobacter coli* | 8 | 0 | azithromycin (389), ciprofloxacin (391), clindamycin (389), erythromycin (391), gentamicin (391), nalidixic acid (375), telithromycin (374), tetracycline (391) |
| *Campylobacter jejuni* | 4 | 0 | ciprofloxacin (508), clindamycin (508), nalidixic acid (506), tetracycline (508) |
| *Clostridioides difficile* | 4 | 0 | azithromycin (203), ceftriaxone (136), clarithromycin (203), moxifloxacin (204) |
| *Corynebacterium diphtheriae* | 3 | 0 | sulfamethoxazole/trimethoprim (236), tetracycline (243), trimethoprim (229) |
| *Enterobacter cloacae* | 9 | 8 | cefepime (103), ciprofloxacin (110), gentamicin (104), imipenem (86), levofloxacin (104), meropenem (99), sulfamethoxazole/trimethoprim (80), tobramycin (99), trimethoprim (124) |
| *Enterococcus faecium* | 11 | 0 | ampicillin (1,680), clindamycin (901), dalfopristin/quinupristin (1,118), doxycycline (345), erythromycin (1,151), gentamicin (1,345), kanamycin (944), streptomycin (1,239), teicoplanin (1,192), tetracycline (958), vancomycin (1,741) |
| *Escherichia coli* | 29 | 15 | amikacin (1,125), ampicillin/sulbactam (283), aztreonam (699), cefazolin (414), cefepime (1,124), cefotaxime (1,149), cefoxitin (1,066), cefpodoxime (175), cefpodoxime/clavulanate (136), ceftazidime (1,466), ceftiofur (203), ceftriaxone (923), cefuroxime (548), cephalothin (472), chloramphenicol (709), colistin (239), doxycycline (227), ertapenem (820), erythromycin (205), imipenem (614), levofloxacin (403), nalidixic acid (406), nitrofurantoin (482), piperacillin/tazobactam (979), streptomycin (386), sulfamethoxazole (341), tetracycline (782), tobramycin (1,070), trimethoprim (861) |
| *Klebsiella aerogenes* | 3 | 0 | ertapenem (83), imipenem (68), meropenem (78) |
| *Klebsiella oxytoca* | 4 | 0 | cefotaxime (98), cefoxitin (99), cefpodoxime (74), ertapenem (112) |
| *Klebsiella pneumoniae* | 32 | 18 | ampicillin/sulbactam (744), aztreonam (1,050), cefazolin (1,141), cefepime (1,233), cefotetan (61), cefoxitin (1,181), cefpodoxime (221), cefpodoxime/ceftriaxone (56), cefpodoxime/clavulanate (207), ceftazidime (1,310), ceftazidime/avibactam (259), ceftriaxone (1,208), cefuroxime (565), chloramphenicol (409), ciprofloxacin (1,706), colistin (336), doripenem (243), doxycycline (141), ertapenem (952), fosfomycin (294), imipenem (1,158), levofloxacin (871), nitrofurantoin (562), ofloxacin (71), polymyxin b (100), streptomycin (80), sulfamethoxazole/trimethoprim (1,546), tetracycline (741), ticarcillin/clavulanate (250), tigecycline (924), tobramycin (1,230), trimethoprim (258) |
| *Klebsiella quasipneumoniae* | 3 | 0 | cefotaxime (101), ceftazidime (102), sulfamethoxazole/trimethoprim (101) |
| *Mycobacterium tuberculosis* | 20 | 1 | amikacin (1,813), amoxicillin (57), capreomycin (670), ciprofloxacin (101), cycloserine (401), delamanid (1,248), ethambutol (1,947), ethionamide (1,774), kanamycin (1,800), levofloxacin (1,248), linezolid (1,248), moxifloxacin (1,524), nicotinamide (166), ofloxacin (643), para-aminosalicylic acid (1,650), prothionamide (170), pyrazinamide (492), rifabutin (1,317), rifampin (1,967), streptomycin (698) |
| *Neisseria gonorrhoeae* | 6 | 0 | cefixime (600), ceftriaxone (1,926), ciprofloxacin (1,746), erythromycin (179), penicillin (1,354), tetracycline (1,371) |
| *Pseudomonas aeruginosa* | 12 | 6 | amikacin (708), aztreonam (235), cefepime (245), ceftazidime (833), ceftazidime/avibactam (309), ceftolozane/tazobactam (259), gentamicin (346), imipenem (447), levofloxacin (702), meropenem (883), piperacillin/tazobactam (210), tobramycin (645) |
| *Salmonella enterica* | 24 | 1 | amoxicillin/clavulanate (799), ampicillin (1,997), ampicillin/sulbactam (183), azithromycin (1,355), cefazolin (185), cefotaxime (941), cefoxitin (1,080), ceftazidime (783), ceftiofur (689), ceftriaxone (1,155), chloramphenicol (1,989), ciprofloxacin (1,503), gentamicin (1,675), kanamycin (584), nalidixic acid (1,599), pefloxacin (348), spectinomycin (396), streptomycin (1,269), sulfamethoxazole (275), sulfamethoxazole/trimethoprim (1,145), sulfathiazole (225), sulfisoxazole (544), tetracycline (1,869), trimethoprim (634) |
| *Shigella flexneri* | 4 | 0 | ampicillin (358), azithromycin (361), ciprofloxacin (361), trimethoprim (362) |
| *Shigella sonnei* | 20 | 0 | amikacin (580), ampicillin (1,674), azithromycin (1,307), cefoperazone (257), cefotaxime (593), cefoxitin (593), cefpirome (333), ceftazidime (582), ceftriaxone (1,163), ciprofloxacin (1,397), colistin (268), ertapenem (338), gentamicin (1,288), nalidixic acid (854), streptomycin (595), sulfamethoxazole/trimethoprim (617), temocillin (333), tetracycline (878), tobramycin (541), trimethoprim (1,104) |
| *Staphylococcus aureus* | 13 | 4 | cefoxitin (1,407), chloramphenicol (1,038), ciprofloxacin (1,624), clindamycin (1,543), daptomycin (904), erythromycin (1,957), fusidic acid (1,101), gentamicin (1,913), mupirocin (736), penicillin (1,508), rifampin (931), sulfamethoxazole/trimethoprim (1,389), tetracycline (1,647) |
| *Staphylococcus haemolyticus* | 3 | 0 | ciprofloxacin (124), fusidic acid (114), tetracycline (104) |
| *Streptococcus agalactiae* | 3 | 0 | clindamycin (141), erythromycin (141), tetracycline (141) |
| *Streptococcus pneumoniae* | 8 | 0 | azithromycin (1,297), cefuroxime (1,732), chloramphenicol (1,859), clindamycin (1,894), erythromycin (1,991), meropenem (1,808), sulfamethoxazole/trimethoprim (1,777), tetracycline (1,971) |
| *Vibrio cholerae* | 3 | 0 | nalidixic acid (91), sulfamethoxazole (74), trimethoprim (74) |
| All | 233 | 55 | |
