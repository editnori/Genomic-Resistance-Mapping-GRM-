# Outside learners: XGBoost, LightGBM and an elastic net on GRM's own folds

Plan block C7 (`docs/assumption-plan.md`) asks whether learners from outside GRM do better than GRM's own on the same folds. Three arms stand against the release's rule as plan block C3's run refitted it, read on each species' C3 machine: XGBoost 3.0.5 and LightGBM 4.6.0 on the 4,000 rows `grm export` screens for each outer fold, read against GRM's boosted trees (`gbt`), and a scikit-learn 1.7.2 elastic net on the unitig build of plan block C3, read against GRM's lasso (`lasso`). The arms are scripts under `scripts/c7-outside/`. The one change under `src/` is `grm export --seed` (below): with no flag, and at `--seed 42`, every byte the export writes is what it wrote before, and no other command changes. No comparison result has been read. The lead ruled the design's 14 questions at 09:05Z on 23 Nov 2026; the "C7 design" entry of `docs/assumption-plan-corrections.md` lists each ruling.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation.

## The arms

| Arm | Library and call | Input per outer fold | Read against |
|---|---|---|---|
| default | GRM's five learners, `grm learn --exp-pick DIR --exp-pick-inner 0`, the C3 run's at seeds 42, 43 and 44 | K31 | itself |
| `xgboost` | XGBoost 3.0.5, `xgboost.train`, `tree_method` exact | the export's 4,000 screened rows of that fold | GRM's boosted trees |
| `lightgbm` | LightGBM 4.6.0, `lightgbm.train` | the same rows | GRM's boosted trees |
| `elasticnet` | scikit-learn 1.7.2, `LogisticRegression(penalty="elasticnet", solver="saga")` | the C3 unitig build projected on the cohort, then the chi-square top 4,000 on each fit's own genomes | GRM's lasso |
| twin | the default at seeds 45, 46 and 47, the C3 run's | K31 | the default |

Each wheel is pinned by its sha256 on pypi.org in `scripts/c7-outside/requirements.txt`, which takes the learner cross-check's file (numpy 1.26.4, scipy 1.15.3, scikit-learn 1.7.2, joblib 1.5.2, threadpoolctl 3.6.0, six 1.17.0, XGBoost 3.0.5) and adds LightGBM 4.6.0 (`sha256:cb19b5af…1c8d`). The runner reads each library's version back and refuses another.

**Boosters.** Both take GRM's settings where GRM has them (`GbtParams::default`): 400 rounds at most, learning rate 0.1, depth 4, L2 weight 1, a hessian floor of 1 per leaf, half the screened rows for each tree, and early stopping after 25 rounds with no gain in log loss on the plan's 15% boosting holdout. XGBoost takes the rest from the learner cross-check's map (`scripts/learner-crosscheck/check_gbt.py`, `PARAMS`), which the runner holds to GRM's values. LightGBM keeps its own defaults where GRM has no knob, with 16 leaves, no start from the class share, `deterministic` and one thread, except one: `min_data_in_leaf` is 1, not the library's 20, since cohorts go down to 25 genomes a class and a fitting side of about 34 genomes, where a floor of 20 allows no split. `min_data_in_bin` stays at the library's 3: a SYNTHETIC test shows a 0/1 row whose fitting-side support is 1 or 2 is still split on (`test_runner.Floors`). Each fit's sample weights are GRM's: balanced on the fitting genomes. The seed is the fold's seed (`seed + 101 × fold`) mod 2^31. A booster that keeps one tree is scored as the library leaves it and counted in `one_tree_folds`.

**Elastic net.** The class-weighted mean log loss plus lambda times ((1 − mix) / 2 × |b|² + mix × |b|₁), each column divided by its standard deviation on the fitting genomes, as GRM's lasso scales its penalty. saga is given the scaled columns centred, which moves only the unpenalized intercept, and the intercept is moved back. For each mix of 0.25, 0.5, 0.75 and 1, 10 penalties run log-spaced from lambda_max(mix) down to 1e-3 of it: 40 points from the outer training genomes. lambda_max is the smallest penalty that keeps every weight at 0; there the fit is written down as all zeros, since saga's relative stopping rule never ends on it. The point with the lowest class-balanced log loss pooled over the plan's 3 inner folds is refitted on every training genome; an unconverged point is left out, and a tie goes to the larger penalty, then the larger mix. saga stops at tol 1e-4 or 2,000 epochs. A refit that reaches 2,000 is scored and counted in `unconverged_refits`; a refit with no nonzero weight is GRM's null fit and is counted in `null_fit_folds`, and, scored at the training share, in `constant_folds` too. A fold where none of the 40 grid points converged has no pick and no refit: it is written as a null fit, scored at the training share and counted in `null_fit_folds`, `constant_folds` and `unconverged_refits` alike, so the count's no-unconverged table leaves its dataset out.

## Inside the folds

The folds are GRM's. On the lineage and random splits the runner reads the plan `grm export --seed s` wrote for the draw's seed (`dataset.json`, `splits[].plan`): each outer fold's training and test genomes, the boosting holdout and the 3 inner folds, with the plan's sha256 checked over its canonical JSON. On the study split it reads the reviewed `study-whole` fold file's outer folds and draws each fold's holdout and inner folds by GRM's rule, duplicate groups whole, with numpy's generator seeded from the sha256 of the dataset, the scheme and the fold. It is GRM's rule there, not GRM's genomes.

A booster reads the rows the export screened for that fold (`screens[f].screened_rows`): GRM's forest screen on all of the fold's training genomes. On the study split, where the export writes no screen, the runner computes the same screen by the export's rule. The elastic net projects the unitig build on the cohort at the dataset's floor, as `Dataset::project` projects K31, then takes the candidates on the outer training genomes and the chi-square top 4,000 on each fit's own genomes: each inner fold's training part for the grid, every training genome for the refit. No outer test genome enters a screen, a holdout, an inner fold or a pick, and a test holds each of those apart.

A fold whose training genomes give no screened row or no candidate scores every test genome at the training share of resistant genomes, GRM's null fit and ruling 16's constant model, counted in `constant_folds`.

## Deviations and limits

The lead ruled these, and each is a departure from GRM's own fit:

- **Study plans.** On the study split the outer folds are the reviewed file's, and the inner folds and the boosting holdout are drawn by GRM's rule: they are other genomes than GRM's own study fit drew.
- **The boosters' rows.** The rows are the export's forest screen on all of a fold's training genomes. The boosters' stopping holdout took part in choosing those 4,000 rows, and GRM's `gbt` holdout did not take part in its own screen, so a difference between a booster and `gbt` reads the screen and the library together.
- **LightGBM's floors.** `min_data_in_leaf` 1 in place of the library's 20 (above). `min_data_in_bin` is the library's 3; the SYNTHETIC test shows it does not block a rare row.
- **Default and twin.** C7 reads the default and the twin of the C3 run on each species' C3 machine, and does not fit them again. `grm` only exports in C7; the C3 run fitted the default and the twin on the same 233 datasets and held seed 42 to the release; fitting them again would cost 33.45 box-hours central and 93.11 upper. The run checks and refuses on each of: every export's lineage plan at seed s equals the default's folds at seed s by class, and by genome (each fold's genome count and the sha256 of its sorted genome ids, against the grm learn report the C3 run wrote at seed s); the export's input fingerprint at seed 42 equals the release report's; the species' C3 run ended `DONE`, built at `8bfd4ca`, with its `SHA256SUMS` holding. C7 opens nothing under a C3 run's `score/`.

One limit is measured and not ruled. saga's tol of 1e-4 stops on the change in the weights and bounds no distance to the optimum. On the 12 two-class panels of `tests/data/lasso-stopping-fixtures.txt`, saga at mix 1, given the same columns and GRM's lasso penalties, came within 1e-4 of GRM's lasso objective at 49 of the 76 points where both report convergence. It came up to 10.4% above it at the smallest penalties of the Streptococcus agalactiae clindamycin, Enterobacter cloacae tobramycin and Staphylococcus haemolyticus fusidic acid panels, and never below it. A difference between the elastic net and `lasso` reads saga's stop as well as the learner.

The C3 runs were built at `8bfd4ca`, the export at this branch's head. `git diff 8bfd4ca..HEAD -- src/` changes:

| Change | Files | Proof on main that a run without the flag writes the same bytes |
|---|---|---|
| C4's `grm learn --exp-screen` and `--exp-floor` | `learn/screen_arm.rs`, `learn/eval.rs`, `analysis.rs`, `main.rs` | "C4 design", unflagged bytes: on the pilot's root 360 files, 317 equal in bytes and 43 equal but for the run clocks, none differing |
| C5's `grm learn --exp-settings` and `--exp-settings-folds`, and `cutoff::calls_at` | `learn/settings.rs`, `cutoff.rs`, `learn/eval.rs`, `main.rs` | "C5 design", bytes: 6 of 6 files equal once the recorded seconds are set to 0, the seed-42 report equal to the release's in 12 of 12 learner rows |
| C7's `grm export --seed` | `main.rs` | a test pins the export of a SYNTHETIC species to the bytes main wrote at `0d987ba`; on the pilot the head's export at `--seed 42` and with no flag equals `8bfd4ca`'s in 12 of 12 files |

## The files

`scripts/c7-outside/outside.py fit` writes one report per arm, dataset, split and seed, `results/outside/<arm>/<dataset>.<scheme>-s<seed>.outside.json`, format `grm-c7-outside-1`. It holds the revision, the dataset, split, seed, arm, library and version, the settings above, the plan's sha256 and the export's input fingerprint and matrix sha256, the cohort and tested genomes by class, and per outer fold: the training, test, fitting and holdout genomes by class, the rows offered and screened, the rounds run and kept or the grid point picked (its mix and its index), the nonzero weights as a count, convergence counts, and the calls on its test genomes (`n`, `tp`, `fp`, `tn`, `fn`, `loss_sum`) at the 0.5 cutoff, with log loss clipped to [1e-4, 1 − 1e-4] as `src/cutoff.rs` clips it. Every key is on a fixed list, and the runner refuses a report with another. A report holds no genome id, no row per genome, no row index, no feature, k-mer, unitig or sequence, no model, and no weight, gain or importance. Its time goes to a file beside it, so two runs write the same bytes.

## The scorer

`scripts/c7-outside/score.py` checks every result `scripts/c7-outside/expect.py` lists and stops before scoring on any failure. The default and the twin are the C3 stage's results, held as `scripts/c3-inputs/score.py` holds them, the default's seed-42 lineage figures to the release report's. Each arm's report must hold its listed keys only; its format, revision, dataset, split, seed, arm, library and version; the GRM learner its pair names; the default's folds at the same dataset, split and seed, by class; and each fold's calls equal to its test genomes. A missing draw is a failed check, never a zero.

It scores each row on its outer calls summed over the folds: balanced accuracy, recall, specificity, the very major and the major error rate, each with its count, denominator and Wilson 95% interval, and log loss (`scripts/c9-pick/score.py`, `measures`). `datasets.tsv` gives each dataset's mean over its draws, and the default again as `k31-at-arm-seeds`, the mean over the draws at the seeds the arms ran, which each arm's pair reads.

`score.py --count` takes several species' tables once all three arms and the twin have lineage rows, or an arm is named with `--withdrawn`, and refuses before that. Five pairs make the family, withdrawn or not:

| Pair | Arm's row | Against |
|---|---|---|
| `xgboost-gbt` | `xgboost` | the default's `gbt`, at the arm's seeds |
| `lightgbm-gbt` | `lightgbm` | the default's `gbt`, at the arm's seeds |
| `elasticnet-lasso` | `elasticnet` | the default's `lasso`, at the arm's seeds |
| `twin-gbt` | the twin's `gbt` | the default's `gbt` |
| `twin-lasso` | the twin's `lasso` | the default's `lasso` |

Per pair it counts the datasets better, worse and the same by 0.02, the arm's value less the default's at 12 decimals, lower being better for log loss and both error rates, with the median difference and its interval from 4,000 species resamples and the one-sided sign test. Holm runs over the 5 lineage pairs on balanced accuracy. Each arm touches one learner, so an arm carries the change when its pair meets condition 1 (interval wholly above 0) and condition 2 (better beats worse at p < 0.05 after Holm). Condition 3 holds when the twin meets 1 in neither pair and 2 in neither pair. Condition 4: the arm's recall interval does not lie wholly below −0.02. Condition 5: its study median has its lineage median's sign. Condition 6 is the maintainer's. `count-no-constant.tsv` is the count without the datasets that hold a constant fold, with Holm over its own rows, and decides nothing; `count-no-unconverged.tsv` is the same without the datasets whose elastic net draws hold a refit that reached 2,000 epochs. Beside condition 3, `conditions.tsv` reads each arm against its own learner's twin pair alone (`decides` 0), which decides nothing.

An arm withdrawn before it runs, `--withdrawn elasticnet` say, keeps its pair in the family at p = 1 and is named in `count.tsv`'s `withdrawn` column. The count refuses a withdrawn arm that has results, and refuses the twin or the default named as withdrawn.

## Datasets

`scripts/c7-outside/lists.py` writes `scripts/c7-outside/lists.tsv` and `species.tsv` from C3's list, `scripts/c3-inputs/lists.tsv`, with each species and antibiotic in full.

| Set | Datasets | Species | Less M, the first run | First run in S |
|---|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 |

The first run is L less M; its 55 datasets in S also get one study run per arm at seed 42, with no twin.

## Running

`scripts/c7-outside/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARM_SEEDS] [ENET_JOBS] [ENET]` runs one species with no prompt, from a clean checkout. On the species' C3 machine `OUT/c3` is a link to the C3 run's folder: nothing of stage 1 runs again, and the run stops unless that run ended `DONE`, built at `8bfd4ca`, with the `SHA256SUMS` its `DONE` names and every line of it outside `score/` holding. On a new box (`OUT/c3` not a link) the first stage is C3's run as it is, `scripts/c3-inputs/run_species.sh KIT SPECIES DATASETS OUT/c3 JOBS SEEDS TWIN unitigs "" grm`, kept for a lost C3 machine.

It then installs the pinned wheels in a venv by sha256, builds `grm` at the branch head for the export alone (the C3 run's binary has no `--seed`), and runs `grm export --seed s` with no network for each seed of `ARM_SEEDS`, 42, 43 and 44 by default. `scripts/c7-outside/plans.py` checks each export's lineage plan at seed s against the default's folds at seed s by class and by genome (the C3 run's report at seed s, under `results/` at 42 and `draws/seed-<s>/` otherwise; only a pass or a fail is printed), and the seed-42 input fingerprint against the release report's, and the run stops before any fit if one check fails. It fits each arm on every listed dataset, split and seed with one thread under GNU time (the boosters `JOBS` at once, then the elastic net `ENET_JOBS` at once), scores into `score/` (the scorer's output too, as `score/score.log`), scans `results/`, `score/` and `logs/` for any genome id of the species and writes `DONE` with the sha256 of its list of output files. Each fit's log comes back with every genome id written as `[genome id]`. Run again, it skips finished steps and finished results. Any failure writes `OUT/FAILED` naming the step.

`ENET` is `run` by default; `leave-out` lists and fits no elastic net, and the count then takes `--withdrawn elasticnet`. A run started with `leave-out` restarts on the same `OUT` with `run`: it keeps the finished steps and the boosters' reports, fits the elastic net, scores again and writes a new `DONE`. Any other change to the datasets, the seeds or the arms after the first start, `run` back to `leave-out` among them, is refused and needs a fresh `OUT`.

An elastic net fit holds its cohort's projection of the unitig build, which peaks at about 5.2 times the build's `patterns.bin`: 4,364,420 KB for a SYNTHETIC build of Escherichia coli's size, 3,351,087 patterns over 2,000 genomes. Unless `ENET_JOBS` is given, the script runs as many at once as three quarters of the box's memory holds at 6 times `patterns.bin` each, from 1 to `JOBS`: 4 on a 32 GB box at Escherichia coli's size.

The pilot, on 23 Nov 2026 on one 16-CPU, 32 GB box: Klebsiella aerogenes with ertapenem (83 genomes, 48 resistant), imipenem (68, 26) and meropenem (78, 26), the lineage split, one thread a fit. It ran first at seed 42: the default (9.50 s wall at 16 threads), the export (3.33 s), the unitig build (8.90 s) and each arm. It ran again after the rulings at seeds 42, 43 and 44: `grm` built at the branch head (42.38 s), the export at each seed (3.54, 3.41 and 3.55 s), the default at seeds 43 and 44 as the C3 run ran it (13.26 and 11.72 s at 16 threads), `plans.py` (21 of 21 checks pass over 9 dataset seeds), the boosters at the three seeds with LightGBM's floor at 1, and the elastic net at 43 and 44. The head's export at `--seed 42` and with no flag equals `8bfd4ca`'s in 12 of 12 files; at seeds 43 and 44 each `dataset.json` differs from seed 42's (6 of 6) and each `matrix.bin` equals it (6 of 6). The reports went to sealed folders that no one opened; what was read is that each fit ran, the input counts, the genome-id scan (0 ids in the 47 files of the first run, 0 in the elastic net rerun's 12, and 0 in the 119 of the three-seed run) and the time. The runner's screen equals the export's on 5 of 5 lineage folds in each dataset at seed 42, candidates and screened rows alike. The export held 539,470, 413,712 and 501,572 rows, and the unitig build 577,629 patterns, of them 539,376, 413,610 and 501,500 on each cohort.

CPU seconds a draw at seeds 42 / 43 / 44, the elastic net's at 42 from the first run's rerun:

| Arm | ertapenem | imipenem | meropenem | Peak memory (KB) |
|---|---|---|---|---|
| `xgboost` | 0.72 / 1.07 / 0.99 | 0.87 / 0.90 / 0.76 | 0.70 / 0.94 / 0.64 | 129,552 to 133,536 |
| `lightgbm`, `min_data_in_leaf` 1 | 2.05 / 2.57 / 2.17 | 2.10 / 1.90 / 2.14 | 2.17 / 2.29 / 2.13 | 141,232 to 143,516 |
| `elasticnet` | 667.37 / 600.63 / 626.61 | 600.08 / 576.22 / 608.25 | 686.13 / 723.63 / 712.24 | 243,080 to 318,480 |

The elastic net carries the cost. At 2,000 genomes a SYNTHETIC fold over 2,000,000 patterns took 603.8 s, its screen 13.3 s in it; a power law through it and the 9 pilot fits' 128.91 s a fold at a mean of 76.33 genomes has exponent ln(603.8 / 128.91) / ln(2000 / 76.33) = 0.4728, and a draw of 5 folds at 2,000 genomes is 3,019.0 s, plus 82.3 s to project Escherichia coli's build on 2,000 genomes. The boosters stay small: 13.40 CPU s a draw for XGBoost and 7.85 for LightGBM at 2,000 genomes, LightGBM's from a SYNTHETIC fold at the old floor of 20.

Over the 233 datasets at three lineage draws and the 55 study runs, on the species' C3 machines, C7's steps take 74.25 box-hours central and 250.54 upper ($14.85 and $50.11 at xlarge's $0.20 an hour); the elastic net is 70.44 and 244.17 of that, and without it the boosters, the exports, the build and the libraries take 3.81 and 6.38 ($0.76 and $1.28). On a new box a species, with stage 1 there too, the run comes to 107.69 central and 343.65 upper, or 37.26 and 99.49 without the elastic net. The plan gives 15 xlarge-hours plus the elastic net after a pilot.

## Tests and planted faults

`tests/c7_outside` holds 83 SYNTHETIC tests: 32 of the runner, 30 of the scorer, 17 of the run script and the plans check, and 4 of the lists and the scan. Two Rust tests in `src/main.rs` hold `grm export --seed`: with no flag and at `--seed 42` the files are the bytes main wrote at `0d987ba`, pinned by sha256, and at `--seed 43` the plans and screens are the ones `grm learn --seed 43` scores, fold by fold, with no key added. XGBoost and LightGBM are each held to a plain booster written apart in `tests/c7_outside/plain.py`, with each library's own rules (XGBoost's exact gain, LightGBM's leaf count read from the hessian and its pre-filter of rare rows), and saga to the stated objective solved apart by FISTA. The suite needs the pinned wheels, and none of its tests is skipped without them. `scripts/c7-outside/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must pass on the clean tree and fail with the fault in.

| Fault | Caught by |
|---|---|
| 1. A booster that fits on the outer test genomes too | `Arms.test_no_outer_test_label_moves_an_xgboost_score` |
| 2. A study fold's boosting holdout drawn over its outer test genomes too | `Plans.test_study_holdouts_and_inner_folds_lie_inside_each_folds_training_genomes` |
| 3. Early stopping read on the outer test fold | `Arms.test_no_outer_test_label_moves_a_lightgbm_score` |
| 4. The elastic net's mix and penalty picked on the outer test fold | `Arms.test_no_outer_test_label_moves_an_elastic_net_score` |
| 5. A study report read where a lineage one belongs | `Checks.test_a_study_report_read_as_a_lineage_one_fails` |
| 6. The elastic net read against GRM's boosted trees, not its lasso | `Checks.test_each_arm_pairs_with_its_own_grm_learner` |
| 7. The twin scored as an arm | `Count.test_conditions_are_read_for_the_three_arms_and_never_for_the_twin` |
| 8. A missing draw read as a zero | `Checks.test_a_missing_draw_fails_and_is_never_read_as_a_zero` |
| 9. A report that carries feature names | `Arms.test_the_report_holds_only_its_listed_keys_and_no_genome` |
| 10. `grm export --seed S` writing seed 42's plans | `tests::export_at_a_seed_writes_the_plans_learn_scores_at_that_seed` (Rust) |
| 11. An export whose folds differ from the C3 default's by class let through | `PlansCheck.test_a_plan_that_differs_from_the_defaults_by_class_is_refused` |
| 12. A withdrawn arm with results counted | `Count.test_a_withdrawn_arm_with_results_is_refused` |
| 13. LightGBM's floor of 20 genomes a leaf back | `Floors.test_a_fitting_side_of_34_genomes_is_split` |

## What it does not do

It reads no sealed, held-out, clinical or spent test set: neither the outside-study genomes nor the other-database genomes, nor the later-records file. It finds, names and ranks no resistance marker: no feature, k-mer, unitig, model, weight, gain or importance leaves a fit. It changes no report and no release rule, and under `src/` only adds `grm export --seed`.
