# Settings arms: a nested search of each learner's settings, and a call cutoff chosen on inner scores

An experiment, off by default, for plan block C5 (`docs/assumption-plan.md`). The release's rule already chooses some settings inside each fit. C5 asks whether a wider search, made inside the training folds, does better: per learner, 30 random settings from the spaces frozen in `scripts/hyperopt/protocol.py`, picked by inner 3-fold balanced accuracy; and, as a second arm, the call cutoff chosen on inner scores in place of 0.5. This branch builds both arms, the driver that runs them, the scripts that list, run and score one species, and a timing pilot. No comparison has been read.

`grm learn --exp-settings FILE` and `grm learn --exp-settings-folds DIR` are the two new flags. Without them nothing new is called and `grm learn` is unchanged, byte for byte. `EVALUATION_REVISION` stays 14 and the report schema stays 2. Neither flag writes a report or removes one.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation.

## The arms

| Arm | On outer fold f of a split |
|---|---|
| Default | The release's rule at seed s, refitted in the same run with the same binary: `grm learn --exp-pick DIR --exp-pick-inner 0`, each outer fold's calls per learner |
| Search | Per learner, the setting the rule below picks among 31 trials, each scored by inner 3-fold balanced accuracy on fold f's training genomes T_f; then that setting refitted once on the whole of T_f and scored once on fold f |
| Cutoff | The release's own fit on T_f, unchanged, with its calls made at a cutoff t_f chosen on the release setting's inner scores in place of 0.5 |
| Twin | The default at seeds 45, 46 and 47, with no search |

Trial 0 is the release's own setting. Trials 1 to 30 are `scripts/hyperopt/search.py`'s `sample_config` at the protocol's sampler seed 20261030, so every dataset, split, draw and fold has the same 30. A trial fixes every value of its learner's space, as `protocol.py`'s contract has it. The lasso's space is one integer, the length of its penalty path, so its 30 draws hold 28 distinct settings. Trial 0 is the default's own fit: on every outer fold the driver checks that the release setting's calls at 0.5 equal the default's, learner by learner, and stops if they differ.

What the release already chooses inside each fit, so no value is searched twice:

| Learner | The release chooses | A sampled setting fixes |
|---|---|---|
| scm | the kind, p from 7 values and 1 to 6 rules, by inner 3 folds | the kind, p in [0.1, 10] on a log scale and 1 to 10 rules, with no inner search |
| cart | depth 1 to 8 by inner 3 folds, leaf 1 | depth 1 to 30 and leaf 1 to 20, with no depth search |
| rf | nothing but the screen to 4,000 | trees, subsample, depth and leaf; the screen stays |
| gbt | rounds by early stopping on its holdout | learning rate, depth, lambda, child weight and column share; rounds stay the cap early stopping runs under |
| lasso | the penalty from a 30-point path on its holdout | the path's length, 10 to 100 points; the holdout still picks the penalty on it |

The screen's size, the support floor and the inner fold count are not searched: they are block C4's.

## The selection rule

The default-anchored rule of `docs/hyperopt-selection.md` at γ = 0.04: a learner keeps trial 0 unless the best inner mean beats trial 0's by more than 0.04, and then takes the best. A trial's inner mean is the mean of its three inner folds' balanced accuracies, compared at 12 decimals. A tie at the best goes to the lowest trial id. A trial whose fit failed is left out; if trial 0 failed, the fold keeps trial 0 and the choices file says so. The plain best is written beside the pick and never refitted.

The cutoff t_f is a point of the grid 0.05, 0.10, ..., 0.95 (19 points, 0.5 among them). Its objective is the balanced accuracy of trial 0's inner scores on T_f, pooled over the three inner folds, a score at or above t called resistant. A tie goes to the point nearest 0.5, and between two equally near, to the lower one. The fitted model is not touched: the outer run fits trial 0 once and writes its calls at 0.5 and at t_f from the same scores. scm scores 0 or 1, so every point gives it the same calls, the tie gives 0.5, and its cutoff arm equals the default on every fold. The scores do not move, so the cutoff arm's log loss is the default's on every dataset: its log loss rows count the same everywhere and decide nothing.

## Inside the folds

The inner folds of the lineage split are the release's own tuning partition of T_f, `training_folds` over lineage groups at `fold_seed(s, f)`, the partition the release's scm and cart tune on in fold f. Each lineage cluster lies in one inner fold. On the study split the same call runs inside each study training fold. That keeps every lineage whole but can split a BioProject across inner folds, since `grm learn` holds no BioProject per genome.

Each inner fit is a fold file fit, so its candidates come from its own inner training genomes and its chi-square screen to 4,000 from its own fitting genomes. The outer test genomes are listed as unused in the inner fold file: they reach no fit and are not scored. The outer run starts only after the picks and the cutoffs are written. Inside each fit the release's own tuning keeps lineages whole on the lineage split and copy groups whole on the study split, as the release's fits do there. The outer refit runs at the seed the default fitted that fold with, `s + 101 f`.

## The flags

- `--exp-settings FILE` (with `--fold-file`): fits the trials of FILE, format `grm-c5-settings-1`, on the fold file's partition, each trial's learners at their own settings, and writes per trial `<fold-out>/<dataset>.<scheme>.t<trial>.json`, format `grm-c5-trial-1`: the fold file result's fields for the learners fitted, the trial's settings, and each learner's calls per fold at 0.5 and at every cutoff FILE gives. FILE also gives the fits' base seed and the groups the fits' own tuning keeps whole, `lineage` or `copies`. A learner whose fit fails is written as failed with its reason, and the others go on.
- `--exp-settings-folds DIR`: after a dataset's evaluation, writes `DIR/<dataset>.<scheme>-s<seed>.folds.json`, format `grm-c5-folds-1`: each genome's outer fold and, per outer fold, its inner fold inside T_f, with each genome's lineage and class. It changes no fit, score or report; it stands beside `--exp-pick`.

Both need `--force`, so a kept report can never skip them. Each refuses beside the other and beside `--learner-seed`, the `--exp-tie-*` options, `--exp-null-refit`, `--exp-mic`, `--exp-cohort` and `--exp-inputs`; `--exp-settings` also refuses `--exp-pick`. Both files hold genome ids, as the reviewed study plans do, and stay on the machine that made them.

## The driver

`scripts/c5-settings/drive.py` (Python standard library) runs one dataset, split, draw and outer fold, for the five learners, each searched and picked apart:

1. It checks the folds file: every genome of T_f in one inner fold, no other genome in any, no lineage in two inner folds, and both classes in each inner fold's training and test genomes.
2. It writes the inner fold file and, per learner, the settings file of its 31 trials at the inner seed (the first 8 bytes of sha256 of `c5-inner/<s>/<f>`, shifted right by 1), and runs `grm learn --fold-file --exp-settings` once per learner, so each learner's time and memory are its own.
3. Per learner it picks a trial by the rule and a cutoff on trial 0's inner calls, and writes `pick.json`.
4. Only then it writes the outer fold file (test fold f, training T_f) and runs, per learner, the release's setting and the picked one on it.
5. It writes `<dataset>.<scheme>-s<seed>.f<f>.choices.json`, format `grm-c5-choices-1`: per learner the trial picked and its setting, the rule's counts (trials better, worse and the same as trial 0, above γ, tied at the best), the cutoff, and the outer calls of the default, the search and the cutoff arm, as counts. It checks the file names no genome before it writes it.

A restart skips every step whose files exist. Everything under its work folder may name genomes and stays on the machine.

## The scorer

`scripts/c5-settings/score.py` is block C3's scorer (`scripts/c3-inputs/score.py`), loaded as a module of its own and set to C5's family: the search, the cutoff and the twin, each against the default. C3's code checks the default's and the twin's results, the default's seed-42 results against the release report included, and scores every draw. C5's own code checks each choices file before it scores: its format, revision, dataset, split, seed and fold; the sampler seed, γ, grid and 30 trials; the five learners in order; the fold's genomes by class equal to the default's fold; the default's calls equal to the default's own result (the default is never refitted with a searched setting); each arm's calls of the fold's test genomes; the pick by the rule; the cutoff a grid point or 0.5; and scm's cutoff arm equal to the default. An arm at a twin seed fails. A draw with a missing choices file is not scored and its dataset's mean is not taken: a missing draw is never read as zero. Any failed check stops the scorer before it scores.

Per split, arm and learner, against the default on the same datasets and folds, the count gives: better, worse and the same by 0.02, the arm's value less the default's, lower better for log loss and both error rates; the median difference with its interval from 4,000 species resamples; the one-sided sign test of better against worse; and Holm over the 15 lineage rows of the primary measure (search, cutoff and the twin by five learners). A dataset's value is the mean of its draws. The conditions are the Rules' six, read with the lead's C3 rulings 15 and 16 and condition 3 as corrected on 23 Nov 2026: the twin meets condition 1 in fewer than 3 of 5 learners and condition 2 in fewer than 3 of 5, each read apart. SCM's result is read beside the synthetic loss of 0.0126 (0.0088 to 0.0163) the plan names for the default-anchored rule.

## Datasets

`scripts/c5-settings/lists.tsv` holds the 69 datasets in both S and L, 9 species. The first run takes the 55 that are not in M, in 8 species: *Acinetobacter baumannii* 2, *Enterobacter cloacae* 8, *Escherichia coli* 15, *Klebsiella pneumoniae* 18, *Mycobacterium tuberculosis* 1, *Pseudomonas aeruginosa* 6, *Salmonella enterica* 1 and *Staphylococcus aureus* 4. Each of the 55 has a reviewed study plan. The 14 in M wait for the maintainer's label rule; *Enterococcus faecium*'s one dataset in S is among them.

## Running

`scripts/c5-settings/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [DRIVES]` runs one species on a box of its own, with no prompt. It checks its kit by sha256 (the release, the data root's list, the curated table, the study plans and the label count's table), refuses a dataset outside S or in M, and refuses at the kit step a box without numpy or a tool a later step needs. It builds `grm` and the species' data root as `scripts/c3-inputs/run_species.sh` does. With one binary it fits the default at each seed with its folds files, checks the seed-42 reports against the release's with `scripts/refit/refit_check.py`, fits the twin and the default's study runs, then runs the driver on every outer fold of every draw and study result, `DRIVES` at a time (8 by default on 16 CPUs). It scores, writes `logs/timing.tsv` and `logs/drive-timing.tsv` (each driver step's wall, user and system seconds and peak memory), and writes `OUT/DONE` with the sha256 of its list of output files. Each log is copied with every genome id written as `[genome id]`, and the end fails if a file under `results/`, `check/`, `score/` or `logs/` names one. Run again, it skips finished steps and finished drives. Any failure writes `OUT/FAILED` and exits 1.

The pilot ran on one 16-CPU, 32 GB box on 23 Nov 2026: *Pseudomonas aeruginosa* with piperacillin/tazobactam (210 genomes), then with aztreonam (235 genomes), seed 42, both splits, five learners, the default and both arms, 8 drives at once with 2 threads each. The ten outer folds of piperacillin/tazobactam took 8 min 47 s from first start to last end and those of aztreonam 9 min 31 s. A drive took 274.05 s on average (311.68 CPU seconds), and every drive's peak was under 1,890,000 KB. Per outer fold, a learner's 31 inner trials took 40 to 57 s on average and its outer refit 7 to 11 s; scm took longest. Every one of the 200 `grm learn` runs ended 0, the 20 choices files have one shape, and none names a genome. The entry "C5 design" in `docs/assumption-plan-corrections.md` gives the cost of the first run.

The run script was rehearsed on the same box with *Pseudomonas aeruginosa* and piperacillin/tazobactam alone, seed 42 and twin seed 45, from the genomes of the release's refit: 19 min 20 s from start to `DONE`, the steps before the fits (the build and the data root) 7 min 2 s of it and the 10 drives 10 min 33 s. The seed-42 report equalled the release's in 12 of 12 learner rows, the scorer's 323 checks all passed with 7 of 7 results present, and no file it keeps names a genome. Started again on the finished folder, it skipped its 12 finished steps and 10 finished drives and wrote the same list of output files, byte for byte; started with another seed list, it refused at the kit step and left `DONE` in place. No score was opened.

## Tests and planted faults

`src/learn/settings.rs` holds 8 SYNTHETIC tests and `src/main.rs` one more on the flags. `tests/c5_settings` holds 30 SYNTHETIC tests of the driver (against a fake `grm`), the scorer, the lists and the expected results. `scripts/c5-settings/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must fail.

| Fault | Caught by |
|---|---|
| 1. A setting picked on the outer test fold: the inner run tests it as its third inner fold | `Drive.test_no_outer_test_genome_is_read_before_the_pick` |
| 2. An inner fold that splits a lineage cluster: drawn by genome in `grm` (2a), or taken by the driver (2b) | `inner_folds_keep_lineages_whole_inside_the_outer_training_genomes`; `Drive.test_a_folds_file_that_splits_a_lineage_is_refused` |
| 3. The screen fitted on the fold's test genomes too | `no_test_genomes_label_reaches_its_folds_fit` |
| 4. The cutoff chosen on outer scores | `Drive.test_the_cutoff_is_chosen_on_inner_scores_alone` |
| 5. A tie broken by the outer score | `Drive.test_a_tie_goes_to_the_lowest_trial_by_inner_scores_alone` |
| 6. The default refitted with a searched setting, in the driver (6a) or taken by the scorer (6b) | `Drive.test_the_default_arm_is_the_defaults_own_fit`; `Score.test_a_choices_file_whose_default_is_not_the_defaults_fit_fails` |
| 7. The twin scored as an arm, in the list (7a) or the scorer (7b) | `Expect.test_no_arm_runs_at_a_twin_seed`; `Score.test_an_arm_at_a_twin_seed_fails` |
| 8. A missing draw read as a zero | `Score.test_a_missing_draw_is_not_scored_and_never_read_as_zero` |
| 9. Holm over the two arms alone, the twin left out | `Count.test_holm_runs_over_the_fifteen_lineage_rows_of_the_primary_measure` |
| 10. A trial fitted at another seed than the default's fold | `the_release_setting_on_one_lineage_fold_reproduces_the_drawn_fold` |
| 11. `--exp-settings` without `--force` | `learn_takes_the_settings_flags_only_with_force_and_alone` |

All 14 plants were caught on 23 Nov 2026, and each test passes with no fault planted. A first run caught 13: fault 7b's test also failed because its added row had no choices files, so it now passes only when the twin-seed check itself fails.

## What it does not do

It reads no sealed, held-out, clinical or spent test set: neither the outside-study genomes, nor the other-database genomes, nor the later-records set. It finds, names and ranks no resistance marker; a choices file holds settings and counts only. It changes no default, report, release rule or the site.
