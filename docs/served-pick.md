# Served learner: an inner pick against the release's pick

An experiment, off by default, for plan block C9 (`docs/assumption-plan.md`). The release serves, per dataset, the learner with the best random-split balanced accuracy among those that call every isolate, a tie going to the later learner (`report::best_learner`). C9 asks whether that is the right learner to serve. Its other arm picks a learner inside each outer lineage training fold, from inner folds grouped by lineage, and scores the pick on that outer fold's test genomes.

`grm learn --exp-pick DIR` writes one result per dataset, split and seed to `DIR`. The reports it writes equal the release's in every fitted field, as `scripts/refit/refit_check.py` compares them; without the flag `grm learn` is unchanged, byte for byte. `EVALUATION_REVISION` stays 14 and the report schema stays 2. The release keeps its rule, its pick and its cap of 2,000 genomes a species.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation.

## The two arms

| Arm | Serves on outer fold f | From |
|---|---|---|
| Default | The release's pick, the same learner on every fold | The same run's random split, by `report::best_learner` at that run's seed |
| Inner pick | The learner with the best inner balanced accuracy in fold f | Inner folds of fold f's training genomes only |

Both arms are scored with the release's own outer fits: the calls each learner made on fold f's test genomes in the split the result is of, lineage or study, at the cutoff of 0.5 (`predict::CALL_CUTOFF`). A null fit is scored at its training prevalence, as the report scores it. Summed over the outer folds, the default's calls give the report's lineage row for its learner.

The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism; the inner pick reads only each outer fold's training genomes, so its served figure cannot. Wherever the two are compared, a same or a worse for the inner pick is never read as showing that the release's pick is right.

## The inner folds

Outer fold f of the lineage split trains on its genomes T_f. The inner pick:

- **Partition.** `training_folds(y, T_f, lineage groups, 3, fold_seed(seed, f))`: 3 folds of T_f, each lineage whole, drawn at the fold's own seed (`ablate::plan::fold_seed`). At 3 folds this is the partition the release's SCM and CART tune on in fold f. A genome outside T_f gets no inner fold.
- **Fits.** For inner fold j, each of the five learners (`scm`, `cart`, `rf`, `gbt`, `lasso`) is fitted on T_f less fold j as the release fits an outer fold: candidates from its own training genomes at the dataset's support floor, tuning on the split's tuning groups, and the seed FNV-1a 64 of `c9-inner/<seed>/<f>/<j>` for its own inner folds, tuning holdouts and draws. It scores fold j's genomes; a null inner fit scores them at its inner training prevalence.
- **Measure.** A learner's inner figure is the balanced accuracy of its calls summed over the 3 inner folds. The figures are compared at 12 decimals.
- **No pick.** Each inner fold must hold both classes in its training and test genomes, and each inner fit's own tuning folds and holdouts must hold both classes; each inner training part must give a candidate. If a check fails, fold f has no pick and records why. A learner whose inner fit ends in an error is left out of fold f's pick; with none left, the fold has no pick. A fold with no pick serves the default's learner and is counted as such (question 4, ruled as taken).
- **Ties.** The learners tied at the best figure are listed. The pick is the later in pipeline order, as `best_learner` breaks a tie, or the earlier with `--exp-pick-tie earlier`.

`--exp-pick-inner 0` fits no inner fold: the run writes the default's pick and the outer calls only, which is what the twin needs. `--exp-pick` is refused beside `--learner-seed`, the `--exp-tie-*` options, `--exp-null-refit`, `--exp-mic` and `--exp-cohort`.

## The result file

`DIR/<dataset>.<scheme>-s<seed>.pick.json`, format `grm-c9-pick-1`:

- the revision, dataset, scheme, seed, inner folds, tie rule, support floor, samples and tested genomes by class;
- the split's reason it is unavailable, if it is;
- on the lineage split, the default's learner, its random-split balanced accuracy and every learner tied with it there; a study result has no default of its own;
- each learner's pooled balanced accuracy on the split, as its report row has it;
- per outer fold: its training and test genomes by class, each learner's calls on its test genomes (`n`, `tp`, `fp`, `tn`, `fn`, log loss sum), and the inner pick: its status, its reason, each inner fold's genomes by class, each learner's calls per inner fold and summed, its inner figure, its null fits and its error, the tied learners and the pick;
- the seconds the inner pick took.

It holds counts and calls only: no genome id and no row per genome.

## Draws, twin and study

- **Draws.** The default and the inner pick come from the same runs: `grm learn --seed s --exp-pick DIR` at s = 42, 43 and 44. A dataset's value is the mean of its 3 draws.
- **Twin.** The default at seeds 45, 46 and 47, from runs with `--exp-pick-inner 0`. Its difference is the mean of its 3 draws less the default's. It is the count with no change.
- **Study.** On the datasets in S: `grm learn --fold-file <dataset>.study-whole.json --exp-pick DIR`, the reviewed study plan. The inner folds are grouped by lineage inside each study training fold. The default is the release's pick from the seed-42 lineage run. `grm learn` refuses `--seed` beside `--fold-file`, so the study split has one run per dataset, at the release's seeds, and no twin. That departs from the plan's Rules, Draws, which give each arm seeds 42, 43 and 44; the entry "C9 design" of `docs/assumption-plan-corrections.md` names it as a deviation (question 6, ruled as taken).

## The scorer

`scripts/c9-pick/score.py EXPECT PICKS RELEASE OUT` reads the results that `expect.py` lists and checks each one:

- it exists, with its format, revision, dataset, scheme, seed and inner folds, and the five learners in order;
- each learner's outer calls, summed over the folds, give the balanced accuracy of its split row, and each fold's calls are its test genomes';
- each inner pick's folds cover its outer fold's training genomes exactly, its summed inner calls give its inner figures, and its pick and tied set follow from those figures under the file's tie rule;
- a twin fits no inner fold, and a fold with no pick says why;
- the default's learner is the later of its tied set, as the release breaks a tie;
- at seed 42, each learner's lineage figure and the default's random-split figure equal the release report's, and the default's learner is the one the release serves (`best_learner` in the release's `index.json`).

A failed check stops the scorer before it scores anything. It writes `checks.tsv`, `draws.tsv`, `datasets.tsv` (species and antibiotic in full), `picks.tsv` and `options.txt`, and prints its check counts only.

Per dataset, split, seed and arm, the served calls are summed over the outer folds, each fold taking the calls of the learner it serves. From them: balanced accuracy (the primary measure), recall, specificity, the very major error rate (resistant called susceptible) and the major error rate (susceptible called resistant), each with its count, denominator and Wilson 95% interval, and log loss with scores clipped at 1e-4. Pick counts beside them: folds per served learner, ties, folds with no pick, and folds whose pick is the default's learner.

Three options carry the design's questions 3 to 5: `--tie later|earlier|default` (`default` sends a tie to the default's learner when it is tied), `--no-pick default|drop` (serve the default's learner, or leave the draw unscored), and `--default per-draw|seed42` (each draw's own random-split pick, or the seed-42 pick at every draw). The lead ruled them at 20:36Z on 22 Nov 2026, before any C9 score was read. `--tie later` is read; `earlier` and `default` are scored beside it and decide nothing. `--no-pick default` is the answer taken. `--default per-draw` is read; the seed-42 pick is written beside it and decides nothing.

`score.py --count OUT <datasets.tsv>...` counts over species, as the plan's Rules do. Per split, pair (the inner pick and the twin, each against the default) and measure, a dataset is better, worse or the same by more than 0.02, compared at 12 decimals; lower is better for log loss and both error rates. Beside the count: the median difference, which for every measure is the arm's value less the default's (log loss and both error rates included), its 95% interval from 4,000 species resamples (the 101st and 3,900th sorted values, each row's generator seeded from the sha256 of its split, pair and measure, every median rounded to 12 decimals), and the one-sided sign test of better against worse. Holm runs over the lineage split's primary measure for the inner pick and the twin. `conditions.tsv` reads the plan's change conditions 1 to 5 off the count; condition 6 is the maintainer's. The count compares the two arms, so the reading under The two arms holds here: the default's served figure can carry its pick's optimism and the inner pick's cannot, and a same or a worse is never read as showing that the release's pick is right.

## Datasets

`scripts/c9-pick/lists.py` writes the lists from the release reports, the label count's table and the reviewed study plans. `scripts/c9-pick/l.tsv` is its list of L as it wrote it: one row per dataset with its species and antibiotic in full, its genomes by class, its lineages, whether it is in M, in S and in the first run, its genome pairs of S and R calls, and the label rules its measures depend on.

| Set | Datasets | Species |
|---|---:|---:|
| L: reports whose lineage split ran | 268 of 274 | 23 |
| L not in M: the first run | 233 | 23 |
| L in M: waits for the label rule | 35 | 9 |
| S in L | 69 of 70 | |
| First run in S: the first study run | 55 of 233 | |
| M in S | 14 of 35 | |

M marks a dataset with a genome pair of S and I calls or of I and R calls. On the 233 datasets not in M, the label rules `r-or-s-wins` (the release's), `intermediate-aside` and `intermediate-class` give the same labels, so no measure there depends on which of the three the maintainer takes. 32 of the 233 hold a genome pair of S and R calls, which the most resistant call labels R; their measures depend on that rule against the other three, and the other 201 depend on none of the four. Every measure on the 35 in M depends on the maintainer's rule.

## Running

`scripts/c9-pick/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN]` runs one species on a box of its own. It checks its kit by sha256 (the release, the data root's list, the curated table, the study plans and the label count's table), refuses a dataset in M, and refuses at the kit step a box without numpy, `objcopy`, `unshare`, `flock` or `/usr/bin/time`. It builds `grm` and records its commit, its sha256 and the sha256 of its code alone, then builds the species' data root from public files: the catalog equal to a pinned list, and each genome with the release's length, contig count and N50. It fits every seed, checks the seed-42 reports against the release's with `scripts/refit/refit_check.py` before any other seed runs, fits the study plans, scores, and writes `OUT/DONE` with the sha256 of its list of output files.

Each step's log is copied with every genome id of the species written as `[genome id]`, and the end fails if any file under `results/`, `check/`, `score/` or `logs/` names one. Run again, the script skips finished steps and finished results. On a later checkout it goes on only when the binary is unchanged and the commits differ only under `scripts/c9-pick/`, `docs/`, `tests/` and `.github/`. Any failure writes `OUT/FAILED` and exits 1.

## Tests and planted faults

`learn::inner_pick` holds 10 SYNTHETIC tests on a fixture of 150 genomes in 30 lineages, and `src/main.rs` one more on the flags. `tests/c9_pick` holds 40 SYNTHETIC tests of the scorer, the lists and the expected results. `scripts/c9-pick/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must fail.

| Fault | Caught by |
|---|---|
| 1. An outer test genome in an inner fold | `no_outer_test_genome_reaches_an_inner_fit`: genomes outside the outer training fold, every feature and label flipped, change no inner pick |
| 2. An inner split that ignores lineage groups | `inner_folds_keep_each_lineage_whole_and_are_the_release_tuning_partition` |
| 3. A pick that ignores its seed: the inner folds (3a) or the inner fits (3b) | `the_seed_keys_every_inner_fit_and_partition` |
| 4. A tie sent the wrong way: in `grm` (4a) or in the scorer's re-pick (4b) | `a_tie_goes_where_the_rule_says_and_is_listed`; `Served.test_a_tie_is_picked_again_by_the_scorers_rule` |
| 5. The default pick read from the lineage split | `the_default_pick_reads_the_random_split_as_the_release_does` |
| 6. The picked learner scored on its inner folds in place of the outer test fold | `Served.test_the_served_figure_is_each_folds_pick_on_its_outer_test_fold` |

## What it does not do

It reads no sealed, held-out, clinical or spent test set: neither the outside-study genomes nor the other-database genomes. It finds, names and ranks no resistance marker. It changes neither the release's pick nor any report; a change to the served learner would need the plan's six conditions and the maintainer's ruling, since the site names the served learner.
