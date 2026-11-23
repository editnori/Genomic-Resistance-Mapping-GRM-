# Input arms: other inputs beside 31-mer presence

An experiment, off by default, for plan block C3 (`docs/assumption-plan.md`). The release fits its five learners on the presence of 31-mers in each genome (K31). C3 asks whether other inputs separate the classes better: unitigs, K31 plus k-mer counts, K31 plus amino-acid k-mers of called genes, the AMRFinderPlus program's calls alone, and K31 plus the program's calls. This branch builds the two program-call arms and designs the other three; they wait for the lead's ruling on their tools.

`grm learn --exp-inputs ARM` fits one arm in place of K31. Without the flag nothing new is called and `grm learn` is unchanged, byte for byte. `EVALUATION_REVISION` stays 14 and the report schema stays 2. An arm writes no report and removes none.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation. GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalog; "the program" on this page is the AMRFinderPlus program itself, a separate thing.

## The arms

| Arm | Columns | State |
|---|---|---|
| K31, the default | the release's, `grm learn` with no `--exp-inputs` | runs today |
| `program` | one per element symbol of the program's Type AMR rows | built |
| `k31-program` | the program's columns, then the dataset's K31 rows less those whose cohort pattern a program column has | built |
| `k31-counts` | K31 rows, then [count ≥ 2] of each 31-mer in a genome | designed, waits |
| `unitigs` | one per distinct unitig pattern, in place of K31 | designed, waits |
| `k31-aa` | K31 rows, then amino-acid 10-mers of called genes | designed, waits |

The program is AMRFinderPlus 4.2.7 with database 2026-08-07.1 and `--plus`, run from the AMRFinderPlus kit as `scripts/refit/run_species.sh` runs it. Each table must equal the sha256 the release recorded for it before any fit reads one. A genome has a column when its table has a row of Type AMR with that element symbol, either scope and any method. The program arms run on A, 237 datasets in 21 species: database 2026-08-07.1 has no organism option for *Mycobacterium tuberculosis* or *Staphylococcus haemolyticus*.

## Inside the folds

An arm hands `evaluate` a new matrix and nothing else. `analysis::project_cohort`, which both the drawn splits and `--fold-file` use, passes the dataset's K31 projection and its genome ids to `input_arms::apply`. The release's per-fit code then runs over the arm's matrix as it runs over K31: the support floor of max(3, n/200) and the collapse of identical patterns on each fit's training genomes, and the chi-square screen to 4,000 for forest, boosted trees and lasso on each fit's own fitting genomes.

What an arm builds once per dataset reads no label. Its one cohort-wide filter is the one `Dataset::project` applies to K31: a program column is kept when its support over the cohort lies in [floor, n − floor], and identical cohort patterns collapse to the first in byte order of the symbol. Every training subset implies that range, so the filter can neither add nor drop a training candidate, whatever the outer test genomes hold. A cohort genome with no table stops the run on that dataset; no row is ever filled.

The folds come from the labels, the lineages and the seed, never from a column, so every arm at seed s is scored on the default's folds at seed s. The scorer checks that.

## The files

`grm learn --exp-inputs ARM --exp-program-tables DIR --exp-pick OUT --exp-pick-inner 0 --force` reads `DIR/<genome id>.tsv` and writes, per dataset, split and seed, C9's result file at no inner fold (`docs/served-pick.md`, The result file): each outer fold's genomes by class and each learner's calls, with no genome id. Beside it, per dataset, it writes `OUT/<dataset>.<tag>.inputs.json`, format `grm-c3-inputs-1`, where the tag is `drawn` for the drawn splits and the fold file's scheme for a study run. That file holds counts only: the tables read, the program names seen, out of range, collapsed and kept, the K31 rows, those left out as a program column's pattern, and the arm's rows. It holds no column name, genome id, k-mer, sequence or position.

`--exp-program-column` takes `symbol` (the default), `subclass` or `class`, and `--exp-program-rows` takes `amr` (the default) or `all`, so the lead's other answers to questions 7 and 8 need no new code. The flag needs `--force`, so a kept report can never skip an arm, and it refuses beside `--learner-seed`, the `--exp-tie-*` options, `--exp-null-refit`, `--exp-mic` and `--exp-cohort`.

## The scorer

`scripts/c3-inputs/score.py` checks every result `expect.py` lists and stops before scoring on any failure:

- C9's checks at no inner fold: format, revision, dataset, scheme, seed, the five learners in order, each fold's calls equal to its test genomes, and the calls summed to the split's balanced accuracy;
- each arm's folds equal to the default's at the same dataset, split and seed, by class;
- the default's seed-42 lineage figures and served learner equal to the release report's;
- the twin at seeds of its own, and neither the default nor the twin with an inputs file;
- each inputs file of its arm, dataset and tag, with a table for every cohort genome and its rows adding up.

It scores each learner on its outer calls summed over the folds: balanced accuracy, recall, specificity, the very major and the major error rate, each with its count, denominator and Wilson 95% interval, and log loss. A dataset's value is the mean of its draws.

`score.py --count` takes several species' tables once every arm of the family has results, and refuses before that. Per arm and learner it counts the datasets better, worse and the same by 0.02 against K31 on the arm's own datasets. Every difference is the arm's value less the default's, at 12 decimals; lower is better for log loss and both error rates. Beside the count stand the median difference with its interval from 4,000 species resamples and the one-sided sign test. Holm runs over 30 lineage pairs: 6 arm rows (unitigs, counts, amino-acid, program, K31 plus program and the twin) by 5 learners. An input arm touches all five learners, so conditions 1 and 2 are read as met for at least 3 of 5.

## Datasets

`scripts/c3-inputs/lists.py` writes `scripts/c3-inputs/lists.tsv` from `scripts/c9-pick/l.tsv` and the AMRFinderPlus comparison's record, with each species and antibiotic in full.

| Set | Datasets | Species | Less M, the first run | First run in S |
|---|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 |
| A | 237 | 21 | 204 | 54 |

On the datasets not in M, three of the four label rules give the same labels. Of L's 233, 32 hold an S and R pair and depend on the most resistant call against the other three; of A's 204, 21 do. The rest depend on none of the four. Every measure on M waits for the maintainer's label rule.

## Running

`scripts/c3-inputs/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARMS]` runs one species on a box of its own, with no prompt. It checks its kit by sha256 (the release, the data root's list, the curated table, the study plans, the label count's table, the AMRFinderPlus kit and the recorded tables' sha256 list). It refuses a dataset in M, and refuses at the kit step a box without numpy, `cargo`, `objcopy`, `unshare`, `flock`, `/usr/bin/time`, `curl`, `tar`, `xz` or `zstd`. It builds `grm` and the species' data root as `scripts/c9-pick/run_species.sh` does, installs the program and proves it on the kit's fixture, runs it on every genome of the species with a recorded table, and holds each table to its recorded sha256.

It fits the default at each seed and checks its seed-42 reports against the release's with `scripts/refit/refit_check.py` before any arm runs. It then fits each arm at each seed, the twin, and the study plans, scores, writes `logs/timing.tsv` of each step's wall time and peak memory, and writes `OUT/DONE` with the sha256 of its list of output files. Each log is copied with every genome id written as `[genome id]`, and the end fails if a file under `results/`, `check/`, `score/` or `logs/` names one. Run again, it skips finished steps and finished results. Any failure writes `OUT/FAILED` and exits 1.

## Tests and planted faults

`input_arms` holds 8 SYNTHETIC tests and `src/main.rs` one more on the flags. `tests/c3_inputs` holds 33 SYNTHETIC tests of the scorer, the lists and the expected results. `scripts/c3-inputs/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must fail.

| Fault | Caught by |
|---|---|
| 1. A cohort filter that reads the outer test genomes' calls | `no_outer_test_genome_changes_a_fits_candidates`: emptying or filling the outer test genomes' calls changes no fit's candidates |
| 2. A joined arm that drops K31 (2a), holds the program block twice (2b) or keeps a K31 copy of a program column (2d), or a scorer that passes one (2c) | `the_joined_arm_holds_each_block_once_the_program_first`; `a_k31_row_equal_to_a_program_column_enters_once_as_the_programs`; `Checks.test_a_joined_arm_that_drops_or_doubles_a_block_fails` |
| 3. A program column built from the label | `an_arm_matrix_does_not_change_with_the_labels` |
| 4. An arm scored on other folds (4a) or paired with the default at another seed (4b) | `Checks.test_an_arm_on_other_folds_fails`; `Checks.test_an_arm_is_held_to_the_default_at_its_own_seed` |
| 5. A twin at a draw's seed, in the scorer (5a) or the list (5b), or fitted with an arm (5c) | `Checks.test_a_twin_at_a_draw_seed_fails`; `Expect.test_a_twin_at_a_draw_seed_is_refused`; `Checks.test_a_twin_with_an_inputs_file_fails` |
| 6. A default not held to the release | `Checks.test_a_default_that_differs_from_the_release_fails` |
| 7. A difference taken as the default's value less the arm's | `Count.test_every_difference_is_the_arms_value_less_the_defaults` |
| 8. Holm over one arm's five learners | `Count.test_holm_runs_over_the_thirty_lineage_pairs_of_the_primary_measure` |
| 9. A count taken on part of the family | `Count.test_the_count_refuses_part_of_the_family` |
| 10. An arm run without `--force` | `learn_takes_the_c3_input_arm_only_when_asked` |

## What it does not do

It reads no sealed, held-out, clinical or spent test set: neither the outside-study genomes nor the other-database genomes. It finds, names and ranks no resistance marker; a program column is named only as the program's table names it, and no column name leaves a fit. It changes no report and no release rule. A winning program-call arm would enter only with the program's version pinned, and the maintainer would rule how GRM may name its calls.
