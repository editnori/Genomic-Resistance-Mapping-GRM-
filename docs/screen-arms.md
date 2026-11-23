# Screen arms: the screen and the floor

An experiment, off by default, for plan block C4 (`docs/assumption-plan.md`). The release keeps, for each fit, the 31-mer rows whose support in the fit's training genomes lies in [floor, n − floor] with floor = max(3, n/200), collapses identical patterns, and lets forest (`rf`), boosted trees (`gbt`) and the lasso keep the top 4,000 of those by the 2x2 Pearson chi-square on their own fitting genomes; SCM and CART read every candidate. C4 asks whether that screen or that floor removes the inputs that separate the classes, or keeps lineage markers. D0 put 15 of 177 miss pairs in Removed, all on `rf`, `gbt` and `lasso` (`docs/refit-result.md`).

This branch builds three arms behind two flags. Without them nothing new is called and `grm learn` is unchanged, byte for byte (The unflagged bytes, below). `EVALUATION_REVISION` stays 14 and the report schema stays 2. An arm writes no report and removes none. No comparison has been read.

Any split here is a retrospective split inside one public database. A split by lineage or by study is not external validation. GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalog; nothing in this block reads the AMRFinderPlus program.

## The arms

| Arm | `grm learn` flags beside `--exp-pick DIR --exp-pick-inner 0` | Screen | Budget | Floor |
|---|---|---|---|---|
| `k31`, the default | none | chi-square | 4,000 | max(3, n/200) |
| `cmh-4000` | `--exp-screen cmh --force` | Cochran-Mantel-Haenszel over lineage strata | 4,000 | max(3, n/200) |
| `cmh-16000` | `--exp-screen cmh --screened-features 16000 --force` | Cochran-Mantel-Haenszel over lineage strata | 16,000 | max(3, n/200) |
| `floor-3` | `--exp-floor 3 --force` | chi-square | 4,000 | 3 |
| `twin` | none, at seeds 45, 46 and 47 | chi-square | 4,000 | max(3, n/200) |

The two screen arms take the default's candidates and change only how `rf`, `gbt` and `lasso` rank them. `scm` and `cart` read no screen, so under a screen arm they are refitted and must make the default's calls in every fold. The floor arm changes the floor alone, in each fit and in the dataset's projection (`Dataset::project`), since a fit at floor 3 can need any row a projection at 3 holds. Where n < 800, max(3, n/200) is 3 already and the floor arm is the default. The twin is the default at seeds of its own.

`--exp-screen` and `--exp-floor` each need `--force`, so a kept report can never skip an arm, and `--exp-pick` with `--exp-pick-inner 0`. Each is refused beside the other, `--learner-seed`, `--exp-tie-rule`, `--exp-tie-trace`, `--exp-tie-out`, `--exp-null-refit`, `--exp-mic`, `--exp-cohort` and `--exp-inputs`. `--exp-floor` is refused beside a budget other than 4,000, so the floor arm cannot change the budget.

## The statistic

For one candidate and one fit, group the fit's fitting genomes by lineage cluster. Stratum k holds n_k genomes, r_k of them resistant; the candidate is present in m_k of them, a_k of those resistant. A stratum is informative when n_k ≥ 2 and 0 < r_k < n_k.

    D = Σ_k (a_k − m_k r_k / n_k)
    V = Σ_k m_k (n_k − m_k) r_k (n_k − r_k) / (n_k² (n_k − 1))
    statistic = D² / V, or 0 when V = 0

The sums run over informative strata in ascending order of the cluster id, in f64 as the default's chi-square is. There is no continuity correction, and the statistic is two-sided. A stratum of one genome or one class adds 0 to D and V, so a candidate constant inside every informative stratum, a pure lineage marker, scores 0. With one stratum holding every fitting genome the statistic is (n − 1)/n times the default's chi-square, so the two screens rank alike and differ by the strata alone. The top k keep the highest statistics, ties to the lower candidate row. When fewer than k candidates score above 0, the budget fills with candidates at 0 by lower row, as the default fills with its own zeros; the inputs file counts them.

## Inside the folds

The strata are the release's lineage clusters (`structure.json`, Mash distance at the release's threshold, the report's `samples[].lineage`), which read no label. Each screen restricts them to the genomes it fits: `rf` the outer fold's training genomes, `gbt` those less its 15% early-stopping holdout, `lasso` those less its 20% penalty holdout, as the default screens each. The final models screen once each, `rf` on every genome of the dataset and `gbt` and `lasso` on those their holdouts leave. The same rule holds on the random, lineage and study splits. On the lineage split a test lineage is wholly out of the fit. On the study split a cluster can hold training and test genomes, and only its fitting genomes enter. The labels a screen reads are those of its fitting genomes, so no outer test genome can change a fit's candidates or its screen. The folds come from the labels, the lineages and the seed, never from a column, the screen or the floor, so every arm at seed s is scored on the default's folds at seed s.

## The files

With either flag, `grm learn` writes C9's result file at no inner fold (`docs/served-pick.md`, The result file) per dataset, split and seed, and beside it `DIR/<dataset>.<tag>-s<seed>.inputs.json`, format `grm-c4-screen-1`, where the tag is `drawn` or the fold file's scheme. It holds counts only:

- the arm, the rule (`cmh` or `chi-square`), the budget, the floor, the dataset's genomes and its rows (the candidates in);
- per split and outer fold, the training genomes and the candidates past the floor;
- per screen (`rf`, `gbt` and `lasso` in each outer fold, and the final models), its rule, the fitting genomes, the candidates, the columns kept and those kept at 0, and for the adjusted screen the strata and informative strata among the fitting genomes;
- per outer fold and learner, the fit's wall and CPU seconds and the process's peak memory.

It holds no genome id, feature, k-mer, sequence or position. The default and the twin write none.

## The scorer

`scripts/c4-screen/score.py EXPECT PICKS RELEASE OUT` checks every result `expect.py` lists and stops before scoring on any failure. It imports C3's scorer (`scripts/c3-inputs/score.py`) and adds what this block needs:

- the default and the twin by C3's checks: C9's checks at no inner fold, the default's seed-42 lineage figures and served learner equal to the release report's, the twin at seeds of its own, and neither with an inputs file;
- every arm on the default's folds at the same dataset, split and seed, and each result's scheme its split's;
- a screen arm's `scm` and `cart` calls equal to the default's in every fold;
- below 800 genomes, the floor arm's result equal to the default's, timing apart;
- each inputs file of its arm, dataset, tag and seed: its rule and budget, its floor (3, or the default's for n), the result's outer folds and training genomes, every screen by the arm's rule, each keeping the budget or every candidate if fewer, `rf`, `gbt` and `lasso` once in each outer fold, `rf` on its training genomes and `gbt` and `lasso` on fewer, the final models once each (none in a fold file's run), the adjusted screen's strata, and a fit time for each learner of each outer fold.

It scores as C3 does: per learner, the outer calls summed over the folds give balanced accuracy, recall, specificity, the very major and the major error rate, each with its count, denominator and Wilson 95% interval, and log loss. A dataset's value is the mean of its draws, and a dataset with a draw missing is not scored.

`score.py --count OUT [--withdrawn ARM] DATASETS.tsv ...` takes every species' table once every arm has results or is withdrawn. Per arm and learner it counts the datasets better, worse and the same by 0.02 against the default, the arm's value less the default's at 12 decimals, with the median difference, its interval from 4,000 species resamples and the one-sided sign test. A screen arm is counted on the learners it touches, `rf`, `gbt` and `lasso`, over every dataset; the floor arm on all five over the datasets of 800 genomes or more, the ones its floor can move; the twin on all five over every dataset and again over those of 800 genomes or more. Holm runs over those 21 lineage pairs (the lead's rulings 10 and 12 of 23 Nov 2026). `count.tsv`'s column `set` names each row's datasets, `all` or `800-or-more`. `conditions.tsv` reads the plan's conditions per arm:

- a screen arm: 1, 2 and 5 for each of `rf`, `gbt` and `lasso`, and 4 for the three;
- the floor arm, over the datasets of 800 genomes or more: C3's ruling 15, a learner carrying the change when it meets 1 and 2, at least 3 of 5 carrying, 4 for all five and 5 for each that carries;
- condition 3 for every arm as corrected for C3 on 23 Nov 2026, with the twin read on the arm's datasets: the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart.

`count-no-constant.tsv` leaves out the constant-model results of C3's ruling 16, which the screen arms cannot have. `count-floor-all.tsv` counts the floor arm on every dataset, those under 800 genomes entering as a difference of 0, with no Holm, and decides nothing.

## Datasets

`scripts/c4-screen/lists.py` writes `scripts/c4-screen/lists.tsv` from C3's L (`scripts/c3-inputs/lists.tsv`), the first run's datasets first, with each species and antibiotic in full, the default floor and whether the floor arm moves the dataset.

| Set | Datasets | Species | Less M, the first run | First run in S | Floor arm moves |
|---|---:|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 | 92 in 11 species, 24 in S |

The floor arm moves no dataset under 800 genomes. On the other 141 of the first run its difference is 0 by construction, so over all 233 its median would be 0 whatever the floor did. The lead ruled at 05:52Z on 23 Nov 2026, before any score, that the floor arm is read over the 92 it can move, with the count over all 233 beside it.

## Running

`scripts/c4-screen/run_species.sh KIT SPECIES DATASETS OUT [JOBS] [SEEDS] [TWIN] [ARMS]` runs one species on a box of its own, with no prompt. Its steps up to the fits are C3's, read from `scripts/c3-inputs/run_species.sh`: the kit held to its sha256, `grm` built from this checkout, the data root, the count, the lineage clusters, the matrix and the reviewed study plans. It fits the default at each seed and checks its seed-42 reports against the release's before any arm runs, then each arm at each seed, the twin, and one study run per arm and dataset in S at the release's seeds. It scores, writes `logs/timing.tsv` and `OUT/DONE` with the sha256 of its list of files, and fails when a file under `results/`, `check/`, `score/` or `logs/` names a genome id. Run again, it skips finished steps and results; it refuses to go on if the binary changed or a file outside `scripts/c4-screen/`, `docs/`, `tests/` or `.github/` did.

The pilot, one draw of *Enterobacter cloacae* (9 datasets of 80 to 124 genomes, 8 in S) at seed 42 with no twin, ran on a 16-CPU, 32 GB box on 23 Nov 2026: 68 of 68 results, 2,231 of 2,231 checks, 0 genome ids in what it keeps. The fit steps took 54.27 s for the default, which also writes the reports, 31.67 s for `cmh-4000`, 35.59 s for `cmh-16000` and 26.36 s for `floor-3`. Summed over 650 fits, the screened learners took 5.66 s of wall under `floor-3` (the default's computation here), 6.56 s under `cmh-4000` and 16.55 s under `cmh-16000`, the lasso most of the last (11.92 s against 2.87 s). Each fit process peaked near 1.3 GB, the count at 2.0 GB. No pilot dataset reaches 800 genomes, so the pilot shows the floor arm equal to the default and cannot time a moved floor.

## The unflagged bytes

`scripts/c4-screen/unflagged.sh` runs `grm` built at main and at this branch, with no C4 flag and no network, each on its own copy of a data root: `grm learn` for the reports, the same with `--exp-pick` for the default's results, and `--fold-file` on each study plan. It does so on the pilot's data root, and on a root of SYNTHETIC genomes made from the contigs of the AMRFinderPlus program's bundled `test_dna.fa` and drawn filler, which `grm count`, `structure` and `matrix` build first. Every file of the two copies must be equal, the run clocks of JSON files apart. On 23 Nov 2026 both held (results in the entry "C4 design" of `docs/assumption-plan-corrections.md`). The Rust test `the_c4_arms_leave_unflagged_fits_alone_and_screen_every_fit_by_their_rule` holds the same on SYNTHETIC fits in CI.

## Negative control

A SYNTHETIC cohort of 400 genomes in 20 strata of 20, the resistant share of the strata from 0.1 to 0.9, with 2,000 candidate rows in a shuffled order: 1,000 lineage markers, present at 0.8 in a random half of the strata and 0.2 in the other half, and 1,000 rows present at 0.5 everywhere. Each of 200 seeded permutations permutes the labels within each stratum, which keeps every stratum's class counts and breaks any link inside a stratum, and each screen keeps its top 100 on all 400 genomes. A permutation rejects when the top 100 hold more lineage markers than a one-sided hypergeometric test at 0.05 allows against 100 draws from the 2,000. The adjusted screen passes at 15 or fewer of 200, about the 95th percentile of Binomial(200, 0.05). The test and its cut were fixed before it ran.

The adjusted screen rejects in 7 of 200; the chi-square, beside it and deciding nothing, in 200 of 200, so the test can see lineage signal. The test is `negative_control_permuted_labels_within_strata_carry_no_lineage_signal`.

## Tests and planted faults

`learn::screen_arm` holds 9 SYNTHETIC tests, among them the statistic against `scripts/c4-screen/cmh_plain.py`, a second, plain implementation in exact rationals, on the committed SYNTHETIC tables of `tests/c4_screen/cmh_cases.json`; `learn::eval` 2 more and `src/main.rs` one on the flags. `tests/c4_screen` tests the scorer, the lists, the expected results, the plain implementation and the unflagged check on SYNTHETIC files. `scripts/c4-screen/planted_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must pass with no fault and fail with it.

| Fault | Caught by |
|---|---|
| 1. A screen that reads an outer test genome | `no_sample_outside_the_fit_changes_the_screen` |
| 2. Strata taken from labels | `the_strata_do_not_change_with_the_labels` |
| 3. A stratum of one class counted | `a_stratum_of_one_class_or_one_sample_changes_no_statistic` |
| 4. The floor arm at another budget, in `grm` (4a) or the scorer (4b) | `learn_takes_the_c4_screen_and_floor_arms_only_when_asked`; `Checks.test_a_floor_arm_at_another_budget_fails` |
| 5. An arm fitted with the default's screen (5a), or passed by the scorer (5b) | `a_fit_under_the_adjusted_arm_takes_its_screen_and_with_no_arm_the_chi_square`; `Checks.test_an_arm_screened_by_the_defaults_rule_fails` |
| 6. The twin scored as an arm (6a) or listed with an inputs file (6b) | `Checks.test_the_twin_passes_as_the_default_at_its_own_seeds_and_fails_with_an_inputs_file`; `Expect.test_only_an_arm_names_an_inputs_file` |
| 7. A missing draw read as a zero | `Scores.test_a_dataset_with_a_missing_draw_is_not_scored` |
| 8. A study fold read as a lineage fold, in the list (8a) or the scorer (8b) | `Expect.test_each_result_is_on_its_own_split`; `Checks.test_a_study_result_read_as_a_lineage_result_fails` |
| 9. Holm over one arm's learners | `Count.test_holm_over_twenty_one_is_held_to_the_plain_step_down` |
| 10. A screen arm's `scm` and `cart` counted | `Count.test_a_screen_arms_scm_and_cart_are_not_counted` |
| 11. A floor arm unlike the default below 800 genomes passed | `Checks.test_a_floor_arm_below_800_genomes_unlike_the_default_fails_and_at_800_may_differ` |
| 12. A screen arm whose `scm` or `cart` differs from the default's passed | `Checks.test_a_screen_arm_whose_scm_or_cart_differs_from_the_default_fails` |
| 13. The floor arm run without `--force` | `learn_takes_the_c4_screen_and_floor_arms_only_when_asked` |

On 23 Nov 2026 the runner caught 17 of 17.

## What it does not do

It reads no sealed, held-out, clinical or spent test set: neither the outside-study genomes nor the other-database genomes. It finds, names and ranks no resistance marker, and no feature, k-mer or genome id leaves a fit. It changes no report, no default and no release rule; the release keeps its cap of 2,000 genomes a species. It cuts no new lineage clusters: the strata are the release's.
