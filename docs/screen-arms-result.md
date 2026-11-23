# Screen and floor result

The release ranks the candidates of the forest (`rf`), boosted trees (`gbt`) and the lasso (`lasso`) by the chi-square statistic and keeps the top 4,000, after a candidate floor of max(3, n/200) for a dataset of n genomes. This page asks whether the Cochran-Mantel-Haenszel screen over lineage strata, at budgets of 4,000 and 16,000 (`cmh-4000`, `cmh-16000`), or the floor fixed at 3 (`floor-3`), serves better on the lineage split than the release's rule. It reads block C4 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/screen-arms.md`.

The reading's rules were fixed at 04:03Z on 23 Nov 2026 and ruled at 05:52Z, before any C4 score was read (`docs/assumption-plan-corrections.md`, entry "C4 design").

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

The Cochran-Mantel-Haenszel screen adjusts the ranking of the candidates for lineage strata and names no marker. It changes which candidates a learner reads, and this page counts only how the learners then score.

"No difference shown" is the rule's reading: an arm that does not meet the change conditions leaves the release as it is. It is no statement that two screens are equal.

## Reading

**Conditions.** A screen arm is read on `rf`, `gbt` and `lasso`, the learners its screen reaches, over every dataset of the first run. It meets condition 1 when, for each of the three, the lineage median difference's interval lies wholly above 0, and condition 2 when, for each of the three, better beats worse by the one-sided sign test at p < 0.05 after Holm over the 21 lineage pairs. Condition 3 holds when the twin, over every dataset, meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. Condition 4 holds when no recall interval of the three lies wholly below -0.02, and condition 5 when each of the three has a study median of its lineage median's sign.

The floor arm is read on all five learners over the datasets of 800 genomes or more, by C3's ruling 15: a learner carries the change when it meets 1 and 2, and at least 3 of 5 must carry it. Condition 3 reads the twin over the same datasets, condition 4 holds for all five, and condition 5 for each learner that carries. The stricter reading of C3 stands beside it and decides nothing. Condition 6, for every arm, is the maintainer's.

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

## Count

The primary measure is each learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. The median difference is the arm's value less the default's, with a 95% interval from 4,000 species resamples, as C3's count takes it.

The Holm family holds 21 pairs on the lineage primary measure: `cmh-4000` and `cmh-16000` each on `rf`, `gbt` and `lasso` over every dataset; `floor-3` on all five learners over the datasets of 800 genomes or more; the twin on all five over every dataset, and again over those of 800 genomes or more. An arm withdrawn before it ran keeps its pairs in the family at p = 1.

A screen arm's `scm` and `cart` read no screen. They are checked equal to the default's in every fold and are not counted.

Recall stands beside log loss wherever log loss is given.

## Per arm

Each arm is read on its own conditions, against the same default and the same twin.

**`cmh-4000`.** The Cochran-Mantel-Haenszel screen at the release's budget of 4,000.

**`cmh-16000`.** The Cochran-Mantel-Haenszel screen at a budget of 16,000.

**`floor-3`.** The release's chi-square screen at 4,000, with the floor fixed at 3. It is read over the datasets of 800 genomes or more.

## Beside it

These counts stand beside the count that decides, and decide nothing.

**Floor arm over every dataset.** `count-floor-all.tsv` counts the floor arm over the whole first run, the 141 datasets of under 800 genomes entering at a difference of 0. With 141 of 233 differences fixed at 0 its median is 0 by construction.

**No constant model.** `count-no-constant.tsv` leaves out every dataset with a study run refused for no candidate and scored as the constant model (ruling 16). Ruling 16 cannot arise in C4: every arm screens the default's candidates or more.

**Stricter reading.** C3's stricter reading of the floor arm's conditions is written beside ruling 15.

**Species.** The counts per species are dataset and learner pairs over each species' own datasets.

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

## Limits

**Study split.** One study run per dataset at the release's seeds, and no study twin. That departs from the plan's Rules, Draws, which give each arm seeds 42, 43 and 44.

**Binaries.** `reference/screen-arms/page-binaries.tsv` names each species' commit, binary and code sha256, and the restarts its run made.

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

**Reading.** The answers read are the run's own: each species' `score/` tables, copied as its run wrote them.

```
python3 -I scripts/c4-screen/run_checks.py check <copies> scripts/c4-screen/lists.tsv <release> read/checks.txt
python3 -I scripts/c4-screen/run_checks.py gather <copies> scripts/c4-screen/lists.tsv read/checks.txt read
python3 -I scripts/c4-screen/score.py --count read/count read/score/<species>/datasets.tsv   # 23 files
python3 -I scripts/c4-screen/page_tables.py write read scripts/c4-screen/lists.tsv read/checks.txt reference/screen-arms
python3 -I scripts/c4-screen/page_tables.py print reference/screen-arms
```

`read/count/withdrawn.txt` names each arm passed to `score.py --count` with `--withdrawn`, one a line, and is empty when none is. `reading.txt` lists every command with the sha256 of each file it wrote.

**Tables.** `scripts/c4-screen/page_tables.py write` checks the reading against 16 gates and writes `reference/screen-arms/` only when all pass. Where `score.py`'s tables already are the page's, it checks and copies them.

| Gate | What must hold |
|---|---|
| 1. files | the reading holds reading.txt, count/ with score.py --count's four tables and withdrawn.txt, and score/; no score.log is in it; the run's checks.txt says every check passed |
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

**Tests.** `tests/c4_screen/test_page_tables.py` holds SYNTHETIC tests of the table script and `tests/c4_screen/test_run_checks.py` of the check script, on three made-up species. `scripts/c4-screen/page_faults.py` plants each fault below in memory and runs the test meant to catch it, which must pass as the script stands and fail with the fault; `tests/c4_screen/test_page_faults.py` holds every fault to that, and each fault's text to one place in its file.

| Fault | Planted | Gate | Caught by |
|---|---|---|---|
| P1 | a score.log in the reading passes | 1 files | `Gates.test_a_score_log_in_the_reading_is_refused` |
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
