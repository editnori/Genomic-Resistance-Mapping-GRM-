# Input-arm result

The release fits its five learners on the presence of 31-mers in each genome (K31). This page asks whether five other inputs separate the classes better: unitigs, K31 plus k-mer counts, K31 plus amino-acid k-mers of called genes, the AMRFinderPlus program's calls alone, and K31 plus the program's calls. It reads block C3 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/input-arms.md`.

<!-- part 2: the time the reading's rules were fixed, and the time the count was read once. -->
The reading's rules were fixed before any C3 score was read (`docs/assumption-plan-corrections.md`, entries "C3 design" and "C3 arms"). The count is read once, with no other option, subset or measure.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalog. "The program" on this page is the AMRFinderPlus program itself, a separate thing.

## Reading

<!-- part 2: Result, per arm by ruling 15: whether at least 3 of 5 learners carry the change, and each condition read. -->

**Ruling 15.** An input arm touches all five learners, so the plan's conditions are read per learner. A learner carries the change when it meets condition 1 and condition 2 both, and at least 3 of the 5 learners must carry it.

Condition 3 follows the plan's words, the twin meets neither 1 nor 2: it holds when the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. Condition 4 holds for all five learners. Condition 5 holds for each learner that carries. Condition 6 is the maintainer's.

| Condition | Reads by ruling 15, which decides |
|---|---|
| 1 and 2 | a learner carries when its lineage median difference's interval lies wholly above 0 and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry |
| 3 | the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart |
| 4 | no learner's recall interval lies wholly below -0.02, for all five |
| 5 | each learner that carries has a study median of its lineage median's sign |
| 6 | the maintainer rules the change |

<!-- part 2: page_tables.py print, [conditions]: each arm's conditions by ruling 15. -->

**The stricter reading.** Beside ruling 15 stands a stricter reading: condition 1 for at least 3 of 5 learners, 2 and 5 for all five, 3 the twin meeting neither, 4 for all five. It decides nothing, and it stands under Beside it.

"No difference shown" is the rule's reading and no statement that two inputs are equal.

The release keeps its input whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run is L's datasets less those in M, which wait for the maintainer's label rule: 233 datasets in 23 species, 55 of them in S with a study split (`docs/input-arms.md`, Datasets). The two program arms run on A's first run, 204 datasets in 21 species, 54 of them in S: database 2026-08-07.1 of the program has no organism option for *Mycobacterium tuberculosis* or *Staphylococcus haemolyticus*. The section Datasets names every species.

<!-- part 2: the datasets and species that entered, and any species set aside, named with its reason. -->

**Arms.** Each arm is fitted in place of K31 with the release's five learners and scored against K31 refitted in the same run, on the same folds, with the same binary (ruling 17). Each arm is counted on its own datasets only: the three matrix arms, the default and the twin on L's first run, the two program arms on A's.

**Draws.** The default and each arm come from the same runs at fold seeds 42, 43 and 44 on the lineage split. A dataset's value is the mean of its draws. The twin is the default at seeds 45, 46 and 47. It is never read as an arm.

**Study split.** Each dataset in S gets one study run per arm and for the default, on its reviewed study plan, at the release's seeds, with no twin. A study row is never counted on the lineage split.

**Refused study runs.** Only a study run can be refused because an outer fold's training genomes give no candidate feature. By ruling 16, an arm's study run refused this way, where the default's is not, is scored as the constant model on the default's folds: every test genome of a fold gets the majority class of the fold's training genomes, R on a tie. That dataset is kept, enters the count, and is marked in `constant`.

The count with those datasets left out stands under Beside it and decides nothing.

<!-- part 2: the study runs refused, per arm. -->

**Not scored.** A dataset with a draw not scored is not scored. A missing draw is never read as zero.

<!-- part 2: each dataset not scored, named in full, with the reason its run recorded. -->

**Runs.** Each species ran on a machine of its own, one run with one binary for the default, the twin and all five arms, and wrote its scorer's tables there. The count is the scorer's count over those tables, copied as the runs wrote them.

## Count

The primary measure is each learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. Every difference is the arm's value less the default's.

The median difference comes with a 95% interval from 4,000 species resamples. The one-sided sign test of better against worse runs with Holm over 30 lineage pairs on the primary measure: 6 arm rows (unitigs, K31 plus counts, K31 plus amino-acid k-mers, the program alone, K31 plus the program and the twin) by 5 learners. An arm withdrawn before it ran keeps its five pairs at p = 1.

<!-- part 2: page_tables.py print, [count]: per arm and learner on the lineage split, conditions 1 and 2. -->

## Per arm

**Unitigs.** One column per distinct unitig pattern of the graph over every 31-mer of the species' genomes, in place of K31. It runs on L.

**K31 plus k-mer counts.** [count ≥ 2] of each 31-mer in a genome, so a genome has the row when the 31-mer occurs in it at least twice, then the K31 rows less those whose cohort pattern a count row has. It runs on L.

**K31 plus amino-acid k-mers.** Amino-acid 10-mers of the genes Prodigal 2.6.3 calls, then the K31 rows less those whose cohort pattern an amino-acid row has. It runs on L.

**The AMRFinderPlus program's calls alone.** One column per element symbol of the program's Type AMR rows, as the program's table names it. It runs on A. The program is AMRFinderPlus 4.2.7 with database 2026-08-07.1 and `--plus`.

**K31 plus the AMRFinderPlus program's calls.** The program's columns, then the dataset's K31 rows less those whose cohort pattern a program column has. It runs on A.

A program-call arm would enter only with the program's version pinned, and the maintainer would rule how GRM may name its calls.

<!-- part 2: per arm, the reading by ruling 15, its datasets and species, and its builds' counts (page_tables.py print, [arms] and [builds]). -->

## Beside it

These tables are scored beside the reading. They decide nothing.

**The stricter reading.** <!-- part 2: page_tables.py print, [stricter]. -->

**Without the constant model.** The count with the constant-model datasets left out, Holm over its own 30 rows. <!-- part 2: page_tables.py print, [no-constant]. -->

**Other measures.** Recall, specificity, the very major and the major error rate and log loss on the lineage split, and the primary measure on the study split. Lower is better for log loss and both error rates. Recall stands beside log loss wherever log loss is given. <!-- part 2: page_tables.py print, [beside]. -->

## Seen after

These points are seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

<!-- part 2 -->

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Labels.** The 35 datasets of M wait for the maintainer's label rule. The count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged.

**Study split.** One study run per dataset at the release's seeds, and no study twin.

**Constant model.** A study run refused for no candidate is scored as the constant model (ruling 16). It enters the count and is marked. The count without it decides nothing.

**The program.** A program-call arm would enter only with the program's version pinned. A program column is named only as the AMRFinderPlus program's table names it.

**Binaries.** <!-- part 2: one binary by binary.txt's code sha256 in every species, or each binary named with its species. -->

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits, arms and learners only. No file holds a genome, a k-mer, a unitig, a sequence, a position or a column name of a matrix.

**The release.** The release keeps its input whatever the count says, until the maintainer rules.

## Method

**Checks of the run.** Before any score is read, `scripts/c3-inputs/run_checks.py` reads each species' copy and prints counts only. Each check can fail, and the count is taken only when all six pass.

| Check | Holds |
|---|---|
| species | the copies are exactly the species of `scripts/c3-inputs/lists.tsv`'s first run |
| `SHA256SUMS` | every line of each species' `SHA256SUMS` holds, the files under `score/` read in `score.tar` when it is not unpacked, `DONE` names its sha256, and the files under `results/` and `check/` are its lines |
| `DONE` and `binary.txt` | each `DONE` names its species, binary, code, result count and `SHA256SUMS` as `binary.txt` does; no `FAILED` but one named, set aside and older than its species' `DONE` |
| `expect.py` | each `expect.tsv` is `expect.py`'s from the first run, `lists.txt` names the design's seeds, twin seeds and arms, and the result files are exactly its rows |
| refusals | each refusal record is an arm's study run refused for no candidate, with five constant-model rows in `draws.tsv`, one per learner, and no others |
| genome ids | no genome id of the release's species lists in any file of the copies that is read or published |

<!-- part 2: run_checks.py's printout, and each check's result. -->

**Reading.** The answers read are the run's own: each species' `score/` tables, copied as its run wrote them. The count runs once over them:

```
python3 -I scripts/c3-inputs/run_checks.py <back dir> scripts/c3-inputs/lists.tsv <release dir>
python3 -I scripts/c3-inputs/score.py --count read/count read/score/<species>/datasets.tsv   # each species once
python3 -I scripts/c3-inputs/page_tables.py write read scripts/c3-inputs/lists.tsv reference/input-arms
python3 -I scripts/c3-inputs/page_tables.py print reference/input-arms
```

`reading.txt` lists every command with the sha256 of each file it wrote.

<!-- part 2: what the count printed. -->

**Tables.** `scripts/c3-inputs/page_tables.py write` checks the reading against thirteen gates and stops at the first that fails, naming it. The count files are the scorer's own tables: the program counts them again from the species' `datasets.tsv` files (gate 13), copies them as they are, and sums them with every other file in `SHA256SUMS`.

| Gate | Passes when |
|---|---|
| 1. reading | the reading holds `reading.txt`, `count/` with `count.tsv`, `count-no-constant.tsv` and `conditions.tsv`, `score/` and `run/`, and nothing else; an arm `count.tsv` marks withdrawn is one of the five |
| 2. species | `score/` and `run/` each hold exactly the species of the first run, each once, with its files, and each species' rows name it alone |
| 3. twin | the twin's rows are lineage draws at seeds 45, 46 and 47 only, no arm's row is at a twin seed, and `conditions.tsv` reads the five arms and never the twin |
| 4. study | a study row is a study run at seed 42 on the study split, and every lineage row is a draw at 42, 43 or 44 or the twin's |
| 5. M | no dataset of M, or outside the first run, stands in any table; species and antibiotic are `lists.tsv`'s |
| 6. own datasets | each arm's rows are on its own set's datasets only: the default, the twin and the three matrix arms on L's first run, the two program arms on A's, the study split on those in S, and a withdrawn arm on none |
| 7. draws | `datasets.tsv` is the mean of `draws.tsv`'s draws, and a dataset with a draw not scored is not scored, never read as zero |
| 8. constant | a constant-model draw is an arm's study run refused for no candidate, marked in `draws.tsv` and counted in `datasets.tsv`, and its dataset is kept |
| 9. checks | each species' `checks.tsv` holds checks and every one passed |
| 10. builds | each species' builds are `grm-c3-build-1` records of its species, one per kind of the matrix arms run, holding counts only |
| 11. binary | each species' `binary.txt` and `DONE` name the same binary and code; one code sha256 in every species, or the page names each binary with its species |
| 12. conditions | ruling 15's rows decide and the stricter reading's decide nothing, and conditions 1 and 2 as read here off `count.tsv` equal `conditions.tsv`'s |
| 13. count | `count.tsv`, `count-no-constant.tsv` and `conditions.tsv` are what `score.py --count` writes from the `datasets.tsv` files, with the arms `count.tsv` marks withdrawn |

**Tests.** `tests/c3_inputs/test_page_tables.py` holds 40 SYNTHETIC tests of the table script and `tests/c3_inputs/test_run_checks.py` 16 of the run's checks, on three made-up species, and `tests/c3_inputs/test_page_faults.py` 4 of the fault runner. `scripts/c3-inputs/page_faults.py` plants each fault below in a copy of the tree and runs the test meant to catch it, which must fail; `tests/c3_inputs/test_page_faults.py` checks that each fault's text occurs once in its file and that its test exists. All 28 were caught.

| Fault | Caught by |
|---|---|
| 1. A folder outside the reading passed | `Gates.test_a_folder_outside_the_reading_is_refused` |
| 2. A species read twice | `Gates.test_a_species_twice_is_refused` |
| 3. The twin read as an arm in `conditions.tsv` | `Gates.test_conditions_that_read_the_twin_as_an_arm_are_refused` |
| 4. A study row counted on the lineage split, by the gate (4a) or in the species table (4b), or conditions 1 and 2 read on the study rows (4c) | `Gates.test_a_study_run_on_the_lineage_split_is_refused`; `Tables.test_a_study_row_stays_on_the_study_split`; `Tables.test_the_count_rows_carry_conditions_one_and_two_on_the_lineage_primary_rows` |
| 5. A dataset of M counted | `Gates.test_a_dataset_of_m_is_refused` |
| 6. An arm's rows passed on datasets not its own (6a), or a program arm counted on L's datasets in the arms table (6b) | `Gates.test_a_program_arm_outside_a_is_refused`; `Tables.test_each_arm_stands_on_its_own_datasets` |
| 7. A value passed for a dataset with a draw not scored (7a), or a missing draw read as zero (7b) | `Gates.test_a_value_for_a_dataset_with_a_draw_not_scored_is_refused`; `Tables.test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero` |
| 8. A constant-model draw passed without its reason (8a) or off an arm's study run (8b), or a study run refused for no candidate left unmarked (8c) | `Gates.test_a_constant_model_draw_without_its_reason_is_refused`; `Gates.test_a_constant_model_draw_off_an_arms_study_run_is_refused`; `Tables.test_a_constant_model_row_is_kept_and_marked` |
| 9. A failed check of the scorer passed | `Gates.test_a_failed_check_is_refused` |
| 10. A build record holding more than its counts passed | `Gates.test_a_build_record_beyond_its_counts_is_refused` |
| 11. A `DONE` line naming another code passed (11a), a page that leaves a species' binary out passed (11b), or two binaries written as one (11c) | `Gates.test_a_done_line_naming_another_code_is_refused`; `Gates.test_a_page_that_leaves_a_species_binary_out_is_refused`; `Tables.test_two_binaries_are_each_named_with_their_species` |
| 12. The stricter reading written as deciding (12a), or `conditions.tsv` passed without reading it again off `count.tsv` (12b) | `Gates.test_a_page_that_writes_the_stricter_reading_as_deciding_is_refused`; `Gates.test_conditions_that_differ_from_the_count_are_refused` |
| 13. A count passed that is not the scorer's | `Gates.test_a_count_that_is_not_the_scorers_is_refused` |
| R1. A file changed since its sum passed | `RunChecks.test_a_file_changed_since_its_sum_fails` |
| R2. A `FAILED` file not set aside passed | `RunChecks.test_a_failed_file_not_set_aside_fails_and_one_set_aside_is_named` |
| R3. A result file more than `expect.py`'s rows passed | `RunChecks.test_a_result_file_more_than_its_rows_fails` |
| R4. A refused default passed | `RunChecks.test_a_refused_default_fails` |
| R5. A refusal for another reason passed | `RunChecks.test_a_refusal_for_another_reason_fails` |
| R6. `score.tar` left out of the genome-id count | `RunChecks.test_a_genome_id_in_score_tar_is_counted` |

**Files.** The reading's tables are in `reference/input-arms/`. `scripts/c3-inputs/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets

Per species, named in full: the datasets of the first run, those in A and in S, those of M waiting, and better, worse and same on the primary measure by arm on the lineage split. A species' median decides nothing.

<!-- part 2: page_tables.py print, [species]. -->
