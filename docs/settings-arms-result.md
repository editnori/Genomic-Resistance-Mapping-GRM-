# Settings-arm result

The release fits each learner with its own settings, some of them chosen inside each fit, and calls a genome resistant when its score is at least 0.5. This page asks whether a search of each learner's settings, made inside each outer training fold, or a call cutoff chosen on inner scores, serves better on the outer test folds. It reads block C5 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/settings-arms.md`.

The design was fixed at 04:08Z on 23 Nov 2026 and the lead's rulings on its questions at 05:52Z, before any C5 score was read (`docs/assumption-plan-corrections.md`, entry "C5 design"). The table script, its gates and this page's form were committed before any score of the run was opened.

<!-- part 2: the clock of the reading, as served-pick-result.md gives it. -->

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

"No difference shown" is the rule's reading and no statement that two rules are equal. Recall stands beside log loss wherever log loss is given. The cutoff arm moves no score, so its log loss equals the default's by construction.

## Reading

<!-- part 2: the result per arm, conditions 1 to 5 with their figures, from page-conditions.tsv and page-count.tsv. -->

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

<!-- part B: the fold <k> as the result gives it, from the rescored checks.tsv detail. -->

<!-- part 2: results, choices files and datasets not scored, with their counts. -->

## Count

The primary measure is each learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. Every difference is the arm's value less the default's; for the very major error rate, the major error rate and log loss, lower is better, so better and worse turn round. The interval is a 95% interval of the median over 4,000 resamples of species. The sign test is one-sided, and Holm runs over the 15 lineage rows of the primary measure: the search, the cutoff arm and the twin, each by five learners.

<!-- part 2: the family's 15 rows, from page_tables.py print [family]. -->

## Per arm

<!-- part 2: per arm, split and learner, each measure's median with better and worse, recall beside log loss, from page_tables.py print [measures] and [count]. -->

**Search.** The search's rows decide with the cutoff arm's, as the conditions read them.

**Cutoff.** The cutoff arm's log loss rows are equal by construction and decide nothing. Its other rows decide.

**What the search chose.** Counts only: per learner, the folds that kept trial 0 of the folds searched, the folds where the plain best differs from the pick, the cutoffs chosen as a count per grid point, and the folds where trial 0 or every trial failed inside. A trial is named by its id alone, and no trial's settings are given.

<!-- part 2: page_tables.py print [choices] and [cutoffs]. -->

## Beside it

These counts stand beside the count. They are in no Holm family and decide nothing.

**Plain best, no anchor.** The plain best against the default, counted as an arm is counted, with no Holm figure.

**Without constant models.** The count with every dataset that holds a constant-model draw left out.

**SCM.** On synthetic cohorts the default-anchored rule lost 0.0126 (0.0088 to 0.0163) for SCM on the lineage split (`docs/hyperopt-selection.md`), so SCM's search result is read beside that loss. The synthetic figure decides nothing.

<!-- part 2: page_tables.py print [beside] and [scm]. -->

## Seen after

These points are seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

**Species.** Per species, named in full, the rows are dataset and learner pairs, with their denominator. They decide nothing.

<!-- part 2: page_tables.py print [species], and what else is seen after the reading. -->

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Study split.** One study run per dataset at the release's seeds, and no study twin.

**Labels.** The 14 datasets of M wait for the label rule. The count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged.

**The cutoff arm.** The cutoff arm moves no score, so its log loss equals the default's by construction. A same on log loss for the cutoff arm says nothing about its calls.

**Times.** A time can equal a genome id as text, so the id check counts apart, by file kind and key or label, a token that is the whole value of `seconds_inner` or `seconds_outer` in a choices file or the value on a `User time (seconds)`, `System time (seconds)` or `Elapsed (wall clock) time` line of a drive's `.time` file, and every other hit fails it.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits, arms, folds and learners only.

<!-- part 2: binaries, from page_tables.py print [binaries]. -->

## Method

**Checks of the run.** Before any score was read, `scripts/c5-settings/run_checks.py check` read each species' copy and printed counts only: the first run's species, each once; each species ended `DONE`, a `FAILED` older than its `DONE` named and set aside; each species' `SHA256SUMS`; `score/score.log` kept out of it; each `DONE` against its `binary.txt`, its results and its choices files; one commit and one code sha256; each `lists.txt`; each `expect.tsv` against `expect.py`; every default and twin result and one choices file per outer fold, and no other; the refusal records and the constant-model rows, counted; and the genome ids over every file that will be read or published and every log. `score/score.log` stays out of the record.

*Klebsiella pneumoniae*'s run stopped at its score step and wrote no `DONE` and no `SHA256SUMS`. `run_checks.py check --failed-at-score klebsiella_pneumoniae` holds its record to four conditions instead, each printed with its counts: its `FAILED` is the at score line in the run's time form, with no `DONE` and no `SHA256SUMS`; `logs/run.log` names as many drives ended ok as it had to run, and none ended otherwise; every line of `SHA256SUMS.lead` holds, its own sha256 is `bd0d9b38edcc6b2cfd0a9c1e80827699d015308ddca1177dac094fadeae2cf7c`, and it names no `score/score.log`; every file of the copy, wherever it lies, is named in it, other than `score/score.log`, `back.tar` and the list itself, and every name in it is a file of the copy; and its `binary.txt` names the run's commit and code sha256. Without the allowance, a `FAILED` with no newer `DONE` fails.

**The scorer after the run.** `scripts/c5-settings/score.py` changed after the run and before any score was opened. The run's scorer, at commit `0d987ba`, stopped *Klebsiella pneumoniae* at its score step on one failed check, "the default's result at the same dataset, split and seed is there and available", for polymyxin B on the lineage split at seed 43, whose default result is unavailable: it read as a failed check a draw that C3's scorer reads as not scored. The changed scorer makes each arm's row there not scored, for the default's reason, by handing C3's scorer the same unavailable result it reads for the default, and a choices file beside such a draw fails a check. For a species with no unavailable default draw it writes `checks.tsv`, `draws.tsv` and `datasets.tsv` byte for byte as `0d987ba`'s scorer does, as a SYNTHETIC test holds on made-up species.

<!-- part B: the 7 finished species' three tables, scored again by the changed scorer, against their own SHA256SUMS lines: n of 21 the same byte for byte. -->

<!-- part 2: the checks' printout, from reference/settings-arms/checks.txt. -->

**Reading.** The answers read are each species' three tables as the changed scorer writes them from its copy, the 7 finished species' tables held first to their own `SHA256SUMS` lines, and each species' choices files as its run wrote them. The count runs once over them:

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

<!-- part 2: the faults' printout and the tests' counts. -->

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

<!-- part 2: per dataset, each arm's value and the default's for each learner, the difference, the draws scored, the constant-model draws and a study run refused, from page_tables.py print [datasets]. -->
