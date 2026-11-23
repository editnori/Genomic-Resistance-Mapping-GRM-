# Cap arms result

The release fits each species on 2,000 genomes picked by its selection rule. This page asks whether a model fitted on the top 4,096 genomes a species by the same rule, `cap-all`, does better than the release's 2,000, `cap-2000`, on the default's held-out genomes, beside a twin of the default. It reads block C2 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/cap-arms.md`, on 8 species: *Enterococcus faecium*, *Escherichia coli*, *Klebsiella pneumoniae*, *Mycobacterium tuberculosis*, *Neisseria gonorrhoeae*, *Salmonella enterica*, *Staphylococcus aureus* and *Streptococcus pneumoniae*.

The reading's tables and its gates were committed before any C2 score was read (`docs/assumption-plan-corrections.md`, entries "C2 design", "C2 rulings" and "C2 review").

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

An added genome only trains and is never scored, so every arm is read on the default's held-out genomes. `cap-4000`, the top 4,000, was withdrawn before any fit and holds no result. Every arm, the default included, is fitted here by `grm arm`, so the default's figures here are not the release report's figures and are never set beside them as if they were.

## Reading

"No difference shown" is the rule's reading. It is not a statement that two caps are equal.

Conditions 1 and 2 are read under the lead's ruling 15: a learner carries the change when it meets conditions 1 and 2 both, and at least 3 of the 5 learners must carry it. Condition 3 holds when the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. Condition 4 reads recall. `conditions.tsv` also writes condition 5, which reads the study split's median; the study split decides nothing in C2, so condition 5 is shown and not read. Condition 6 is the maintainer's. The stricter reading `conditions.tsv` writes beside decides nothing.

The release keeps its cap of 2,000 genomes a species and its pipeline's rule whatever the count says, until the maintainer rules. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run: the 143 CAP datasets of the 8 species that are not in M, 39 of them in S (`scripts/c2-cap/lists.tsv`). The datasets of M wait for the label rule. A dataset's genomes are its genomes in the release's lists.

**Arms.** `cap-2000` is the release's 2,000 genomes a species. `cap-all` is the top 4,096 by the release's selection rule on the rebuilt revision 14 table, or every genome the rule ranks when fewer. The twin is `cap-2000` again on slots of its own, scored as an arm of no change. `cap-4000` was withdrawn before any fit.

**Ties at the cut.** Genomes of equal score sort by id, so where the 4,096 cut falls inside a tie of the selection rule, which of the tied genomes are added rests on their ids and not on their calls.

**None-added.** A dataset with no added genome labelled under `cap-all` is none-added: not fitted, and out of `cap-all`'s count, its denominator included, on both splits. The section Datasets names each.

**Draws.** The default's and `cap-all`'s draws are `grm arm` slots 0 to 2 at evaluation seed 42, and the twin's are the default's slots 3 to 5. A dataset's value is the mean of its draws. A missing draw stops the scorer and is never read as zero.

**Unavailable lineage split.** A dataset whose lineage draws `grm arm` refuses past both replacements has no lineage test for any arm. It leaves the lineage count of every arm and the twin alike, and the section Datasets names it with its reason. Its study rows are scored as before.

**Study split.** The datasets in S get one study fit a dataset and arm, `grm learn --fold-file` at seed 42 on the default's `study-whole` test folds. An added genome trains there in every fold or in none, and only when its BioProject is known and no tested genome has it or its k-mer digest.

## Count

The primary measure is the pooled lineage balanced accuracy on the default's outer test folds, per learner. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. The family is `cap-4000`, `cap-all` and the twin against `cap-2000`, Holm over 15 lineage pairs, 3 arms by 5 learners. `cap-4000`'s five pairs stay in the family at p = 1, named withdrawn.

The datasets that add fewer than 20 labelled genomes under `cap-all` stay in the count: the count was fixed with no threshold before the cohorts were listed (the lead's ruling of 12:19Z on 23 Nov 2026).

## Per arm

Recall stands beside log loss wherever log loss is given.

## Beside it

These counts decide nothing.

**Without the few added.** The count over `cap-all`'s datasets less those that add fewer than 20 labelled genomes, made with the scorer's own counting code and Holm over its own rows.

**Study split.** The study split decides nothing. Its fits, one a dataset and arm, are counted beside the lineage count.

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

## Limits

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**No difference shown.** "No difference shown" is the rule's reading. It is not a statement that two caps are equal.

**Added genomes.** An added genome only trains and is never scored, so every arm is read on the default's held-out genomes.

**Withdrawn.** `cap-4000` was withdrawn before any fit. It holds no result, and its five lineage pairs stand in the family at p = 1.

**Ties at the cut.** The 4,096 cut can fall inside a tie of the selection rule, and the tied genomes added are then picked by id.

**The study split decides nothing.** No condition that decides reads it.

**Labels.** The datasets of M wait for the label rule.

**No final model.** No final model is built here. `grm arm` saves no final model, so a cap arm that passed would need a final fit on its cohort per dataset before block T1 could score it.

**The release.** The release keeps its cap of 2,000 genomes a species and its pipeline's rule whatever the count says, until the maintainer rules.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, splits, arms and learners only.

## Method

**Checks of the run.** Before any score was read, `scripts/c2-cap/run_checks.py` read each species' copy and printed counts only: the species, every `SHA256SUMS`, `DONE` and `binary.txt`, one binary, `expect.py`'s rows, the cohorts against the selection rule's counts, genome ids and BioProject accessions, and the refused lineage draws. It reads the tables of each species' `score.tar` in place, without unpacking it, and `score/score.log` is never read or hashed by it. A copy that lacks `lists.txt` takes it rebuilt from `lists.tsv` and the driver's arguments, and the rebuilt text must match its `SHA256SUMS` line.

**Tables.** `scripts/c2-cap/page_tables.py write` checks the reading against its gates and writes `reference/cap-arms/`. `scripts/c2-cap/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

**Tests.** `tests/c2_cap/test_page_tables.py`, `tests/c2_cap/test_run_checks.py`, `tests/c2_cap/test_review_checks.py` and `tests/c2_cap/test_review_faults.py` hold SYNTHETIC tests of the two scripts. `scripts/c2-cap/page_faults.py` plants each fault of the table script in a copy of the tree and runs the test meant to catch it.

A column that reads 0 by construction stays off the printed page: C2's `score.py` writes `constant` 0 on every draw, so that column is not printed, and its count writes no count without constant learners.

## Datasets

Per species, the datasets of the first run, those in S, the cohort sizes, the added genomes, and per dataset the default's and `cap-all`'s figures, with each dataset not fitted named with its reason.
