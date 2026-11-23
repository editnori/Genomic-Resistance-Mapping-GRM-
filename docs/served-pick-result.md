# Served-learner result

The release serves, per dataset, the learner with the best random-split balanced accuracy, a tie going to the later learner. This page asks whether a learner picked inside each outer training fold, from inner folds grouped by lineage, serves better on the outer test folds. It reads block C9 of the assumption plan (`docs/assumption-plan.md`), with the method of `docs/served-pick.md`.

The reading's rules were fixed at 20:36Z on 22 Nov 2026, before any C9 score was read (`docs/assumption-plan-corrections.md`, entry "C9 design"). The count was read once, with no other option, subset or measure.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. The margin of 0.02 is an adapted research screen. It is no applicable FDA validation criterion.

The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism, and the inner pick's cannot. Wherever the two are compared on this page, a same or a worse for the inner pick is never read as showing that the release's pick is right (ruling 9).

## Reading

[AFTER THE READING: the answer first, condition by condition, for the inner pick and for the twin, from `page-conditions.tsv` and `page-count.tsv`.]

The release keeps its rule whatever the count says, until the maintainer rules. Condition 6 is the maintainer's. Nothing changes in the release or on the site.

## Data

**Datasets.** The first run: the datasets of L, the reports whose lineage split ran, less those in M, which wait for the maintainer's label rule. The section Datasets names them.

**Arms.** The default serves the release's pick on every outer lineage fold, read from the same run's random split at that run's seed. The inner pick serves, per outer fold, the learner with the best balanced accuracy over 3 inner folds of that fold's training genomes, each lineage whole. Both arms are scored with the release's own outer fits.

**Draws.** The default and the inner pick come from the same runs at fold seeds 42, 43 and 44. A dataset's value is the mean of its draws. The twin is the default at seeds 45, 46 and 47.

**Study split.** The first-run datasets in S get one study run each, on their reviewed study plans, at the release's seeds, with no twin. Their default is the seed-42 lineage run's pick. One study run and no twin is a deviation from the plan's Rules, Draws.

**Runs.** Each species ran on a machine of its own and wrote its scorer's tables there, under the answers read.

## Count

The primary measure is the served learner's pooled lineage balanced accuracy on the outer test folds. A dataset is better, worse or the same by more than 0.02, compared at 12 decimals. The median difference is the inner pick's value less the default's, with a 95% interval from 4,000 species resamples read at the 101st and 3,900th values. The one-sided sign test of better against worse runs with Holm over two pairs on the lineage primary measure: the inner pick and the twin, each against the default.

[AFTER THE READING: the count's tables.]

## Beside it

Three answers of the design's questions 3 and 5 were scored beside the read ones, from the same results, each with its own count. They decide nothing.

| Answer | Question | Changes |
|---|---|---|
| `--tie earlier` | 3. Ties | an inner tie goes to the earlier learner |
| `--tie default` | 3. Ties | an inner tie goes to the default's learner when it is tied |
| `--default seed42` | 5. Default | the default serves the seed-42 pick at every draw |

[AFTER THE READING: the beside counts and the pick tables.]

## Seen after

These points were seen in the tables after the reading. They change no rule and no reading above, and decide nothing.

[AFTER THE READING.]

## Limits

**Ruling 9.** The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism, and the inner pick's cannot. A same or a worse for the inner pick is never read as showing that the release's pick is right.

**Study split.** One study run per dataset at the release's seeds, and no study twin. That departs from the plan's Rules, Draws, which give each arm seeds 42, 43 and 44.

**Binaries.** [AFTER THE CHECKS: the species per binary, and Klebsiella pneumoniae's seed-42 equality with its first run.]

**Labels.** The datasets of M wait for the maintainer's label rule, and the count over them is taken when that rule lets them run, with these rules, this scorer and these arguments unchanged. On the first run three of the four label rules give the same labels.

**One database.** Every split is a retrospective split inside one public database. A split by lineage or by study is not external validation.

**Conditions.** No marker is found, named or ranked.

## Method

[AFTER THE READING: the checks of the run, the reading's commands, and the gates.]

**Files.** The reading's tables are in `reference/served-pick/`. `scripts/c9-pick/page_tables.py write` wrote them from the run's own scorer tables and the three scorings beside them, after eight gates passed. `scripts/c9-pick/page_tables.py print` prints this page's tables from those files, rounded to 3 decimals, and `tests/c9_pick/test_page.py` checks that each printed row stands here.

## Datasets

[AFTER THE READING: per species and per dataset, every species and antibiotic in full.]
