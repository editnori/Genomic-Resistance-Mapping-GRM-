# Learner overlap on the lineage split (revision 14)

## Plan

This plan was fixed on 24 November 2026, before the table script was written and before any report was read for this question. It fits nothing, trains nothing and scores no genome again. Every figure comes from the held-out scores and labels that the 274 revision 14 release reports already store.

### Question

Each revision 14 report stores, for its lineage split, every learner's held-out score for every genome. When the learners miss a resistant genome, do they miss the same one? If the five learners miss largely the same resistant genomes, a sixth learner of any shape has little room to lift recall on them, and for those genomes the data limits recall. If they miss different genomes, a combination of the five has room. The lineage split is a retrospective split inside one database, BV-BRC. It is not external validation.

### Data

- The 274 report files of the revision 14 release's `datasets/` directory, one per species and antibiotic.
- Per report: `samples[].resistant` (the label), `samples[].catalog_resistant` (GRM's own reader of the AMRFinderPlus catalogue, the known-gene rule) and, in the split whose `scheme` is `lineage`, each learner entry's `learner`, `scores`, `metrics`, `no_call` and `fold_complexity`, and the split's `test_fold`. `counts.samples` and `counts.resistant` are checked against `samples`.
- The random split is not read. A `samples` entry also holds a genome id. The table script drops every key named `id` as a file is parsed, so no genome id is held, printed or written.
- The cutoff is the release's own. A score of 0.5 or more is a resistant call and a score under 0.5 is not. No other cutoff is tried.

### Learners

The learners are the lineage split's entries other than `catalog`, which is GRM's reader of the AMRFinderPlus catalogue and not a learner. A first look at one report's keys found `scm`, `cart`, `rf`, `gbt`, `lasso` and `catalog`. No combined or mean-vote entry stands in the split. The read prints the learner names it finds on every report. The five are the name set found on the most reports. A report with another set is listed by species and antibiotic, keeps its row and stays out of the pooled rows. `catalog` enters only through `catalog_resistant` and through check 9.

### Table `reference/learner-overlap/datasets.tsv`

One row per report used, an aggregate that keeps its genome count. Columns, in order:

- `species`, `antibiotic`: as the report names them, then `genomes`, `resistant` and `susceptible`;
- `learners`: the learner names used, comma-joined, and `learner_count`;
- `miss_0` to `miss_5`: resistant genomes scored under 0.5 by exactly k of the learners, which add to `resistant`;
- `wrong_0` to `wrong_5`: susceptible genomes scored 0.5 or more by exactly k of the learners, which add to `susceptible`;
- `miss_scm`, `miss_cart`, `miss_rf`, `miss_gbt`, `miss_lasso`: each learner's resistant genomes under 0.5;
- `miss_all5`: resistant genomes missed by all five (`miss_5`);
- `miss_all5_no_catalog_call`: of those, the ones `catalog_resistant` leaves without a resistant call;
- `miss_all5_in_null_fold`: of those, the ones in an outer fold where at least one of the five had a null fit (`fold_complexity` 0). Such a fold scores every genome at its training prevalence (docs/equal-scoring.md);
- `null_fold_pairs`: learner and fold pairs with a null fit;
- `expected_all5`: the yardstick, `resistant` times the product of the five miss rates `miss_<learner> / resistant`. It is the count the five would miss together if they missed independently. It is a yardstick. It is neither a model nor a test;
- `ratio_all5`: `miss_all5 / expected_all5`, blank when `expected_all5` is 0.

For a report with fewer than five learners, the `miss_k` and `wrong_k` columns above its learner count, `miss_all5` and the columns after it are blank. Every decimal is written with 12 places.

### Table `reference/learner-overlap/pooled.tsv`

The same counts summed over the pooled reports: those used, carrying the five, with 10 or more resistant genomes. One row with `scope` `all`, then one row per species with `scope` `species`. Beside the sums:

- `datasets`, `datasets_with_ratio` (`expected_all5` over 0) and `datasets_expected_zero`, the reports where one learner misses no resistant genome, reported apart and never dropped silently;
- `ratio_median`, `ratio_q1`, `ratio_q3`: the median and quartiles of `ratio_all5` over the reports with a ratio, each report weighted once, none by its genome count (quartiles by Python's `statistics.quantiles` with `n=4` and `method="inclusive"`; with one report the three are its ratio);
- `datasets_ratio_over_1`: the reports whose ratio is over 1;
- `share_all5_median`: the median over the pooled reports of `miss_all5 / resistant`;
- `share_some_median`: the median over the pooled reports of `(miss_1 + miss_2 + miss_3 + miss_4) / resistant`, the resistant genomes at least one of the five calls and at least one misses.

### Checks

Each prints counts only. A check marked stop ends the read with no table written.

1. Every report file parses, and the files read equal the files found (stop). The files found must be 274 (stop).
2. For every learner entry of every report used, `len(scores) == len(samples)`, each score is a number from 0 to 1, and `len(test_fold) == len(samples)` (stop).
3. `counts.resistant` equals the labels counted in `samples`, and `counts.samples` equals its length (stop).
4. The learner name set is the same on every report used, or the reports that differ are listed by species and antibiotic and counted.
5. Reports with no lineage split, or with fewer than 2 learners, are skipped and counted, and the skipped plus the used equal the files found (stop). A report with two lineage splits stops the read.
6. No genome id stands in any written file: every file the read and the write make is scanned for a token of 2 to 7 digits, a point and 1 to 9 digits, standing alone (the shape of a BV-BRC genome id, a taxon number and a genome number), and for a key named `id` or `samples`. It must find 0 (stop; the files stay unwritten).
7. Each learner's counts at 0.5, recounted from its scores, equal its stored `metrics` `tp`, `fn`, `fp` and `tn` (stop). This proves the scores line up with `samples` and the cutoff reads as the release reads it.
8. No species and antibiotic pair stands on two reports (stop).
9. Where a `catalog` entry stands, its calls at 0.5 equal `catalog_resistant` on every genome (stop).
10. Each learner's genomes in its null-fit folds, by class, are compared with its `no_call`. The rows that differ are counted and the read goes on.

### Reading rule

The statistic is the median over the pooled reports of `ratio_all5`. The reading rule, fixed here:

- a median of 10 or more: the misses stack. The five learners miss largely the same resistant genomes, and for those genomes the data limits recall.
- a median over 1 and under 10: the five share misses more often than the yardstick, by under tenfold.
- a median of 1 or less: no stacking above the yardstick is shown.

Beside it, `share_some_median` against `share_all5_median` says whether more of a typical report's resistant genomes sit where some of the five call them, where a combination has room, or where all five miss. No p-value is computed and nothing here is a test of a null hypothesis. The page will not say that any model would fail, and it names no gene, variant or position.
