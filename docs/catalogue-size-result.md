# Catalogue size result

Is the set of ways *Mycobacterium tuberculosis* can resist an antibiotic small for some antibiotics and large for others, and does the WHO catalogue's size show it? The lineage curve's page (`docs/lineage-curve-result.md`) counts, per antibiotic, the distinct graded variants of the WHO catalogue among the resistant genomes (`catalogue_variants` in `reference/lineage-curve/page-endpoints.tsv`). An antibiotic studied for decades may have both a large catalogue and many resistant genomes on deposit, so that count may record attention rather than how many ways the antibiotic can fail. This page asks whether the count tracks the resistant share of an antibiotic's genomes or only how many resistant genomes there were to study, and says so plainly where the numbers cannot separate the two.

It fits nothing, trains nothing and scores no genome again. It reads the lineage curve's committed table beside three counts of each release report: its genomes, its resistant genomes and its susceptible genomes.

## Rule

This rule was written and committed before the release's counts were read, and it is not changed after.

- **Antibiotics.** The rows of `reference/lineage-curve/page-endpoints.tsv` whose `catalogue_variants` is 1 or more and whose `in_M` is 0. A row with `in_M` 1 waits on the maintainer's label rule: it stands in the table and is left out of every correlation. A row the tool does not read, or with no graded variant, stands in the table and is left out too. Every correlation of a set uses the same antibiotics.
- **Pairs.** The Spearman rank correlation of `catalogue_variants` against each of: the resistant share (`counts.resistant` over `counts.samples` of the antibiotic's release report), the resistant genomes (`counts.resistant`) and the genomes (`counts.samples`). Tied values take the average of their ranks, and the correlation is Pearson's over the ranks.
- **Leave one out.** For each pair, the correlation over the set with one antibiotic dropped, for each antibiotic in turn: the smallest and the largest, each with the antibiotic dropped.
- **Decision.** Let R be the correlation with the resistant genomes and S the one with the resistant share, each in absolute value and unrounded. If R is as large as S or larger, the answer is that catalogue size tracks how many resistant genomes there were to study, and target size is not separable here. If S is larger than R by 0.15 or more, the answer is that the number of graded variants tracks the resistant share, in one species, one catalogue and one database, and no causal reading follows. If S is larger than R by less than 0.15, the answer is that these numbers cannot separate the two. A difference within 0.000000001 of a bound counts as meeting it.
- **Beside.** The correlation with the genomes is reported and decides nothing. The same three pairs over the antibiotics of the set that the tool reads by their own drug (kind `primary`) are reported and decide nothing: a proxy counts another drug's graded variants.
- **Not done.** No p-value, no test, no fit, no fitted line and no model.

## Plan

**Table.** `reference/catalogue-size/antibiotics.tsv`, one row per report of the species in `page-endpoints.tsv`, species and antibiotic in full as that file names them: `species`, `antibiotic`, `in_M`, `kind`, `genomes`, `resistant`, `susceptible`, `resistant_share`, `resistant_with_result`, `catalogue_variants`, `genes_named`, `without_determinant`, `share_without`, `in_correlation`, `left_out_because`. `genomes`, `resistant` and `susceptible` are the report's `counts.samples`, `counts.resistant` and `counts.susceptible`; `resistant_share` is `resistant` over `genomes`. `in_M`, `kind` and the five columns from `resistant_with_result` to `share_without` are carried from `page-endpoints.tsv` unchanged.

**Correlations.** `reference/catalogue-size/correlations.tsv`, one row per set and pair: `set` (`rule` or `primary only`), `x`, `y`, `n`, `rho`, `loo_min`, `loo_min_without`, `loo_max`, `loo_max_without`, `decides`, `comparison_gives` and `antibiotics`, the set listed.

**Guards.** Each prints counts only and stops the run when it fails:

1. every antibiotic of `page-endpoints.tsv` is matched to exactly one release report by its dataset name; the unmatched on either side are listed and counted;
2. each report's `counts.resistant` and `counts.susceptible` add to its `counts.samples`;
3. the antibiotics of every correlation of a set are the same, listed, with n printed beside each correlation;
4. the rank correlation gives the known answers of three hand-made cases, one with ties, before it touches a real number;
5. no genome id of the release stands in any written file: a scan prints 0 hits;
6. for each report the tool reads, `counts.resistant` equals the page's `resistant_with_result` plus `resistant_no_result`.

A set of under 3 antibiotics, or a column whose values are all tied, stops the run.
