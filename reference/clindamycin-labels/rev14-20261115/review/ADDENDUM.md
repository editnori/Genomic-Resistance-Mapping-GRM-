# revclinda addendum: fold plan spread for `staphylococcus_aureus__clindamycin` (brief point 4)

Exploratory. Written before any run of this point and before any score of a seed other than 42 was seen. What I had seen when writing it: lane clinda's own results at seed 42 (its report and `results/`), and my recompute of them (`results/refit/`).

## What I will run

- Binary: `grm` built from `c1f36fc` (branch `owner/main-c1f36fc`) in `/home/user/scratch/revclinda/target-c1f`, Rust 1.98.1, `--locked --offline`. Every call under `unshare -rn`, `RAYON_NUM_THREADS` at most 6.
- Roots, in `/home/user/scratch/revclinda/roots/`: `with-sa` and `without-sa`. Both hold one link per built file of lane clinda's `root-with/species/staphylococcus_aureus/` (counted.json, genomes.txt, matrix, structure.json) and a real empty `results/`; symbolic-link trees of the public root for bindings, bvbrc, catalog, fasta and species-check; and a real copy of the public `curated/`. In `without-sa`, `phenotypes.jsonl` loses the clindamycin rows of the 659 genomes of PMID 34100643 as found by my own recount, and nothing else.
- Check first: `grm learn --species staphylococcus_aureus --only staphylococcus_aureus__clindamycin` at the default seed on both roots. The point goes on only if every out-of-fold score equals lane clinda's "with" and "without" report. If not, I report the difference and stop the point.
- Then the same command with `--seed s` for s = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 on both roots: 11 fold plans per arm with the default (42).

## Measures

Per split (random, lineage) and row (scm, cart, rf, gbt, lasso, known-gene rule), from each run's out-of-fold scores at the 0.5 cutoff: balanced accuracy and recall (specificity shown beside them).

- V1(s): arm "with", all 1,543 samples, fold plan s. V2(s): the same run scored on the 884 samples outside the study. V3(s): arm "without", fold plan s.
- Spread within an arm: minimum, maximum and range (maximum minus minimum) over the 11 fold plans, for V1, V2 and V3.
- Lane's change: V3(42) - V1(42), as in its report.

## What counts as "the change exceeds the spread between fold plans"

For one split, row and measure:

1. Total change (primary): the ranges of V1 and of V3 over the 11 fold plans do not overlap, and V3's range lies on the side the lane's change points to.
2. Training part (secondary): the same test between the ranges of V2 and V3.

The answer per split is the number of the lane's passing rows (`results/refit/rescore-staphylococcus_aureus/rule.tsv`, measures balanced accuracy and recall) that meet test 1, and the number that meet test 2. A split's changes "stand outside the spread between fold plans" when every one of its passing rows meets test 1. Rows that meet test 1 but not test 2 are read as: the total change stands, its training part does not stand clear of the fold plan.

Not decided in advance and so exploratory if done: any other seed count, any other dataset.
