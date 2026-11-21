# Plan: what an exact tie rule changes in cart, rf and gbt

Lane `tierule`, fixed at 2026-11-15 12:24Z, before any run of the branch, any flag-on fit and the reproduction check. Its sha256 is in `PLAN.sha256` and in `REPORT.md`. Anything decided after a result is read is labeled exploratory in `REPORT.md`.

## Question

Does any released figure of `rev14-274-20261112` (evaluation revision 14, built from commit `469a4f9`, whose `src`, `Cargo.toml`, `Cargo.lock`, `.cargo` and `build.rs` equal `main` at `852ebc9`) depend on how the decision tree (cart), the forest (rf) and the boosted trees (gbt) break ties between splits of equal gain?

## The two rules

Present rule. Each candidate's gain is a 64-bit float, and the larger gain wins, with bit-equal gains going to the lowest-numbered feature. One side's class sums are counted and the other's are taken by subtraction, so two splits of exactly equal gain can differ in the last bits, and rounding then decides which one is kept.

Exact rule, cart and rf. A candidate's class-weighted Gini impurity of its two sides is the exact fraction (a·b·D2 + c·d·D1)/(D1·D2), with (a, b) the present side's resistant and susceptible counts, (c, d) the absent side's, P and N the fit's class counts, D1 = a·N + b·P and D2 = c·N + d·P. The lower impurity wins, compared as 128-bit integers with a 256-bit cross product, and equal fractions go to the lowest-numbered feature. Which candidates are admitted (64-bit gain above 1e-12) is unchanged.

Exact rule, gbt. Two gbt candidates are equal when they divide the node's fitting samples into the same two groups, with the same present side or the sides swapped. The gbt gain is a sum of floating gradients, so its exact value is not a count fraction, and this partition test is the definition of "equal" used here. The exact rule takes the present pick and moves the split to the lowest-numbered of the tree's candidates with the same partition, recomputing that feature's gain.

Feature numbers are the learner's own: cart uses the training-candidate matrix row, rf the position in the chi-square screen, gbt the candidate row id. No feature is named.

## The flag

`grm learn --exp-tie-rule float|exact`, default `float`. `--exp-tie-trace FILE` writes one JSON line per cart tree, rf forest and gbt model with its tally (nodes searched, nodes with a tie at the best exact gain, nodes where the picks differ) and a digest of each tree. A traced fit splits on the configured rule, so it is the same model as the untraced fit. Off, the report and its input fingerprint are unchanged. On, the report records `eval_config.tie_rule` and the fingerprint hashes it. `EVALUATION_REVISION` stays 14.

Branch `exp/tie-rule-exact` at `88af130`. Binaries, built with Rust 1.98.1 under `CARGO_TARGET_DIR=/home/user/scratch/tierule/target`:

- main `852ebc9`: `/home/user/scratch/tierule/bin/grm-main-852ebc9`, sha256 `8ab7cbff1145cc0f4c4cc51dfa659fdf00bf51d6bf127a9d947d3c29c71fb5e9`
- branch `88af130`: `/home/user/scratch/tierule/bin/grm-branch-88af130`, sha256 `f2423e4f2186ef7a38602e79c5e6d484ad1373c03c85bb5577621f3afafaf8b5`

If the branch changes before the runs end, the runs restart on the new binary and both binaries are listed.

## Datasets, in a fixed order

Rule: species in ascending order of the release's summed per-dataset `seconds` (ties by slug), and datasets in alphabetical order of id within a species.

Tier 1, run in full: the 12 species with the smallest sums, 51 datasets.

1. streptococcus_agalactiae (3)
2. vibrio_cholerae (3)
3. staphylococcus_haemolyticus (3)
4. klebsiella_quasipneumoniae (3)
5. klebsiella_oxytoca (4)
6. shigella_flexneri (4)
7. clostridioides_difficile (4)
8. corynebacterium_diphtheriae (3)
9. klebsiella_aerogenes (3)
10. campylobacter_coli (8)
11. campylobacter_jejuni (4)
12. enterobacter_cloacae (9)

Tier 2, in this order, as time allows: shigella_sonnei (20), neisseria_gonorrhoeae (7), enterococcus_faecium (14), staphylococcus_aureus (18), streptococcus_pneumoniae (12), mycobacterium_tuberculosis (22), salmonella_enterica (25), pseudomonas_aeruginosa (13), acinetobacter_baumannii (18), klebsiella_pneumoniae (38), escherichia_coli (36). A Tier 2 species starts only after Tier 1 is done and before 17:00Z. No new dataset starts after 18:30Z. A dataset counts only when all five arms below finished. The denominator of every share is the number of datasets that counted, stated with each figure. The finding covers those datasets only.

## Inputs, folds and seeds

The release's inputs: the same FASTA, genome lists (set and order equal to the release's species files) and BV-BRC feature cache, rebuilt with the main binary by `grm count`, `grm structure` and `grm matrix` (k = 31, minimum support 2). Folds and seeds come from the release's evaluation settings (seed 42), so every arm uses the same matrix, folds and seeds. Every run is under `unshare -rn` (no network), with `RAYON_NUM_THREADS=8`, or 5 while both other lanes build or fit.

## Arms, per dataset, in this order

Each arm is `grm learn --species <slug> --only <id> --force` on one root, and its report is moved to `results/runs/<arm>/<id>.json` afterwards.

- A: main binary.
- B: branch binary, flag off.
- C: branch binary, flag off, traced.
- D: branch binary, `--exp-tie-rule exact`.
- E: branch binary, `--exp-tie-rule exact`, traced.

## Reproduction, before any exact outcome is read

Reports carry wall-clock `seconds`, so equality is byte equality after every `"seconds":<number>` is masked. The first dataset (streptococcus_agalactiae__clindamycin) runs arms A and B alone first. If A or B differs from the released report, the work stops and `REPORT.md` says what differs. For every dataset: A must equal the release, B must equal A, C must equal B and E must equal D, all masked. A dataset whose A differs from the release is left out of the measures and listed. Any B that differs from A ends the flag-off claim and is reported.

## Counts, per learner, per split (lineage first)

From the traces of C (the present trees) and E (the exact trees):

- nodes searched;
- nodes where the exact pick differs from the present pick: for cart and rf split into equal exact gain (a true tie decided by rounding or feature number) and lower exact gain (rounding put a worse split first); for gbt, nodes where a lower-numbered candidate gives the same two groups;
- trees that differ (C against E digests, position by position);
- fold models that differ: the cart tree, rf forest and gbt model scored on each outer test fold, counted apart from the inner depth-tuning trees, the MIC models and the final models;
- held-out calls that change (score at or above 0.5 against below), from the `scores` of B and D, per learner and split.

A fit record is labeled by its place in the fixed order of `grm learn` (split, fold, learner, inner or outer, then MIC, then final), checked against its learner, seed and fit-sample count. scm and lasso are controls: they must not change.

## Released measures, both ways (B against D)

For each of cart, rf and gbt apart, for each split with the lineage split first:

- recall (`sensitivity`) beside log loss, the mean of −[y ln p + (1 − y) ln(1 − p)] over held-out genomes with p the score clamped to [1e-12, 1 − 1e-12];
- the rest of `metrics`: balanced accuracy, specificity, precision, F1, MCC, VME, ME, AUROC, AUPRC, accuracy and the four counts;
- `fold_balanced_accuracy`, `no_call` and `called_metrics` where present;
- the report's `mic` block and the final models (`models.cart`, `models.cart_depth`, `models.rf_importance`, `models.gbt_importance`), equal or not.

A dataset changes a released measure for a learner if any of these fields differ in either split. For each learner and measure: the number of datasets that change, out of the datasets run; the median and largest absolute change over all datasets run (zeros included) and over the changed datasets alone. Per learner and in total.

## Time cost

The learners' `seconds` of B against D per learner, summed over datasets, with the median per-dataset ratio, and the wall-clock of each untraced arm. Traced arms are not timed. The machine is shared, so the cost is read as rough.

## What the result means

"The release's figures do not depend on the tie rule" means: for cart, rf and gbt each, zero of the datasets run change any released measure above in either split. If trees or fold models differ but no measure changes, the reading is that the trees depend on the rule and the released figures do not. If any measure changes, the reading is a dependence, sized by the counts above. A result on the smaller species does not speak for the species left unrun. Either way the default stays as it is: the decision belongs to the project owner.
