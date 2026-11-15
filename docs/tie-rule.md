# How the tree learners break ties between equal splits

When two candidate splits are equally good, the decision tree (cart), the random forest (rf) and the boosted trees (gbt) keep one of them by a rule. This page says what the rule is, and how much the figures of release `rev14-274-20261112` move when that rule alone changes. The present rule stays the production rule of evaluation revision 14. Whether a later revision adopts another rule is Layth's decision: the options are at the end, and this page recommends none.

A change under the exact rule is a sensitivity of a figure to the tie rule. It does not show that either rule is more accurate. Every figure is cross-validated on public BV-BRC genomes, and none is a clinical validation. A feature is a row of the matrix, and this page names none.

The files this page quotes are in [`reference/tie-rule/`](../reference/tie-rule/README.md), and a name such as `measures.tsv` below is a file of that folder. The experiment ran on 15 November 2026. An independent review recomputed its figures from the reports with its own code on the same day and corrected parts of its reading; this page uses the corrected reading. `scripts/tie-rule/record.py` wrote the folder's figures from the reports again.

## The two rules

**The present rule**, in production. A candidate split's gain is a 64-bit float. The larger gain wins, and two bit-equal gains go to the lower feature number. One side's class sums are counted and the other's are taken by subtraction, so two splits of exactly equal gain can differ in the last bits, and rounding then decides which one is kept. In the code at commit `3e5b9a2`:

- cart: `fit_tree` in `src/learn/tree.rs` (line 72), whose closure `better` (line 126) compares two candidates.
- rf: `fit_forest` (`src/learn/tree.rs`, line 308) calls `fit_tree` with the features it draws at each node (`max_features`), so the same closure decides among them.
- gbt: `grow` in `src/learn/gbt.rs` (line 228), whose `reduce_with` (line 259) compares two candidates the same way.

**The exact rule**, an experiment only. For cart and rf, a candidate's class-weighted Gini impurity of its two sides is an exact fraction of counts. The lower impurity wins, compared in integers with no rounding, and equal fractions go to the lower feature number. For gbt the gain is a sum of real-valued gradients with no exact form, so the rule takes the present pick and moves the split to the lowest-numbered candidate that divides the node's fitting samples into the same two groups, with the sides the same or swapped. Under both rules the same 64-bit test of the gain decides which candidates may be picked.

The exact rule is the flag `--exp-tie-rule exact` of branch `exp/tie-rule-exact`, which is not on `main`. Arm A ran `main`'s binary at `852ebc9`, and arms B to E the branch's binary at `88af130`. The branch's later commits, up to its head `35e1ae2`, came during the runs and change tests only; the plan asked for a restart on a new binary, and the review's reruns of the exact arm from `35e1ae2` equal these reports on every dataset it reran.

Feature numbers are each learner's own: cart's is the row of the training matrix, rf's the position in the chi-square screen, and gbt's the candidate's row.

## What was run

Each dataset ran in five arms on the release's matrices, folds and seeds: A, `main`'s binary; B, the branch's binary with the flag off, which is the present rule; C, as B with a trace of every fit; D, the exact rule; E, as D with a trace. With the `seconds` fields masked, A equals the released report, B equals A, C equals B and E equals D on 169 of 169 datasets (`reproduction.tsv`). The plan was fixed before any run (`PLAN.md`).

169 of the 274 datasets and 19 of the 23 species were run. The plan takes the species in ascending order of the release's summed fitting time, and its stop rule left the 4 species with the largest sums unrun: `acinetobacter_baumannii`, `escherichia_coli`, `klebsiella_pneumoniae` and `pseudomonas_aeruginosa`, 105 datasets. They differ from the datasets run: at the median, 708 genomes against 600, and 1,285,698 features against 235,169, 5.5 times as many (`datasets.tsv`). Nothing on this page speaks for them.

The lineage split exists for 163 of the 169 datasets, so its denominators are 163, and the random split's are 169. A held-out call is a score at or above 0.5. Each learner scores 137,303 held-out genomes on the lineage split and 139,498 on the random split.

## The effect per learner

Each tree learner's released figures depend on the tie rule on part of the 169 datasets run (cart 61, rf 166, gbt 149), and some tree learner's log loss changes on all 169.

A metric changes when any field of the learner's `metrics` differs between the two rules. A call changes when the score crosses the call threshold one way or the other. From `measures.tsv` and `summary.tsv`:

| Learner | Datasets with a changed metric, lineage / random / either split | Datasets with a changed log loss | Calls changed, lineage / random | Largest change: recall / balanced accuracy / log loss |
|---|---|---|---|---|
| scm | 0 of 163 / 0 of 169 / 0 of 169 | 0 of 169 | 0 / 0 | 0.000 / 0.000 / 0.000 |
| cart | 31 of 163 / 44 of 169 / 61 of 169 | 61 of 169 | 519 / 233 | 0.387 / 0.141 / 6.567 |
| rf | 161 of 163 / 164 of 169 / 166 of 169 | 169 of 169 | 665 / 360 | 0.311 / 0.156 / 0.110 |
| gbt | 130 of 163 / 124 of 169 / 149 of 169 | 165 of 169 | 155 / 33 | 0.090 / 0.045 / 0.088 |
| lasso | 0 of 163 / 0 of 169 / 0 of 169 | 0 of 169 | 0 / 0 | 0.000 / 0.000 / 0.000 |

Over the three tree learners and both splits, 1,965 held-out calls change. scm and lasso, which grow no trees, change on no dataset.

Datasets whose change, in absolute value, is over 0.01 / over 0.05, out of 163 on the lineage split and 169 on the random split:

| Learner, split | Recall | Log loss | Balanced accuracy |
|---|---|---|---|
| cart, lineage | 11 / 5 | 23 / 13 | 12 / 4 |
| cart, random | 12 / 2 | 35 / 23 | 6 / 0 |
| rf, lineage | 22 / 5 | 14 / 1 | 21 / 1 |
| rf, random | 16 / 0 | 10 / 0 | 7 / 0 |
| gbt, lineage | 9 / 2 | 11 / 2 | 7 / 0 |
| gbt, random | 7 / 0 | 1 / 0 | 1 / 0 |

## The paired test

For each tree learner and split, the change of each dataset, exact minus present, is averaged over the datasets, and its 95% interval comes from a bootstrap over species (19 species, 4,000 draws, seed 20261115; `paired.tsv`).

| Learner, split | Recall | Log loss | Balanced accuracy | Specificity |
|---|---|---|---|---|
| cart, lineage | +0.0064 [+0.0013, +0.0125] | -0.0271 [-0.1203, +0.0534] | +0.0004 [-0.0018, +0.0032] | -0.0055 [-0.0112, -0.0004] |
| cart, random | -0.0009 [-0.0020, -0.00002] | -0.0015 [-0.0330, +0.0424] | -0.0004 [-0.0013, +0.0003] | +0.0001 [-0.0015, +0.0018] |
| rf, lineage | +0.0006 [-0.0029, +0.0044] | +0.0006 [-0.0005, +0.0026] | +0.0004 [-0.0014, +0.0024] | +0.0001 [-0.0015, +0.0014] |
| rf, random | -0.00004 [-0.0009, +0.0008] | +0.0002 [-0.0004, +0.0010] | +0.0001 [-0.0004, +0.0007] | +0.0003 [-0.0003, +0.0008] |
| gbt, lineage | +0.0006 [-0.0004, +0.0022] | +0.0001 [-0.0013, +0.0013] | +0.0004 [-0.0003, +0.0014] | +0.0002 [-0.0015, +0.0018] |
| gbt, random | -0.0003 [-0.0009, +0.0001] | +0.0002 [-0.00002, +0.0005] | -0.0002 [-0.0007, +0.0001] | -0.0001 [-0.0008, +0.0006] |

A paired test over datasets, with a species bootstrap, finds no difference in balanced accuracy or log loss for any tree learner on either split. For cart on the lineage split, the exact rule moves recall up by 0.006 [0.001, 0.012] and specificity down by 0.006 [0.0004, 0.011], beside log loss -0.027 [-0.120, +0.053].

3 of the 24 intervals exclude 0: cart's recall and specificity on the lineage split, and cart's recall on the random split. About 1.2 would at 95% by chance. The test shows no difference in balanced accuracy or log loss; it does not show that the two rules are equally accurate.

## The learner order

On 7 of the 169 datasets run, the learner with the highest random-split balanced accuracy changes when only the tie rule changes; the order of any two learners by balanced accuracy changes on at most 8 of 163 datasets (lineage split) and 9 of 169 (random split).

That learner is the index's best learner (`src/report.rs`, `best_learner`). Datasets where the order of two learners by balanced accuracy changes with the tie rule alone (`order.tsv`):

| Pair of learners | Lineage split, of 163 | Random split, of 169 |
|---|---|---|
| scm, cart | 4 | 8 |
| scm, rf | 2 | 3 |
| scm, gbt | 2 | 3 |
| scm, lasso | 0 | 0 |
| cart, rf | 7 | 5 |
| cart, gbt | 8 | 7 |
| cart, lasso | 7 | 9 |
| rf, gbt | 3 | 5 |
| rf, lasso | 1 | 3 |
| gbt, lasso | 1 | 4 |

The tie effect of a dataset and split is the largest change in balanced accuracy over cart, rf and gbt. Under the present rule, the gap from the best learner to the second has a median of 0.0070 on the lineage split and 0.0047 on the random split, and the tie effect a median of 0.0011 and 0.0009. The tie effect is above 0 and at least the gap on 48 of 163 datasets (lineage split) and 46 of 169 (random split) (`splits.tsv`). There, the difference between the two best learners is no larger than what the tie rule alone moves.

## Feature numbering

At 52,749 of the 71,519 cart nodes searched (73.8%) and 6,156,787 of the 18,639,382 rf nodes (33.0%), more than one candidate shares the best exact impurity. There the exact rule takes the lowest feature number, and the present rule takes it too unless rounding departs, so under both rules the split depends on the order of the matrix rows. Rounding departs at 1,761 cart nodes (2.5%) and 244,659 rf nodes (1.3%). Each of them is an exact tie, and at none did rounding put a strictly worse split first (`traces.tsv`, the present trees).

So the experiment measures where rounding departs from the exact order. It does not measure how much the figures depend on the order of the matrix rows.

## gbt

For gbt, the exact rule moves a split only to a lower feature with the same two groups. It moves 10,602 of the 997,991 nodes searched in the logistic models, in 1,593 of 1,829 fits, and 7,044 of 343,562 in the MIC models, in 364 of 365 fits (`traces.tsv`).

Two gbt candidates with equal gains and different groups still follow rounding. In the first round of a logistic fit every sample of a class has the same gradient, so two candidates on equal numbers of resistant and susceptible samples have equal gains whatever samples they hold. The review built such a case, and both rules split on the higher feature number. So for gbt the exact rule does not make the stated rule true.

## The options

The review gives five options, the last two of them its own. This page recommends none.

- **Keep the present rule and merge the flag.** Needs a Methods sentence that ties are broken on 64-bit gains, so the lower feature number wins only when two gains are equal to the last bit; tests for the code that passes the rule and the counts to the fits, where the review planted faults that every test passed and that change real reports; exact reports kept out of results directories, because the site and the index read neither `tie_rule` nor the input fingerprint; and the branch binary's lasso time solved first. With the flag off that binary writes `main`'s reports, but took a median of 3.02 times `main`'s lasso seconds per dataset, and 1.17 times its cart seconds, on the run's shared machine (`seconds.tsv`). The review measured the lasso slowdown again on an idle machine.
- **Keep the present rule and leave the branch unmerged.** Needs the Methods sentence only. The measurement stays in this record.
- **Adopt the exact rule at a later revision.** Needs `EVALUATION_REVISION` raised and all 274 datasets fitted again, the 105 not run included; a decision on gbt candidates with equal gains and different groups; the site and the index to show `tie_rule` wherever exact and present reports can be read side by side; and the time cost measured on an idle machine. On the run's shared machine the exact rule took 1.18 (cart), 1.11 (rf) and 1.25 (gbt) times the present rule's seconds, median per dataset.
- **Keep the present rule and report the tie effect as a stated sensitivity**, beside the gap between the two best learners. Needs a Methods paragraph and a field or table per dataset (`splits.tsv` holds one for the datasets run), and no refit.
- **A rule that does not depend on feature numbers.** Both rules take the lowest number at 73.8% of cart and 33.0% of rf nodes searched. Needs its own design and experiment.

## Limits

- 169 of the 274 datasets were run. The 105 datasets of the 4 largest species were not, and they have more features.
- A change under the exact rule is a sensitivity to the tie rule, not a difference in accuracy, and the figures are internal cross-validation on public BV-BRC genomes.
- For gbt, two candidates are equal only when they make the same two groups, as above.
- Which candidates may be picked is the same 64-bit test under both rules.
- The traces count nodes per fit and hold no record per node, so a single node cannot be looked up.
- The seconds were measured on a machine shared with other work, and are rough.
- The reports are not in the repository. `inputs.sha256` gives the sha256 of each report and trace the record was written from.

## Files and commands

- [`reference/tie-rule/README.md`](../reference/tie-rule/README.md) lists each file of the record, what it holds and how it was written.
- `python3 -I scripts/tie-rule/record.py RUNS RELEASE PLAN_DIR --check reference/tie-rule` checks the record against the reports where they are unpacked. CI does not have them.
- `tests/tie_rule` runs the script on SYNTHETIC reports written in the test, holds every figure of this page to the record, and checks the record's files against each other and against `SHA256SUMS`.
