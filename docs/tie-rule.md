# Tie rule experiment

`grm learn --exp-tie-rule exact` changes how the decision tree (cart), the forest (rf) and the boosted trees (gbt) break ties between equally good features. It is off by default, and a run without it writes what `main` writes. The code is `src/learn/tie.rs`.

## The two rules

Production's rule (`float`): among candidates of equal 64-bit gain, a node splits on the lowest-numbered feature. Two splits that are equal on paper can differ in the last bit of their 64-bit gains, so the pick can follow rounding.

The exact rule (`exact`) makes the stated tie rule true for cart and rf; for gbt, only between splits with the same two groups.

- **cart and rf.** A candidate's gain is a fraction of the node's class counts, compared exactly; equal fractions go to the lowest-numbered feature. Which candidates may split is still decided on the 64-bit gain, as in production.
- **gbt.** Gains are sums of gradients, with no exact form. The rule moves a split only to a lower feature with the same two groups: the node takes production's pick, then the lowest-numbered candidate that divides the node's samples into the same two groups. Candidates with equal gains and different groups still follow rounding. In the first round of a logistic fit the gradients are equal within a class, so gbt candidates with equal counts have equal gains; this rule does not treat them as equal. The test `equal_counts_on_other_samples_follow_the_rounded_gain_under_both_rules` in `src/learn/gbt.rs` pins such a case.

A fit under the exact rule checks once, before its first node, that it has fewer than 2^31 resistant and 2^31 susceptible fitting samples, and stops with a message if not (`tie::check_counts`).

## What the review found

From the review of 15 Nov 2026, on the experiment lane's 169 datasets: for cart and rf the exact rule makes the stated rule true and never chose a worse split: all 246,420 cart and rf differences were exact ties. For gbt it moves a split only to a lower feature with the same two groups; candidates with equal gains and different groups still follow rounding. Balanced accuracy and log loss show no difference in a paired test.

## Its reports

With `--exp-tie-rule exact` a report adds `eval_config.tie_rule: "exact"` and the input fingerprint hashes it, so `grm learn` refits it and the metric audit flags it. The site does not read either field, so exact reports are kept out of results directories:

- At the release's seeds `grm learn` writes to `species/<slug>/results/`, the folder `grm report` indexes. An exact run refuses that folder unless `--exp-tie-out DIR` names another one.
- `grm report` stops at a report whose `eval_config.tie_rule` is not production's, and names its file.

`--exp-tie-trace FILE` writes one JSON line per tree-learner fit: its nodes searched, the nodes where the two rules pick different features, and a digest of each tree. The fits are the same as without it.

`--fold-file` runs at the release's seeds and production's rule; it refuses the experiment's options.
