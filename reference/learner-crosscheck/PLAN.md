# Lane crosscheck: plan, fixed before any comparison is read

Written on 14 Nov 2026 before any learner was fitted on a panel. Its sha256 is in `PLAN.sha256` beside it.

## Question

Does each of GRM's five learners, as built from `main` at `c435000` (evaluation revision 14), give the result a standard version of the same method gives on the 16 saved test panels?

This is a check of software correctness and reproducibility on small test panels. It says nothing about clinical performance and nothing about how well any learner predicts resistance.

The counts recorded for October 2026 (465 of 480, 5,915 of 6,000 and 777 of 800 trees) are history. They are no target, and no rule below was chosen to reach them.

## What was read and done before this plan

- Read: `src/learn/{mod,scm,tree,gbt,lasso,eval}.rs`, `docs/lasso-stopping.md`, `docs/ci.md`, the header and the `fixture` lines of `tests/data/lasso-stopping-fixtures.txt`.
- Built `grm` at `c435000` with Rust 1.98.1 (release, rc 0).
- Installed the reference libraries from hash-checked files (below) and read the source of pyscm-ml 1.1.2.
- Tried each library's interface on made-up arrays (`toy/api_probe.py` in the scratch). No panel was decoded into a learner and no learner was fitted on a panel.

## Panels

The 16 panels of `tests/data/lasso-stopping-fixtures.txt` (sha256 `7f354884a5dec09bb7e87672d56842f0c3ebd1aa5ef785e4479a83c9a1f99550`), all of them:

| Panel | Samples | Features |
|---|---:|---:|
| `campylobacter_jejuni__ciprofloxacin_lineage_fold0` | 508 | 128 |
| `campylobacter_jejuni__ciprofloxacin_random_fold0` | 508 | 128 |
| `constant_balanced` | 80 | 4 |
| `constant_imbalanced` | 80 | 4 |
| `enterobacter_cloacae__tobramycin_lineage_fold0` | 99 | 128 |
| `enterobacter_cloacae__tobramycin_random_fold0` | 99 | 128 |
| `imbalanced_5pct` | 200 | 16 |
| `marker_dup_constant` | 120 | 12 |
| `one_class_train` | 60 | 4 |
| `planted_logistic` | 300 | 40 |
| `staphylococcus_haemolyticus__fusidic_acid_lineage_fold0` | 114 | 128 |
| `staphylococcus_haemolyticus__fusidic_acid_random_fold0` | 114 | 128 |
| `streptococcus_agalactiae__clindamycin_lineage_fold0` | 141 | 128 |
| `streptococcus_agalactiae__clindamycin_random_fold0` | 141 | 128 |
| `tied_counts` | 100 | 3 |
| `xor_balanced` | 80 | 2 |

Each panel gives a label per sample (`y`, resistant = 1), a training mask (`fit`) and one presence row per feature. Every learner is fitted on the samples of `fit`, called T below, and scored on all samples. The file's `point` lines (the saved lasso paths) are not used.

Both sides decode the file on their own: the Rust harness and a Python decoder. The run stops if their sample counts, label counts, mask counts or per-feature presence counts differ.

## GRM side

A Cargo example, `learner-crosscheck`, calls the public functions of `grm::learn` and writes every fitted model as JSON. `src/` is not changed, so the learners are those of `c435000`. Settings are the production defaults of `EvalConfig`, `ForestParams` and `GbtParams`:

| Learner | Call | Settings |
|---|---|---|
| Set covering machine | `scm::fit` | both kinds, each p of `P_GRID` (0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0), `max_rules` 6: 14 fits per panel, 224 in all |
| Decision tree | `tree::fit_tree` | `max_depth` 8, `min_leaf` 1, every feature at every node; scores at each depth 1 to 8 |
| Random forest | `tree::fit_forest` on the features in the order `screen` gives | 200 trees, subsample 0.632 without replacement, `max_depth` 30, `min_leaf` 1, ceil(sqrt(F)) features drawn per node, seed 17 |
| Boosted trees, tree by tree | `gbt::fit`, logistic | 400 rounds, learning rate 0.1, `max_depth` 4, lambda 1, `min_child_weight` 1, `colsample` 0.5, seed 23; fitted on T with T as the validation mask and patience 400, so early stopping drops only trailing rounds that no longer lower the training loss |
| Boosted trees, end to end | `tuning_holdout(T, 0.15, seed 42)`, then `gbt::fit` | as production, with early stopping (patience 25), at `colsample` 1.0 and at `colsample` 0.5 |
| L1 logistic regression | not refitted here | the lasso suite's result at `c435000` is cited |

The harness also replays the random draws GRM makes, from the documented seeds, and writes them: each forest tree's subsample and each node's drawn features, and each boosting round's drawn features. It stops if a tree splits on a feature outside its draw. The draws are data handed to the reference; no reference library can be given GRM's random stream.

The harness records, without comparing them, the settings production would select on each panel (`scm::select` and `tree::select_depth` on `training_folds(T, 3, seed 42)`), or the reason no selection is possible.

The harness runs twice, with 1 and with 8 rayon threads. The two outputs must be byte-identical.

## References

Installed into an environment of their own, from pypi.org files named by sha256 (`pip install --no-deps --require-hashes`), on Python 3.12.3:

| Library | Version | Used for |
|---|---|---|
| scikit-learn | 1.7.2 | decision tree, forest |
| XGBoost | 3.0.5, manylinux2014 wheel (the CPU build of that release) | boosted trees |
| pyscm-ml | 1.1.2, built from its source archive with setuptools 80.9.0 and wheel 0.45.1 | set covering machine |
| numpy 1.26.4, scipy 1.15.3, joblib 1.5.2, threadpoolctl 3.6.0, six 1.17.0 | | dependencies |

A pinned implementation of the set covering machine exists (pyscm-ml), so it is used. Its rules are not GRM's (below), so a small independent implementation is written and tested as well.

### Set covering machine

Reference A, pyscm-ml 1.1.2. `SetCoveringMachineClassifier(model_type, p, max_rules=6)`. From its source:

- Its rules are stumps `x > t` and `x <= t`, with t at each value a feature takes among the remaining samples. On 0/1 features the stumps at t = 0 are GRM's presence and absence rules. The stumps at the largest value reject every sample (`x > max`) or accept every sample (`x <= max`). GRM offers such a rule only through a feature that is constant on the remaining samples. The library cannot leave these two stumps out.
- It takes utilities within 1e-7 as equal; GRM compares the float64 values.
- Among equal utilities it takes the largest covered minus lost, then the first found. GRM takes the fewest lost, then the lowest feature, presence before absence.
- It adds rules until no negative remains or `max_rules` is reached. GRM also stops when its best rule covers no negative.
- It refuses a training set with one class.

Reference B, an independent implementation in plain Python, written from the method's definition (Marchand and Shawe-Taylor 2002: greedy choice of the rule with the largest covered minus p times lost) with GRM's rule set, stop rule and tie order as `src/learn/scm.rs` documents them. It is tested on hand-worked cases before it is used.

Compared, for each of the 224 fits:

1. Reference B's whole path against GRM's: the rules in order, each step's covered, lost and utility, the stop, and the calls on every sample for every prefix of the path.
2. Reference A step by step: at each step GRM takes, the library's own utility search (`pyscm._scm_utility.find_max`) runs on the samples that remain after GRM's earlier rules. A step is **certified** when GRM's rule is among the library's best rules. It is **dominated** when every one of the library's best rules is a stump at a feature's largest value, with a utility more than 1e-7 above GRM's.
3. Reference A run freely (`fit` with `max_rules=6`): whether it gives the same rules and the same calls. Reported, not judged, because it follows from 2 and the rule differences above.

A fit **matched** when Reference B reproduces it exactly and every step is certified. A fit has a **difference with a shown cause** when Reference B reproduces it exactly, at least one step is dominated and no step is anything else; the cause is the library's two extra stumps per feature, shown by the best rules' thresholds and counts. Anything else is a **difference without a shown cause**. A fit the library refuses (one class) is judged on Reference B alone and marked so.

A panel matched when all 14 of its fits matched.

### Decision tree

Reference: scikit-learn 1.7.2 `DecisionTreeClassifier(criterion="gini", splitter="best", min_samples_split=2, min_samples_leaf=1, max_features=None, random_state=0)`.

| GRM | scikit-learn |
|---|---|
| class-balanced weights n / (2 n_class) on T | `class_weight="balanced"`, the same formula |
| impurity: weighted Gini | `criterion="gini"` |
| a node splits only on gain above 1e-12 (weighted count times impurity fall) | `min_impurity_decrease = 1e-12 / W`, W the weighted sample count of the fit; the default 0 splits on zero gain |
| `min_leaf` 1, `max_depth` 8 | `min_samples_leaf=1`, `max_depth=8` |
| equal gains: the lowest feature | features are visited in a seeded random order and the first best is kept. No setting maps this |
| node value: weighted resistant share | `tree_.value`, the same quantity |

Because scikit-learn breaks ties at random, a tied node has no single reference answer. Exact arithmetic settles it: for each node the gain of every candidate feature is computed as a fraction, giving the set of features with the largest gain.

Compared:

1. Node by node. For each node of GRM's tree, scikit-learn fits a one-split tree (`max_depth=1`) on that node's training samples with GRM's class weights as `sample_weight`. A node **agrees** when both make it a leaf, or both split and (a) scikit-learn's feature is in the exact best set, (b) GRM's feature is the lowest feature of the exact best set, and (c) GRM's two child values equal scikit-learn's within 1e-9 when the two features split the node's samples alike. The count of nodes where the two picked the same feature is reported beside the count where they picked different members of a tie.
2. Stop rules: a node at depth 8, a pure node or a node with one sample must be a leaf, and an impure node above depth 8 is a leaf only if the exact best gain is not positive.
3. Scores: GRM's scores at each depth 1 to 8, on every sample, equal the leaf values of GRM's tree walked in Python, within 1e-12.
4. The whole tree grown freely by scikit-learn (`max_depth=8`, `class_weight="balanced"`), walked beside GRM's by the training samples each node holds. Reported: the same tree, or the first node where they part and whether the two splits there have exactly equal gain.
5. The same free tree at scikit-learn's default `min_impurity_decrease=0`. Reported, not judged.

A panel **matched** when every node agrees (1), the stop rules hold (2), the scores agree (3) and the free tree (4) is the same tree or parts from GRM's only at exact ties. A **difference with a shown cause** is a node where GRM's feature is in the exact best set but is not its lowest feature while the gains GRM computes in float64 differ (rounding broke the tie order). Anything else is a **difference without a shown cause**.

### Random forest

No library can be handed GRM's random draws, so no library forest can equal GRM's tree for tree. The reference is scikit-learn's split search, applied node by node to the samples and the features GRM drew, and the standard definition of a forest score.

Compared, for each of the 200 trees of each panel (3,200 trees):

1. The subsample: round(0.632 |T|) distinct samples of T.
2. Each node, as in the decision tree's item 1, with the candidate features limited to that node's draw of ceil(sqrt(F)) distinct features, and class weights taken on the tree's subsample. A node whose draw holds no feature with positive gain must be a leaf.
3. Stop rules as in the decision tree's item 2, at depth 30.
4. The forest's score for every sample equals the mean of its trees' leaf values, computed in Python, within 1e-12.

Also reported, not judged: scikit-learn's own `RandomForestClassifier` (200 trees, `max_samples` round(0.632 |T|) with replacement, `max_features` ceil(sqrt(F)), `max_depth=30`, `class_weight="balanced_subsample"`, `min_impurity_decrease` mapped as above, `random_state` 0 and 1): the mean absolute difference of its scores from GRM's and the count of samples called differently at 0.5, beside the same two figures between the two scikit-learn seeds. It samples with replacement and draws again when no drawn feature splits, so equality is not expected.

A panel **matched** when all 200 trees pass 1 to 3 and the score passes 4. Differences are classed as for the decision tree.

### Boosted trees

Reference: XGBoost 3.0.5, `tree_method="exact"`, `objective="binary:logistic"`, one thread.

| GRM | XGBoost |
|---|---|
| learning rate 0.1, `max_depth` 4, lambda 1, `min_child_weight` 1 | `eta=0.1`, `max_depth=4`, `reg_lambda=1`, `min_child_weight=1`, `reg_alpha=0`, `gamma=0`, `max_delta_step=0`, `subsample=1` |
| class-balanced weights n / (2 n_class) on the fitting samples | the same values as `weight` |
| start at margin 0 | `base_score=0.5`, or `base_margin` given |
| split gain half of GL²/(HL+λ) + GR²/(HR+λ) − G²/(H+λ); leaf −G/(H+λ) times the rate | the same formulas, gain not halved |
| a node splits on gain above 1e-9 | on loss change above 1e-6, which is 5e-7 on GRM's scale. No setting maps this |
| equal gains: the lowest feature | the lowest feature |
| float64 throughout; Hessian floor 1e-6 | float32 gradients, margins and leaves; Hessian floor 1e-16 |
| a root below `min_child_weight` keeps its leaf value | a leaf below `min_child_weight` is set to 0 |
| `colsample` 0.5: ceil(0.5 F) features drawn per tree | its own sampler; GRM's draw is handed over as the columns of that round's matrix, in ascending feature order |
| early stopping: an improvement must exceed 1e-9, counted from the loss with no tree, and a model can keep no tree | any improvement counts, from the first round, and at least one tree is kept |

**Tree by tree.** For each tree GRM grew in the 400-round fit on T, XGBoost grows one tree from the same state: `base_margin` set to GRM's margin before that round, the same weights, the round's drawn features. Up to 6,400 trees; the count GRM returns per panel is reported.

A tree is **identical** when both trees have the same shape and the same feature at every split, and every leaf value agrees within 1e-5.

A tree that is not identical gets a cause at the first node where the two part, walking both from the root by the fitting samples each node holds. G and H are recomputed in float64 from GRM's margin. ν is the larger change in the two candidates' loss changes when the same quantities are recomputed from float32-rounded margins, gradients and Hessians.

| Cause | Shown when |
|---|---|
| equal split | the two features split the node's fitting samples alike (the same samples or the complement) |
| float32 rounding | GRM's feature has the largest float64 gain at the node, and the float64 loss change of XGBoost's feature is below it by at most 4ν |
| minimum split gain | GRM splits and XGBoost makes a leaf, and GRM's loss change is at most 1e-6 + 4ν |
| leaf set to 0 | the shapes agree, XGBoost's leaf is 0 and the node's H is below 1 |
| Hessian floor | the shapes agree and the leaf recomputed in float64 with XGBoost's floor agrees with XGBoost's within 1e-5 |
| none | anything else |

A panel **matched** when all its trees are identical. It has a **difference with a shown cause** when every tree that is not identical has a cause other than none.

**End to end.** Where `tuning_holdout` gives a fitting and a validation part with both classes, GRM's fit at `colsample` 1.0 is set against `xgboost.train` with 400 rounds, `early_stopping_rounds=25`, `eval_metric="logloss"` on the validation part with the same weights, and scores from `iteration_range=(0, best_iteration + 1)`. `colsample` 1.0 is used because XGBoost's own loop cannot be given GRM's column draws. GRM's fit at `colsample` 0.5 is recorded for its tree count only.

A panel **matched** end to end when both keep the same number of trees and every sample's probability agrees within 1e-4. Otherwise the cause is, in this order:

| Cause | Shown when |
|---|---|
| no tree kept | GRM keeps no tree, and XGBoost's validation loss in its first 25 rounds is never more than 1e-6 below ln 2, the loss with no tree |
| a tree differs | the first tree that is not identical gets a cause from the table above, computed on GRM's margin before that round; later trees follow from it |
| stopping rule | all trees both keep are identical, and GRM's rule (1e-9, patience 25, counted from ln 2) replayed on XGBoost's validation losses gives GRM's tree count |
| none | anything else |

A panel `tuning_holdout` refuses is not applicable end to end, with the error recorded.

### L1 logistic regression

Not refitted. `docs/lasso-stopping.md` describes what the lasso suite holds on these 16 panels. Its result at `c435000` is taken from `cargo test --release --locked --lib learn::lasso` and cited with its counts.

## Outcome per learner

For each of the four learners, over the 16 panels: the count that matched by the rule above, the count with a difference whose cause is shown, the count with a difference without a shown cause, and the count not applicable. A difference without a shown cause is a finding and is reported as one, with the panel, the learner and what differs. For the boosted trees the tree-by-tree and the end-to-end counts are given apart.

There is no pass mark on these counts. The lane's result is the table.

## Not compared

- The inner-fold selection of the set covering machine's kind, p and rule count and of the tree's depth, beyond recording what it selects.
- The chi-square screen: with at most 128 features it only orders them.
- The boosted trees' censored MIC objective.
- Data of full size, thread counts other than 1 and 8, and platforms other than Linux x86-64.

## Kept

The harness, the Python kit, the pinned requirements and a result table go into the repository on branch `exp/learner-crosscheck`, with a test that needs no reference library.

Anything decided after a comparison has been read is labelled exploratory in the report, and so is any change to the kit's logic made after its first run on the panels.
