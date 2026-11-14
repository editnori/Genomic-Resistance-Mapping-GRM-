# Learner cross-check: the result of 14 November 2026

GRM runs its own code for each of its five learners. On 14 November 2026 each learner was fitted on 16 small test panels, once with GRM and once with a standard version of the same method, and every difference was traced. The kit is in [`scripts/learner-crosscheck`](../../scripts/learner-crosscheck/README.md). This folder holds what it wrote.

This is a software check on small test panels. It says nothing about clinical performance and nothing about how well a learner predicts resistance. No genome data was read.

**No difference was left without a shown cause.** The differences that have one are below. One of them is a property of GRM's code that its authors may want to change: where several features have exactly equal gain, the trees do not always split on the lowest-numbered one.

## How it was made

| | |
|---|---|
| GRM | `main` at commit `c4350006b1406bf480e5885cf9501b7dd3584324`, evaluation revision 14, `src/` tree `8858b8eb434e5ca0a2ba02457ef1bcc12ed3325f`, unchanged by the kit |
| Compiler | rustc 1.98.1 (48a229cea 2026-09-01), release build |
| Panels | the 16 of `tests/data/lasso-stopping-fixtures.txt`, sha256 `7f354884a5dec09bb7e87672d56842f0c3ebd1aa5ef785e4479a83c9a1f99550` |
| References | Python 3.12.3 with scikit-learn 1.7.2, XGBoost 3.0.5 (the CPU wheel) and pyscm-ml 1.1.2, on numpy 1.26.4 and scipy 1.15.3, each installed from a file named by its sha256 (`environment.tsv`, `versions.json`) |
| Command | `scripts/learner-crosscheck/run.sh ENV OUT`, from 17:32:43 to 17:33:01 UTC on 14 November 2026 (`run.json`) |

The harness wrote the same bytes with 1 and with 8 threads. A first run at 17:28 UTC, in another environment built from the same pinned files, gave the same bytes in `summary.tsv`, `differences.tsv`, `totals.json`, `detail.json` and `versions.json`.

The harness file, `main.rs`, was reformatted with `cargo fmt` after that run: white space, trailing commas and one pair of braces. With the kit committed, the whole check ran again from the committed files, from 17:53:09 to 17:53:25 UTC, and wrote the same bytes in every file here that `run.sh` writes except `run.json`, which holds the commit and the times.

## The plan

`PLAN.md` was fixed at 17:09:18 UTC, before any learner was fitted on a panel. It names the panels, the settings, each reference, the rule by which a unit matches, and the causes a difference may be given. It sets no pass mark. `PLAN-clarifications.md` was fixed at 17:28:21 UTC, before the first run, and settles 16 places where the plan's words left a choice. Both are here as they were hashed; their sha256 values are at the end of this page. `exploratory.py` and its counts came after the first comparison was read. They are labelled exploratory and no outcome depends on them.

## Result

| Learner and check | Unit | Compared | Matched | Shown cause | No shown cause | Panels matched | Panels with a shown cause | Panels with none shown | Panels not run |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Set covering machine | fits | 224 | 194 | 30 | 0 | 12 | 4 | 0 | 0 |
| Decision tree | nodes | 400 | 395 | 5 | 0 | 12 | 4 | 0 | 0 |
| Random forest | trees | 3,200 | 2,774 | 426 | 0 | 3 | 13 | 0 | 0 |
| Boosted trees, tree by tree | trees | 5,200 | 4,766 | 434 | 0 | 7 | 9 | 0 | 0 |
| Boosted trees, end to end | panels | 15 | 11 | 4 | 0 | 11 | 4 | 0 | 1 |

A unit is what the check counts: a fit of the set covering machine, a searched node of the decision tree, a tree of the forest or of the boosted model, or a panel. A difference with a shown cause is one that the plan's rules class and that the numbers in `differences.tsv` show. A difference with no shown cause would be a finding of its own; there is none.

By panel, the units that matched of those compared:

| Panel | Samples | Features | Set covering machine | Decision tree | Random forest | Boosted trees, tree by tree | Boosted trees, end to end |
|---|---:|---:|---:|---:|---:|---:|---:|
| `campylobacter_jejuni__ciprofloxacin_lineage_fold0` | 508 | 128 | 14 of 14 | 50 of 51 | 109 of 200 | 383 of 400 | 1 of 1 |
| `campylobacter_jejuni__ciprofloxacin_random_fold0` | 508 | 128 | 14 of 14 | 65 of 67 | 133 of 200 | 385 of 400 | 0 of 1 |
| `constant_balanced` | 80 | 4 | 14 of 14 | 1 of 1 | 200 of 200 | 0 of 0 | 0 of 1 |
| `constant_imbalanced` | 80 | 4 | 14 of 14 | 1 of 1 | 200 of 200 | 0 of 0 | 0 of 1 |
| `enterobacter_cloacae__tobramycin_lineage_fold0` | 99 | 128 | 14 of 14 | 24 of 25 | 180 of 200 | 400 of 400 | 1 of 1 |
| `enterobacter_cloacae__tobramycin_random_fold0` | 99 | 128 | 14 of 14 | 27 of 27 | 178 of 200 | 396 of 400 | 1 of 1 |
| `imbalanced_5pct` | 200 | 16 | 8 of 14 | 27 of 27 | 193 of 200 | 400 of 400 | 1 of 1 |
| `marker_dup_constant` | 120 | 12 | 14 of 14 | 21 of 21 | 195 of 200 | 400 of 400 | 1 of 1 |
| `one_class_train` | 60 | 4 | 14 of 14 | 1 of 1 | 200 of 200 | 32 of 400 | not run |
| `planted_logistic` | 300 | 40 | 6 of 14 | 75 of 75 | 134 of 200 | 399 of 400 | 1 of 1 |
| `staphylococcus_haemolyticus__fusidic_acid_lineage_fold0` | 114 | 128 | 14 of 14 | 31 of 31 | 156 of 200 | 386 of 400 | 1 of 1 |
| `staphylococcus_haemolyticus__fusidic_acid_random_fold0` | 114 | 128 | 14 of 14 | 31 of 31 | 153 of 200 | 388 of 400 | 1 of 1 |
| `streptococcus_agalactiae__clindamycin_lineage_fold0` | 141 | 128 | 14 of 14 | 14 of 15 | 175 of 200 | 399 of 400 | 1 of 1 |
| `streptococcus_agalactiae__clindamycin_random_fold0` | 141 | 128 | 14 of 14 | 15 of 15 | 172 of 200 | 398 of 400 | 1 of 1 |
| `tied_counts` | 100 | 3 | 10 of 14 | 11 of 11 | 199 of 200 | 400 of 400 | 0 of 1 |
| `xor_balanced` | 80 | 2 | 2 of 14 | 1 of 1 | 197 of 200 | 0 of 0 | 1 of 1 |

`0 of 0` means GRM kept no boosted tree on the panel, so there was no tree to compare; these three panels count as matched in the table above. `not run` is `one_class_train`: its training samples have one class, and production's two-class gate on the tuning holdout refuses the fit.

## The differences and their causes

| Learner and check | Cause | Difference rows | Panels |
|---|---|---:|---:|
| Set covering machine | the library's stump at a feature's largest value, which rules out every remaining sample or none, has the higher utility | 30 | 4 |
| Decision tree | float64 rounding broke an exact tie: GRM's own gains for the two features differ | 5 | 4 |
| Random forest | float64 rounding broke an exact tie: GRM's own gains for the two features differ | 507 | 13 |
| Boosted trees, tree by tree | equal split | 64 | 7 |
| Boosted trees, tree by tree | float32 rounding | 2 | 2 |
| Boosted trees, tree by tree | leaf set to 0 | 368 | 1 |
| Boosted trees, end to end | a tree differs (equal split) | 1 | 1 |
| Boosted trees, end to end | a tree differs (float32 rounding) | 1 | 1 |
| Boosted trees, end to end | no tree kept | 2 | 2 |

**Set covering machine.** An independent implementation, written from the method's definition in exact arithmetic (`scm_reference.py`), gave GRM's rules, counts, utilities, stop and calls on all 224 fits. pyscm-ml's own search for the best rule was run at each step GRM took. GRM's rule was among the library's best at 405 of 435 steps. At the other 30, all of them a first step, the library's best rule was one GRM does not have: a threshold at a feature's largest value, which rejects every remaining sample or none of them. So 30 fits differ from what the library would fit, and the cause is the library's rule set. The 14 fits of `one_class_train` could not be read against the library, which refuses one class; they rest on the independent implementation alone.

**Decision tree and random forest.** scikit-learn searched one split on the training samples of every node GRM searched, and exact fractions gave the set of features with the largest gain. On 395 of 400 decision tree nodes and 108,575 of 109,082 forest nodes, both made a leaf or both split on a feature of the largest gain, GRM on the lowest-numbered such feature. At the other 5 and 507 nodes, GRM split on a feature of the largest gain that was not the lowest-numbered. In each case the gains of the two features are equal as fractions and differ in the last digits as GRM computes them in 64-bit floating point. GRM's code gives equal gains to the lowest-numbered feature, but it compares the rounded gains. The 507 nodes are in 426 trees. Every score GRM gave equals the score its own trees give when walked in Python.

**Boosted trees, one tree at a time.** XGBoost grew each round's tree from GRM's own margin, on the features GRM drew for the round. Of 5,200 kept trees, 4,766 were identical. In 368, all on `one_class_train`, XGBoost set a leaf to 0 because its weight was under the minimum child weight, and GRM gave the leaf its value. In 64, the two sides split a node on two features that divide its samples into the same two groups. In 2, the two candidate features differ in 64-bit arithmetic and are equal in XGBoost's 32-bit arithmetic.

**Boosted trees, end to end.** With XGBoost running by itself, early stopping included, 11 of 15 panels gave the same number of trees and probabilities within 0.0001. On `constant_balanced` and `constant_imbalanced` GRM kept no tree and XGBoost kept one, and the probabilities were equal. GRM keeps no tree when none lowers the validation loss, and counts such a model as a no-call (`src/learn/gbt.rs`). On `campylobacter_jejuni__ciprofloxacin_random_fold0` both kept 15 trees, the first tree differed by an equal split, and the probabilities differed by up to 0.0843. On `tied_counts` GRM kept 399 trees and XGBoost 400, the tree of round 1 differed by 32-bit rounding, and the probabilities differed by up to 0.000234.

## What this says about GRM

1. **Equal gains.** Among features of exactly equal gain, the decision tree, the forest and the boosted trees do not always split on the lowest-numbered one. GRM computes one side of a split directly and the other by subtraction, so two features that are equally good can get gains that differ in the last digits, and the larger rounded gain wins. The fits stay the same from run to run and from 1 to 8 threads. The choice is between features that fit the training samples equally well, but it can change a call on a sample outside training. Whether the code should compare exact counts is a decision for its authors. This check changed no code in `src/`.
2. **Set covering machine.** GRM's rule set has no rule that rejects every sample or none. Its search, tie order and stop are what the independent implementation gives.
3. **Boosted trees.** A model that keeps no tree is kept on purpose. The two other kinds of difference from XGBoost are equal splits and 32-bit rounding.

## Further counts

| Count, over the 16 panels | Value |
|---|---:|
| Set covering machine: fits the independent implementation gave rule for rule | 224 |
| Set covering machine: rules in GRM's fits | 449 |
| Set covering machine: steps read against pyscm-ml | 435 |
| Set covering machine: steps where GRM's rule was among the library's best | 405 |
| Set covering machine: of those, by the same effect only | 7 |
| Set covering machine: steps where the library preferred a rule GRM does not have | 30 |
| Set covering machine: steps of any other kind | 0 |
| Set covering machine: fits pyscm-ml refused (one class) | 14 |
| Set covering machine: fits pyscm-ml ran freely | 210 |
| Set covering machine: free fits with GRM's rules | 86 |
| Set covering machine: free fits with GRM's calls on every sample | 176 |
| Decision tree: nodes | 452 |
| Decision tree: nodes searched for a split | 400 |
| Decision tree: both made a leaf | 182 |
| Decision tree: both split on the same feature | 158 |
| Decision tree: they split on two features of one exact tie, GRM on the lowest | 55 |
| Decision tree: nodes with an exact tie for the largest gain | 89 |
| Decision tree: stop rules broken | 0 |
| Decision tree: node records wrong | 0 |
| Decision tree: panels where scikit-learn's own tree was the same tree | 7 |
| Decision tree: places where its own tree parted from GRM's at an exact tie | 18 |
| Decision tree: places where it parted anywhere else | 0 |
| Decision tree: panels with the same tree at scikit-learn's default floor of 0 | 4 |
| Random forest: nodes searched for a split | 109,082 |
| Random forest: nodes that agreed | 108,575 |
| Random forest: both made a leaf | 56,141 |
| Random forest: both split on the same feature | 45,563 |
| Random forest: they split on two features of one exact tie, GRM on the lowest | 6,871 |
| Random forest: nodes with an exact tie for the largest gain | 11,236 |
| Boosted trees, tree by tree: rounds grown | 6,400 |
| Boosted trees, tree by tree: trees kept | 5,200 |
| Boosted trees, end to end: trees GRM kept | 1,026 |
| Boosted trees, end to end: trees XGBoost kept | 1,029 |
| Boosted trees, end to end: trees GRM kept at the production column share | 1,260 |

"By the same effect only" counts steps on a feature that is present on every remaining sample. The library has no rule at 0 there, and its rule at the largest value rules out the same samples as GRM's rule.

scikit-learn's own forest cannot equal GRM's tree for tree: it samples with replacement and draws features again when none splits, and GRM does neither. `detail.json` gives, per panel, how far its scores were from GRM's beside how far two of its own seeds were from each other. Those figures are reported and not judged.

Exploratory, written after the first comparison was read:

| Exploratory count, over the 16 panels | Value |
|---|---:|
| Decision tree: split nodes | 218 |
| Decision tree: split nodes whose stored gain is the gain the kit recomputes | 218 |
| Random forest: split nodes | 52,941 |
| Random forest: split nodes whose stored gain is the gain the kit recomputes | 52,941 |
| Decision tree: nodes where rounding broke an exact tie | 5 |
| Decision tree: of those, the two features are present on opposite samples of the node | 2 |
| Random forest: nodes where rounding broke an exact tie | 507 |
| Random forest: of those, the two features are present on opposite samples of the node | 246 |
| Boosted trees: equal splits where the two features are present on the same samples | 0 |
| Boosted trees: equal splits where they are present on opposite samples | 64 |
| Boosted trees: equal splits where GRM names the higher-numbered feature | 64 |
| Boosted trees, end to end: trees both sides keep | 1,026 |
| Boosted trees, end to end: of those, identical to XGBoost's tree of the same round | 624 |

The last two rows show what a match end to end means. It is the same number of trees and probabilities within 0.0001, which is the plan's rule. It is not the same trees.

## The L1 logistic regression

It was not compared again. Its own suite ran at the same commit, `cargo test --locked --offline --release --lib learn::lasso`: 28 passed, 0 failed, 1 ignored (`lasso-suite.txt`). The ignored test, `mtb_delamanid_fold0_reaches_the_certified_optima`, needs a local folder that is not in the repository and was not on the machine, so it did not run. `docs/lasso-stopping.md` says what the suite requires.

## The check of October 2026

An earlier check of the same kind ran in October 2026 under revision 7. Its files and scripts are lost; the 16 panels are what is left of it. The site said of it that the set covering machine, the decision tree and the forest matched their reference on all 16 cases. That statement has no file behind it, and its rules are not known, so it cannot be set beside the table above. Under the rules fixed here the three learners do not match on every panel, for the causes given.

## Limits

- Sixteen small panels. A learner that agrees here can still differ on other inputs.
- The harness calls the learners with production's default settings. It does not run `grm evaluate`, so the settings production selects inside each fit are recorded in `detail.json` and not compared.
- At a node where features tie exactly, no library gives one answer. The reference there is the exact arithmetic of `exact.py`, which is part of this kit and not a library.
- The fits GRM made are not committed, because they take 12 MB. `versions.json` holds the sha256 of each panel's fits, and `run.sh` writes them again.
- `tests/learner_crosscheck` checks that these files agree with each other, with this page and with the Methods text of the site. It cannot fit a learner, so it cannot tell whether the tables are still true of a later commit.

## Files

The sha256 of each file, in the form `sha256sum -c` reads from this folder:

```
cc839abd2816fa4360408d0e0759ea71108a800600ad04092f98455d251ac75d  PLAN-clarifications.md
7c5cf3fabcd1c91d12dcf9dc80dbfeeaee8122b6aa0dabcb986dfa4a7c8909e1  PLAN.md
882bdb833536a0e092aa8b8840d6d1469fcde310893b9f75cbd3c10638023952  detail.json
19eae122edbc1c998ead379aa80002cdac8883c483c1caf5f3ea3efc9ecc4f4b  differences.tsv
1f32f76a24d3b94f77b9e0581295ed1e34ff25e4d6dc7e75ccf2432a4ff0da6f  environment.tsv
f18c073c8e12328550cf7a8292fedbaf740bbe4f6738fd408402deefd902f74b  exploratory.json
a7228320cb88696be1d7d1e7f0b6400326b05efcc971a5bb98ea089a538df402  lasso-suite.txt
4b8449284667b1c9f3d403c4ed9064b165bf20f44c8b52644c930583e400623b  run.json
3b2bc04e0d7d4c1398fa52a1f6215d6b3801482eca40357568dc45abab77b73d  summary.tsv
072f98d19274191c0f5a7312dc27791f872ec1e478ae076442d741249b7ac086  totals.json
0799875088d4286d3750a10cf7f92899e046dd9e41428ee592ea79d4435c669b  versions.json
```
