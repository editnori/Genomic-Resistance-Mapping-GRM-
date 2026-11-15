# Lane crosscheck: clarifications to the plan, fixed before any saved panel was fitted

`PLAN.md` (sha256 `7c5cf3fabcd1c91d12dcf9dc80dbfeeaee8122b6aa0dabcb986dfa4a7c8909e1`) is unchanged. Writing the kit showed places where the plan's words leave a choice. Each choice below was made before the harness or the kit ran on a saved panel. Until then the kit ran only on six made-up panels (`toy/toy-fixtures.txt` in the scratch), and on copies of their harness output with one value changed, to see that the kit reports a changed fit. Nothing from those made-up panels is evidence about the saved panels.

## Set covering machine

1. A step is certified when the library's best rules hold GRM's rule: the stump at 0 on the same feature, `x > 0` for presence and `x <= 0` for absence. On a feature that is 1 on every remaining sample the library has no stump at 0. Its stump `x > 1` on that feature then rules out the same remaining samples as GRM's absence rule, and the step is certified too. The count of steps certified only in this second way is reported.
2. On a training set with one class the library's fit is refused, and its utility search is not run either. Those fits are judged on the independent implementation alone, as the plan says.
3. The independent implementation compares utilities as exact fractions of the float p. A step's utility agrees when GRM's value is the float64 expression `covered - p * lost` and is within 1e-12 of the exact one.
4. The stop is judged on the independent implementation. What the library's own fit does after GRM stops is part of the free fit, reported and not judged.

## Decision tree and forest

5. scikit-learn's one-split tree runs on every node above the depth limit, pure nodes included. A node at the depth limit is a leaf by the limit and no split search, so it is covered by the stop rules and by the free tree. In the forest a pure node drew nothing, and scikit-learn is given every feature there.
6. A node's record is also checked: its sample count, its resistant count and its value against the exact class-weighted share, within 1e-12.
7. The subsample size is round(0.632 |T|) with halves rounded away from zero, as Rust rounds.
8. Units: the decision tree is counted in searched nodes, the forest in trees. A tree of the forest matched when its subsample, every node and every stop rule pass.

## Boosted trees

9. The tree-by-tree fit offers the features in the order `screen` gives on T, as production does. The features a round drew are handed to XGBoost in ascending order.
10. Two trees are walked breadth first, the absent child before the present one. A difference of shape or feature is classed at the first such node. When shapes and features agree, every leaf that differs by more than 1e-5 needs a cause of its own, and the tree has a shown cause only if each has one.
11. ν is computed at the node from the two candidates' loss changes, each once from float64 quantities and once from float32-rounded margin, gradients and Hessians with XGBoost's Hessian floor. "GRM's feature has the largest float64 gain" allows a relative slack of 1e-9 for the kit's own summation order, which is not GRM's.
12. "Leaf set to 0" reads the node's Hessian sum as XGBoost holds it, from float32 gradients.
13. End to end, the cause "stopping rule" applies only when the two keep different numbers of trees. Equal counts, identical trees and a probability gap above 1e-4 have no shown cause.
14. The end-to-end fit at the production column share uses seed 42 for its draws, the seed production passes to both the holdout and the booster. Only its tree count is recorded.
15. If the kit's own margin, summed from GRM's trees, differs from GRM's by more than 1e-9, the kit stops: it would not be reading GRM's trees as GRM does.

## Outcome per panel

16. A panel's outcome for a check is matched when no difference row is written for it, explained when every row has a cause, and unexplained when a row has none.
