# Hyperparameter selection: why tuned settings can lose to a default, what a rule must meet, and the fixed default

**What this covers.** This note covers selecting **learner hyperparameters** only: the per-family knobs in `scripts/hyperopt/protocol.py` (scm, cart, rf, gbt, lasso). It does not cover selecting the **input-data configuration**: genome and cohort count, which genomes are included, k-mer length, feature families, or the retained-feature budget. That selection is tracked separately as the primary open capability, and nothing here applies to it.

**Status.** `select.py` selects each family's **fixed default**: no tuned configuration replaces production's default (section 4, "Adopted: the fixed default"). Main decided this on 5 Nov 2026 under Layth's delegation of that day (DECISIONS F4, UD-8). The criteria A1 to A6 below were ratified as written, and the default-anchored rule was evaluated against them: keep the default unless the best inner mean beats the default's by more than γ = 0.04. It met A1 to A5 on the toy suite with the bounds as written. A6's paired outer comparison on the synthetic harness, with fresh cohort seeds, found a mean loss of 0.0099 against the default for SCM on the lineage split, above A2's 0.005 (SE 0.0074). DECISIONS F4 sets the fallback for that case: selection returns the fixed default. The one-standard-error rule with the complexity preference, which `select.py` held before, is described in sections 1 to 3. It answered finding F4 of the 2026-11-01 hyperparameter review: on one synthetic cohort it chose CART and GBT configurations that scored worse on untouched outer folds than production's defaults. `tests/hyperopt/test_selection_toy.py` encodes the criteria on generic toy data for all three rules: the fixed default, the default-anchored candidate and the one-SE baseline. It pins each value and verdict, and where a rule misses a bound, the test asserts the miss and its size. This note covers learner hyperparameters only and changes no search space, protocol or production setting.

## 1. The rule before 5 Nov 2026: one SE plus complexity

For each family and outer training fold, `select.py` scored each feasible trial *c* by the mean *m_c* and the standard error *ŝ_c* = SD(n−1)/√k of its *k* inner-fold scores (*k* = 3). The best trial *b* had the highest mean. The candidates were the trials with *m_c* ≥ *m_b* − *ŝ_b*. The rule then took the candidate with the lowest reported complexity, then the default, then the lowest trial id. Complexity is whatever the objective reports: SCM rules, CART leaves, RF trees, GBT trees kept, lasso nonzero coefficients.

Notation for section 3: trial *c*'s score on inner fold *j* is *x_cj* = *ν_c* + *f_j* + *e_cj*. Here *f_j* ~ N(0, τ²) is a fold effect that every trial shares, *e_cj* ~ N(0, σ²) is the trial's own noise, *s*² = τ² + σ², and ρ = τ²/*s*² is the shared share. *ν_c* is what the inner folds measure, and *μ_c* is the trial's true held-out score on an outer fold. They differ when the inner and outer training sets differ in size (section 3.6).

## 2. What the archived synthetic run shows

These figures were recomputed from the archived run A outputs (`selection.json`, `outer.json`, `history.jsonl`, `diagnose.json`, `summary.json`), with no refitting. All 25 selections match `select.py`'s, and every pooled balanced accuracy (BA) below equals the archived value exactly.

| family | default | one-SE selection | best inner mean (no one-SE) | static |
|---|---|---|---|---|
| CART | 0.8742 | 0.7948 | 0.8518 | 0.8962 |
| GBT | 0.8404 | 0.8285 | 0.8665 | |
| SCM | 0.8358 | 0.8895 | 0.9033 | 0.8220 |
| RF | 0.6613 | 0.7338 | 0.7327 | |

These are pooled outer BAs over 315 samples, from one cohort with one sampler seed.

- **CART.** Picking the best inner mean already loses 0.022 to the default. The preference loses another 0.057.
  - In folds 0, 1 and 4 the default was inside the window and lost on complexity: 12.3 leaves against 8.7, 9.7 against 9.3, and 11.7 against 4.0. In folds 2 and 3 it sat just under the threshold (0.690 against 0.708, and 0.813 against 0.816).
  - The fold 4 example in F4 holds. The candidates' inner means ran from 0.7295 to 0.7885, with best SE 0.069. The rule picked `{"max_depth":3,"min_leaf":1}`: 4 leaves using 3 of the 8 planted markers, outer BA 0.607 against the default's 0.797.
- **GBT. The loss has two different sources, and neither is "fewer trees than the default."**
  - **Fold 0 is selection noise.** The selected trial 5 is the highest-inner-mean trial, so the one-SE preference played no part there. Its inner mean was 0.839, and its outer BA fell to 0.754 against the default's 0.826.
  - **Folds 1 to 4 are the preference.** The rule took a lower-complexity candidate than the best-mean trial: 46 trees kept against 70, 93 against 116, 283 against 293, and 118 against 147. That cost 0.0152, 0.0525, 0.0615 and 0.0714 outer BA against the best-mean trial, which beat the default in each of those folds.
  - **The default was never a candidate.** Its inner mean was 0.06 to 0.15 below the best.
  - **The selected models were larger than the default in every fold**, both in the inner complexity proxy (68, 46, 93, 283, 118 against 43, 25, 35, 42, 33) and in the final outer refits (69, 68, 80, 324, 184 trees against 57, 37, 42, 118, 60). The inner complexity proxy, parameter interactions and final capacity are three different things (section 3.4).
- **Window size.** The window admitted 8 to 16 of 32 CART trials, 2 to 4 GBT, 1 to 4 SCM and 1 to 3 RF. "Most of the search" overstates it outside CART.
- **The inner fold scores are informative.** Across the 32 configurations, the Spearman correlation between inner mean and outer BA was 0.40 to 0.94 per fold and family.
- **On the lineage split,** the GBT selection beat its default (0.856 against 0.820), and CART still lost (0.839 against 0.857). The GBT direction is not stable across split schemes.

## 3. Why the rule can choose a worse held-out result

The toy numbers below come from section 4's suite: 32 configurations, 3 folds, *s* = 0.045, ρ = 0.3 unless stated, and 100,000 draws. The toy has no genomes, learners or real data.

### 3.1 The winner's curse in inner means

The shared fold effect shifts every trial equally, so it cannot change which trial is best. The trial's own noise can. With *K* equal trials, the best mean overstates its trial by E[*m_b* − *ν_b*] = (σ/√k)·*z_K*, where *z_K* is the expected maximum of *K* standard normals. *z_32* = 2.07, so the bias is 0.045 at σ = 0.038 and k = 3. The toy measured 0.0450 against a prediction of 0.0450.

The bias grows with *K* (0.031 at *K* = 8, 0.056 at 128) and shrinks slowly with *k* (0.035 at k = 5, 0.025 at k = 10, for *K* = 32). The one-SE pick inherits most of it: 0.0340 on the Gaussian null and 0.0405 under binomial noise on 84 samples.

In the archive, the selected trial's inner-minus-outer gap exceeded the fold's average gap over all 32 configurations by +0.117 (CART fold 0), +0.141 (CART fold 4) and +0.103 (GBT fold 0). Over five folds the excess averaged +0.056 for CART and +0.077 for GBT. GBT fold 0 is this mechanism alone.

### 3.2 One noisy SE sets the window

*ŝ_b* comes from three numbers, so it has 2 degrees of freedom. For Gaussian folds, *ŝ*/(*s*/√k) has mean 0.886 and SD 0.463. It falls below 0.5 with probability 1 − e^(−0.25) = 0.221, and exceeds 1.5 with probability e^(−2.25) = 0.105. So in about one draw in five the window is under half its typical width, and in one in ten it is over one and a half times as wide. At k = 10 those probabilities are 0.013 and 0.016. The candidate set, and so the pick, swings with this noise.

### 3.3 The window uses the wrong SE for a paired comparison

Every trial is scored on the same folds, so a difference between two trials cancels the fold effect: Var(*m_c* − *m_b*) = 2σ²/k. The window instead uses one trial's SE, whose expected square is (τ² + σ²)/k. The ratio of the window to the paired SE of a difference is 1/√(2(1 − ρ)): 0.71 at ρ = 0, 0.85 at 0.3, 1.00 at 0.5, 1.29 at 0.7, 1.58 at 0.8 and 2.24 at 0.9.

The more the folds share, the more clearly worse trials the window admits. On the toy null, the mean candidate count is 5.3 of 32 at ρ = 0, 6.7 at 0.3, 8.4 at 0.5 and 14.9 at 0.8. In the archive, a two-way decomposition put the shared share at 0.04 to 0.67 for CART and 0.05 to 0.48 for GBT. 18 of the 54 non-best CART candidates trailed the best by more than one paired SE (5 of 7, 1 of 15, 0 of 10, 8 of 8 and 4 of 14 by fold).

Overlapping training sets push the other way. With k = 3, any two inner training sets share half their samples, which correlates a trial's fold scores. As a measure of how far *m_c* is from *ν_c*, the naive SE is therefore too small. The common heuristic correction (Nadeau and Bengio, 2003) multiplies the variance by (1/k + n_test/n_train)/(1/k) = 2.5 at k = 3. The SE is calibrated neither for ranking nor for uncertainty. `select.SE_LABEL` already calls it descriptive.

### 3.4 The complexity preference picks an edge of the window, and the direction depends on the learner

Within the window the rule takes the lowest complexity. The result depends on how the true score moves with the proxy:

- If the score rises with the proxy over the window (an underfitting side), the rule picks the worst candidate.
- If it falls (overfitting), the rule picks the best one. This is the regime the one-SE heuristic assumes, as in pruning a single tree.
- If the proxy is unrelated to the score, the pick is close to a uniform draw from the window.

On the toy, the pick's cost against a uniform pick from the same candidates was +0.0084 on `underfit`, +0.0103 on `shared_folds`, +0.0039 on `clear_winner` (whose best configuration is the more complex one), +0.0007 on `unrelated_proxy`, −0.0002 on `overfit` and −0.0047 on `default_best`.

In the archive, the Spearman correlation between reported complexity and outer BA across the 32 configurations was positive in every fold: +0.46 to +0.82 for CART and +0.23 to +0.62 for GBT. On this cohort both families sat on the underfitting side of their own proxies.

For GBT the proxy is not a capacity order. Trees kept under early stopping fall as the learning rate rises (Spearman between learning rate and trees kept was −0.14 to −0.54 per fold), and depth, rounds and regularisation interact with it. So "fewest trees kept" partly means "highest learning rate", and comparing the inner proxy says nothing about the selected model's size against the default's (section 2).

### 3.5 The default is not protected

Write *μ\** for the best true score. Loss against the default is *μ_0* − *μ_sel* = (*μ\** − *μ_sel*) − (*μ\** − *μ_0*). When the default is near the best, any rule that never tests a candidate against the default turns its own oracle regret into a loss. The default wins only from inside the window and at the lowest complexity, which is why the CART default lost in folds 0, 1 and 4.

On the toy `default_best`, the loss against the default was 0.0080 for the one-SE rule and 0.0097 for picking the best mean. Dropping the one-SE preference does not fix this.

### 3.6 Inner-to-outer shift

Inner fits train on (k − 1)/k of the outer training set: 2/3 here, about 168 samples against 252. So the inner folds measure *ν_c* = *μ_c*(*n_in*), while the target is *μ_c*(*n_out*). Take a learning curve *μ_c*(*n*) ≈ *a_c* − *b_c*·*n*^(−γ). Higher-capacity configurations usually have larger *a_c* and *b_c*, so they gain more from the extra data. Inner selection is then biased toward low capacity, and the complexity preference pushes the same way.

In the archive, CART and GBT inner means sat below the same configuration's outer BA on average in every fold, by 0.02 to 0.15, and the gap varied by configuration. CART's default is an adaptive procedure: it reruns its depth search when refit on the outer training set. Its gap was larger than the fold average in 4 of 5 folds (−0.065 against −0.025, −0.120 against −0.076, −0.136 against −0.082, and −0.049 against −0.018). The default beat even the best-mean CART pick on outer folds 0, 1, 3 and 4.

On the toy `shift_capacity` scenario, where complex configurations lose up to 0.05 at the inner size, the one-SE rule lost 0.018 to the default, against 0.003 for the best mean. Inner data cannot reveal a shift like this, so section 4 reports it rather than gating on it, and A6 relies on outer folds.

### 3.7 A retrospective comparison does not validate a replacement

`diagnose.json`'s best-mean figures (CART 0.852 against 0.795, GBT 0.866 against 0.828) were computed after the fact, on the outer folds that F4 was judged on. Choosing a rule because it scored better on those folds turns them into a selection set again.

For two rules with equal expected performance, the observed difference *D* ~ N(0, *s_D*²) still gives the rule that is picked an expected advantage of E|*D*| = *s_D*·√(2/π) ≈ 0.8·*s_D*. The best-mean rule also fails A2 and A5 on the toy, still lost to the CART default, and was never scored on the lineage split. Nothing here supports adopting it.

## 4. Acceptance criteria

These hold for any learner-hyperparameter selection rule. Fix the bounds before evaluating a candidate. Main ratified them as written on 5 Nov 2026 (DECISIONS F4). They are set on the toy's scale: per-fold SD 0.045, so the SE of a 3-fold mean is 0.026. For an objective with a different per-fold SD, scale them in proportion.

**The toy suite.** Every scenario has 32 configurations, 3 folds, configuration 0 as the default, Gaussian noise with *s* = 0.045 and ρ = 0.3, and *ν* = *μ*, unless stated otherwise. "Capacity" *h* runs from 2 to 40 over configurations 1 to 31, and the default sits at 16 unless noted. Each test uses 4,000 draws from a fixed seed.

| scenario | true scores *μ* | reported complexity | what it probes |
|---|---|---|---|
| `null` | all 0.80 | capacity | winner's curse, calibration |
| `null_binomial` | all 0.80; binomial on 84 samples | capacity | the same under binomial noise |
| `default_best` | default 0.85, others 0.85 − 0.003·i | a fixed shuffle of capacity | protecting a good default |
| `underfit` | 0.88 − 0.30·e^(−h/8), default at h = 28 | capacity | a proxy that runs against the score |
| `overfit` | 0.85 − 0.15·e^(−h/3) − 0.0015·h, default at h = 30 | capacity | the regime the one-SE rule assumes |
| `unrelated_proxy` | as `underfit` | a fixed shuffle of capacity | an uninformative proxy |
| `shared_folds` | as `underfit`, ρ = 0.8 | capacity | a window inflated by the fold effect |
| `clear_winner` | default 0.80, one configuration 0.88 (h ≈ 26), others 0.80 − 0.10·i/31 | capacity | power |

**The criteria:**

- **A1. Untouched outer evaluation.** The rule is a deterministic function of the inner history, the protocol and the fold plan it pins. No outer score, outer sample or post-hoc outer comparison enters it: the fold plan carries inner fold labels only, with outer test samples null, and a record with any field outside the schema is refused. Record order does not change a flat history's selection, and a Task B history out of run order is refused (`off_plan_order`) rather than re-selected. Outer folds score each arm once, after selection, and are never used to choose or tune the rule.
- **A2. No harm against the fixed default.** On every suite scenario, E[*μ_0* − *μ_sel*] ≤ 0.005.
- **A3. Power.** On `clear_winner`, E[*μ_sel* − *μ_0*] ≥ 0.5 × 0.08 = 0.04.
- **A4. A preference must point the right way.** Among the candidates the rule treats as indistinguishable, its pick is no worse than a uniform pick from them: E[mean *μ* of candidates − *μ_sel*] ≤ 0.002 on every suite scenario.
- **A5. Calibrated reported scores.** On every suite scenario, any score reported for the selected configuration has |E[reported − *μ_sel*]| ≤ 0.01. A raw inner mean of the selected trial may be shown only labelled as a selection-biased tuning score, next to a score from data the selection did not see. F3 covers the labelling.
- **A6. Adoption evidence.** This is a process criterion, not a toy test. A new rule is adopted only if all of the following hold:
  - it passes A1 to A5 here, with the bounds fixed in advance;
  - its paired outer comparison with the fixed default, on the GRM synthetic harness, uses outer folds or cohort seeds that played no part in designing or choosing it;
  - each arm is scored once per outer fold;
  - results are reported per family and per split scheme, with the loss against the default and its uncertainty.

  A retrospective comparison on outer folds already seen (section 3.7) is a hypothesis, not validation.

### Adopted: the fixed default

For each family and outer training fold, `select.py` selects trial 0, the family default, whatever the inner scores. The history is still checked, refused or accepted as before (A1), every trial is still scored and listed, and the best feasible trial (ties to the default, then the lowest trial id) and `default_gain`, how far its mean beats the default's, are reported for the record only. A default that is not a feasible trial (failed, unconverged, off the history's folds or outside a constraint) is still the selection, and the family (`default_feasible: false`) and the run (`warnings`, code `default_not_feasible`) say so. A history with no trial 0 selects nothing (`no_default`). The report's schema is `grm-hyperopt-report/3`.

The selected entry's inner mean is labelled a selection-biased tuning score (`select.SELECTED_SCORE_LABEL`; `SCORE_LABEL` says so too, and the report heads its column "inner-tuning mean (selection-biased)"). The score to report for a selection is its outer-fold score: the configuration refitted on the outer training fold and scored once on the outer test fold. That refit is still outside this package (`select.OUTER_EVALUATION`).

**The fixed default on the toy** (the `A*` classes): A1 holds (6 tests: the signature, flat record order, Task B run order refused, an outer-score field refused, the pinned fold plan required, and every selection the default). A2's loss and A4's cost are 0 on every suite scenario. A5's outer-fold score has a bias within 0.0013 of zero on every scenario (the same outer draws as the candidate's, below). A3 is missed by design: the gain on `clear_winner` is 0 against a bound of 0.04. That is the price DECISIONS F4 set for a candidate that fails A6. `RuleCopy` checks the vectorised fixed-default rule against `select.select_history` on 175 flat histories and on 21 Task B histories written by `search.run_trial`.

### Evaluated and not adopted: the default-anchored rule

**The rule.** Among the feasible trials, keep the default unless the best inner mean *m_b* beats the default's *m_0* by more than γ; then take the best trial. The comparison is paired: every trial is scored on the same inner folds, so the shared fold effect cancels in *m_b* − *m_0* (section 3.3), and the SE plays no part. Complexity plays no part either. The candidates it treats as indistinguishable are the trials within γ of the best mean.

**γ = 0.04.** On the 100,000-draw table of the review lane `grm-hpo-f4-20261102`, the rule met A2 and A3 on every suite scenario for γ = 0.03, 0.04 and 0.05: its `default_best` losses were 0.0040, 0.0025 and 0.0013 and its `clear_winner` gains 0.075, 0.072 and 0.067; at γ = 0.02 the `default_best` loss was 0.0060. γ = 0.04 sits mid-range. At the toy's noise it is about 1.3 paired SEs of a two-trial difference (0.031), and close to the winner's-curse margin the best of 32 equal trials has over an equal default (0.045, section 3.1).

**On the toy** (the `Candidate*` classes; 4,000 draws, seed 20261102; each value pinned to within 3 Monte Carlo SEs, each verdict clearing its bound by at least six):

| criterion | passes | fails |
|---|---|---|
| A1 | yes, as the fixed default (the rule reads only the inner means) | |
| A2 (≤ 0.005) | `null` 0, `default_best` 0.00243, `underfit` −0.00130, `overfit` −0.01440, `unrelated_proxy` −0.00137, `shared_folds` −0.00044, `clear_winner` −0.07214 | |
| A3 (≥ 0.04) | gain 0.07214 | |
| A4 (≤ 0.002) | `null` 0, `default_best` −0.01418, `underfit` −0.00136, `overfit` 0.00011, `unrelated_proxy` −0.00134, `shared_folds` −0.00094, `clear_winner` −0.00511 | |
| A5 (\|bias\| ≤ 0.01; the outer-fold score) | `null` 0.00073, `null_binomial` −0.00016, `default_best` −0.00009, `underfit` −0.00014, `overfit` 0.00128, `unrelated_proxy` 0.00045, `shared_folds` 0.00097, `clear_winner` −0.00078 | |

The toy's outer-fold score is the selected configuration's true score plus independent noise with the per-fold SD, so its bias is zero up to Monte Carlo error (SE about 0.0007). The inner mean of the same selection overstates it by 0.03501 on `null`, which is why it may be shown only labelled. The rule kept the default in 42% of `null` draws, 79% of `default_best`, 50% of `underfit`, 90% of `shared_folds` and 9% of `clear_winner`.

**A6: the paired outer comparison.** The lane `grm-method-20261105` ran it on 5 Nov 2026 with a design fixed before any run (its `NOTES.md`):

- **Cohorts.** L's generator (`scripts/ablation/synthetic.py`) at seeds 900001 to 900010, which played no part in designing the rule. One override, `n_lineages` 12, gives the lineage split the 2 × 5 lineages evaluation requires. Each cohort has 400 genomes, 123 binary features (Count and Categorical columns left out) and 20 duplicate pairs.
- **Splits.** Evaluation's own split plan, seed 42, 5 outer and 3 inner folds: random, with duplicate isolation and duplicate-grouped tuning; and lineage, with lineage tuning.
- **Search and fits.** Per outer fold, Task B's protocol (32 trials per family, all five families) ran the search, scored by the real learners through a review adapter that links this tree's crate. The candidate's selection code (this lane's first F4 tree, γ = 0.04) picked the configuration. Then each arm, the fixed default and the rule's pick, was refitted once on the outer training fold and scored once on the outer test fold (16,000 trials, all ok).
- **Loss.** BA(default) − BA(selected) per outer fold; a null fit is scored 0.5, which happened once (lineage GBT, default arm). The cohort is the unit: the mean of its 5 folds, then the mean over the 10 cohorts, with SE = SD/√10 and a 95% t interval.

| scheme | family | cohorts | folds | default kept | folds differing | BA default | BA selected | loss vs default (mean of cohort means) | SE | 95% CI | <= 0.005 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lineage | cart | 10 | 50 | 40 | 9 | 0.7285 | 0.7335 | -0.0050 | 0.0040 | [-0.0141, +0.0041] | yes |
| lineage | gbt | 10 | 50 | 29 | 21 | 0.7593 | 0.7703 | -0.0110 | 0.0064 | [-0.0255, +0.0035] | yes |
| lineage | lasso | 10 | 50 | 50 | 0 | 0.7616 | 0.7616 | +0.0000 | 0.0000 | [+0.0000, +0.0000] | yes |
| lineage | rf | 10 | 50 | 35 | 15 | 0.7805 | 0.7840 | -0.0035 | 0.0028 | [-0.0097, +0.0027] | yes |
| lineage | scm | 10 | 50 | 36 | 14 | 0.7223 | 0.7124 | +0.0099 | 0.0074 | [-0.0068, +0.0266] | NO |
| random | cart | 10 | 50 | 44 | 4 | 0.7262 | 0.7277 | -0.0015 | 0.0029 | [-0.0080, +0.0050] | yes |
| random | gbt | 10 | 50 | 33 | 16 | 0.8041 | 0.8044 | -0.0004 | 0.0027 | [-0.0065, +0.0057] | yes |
| random | lasso | 10 | 50 | 50 | 0 | 0.8259 | 0.8259 | +0.0000 | 0.0000 | [+0.0000, +0.0000] | yes |
| random | rf | 10 | 50 | 43 | 6 | 0.8172 | 0.8186 | -0.0014 | 0.0023 | [-0.0066, +0.0038] | yes |
| random | scm | 10 | 50 | 40 | 5 | 0.7309 | 0.7292 | +0.0016 | 0.0012 | [-0.0011, +0.0044] | yes |

The rule lost more than 0.005 to the default in one cell, SCM on the lineage split (+0.0099). The interval there includes zero, and a few folds carry the loss: +0.252 and +0.198 BA (cohorts 900004 and 900003), then +0.059 and +0.055. The pre-stated rule reads the mean, though, and the mean is over the bound, so the rule is not adopted. Lasso kept its default in all 100 folds; its only searched knob is `n_lambdas`. Receipts: `grm-method-20261105/logs/meth-a6/grm-logs/a6/` (`summary.json`, `summary.md`, `runs.tgz` with every history and selection, `outputs.sha256`).

### The replaced rule, recorded as a baseline

**The one-SE rule on the toy** (4,000 draws; the `Baseline*` classes, each value pinned to within 3 Monte Carlo SEs and its verdict asserted):

| criterion | passes | fails |
|---|---|---|
| A1 | yes (5 tests: the signature, flat record order, Task B run order refused, an outer-score field refused, the pinned fold plan required) | |
| A2 | `null` 0, `overfit` −0.0169, `unrelated_proxy` −0.0006, `clear_winner` −0.0704 | `default_best` 0.0079, `underfit` 0.0069, `shared_folds` 0.0083 |
| A3 | gain 0.0704 | |
| A4 | `null` 0, `overfit` −0.0002, `default_best` −0.0048, `unrelated_proxy` 0.0007 | `underfit` 0.0082, `shared_folds` 0.0103, `clear_winner` 0.0039 |
| A5 | `clear_winner` 0.0073 | `null` 0.0335, `null_binomial` 0.0407, `default_best` 0.0229 |

Every verdict clears its bound by at least six Monte Carlo SEs. `RuleCopy` checks the vectorised rule the Monte Carlo uses against `select.select_history` on 175 flat histories and on 21 Task B histories written by `search.run_trial`.

The baseline classes test a vectorised copy of the replaced rule; `select.py` no longer holds it. A future rule replaces the fixed default only through the same steps: A1 to A5 here with the bounds as written, then A6 on fresh cohort seeds.

## 5. Reproduce

```bash
python3 -I -B -W error -m unittest discover -s tests/hyperopt -t . -p test_selection_toy.py -v
```

This runs the 59 toy tests in about 2 s with numpy 1.26.4. All pass, and the run ends in a plain `OK`, as CI's last-line check requires: a test for a bound a rule misses asserts the miss. A6's driver and adapter are review tools, kept with the lane's receipts (`grm-method-20261105/tools/a6/`), not in the repository. The archive recomputation and the 100,000-draw tables come from the review lane `grm-hpo-f4-20261102` (`tools/repro_f4.py`, `tools/toy_sim.py`). They read the 2026-11-01 review's saved JSON and the toy problems above only.
