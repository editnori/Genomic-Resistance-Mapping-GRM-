# What an ablation arm runner must do

No code in the tree runs ablation arms yet. The parts exist: G's fold plans (`src/ablate/plan.rs`), H's blocks (`blocks.rs`; 4 of 6 blocks available: N and X have adapter types but no source), I's encoders (`encode.rs`) and F's group-sampling plan (`src/curate.rs`), each tested on its own. The order they run in and the samples each one sees are the caller's job, and a mistake in either gives a wrong result with no error. This page lists what the runner owes, and what the curation code that calls `plan_group_sampling` owes.

Section numbers are those of the GRM ablation data contract v1.14 (main orchestration ruling, 2026-10-31, sha256 `2300fbdbf920f847518a3f1628843826360f709b170c551f38db50d1f13deb76`). Where this page or the code disagrees with the contract, the contract wins.

## Each outer fold, in this order

G's `fold_plan` gives one `FoldPlan` per dataset, split and repeat. Every arm of that dataset, split and repeat uses it, so arms differ only in their blocks. For outer fold `f`:

1. **Masks (§16).** Take every mask from the plan. The outer training samples have `outer[i] != f`. An inner CV training fold is the outer training samples outside that inner fold of `inner[f]`. A tuning fit (the `fit_m` part of `tuning_holdout`) is the outer training samples outside `holdouts[f].gbt` for boosting, or `holdouts[f].lasso` for the lasso. Draw no other folds or holdouts.
2. **Encode (§16, M-2).** Call `fit_encoder` on the outer training samples, then `transform`. A block that can't be encoded, such as one with no row for any training sample, returns an error with its reason. Record the block as unavailable with that reason, and never drop it from the arm silently.
3. **`block_candidates` on the outer training mask (§16).** The support filter and within-block collapse read no labels.
4. **`screen_block` once per fit, on that fit's rows (§16, M-3).** Screening reads labels. Each inner CV training fold, each tuning fit and the final outer fit get their own call with their own mask from step 1. A block screened on the outer training mask must never reach inner tuning, or validation labels choose its columns.
5. **`assemble` each fit's screened blocks, then the learners (§18).** Call `assemble` before any learner's own duplicate collapse (`training_candidates`), which would put the representative's placeholder 0 back over a lower-priority twin's observed value. Pass `assemble` the outer training mask. That is the mask the base learners' `training_candidates` runs on (`src/learn/eval.rs:831`), so `assemble`'s collapse matches the learners' own. I's tests pass the same mask.

So the order on each fit is `block_candidates` (once per outer fold, on the outer training mask), then `screen_block` (on that fit's rows), then `assemble` (on the outer training mask), then the learners.

## After `assemble`

**Never re-split a coalesced column (§18).** A collapsed column takes, on each sample, the value of the highest-priority source that observes it. Splitting it back into per-block columns, for example to report each block's importance, puts the representative's placeholder back on every sample that only a lower-priority twin observes. That is what §18 exists to prevent. To attribute a collapsed column to blocks, read its provenance instead: every source in priority order, and how many samples each supplied (`ColumnProv::sources`, `ColumnSource::supplied`).

**Within-block twins don't coalesce (§19).** `block_candidates` keeps the first of a block's training-identical columns, so a dropped twin's observation is lost where the kept column isn't assessable. The loss is the same in every arm that contains the block, so it can't bias a leave-one-block-out contrast. Coalescing within a block is the candidate refinement if measured blocks show the loss matters. It would break the K-only equivalence with the base pipeline's `training_candidates`, so it needs a contract change first, and a runner must not do it on its own.

## Inputs and outputs

No runner function exists, so it has no signature yet. The data at its edges is set by types that exist. Line numbers are at commit `b90be1ca`.

**In:**
- one `FoldPlan` per (dataset, split, repeat), from `plan::fold_plan` (`src/ablate/plan.rs:388`), that passes `FoldPlan::verify` (`plan.rs:360`). A plan that can't be scored (`is_available`, `plan.rs:262`) has its `unavailable_reason`, and K's run manifest records such a run as `unavailable` with that reason ([ablation-report.md](ablation-report.md#run-manifest)).
- the arm's blocks, each admitted with `blocks::admit` (`src/ablate/blocks.rs:843`) over the dataset's sample ids. 4 of 6 blocks are available: N and X have adapter types but no source, and `admit` refuses them with their reasons (`blocks.rs:204-205`, `:848-850`), and an `ok` run uses no unavailable block.
- the labels, in the plan's sample order, and an `EvalConfig` (`src/learn/eval.rs:54-66`) for the learners. G's plan for the pipeline's own folds takes its `k`, `inner_k` and `seed` from `EvalConfig`'s `folds`, `inner_folds` and `seed` (`plan.rs:237-242`, `:370-375`). The encoder takes an `EncodeConfig` (`src/ablate/encode.rs:48-53`).

**Out:** contract §4 prediction rows, one per (dataset, arm, learner, split, repeat, sample_id), with J's field names (`scripts/ablation/stats.py:150-151`) and K's `PREDICTION` record (`scripts/ablation/schema.py:310-314`). A row's `fold` is the sample's outer test fold in the plan, and its `fold_plan_sha256` is the plan's `sha256`, the same across arms and learners. A no-call row has `called` false and scores the training prevalence.

The learner side doesn't fit this yet. `learn::eval::evaluate` (`eval.rs:790`) takes no `FoldPlan`: it draws its own outer folds (`eval.rs:941-952`) and returns an `EvalOutput` (`eval.rs:308`) with no §4 rows. `fit_fold` (`eval.rs:608`), which fits one learner on one outer fold, is private.

## Group sampling: one species per call

Call `curate::plan_group_sampling` once per species, with only that species' eligible genomes. Group ids aren't tied to a species, and two species' metadata can both use `g1`. In one call, those genomes form one group. That group gets one round-robin draw per round instead of one per species, and every `probability` in it is computed for the mixed group. Nothing fails, because `EligibleGenome` carries no species. A mixed call is wrong even when the species' group ids are disjoint: one call has one `cap` and one always-included stratum, so the species share them instead of each having its own. The curation code that calls the planner should check that every genome in a call is of one species.

The docs of `screen_block`, `block_candidates`, `assemble` and `plan_group_sampling` state their parts of this, and the last two point here.
