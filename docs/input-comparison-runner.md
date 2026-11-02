# Full-pipeline input comparison: the requirement and what exists

Layth specified a full-pipeline comparison over the input data. It covers five inputs: genome and cohort count, which genomes are included, k-mer length, feature families and the retained-feature budget. Each input is a named, versioned setting, with today's value as its default. Each configuration of those settings goes through the whole pipeline. Choices are made on inner folds only, every configuration is scored once on the same outer folds, and the outer predictions are compared in pairs, to see whether the choice of inputs makes a difference.

**None of this is built, integrated or run.** No code in the tree runs more than one input configuration through the pipeline and compares the results, and no difference has been measured. The components below exist and are tested on their own, but no `grm` command reaches most of them.

This page records the requirement, what each existing component takes and returns, the rules recorded for them, and what is missing. It chooses no values or ranges; those are Layth's. [input-configuration.md](input-configuration.md) traces each input as the code has it today. Line numbers are at commit `b90be1ca`.

## The five settings

| Setting | Today's value, the default | Fixed at | What records it today |
|---|---|---|---|
| Genome and cohort count | cap 2,000 per species; class minimum 25 | `grm curate --cap`, `--min-class` (`src/main.rs:35-41`) | the class minimum: `qc.json`, `datasets.json` and each report (`src/analysis.rs:746`). The cap: nothing |
| Which genomes are included | `select_genomes`, the legacy label-weighted ranking | `main.rs:611` | no version name. F's alternative records `group-round-robin-v1` (`src/curate.rs:675`) |
| k-mer length | 31 | `PIPELINE_K` (`main.rs:259`) | `counted.json`, structure, matrix and each report (`analysis.rs:743`) |
| Feature families | the K block only | `analysis.rs:793` | nothing names the family. A `BlockSpec` records a block's `version` and `source_sha256` (`src/ablate/blocks.rs:77-88`), but no production path builds one |
| Retained-feature budget | `screened_features` 4,000 per fit, with the support floors in [input-configuration.md](input-configuration.md#retained-feature-budget) | `src/learn/eval.rs:77` | each report's `eval_config` (`analysis.rs:192`, `:727`) and the input fingerprint (test `fingerprint_tracks_every_evaluation_config_field`, `analysis.rs:2255`) |

## Stages and the components that exist

The requirement names the stages a configuration goes through: curate, count, blocks and encoder per fit, the retained-feature screen, `evaluate`, then the paired comparison and report. The table gives the component that exists for each, what it takes and returns, and whether a production caller reaches it.

| Stage | Component | Takes | Returns | Production caller |
|---|---|---|---|---|
| Curate | `curate::select_genomes` (`curate.rs:629`) | phenotypes, datasets, one species, `cap` | the selected genome ids | yes: `grm curate` (`main.rs:611`) |
| Curate, opt-in | `curate::plan_group_sampling` (`curate.rs:798`) | eligible genomes with group ids, an always-included stratum, `cap`, `seed` | a `GroupSamplingPlan` (`curate.rs:724-742`) | no |
| Count | `grm count --k` (`main.rs:64-78`) | species, k from 5 to 32 | ksets, gene calls and `counted.json` | yes, one k per run; `grm pipeline` counts at 31 only |
| Fold plan | `ablate::plan::fold_plan` (`src/ablate/plan.rs:388`) | labels, `GroupIds` (lineages, duplicate groups, an optional source), scheme, k, inner_k, seed, repeat | a `FoldPlan` (`plan.rs:232-257`): outer folds, and per outer fold its inner folds and tuning holdouts, with a `sha256` | no |
| Blocks | `ablate::blocks::admit` (`blocks.rs:843`) over one of H's adapters | the adapter, sample ids | a `PredictorBlock`, or an error naming why the block is unavailable | no |
| Encode | `encode::fit_encoder` (`src/ablate/encode.rs:334`), `transform` (`:403`) | a raw block, sample ids, the outer training mask, an `EncodeConfig` | an `EncodedBlock` (`encode.rs:99-105`) | no |
| Candidates | `encode::block_candidates` (`encode.rs:528`) | an encoded block, the outer training mask | the block's columns past the support filter, within-block twins dropped | no |
| Screen | `encode::screen_block` (`encode.rs:562`) | a candidate block, labels, one fit's rows, the block's `screen_cap` | at most `screen_cap` columns, by chi-square on that fit's rows | no |
| Assemble | `encode::assemble` (`encode.rs:611`) | one fit's screened blocks, the outer training mask | one arm matrix, cross-block twins coalesced, with each column's sources | no |
| Evaluate | `learn::eval::evaluate` (`eval.rs:790`) | a feature matrix, labels, lineage groups, duplicate groups, `min_support`, MIC input, an `EvalConfig` | an `EvalOutput` (`eval.rs:308`): per split and learner, out-of-fold scores and metrics, and the final models | yes: `grm learn` (`analysis.rs:542`), on K's patterns only |
| Pair | J's `pair_arms`, `compare`, `pool`, `finalize` (`scripts/ablation/stats.py:235`, `:547`, `:627`, `:746`) | prediction rows (ablation contract §4) for two arms | contrast and pooled records | no |
| Report | K's `validate_predictions` (`scripts/ablation/schema.py:328`), `build_report` (`scripts/ablation/report.py:303`) | the run manifest, the records, the planned ids, the prediction files | `ablation-report.json` and `ablation-report.md` | no |

## Rules recorded for these components

- **Choices on inner folds only; every configuration scored once on the same outer folds.** G's `FoldPlan` fixes one (dataset, split, repeat)'s outer folds, inner folds and tuning holdouts, and every arm of it uses the same plan (`plan.rs:3-9`). Inside `evaluate`, hyperparameters are chosen in each outer training fold, and screening sees only fitting samples (`eval.rs:22-26`).
- **Label-dependent screening uses the rows of the fit it serves** (ablation contract §16). Each inner fit gets its own `screen_block` call; only the final outer fit uses the outer training mask (`encode.rs:551-561`).
- **`assemble` comes before any learner's own duplicate collapse, on the outer training mask, and a coalesced column is never split again** (contract §18; [ablation-runner.md](ablation-runner.md)).
- **One species per `plan_group_sampling` call** (`curate.rs:748-753`).
- **Pairs share samples and folds.** J pairs two arms only when they have the same samples (`stats.py:258`), the same fold plan in each repeat (`:262`), and one label, cluster and fold per sample (`:267-277`).
- **Claim levels.** K's levels are `predictive_increment`, `association_stable`, `mechanism_candidate` and `causal_supported` (`schema.py:20-21`). A contrast always has `predictive_increment` (`schema.py:23`, `:605-606`; contract §1 and §8), and the report prints a claim level only beside an `improves` decision (`report.py:81`). Features get at most `association_stable` from pipeline output; the two higher levels come only from a named manual annotation with evidence references (`schema.py:24-27`, `report.py:82-83`).
- **Decisions.** A decision is `improves`, `worse` or `inconclusive`; there is no equivalence decision (contract §15). So the most a comparison of configurations could report is a predictive increment.

## What doesn't exist

- **The caller.** No function or `grm` command takes an input configuration through curate, count, blocks, screen and `evaluate`, or passes outer predictions to J and K. `ablate::*` and `plan_group_sampling` have no caller outside their own modules, and there is no arm runner ([ablation-runner.md](ablation-runner.md)).
- **A configurable k in the pipeline.** `grm pipeline` counts at `PIPELINE_K` = 31 and recounts stages built at any other k ([input-configuration.md](input-configuration.md#k-mer-length)).
- **A budget setting.** No flag reaches `screened_features`; `grm learn` takes `EvalConfig::default()` (`main.rs:385`).
- **N and X.** 4 of 6 adapter types present; N and X unavailable (`blocks.rs:204-205`).
- **A selection-rule setting.** No `grm curate` flag chooses the rule. `run_curate` has no seed, no stratum and no group ids at curate time, and it doesn't record the cap ([input-configuration.md](input-configuration.md#which-genomes-are-included)).
- **Learners on a given fold plan.** `evaluate` draws its own outer folds (`outer_splits`, `eval.rs:941-952`) and takes no `FoldPlan`. G's plan gives the same outer folds only at repeat 0 with `EvalConfig::default()` (`plan.rs:370-375`). `fit_fold`, which fits one learner on one outer fold, is private (`eval.rs:608`) and draws its own inner folds (`eval.rs:619`).
- **Prediction rows from `grm`.** `EvalOutput` gives each learner's out-of-fold score per sample, NaN for a sample without a call (`eval.rs:196-199`), and each sample's outer fold (`eval.rs:223-226`). Contract §4 rows also need `arm`, `cluster`, `called` and `fold_plan_sha256`, and a no-call row scores the training prevalence. Nothing converts one into the other.
- **A pairing across cohorts or fold plans.** J refuses two arms with different samples or different fold plans (`stats.py:258`, `:262`). A cohort setting (the cap, the class minimum or the selection rule) changes the samples. A change of k can change the lineage labels and duplicate groups that the fold plan reads ([input-configuration.md](input-configuration.md#k-mer-length)). J and K have no rule for either comparison.
- **A link between the Python packages and `grm`.** `scripts/hyperopt` and `scripts/ablation` call no `grm` binary, and no `grm` command calls them. `scripts/hyperopt` searches learner settings only (`scripts/hyperopt/protocol.py:141-196`).

## Status

No result exists, and none is claimed. Running a comparison on real cohorts is held. The real-cohort k comparison was refused and isn't rerouted (`docs/work-program-20261027.md:148`, `:339`). The values and ranges of each setting are Layth's to choose.
