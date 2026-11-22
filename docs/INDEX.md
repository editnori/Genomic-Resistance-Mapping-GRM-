# Documentation index

Every document under `docs/`, grouped by purpose. For the headline Revision 14 numbers and links to the live dashboard, start with [README.md](../README.md).

## Findings and evaluations

| Document | Summary |
| --- | --- |
| [release-summary.md](release-summary.md) | Main Revision 12 and 14 figures: per-learner medians with species bootstrap intervals, random-to-lineage drops, error-screen counts, FDA 2009 guidance wording, and known-gene rule bands. |
| [study-split.md](study-split.md) | Retrospective split by study (BioProject) inside BV-BRC across 70 datasets in 9 species, compared with random and matched random folds. |
| [time-split.md](time-split.md) | Retrospective first-public-day split inside BV-BRC across 24 datasets in 4 species, compared with size-matched random draws and study splits. |
| [outside-study-results.md](outside-study-results.md) | Retrospective test of frozen models on newly deposited BV-BRC studies (7,258 primary genomes and 3,155 same-study genomes), including Revision 14 reading. |
| [outside-study-results-by-species.md](outside-study-results-by-species.md) | Species-by-species breakdown of the outside-study test. |
| [outside-study-test.md](outside-study-test.md) | Frozen evaluation plan for the outside-study test. |
| [other-database-test.md](other-database-test.md) | Released Revision 14 models scored on 13,377 NCBI Pathogen Detection genomes outside BV-BRC's labelled set (9 species). |
| [kover-comparison.md](kover-comparison.md) | Kover 2.0.7 compared with GRM's set covering machine on the 12 smallest lineage-split datasets (1,118 genome-dataset entries). |
| [dense-network.md](dense-network.md) | One-setting dense neural network on 4,000 screened features compared with the five classical learners across 24 datasets in 12 species. |
| [amrfinderplus-comparison.md](amrfinderplus-comparison.md) | Full AMRFinderPlus 4.2.7 program compared with GRM's known-gene rule and five learners across 235 datasets in 21 species. |
| [tbprofiler-comparison.md](tbprofiler-comparison.md) | TBProfiler 6.6.3 run on 2,000 *Mycobacterium tuberculosis* assemblies across 16 named endpoints (14 graded) against GRM's learners and known-gene rule. |
| [baselines.md](baselines.md) | No-genome baselines (`grm baselines`): constant training-fold MIC and lineage-identity baseline. |
| [clindamycin-labels.md](clindamycin-labels.md) | Single-study label concentration and sensitivity refits across the six clindamycin datasets. |
| [input-configuration.md](input-configuration.md) | Default input settings, feature families, label rules, and the planned input comparison. |
| [input-comparison-runner.md](input-comparison-runner.md) | Runner, pairing rules, and second reading for the 4-species input comparison (H-1). |
| [tie-rule.md](tie-rule.md) | Decision tree tie-breaking evaluation record. |
| [hyperopt-selection.md](hyperopt-selection.md) | Hyperparameter tuning evaluation and why fixed defaults remain in place. |
| [lasso-stopping.md](lasso-stopping.md) | Lasso convergence rule and the 18-dataset Revision 14 scan. |
| [ablation-report.md](ablation-report.md) | Schema and validator contract (`ablation-report` 1) for feature-block ablations. |
| [ablation-runner.md](ablation-runner.md) | Runner contract for feature-block ablations. |
| [association-plan.md](association-plan.md) | Planned association analysis specification. |
| [assumption-plan.md](assumption-plan.md) | Fixed plan, written on 21 November 2026 before any fit, for testing each assumption a GRM result depends on: one diagnostic, nine comparisons (cutoff, cap, inputs, screen, settings, labels and MIC, outside learners, stability, pick), and a later-records test. |
| [assumption-plan-later-records.md](assumption-plan-later-records.md) | Counts behind the plan's later-records test: 252,148 BV-BRC lab records set aside on 21 November 2026, counted with no fit. |
| [assumption-plan-corrections.md](assumption-plan-corrections.md) | Dated corrections to the fixed plan: the study-split row, the release's missing fold models, and the design of the MIC comparison. |

## Methodology, pipeline, and contracts

| Document | Summary |
| --- | --- |
| [duplicate-isolation.md](duplicate-isolation.md) | Keeping k-mer-identical duplicate genomes together in outer, inner, and final tuning folds (revisions 10 to 12). |
| [equal-scoring.md](equal-scoring.md) | Scoring null fits at training prevalence so every learner covers every isolate (revision 14). |
| [species-check.md](species-check.md) | Species-identity check that drops genomes whose sequence does not fit their labelled species (revision 10). |
| [fold-draws.md](fold-draws.md) | Repeated fold-draw design and stability checks. |
| [fold-file.md](fold-file.md) | Supplied fold-plan file format for `grm learn --fold-file`. |
| [mic-plan.md](mic-plan.md) | Experimental MIC models on a fold plan (`grm learn --exp-mic`, off by default), and the scripts that run and score one species. |
| [null-fold-refit.md](null-fold-refit.md) | Experimental null-fold refit flag (off by default). |
| [matrix-scaling.md](matrix-scaling.md) | Out-of-core partial matrix building (`--disk-sparse`) and the `REVIEWED_COLUMNS` (4,096) review gate. |
| [lasso-layout.md](lasso-layout.md) | Machine-code layout check for `grm::learn::lasso::column_sum`. |
| [prediction-contract.md](prediction-contract.md) | `grm predict` input/output contract. |
| [dataset-export.md](dataset-export.md) | `grm export` binary matrix and metadata contract for external learners. |
| [provenance-fields.md](provenance-fields.md) | `grm fetch genome-dates` and `grm provenance` sidecar fields. |
| [report-rule-fields.md](report-rule-fields.md) | Optional `fold_models` and `rule_groups` report fields. |
| [rule-words.md](rule-words.md) | Plain-language wording rules for Set Covering Machine rules on the site. |
| [synthetic-assurance.md](synthetic-assurance.md) | Extended synthetic fold-mask and k-mer assurance suites. |
| [assurance-open-proofs.md](assurance-open-proofs.md) | Open equivalence targets in the Bend2 fold-mask proofs. |
| [givens-export.md](givens-export.md) | Export contract for `reference/givens-export/grm.json`. |
| [givens-terms.md](givens-terms.md) | Generated readable summary of GRM's Givens definitions, assumptions, and counts. |

## Operations and hosting

| Document | Summary |
| --- | --- |
| [operations-reference.md](operations-reference.md) | Detailed pipeline reuse checks, `grm curate` reruns, `grm intervals`, `grm baselines`, `grm provenance`, deployment test skip rules, and result checker options moved out of `README.md`. |
| [static-hosting.md](static-hosting.md) | Current `genomicresistance.com` hosting on Cloudflare Worker + R2 (`ops/cloudflare/`). |
| [deployment.md](deployment.md) | Host-based `grm serve` + Cloudflare Tunnel deployment toolkit (`scripts/deployment/`), retained for rollback/host use and dry-run publication checks. |
| [ci.md](ci.md) | Continuous integration jobs, exact test counts, pinned toolchains, and checks outside CI. |
| [data-collection.md](data-collection.md) | Archiving and verifying project files with `scripts/collect/collect.py`. |

## Historical snapshots

| Document | Summary |
| --- | --- |
| [revision-history.md](revision-history.md) | Revision 5 through 14 release history, retired Vast host notes, revision-5 `research.json` packaging details, and legacy archive notes moved verbatim from `README.md`. |
| [research/genotype-to-phenotype.md](research/genotype-to-phenotype.md) | **Historical snapshot** (October 2026, evaluation revision 7). Literature review and early diagnostic/surveillance assessment; see [release-summary.md](release-summary.md) and the newer evaluation pages above for Revision 14 figures. |
