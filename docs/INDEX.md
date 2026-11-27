# Documentation index

Every document under `docs/`, grouped by purpose. For the headline Revision 14 numbers and links to the live dashboard, start with [README.md](../README.md).

## Release summary and validation splits

| Document | Summary |
| --- | --- |
| [release-summary.md](release-summary.md) | Main Revision 12 and 14 figures: per-learner medians with species bootstrap intervals, random-to-lineage drops, error-screen counts, FDA 2009 guidance wording, known-gene rule bands, and the path from lab records to datasets (revision 14). |
| [study-split.md](study-split.md) | Retrospective split by study (BioProject) inside BV-BRC across 70 datasets in 9 species, compared with random and matched random folds. |
| [time-split.md](time-split.md) | Retrospective first-public-day split inside BV-BRC across 24 datasets in 4 species, compared with size-matched random draws and study splits. |
| [time-split-stagger.md](time-split-stagger.md) | Frozen plan for a staggered split by time inside BV-BRC, on the public day and the collection year: cut points, windows, contrasts, rule and predictions. Not run yet. |
| [outside-study-results.md](outside-study-results.md) | Retrospective test of frozen models on newly deposited BV-BRC studies (7,258 primary genomes and 3,155 same-study genomes), including Revision 14 reading. |
| [outside-study-results-by-species.md](outside-study-results-by-species.md) | Species-by-species breakdown of the outside-study test. |
| [outside-study-test.md](outside-study-test.md) | Frozen evaluation plan for the outside-study test. |
| [other-database-test.md](other-database-test.md) | Released Revision 14 models scored on 13,377 NCBI Pathogen Detection genomes outside BV-BRC's labelled set (9 species). |
| [clindamycin-labels.md](clindamycin-labels.md) | Single-study label concentration and sensitivity refits across the six clindamycin datasets. |

## Assumption plan and block results

| Document | Summary |
| --- | --- |
| [assumption-plan.md](assumption-plan.md) | Fixed plan, written on 21 November 2026 before any fit, for testing each assumption a GRM result depends on: one diagnostic, nine comparisons (cutoff, cap, inputs, screen, settings, labels and MIC, outside learners, stability, pick), and a later-records test. |
| [assumption-plan-later-records.md](assumption-plan-later-records.md) | Counts behind the plan's later-records test: 252,148 BV-BRC lab records set aside on 21 November 2026, counted with no fit. |
| [assumption-plan-corrections.md](assumption-plan-corrections.md) | Dated corrections to the fixed plan: the study-split row, the release's missing fold models, the design of the MIC comparison, four changes to its checks, the readings of D0's causes, and the lead's rulings on the label rule, condition 6, the rule for the one release change and the run on the datasets of M; where the site's label text now sits; and which of the plan's counts are the rebuilt table's. |
| [assumption-refit.md](assumption-refit.md) | Refit of the release's rule with fold models saved, its check against the release, and the plan's D0 and C8 readers run on it, one species a box. |
| [refit-result.md](refit-result.md) | Result of the plan's D0 and C8 on the refit of 238 datasets in 23 species, read once by rules fixed before any count: Absent holds 107 of 177 miss pairs, so C3 runs first; the median overlap of chosen features across lineage folds is 0.050 to 0.183 per learner. |
| [cohort-draw.md](cohort-draw.md) | Experimental small cohorts drawn from large datasets (`grm learn --exp-cohort`, off by default), and the scripts that draw, run and score one species. |
| [cohort-result.md](cohort-result.md) | Small-cohort result of plan block C1 part B: the estimate's error at 10, 15 and 20 genomes a class against 25, on 89 datasets in 10 species, read once by rules fixed before any score. No lower line is a candidate. |
| [cap-arms.md](cap-arms.md) | Cap arms beside the release's 2,000 genomes a species for plan block C2: the top 4,096 by the release's selection rule (the top 4,000 withdrawn before any fit), scored on the default's held-out genomes, how added genomes are kept out of held-out lineages and BioProjects (fetched by id from BV-BRC's genome table), the datasets where no added genome is labelled, and the scripts that list, fit and score one species, up to four datasets at once. |
| [cap-arms-result.md](cap-arms-result.md) | Plan block C2 read once on the first run's 143 datasets in 8 species: the top 4,096 genomes a species by the release's selection rule against the release's 2,000, on the default's held-out genomes, with the twin; no learner carries conditions 1 and 2, so the rule reads no difference shown, and the release keeps its cap of 2,000 until the maintainer rules; every species and antibiotic named. |
| [input-arms.md](input-arms.md) | Experimental input arms beside 31-mer presence (`grm learn --exp-inputs`, off by default) for plan block C3: the AMRFinderPlus program's calls alone and joined to K31, the designs of the other three arms, and the scripts that list, run and score one species. |
| [input-arms-m-count.md](input-arms-m-count.md) | The count over the datasets of M: each arm against the conditions of ruling 15, the counts per pair and learner, the size table, the conditions table and the floor; written whole by `scripts/c3-inputs/m_page.py page` from `reference/input-arms-m/`. |
| [input-arms-m-left-out.md](input-arms-m-left-out.md) | The datasets of M left out of the run's second reading, each with its cohort, listed and left counts and the floor; written whole by `scripts/c3-inputs/m_page.py left` from `reference/input-arms-m/`. |
| [input-arms-result.md](input-arms-result.md) | Plan block C3 read once on the first run's 233 datasets in 23 species: five inputs against 31-mer presence refitted in the same run, with the twin; by ruling 15 the AMRFinderPlus program's calls alone meet conditions 1 to 5 on 203 of A's 204 datasets in 21 species, the other four arms do not, and the release keeps its input until the maintainer rules. |
| [screen-arms.md](screen-arms.md) | Experimental screen and floor arms (`grm learn --exp-screen cmh`, `--exp-floor`, off by default) for plan block C4: a Cochran-Mantel-Haenszel screen over the release's lineage clusters at budgets of 4,000 and 16,000, a candidate floor fixed at 3, the scripts that list, run and score one species, the unflagged-bytes check and the negative control. |
| [screen-arms-result.md](screen-arms-result.md) | Plan block C4 read once on the first run's 233 datasets in 23 species: the Cochran-Mantel-Haenszel screen at budgets of 4,000 and 16,000 and the candidate floor fixed at 3 against the release's chi-square screen, condition by condition, with the twin, the counts beside and every species and antibiotic; the release's screen, budget and floor do not change. |
| [settings-arms.md](settings-arms.md) | Experimental settings arms (`grm learn --exp-settings` and `--exp-settings-folds`, off by default) for plan block C5: a nested search of 30 random settings per learner inside each outer fold's training genomes, a call cutoff chosen on inner scores, the driver and the scripts that list, run and score one species. |
| [settings-arms-result.md](settings-arms-result.md) | Plan block C5 read once on the first run's 55 datasets in 8 species: a search of each learner's settings inside each outer training fold and a call cutoff chosen on inner scores, each against the release's rule refitted in the same run, with the twin and the plain best beside; no learner carries conditions 1 and 2 for either arm, so the rule reads no difference shown, and the release keeps its settings and its 0.5 cutoff until the maintainer rules; every species and antibiotic named. |
| [mic-plan.md](mic-plan.md) | Experimental MIC models on a fold plan (`grm learn --exp-mic`, off by default), and the scripts that run and score one species. |
| [mic-comparison.md](mic-comparison.md) | MIC comparison of plan block C6 part B: the released MIC model, an ordered model and the training median on lineage and study folds, 100 and 28 datasets in 12 species, read by rules fixed before any score. |
| [outside-learners.md](outside-learners.md) | Outside learners for plan block C7: XGBoost 3.0.5 and LightGBM 4.6.0 on the 4,000 rows `grm export` screens, read against GRM's boosted trees, and a scikit-learn 1.7.2 elastic net on unitigs, read against GRM's lasso; the pins, the folds, the scorer, the datasets, the pilot on Klebsiella aerogenes and the planted faults. |
| [outside-learners-result.md](outside-learners-result.md) | Plan block C7 read once on the first run's 233 datasets in 23 species: XGBoost 3.0.5 and LightGBM 4.6.0 against GRM's boosted trees and a scikit-learn 1.7.2 elastic net against GRM's lasso, condition by condition, with the twin, the counts beside and every species and antibiotic; LightGBM and the elastic net meet conditions 1 to 5 and XGBoost meets 1 to 4, and the release keeps its learners until the maintainer rules. |
| [served-pick.md](served-pick.md) | Experimental pick of the served learner inside each outer training fold (`grm learn --exp-pick`, off by default) for plan block C9, and the scripts that list, run and score one species. |
| [served-pick-result.md](served-pick-result.md) | Plan block C9 read once on the first run's 233 datasets in 23 species: the learner picked inside each outer training fold against the release's pick, condition by condition, with the twin, the three answers beside and every species and antibiotic; the release's pick does not change. |
| [mechanism-split.md](mechanism-split.md) | Plan blocks M1 and M2 on the AMRFinderPlus program's tables, fitting nothing: each learner's recall on resistant genomes whose catalogue calls are point mutations, acquired genes, both or neither, and the mix of catalogue calls over collection years with one vote per BioProject; M3 not built; M4, the lineage curve for Mycobacterium tuberculosis, read once in `lineage-curve-result.md`; the scripts that run one species and combine them. |
| [mechanism-split-result.md](mechanism-split-result.md) | Plan blocks M1 and M2 of the mechanism split read once on the 237 datasets of 21 species by the method page's fixed rules: on the 19 datasets with 10 or more resistant genomes of each kind, each learner's recall on resistant genomes with only an acquired gene against those with only a point mutation reads no difference shown, with GRM's known-gene rule and the AMRFinderPlus program beside; 110 of 285 species and drug class cells pass the BioProject vote over collection years and none moves; every species and antibiotic named. |
| [lineage-curve-result.md](lineage-curve-result.md) | Plan block M4 of the mechanism split, the lineage curve for Mycobacterium tuberculosis, read once by the method page's fixed rules over TBProfiler 6.6.3 with `tbdb` commit 3ed1cb99 run again on the 2,000 genomes of the 22 revision 14 reports: per antibiotic, the distinct graded variants of the WHO catalogue among resistant genomes as lineages are added, the resistant genomes with no graded variant, and the variants seen in one lineage only. |
| [learner-overlap-result.md](learner-overlap-result.md) | Learner overlap on the lineage split of the 274 revision 14 reports, read by a plan fixed before the script and fitting nothing, with a declared second pass: on the 268 reports with five learners, the resistant genomes each number of the five misses at the release's 0.5 cutoff, beside the count expected if they missed independently. All five miss 6,316 of 93,786 resistant genomes where 656.85 are expected, a median ratio of 488.8 over 259 reports (535.0 over 253 with null-fit folds left out), and 23,859 more are missed by 1 to 4 of the five. The five-learner mean misses 14,383 of the 93,786. |
| [catalogue-size-result.md](catalogue-size-result.md) | Whether the number of graded variants of the WHO catalogue among the resistant genomes of each antibiotic of Mycobacterium tuberculosis, from `lineage-curve-result.md`, tracks the resistant share of the antibiotic's revision 14 report or only how many resistant genomes there were to study: Spearman rank correlations over 16 antibiotics, read once by a rule committed before the read; the share's correlation is the larger by less than the rule's margin of 0.15, so these numbers cannot separate the two. |
| [mechanism-year-result.md](mechanism-year-result.md) | The mix by collection year over the primary set's 204 datasets, read once from the 274 revision 14 reports alone, fitting nothing: M1's four strata were not readable from the reports, so it gives, year by year, the share of resistant genomes GRM's known-gene rule makes no call on and the five learners' recall within each group, with the plan fixed before the read, and relates that no call to M1's strata and to the AMRFinderPlus program's own call from M1's committed per-dataset table. |
| [combiner-result.md](combiner-result.md) | Plan, gates and planted faults for FREEZE-COMBINE, committed before any reading: whether a weighted mean of the five learners' stored held-out scores, or a logistic combiner on them, fitted on the other folds of the release's lineage split, calls resistant genomes their mean misses on the 268 reports the learner overlap pooled. Not read yet. |
| [matched-accuracy.md](matched-accuracy.md) | The matched accuracy protocol: stage 1, the script that counts eligible genome rows and freezes them, opening no score, which ran once; and stage 2's script, which compares the five learners' scores inside each matched group; stage 2's real reading is not authorised. Not run yet. |
| [four-track-plan.md](four-track-plan.md) | Frozen pre-run specification for four parallel post-Revision-14 evaluation tracks (Track A published AMR tool benchmark, Track B point-mutation representations, Track C cross-database deduplication CLI, Track D nested selective calling and abstention), keeping Revision 14 artifacts unchanged. Not run yet. |
| [point-mutation-plan.md](point-mutation-plan.md) | Track B protocol for comparing 31-mer baselines against protein-aligned known-target variants, gene and active-site burden, and audited foundation-model representations (`esm2`, `evo2`, `prokbert`, `bacformer`; AlphaGenome 2 excluded) under a fixed downstream learner on `lineage` and `study-whole` folds. Not run yet. |

## External comparisons, baselines, and earlier studies

| Document | Summary |
| --- | --- |
| [kover-comparison.md](kover-comparison.md) | Kover 2.0.7 compared with GRM's set covering machine on the 12 smallest lineage-split datasets (1,118 genome-dataset entries). |
| [dense-network.md](dense-network.md) | One-setting dense neural network on 4,000 screened features compared with the five classical learners across 24 datasets in 12 species. |
| [amrfinderplus-comparison.md](amrfinderplus-comparison.md) | Full AMRFinderPlus 4.2.7 program compared with GRM's known-gene rule and five learners across 235 datasets in 21 species. |
| [tbprofiler-comparison.md](tbprofiler-comparison.md) | TBProfiler 6.6.3 run on 2,000 *Mycobacterium tuberculosis* assemblies across 16 named endpoints (14 graded) against GRM's learners and known-gene rule. |
| [baselines.md](baselines.md) | No-genome baselines (`grm baselines`): constant training-fold MIC and lineage-identity baseline. |
| [input-configuration.md](input-configuration.md) | Default input settings, feature families, label rules, and the input comparison, since run on four species. |
| [input-comparison-runner.md](input-comparison-runner.md) | Runner, pairing rules, and second reading for the 4-species input comparison (H-1). |
| [tie-rule.md](tie-rule.md) | Decision tree tie-breaking evaluation record. |
| [hyperopt-selection.md](hyperopt-selection.md) | Hyperparameter tuning evaluation and why fixed defaults remain in place. |
| [lasso-stopping.md](lasso-stopping.md) | Lasso convergence rule and the 18-dataset Revision 14 scan. |
| [ablation-report.md](ablation-report.md) | Schema and validator contract (`ablation-report` 1) for feature-block ablations. |
| [ablation-runner.md](ablation-runner.md) | Runner contract for feature-block ablations. |
| [association-plan.md](association-plan.md) | Planned association analysis specification. |

## Methodology, pipeline, and contracts

| Document | Summary |
| --- | --- |
| [duplicate-isolation.md](duplicate-isolation.md) | Keeping k-mer-identical duplicate genomes together in outer, inner, and final tuning folds (revisions 10 to 12). |
| [equal-scoring.md](equal-scoring.md) | Scoring null fits at training prevalence so every learner covers every isolate (revision 14). |
| [species-check.md](species-check.md) | Species-identity check that drops genomes whose sequence does not fit their labelled species (revision 10). |
| [fold-draws.md](fold-draws.md) | Repeated fold-draw design and stability checks. |
| [fold-file.md](fold-file.md) | Supplied fold-plan file format for `grm learn --fold-file`. |
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
