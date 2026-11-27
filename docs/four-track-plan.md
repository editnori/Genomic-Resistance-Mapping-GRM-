# Four-track evaluation plan

This specification defines four parallel, pre-registered evaluation tracks extending GRM beyond Revision 14 without modifying any released artifact, default pipeline behavior, or committed reference record.

All four tracks are scaffolding and protocol specifications. No real held-out or external cohort results have been read under these tracks, and no released default, model, or site figure changes.

## Shared Revision 14 invariants

Every track adheres to the following frozen invariants:

1. **Revision 14 immutability.** `EVALUATION_REVISION` remains `14` and `REPORT_SCHEMA` remains `2`. Every Revision 14 dataset report retains `LearnerResult.conformal == []`. Pinned input files (`genomes.jsonl`, curated files, and released reports), `scripts/check-results.py`, and all committed records under `reference/` remain untouched.
2. **Sidecars and off-by-default flags.** All track machinery lives in standalone scripts, sidecar JSON/JSONL artifacts, or off-by-default experimental flags. Running `grm` without experimental flags produces byte-identical Revision 14 artifacts.
3. **Retrospective scope labeling.** Every split inside a single public database (`random`, `lineage`, `study-whole`, and insertion-date or collection-year time splits) is retrospective and is never presented as prospective or clinical validation.
4. **Output hygiene.** Summary records and comparison reports hold aggregate counts, metrics, and schema provenance only: no genome identifiers, no BioProject accessions, no per-genome rows or scores, no raw sequences, and no unapproved variant or position tokens. Reports are deterministic across hash seeds and pass automated genome-ID scanning (`scripts/c7-outside/scan.py` rules).
5. **Label-rule separation.** Datasets in set `M` (where intermediate or conflicting calls depend on the pending label rule) remain explicitly marked and separated from primary non-`M` readings.

## Track A: Published AMR tool benchmark harness

Track A establishes a uniform, auditable benchmark comparing published rule-based, species-specific, and trainable antimicrobial resistance tools against GRM's five built-in learners (`scm`, `cart`, `rf`, `gbt`, `lasso`) and GRM's known-gene rule (kept strictly distinct from the standalone `amrfinderplus` program).

### Benchmark tool panel (11 tools)

| Category | Tools |
|---|---|
| Published catalogue and species-specific callers (7) | `amrfinderplus`, `resfinder_pointfinder`, `card_rgi`, `staramr`, `kleborate`, `tbprofiler`, `mykrobe` |
| Trainable sequence and tabular ML baselines (4) | `kover`, `xgboost`, `lightgbm`, `elasticnet` |

### Locked Track A rules

- **Four-tier species-drug support registry.** Before any phenotype is scored, every `(tool, species, antibiotic)` pair is assigned to one of four frozen tiers:
  - **Tier 1 (direct native support):** tool ships explicit curated or graded rules/models for that exact species and antibiotic.
  - **Tier 2 (class-level or proxy mapping):** tool emits a class-level call or a documented proxy drug call mapped by a frozen table (never inferred post hoc).
  - **Tier 3 (partial or ungraded locus support):** tool covers only a subset of target loci or reports ungraded/uncertain hits for that species-drug pair.
  - **Tier 4 (unsupported):** species or antibiotic is outside the tool's scope; scored as unsupported rather than susceptible.
- **Separate non-binary outcome accounting.** Four non-definitive outcomes are counted in separate columns on every dataset and never collapsed into susceptible (`S`) or silently dropped:
  1. `unsupported` (Tier 4 species-drug cell),
  2. `no_call` (tool ran on a supported cell but emitted no definitive call, e.g., missing target locus or low coverage),
  3. `intermediate` (`I` call emitted by the tool),
  4. `failure` (tool crash, timeout, or unparseable output).
- **Pairwise common-support and full-denominator coverage.** Every pairwise comparison between callers reports two complementary views side by side:
  - **Common-support performance:** balanced accuracy, sensitivity, specificity, VME, and ME computed strictly on the intersection of isolates where both callers emit a definitive R/S call.
  - **Full-denominator coverage:** call rate, `unsupported`, `no_call`, `intermediate`, and `failure` shares over the full cohort denominator, accompanied by conservative sensitivity/specificity bounds counting non-calls as wrong or right.
- **Retrospective labeling of modern databases on insertion-date time splits.** When any modern reference catalogue or pretrained database is evaluated on an insertion-date or collection-year temporal split whose cutoff predates the database build, those scores are explicitly labeled **retrospective** because later isolates and resistance determinants may already have entered the modern database. They are never reported as prospective temporal performance.
- **Zero external-cohort winner selection.** No tool, reading mode, threshold, or hyperparameter setting is selected on any external or sealed cohort. All trainable tools tune strictly inside outer training folds on internal retrospective splits; external cohorts are evaluated once under frozen pre-registered rules.

## Track B: Point-mutation representation comparison

Track B tests whether biologically structured representations of point-mutation resistance improve generalization on held-out lineages and held-out studies relative to GRM's 31-mer baseline under identical learner and tuning controls. Full protocol details are specified in [point-mutation-plan.md](point-mutation-plan.md).

### Locked Track B rules

- **Fixed known-catalogue mechanism categories.** Mechanism strata (`point only`, `gene only`, `both`, `neither`) and known target-gene panels are frozen strictly from already-published catalogues (AMRFinderPlus core `AMR`/`POINT`/`POINT_DISRUPT` categories from [mechanism-split.md](mechanism-split.md) and WHO *Mycobacterium tuberculosis* catalogue v2 graded loci from [lineage-curve-result.md](lineage-curve-result.md)).
- **Fixed splits, fixed downstream learner, and fixed tuning budget.** Every representation arm is evaluated on identical `lineage` (`grm-dataset-export-v1`) and `study-whole` (`grm-fold-file-v1`) outer folds, identical inner folds, identical cohort sample support, and the exact same downstream learner (`L2`-regularized class-balanced logistic regression with deterministic inner-fold regularization selection over a frozen log-spaced `lambda` grid).
- **Representation arms compared:**
  1. 31-mer baseline from `grm export` (`kmer31`),
  2. Protein-aligned known-target variants (`ProteinVariantEncoder`) with explicit `unseen_variant_count` side-feature,
  3. Gene and predefined active-site/region burden summaries (`BurdenEncoder`) across synonymous, nonsynonymous, truncating, and indel consequences,
  4. Pretrained foundation-model zero-shot scores and embeddings (`FoundationScoreAdapter`) from **ESM-2**, **Evo 2**, **ProkBERT**, and **Bacformer**,
  5. Concatenated hybrid arms combining 31-mers with structured or foundation-model representations.
- **Explicit exclusion of AlphaGenome 2.** AlphaGenome 2 is excluded for bacterial genomes. Its architecture and training corpus target mammalian eukaryotic regulatory genomics (chromatin accessibility, histone modification, transcription-factor binding, and eukaryotic splicing), which are biologically incompatible with prokaryotic operon structure, intronless coding sequences, and chromosomal point-mutation resistance mechanisms. The exposure validator (`scripts/point-mutations/exposure_audit.py`) hard-refuses AlphaGenome 2 manifests.
- **Three-way pretraining exposure audit (`grm-rep-exposure-v1`).** Every foundation-model arm must carry a validated audit record pinning checkpoint ID, checkpoint SHA-256, training corpus provenance, and three independent exposure axes: `sequence_overlap_status`, `phenotype_annotation_status`, and `eval_tuning_status`. Claim ceilings derive deterministically from those three axes: `"unknown"` sequence overlap caps claims at `"retrospective_only"` and can never obtain `"certified_unseen_sequence"`.
- **Descriptive-only DMS context.** Published deep mutational scanning (DMS) assays serve solely as descriptive biophysical context for zero-shot protein/sequence scores. DMS correlations are never used to select checkpoints, filter positions, or tune thresholds on AMR evaluation cohorts.
- **Zero held-out pocket, position, or variant selection.** All target genes and active-site/region coordinates are frozen from public annotations prior to fitting. Variant vocabularies, standardization statistics, missing-value imputers, and PCA projections are fit strictly inside each outer (and inner) training fold; mutating outer test labels or features must leave every fitted encoder state, selected hyperparameter, and model weight vector bit-for-bit identical.

## Track C: Reusable cross-database deduplication CLI and external corpus ingestion

Track C replaces ad-hoc per-experiment overlap filters (such as the one-off checks documented in [other-database-test.md](other-database-test.md)) with a deterministic, reusable multi-layer deduplication CLI and provenance pipeline for external bacterial AMR corpora.

### Locked Track C rules

- **Multi-layer cross-database deduplication.** Every candidate external isolate is checked against training and previously spent cohorts across four ordered layers:
  1. **INSDC/ENA/NCBI alias map:** resolves cross-database accession equivalences across BioSample (`SAMN`/`SAMEA`/`SAMD`), BioProject (`PRJNA`/`PRJEB`/`PRJDB`), run (`SRR`/`ERR`/`DRR`), GenBank/RefSeq assembly (`GCA`/`GCF`, exact and version-stripped), and BV-BRC genome identifiers.
  2. **Sidecar provenance JSONL:** records source database, accession aliases, collection/insertion metadata, and AST provenance in sidecar JSONL files keyed by content digest without altering `genomes.jsonl` or Revision 14 reports.
  3. **Normalized assembly and canonical k-mer set identity:** computes canonicalized contig-normalized assembly digests alongside dual canonical 31-mer set digests combining `fnv64_k31` (for compatibility with GRM's existing kset digest) and `sha256_k31` (cryptographic collision-resistant canonical 31-mer set digest).
  4. **Species-specific Mash distance thresholds:** applies frozen species-specific Mash distance cutoffs that explicitly separate biological near-relatedness (shared clonal complex, outbreak, or lineage cluster) from same-specimen / re-sequenced / re-assembled technical identity.
- **Unresolved overlaps excluded from certified external cohorts.** Isolates with missing accessions, unresolvable alias chains, or ambiguous sequence-identity status against training or spent cohorts are excluded from certified external cohorts rather than assumed disjoint.
- **Sealed records non-queryable for membership by default.** Sealed holdout cohorts are non-queryable for membership by default so iterative cohort construction cannot probe sealed contents. Before formal unsealing under a frozen analysis plan, any overlap with a sealed cohort is declared `unknown`.

## Track D: Nested selective calling and abstention

Track D evaluates whether fold-safe selective calling (abstaining on uncertain or genetically distant isolates) reduces clinical-style error rates among called isolates without hiding resistant cases or degrading calibration, using strictly nested calibration outside Revision 14 reports.

### Locked Track D rules

- **Nested group-separated inner-fold calibration.** Abstention thresholds and combined uncertainty-distance rules are calibrated strictly inside each outer training fold via group-separated inner folds (following the `src/learn/inner_pick.rs` pattern: lineage-grouped inner folds on `lineage` and group-separated inner folds on `study-whole`, keeping k-mer-identical duplicates whole). Outer test scores and labels never enter calibration.
- **Fold-safe Mash distances (`d_inner` / `d_outer`).** Nearest-training-genome Mash distances are computed strictly within each training partition: inner validation isolates measure distance exclusively to their inner training partition (`d_inner`), and outer test isolates measure distance exclusively to their outer training partition (`d_outer`). No test-to-test or held-out-to-train distance leakage is permitted.
- **Frozen primary coverage and full reporting ladder.**
  - Primary target coverage: `c = 0.85` (frozen pre-run).
  - Full secondary coverage ladder: `{1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50}`.
- **Frozen tie rule and sparse-denominator policy.**
  - Ties at any selective-calling threshold are resolved by a deterministic, outcome-blind tie rule frozen pre-run (never consulting true labels).
  - Sparse-denominator policy reports raw confusion counts on every cell and enforces a frozen minimum called-class denominator before computing conditional error rates (called VME, called ME, called NPV); cells below the floor are marked unestimable rather than coerced to zero or silently omitted.
- **Dual operating-point reporting.** Every dataset reports two distinct operating views separately:
  1. **Training-calibrated operating point:** threshold chosen inside the outer training fold for target `c` and applied unchanged to the outer test fold, reporting realized coverage honestly under lineage or study shift.
  2. **Exact `floor(c * n)` matched-coverage comparison:** ranks test isolates by abstention score to exact count `floor(c * n)` to isolate ranking quality across two pre-registered contrasts:
     - probability margin vs outcome-blind random thinning, and
     - margin + fold-safe Mash distance vs margin alone.
- **Asymmetric clinical-style metrics and cross-split species sign test.** Every operating point reports `coverage_r` (coverage among true resistant isolates), `coverage_s` (coverage among true susceptible isolates), called VME, called ME, called NPV, and full-denominator `FN / n_r` (preventing abstention from masking missed resistant isolates). Primary statistical inference uses a species-level one-sided exact binomial sign test (direction only, aggregating datasets within each species before testing across species), and a claim requires significance on **both** `lineage` and `study-whole` splits.

## Cross-track integration boundaries

| Track | Artifact / code boundary | Revision 14 reports changed? | Real scores read? |
|---|---|---|---|
| Track A | External tool adapters, 4-tier registry, common-support scorer | No | Not under this plan |
| Track B | `scripts/point-mutations/`, `tests/point_mutations/`, `docs/point-mutation-plan.md` | No | Not run yet |
| Track C | Cross-database deduplication CLI and sidecar provenance JSONL | No | Not run yet |
| Track D | Nested selective calling sidecar evaluator (`c = 0.85` ladder) | No (`conformal == []` unchanged) | Not run yet |
