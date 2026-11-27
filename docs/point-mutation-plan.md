# Point-mutation representation comparison plan (Track B)

Not run yet. This document fixes the pre-run protocol, pilot species-drug panel, three-exposure pretraining audit rules, descriptive deep mutational scanning (DMS) boundaries, and fixed-learner comparison contract for Track B of [four-track-plan.md](four-track-plan.md).

The scaffolding lives in `scripts/point-mutations/` and synthetic unit plus planted-fault leakage tests live in `tests/point_mutations/`. No Revision 14 report, `src/` file, `scripts/check-results.py` check, or committed `reference/` record is modified. Every split inside BV-BRC (`lineage`, `study-whole`, and reference `random`) is a retrospective split inside one public database and is never external or clinical validation.

## Background and terminology alignment

Track B builds directly on two completed Revision 14 diagnostic readings and two frozen split contracts:

- **Mechanism split ([mechanism-split.md](mechanism-split.md), [mechanism-split-result.md](mechanism-split-result.md)).** Across 237 datasets in 21 species, AMRFinderPlus 4.2.7 (`database 2026-08-07.1`) stratifies resistant genomes into `point only` (`POINT` or `POINT_DISRUPT` rows and no `AMR` row), `gene only` (`AMR` row and no point row), `both`, and `neither`. In the primary non-`M` set, 47 of 204 datasets hold 10 or more `point only` resistant genomes, and 19 datasets in 5 species hold 10 or more in both `gene only` and `point only`. On those 19 datasets, all five 31-mer learners trail on `point only` relative to `gene only` across cephalosporin, carbapenem, penicillin-inhibitor, and aminoglycoside datasets while leading on quinolone datasets. GRM's known-gene rule is kept strictly separate from the AMRFinderPlus program and is never equated with it.
- ***Mycobacterium tuberculosis* lineage curve ([lineage-curve-result.md](lineage-curve-result.md), [tbprofiler-comparison.md](tbprofiler-comparison.md)).** Over 2,000 BV-BRC *M. tuberculosis* assemblies scored with TBProfiler 6.6.3 (`tbdb` commit `3ed1cb99`, WHO catalogue v2), chromosomal point mutations and coding disruptions exhibit markedly different lineage dispersion across drugs: 78 of 94 graded *pncA* variants for pyrazinamide and 121 of 136 graded variants (predominantly *ethA*) for ethionamide occur in only one main lineage, creating severe rare-variant generalization pressure under lineage holdout where exact 31-mer matching misses unseen alleles.
- **Export and split contracts ([dataset-export.md](dataset-export.md), [fold-file.md](fold-file.md), [outside-learners.md](outside-learners.md), [study-split.md](study-split.md)).**
  - `lineage` (and optional reference `random`, always labeled "training and test share lineages") reads outer folds, inner folds, sample ordering, labels, and the 31-mer binary matrix from `grm export` (`format: grm-dataset-export-v1`, `evaluation_revision: 14`), verifying `matrix_sha256` and canonical `plan.sha256`.
  - `study-whole` reads outer BioProject-grouped folds from a reviewed `grm-fold-file-v1` file, keeping k-mer-identical duplicate groups intact and drawing deterministic group-preserving 3-fold inner partitions inside each outer training fold (matching the `scripts/c7-outside/outside.py` discipline).

## Pilot point-mutation-heavy species-drug panel

All pilot datasets and target loci are fixed strictly from already-published resistance mechanism catalogues (WHO *M. tuberculosis* catalogue v2 graded loci and AMRFinderPlus core `POINT`/`POINT_DISRUPT` organism panels). No target gene, residue coordinate, or active-site boundary is selected from held-out evaluation performance. Datasets in set `M` remain explicitly marked and separated from the primary non-`M` summary.

### Group 1: *Mycobacterium tuberculosis* chromosomal point-mutation panel

| Dataset (antibiotic) | Label set | Published target loci (WHO v2) | Mechanism profile from M4 |
|---|---|---|---|
| `rifampin` | primary (`not M`) | *rpoB* (RRDR codons 426-452 and coding region) | Concentrated hotspot + 20/40 single-lineage variants |
| `ethambutol` | primary (`not M`) | *embB*, *embA* promoter | Recurrent hotspot substitutions across lineages |
| `pyrazinamide` | primary (`not M`) | *pncA* (full coding sequence and active-site pocket) | Dispersed loss-of-function (78/94 in one main lineage) |
| `streptomycin` | primary (`not M`) | *rpsL*, *rrs*, *gid* | Mixed canonical substitutions + *gid* loss-of-function |
| `levofloxacin` | primary (`not M`) | *gyrA*, *gyrB* (QRDR) | Topoisomerase QRDR point substitutions |
| `moxifloxacin` | primary (`not M`) | *gyrA*, *gyrB* (QRDR) | Topoisomerase QRDR point substitutions |
| `ethionamide` | primary (`not M`) | *ethA*, *inhA* | High-diversity *ethA* disruption (121/136 in one main lineage) |
| `isoniazid` | `M` (reported apart) | *katG*, *inhA* | Canonical *katG* S315 + dispersed *katG* loss-of-function |

### Group 2: Gram-negative fluoroquinolone target panel (point + acquired coexistence)

| Species | Antibiotics | Label set | Published target loci (AMRFinderPlus core) |
|---|---|---|---|
| *Escherichia coli* | `ciprofloxacin`, `levofloxacin`, `nalidixic acid` | primary (`not M`) | *gyrA*, *gyrB*, *parC*, *parE* (QRDR); efflux regulators *acrR*, *marR*, *soxR* |
| *Klebsiella pneumoniae* | `ciprofloxacin`, `levofloxacin` | primary (`not M`) | *gyrA*, *parC* (QRDR); *oqxR*, *ramR* |
| *Salmonella enterica* | `ciprofloxacin`, `nalidixic acid` | primary (`not M`) | *gyrA*, *gyrB*, *parC*, *parE* (QRDR) |
| *Shigella sonnei* | `ciprofloxacin`, `nalidixic acid` | primary (`not M`) | *gyrA*, *parC* (QRDR) |
| *Pseudomonas aeruginosa* | `ciprofloxacin`, `levofloxacin` | primary (`not M`) | *gyrA*, *gyrB*, *parC*, *parE*; *mexR*, *nfxB* |
| *Acinetobacter baumannii* | `ciprofloxacin`, `levofloxacin` | primary (`not M`) | *gyrA*, *parC* (QRDR) |

### Group 3: Gram-negative beta-lactam, porin, and beta-lactamase target panel

| Species | Antibiotics | Label set | Published target loci (AMRFinderPlus core) |
|---|---|---|---|
| *Klebsiella pneumoniae* | `cefuroxime`, `ceftriaxone`, `ceftazidime`, `cefoxitin`, `ertapenem`, `doripenem`, `imipenem`, `ampicillin/sulbactam` | primary (`not M`) | Porins *ompK35*, *ompK36* (`POINT` / `POINT_DISRUPT`); active-site variants in *blaSHV*, *blaTEM*, *blaCTX-M*, *blaKPC* |
| *Escherichia coli* | `cefuroxime`, `cephalothin`, `ertapenem` | primary (`not M`) | *ampC* promoter/attenuator, *ftsI* (PBP3), *ompC*, *ompF*; *blaTEM*, *blaCTX-M* |
| *Klebsiella pneumoniae* | `meropenem`, `piperacillin/tazobactam` | `M` (reported apart) | *ompK35*, *ompK36*, *blaKPC*, *blaSHV* |
| *Escherichia coli* | `amoxicillin/clavulanate`, `ampicillin`, `meropenem` | `M` (reported apart) | *ampC*, *ftsI*, *ompC*, *ompF*, *blaTEM* |

### Group 4: Secondary point-target stress controls (exploratory beside primary panel)

- *Pseudomonas aeruginosa* `amikacin` and `tobramycin` (*fusA1*, *mexZ* point mutations beside acquired aminoglycoside-modifying enzymes).
- *Staphylococcus aureus* `ciprofloxacin`, `rifampin`, and `fusidic acid` (*gyrA*, *grlA*, *rpoB*, *fusA*).
- *Enterococcus faecium* beta-lactam and lipopeptide point targets (*pbp5*, *liaFSR*) where available in `S` or `L`.

Datasets enter `lineage` when their Revision 14 export has an available lineage plan, and enter `study-whole` when they appear in the reviewed study-split set `S` (`docs/study-split.md`). No dataset is dropped post hoc after reading any score.

## Representation arms and fold-leak-free extraction

All feature builders in `scripts/point-mutations/features.py` follow a strict `fit(train)` / `transform(eval)` separation inside every outer training fold and every inner training fold:

| Arm | Builder | Feature definition | Fold-fitted state (training fold only) |
|---|---|---|---|
| `kmer31` (baseline) | `Kmer31FoldEncoder` | `grm export` binary 31-mer presence matrix, candidate support filter `[floor, n_train - floor]`, deduplicated training patterns, chi-square top-4,000 | Candidate set, top-4,000 columns, column mean and std |
| `protein_variants` | `ProteinVariantEncoder` | Presence/absence of protein-aligned amino-acid substitutions, indels, and truncations on frozen published target genes, plus explicit `unseen_variant_count` side-feature counting distinct eligible variants absent from the fitted vocabulary | Frozen sorted vocabulary (`train_only`, `catalogue_only`, or `train_union_catalogue`) |
| `gene_burden` | `BurdenEncoder` | Per-gene and per-predefined-region (e.g., RRDR, QRDR, catalytic pocket) counts of distinct `synonymous`, `nonsynonymous`, `truncating`, and `indel` variants from frozen public coordinates | Training-fold column mean and std (zero-variance columns mapped to `0.0`) |
| `foundation_<family>` | `FoundationScoreAdapter` | Precomputed per-isolate zero-shot log-likelihood ratios (LLRs) and/or embeddings from `esm2`, `evo2`, `prokbert`, or `bacformer` aligned strictly to `dataset.json` or `grm-fold-file-v1` row manifests | Optional train-fold mean imputer, train-fold column mean and std, optional sign-canonical SVD/PCA projection |
| `hybrid_*` | `make_concat_arm` | Column concatenation of `kmer31` with `protein_variants`, `gene_burden`, and/or `foundation_<family>` | Independent fold-safe fit of each component arm strictly on training indices |

### Explicit model-family exclusions

- **AlphaGenome 2 (`alphagenome2`) is explicitly excluded for bacterial genomes.** AlphaGenome 2 models mammalian eukaryotic regulatory tracks (chromatin accessibility, histone states, transcription-factor binding, and eukaryotic RNA splicing). Bacterial AMR point mutations operate primarily in intronless prokaryotic coding regions and bacterial promoters/ribosomal RNAs; applying a mammalian regulatory sequence model is biologically invalid. `exposure_audit.py` hard-refuses any manifest naming `alphagenome2` (or variants of that name) or declaring an `eukaryotic_only` / `mammalian_only` taxonomic scope.

## Three-exposure pretraining audit protocol (`grm-rep-exposure-v1`)

Every foundation-model manifest (`grm-rep-foundation-v1`) MUST embed a valid `grm-rep-exposure-v1` record validated by `scripts/point-mutations/exposure_audit.py`. Manifests missing an audit record, failing schema validation, or exceeding their derived claim ceiling are refused before any model fit begins.

### Required fields

- `format`: `"grm-rep-exposure-v1"`
- `model_family`: one of `esm2`, `evo2`, `prokbert`, `bacformer`
- `checkpoint_id`: pinned checkpoint identifier string
- `checkpoint_sha256`: 64-character lowercase hexadecimal SHA-256 of checkpoint weights/archive
- `training_corpus`: object with `name`, `version_or_cutoff`, non-empty `sources` list, and `taxonomic_scope`
- `sequence_overlap_status`: `verified_disjoint` | `partial_overlap` | `confirmed_overlap` | `unknown`
- `phenotype_annotation_status`: `none_verified` | `indirect_catalogue` | `supervised_phenotype` | `unknown`
- `eval_tuning_status`: `frozen_preregistered` | `outer_train_only` | `eval_touched` | `unknown`
- `audit_evidence`: object with `method`, `checked_at_utc` (ISO-8601 UTC or `"unverified"`, with timestamp required for `verified_disjoint`), and `notes`
- `claimed_ceiling` (optional): if provided, must not exceed the automatically derived ceiling

### Automatic claim-ceiling derivation

| `eval_tuning_status` | `phenotype_annotation_status` | `sequence_overlap_status` | Derived `claim_ceiling` |
|---|---|---|---|
| `eval_touched` or `unknown` | any | any | `contaminated_exploratory` |
| `frozen_preregistered` or `outer_train_only` | `indirect_catalogue`, `supervised_phenotype`, or `unknown` | any | `non_independent_retrospective` |
| `frozen_preregistered` or `outer_train_only` | `none_verified` | `unknown`, `partial_overlap`, or `confirmed_overlap` | `retrospective_only` |
| `frozen_preregistered` or `outer_train_only` | `none_verified` | `verified_disjoint` | `certified_unseen_sequence` |

In particular, an unaudited pretraining corpus (`sequence_overlap_status = "unknown"`) always caps claims at `"retrospective_only"` and can never claim `"certified_unseen_sequence"`. Concatenated hybrid arms inherit the strictest (lowest-rank) claim ceiling across their component foundation arms.

## Descriptive-only DMS alignment rules

Published deep mutational scanning (DMS) fitness datasets for antibiotic target proteins (e.g., TEM-1 beta-lactamase, FolA/DHFR, or other experimental enzyme scans) may be aligned against zero-shot foundation-model log-likelihood ratios solely as **descriptive biophysical context**:

1. **Permitted use:** reporting descriptive rank concordance (Spearman correlation) between frozen zero-shot variant scores and experimental fitness measurements to document whether a checkpoint captures bacterial protein biophysics.
2. **Prohibited uses:**
   - Selecting checkpoints, transformer layers, context windows, or prompt templates based on DMS scores or evaluation outcomes.
   - Deriving variant inclusion thresholds, pocket masks, or feature weights from DMS assays tuned against AMR phenotype splits.
   - Substituting DMS concordance for held-out lineage or study-split AMR evaluation.

## Fixed-learner experimental controls (`grm-rep-compare-v1`)

To attribute any performance difference solely to the genome representation rather than learner capacity or tuning budget, `scripts/point-mutations/compare.py` holds every non-representation factor strictly fixed across arms:

- **Downstream learner:** deterministic class-balanced L2-regularized logistic regression solved to exact KKT tolerance (`tol = 1e-9`, `max_iter = 60`) in pure NumPy using primal Newton-IRLS (`p <= n`) and exact Woodbury SPD Newton-IRLS (`p > n`).
- **Tuning budget and tie rule:** identical log-spaced regularization ladder `lambda in {1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0}` selected by lowest pooled inner-fold class-balanced log loss (rounded to 12 decimals), breaking ties toward larger `lambda`.
- **Null-fit parity:** any fold yielding zero features or exact zero gradient falls back to the outer training fold's resistant prevalence and is counted in `null_fit_folds` and `constant_folds`, matching GRM Revision 14 equal-scoring rules.
- **Metrics and cutoff:** resistant call at `score >= 0.5`, clipped log loss on `[1e-4, 1 - 1e-4]`, reporting per-fold counts (`n`, `tp`, `fp`, `tn`, `fn`, `loss_sum`), pooled `balanced_accuracy`, `recall`, `specificity`, `vme`, `me`, `log_loss`, and exact paired arm-minus-baseline differences on identical tested support (`common_tested_n`).
- **Zero-ID report hygiene:** `grm-rep-compare-v1` enforces an exact key allowlist with zero genome IDs, zero per-genome predictions, zero feature/variant names, zero weights, and zero timestamps, verified against `scripts/c7-outside/scan.py`.

## Verification and planted-fault leakage guards

`tests/point_mutations/` verifies all components on purely synthetic fixtures (no real genomes, no network, standard library + NumPy 1.26.4 only):

1. **Exposure audit tests (`test_exposure_audit.py`):** verifies all claim-ceiling combinations, `"unknown"` sequence overlap capping at `"retrospective_only"`, AlphaGenome 2 refusal, over-claim rejection, and malformed/missing audit rejection.
2. **Feature builder tests (`test_features.py`):** verifies `ProteinVariantEncoder` (including `unseen_variant_count`), `BurdenEncoder` (frozen gene/region consequence aggregation and train-only scaling), and `FoundationScoreAdapter` (strict manifest alignment to `grm-dataset-export-v1` and `grm-fold-file-v1`, audit gating, train-only imputation/scaling/PCA).
3. **Fixed-learner comparator tests (`test_compare.py`):** verifies exact KKT optimality of both primal and Woodbury L2 logistic solvers, end-to-end multi-arm comparisons on synthetic `lineage` and `study-whole` fixtures, deterministic JSON bytes across hash seeds, and zero genome-ID hits under `scripts/c7-outside/scan.py`.
4. **Planted-fault leakage tests (`test_leakage.py`):** mutates every outer test label, variant list, burden count, 31-mer bit, and foundation embedding vector, and asserts that every fitted vocabulary, burden scaler, PCA projection matrix, selected inner `lambda`, and training-fold weight vector + intercept is bit-for-bit identical; then plants five concrete leakage faults (full-cohort vocabulary, full-cohort burden scaling, full-cohort PCA, outer-test `lambda` selection, and outer-test refit inclusion) and confirms every planted fault is caught.
