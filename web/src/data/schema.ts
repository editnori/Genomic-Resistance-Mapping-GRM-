// The dashboard's data contract. Every type mirrors JSON written by the Rust
// side (`grm report` in src/report.rs, per-dataset reports in
// src/analysis.rs). Change both sides together; `schema` versions the files.
//
// Serialization notes from serde: f64 NaN becomes null; Option<T> is T | null;
// tuples are arrays; enums without a rename serialize as their variant name.

/** The schema `grm` writes. Schema 2 lets a learner row hold null scores (no
 *  call), null fold balanced accuracy, and null metrics or rates when it
 *  called no isolate, and lets a catalog row have no best learner. A
 *  dashboard that reads only schema 1 stops at the version check. */
export const SCHEMA = 2;

/** Schema 1 files, from revisions 5 and 6, have none of those nulls, so this
 *  dashboard reads them unchanged. */
export const READ_SCHEMAS: readonly number[] = [1, 2];

export type Learner = "scm" | "cart" | "rf" | "gbt" | "lasso" | "catalog";
export type Scheme = "random" | "lineage";

/** Clinical agreement on one set of predictions. VME = resistant called
 *  susceptible over all resistant; ME = susceptible called resistant over all
 *  susceptible. */
export interface BinaryMetrics {
  n: number;
  tp: number;
  fp: number;
  tn: number;
  fn: number;
  accuracy: number;
  balanced_accuracy: number;
  sensitivity: number;
  specificity: number;
  precision: number;
  f1: number;
  mcc: number;
  vme: number;
  me: number;
  auroc: number | null;
  auprc: number | null;
}

/** The summary subset used in catalog tables. */
export interface Lite {
  balanced_accuracy: number;
  sensitivity: number;
  specificity: number;
  vme: number;
  me: number;
  auroc: number | null;
}

/** A learner's Lite on one split. Its rates are null when it called no
 *  isolate; otherwise they cover the isolates it called. */
export interface LearnerLite {
  balanced_accuracy: number | null;
  sensitivity: number | null;
  specificity: number | null;
  vme: number | null;
  me: number | null;
  auroc: number | null;
  /** Isolates without a call. Absent on rows written before no-calls. */
  no_call?: NoCall;
}

export interface ConformalPoint {
  alpha: number;
  /** Share of isolates that get a single call. */
  coverage: number;
  vme: number;
  me: number;
  /** Share that get both labels: "cannot call". */
  ambiguous: number;
}

export interface LearnerResult {
  learner: Learner;
  /** Pooled out-of-fold calls; isolates without a call are left out. Null
   *  when the learner called no isolate. */
  metrics: BinaryMetrics | null;
  /** One per outer fold; null for a fold whose fit made no call. */
  fold_balanced_accuracy: (number | null)[];
  /** ROC points [fpr, tpr]; empty for 0/1 learners (SCM, catalog). */
  roc: [number, number][];
  /** Out-of-fold resistant score per sample, in `samples` order; null for
   *  an isolate without a call. */
  scores: (number | null)[];
  /** Empty while the pipeline fits no separate calibration set. */
  conformal: ConformalPoint[];
  /** Per-fold model size: SCM rules, CART leaves, forest trees, boosted
   *  trees kept by early stopping, or lasso non-zero weights. Zero for
   *  boosting or the lasso is a null fit: rows without `no_call` scored
   *  that fold 0.5, which the 0.5 cutoff calls resistant; rows with it
   *  leave the fold uncalled. */
  fold_complexity: number[];
  seconds: number;
  /** Isolates in folds whose fit was null, by class. Absent on catalog rows
   *  and on reports written before null fits became no-calls. */
  no_call?: NoCall;
}

export interface SplitResult {
  scheme: Scheme;
  folds: number;
  n_groups: number;
  /** Outer fold that held out each sample, in `samples` order. Absent when
   *  the split is unavailable and in reports written before it was saved. */
  test_fold?: number[];
  /** Empty when the split could not be evaluated (see unavailable_reason). */
  learners: LearnerResult[];
  /** Why the planned partitions could not support evaluation; present only
   *  then, with no scores and no catalog row. */
  unavailable_reason?: string;
}

export type Censor = "Exact" | "AtMost" | "Above";

/** log2 dilution step; 0 = 1 mg/L. */
export interface Mic {
  log2: number;
  censor: Censor;
}

export interface Breakpoint {
  /** Resistant when MIC > 2^r_above mg/L. */
  r_above: number;
  agreement: number;
  n: number;
}

export interface Sample {
  id: string;
  resistant: boolean;
  /** Whether the label was inferred from MIC; false in audited reports,
   *  whose samples carry reported calls only. */
  from_mic: boolean;
  mic: Mic | null;
  lineage: number;
  country: string | null;
  year: string | null;
  host: string | null;
  mlst: string | null;
  /** Indices into DatasetReport.catalog_groups. */
  groups: number[];
  catalog_resistant: boolean | null;
}

export interface GroupRow {
  group: string;
  class: string;
  subclass: string;
  counts_for_rule: boolean;
  n_resistant: number;
  n_susceptible: number;
}

export interface CatalogBaseline {
  version: string;
  /** AND of groups, each an OR of "C:<class>" / "S:<subclass>" tokens. */
  rule: string[][];
  metrics: BinaryMetrics;
  groups: GroupRow[];
}

export interface PanelResult {
  id: string;
  name: string;
  mirrors: string;
  targets: string[];
  metrics: BinaryMetrics;
  /** [group, resistant carriers, susceptible carriers] */
  detected: [string, number, number][];
  no_call?: NoCall;
}

export interface NoCall {
  resistant: number;
  susceptible: number;
}

export interface UnsupportedAssay {
  id: string;
  name: string;
  mirrors: string;
  reason: string;
}

export interface MicResult {
  n: number;
  n_exact: number;
  essential_agreement: number;
  essential_agreement_all: number;
  /** Calls from the predicted MIC at each training fold's own threshold,
   *  scored against the reported R/S labels; null when a fold has none. */
  categorical: BinaryMetrics | null;
  /** The cohort breakpoint (labels isolates with an MIC but no call); not the
   *  threshold the categorical metrics use. */
  breakpoint: Breakpoint | null;
  /** One empirical threshold per training fold; absent in reports written
   *  before fold thresholds. */
  fold_breakpoints?: (Breakpoint | null)[];
  /** [observed log2, predicted log2, count] over exact MICs. */
  confusion: [number, number, number][];
  /** Out-of-fold predicted log2 MIC per sample; null without an MIC. */
  predicted: (number | null)[];
}

export interface Rule {
  feature: number;
  presence: boolean;
}

export interface ScmStep {
  covered: number;
  lost: number;
  utility: number;
  equivalent: Rule[];
  n_equivalent: number;
}

export interface FinalScm {
  kind: "conjunction" | "disjunction";
  p: number;
  cv_error: number;
  rules: Rule[];
  steps: ScmStep[];
  /** Increase in training balanced error when each rule is dropped. */
  importance: number[];
}

export interface TreeNode {
  value: number;
  n: number;
  n_resistant: number;
  depth: number;
  split: { feature: number; gain: number; absent: number; present: number } | null;
}

export interface WeightedFeature {
  feature: number;
  weight: number;
}

export interface FinalModels {
  scm: FinalScm;
  cart: { nodes: TreeNode[] };
  cart_depth: number;
  lasso: { intercept: number; coef: [number, number][]; lambda: number } | null;
  rf_importance: WeightedFeature[];
  gbt_importance: WeightedFeature[];
}

export interface Locus {
  genome_id: string;
  contig: string;
  start: number;
  end: number;
  kmers_placed: number;
  gene: string | null;
  product: string | null;
  patric_id: string | null;
  gene_strand: string | null;
  codons: [number, number] | null;
  /** [gene or product, signed distance in bp] when intergenic. */
  nearest: [string, number] | null;
  /** Contig sequence around the span; the span itself is uppercase. */
  context: string;
}

export interface FeatureCard {
  feature: number;
  n_kmers: number;
  kmer: string;
  catalog_families: string[];
  catalog_products: string[];
  locus: Locus | null;
  n_resistant: number;
  n_susceptible: number;
  /** Presence over samples: base64 of little-endian u64 words. */
  bits: string;
}

/** Whether a report's metrics pass the internal evaluation audit, and what
 *  that covers. Absent in reports written before the audit (read as not valid). */
export interface EvaluationAudit {
  revision: number;
  valid: boolean;
  /** Why the metrics are not valid; empty when they are. */
  reasons: string[];
  /** What a valid audit covers (internal cross-validation, not external or
   *  clinical validation). Absent before audit revision 3. */
  scope?: string;
}

/** The cross-validation and learner settings behind a report. Seeds that
 *  each fit derives from `seed` and its outer fold are not recorded. */
export interface EvalConfig {
  folds: number;
  inner_folds: number;
  scm_max_rules: number;
  cart_max_depth: number;
  forest: { n_trees: number; subsample: number; max_depth: number; min_leaf: number };
  gbt: { rounds: number; learning_rate: number; max_depth: number; lambda: number; min_child_weight: number; colsample: number; patience: number };
  screened_features: number;
  lasso_lambdas: number;
  seed: number;
}

export interface DatasetReport {
  schema: number;
  id: string;
  evaluation?: EvaluationAudit;
  /** Absent in reports written before it was recorded. */
  eval_config?: EvalConfig;
  species: string;
  species_slug: string;
  antibiotic: string;
  drug_class: string;
  breakpoint: Breakpoint | null;
  counts: {
    samples: number;
    resistant: number;
    susceptible: number;
    /** Samples whose label was inferred from MIC; 0 in audited reports. */
    labels_from_mic: number;
    /** Matrix genomes whose resistant or susceptible label for this drug was
     *  inferred from MIC, left out of training and testing; out of samples +
     *  this. Absent in reports written before it was counted. */
    labels_from_mic_excluded?: number;
    features: number;
    patterns: number;
    kmers: number;
    /** k of the species matrix; absent in reports written before it was recorded. */
    k?: number;
    lineages: number;
    min_support: number;
    /** Fewest resistant and fewest susceptible samples the dataset needed.
     *  Absent in reports written before it was recorded, when `grm learn`
     *  gated every dataset at 25; those read as 25. */
    min_class?: number;
  };
  lineage_threshold: number;
  catalog_groups: string[];
  samples: Sample[];
  splits: SplitResult[];
  catalog: CatalogBaseline | null;
  panels: PanelResult[];
  unsupported_assays?: UnsupportedAssay[];
  mic: MicResult | null;
  models: FinalModels;
  features: FeatureCard[];
  seconds: number;
}

export interface Funnel {
  records: number;
  records_class_label: number;
  records_not_antibiotic: number;
  records_unresolved_name: number;
  records_no_measure: number;
  genomes_with_records: number;
  genomes_without_metadata: number;
  genomes_failed_qc: Record<string, number>;
  genomes_passing_qc: number;
  pairs: number;
  pairs_conflict: number;
  pairs_intermediate_only: number;
  pairs_reported: number;
  pairs_labeled_from_mic: number;
  pairs_unlabeled: number;
}

export interface SpeciesRow {
  species: string;
  slug: string;
  datasets_curated: number;
  datasets_reported: number;
  genomes_labeled: number;
  genomes_analyzed: number;
  lineages: number | null;
  lineage_threshold: number | null;
}

export interface CatalogRow {
  id: string;
  species: string;
  species_slug: string;
  antibiotic: string;
  drug_class: string;
  /** The dataset report's audit; report.rs writes revision 0 for legacy
   *  reports. Absent in manifests written before the audit. */
  evaluation?: EvaluationAudit;
  n: number;
  n_resistant: number;
  labels_from_mic: number;
  /** The report's counts.labels_from_mic_excluded, when it records one. */
  labels_from_mic_excluded?: number;
  /** The report's counts.min_class; report.rs writes 25 for a report that
   *  records none. Absent in manifests built before report.rs wrote it,
   *  which read as 25 too. */
  min_class?: number;
  features: number;
  lineages: number;
  best_learner: Learner | null;
  best_random: number | null;
  best_lineage: number | null;
  lineage_unavailable_reason?: string;
  learners: { learner: Learner; random: LearnerLite | null; lineage: LearnerLite | null }[];
  catalog: Lite | null;
  mic_essential_agreement: number | null;
  /** Plain-words SCM rules, e.g. "presence of blaKPC". */
  scm_rules: string[];
  /** [panel id, metrics] */
  panels: [string, Lite][];
  /** Published Kover result for the same species and drug, if any. */
  kover: { sensitivity: number; specificity: number; balanced_accuracy: number; n_rules: number; n_examples: number } | null;
}

export type Confusion = Pick<BinaryMetrics, "n" | "tp" | "fp" | "tn" | "fn">;

export interface PanelSummaryRow {
  panel: string;
  name: string;
  dataset: string;
  species: string;
  antibiotic: string;
  metrics: Lite;
  best_model: Lite;
  best_learner: Learner;
  counts?: Confusion;
  best_model_counts?: Confusion;
  no_call?: NoCall;
}

export interface UnsupportedPanelRow {
  panel: string;
  name: string;
  dataset: string;
  species: string;
  antibiotic: string;
  reason: string;
}

export interface Manifest {
  schema: number;
  generated_at: string;
  evaluation_revision?: number;
  /** The counts.k every published report records; null when two differ, one
   *  predates the field, or none is published. Older manifests hard-code 31. */
  k: number | null;
  amr_records: number;
  amrfinder_version: string;
  name_vocabulary: string | null;
  funnel: Funnel;
  species: SpeciesRow[];
  datasets: CatalogRow[];
  panels: PanelSummaryRow[];
  unsupported_panels?: UnsupportedPanelRow[];
}

export interface Merge {
  a: number;
  b: number;
  height: number;
  size: number;
}

export interface GenomeRow {
  id: string;
  lineage: number;
  x: number;
  y: number;
  /** Mash distance to the closest other genome. */
  nearest: number;
  country: string | null;
  year: string | null;
  host: string | null;
  mlst: string | null;
  length: number;
  contigs: number;
  n50: number;
  groups: number[];
  /** One character per drug in SpeciesReport.drugs: R, S, I, or "." (untested). */
  profile: string;
}

export interface SpeciesReport {
  schema: number;
  species: string;
  slug: string;
  drugs: string[];
  drug_classes: string[];
  catalog_groups: string[];
  lineage_threshold: number;
  n_lineages: number;
  /** Dendrogram leaf order (genome indices); merges use scipy-style ids. */
  leaf_order: number[];
  merges: Merge[];
  genomes: GenomeRow[];
  cross: { a: string; b: string; n: number; rr: number; rs: number; sr: number; ss: number; phi: number }[];
}

export interface DrugNameRow {
  records: number;
  class: string;
  raw: string;
  canonical: string;
  kind: "drug" | "class" | "excluded" | "unresolved";
  components: { text: string; canonical: string; source: string }[];
}

export interface Corpus {
  schema: number;
  by_species: { species: string; pairs_resistant: number; pairs_susceptible: number; genomes: number }[];
  years: [string, number][];
  countries: [string, number][];
  methods: [string, number][];
  standards: [string, number][];
  breakpoints: {
    species: string;
    antibiotic: string;
    r_above_log2: number;
    r_above_mg_l: number;
    agreement: number;
    n: number;
    labels_from_mic: number;
  }[];
  drug_names: DrugNameRow[];
}

// research.json, written by `grm research` (src/research.rs): revision-5
// research pinned to the revision-5 cohorts it describes. An entry applies to
// a loaded dataset report only through researchFor in research.ts, which
// checks the entry's guard against that report.

export const RESEARCH_SCHEMA = 1;
/** sha256 of RESEARCH-SNAPSHOT-SHA256SUMS, the revision-5 snapshot both
 *  research lanes ran on. */
export const REVISION5_SNAPSHOT_SHA256 = "75a724c1c362e6f5c9555886788c57cbf99f617d25f8a91dc8beb4c54b0b92fe";

/** SHA-256 digests of the report values a lane read (encoding in research.ts). */
export interface CohortGuard {
  samples: number;
  sample_ids: string;
  labels: string;
  lineages: string;
}

export interface MicGuard extends CohortGuard {
  mic_values: string;
  /** null when the revision-5 report has no MIC model. */
  mic_predicted: string | null;
}

export interface SelectiveGuard extends CohortGuard {
  /** Every split's held-out scores: the lineage split for the verdicts, the
   *  random split for category C. */
  scores: string;
}

/** One prediction set, as in the MIC lane's prediction_quality.csv. A
 *  baseline run scores the same isolates as the model run on its split, so
 *  the two share n_predicted and n_exact. */
export interface MicRun {
  n_predicted: number;
  n_exact: number;
  /** null when no exact MIC was scored. */
  ea_exact: number | null;
  /** Lineage-split model runs: folds whose early-stopping holdout fell back
   *  to a random holdout. */
  holdout_random_folds?: number;
}

export interface MicDataset {
  id: string;
  guard: MicGuard;
  zone_diameters: {
    /** Report MICs the rerun dropped because they were zone diameters. */
    removed: number;
    /** Scored MICs that are zone diameters; always 0, because a model that
     *  scores any is in mic.withheld instead. */
    scored: number;
  };
  /** A rerun has both splits; a saved report model has the random split only. */
  model: { source: "rerun"; random: MicRun; lineage: MicRun } | { source: "saved_report"; random: MicRun; lineage: null };
  /** The no-genome baseline on the rerun's folds; null for a saved report model. */
  baseline: { random: MicRun; lineage: MicRun } | null;
}

/** A MIC model that scores zone diameters as MICs: no numbers, only why.
 *  The reason holds no revision-5 count, so it reads the same on every
 *  report revision. */
export interface MicWithheld {
  id: string;
  reason: string;
}

export type SelectiveCategory = "A1" | "A2" | "B" | "C" | "D1" | "D2" | "no_lineage";

/** Genome-only calls on the lineage split. VME and ME are among called
 *  isolates, with two-sided 95% Clopper-Pearson bounds; per-100 values count
 *  every resistant or susceptible isolate in the dataset. */
export interface SelectiveResult {
  n: number;
  n_r: number;
  n_s: number;
  tp: number;
  fn: number;
  fp: number;
  tn: number;
  coverage: number;
  coverage_r: number;
  coverage_s: number;
  vme: number;
  vme_lo: number;
  vme_hi: number;
  me: number;
  me_lo: number;
  me_hi: number;
  ca: number;
  per100r_caught: number;
  per100r_missed: number;
  per100r_lab: number;
  per100s_cleared: number;
  per100s_overcalled: number;
  per100s_lab: number;
}

interface SelectiveBase {
  id: string;
  guard: SelectiveGuard;
}

/** Only A1 and A2 carry numbers. */
export type SelectiveDataset =
  | (SelectiveBase & { category: "A1"; result: SelectiveResult })
  | (SelectiveBase & { category: "A2"; result: SelectiveResult & { learner: "gbt" | "rf" } })
  | (SelectiveBase & { category: "B" | "C" | "D1" | "D2" | "no_lineage" });

export interface Research {
  schema: number;
  snapshot: { revision: number; sha256: string; datasets: number };
  /** Every lane file read, relative to data/research. */
  sources: { path: string; sha256: string; target?: string }[];
  mic: {
    /** zone_fix holds {removed}, where a view writes the entry's removed
     *  zone diameters. */
    limits: { scope: string; baseline: string; exact_minimum: string; saved_only: string; zone_fix: string };
    datasets: MicDataset[];
    withheld: MicWithheld[];
  };
  selective: {
    /** Every category in a fixed order, with what it means in full
     *  sentences. A2's label holds {learner}, where a view writes the refit
     *  learner. */
    categories: { code: SelectiveCategory; label: string }[];
    /** bar is the pass rule every category refers to. */
    limits: { scope: string; bar: string; no_call: string; optimism: string; procedure: string; per_100: string };
    datasets: SelectiveDataset[];
  };
}
