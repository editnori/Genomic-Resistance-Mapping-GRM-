// Bundled summaries of the five additional revision 14 studies:
// first-public-day (time) split, Kover 2.0.7 comparison, dense neural
// network, AMRFinderPlus 4.2.7 program comparison, and TBProfiler 6.6.3
// comparison on M. tuberculosis assemblies. Generated from the committed
// reference/ TSVs by scripts/release-summary/studies_site_data.py.

import raw from "../additional-studies.json";

export type ComparisonLearner = "scm" | "cart" | "rf" | "gbt" | "lasso";
export const COMPARISON_LEARNERS: readonly ComparisonLearner[] = ["scm", "cart", "rf", "gbt", "lasso"];

export interface TimeMeasure {
  median: number;
  low: number;
  high: number;
  reading: string;
  worse: number;
  better: number;
}

export interface DenseCell {
  diff: number;
  low: number;
  high: number;
  positive: number;
  negative: number;
}

export interface AmrLearnerRow {
  ba_median: number;
  ba_lo: number;
  ba_hi: number;
  ba_ahead: number;
  ba_level: number;
  ba_behind: number;
  sens_median: number;
  spec_median: number;
}

export interface TbLearnerRow {
  tool_higher: number;
  tool_lower: number;
  equal: number;
  median_learner_minus_tool: number;
}

export interface AdditionalStudiesData {
  source: string;
  evaluation_revision: number;
  time_split: {
    datasets: number;
    species: number;
    tested_pairs: number;
    resistant_pairs: number;
    susceptible_pairs: number;
    answer: string;
    learners_with_cost: number;
    total_log_loss_worse: number;
    total_balanced_accuracy_worse: number;
    total_recall_worse: number;
    total_specificity_worse: number;
    by_learner: Record<ComparisonLearner, { recall: TimeMeasure; log_loss: TimeMeasure; specificity: TimeMeasure; balanced_accuracy: TimeMeasure }>;
  };
  kover: {
    datasets: number;
    species: number;
    genomes: number;
    folds: number;
    kover_one_rule_folds: number;
    balanced_accuracy: { value: number; low: number; high: number; higher: number; equal: number; lower: number };
    recall: { value: number; low: number; high: number };
    specificity: { value: number; low: number; high: number };
    log_loss: { value: number; low: number; high: number };
    ties_equal_ba: { value: number; low: number; high: number };
  };
  dense_network: {
    datasets: number;
    species: number;
    outcome: string;
    vs_best5: { balanced_accuracy: DenseCell; recall: DenseCell; log_loss: DenseCell };
    by_learner: Record<ComparisonLearner, { balanced_accuracy: DenseCell; recall: DenseCell; log_loss: DenseCell }>;
  };
  amrfinderplus: {
    datasets: number;
    species: number;
    program_vs_known_gene_rule: {
      ahead: number;
      level: number;
      behind: number;
      median: number;
      mean: number;
      ahead_by_10pts: number;
      ahead_by_10pts_quinolone: number;
    };
    by_learner_lineage: Record<ComparisonLearner, AmrLearnerRow>;
  };
  tbprofiler: {
    by_name_endpoints: number;
    by_name_rows: number;
    by_name_rows_with_call: number;
    by_name_rows_without_call: number;
    graded_endpoints: number;
    graded_rows: number;
    graded_rows_with_call: number;
    graded_rows_without_call: number;
    by_name_sensitivity_median: number;
    by_name_specificity_median: number;
    graded_sensitivity_median: number;
    graded_specificity_median: number;
    graded_vs_learners_lineage: Record<ComparisonLearner, TbLearnerRow>;
  };
}

export const ADDITIONAL_STUDIES = raw as unknown as AdditionalStudiesData;
