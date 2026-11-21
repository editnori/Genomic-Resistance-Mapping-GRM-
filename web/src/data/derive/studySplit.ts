// The retrospective split by study inside BV-BRC (docs/study-split.md):
// 70 datasets across 9 species, scored with folds keeping each BioProject whole
// beside random folds and matched random folds. Bundled in data/study-split.json
// by scripts/study-split/site_data.py from reference/study-split/.

import bundled from "../study-split.json";

export const STUDY_LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"] as const;
export type StudyLearner = (typeof STUDY_LEARNERS)[number];
export type StudyFamily = "SR" | "SM" | "MR";

export interface StudyContrast {
  median: number;
  low: number;
  high: number;
  reading: "worse" | "better" | "not shown";
  datasets: number;
  worse: number;
  equal: number;
  better: number;
}

export interface StudyLineageContrast extends StudyContrast {
  keep: number;
}

export interface FirstEntry {
  recall: StudyContrast;
  specificity: StudyContrast;
  balanced_accuracy: StudyContrast;
  log_loss: StudyContrast;
  error_rate?: StudyContrast;
  rule: string;
  keep: number;
}

export interface SecondFamilyEntry {
  recall: StudyContrast;
  specificity: StudyContrast;
  loss: StudyContrast;
  loss_measure: "error_rate" | "log_loss";
  pairs: number;
  cost: number;
  recall_worse: number;
  loss_worse: number;
  answer: string;
}

export interface StudySplitData {
  source: string;
  evaluation_revision: number;
  datasets: number;
  species: number;
  pairs: number;
  resistant_pairs: number;
  no_study_pairs: number;
  distinct_genomes: number;
  distinct_studies: number;
  lineage_datasets: number;
  species_counts: readonly { slug: string; name: string; datasets: number }[];
  first: Record<StudyLearner, FirstEntry>;
  second: Record<StudyLearner, Record<StudyFamily, SecondFamilyEntry>>;
  vs_lineage: Record<StudyLearner, { recall: StudyLineageContrast; specificity: StudyLineageContrast; log_loss: StudyLineageContrast }>;
}

export const STUDY_SPLIT = bundled as StudySplitData;

export interface StudySplitSummary {
  datasets: number;
  species: number;
  pairs: number;
  costShown: number;
  logLossWorse: number;
  specificityWorse: number;
  recallWorseLearners: readonly StudyLearner[];
  recallNotShownLearners: readonly StudyLearner[];
  lineageLogLossNotShown: number;
}

export const studySplitSummary = (d: StudySplitData = STUDY_SPLIT): StudySplitSummary => ({
  datasets: d.datasets,
  species: d.species,
  pairs: d.pairs,
  costShown: STUDY_LEARNERS.filter((l) => d.first[l].rule === "cost").length,
  logLossWorse: STUDY_LEARNERS.filter((l) => d.first[l].log_loss.reading === "worse").length,
  specificityWorse: STUDY_LEARNERS.filter((l) => d.first[l].specificity.reading === "worse").length,
  recallWorseLearners: STUDY_LEARNERS.filter((l) => d.first[l].recall.reading === "worse"),
  recallNotShownLearners: STUDY_LEARNERS.filter((l) => d.first[l].recall.reading !== "worse"),
  lineageLogLossNotShown: STUDY_LEARNERS.filter((l) => d.vs_lineage[l].log_loss.reading === "not shown").length,
});
