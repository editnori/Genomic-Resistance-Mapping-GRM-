// In-memory research.json pieces for the tests. Each guard is the one
// reportGuard gives the report the entry is built from, so the entry applies
// to that report until a test changes either. Nothing here is a measurement:
// every limit and label is a placeholder a test can find in a view.

import { reportGuard } from "../../data/research";
import {
  RESEARCH_SCHEMA,
  REVISION5_SNAPSHOT_SHA256,
  type DatasetReport,
  type MicDataset,
  type MicRun,
  type MicWithheld,
  type Research,
  type SelectiveCategory,
  type SelectiveDataset,
  type SelectiveResult,
} from "../../data/schema";

export const SELECTIVE_CODES: readonly SelectiveCategory[] = ["A1", "A2", "B", "C", "D1", "D2", "no_lineage"];

/** One prediction set, with enough exact MICs to need no exact-minimum limit. */
export const micRun = (ea: number | null, n_exact = 60): MicRun => ({ n_predicted: 80, n_exact, ea_exact: ea });

/** Genome-only calls for 32 of 40 resistant and 54 of 60 susceptible isolates. */
export const SELECTIVE_RESULT: SelectiveResult = {
  n: 100, n_r: 40, n_s: 60, tp: 30, fn: 2, fp: 4, tn: 50,
  coverage: 0.86, coverage_r: 0.8, coverage_s: 0.9, vme: 2 / 32, vme_lo: 0.0077, vme_hi: 0.2081, me: 4 / 54, me_lo: 0.0212, me_hi: 0.1789, ca: 80 / 86,
  per100r_caught: 75, per100r_missed: 5, per100r_lab: 20, per100s_cleared: 250 / 3, per100s_overcalled: 20 / 3, per100s_lab: 10,
};

/** The reason src/research.rs writes for a withheld saved model. */
export const WITHHELD_REASON = "The revision-5 saved MIC model scored disk-diffusion zone diameters as MICs, and there is no rerun or no-genome baseline without them.";

export const MIC_LIMITS: Research["mic"]["limits"] = {
  scope: "MIC scope limit.",
  baseline: "MIC baseline limit.",
  exact_minimum: "MIC exact-minimum limit.",
  saved_only: "MIC saved-only limit.",
  zone_fix: "MIC zone-fix limit for {removed}.",
};

export const SELECTIVE_LIMITS: Research["selective"]["limits"] = {
  scope: "Selective scope limit.",
  bar: "Selective bar.",
  no_call: "Selective no-call limit.",
  optimism: "Selective optimism limit.",
  procedure: "Selective procedure limit.",
  per_100: "Selective per-100 limit.",
};

/** A placeholder sentence per category; A2's names its learner. */
export const categoryLabel = (code: SelectiveCategory): string => (code === "A2" ? "Category A2 with {learner}." : `Category ${code}.`);

/** A MIC entry guarding `r`: a rerun with both splits and a baseline, or a
 *  saved report model with neither. */
export const micEntry = (r: DatasetReport, saved = false): MicDataset => {
  const { scores: _scores, ...guard } = reportGuard(r);
  return {
    id: r.id,
    guard,
    zone_diameters: { removed: 0, scored: 0 },
    model: saved
      ? { source: "saved_report", random: micRun(0.8), lineage: null }
      : { source: "rerun", random: micRun(0.8), lineage: { ...micRun(0.6), holdout_random_folds: 0 } },
    baseline: saved ? null : { random: micRun(0.7), lineage: micRun(0.5) },
  };
};

export const selectiveEntry = (r: DatasetReport, category: SelectiveCategory): SelectiveDataset => {
  const { samples, sample_ids, labels, lineages, scores } = reportGuard(r);
  const base = { id: r.id, guard: { samples, sample_ids, labels, lineages, scores } };
  if (category === "A1") return { ...base, category, result: SELECTIVE_RESULT };
  if (category === "A2") return { ...base, category, result: { ...SELECTIVE_RESULT, learner: "gbt" } };
  return { ...base, category };
};

export const researchWith = (mic: MicDataset[], selective: SelectiveDataset[], withheld: MicWithheld[] = []): Research => ({
  schema: RESEARCH_SCHEMA,
  snapshot: { revision: 5, sha256: REVISION5_SNAPSHOT_SHA256, datasets: 234 },
  sources: [],
  mic: { limits: MIC_LIMITS, datasets: mic, withheld },
  selective: {
    categories: SELECTIVE_CODES.map((code) => ({ code, label: categoryLabel(code) })),
    limits: SELECTIVE_LIMITS,
    datasets: selective,
  },
});
