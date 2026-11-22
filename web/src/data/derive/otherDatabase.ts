// Test on data from another public database, NCBI Pathogen Detection
// (docs/other-database-test.md): 13,377 genomes across 9 species (E. coli pilot
// and 8 more species). Bundled in data/other-database.json by
// scripts/other-database/site_data.py from reference/other-database/.

import bundled from "../other-database.json";

export const OTHER_DB_LEARNERS = ["scm", "cart", "lasso"] as const;
export type OtherDbLearner = (typeof OTHER_DB_LEARNERS)[number];

export interface OtherDbSpecies {
  slug: string;
  name: string;
  part: "pilot" | "eight";
  scored: number;
  reports_scored: number;
  shared_bioproject: number;
  largest_bioproject: string;
  largest_bioproject_share: number;
  largest_bioproject_training: number;
}

export interface OtherDbTotals {
  species: number;
  plan_meets: number;
  plan_fails: number;
  plan_not_assessable: number;
  estimable: number;
  holds: number;
  recall_below: number;
  log_loss_worse: number;
  specificity_under_half: number;
}

export interface OtherDbSpeciesLearner {
  species: string;
  name: string;
  part: "pilot" | "eight";
  learner: OtherDbLearner;
  estimable: number;
  holds: number;
  recall_below: number;
  log_loss_worse: number;
  specificity_under_half: number;
  lowest_specificity: number | null;
  lowest_specificity_dataset: string | null;
  rule: "meets" | "fails" | "not assessable";
  margin: number | null;
}

export interface OtherDatabaseData {
  source: string;
  evaluation_revision: number;
  species: readonly OtherDbSpecies[];
  totals: Record<OtherDbLearner, OtherDbTotals>;
  by_species: readonly OtherDbSpeciesLearner[];
}

export const OTHER_DATABASE = bundled as OtherDatabaseData;

export interface OtherDatabaseSummary {
  totalGenomes: number;
  totalSpecies: number;
  pilotGenomes: number;
  eightGenomes: number;
  eightSpecies: number;
  totals: Record<OtherDbLearner, OtherDbTotals & { assessable: number }>;
  pilot: Record<OtherDbLearner, OtherDbSpeciesLearner>;
}

export const otherDatabaseSummary = (d: OtherDatabaseData = OTHER_DATABASE): OtherDatabaseSummary => {
  const pilotSp = d.species.filter((s) => s.part === "pilot");
  const eightSp = d.species.filter((s) => s.part === "eight");
  const totals = Object.fromEntries(
    OTHER_DB_LEARNERS.map((l) => {
      const t = d.totals[l];
      return [l, { ...t, assessable: t.species - t.plan_not_assessable }];
    }),
  ) as Record<OtherDbLearner, OtherDbTotals & { assessable: number }>;
  const pilot = Object.fromEntries(
    OTHER_DB_LEARNERS.map((l) => [l, d.by_species.find((r) => r.part === "pilot" && r.learner === l)!]),
  ) as Record<OtherDbLearner, OtherDbSpeciesLearner>;
  return {
    totalGenomes: d.species.reduce((a, s) => a + s.scored, 0),
    totalSpecies: d.species.length,
    pilotGenomes: pilotSp.reduce((a, s) => a + s.scored, 0),
    eightGenomes: eightSp.reduce((a, s) => a + s.scored, 0),
    eightSpecies: eightSp.length,
    totals,
    pilot,
  };
};
