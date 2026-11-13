// The input comparison's second reading across four species
// (docs/input-comparison-runner.md, "A second reading across species"). Each
// changed input is scored against the default on the same datasets, by log
// loss with whole lineages held out. The counts are bundled in
// data/input-comparison.json, which scripts/input-comparison/site_data.py
// writes from reference/input-comparison/second-reading.json; nothing here
// tests or decides anything again.

import bundled from "../input-comparison.json";

/** Datasets where the change's mean log-loss change over the learners is above zero (better), below it (worse) or exactly zero. */
export interface DatasetCounts {
  better: number;
  worse: number;
  unchanged: number;
  /** Datasets with at least one comparison made. */
  datasets: number;
}

export interface InputArm {
  id: string;
  /** The default run again: a measure of how far a score moves when no input changes. */
  twin: boolean;
  /** Over the four species, with the two-sided exact sign test's p and the median change. */
  all: DatasetCounts & { p: number; median: number };
  /** In the order of `species`. */
  species: readonly DatasetCounts[];
  /** Comparisons per dataset and learner over the four species: made, within the margin, at least the margin either way, exactly zero. */
  contrasts: { made: number; inside: number; beyond: number; zero: number };
}

export interface InputComparison {
  source: string;
  margin: number;
  species: readonly { name: string; datasets: number; answered: number }[];
  family: { size: number; tested: number };
  firm: { better: number; worse: number; learners: Record<string, number> };
  predictions: {
    species: readonly string[];
    tests: readonly { arm: string; predicted: "better" | "worse"; better: number; worse: number; unchanged: number; p_one_sided: number }[];
  };
  arms: readonly InputArm[];
}

export const INPUT_COMPARISON = bundled as InputComparison;

const sum = (xs: readonly number[]): number => xs.reduce((a, b) => a + b, 0);

export const armById = (id: string, c: InputComparison = INPUT_COMPARISON): InputArm => {
  const a = c.arms.find((x) => x.id === id);
  if (a === undefined) throw new Error(`no ${id} in the input comparison`);
  return a;
};

export interface InputSummary {
  species: number;
  datasets: number;
  answered: number;
  changed: readonly InputArm[];
  twins: readonly InputArm[];
  firm: number;
  firmBetter: number;
  firmWorse: number;
  firmForest: number;
  tested: number;
}

export const inputSummary = (c: InputComparison = INPUT_COMPARISON): InputSummary => ({
  species: c.species.length,
  datasets: sum(c.species.map((s) => s.datasets)),
  answered: sum(c.species.map((s) => s.answered)),
  changed: c.arms.filter((a) => !a.twin),
  twins: c.arms.filter((a) => a.twin),
  firm: c.firm.better + c.firm.worse,
  firmBetter: c.firm.better,
  firmWorse: c.firm.worse,
  firmForest: c.firm.learners.rf ?? 0,
  tested: c.family.tested,
});

/** A species' place in `species`, by its name. */
export const speciesIndex = (name: string, c: InputComparison = INPUT_COMPARISON): number => {
  const i = c.species.findIndex((s) => s.name === name);
  if (i < 0) throw new Error(`no ${name} in the input comparison`);
  return i;
};

/** The prediction's test for one change, on the species it named. */
export const prediction = (arm: string, c: InputComparison = INPUT_COMPARISON) => {
  const t = c.predictions.tests.find((x) => x.arm === arm);
  if (t === undefined) throw new Error(`no prediction for ${arm}`);
  return { ...t, of: t.better + t.worse };
};
