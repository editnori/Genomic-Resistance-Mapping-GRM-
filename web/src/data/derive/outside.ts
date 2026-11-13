// The outside-study test (docs/outside-study-results.md): each report's final
// model run, unchanged, on newly deposited genomes. The rows are bundled in
// data/outside-study.json, which scripts/outside-study-results/site_data.py
// writes from reference/outside-study with the documents' own reading rules,
// and every count here is made from those rows. The two populations are
// summarised apart and never added together or differenced.

import bundled from "../outside-study.json";

export type Learner = "scm" | "cart" | "lasso" | "known_gene";
export type Measure = "sensitivity" | "specificity";

/** The report's lineage-split figure, the new estimate and its interval. */
export interface Figure {
  report: number | null;
  estimate: number;
  lo: number;
  hi: number;
}

export interface OutsideRow {
  dataset: string;
  learner: Learner;
  /** "below", "consistent" or "above"; anything else is a row the plan leaves unread. */
  sensitivity: string;
  specificity: string;
  /** The set covering machine's rows only. */
  figures?: Record<Measure, Figure>;
  /** The known-gene rule's rows only: its report figure is exactly 0 or 1 on a measure. */
  report_at_edge?: boolean;
}

export interface Population {
  id: "primary" | "same-study";
  species: readonly { slug: string; frozen: number; scored: number }[];
  datasets: readonly {
    id: string;
    species: string;
    antibiotic: string;
    estimates: boolean;
    top_study_resistant_share: number;
    /** The study that holds the most resistant genomes. */
    top_study_resistant: string;
    /** Genomes scored resistant or susceptible. */
    scored: number;
  }[];
  rows: readonly OutsideRow[];
}

export interface OutsideStudy {
  source: string;
  evaluation_revision: number;
  confidence: number;
  populations: readonly Population[];
}

export const OUTSIDE_STUDY = bundled as OutsideStudy;

/** The learners in the documents' order. */
export const OUTSIDE_LEARNERS: readonly Learner[] = ["scm", "cart", "lasso", "known_gene"];

const READ = new Set(["below", "consistent", "above"]);

/** Whether a row holds: neither reading is "below". Null when either reading
 *  is one the plan leaves unread (tables.py, holds). */
export const holds = (r: OutsideRow): boolean | null =>
  READ.has(r.sensitivity) && READ.has(r.specificity) ? r.sensitivity !== "below" && r.specificity !== "below" : null;

export interface PopulationSummary {
  id: Population["id"];
  frozen: number;
  scored: number;
  species: number;
  datasets: number;
  withEstimates: number;
  /** Datasets with estimates where one study holds at least half of the resistant genomes. */
  oneStudyHalf: number;
  /** Datasets with estimates where one study holds at least ONE_STUDY_ALL of the resistant genomes. */
  oneStudyAll: number;
  /** The distinct studies that hold the most resistant genomes on those datasets. */
  oneStudyAllStudies: number;
  /** The median, over the datasets with estimates, of the genomes scored resistant or susceptible. */
  medianScored: number;
  /** Datasets where the known-gene rule holds and its report figure is exactly 0% or 100% on a measure. */
  ruleHoldsAtEdge: number;
  /** Per learner: datasets where it holds, of the datasets read for it. */
  holds: Record<Learner, { holds: number; read: number }>;
  /** The set covering machine's readings over the rows read: each measure
   *  above and below the report's figure, and the two together. */
  scm: { sensitivityAbove: number; sensitivityBelow: number; specificityAbove: number; specificityBelow: number; onlySpecificity: number; onlySensitivity: number; both: number };
}

/** The share above which one study holds a dataset's resistant genomes all but wholly. */
export const ONE_STUDY_ALL = 0.999;

const median = (xs: readonly number[]): number => {
  if (xs.length === 0) return Number.NaN;
  const s = [...xs].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 === 1 ? s[m]! : (s[m - 1]! + s[m]!) / 2;
};

const count = <T,>(xs: readonly T[], f: (x: T) => boolean): number => xs.reduce((n, x) => n + (f(x) ? 1 : 0), 0);

export const populationSummary = (p: Population): PopulationSummary => {
  const learner = (l: Learner) => p.rows.filter((r) => r.learner === l);
  const tally = Object.fromEntries(
    OUTSIDE_LEARNERS.map((l) => {
      const read = learner(l).filter((r) => holds(r) !== null);
      return [l, { holds: count(read, (r) => holds(r) === true), read: read.length }];
    }),
  ) as Record<Learner, { holds: number; read: number }>;
  const scm = learner("scm");
  const read = scm.filter((r) => holds(r) !== null);
  const est = p.datasets.filter((d) => d.estimates);
  return {
    id: p.id,
    frozen: p.species.reduce((a, s) => a + s.frozen, 0),
    scored: p.species.reduce((a, s) => a + s.scored, 0),
    species: p.species.length,
    datasets: p.datasets.length,
    withEstimates: est.length,
    oneStudyHalf: count(est, (d) => d.top_study_resistant_share >= 0.5),
    oneStudyAll: count(est, (d) => d.top_study_resistant_share >= ONE_STUDY_ALL),
    oneStudyAllStudies: new Set(est.filter((d) => d.top_study_resistant_share >= ONE_STUDY_ALL).map((d) => d.top_study_resistant)).size,
    medianScored: median(est.map((d) => d.scored)),
    ruleHoldsAtEdge: count(learner("known_gene"), (r) => holds(r) === true && r.report_at_edge === true),
    holds: tally,
    scm: {
      sensitivityAbove: count(read, (r) => r.sensitivity === "above"),
      sensitivityBelow: count(read, (r) => r.sensitivity === "below"),
      specificityAbove: count(read, (r) => r.specificity === "above"),
      specificityBelow: count(read, (r) => r.specificity === "below"),
      onlySpecificity: count(read, (r) => r.specificity === "below" && r.sensitivity !== "below"),
      onlySensitivity: count(read, (r) => r.sensitivity === "below" && r.specificity !== "below"),
      both: count(read, (r) => r.sensitivity === "below" && r.specificity === "below"),
    },
  };
};

/** One drawn mark: a set covering machine row read on one measure. */
export interface OutsidePoint {
  dataset: string;
  measure: Measure;
  report: number;
  figure: Figure;
  reading: string;
}

/** The set covering machine's rows that can be read on a measure, as marks. */
export const outsidePoints = (p: Population, measure: Measure): OutsidePoint[] =>
  p.rows.flatMap((r) => {
    const f = r.figures?.[measure];
    const reading = r[measure];
    return r.learner === "scm" && f !== undefined && f.report !== null && READ.has(reading) ? [{ dataset: r.dataset, measure, report: f.report, figure: f, reading }] : [];
  });

export const populationById = (id: Population["id"]): Population => {
  const p = OUTSIDE_STUDY.populations.find((x) => x.id === id);
  if (p === undefined) throw new Error(`no ${id} population in the outside-study rows`);
  return p;
};
