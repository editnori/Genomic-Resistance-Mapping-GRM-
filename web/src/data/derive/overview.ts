// Corpus-level findings computed from the manifest. Every summary keeps the
// count it was computed over, so each sentence can print its basis.

import type { CatalogRow, Lite, Manifest } from "../schema";
import { int, listWords, median, plural } from "../format";

/** " 3 datasets with a catalog rule have no learner that calls every isolate
 *  and are left out." Empty when every dataset has a best learner. */
export const noBestWords = (n: number, which = ""): string =>
  n === 0 ? "" : ` ${plural(n, "dataset")}${which} ${n === 1 ? "has" : "have"} no learner that calls every isolate and ${n === 1 ? "is" : "are"} left out.`;

export interface SplitSummary {
  /** Datasets with a best learner, so a random-split score. */
  n: number;
  /** Datasets where no learner calls every isolate, left out of every median. */
  nNoBest: number;
  /** Datasets that also ran the lineage split. */
  nLineage: number;
  medianRandom: number | null;
  medianLineage: number | null;
  /** Median of the per-dataset gaps (random minus lineage), same learner. */
  medianGap: number | null;
  /** Datasets losing more than 5 points when lineages are held out. */
  nLargeGap: number;
}

export const splitSummary = (rows: readonly CatalogRow[]): SplitSummary => {
  const eligible = rows.filter((r) => r.best_random !== null);
  const paired = eligible.filter((r) => r.best_lineage !== null);
  return {
    n: eligible.length,
    nNoBest: rows.length - eligible.length,
    nLineage: paired.length,
    medianRandom: median(paired.map((r) => r.best_random!)),
    medianLineage: median(paired.map((r) => r.best_lineage!)),
    medianGap: median(paired.map((r) => r.best_random! - r.best_lineage!)),
    nLargeGap: paired.filter((r) => r.best_random! - r.best_lineage! > 0.05).length,
  };
};

/** Paired datasets with the largest per-dataset gaps (random minus lineage),
 *  largest first: the figure's direct labels. */
export const largestGaps = (rows: readonly CatalogRow[], k: number): { row: CatalogRow; gap: number }[] =>
  rows
    .flatMap((row) => (row.best_random === null || row.best_lineage === null ? [] : [{ row, gap: row.best_random - row.best_lineage }]))
    .sort((a, b) => b.gap - a.gap || a.row.id.localeCompare(b.row.id))
    .slice(0, k);

/** The catalog rule is fixed in advance, so it is compared with the best
 *  learner's score with lineages held out. */
export const catalogDiff = (r: CatalogRow): number | null =>
  r.catalog === null || r.best_lineage === null ? null : r.catalog.balanced_accuracy - r.best_lineage;

export const CATALOG_CUTOFFS = [0.01, 0.02, 0.05] as const;
/** The cutoff the sentences use; the other two print beside it. */
export const CATALOG_CUTOFF = 0.02;

/** A catalog rule that calls at least 90% of one class the other way scores
 *  near 50% balanced accuracy. Specificity at most 10% means it calls at
 *  least 90% of susceptible isolates resistant; sensitivity at most 10%, at
 *  least 90% of resistant isolates susceptible. Read from the rule's own
 *  sensitivity and specificity; it changes no rule. */
export const UNIFORM_SHARE = 0.1;

export type UniformCall = "resistant" | "susceptible";

export const uniformCall = (l: Lite): UniformCall | null =>
  l.specificity <= UNIFORM_SHARE ? "resistant" : l.sensitivity <= UNIFORM_SHARE ? "susceptible" : null;

export interface CatalogSummary {
  /** Datasets whose drug has an acquired-gene rule. */
  nWithRule: number;
  /** Of those, datasets that also ran the lineage split. */
  nCompared: number;
  /** Datasets where the catalog is within each cutoff of the best learner,
   *  or ahead of it. */
  within: Map<number, number>;
  medianDiff: number | null;
  /** The compared datasets by catalog minus learner, at CATALOG_CUTOFF: ahead
   *  by more, within it either way, behind by more. They sum to nCompared. */
  regions: { ahead: number; within: number; behind: number };
  /** The compared datasets by uniformCall; they sum to nCompared. */
  uniform: { resistant: number; susceptible: number; other: number };
}

export const catalogSummary = (rows: readonly CatalogRow[]): CatalogSummary => {
  const compared = rows.flatMap((r) => {
    const d = catalogDiff(r);
    return d === null || r.catalog === null ? [] : [{ d, call: uniformCall(r.catalog) }];
  });
  const diffs = compared.map((c) => c.d);
  return {
    nWithRule: rows.filter((r) => r.catalog !== null).length,
    nCompared: diffs.length,
    within: new Map(CATALOG_CUTOFFS.map((c) => [c, diffs.filter((d) => d >= -c).length])),
    medianDiff: median(diffs),
    regions: {
      ahead: diffs.filter((d) => d > CATALOG_CUTOFF).length,
      within: diffs.filter((d) => d >= -CATALOG_CUTOFF && d <= CATALOG_CUTOFF).length,
      behind: diffs.filter((d) => d < -CATALOG_CUTOFF).length,
    },
    uniform: {
      resistant: compared.filter((c) => c.call === "resistant").length,
      susceptible: compared.filter((c) => c.call === "susceptible").length,
      other: compared.filter((c) => c.call === null).length,
    },
  };
};

export type ClassVerdict = "catalog" | "beyond" | "no-rule";

export interface ClassRow {
  drugClass: string;
  datasets: number;
  /** Datasets compared (catalog rule and lineage split both present). */
  compared: number;
  medianCatalog: number | null;
  medianLearner: number | null;
  /** Median of best learner (lineage) minus catalog, per dataset. */
  medianGap: number | null;
  verdict: ClassVerdict;
}

/** Per drug class: does the acquired-gene rule match the learners, or do the
 *  learners find resistance the rule does not test? */
export const classRows = (rows: readonly CatalogRow[]): ClassRow[] => {
  const by = new Map<string, CatalogRow[]>();
  rows.forEach((r) => by.set(r.drug_class, [...(by.get(r.drug_class) ?? []), r]));
  return [...by.entries()]
    .map(([drugClass, group]): ClassRow => {
      const compared = group.filter((r) => catalogDiff(r) !== null);
      const medianGap = median(compared.map((r) => -catalogDiff(r)!));
      return {
        drugClass,
        datasets: group.length,
        compared: compared.length,
        medianCatalog: median(compared.map((r) => r.catalog!.balanced_accuracy)),
        medianLearner: median(compared.map((r) => r.best_lineage!)),
        medianGap,
        verdict: medianGap === null ? "no-rule" : medianGap <= CATALOG_CUTOFF ? "catalog" : "beyond",
      };
    })
    .sort((a, b) => (a.medianGap ?? Infinity) - (b.medianGap ?? Infinity) || a.drugClass.localeCompare(b.drugClass));
};

export const speciesCount = (m: Manifest): number => new Set(m.datasets.map((d) => d.species)).size;

/** A report's class minimum. `grm learn` gated every dataset at 25 before
 *  reports recorded one, so a report or manifest row without it reads as 25,
 *  as report.rs writes it. */
export const minClassOf = (recorded: number | undefined): number => recorded ?? 25;

/** The class minimum of the analyzed reports: one value for every report,
 *  each value with its report count when they differ, or none before any
 *  dataset has a report. */
export type ClassMinimum =
  | { kind: "value"; value: number }
  | { kind: "mixed"; reports: [minimum: number, reports: number][] }
  | { kind: "unavailable" };

export const classMinimum = (rows: readonly CatalogRow[]): ClassMinimum => {
  const by = new Map<number, number>();
  rows.forEach((r) => {
    const minimum = minClassOf(r.min_class);
    by.set(minimum, (by.get(minimum) ?? 0) + 1);
  });
  const reports = [...by.entries()].sort((a, b) => a[0] - b[0]);
  if (reports.length === 0) return { kind: "unavailable" };
  return reports.length === 1 ? { kind: "value", value: reports[0]![0] } : { kind: "mixed", reports };
};

/** What makes a species and antibiotic a dataset, completing "A dataset is
 *  one species and one antibiotic with ...". */
export const datasetRule = (c: ClassMinimum): string => {
  const genomes = "genomes that passed QC and have a reported lab call";
  if (c.kind === "value") return `at least ${int(c.value)} resistant and ${int(c.value)} susceptible ${genomes}`;
  if (c.kind === "unavailable") return `enough resistant and susceptible ${genomes}; no dataset report records the minimum yet`;
  const each = listWords(c.reports.map(([minimum, n]) => `at least ${int(minimum)} of each in ${plural(n, "report")}`));
  return `enough resistant and susceptible ${genomes}: ${each}`;
};

export interface MicSummary {
  n: number;
  median: number | null;
  /** Datasets where at least 90% of exact MICs are predicted within one dilution. */
  nAbove90: number;
}

export const micSummary = (rows: readonly CatalogRow[]): MicSummary => {
  const ea = rows.map((r) => r.mic_essential_agreement).filter((v): v is number => v !== null);
  return { n: ea.length, median: median(ea), nAbove90: ea.filter((v) => v >= 0.9).length };
};
