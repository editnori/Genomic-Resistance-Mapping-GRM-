// What the index says beside the best-learner medians: each learner on both
// splits, the two error rates against the lab bars with their uncertainty,
// and the known-gene rule grouped by its own score. Every summary keeps the
// count it was computed over, so each sentence can print its basis.
//
// The bar outcomes rest on the exact (Clopper-Pearson) 95% interval of a
// rate, read as two binomial tails: the interval's upper bound is at or
// under a bar exactly when x or fewer events have probability at most 2.5%
// at the bar, and its lower bound is over the bar exactly when x or more
// events have probability under 2.5% there. The counts come from the index:
// a rate is stored over the genomes its learner called, so the events are
// the rate times that count. From evaluation revision 14 a null fit's genomes
// are scored too (null_fits_scored), so such a row's rates cover every genome.

import type { CatalogRow, Learner, LearnerLite, NoCall, Scheme } from "../schema";
import { median } from "../format";
import { REFERENCE_BARS } from "../metrics";
import { CATALOG_CUTOFF, catalogDiff, uniformCall } from "./overview";

/** The five learners that train on genomes. The known-gene rule does not. */
export const FITTED = ["scm", "cart", "rf", "gbt", "lasso"] as const satisfies readonly Learner[];
export type Fitted = (typeof FITTED)[number];

/** One side of a two-sided 95% interval. */
const TAIL = 0.025;

/** P(X <= x) for X ~ Binomial(n, p) with 0 < p < 1. Each term is carried as
 *  a logarithm, so a cohort of thousands does not underflow to zero. */
export const lowerTail = (x: number, n: number, p: number): number => {
  if (x < 0) return 0;
  if (x >= n) return 1;
  const step = Math.log(p) - Math.log1p(-p);
  let logTerm = n * Math.log1p(-p);
  let total = Math.exp(logTerm);
  for (let k = 0; k < x; k += 1) {
    logTerm += Math.log((n - k) / (k + 1)) + step;
    total += Math.exp(logTerm);
  }
  return Math.min(1, total);
};

/** The exact 95% upper bound of x events in n is at or under the bar. */
export const shownAtMost = (x: number, n: number, bar: number): boolean => n > 0 && lowerTail(x, n, bar) <= TAIL;

/** The exact 95% lower bound of x events in n is over the bar. */
export const shownAbove = (x: number, n: number, bar: number): boolean => n > 0 && x > 0 && 1 - lowerTail(x - 1, n, bar) < TAIL;

export interface RateCount {
  events: number;
  n: number;
}

export interface RateCounts {
  /** Resistant genomes called susceptible, over the resistant genomes called. */
  vme: RateCount;
  /** Susceptible genomes called resistant, over the susceptible genomes called. */
  me: RateCount;
}

const NONE: NoCall = { resistant: 0, susceptible: 0 };

/** The genomes a result's rates leave out, by class. A row written from
 *  evaluation revision 14 still counts its null-fit genomes in no_call, but
 *  it scores them (null_fits_scored), so its rates leave none out. */
export const leftOut = (result: LearnerLite | null): NoCall => (result === null || result.no_call === undefined || result.null_fits_scored === true ? NONE : result.no_call);

/** A result's two error rates as counts. Null when a rate is missing or the
 *  result called no genome of one class. A known-gene rule's Lite fits here
 *  too: it calls every genome, so it carries no no_call. */
export const rateCounts = (row: CatalogRow, result: LearnerLite | null): RateCounts | null => {
  if (result === null || result.vme === null || result.me === null) return null;
  const out = leftOut(result);
  const resistant = row.n_resistant - out.resistant;
  const susceptible = row.n - row.n_resistant - out.susceptible;
  if (resistant <= 0 || susceptible <= 0) return null;
  return {
    vme: { events: Math.round(result.vme * resistant), n: resistant },
    me: { events: Math.round(result.me * susceptible), n: susceptible },
  };
};

/** Shown to meet both bars, shown to miss one, or too few genomes to tell. */
export type BarOutcome = "meets" | "misses" | "open";

export const barOutcome = (c: RateCounts): BarOutcome => {
  if (shownAtMost(c.vme.events, c.vme.n, REFERENCE_BARS.vme) && shownAtMost(c.me.events, c.me.n, REFERENCE_BARS.me)) return "meets";
  if (shownAbove(c.vme.events, c.vme.n, REFERENCE_BARS.vme) || shownAbove(c.me.events, c.me.n, REFERENCE_BARS.me)) return "misses";
  return "open";
};

/** Whose calls a summary counts: each dataset's best learner, one learner
 *  for every dataset, or the known-gene rule. */
export type Who = "best" | Learner;

/** One result of a dataset. The known-gene rule has a single result, which
 *  answers for either split. */
export const resultOf = (row: CatalogRow, who: Who, scheme: Scheme): LearnerLite | null => {
  const learner = who === "best" ? row.best_learner : who;
  if (learner === null) return null;
  if (learner === "catalog") return row.catalog;
  return row.learners.find((l) => l.learner === learner)?.[scheme] ?? null;
};

export interface BarSummary {
  /** Datasets with both error rates for this result. */
  n: number;
  meets: number;
  misses: number;
  open: number;
  /** Of open, datasets whose two rates as printed are at or under their
   *  bars: the interval is what keeps them from "meets". */
  openAtBars: number;
}

export const barSummary = (rows: readonly CatalogRow[], who: Who, scheme: Scheme): BarSummary => {
  const s: BarSummary = { n: 0, meets: 0, misses: 0, open: 0, openAtBars: 0 };
  rows.forEach((row) => {
    const c = rateCounts(row, resultOf(row, who, scheme));
    if (c === null) return;
    const outcome = barOutcome(c);
    s.n += 1;
    s[outcome] += 1;
    if (outcome === "open" && c.vme.events / c.vme.n <= REFERENCE_BARS.vme && c.me.events / c.me.n <= REFERENCE_BARS.me) s.openAtBars += 1;
  });
  return s;
};

/** The fewest resistant genomes at which a result with no missed resistant
 *  genome is shown to meet the VME bar: (1 - bar)^n <= 2.5%. */
export const ZERO_MISS_FLOOR = Math.ceil(Math.log(TAIL) / Math.log1p(-REFERENCE_BARS.vme));

/** Datasets with at least ZERO_MISS_FLOOR resistant genomes. */
export const atZeroMissFloor = (rows: readonly CatalogRow[]): number => rows.filter((r) => r.n_resistant >= ZERO_MISS_FLOOR).length;

export interface LearnerRow {
  learner: Fitted;
  /** Datasets where this learner has a score on both splits. */
  n: number;
  medianRandom: number | null;
  medianLineage: number | null;
  /** Median of the per-dataset gaps (random minus lineage). */
  medianGap: number | null;
  /** Of n, datasets where the learner left a genome without a call on
   *  either split; its scores there cover the genomes it called. */
  nNoCall: number;
  /** Datasets with a best learner and this learner, both scored with
   *  lineages held out. */
  nBest: number;
  /** Of nBest: the best learner scores higher than this learner used alone,
   *  or lower. The rest score the same; on many the two are one learner. */
  bestHigher: number;
  bestLower: number;
  /** Median of the best learner's score minus this learner's, per dataset. */
  medianBestLead: number | null;
}

export interface LearnerSummary {
  rows: LearnerRow[];
  /** Datasets with a best learner scored with lineages held out. */
  nBest: number;
  /** Of nBest, datasets where no learner scored on every genome outscores
   *  the best learner with lineages held out. A tie counts. A learner that
   *  leaves a genome without a call there is scored on fewer genomes than
   *  the best learner, so it is not set against it. */
  bestAlsoTops: number;
  /** Of nBest, datasets where such a learner outscores the best learner by
   *  more than CATALOG_CUTOFF, the margin the rule comparison uses. */
  otherAheadByCutoff: number;
}

const balanced = (row: CatalogRow, learner: Fitted, scheme: Scheme): number | null => resultOf(row, learner, scheme)?.balanced_accuracy ?? null;

const withoutCall = (result: LearnerLite | null): boolean => leftOut(result).resistant + leftOut(result).susceptible > 0;

/** The best learner is picked by its random-split score, among the learners
 *  scored on every genome; the lineage test plays no part in the pick. These
 *  rows set it beside each learner used alone for every dataset, where no
 *  pick is made. */
export const learnerSummary = (rows: readonly CatalogRow[]): LearnerSummary => {
  const perLearner = FITTED.map((learner): LearnerRow => {
    const paired = rows.flatMap((row) => {
      const random = balanced(row, learner, "random");
      const lineage = balanced(row, learner, "lineage");
      return random === null || lineage === null ? [] : [{ row, random, lineage }];
    });
    const leads = paired.flatMap((p) => (p.row.best_lineage === null ? [] : [p.row.best_lineage - p.lineage]));
    return {
      learner,
      n: paired.length,
      medianRandom: median(paired.map((p) => p.random)),
      medianLineage: median(paired.map((p) => p.lineage)),
      medianGap: median(paired.map((p) => p.random - p.lineage)),
      nNoCall: paired.filter((p) => withoutCall(resultOf(p.row, learner, "random")) || withoutCall(resultOf(p.row, learner, "lineage"))).length,
      nBest: leads.length,
      bestHigher: leads.filter((d) => d > 0).length,
      bestLower: leads.filter((d) => d < 0).length,
      medianBestLead: median(leads),
    };
  });
  // How far the top learner scored on every genome sits over the best
  // learner with lineages held out; zero when none does.
  const margins = rows.flatMap((row) => {
    if (row.best_lineage === null) return [];
    const scores = FITTED.flatMap((learner) => {
      const result = resultOf(row, learner, "lineage");
      const score = result?.balanced_accuracy ?? null;
      return score === null || withoutCall(result) ? [] : [score];
    });
    return [Math.max(row.best_lineage, ...scores) - row.best_lineage];
  });
  return {
    rows: perLearner,
    nBest: margins.length,
    bestAlsoTops: margins.filter((d) => d <= 0).length,
    otherAheadByCutoff: margins.filter((d) => d > CATALOG_CUTOFF).length,
  };
};

/** The rule's own balanced accuracy in four groups. A rule that calls every
 *  genome one way scores exactly 50%, a coin toss, so the first cut sits just
 *  above it at 55%. The groups above it are even 15-point steps. */
export const RULE_BAND_EDGES = [0.55, 0.7, 0.85] as const;

export interface RuleBand {
  /** The group covers from <= score < to. */
  from: number;
  to: number;
  /** Datasets whose rule scores in this group. */
  datasets: number;
  /** Of those, datasets that also ran the lineage split. */
  compared: number;
  /** The compared datasets by rule minus best learner, at CATALOG_CUTOFF:
   *  the rule ahead of the learner or within the cutoff, or behind by more.
   *  They sum to compared. */
  ruleHolds: number;
  ruleBehind: number;
  /** Median of best learner (lineage held out) minus rule, per dataset. */
  medianLearnerLead: number | null;
  /** Datasets whose rule calls at least 90% of one class the other way
   *  (uniformCall), by the call it makes. */
  uniform: { resistant: number; susceptible: number };
}

export const ruleBands = (rows: readonly CatalogRow[]): RuleBand[] => {
  const edges = [-Infinity, ...RULE_BAND_EDGES, Infinity];
  return edges.slice(0, -1).map((from, i): RuleBand => {
    const to = edges[i + 1] ?? Infinity;
    const group = rows.filter((r) => r.catalog !== null && r.catalog.balanced_accuracy >= from && r.catalog.balanced_accuracy < to);
    const diffs = group.flatMap((r) => {
      const d = catalogDiff(r);
      return d === null ? [] : [d];
    });
    const calls = group.map((r) => (r.catalog === null ? null : uniformCall(r.catalog)));
    return {
      from,
      to,
      datasets: group.length,
      compared: diffs.length,
      ruleHolds: diffs.filter((d) => d >= -CATALOG_CUTOFF).length,
      ruleBehind: diffs.filter((d) => d < -CATALOG_CUTOFF).length,
      medianLearnerLead: median(diffs.map((d) => -d)),
      uniform: { resistant: calls.filter((c) => c === "resistant").length, susceptible: calls.filter((c) => c === "susceptible").length },
    };
  });
};
