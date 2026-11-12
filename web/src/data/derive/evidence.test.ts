import { describe, expect, it } from "vitest";
import type { CatalogRow, Learner, LearnerLite, Lite, NoCall } from "../schema";
import { FITTED, ZERO_MISS_FLOOR, atZeroMissFloor, barOutcome, barSummary, learnerSummary, leftOut, lowerTail, rateCounts, resultOf, ruleBands, shownAbove, shownAtMost } from "./evidence";

/** A learner's result on one split: balanced accuracy and the two error rates. */
const fit = (ba: number | null, vme: number | null = 0, me: number | null = 0, noCall?: NoCall): LearnerLite => ({
  balanced_accuracy: ba,
  sensitivity: vme === null ? null : 1 - vme,
  specificity: me === null ? null : 1 - me,
  vme,
  me,
  auroc: null,
  ...(noCall === undefined ? {} : { no_call: noCall }),
});

/** The same result as a revision-14 row writes it: the null-fit genomes are
 *  still counted in no_call, and they are scored. */
const scored = (result: LearnerLite): LearnerLite => ({ ...result, null_fits_scored: true });

const rule = (ba: number, sensitivity: number, specificity: number): Lite => ({ balanced_accuracy: ba, sensitivity, specificity, vme: 1 - sensitivity, me: 1 - specificity, auroc: null });

type Fits = Partial<Record<Learner, { random: LearnerLite | null; lineage: LearnerLite | null }>>;

const row = (id: string, n: number, nResistant: number, best: Learner | null, fits: Fits, catalog: Lite | null = null): CatalogRow => {
  const learners = FITTED.flatMap((learner) => {
    const f = fits[learner];
    return f === undefined ? [] : [{ learner, random: f.random, lineage: f.lineage }];
  });
  const own = best === null ? undefined : fits[best];
  return {
    id,
    species: "Klebsiella pneumoniae",
    species_slug: "klebsiella_pneumoniae",
    antibiotic: id,
    drug_class: "x",
    n,
    n_resistant: nResistant,
    labels_from_mic: 0,
    features: 10,
    lineages: 12,
    best_learner: best,
    best_random: own?.random?.balanced_accuracy ?? null,
    best_lineage: own?.lineage?.balanced_accuracy ?? null,
    learners,
    catalog,
    mic_essential_agreement: null,
    scm_rules: [],
    panels: [],
    kover: null,
  };
};

describe("lowerTail", () => {
  it("sums a binomial's lower tail", () => {
    expect(lowerTail(2, 4, 0.5)).toBeCloseTo(11 / 16, 12);
    expect(lowerTail(0, 6, 0.4593)).toBeCloseTo((1 - 0.4593) ** 6, 12);
    expect(lowerTail(5, 6, 0.3)).toBeCloseTo(1 - 0.3 ** 6, 12);
  });

  it("covers the ends without a sum", () => {
    expect(lowerTail(-1, 10, 0.2)).toBe(0);
    expect(lowerTail(10, 10, 0.2)).toBe(1);
  });

  it("keeps its precision on a cohort of thousands", () => {
    // The mean is 300 events, so the tail up to 300 is near one half.
    expect(lowerTail(300, 20000, 0.015)).toBeCloseTo(0.51535, 4);
    expect(lowerTail(250, 20000, 0.015)).toBeCloseTo(0.00157, 4);
    expect(lowerTail(0, 20000, 0.5)).toBe(0);
  });
});

describe("the exact 95% interval against a bar", () => {
  // Bounds from scipy.stats.beta (Clopper-Pearson): 5 of 65 is 2.5448% to
  // 17.0456%, 0 of 6 has an upper bound of 45.9258%, 0 of 1,886 of 0.1954%.
  it("agrees with the interval's upper bound", () => {
    expect([shownAtMost(5, 65, 0.1705), shownAtMost(5, 65, 0.1704)]).toEqual([true, false]);
    expect([shownAtMost(0, 6, 0.4593), shownAtMost(0, 6, 0.4592)]).toEqual([true, false]);
    expect([shownAtMost(0, 1886, 0.001955), shownAtMost(0, 1886, 0.001953)]).toEqual([true, false]);
  });

  it("agrees with the interval's lower bound", () => {
    expect([shownAbove(5, 65, 0.0254), shownAbove(5, 65, 0.0255)]).toEqual([true, false]);
  });

  it("never shows zero events above a bar, or anything on zero genomes", () => {
    expect(shownAbove(0, 10, 0.015)).toBe(false);
    expect([shownAtMost(0, 0, 0.015), shownAbove(0, 0, 0.015)]).toEqual([false, false]);
  });

  it("needs 245 resistant genomes before a result with no miss is shown under the VME bar", () => {
    expect(ZERO_MISS_FLOOR).toBe(245);
    expect([shownAtMost(0, 245, 0.015), shownAtMost(0, 244, 0.015)]).toEqual([true, false]);
    const rows = [row("a", 600, 245, null, {}), row("b", 600, 244, null, {})];
    expect(atZeroMissFloor(rows)).toBe(1);
  });
});

describe("rateCounts", () => {
  const base = row("a", 100, 40, "gbt", {});

  it("counts a rate over the genomes the learner called", () => {
    // 8 resistant genomes have no call: 4 of 32 missed, 3 of 60 false alarms.
    expect(rateCounts(base, fit(0.9, 0.125, 0.05, { resistant: 8, susceptible: 0 }))).toEqual({ vme: { events: 4, n: 32 }, me: { events: 3, n: 60 } });
  });

  it("counts a rate over every genome when the null-fit genomes are scored", () => {
    // Revision 14: the same 8 genomes sit in null-fit folds and are scored, so 5 of 40 and 3 of 60.
    const result = scored(fit(0.9, 0.125, 0.05, { resistant: 8, susceptible: 0 }));
    expect(leftOut(result)).toEqual({ resistant: 0, susceptible: 0 });
    expect(rateCounts(base, result)).toEqual({ vme: { events: 5, n: 40 }, me: { events: 3, n: 60 } });
    // A whole class in null-fit folds still has counts: every genome has a score.
    expect(rateCounts(base, scored(fit(0.5, 0, 0.05, { resistant: 40, susceptible: 0 })))).toEqual({ vme: { events: 0, n: 40 }, me: { events: 3, n: 60 } });
    expect(leftOut(fit(0.9, 0.125, 0.05, { resistant: 8, susceptible: 0 }))).toEqual({ resistant: 8, susceptible: 0 });
    expect(leftOut(null)).toEqual({ resistant: 0, susceptible: 0 });
  });

  it("counts a known-gene rule over every genome", () => {
    expect(rateCounts(base, rule(0.9, 0.9, 0.95))).toEqual({ vme: { events: 4, n: 40 }, me: { events: 3, n: 60 } });
  });

  it("has no counts without both rates, or when a class has no called genome", () => {
    expect(rateCounts(base, null)).toBeNull();
    expect(rateCounts(base, fit(null, null, 0.05))).toBeNull();
    expect(rateCounts(base, fit(0.5, 0, 0.05, { resistant: 40, susceptible: 0 }))).toBeNull();
    expect(rateCounts(base, fit(0.5, 0.1, 0, { resistant: 0, susceptible: 60 }))).toBeNull();
  });
});

describe("barOutcome", () => {
  it("separates shown to meet, shown to miss and too few to tell", () => {
    expect(barOutcome({ vme: { events: 0, n: 300 }, me: { events: 0, n: 300 } })).toBe("meets");
    expect(barOutcome({ vme: { events: 10, n: 100 }, me: { events: 0, n: 300 } })).toBe("misses");
    expect(barOutcome({ vme: { events: 0, n: 300 }, me: { events: 3, n: 60 } })).toBe("open");
    expect(barOutcome({ vme: { events: 0, n: 50 }, me: { events: 0, n: 50 } })).toBe("open");
  });

  it("needs both rates under their bars to meet, and one over to miss", () => {
    // 0 of 300 is shown under the 1.5% bar; 0 of 50 is not shown under either bar.
    expect(barOutcome({ vme: { events: 0, n: 300 }, me: { events: 0, n: 50 } })).toBe("open");
    expect(barOutcome({ vme: { events: 0, n: 50 }, me: { events: 20, n: 100 } })).toBe("misses");
  });
});

describe("barSummary", () => {
  const both = (lineage: LearnerLite | null): Fits => ({ gbt: { random: fit(0.9), lineage } });
  const rows = [
    row("meets", 600, 300, "gbt", both(fit(1, 0, 0))),
    row("misses", 200, 100, "gbt", both(fit(0.95, 0.1, 0))),
    row("open-at-bars", 100, 50, "gbt", both(fit(1, 0, 0))),
    row("open-over-a-bar", 200, 100, "gbt", both(fit(0.99, 0.02, 0))),
    row("no-lineage", 200, 100, "gbt", both(null)),
    row("no-best", 200, 100, null, both(fit(1, 0, 0))),
  ];

  it("sorts the best learner's datasets and keeps its basis", () => {
    expect(barSummary(rows, "best", "lineage")).toEqual({ n: 4, meets: 1, misses: 1, open: 2, openAtBars: 1 });
  });

  it("reads one learner for every dataset", () => {
    expect(barSummary(rows, "gbt", "lineage").n).toBe(5);
    expect(barSummary(rows, "scm", "lineage")).toEqual({ n: 0, meets: 0, misses: 0, open: 0, openAtBars: 0 });
  });

  it("gives the known-gene rule one result for either split", () => {
    const withRule = [row("a", 200, 100, "gbt", both(null), rule(0.9, 0.9, 0.9)), row("b", 200, 100, null, {})];
    expect(resultOf(withRule[0]!, "catalog", "random")).toBe(resultOf(withRule[0]!, "catalog", "lineage"));
    expect(barSummary(withRule, "catalog", "lineage")).toEqual({ n: 1, meets: 0, misses: 1, open: 0, openAtBars: 0 });
  });
});

describe("learnerSummary", () => {
  const rows = [
    row("a", 100, 50, "gbt", { scm: { random: fit(0.9), lineage: fit(0.8) }, gbt: { random: fit(0.95), lineage: fit(0.85) } }),
    // Gradient boosting outscores the best learner once lineages are held out.
    row("b", 100, 50, "scm", { scm: { random: fit(0.92), lineage: fit(0.9) }, gbt: { random: fit(0.9), lineage: fit(0.95) } }),
    row("c", 100, 50, "gbt", { scm: { random: fit(0.7), lineage: null }, gbt: { random: fit(0.8), lineage: null } }),
    row("d", 100, 50, "scm", {
      scm: { random: fit(0.86), lineage: fit(0.75) },
      gbt: { random: fit(0.85), lineage: fit(0.7, 0, 0, { resistant: 1, susceptible: 0 }) },
    }),
  ];
  const s = learnerSummary(rows);
  const of = (learner: Learner) => s.rows.find((r) => r.learner === learner)!;

  it("has a row for each of the five learners, in order", () => {
    expect(s.rows.map((r) => r.learner)).toEqual(["scm", "cart", "rf", "gbt", "lasso"]);
  });

  it("takes each learner's medians over the datasets where it has both scores", () => {
    const scm = of("scm");
    expect(scm.n).toBe(3);
    expect(scm.medianRandom).toBeCloseTo(0.9);
    expect(scm.medianLineage).toBeCloseTo(0.8);
    expect(scm.medianGap).toBeCloseTo(0.1);
    const gbt = of("gbt");
    expect(gbt.n).toBe(3);
    expect(gbt.medianRandom).toBeCloseTo(0.9);
    expect(gbt.medianLineage).toBeCloseTo(0.85);
    expect(gbt.medianGap).toBeCloseTo(0.1);
  });

  it("sets the best learner beside each learner used alone", () => {
    const scm = of("scm");
    expect([scm.nBest, scm.bestHigher, scm.bestLower]).toEqual([3, 1, 0]);
    expect(scm.medianBestLead).toBeCloseTo(0, 6);
    const gbt = of("gbt");
    expect([gbt.nBest, gbt.bestHigher, gbt.bestLower]).toEqual([3, 1, 1]);
    expect(gbt.medianBestLead).toBeCloseTo(0, 6);
  });

  it("counts the datasets where a learner left a genome without a call", () => {
    expect([of("scm").nNoCall, of("gbt").nNoCall]).toEqual([0, 1]);
  });

  it("counts how often no learner that calls every genome outscores the best learner with lineages held out", () => {
    expect([s.nBest, s.bestAlsoTops]).toEqual([3, 2]);
  });

  it("counts the datasets where such a learner is ahead by more than the cutoff", () => {
    // Dataset b: 95% against the best learner's 90%.
    expect(s.otherAheadByCutoff).toBe(1);
    // One point ahead outscores the best learner but stays inside the 2-point cutoff.
    const near = learnerSummary([row("g", 100, 50, "scm", { scm: { random: fit(0.9), lineage: fit(0.8) }, gbt: { random: fit(0.85), lineage: fit(0.81) } })]);
    expect([near.nBest, near.bestAlsoTops, near.otherAheadByCutoff]).toEqual([1, 0, 0]);
  });

  it("does not set a learner that leaves a genome without a call against the best learner", () => {
    const partial = fit(0.99, 0, 0, { resistant: 1, susceptible: 0 });
    const one = learnerSummary([row("e", 100, 50, "scm", { scm: { random: fit(0.9), lineage: fit(0.8) }, gbt: { random: fit(0.85), lineage: partial } })]);
    expect([one.nBest, one.bestAlsoTops, one.otherAheadByCutoff]).toEqual([1, 1, 0]);
    // A revision-14 row scores its null-fit genomes, so the same row counts, and it is no longer a no-call dataset.
    const rescored = learnerSummary([row("h", 100, 50, "scm", { scm: { random: fit(0.9), lineage: fit(0.8) }, gbt: { random: fit(0.85), lineage: scored(partial) } })]);
    expect([rescored.nBest, rescored.bestAlsoTops, rescored.otherAheadByCutoff]).toEqual([1, 0, 1]);
    expect(rescored.rows.find((r) => r.learner === "gbt")!.nNoCall).toBe(0);
    expect(one.rows.find((r) => r.learner === "gbt")!.nNoCall).toBe(1);
    // The same score from a learner that calls every genome does count.
    const whole = learnerSummary([row("f", 100, 50, "scm", { scm: { random: fit(0.9), lineage: fit(0.8) }, gbt: { random: fit(0.85), lineage: fit(0.99) } })]);
    expect([whole.nBest, whole.bestAlsoTops, whole.otherAheadByCutoff]).toEqual([1, 0, 1]);
  });

  it("leaves a learner with no result empty", () => {
    const cart = of("cart");
    expect([cart.n, cart.nBest, cart.medianRandom, cart.medianLineage, cart.medianGap, cart.medianBestLead]).toEqual([0, 0, null, null, null, null]);
  });
});

describe("ruleBands", () => {
  const scored = (id: string, lineage: number | null, catalog: Lite | null): CatalogRow =>
    row(id, 100, 50, "gbt", { gbt: { random: fit(0.99), lineage: lineage === null ? null : fit(lineage) } }, catalog);
  const rows = [
    scored("all-resistant", 0.9, rule(0.5, 1, 0)),
    scored("all-susceptible", null, rule(0.5, 0.05, 0.95)),
    scored("on-the-first-cut", 0.56, rule(0.55, 0.55, 0.55)),
    scored("middle", 0.9, rule(0.8, 0.8, 0.8)),
    scored("on-the-last-cut", 0.8, rule(0.85, 0.85, 0.85)),
    scored("top", 0.99, rule(0.95, 0.95, 0.95)),
    scored("no-rule", 0.9, null),
  ];
  const bands = ruleBands(rows);

  it("puts each rule in one of four groups, a cut belonging to the group above it", () => {
    expect(bands.map((b) => b.datasets)).toEqual([2, 1, 1, 2]);
    expect(bands.map((b) => b.compared)).toEqual([1, 1, 1, 2]);
    expect(bands.map((b) => [b.from, b.to])).toEqual([
      [-Infinity, 0.55],
      [0.55, 0.7],
      [0.7, 0.85],
      [0.85, Infinity],
    ]);
  });

  it("counts the rule as holding when it is ahead or within 2 points", () => {
    expect(bands.map((b) => [b.ruleHolds, b.ruleBehind])).toEqual([
      [0, 1],
      [1, 0],
      [0, 1],
      [1, 1],
    ]);
    expect(bands[0]!.medianLearnerLead).toBeCloseTo(0.4, 6);
    expect(bands[3]!.medianLearnerLead).toBeCloseTo(-0.005, 6);
  });

  it("counts the rules that call nearly every genome one way", () => {
    expect(bands.map((b) => b.uniform)).toEqual([
      { resistant: 1, susceptible: 1 },
      { resistant: 0, susceptible: 0 },
      { resistant: 0, susceptible: 0 },
      { resistant: 0, susceptible: 0 },
    ]);
  });

  it("has four empty groups without a rule", () => {
    expect(ruleBands([scored("no-rule", 0.9, null)]).map((b) => [b.datasets, b.compared, b.medianLearnerLead])).toEqual([
      [0, 0, null],
      [0, 0, null],
      [0, 0, null],
      [0, 0, null],
    ]);
  });
});
