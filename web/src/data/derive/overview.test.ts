import { describe, expect, it } from "vitest";
import type { CatalogRow, Lite } from "../schema";
import { catalogSummary, classMinimum, classRows, datasetRule, noBestWords, splitSummary, uniformCall } from "./overview";

const lite = (ba: number): Lite => ({ balanced_accuracy: ba, sensitivity: ba, specificity: ba, vme: 1 - ba, me: 1 - ba, auroc: null });

const row = (id: string, drugClass: string, random: number, lineage: number | null, catalog: number | null): CatalogRow => ({
  id,
  species: "Klebsiella pneumoniae",
  species_slug: "klebsiella_pneumoniae",
  antibiotic: id,
  drug_class: drugClass,
  n: 100,
  n_resistant: 50,
  labels_from_mic: 0,
  features: 10,
  lineages: 12,
  best_learner: "gbt",
  best_random: random,
  best_lineage: lineage,
  learners: [],
  catalog: catalog === null ? null : lite(catalog),
  mic_essential_agreement: null,
  scm_rules: [],
  panels: [],
  kover: null,
});

describe("splitSummary", () => {
  it("takes the median of per-dataset gaps, not the gap between medians", () => {
    // Gaps 0.01, 0.10, 0.02 -> median 0.02. Medians 0.9 and 0.8 differ by 0.1.
    const rows = [row("a", "x", 0.95, 0.94, null), row("b", "x", 0.9, 0.8, null), row("c", "x", 0.8, 0.78, null)];
    const s = splitSummary(rows);
    expect(s.medianGap).toBeCloseTo(0.02);
    expect(s.medianRandom).toBeCloseTo(0.9);
    expect(s.medianLineage).toBeCloseTo(0.8);
    expect(s.nLargeGap).toBe(1);
  });

  it("uses the same paired datasets for all three medians", () => {
    const s = splitSummary([row("a", "x", 0.1, null, null), row("b", "x", 0.8, 0.7, null)]);
    expect([s.n, s.nLineage]).toEqual([2, 1]);
    expect(s.medianRandom).toBeCloseTo(0.8);
    expect(s.medianLineage).toBeCloseTo(0.7);
    expect(s.medianGap).toBeCloseTo(0.1);
  });

  it("reports no medians when no dataset has both splits", () => {
    const s = splitSummary([row("a", "x", 0.9, null, null)]);
    expect([s.medianRandom, s.medianLineage, s.medianGap]).toEqual([null, null, null]);
  });

  it("leaves out and counts a dataset where no learner calls every isolate", () => {
    const none = { ...row("a", "x", 0.9, 0.8, null), best_learner: null, best_random: null, best_lineage: null };
    const s = splitSummary([none, row("b", "x", 0.8, 0.7, null)]);
    expect([s.n, s.nLineage, s.nNoBest]).toEqual([1, 1, 1]);
    expect(s.medianRandom).toBeCloseTo(0.8);
    expect(s.medianGap).toBeCloseTo(0.1);
    expect(noBestWords(s.nNoBest)).toBe(" 1 dataset has no learner that calls every isolate and is left out.");
    expect(noBestWords(2, " with a simulated panel")).toBe(" 2 datasets with a simulated panel have no learner that calls every isolate and are left out.");
    expect(noBestWords(0)).toBe("");
  });
});

describe("catalogSummary", () => {
  it("counts datasets where the catalog is within each cutoff of the lineage-held-out learner", () => {
    const rows = [
      row("a", "x", 0.99, 0.95, 0.945), // 0.5 points behind
      row("b", "x", 0.99, 0.95, 0.92), // 3 points behind
      row("c", "x", 0.99, 0.95, 0.97), // ahead
      row("d", "x", 0.99, null, 0.9), // no lineage split: not compared
      row("e", "x", 0.99, 0.95, null), // no rule
    ];
    const s = catalogSummary(rows);
    expect([s.nWithRule, s.nCompared]).toEqual([4, 3]);
    expect([s.within.get(0.01), s.within.get(0.02), s.within.get(0.05)]).toEqual([2, 2, 3]);
  });
});

describe("classRows", () => {
  it("marks a class the catalog explains, one it does not, and one with no rule", () => {
    const rows = [
      row("meropenem", "carbapenem", 0.98, 0.97, 0.965),
      row("ciprofloxacin", "quinolone", 0.95, 0.93, 0.7),
      row("rifampicin", "rifamycin", 0.95, 0.93, null),
    ];
    const verdicts = Object.fromEntries(classRows(rows).map((c) => [c.drugClass, c.verdict]));
    expect(verdicts).toEqual({ carbapenem: "catalog", quinolone: "beyond", rifamycin: "no-rule" });
  });
});

describe("classMinimum", () => {
  const at = (id: string, min_class?: number): CatalogRow => ({ ...row(id, "x", 0.9, 0.8, null), min_class });

  it("states the one minimum every report records", () => {
    const c = classMinimum([at("a", 10), at("b", 10)]);
    expect(c).toEqual({ kind: "value", value: 10 });
    expect(datasetRule(c)).toBe("at least 10 resistant and 10 susceptible genomes that passed QC and have a reported lab call");
  });

  it("lists each minimum with its reports", () => {
    const mixed = classMinimum([at("a", 25), at("b", 10), at("c", 25)]);
    expect(mixed).toEqual({ kind: "mixed", reports: [[10, 1], [25, 2]] });
    expect(datasetRule(mixed)).toBe(
      "enough resistant and susceptible genomes that passed QC and have a reported lab call: at least 10 of each in 1 report and at least 25 of each in 2 reports",
    );
  });

  it("reads a report that records no minimum as the old learn gate of 25", () => {
    expect(classMinimum([at("a"), at("b")])).toEqual({ kind: "value", value: 25 });
    expect(classMinimum([at("a", 25), at("b")])).toEqual({ kind: "value", value: 25 });
    const partial = classMinimum([at("a", 10), at("b")]);
    expect(partial).toEqual({ kind: "mixed", reports: [[10, 1], [25, 1]] });
    expect(datasetRule(partial)).toContain("at least 10 of each in 1 report and at least 25 of each in 1 report");
  });

  it("is unavailable before any dataset has a report", () => {
    const c = classMinimum([]);
    expect(c).toEqual({ kind: "unavailable" });
    expect(datasetRule(c)).toBe("enough resistant and susceptible genomes that passed QC and have a reported lab call; no dataset report records the minimum yet");
  });
});


describe("catalog display regions", () => {
  it("partitions only paired catalog comparisons into ahead, within, and behind", () => {
    const s = catalogSummary([row("ahead", "x", 0.9, 0.6, 0.7), row("within", "x", 0.9, 0.7, 0.71), row("behind", "x", 0.9, 0.8, 0.7), row("unpaired", "x", 0.9, null, 0.99), row("no rule", "x", 0.9, 0.8, null)]);
    expect(s.regions).toEqual({ ahead: 1, within: 1, behind: 1 });
    expect(Object.values(s.regions).reduce((a, b) => a + b, 0)).toBe(s.nCompared);
    expect(s.within.get(0.02)).toBe(s.regions.ahead + s.regions.within);
  });
  it("marks the stated ≥90% calls from sensitivity and specificity, including the boundary", () => {
    expect(uniformCall({ ...lite(0.5), specificity: 0.1 })).toBe("resistant");
    expect(uniformCall({ ...lite(0.5), sensitivity: 0.1 })).toBe("susceptible");
    expect(uniformCall({ ...lite(0.5), specificity: 0.100001, sensitivity: 0.100001 })).toBeNull();
  });
});
