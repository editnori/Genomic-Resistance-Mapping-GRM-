import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../validate";
import type { Manifest } from "../schema";
import { gapWords } from "../format";
import { auditOf, LEGACY_AUDIT, manifestAudit, NO_CALL_REVISION, nullFitWords } from "./evaluation";
import { classPools, distinctResistant, joinPanels, notScored, pool, share } from "./tests";
import { withRealProbeRows } from "../../dev/fixtures/probeRows";
import { splitAbsence, splitOf } from "./dataset";

const manifest = parseManifest(fixtureFile("/data/index.json"));

describe("evaluation audit", () => {
  it("words null fits as no-calls only from the revision that ships them", () => {
    const before = nullFitWords(NO_CALL_REVISION - 1);
    const from = nullFitWords(NO_CALL_REVISION);
    expect(before).toContain("scores every isolate of that fold 0.5, which the 0.5 cutoff calls resistant");
    expect(before).not.toContain("no call");
    expect(nullFitWords(undefined)).toBe(before);
    expect(from).toContain(`From evaluation revision ${NO_CALL_REVISION} it makes no call`);
    expect(from).toContain("earlier revisions scored those isolates 0.5");
  });

  it("reads a report without an audit as the legacy, provisional audit", () => {
    expect(auditOf({})).toBe(LEGACY_AUDIT);
    expect(LEGACY_AUDIT.valid).toBe(false);
  });

  it("counts datasets that fail the audit and keeps the scope of those that pass", () => {
    const scope = "Internal cross-validation only.";
    const m: Manifest = {
      ...manifest,
      datasets: manifest.datasets.map((d, i) =>
        i === 0 ? { ...d, evaluation: { revision: 3, valid: true, reasons: [], scope } } : { ...d, evaluation: { revision: 2, valid: false, reasons: ["full cohort"] } },
      ),
    };
    const a = manifestAudit(m);
    expect([a.total, a.invalid]).toEqual([m.datasets.length, m.datasets.length - 1]);
    expect(a.reasons).toEqual(["full cohort"]);
    expect(a.scopes).toEqual([scope]);
  });

  it("marks an audit from an earlier revision provisional against the manifest's current revision", () => {
    const scope = "Internal cross-validation only.";
    const current = { revision: 6, valid: true, reasons: [], scope };
    const earlier = { revision: 5, valid: true, reasons: [], scope };
    expect(auditOf({ evaluation: current }, 6)).toEqual(current);
    expect(auditOf({ evaluation: earlier })).toEqual(earlier);
    expect(auditOf({}, 6)).toBe(LEGACY_AUDIT);
    const stale = auditOf({ evaluation: earlier }, 6);
    expect(stale.valid).toBe(false);
    expect(stale.revision).toBe(5);
    expect(stale.reasons).toHaveLength(1);
    expect(stale.reasons[0]).toContain("revision 5");
    expect(stale.reasons[0]).toContain("current revision is 6");
  });

  it("counts earlier-revision datasets as provisional in a mixed rebuild and keeps the scope of current ones", () => {
    const scope = "Internal cross-validation only.";
    const mixed: Manifest = {
      ...manifest,
      evaluation_revision: 6,
      datasets: manifest.datasets.map((d, i) => ({ ...d, evaluation: { revision: i === 0 ? 6 : 5, valid: true, reasons: [], scope } })),
    };
    const a = manifestAudit(mixed);
    expect([a.total, a.invalid]).toEqual([mixed.datasets.length, mixed.datasets.length - 1]);
    expect(a.reasons).toHaveLength(1);
    expect(a.reasons[0]).toContain("revision 5");
    expect(a.scopes).toEqual([scope]);
    expect(manifestAudit({ ...mixed, evaluation_revision: undefined }).invalid).toBe(0);
  });

  it("rejects an evaluation block of the wrong shape at the edge", () => {
    const raw = fixtureFile(`/data/datasets/${manifest.datasets[0]!.id}.json`) as Record<string, unknown>;
    expect(() => parseDataset({ ...raw, evaluation: { revision: 3, valid: "yes", reasons: [] } })).toThrow(/evaluation/);
    expect(() => parseDataset({ ...raw, evaluation: { revision: 3, valid: true, reasons: [], scope: "internal" } })).not.toThrow();
  });
});

describe("panel pools", () => {
  const rows = joinPanels(manifest);

  it("counts each dataset's resistant isolates once however many panels score them", () => {
    const d = distinctResistant(rows);
    const byDataset = new Map(rows.map((r) => [r.dataset, r.resistant]));
    expect(d.datasets).toBe(byDataset.size);
    expect(d.resistant).toBe([...byDataset.values()].reduce((a, b) => a + b, 0));
    expect(rows.reduce((a, r) => a + r.resistant, 0)).toBeGreaterThan(d.resistant);
  });

  it("picks each dataset's panel by balanced accuracy, not by resistant isolates caught", () => {
    const base = rows.find((r) => rows.filter((x) => x.dataset === r.dataset).length > 1)!;
    const broad = { ...base, panel: "broad", metrics: { ...base.metrics, sensitivity: 1, specificity: 0.4 }, caught: base.resistant, falseCalls: 999 };
    const narrow = { ...base, panel: "narrow", metrics: { ...base.metrics, sensitivity: 0.95, specificity: 0.98 }, caught: Math.round(base.resistant * 0.95), falseCalls: 1 };
    const [pool] = classPools([broad, narrow]);
    expect(pool!.falseCalls).toBe(1);
  });

  it("takes exact counts from a panel row and falls back to rates times n without them", () => {
    const base = manifest.panels[0]!;
    const row = manifest.datasets.find((d) => d.id === base.dataset)!;
    const rates = { ...base, metrics: { ...base.metrics, sensitivity: 0.5, specificity: 0.5 }, best_model: { ...base.best_model, sensitivity: 0.9 } };
    const exact = { ...rates, counts: { n: 20, tp: 7, fp: 3, tn: 9, fn: 1 }, best_model_counts: { n: 20, tp: 8, fp: 2, tn: 10, fn: 0 } };
    const partial = { ...exact, best_model_counts: undefined };
    const [fromRates, fromCounts, fromBoth] = joinPanels({ ...manifest, panels: [rates, exact, partial] });
    expect(fromCounts).toMatchObject({ resistant: 8, caught: 7, modelCaught: 8, susceptible: 12, falseCalls: 3 });
    expect(fromRates).toMatchObject({
      resistant: row.n_resistant,
      caught: Math.round(0.5 * row.n_resistant),
      modelCaught: Math.round(0.9 * row.n_resistant),
      susceptible: row.n - row.n_resistant,
      falseCalls: Math.round(0.5 * (row.n - row.n_resistant)),
    });
    expect(fromBoth).toMatchObject({ resistant: 8, caught: 7, modelCaught: Math.round(0.9 * 8), falseCalls: 3 });
    expect(fromCounts).toMatchObject({ noCallResistant: 0, noCallSusceptible: 0 });
  });

  it("scores real probe rows on the isolates each method calls and keeps no-calls apart", () => {
    const rows = joinPanels(withRealProbeRows(manifest));
    rows.forEach((r) => {
      const d = r.row!;
      expect(r.resistant).toBe(d.n_resistant);
      expect(r.susceptible).toBe(d.n - d.n_resistant);
      expect(r.calledResistant + r.noCallResistant).toBe(r.resistant);
      expect(r.calledSusceptible + r.noCallSusceptible).toBe(r.susceptible);
      expect(r.caught / r.calledResistant).toBeCloseTo(r.metrics.sensitivity, 12);
      expect(1 - r.falseCalls / r.calledSusceptible).toBeCloseTo(r.metrics.specificity, 12);
      expect(r.modelResistant).toBe(r.resistant);
      expect(r.modelCaught / r.modelResistant).toBeCloseTo(r.best_model.sensitivity, 12);
    });
    expect(rows.map((r) => [r.noCallResistant, r.noCallSusceptible])).toEqual([
      [0, 0],
      [4, 2],
      [4, 2],
    ]);
    expect(pool(rows)).toEqual({
      datasets: 3,
      resistant: 1691,
      caught: 1580,
      calledResistant: 1683,
      modelResistant: 1691,
      modelCaught: 1571,
      susceptible: 1659,
      calledSusceptible: 1655,
      falseCalls: 177,
      noCallResistant: 8,
      noCallSusceptible: 4,
    });
  });

  it("lets no abstention move the panel's or the learner's share", () => {
    const m = withRealProbeRows(manifest);
    const [, withNoCalls] = joinPanels(m);
    const [without] = joinPanels({ ...m, panels: [{ ...m.panels[1]!, no_call: undefined }] });
    expect(share(withNoCalls!.caught, withNoCalls!.calledResistant)).toBe(share(without!.caught, without!.calledResistant));
    expect(share(withNoCalls!.falseCalls, withNoCalls!.calledSusceptible)).toBe(share(without!.falseCalls, without!.calledSusceptible));
    expect(share(withNoCalls!.modelCaught, withNoCalls!.modelResistant)).toBe(share(without!.modelCaught, without!.modelResistant));
    expect([withNoCalls!.resistant - without!.resistant, withNoCalls!.susceptible - without!.susceptible]).toEqual([4, 2]);
  });

  it("words an unsupported reason as one sentence", () => {
    expect(notScored("gene rpoB has no verified wild-type sequence")).toBe("Not scored. Gene rpoB has no verified wild-type sequence.");
    expect(notScored("In E. coli, synonymous changes read as mutations.")).toBe("Not scored. In E. coli, synonymous changes read as mutations.");
  });
});

describe("splits", () => {
  const report = parseDataset(fixtureFile(`/data/datasets/${manifest.datasets[0]!.id}.json`));

  it("treats a planned split with no learners as absent and gives its reason", () => {
    const reason = "Outer folds cannot hold both classes.";
    const r = { ...report, splits: [...report.splits.filter((s) => s.scheme === "random"), { scheme: "lineage" as const, folds: 5, n_groups: 7, learners: [], unavailable_reason: reason }] };
    expect(splitOf(r, "lineage")).toBeUndefined();
    expect(splitAbsence(r, "lineage")).toBe(reason);
    expect(splitAbsence(r, "random")).toBeNull();
  });

  it("words a per-dataset difference by its sign", () => {
    expect(gapWords(0.031)).toBe("3.1 points lower with lineages held out");
    expect(gapWords(-0.02)).toBe("2.0 points higher with lineages held out");
  });
});
