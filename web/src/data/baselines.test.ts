// The baselines loader and its pin. Every file and report here is
// synthetic: counts written by the test, not BV-BRC data.
import { describe, expect, it } from "vitest";
import {
  baselinesFor,
  leadOver,
  leadWords,
  micBaselineSummary,
  micConstantFor,
  micConstants,
  micSummaryFor,
  parseBaselines,
  parseMicSummary,
  signedPoints,
} from "./baselines";
import { ShapeError } from "./validate";
import type { CatalogRow, DatasetReport, Manifest } from "./schema";
import { SHA, agreement, file, lineageEntry, micPart, summary, summaryEntry } from "../dev/fixtures/baselines";

const refused = (raw: unknown): string => {
  try {
    parseBaselines(raw);
  } catch (e) {
    expect(e).toBeInstanceOf(ShapeError);
    return (e as Error).message;
  }
  throw new Error("accepted");
};

describe("parseBaselines (synthetic files)", () => {
  it("accepts a file whose figures are its own counts'", () => {
    const b = parseBaselines(file());
    expect(b.mic?.model.essential_agreement).toBe(58 / 96);
    expect(b.lineage.map((s) => s.scheme)).toEqual(["random", "lineage"]);
  });

  it("refuses another schema as a version mismatch", () => {
    try {
      parseBaselines(file({ schema: 2 }));
      throw new Error("accepted");
    } catch (e) {
      expect(e).toBeInstanceOf(ShapeError);
      expect((e as ShapeError).foundSchema).toBe(2);
    }
  });

  it("refuses figures that are not their counts'", () => {
    expect(refused(file({ report_sha256: "abc" }))).toContain("SHA-256");
    const m = micPart(96, 58, 81, 60);
    expect(refused(file({ mic: { ...m, median: { ...m.median, essential_agreement: 0.85 } } }))).toContain("agree / n");
    expect(refused(file({ mic: { ...m, model: agreement(58, 95) } }))).toContain("96 exact MICs");
    expect(refused(file({ mic: { ...m, fold_constants: m.fold_constants.slice(1) } }))).toContain("one fold_constants entry per fold");
    const moved = m.fold_constants.map((f, i) => (i === 0 ? { ...f, test_exact: f.test_exact + 1 } : f));
    expect(refused(file({ mic: { ...m, fold_constants: moved } }))).toContain("add up to n_exact");
    const noTrain = m.fold_constants.map((f, i) => (i === 0 ? { ...f, train_exact: 0 } : f));
    expect(refused(file({ mic: { ...m, fold_constants: noTrain } }))).toContain("null exactly when");
    const random = lineageEntry("random", 120, 50, 10, 45, 8);
    expect(refused(file({ lineage: [{ ...random, balanced_accuracy: random.balanced_accuracy + 0.01 }] }))).toContain("counts' own");
    expect(refused(file({ lineage: [{ ...random, unseen_share: 0.5 }] }))).toContain("unseen / n");
    expect(refused(file({ lineage: [random, random] }))).toContain("one lineage entry per split");
  });

  it("reads the prediction check as the fold check and keeps the threshold check apart", () => {
    const m = micPart(96, 58, 81, 60);
    const check = (groups: number, isolates: number, spanning: number, confirmed: boolean) => file({ mic: { ...m, prediction_check: { confirmed, groups, isolates, spanning, constant_folds: [] } } });
    // Confirmed with no group, or a group of one: refused.
    expect(refused(check(0, 0, 0, true))).toContain("prediction_check");
    expect(refused(check(2, 3, 0, true))).toContain("prediction_check");
    expect(refused(file({ mic: { ...m, prediction_check: { ...m.prediction_check, constant_folds: [3, 1] } } }))).toContain("ascending");
    expect(refused(file({ mic: { ...m, threshold_check: { confirmed: true, checked: 0 } } }))).toContain("threshold_check");
    expect(refused(file({ mic: { ...m, exact_with_training_copy: 97 } }))).toContain("exact_with_training_copy");
    // A file whose thresholds confirm nothing still reads as confirmed by
    // the predictions, and the reverse reads as unconfirmed.
    const thresholdsOnly = parseBaselines(file({ mic: { ...micPart(96, 58, 81, 60, false), threshold_check: { confirmed: true, checked: 5 } } }));
    expect(thresholdsOnly.mic?.prediction_check.confirmed).toBe(false);
    const predictionsOnly = parseBaselines(file({ mic: { ...m, threshold_check: { confirmed: false, checked: 0 } } }));
    expect(predictionsOnly.mic?.prediction_check.confirmed).toBe(true);
  });

  it("counts the groups spanning folds: none to be confirmed, at most the inside ones to be written", () => {
    const m = micPart(96, 58, 81, 60);
    const check = (groups: number, isolates: number, spanning: number, confirmed: boolean) => file({ mic: { ...m, prediction_check: { confirmed, groups, isolates, spanning, constant_folds: [] } } });
    // 1 spanning group beside 3 inside ones: written, unconfirmed.
    expect(parseBaselines(check(3, 6, 1, false)).mic?.prediction_check).toEqual({ confirmed: false, groups: 3, isolates: 6, spanning: 1, constant_folds: [] });
    // As many spanning as inside: still written.
    expect(parseBaselines(check(2, 4, 2, false)).mic?.prediction_check.spanning).toBe(2);
    // Confirmed with a spanning group, or unconfirmed with none: refused.
    expect(refused(check(3, 6, 1, true))).toContain("none spanning folds");
    expect(refused(check(3, 6, 0, false))).toContain("none spanning folds");
    // More spanning than inside: the command refuses the dataset, so a file
    // that says so is refused.
    expect(refused(check(2, 4, 3, false))).toContain("outnumber");
    expect(refused(check(0, 0, 1, false))).toContain("outnumber");
    // A file with no count of the spanning groups: refused.
    const { spanning: _, ...old } = m.prediction_check;
    expect(refused(file({ mic: { ...m, prediction_check: old } }))).toContain("spanning");
  });
});

describe("baselinesFor: the pin (synthetic)", () => {
  const report = (over: Partial<DatasetReport> = {}): DatasetReport =>
    ({
      id: "synthetic__fixturemycin",
      evaluation: { revision: 14, valid: true, reasons: [], scope: "" },
      samples: new Array(120).fill({}),
      splits: [
        { scheme: "random", folds: 5 },
        { scheme: "lineage", folds: 5 },
      ],
      mic: { n_exact: 96, essential_agreement: 58 / 96 },
      ...over,
    }) as unknown as DatasetReport;
  const b = parseBaselines(file());

  it("gives the file when dataset, revision and report hash all match", () => {
    const got = baselinesFor(b, report(), SHA);
    expect(got?.mic?.median?.agree).toBe(81);
    expect(got?.lineage).toHaveLength(2);
  });

  it("gives nothing for an absent file, an unknown digest, or another report", () => {
    expect(baselinesFor(null, report(), SHA)).toBeNull();
    expect(baselinesFor(b, report(), undefined)).toBeNull();
    expect(baselinesFor(b, report(), null)).toBeNull();
    expect(baselinesFor(b, report(), "b".repeat(64))).toBeNull();
    expect(baselinesFor(b, report({ id: "other__drug" }), SHA)).toBeNull();
    expect(baselinesFor(b, report({ evaluation: { revision: 13, valid: true, reasons: [], scope: "" } }), SHA)).toBeNull();
    expect(baselinesFor(b, report({ evaluation: undefined }), SHA)).toBeNull();
  });

  it("keeps the MIC part only for the report's own model, and a split only where the report has it", () => {
    expect(baselinesFor(b, report({ mic: null }), SHA)?.mic).toBeNull();
    expect(baselinesFor(b, report({ mic: { n_exact: 96, essential_agreement: 0.6 } } as Partial<DatasetReport>), SHA)?.mic).toBeNull();
    const one = baselinesFor(b, report({ splits: [{ scheme: "random", folds: 5 }] } as Partial<DatasetReport>), SHA);
    expect(one?.lineage.map((s) => s.scheme)).toEqual(["random"]);
    expect(baselinesFor(b, report({ samples: new Array(119).fill({}) }), SHA)?.lineage).toEqual([]);
  });
});

describe("leads and their sign", () => {
  it("is positive when the score is above its baseline", () => {
    expect(leadOver(0.9453, 0.7724)).toBeCloseTo(0.1729, 10);
    expect(signedPoints(leadOver(0.9453, 0.7724))).toBe("+17.3");
    expect(signedPoints(leadOver(58 / 96, 81 / 96))).toBe("−24.0");
    expect(signedPoints(0)).toBe("0.0");
    expect(leadWords(leadOver(58 / 96, 81 / 96))).toBe("24.0 points below it");
    expect(leadWords(leadOver(0.9, 0.8))).toBe("10.0 points above it");
    expect(leadWords(0)).toBe("level with it");
    // A lead that shows as 0.0 points is level, on either side of 0.
    expect(leadWords(0.0004)).toBe("level with it");
    expect(leadWords(-0.0004)).toBe("level with it");
    expect(leadWords(0.01)).toBe("1.0 point above it");
    expect(leadWords(-0.033, "lineage identity alone")).toBe("3.3 points below lineage identity alone");
    expect(leadWords(0.0004, "lineage identity alone")).toBe("level with lineage identity alone");
  });
});

describe("the summary beside the index (synthetic)", () => {
  const row = (id: string, ea: number | null) => ({ id, mic_essential_agreement: ea }) as CatalogRow;
  const INDEX = "d".repeat(64);
  const manifestOf = (rows: CatalogRow[], revision = 14) => ({ evaluation_revision: revision, datasets: rows }) as unknown as Manifest;
  const rows = [row("a__drug", 0.9), row("b__drug", 0.6), row("c__drug", 0.8), row("d__drug", 0.7), row("e__drug", null)];
  const pinned = parseMicSummary(summary(INDEX, [summaryEntry("a__drug", 0.9, 0.92), summaryEntry("b__drug", 0.6, 0.5, false), summaryEntry("c__drug", 0.8, 0.8)]));

  const notParsed = (raw: unknown): string => {
    try {
      parseMicSummary(raw);
    } catch (e) {
      expect(e).toBeInstanceOf(ShapeError);
      return (e as Error).message;
    }
    throw new Error("accepted");
  };

  it("reads its own schema, a SHA-256 pin and one entry per dataset", () => {
    expect(pinned.datasets.map((d) => d.dataset)).toEqual(["a__drug", "b__drug", "c__drug"]);
    try {
      parseMicSummary({ ...summary(INDEX, []), schema: 2 });
      throw new Error("accepted");
    } catch (e) {
      expect((e as ShapeError).foundSchema).toBe(2);
    }
    expect(notParsed({ ...summary(INDEX, []), index_sha256: "abc" })).toContain("index_sha256");
    expect(notParsed(summary(INDEX, [summaryEntry("a__drug", 0.9, 1.2)]))).toContain("datasets[0]");
    expect(notParsed(summary(INDEX, [{ ...summaryEntry("a__drug", 0.9, 0.8), report_sha256: "x" }]))).toContain("datasets[0]");
    expect(notParsed(summary(INDEX, [summaryEntry("a__drug", 0.9, 0.8), summaryEntry("a__drug", 0.9, 0.8)]))).toContain("one entry per dataset");
  });

  it("is taken whole only for the index bytes it is pinned to, and refused whole otherwise", () => {
    const m = manifestOf(rows);
    expect(micSummaryFor(pinned, m, INDEX)).toBe(pinned);
    // Another index's bytes, an unknown digest, another revision: refused.
    expect(micSummaryFor(pinned, m, "e".repeat(64))).toBeNull();
    expect(micSummaryFor(pinned, m, undefined)).toBeNull();
    expect(micSummaryFor(pinned, m, null)).toBeNull();
    expect(micSummaryFor(pinned, manifestOf(rows, 13), INDEX)).toBeNull();
    expect(micSummaryFor(null, m, INDEX)).toBeNull();
    // One entry the index does not hold, or holds with another model figure
    // or no MIC model, refuses every entry.
    expect(micSummaryFor(parseMicSummary(summary(INDEX, [...pinned.datasets, summaryEntry("z__drug", 0.5, 0.5)])), m, INDEX)).toBeNull();
    expect(micSummaryFor(parseMicSummary(summary(INDEX, [summaryEntry("a__drug", 0.91, 0.92), summaryEntry("c__drug", 0.8, 0.8)])), m, INDEX)).toBeNull();
    expect(micSummaryFor(parseMicSummary(summary(INDEX, [summaryEntry("e__drug", 0.5, 0.5)])), m, INDEX)).toBeNull();
    // A summary that covers some datasets is taken as it is.
    const some = parseMicSummary(summary(INDEX, [summaryEntry("c__drug", 0.8, 0.8)]));
    expect(micSummaryFor(some, m, INDEX)).toBe(some);
  });

  it("gives each covered row its constant, on the summary's revision only", () => {
    expect(micConstantFor(row("a__drug", 0.9), pinned.datasets[0])).toEqual({ model: 0.9, constant: 0.92, confirmed: true });
    expect(micConstantFor(row("b__drug", 0.9), pinned.datasets[0])).toBeNull();
    expect(micConstantFor(row("a__drug", 0.91), pinned.datasets[0])).toBeNull();
    expect(micConstantFor(row("a__drug", 0.9), undefined)).toBeNull();
    expect([...micConstants(rows, pinned, 14).keys()]).toEqual(["a__drug", "b__drug", "c__drug"]);
    expect(micConstants(rows, pinned, 13).size).toBe(0);
    expect(micConstants(rows, null, 14).size).toBe(0);
  });

  it("counts a model level with its constant as not above it", () => {
    const s = micBaselineSummary(rows, micConstants(rows, pinned, 14));
    expect(s).toEqual({
      n: 4,
      covered: 3,
      medianModel: 0.8,
      medianConstant: 0.8,
      medianLead: 0,
      notAbove: 2,
      constantAt90: 1,
      modelAt90: 1,
      bothAt90: 1,
      unconfirmed: 1,
    });
  });
});
