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
  parseBaselines,
  signedPoints,
  type Baselines,
} from "./baselines";
import { ShapeError } from "./validate";
import type { CatalogRow, DatasetReport } from "./schema";
import { SHA, agreement, file, lineageEntry, micPart } from "../dev/fixtures/baselines";

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
    // Confirmed with no group, or a group of one: refused.
    expect(refused(file({ mic: { ...m, prediction_check: { confirmed: true, groups: 0, isolates: 0, constant_folds: [] } } }))).toContain("prediction_check");
    expect(refused(file({ mic: { ...m, prediction_check: { confirmed: true, groups: 2, isolates: 3, constant_folds: [] } } }))).toContain("prediction_check");
    expect(refused(file({ mic: { ...m, prediction_check: { ...m.prediction_check, constant_folds: [3, 1] } } }))).toContain("ascending");
    expect(refused(file({ mic: { ...m, threshold_check: { confirmed: true, checked: 0 } } }))).toContain("threshold_check");
    expect(refused(file({ mic: { ...m, exact_with_training_copy: 97 } }))).toContain("exact_with_training_copy");
    // A file whose thresholds confirm nothing still reads as confirmed by
    // the predictions, and the reverse reads as unconfirmed.
    const thresholdsOnly = parseBaselines(file({ mic: { ...micPart(96, 58, 81, 60, false), threshold_check: { confirmed: true, checked: 5 } } }));
    const row = { id: "synthetic__fixturemycin", mic_essential_agreement: 58 / 96 } as CatalogRow;
    expect(micConstantFor(row, thresholdsOnly, 14)?.confirmed).toBe(false);
    const predictionsOnly = parseBaselines(file({ mic: { ...m, threshold_check: { confirmed: false, checked: 0 } } }));
    expect(micConstantFor(row, predictionsOnly, 14)?.confirmed).toBe(true);
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
  });
});

describe("the constants beside the index (synthetic)", () => {
  const row = (id: string, ea: number | null) => ({ id, mic_essential_agreement: ea }) as CatalogRow;
  const named = (id: string, model: number, median: number, confirmed = true): Baselines =>
    parseBaselines(file({ dataset: id, mic: micPart(100, model, median, median, confirmed) }));

  it("pins on dataset, revision and the model's own figure", () => {
    const f = named("a__drug", 90, 80);
    expect(micConstantFor(row("a__drug", 0.9), f, 14)).toEqual({ model: 0.9, constant: 0.8, mode: 0.8, confirmed: true });
    expect(micConstantFor(row("b__drug", 0.9), f, 14)).toBeNull();
    expect(micConstantFor(row("a__drug", 0.9), f, 13)).toBeNull();
    expect(micConstantFor(row("a__drug", 0.91), f, 14)).toBeNull();
    expect(micConstantFor(row("a__drug", 0.9), undefined, 14)).toBeNull();
  });

  it("counts a model level with its constant as not above it", () => {
    const rows = [row("a__drug", 0.9), row("b__drug", 0.6), row("c__drug", 0.8), row("d__drug", 0.7), row("e__drug", null)];
    const files = new Map([
      ["a__drug", named("a__drug", 90, 92)],
      ["b__drug", named("b__drug", 60, 50, false)],
      ["c__drug", named("c__drug", 80, 80)],
    ]);
    const s = micBaselineSummary(rows, micConstants(rows, files, 14));
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
