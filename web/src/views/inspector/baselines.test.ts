// @vitest-environment happy-dom
// The baselines beside a dataset report's scores. The files are synthetic
// (counts written by the test); the reports are the dev fixtures'.
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { lineageEntry, micPart } from "../../dev/fixtures/baselines";
import { parseManifest } from "../../data/validate";
import { sha256Hex } from "../../data/load";
import { leadOver, micConstants, parseBaselines, signedPoints, type MicConstant } from "../../data/baselines";
import { pct } from "../../data/format";
import type { CatalogRow, DatasetReport, Manifest } from "../../data/schema";
import type { Navigator } from "../../app/router";
import { ThemeProvider } from "../../theme/ThemeProvider";
import { DatasetInspector } from "./DatasetInspector";
import { micFinding, micLead } from "../findings";
import { datasetColumns } from "../DatasetsView";
import { constantAt90Sentence } from "../MethodsView";
import { micBaselineSummary } from "../../data/baselines";
import { micSummary } from "../../data/derive/overview";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const manifest = parseManifest(fixtureFile("/data/index.json"));
const withMic = manifest.datasets.findIndex((d) => d.mic_essential_agreement !== null);
const passing = { revision: 14, valid: true, reasons: [], scope: "Internal cross-validation only." };

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

/** A fixture report under a fresh id (the caches key on the path), with an
 *  MIC model whose agreement is 58 of its exact MICs. */
const reportFor = (id: string): DatasetReport & Record<string, unknown> => {
  const base = fixtureFile(`/data/datasets/${manifest.datasets[withMic]!.id}.json`) as DatasetReport & Record<string, unknown>;
  const mic = base.mic!;
  return { ...base, id, evaluation: passing, mic: { ...mic, essential_agreement: 58 / mic.n_exact } };
};

/** A synthetic baselines file for `report`, pinned to `sha`. */
const fileFor = (report: DatasetReport, sha: string, confirmed = true) => {
  const n = report.samples.length;
  const nExact = report.mic!.n_exact;
  const tp = Math.floor(n * 0.3);
  const fp = Math.floor(n * 0.2);
  const tn = Math.floor(n * 0.3);
  return {
    schema: 1,
    dataset: report.id,
    report_sha256: sha,
    evaluation_revision: 14,
    mic: micPart(nExact, 58, Math.min(nExact, 81), Math.min(nExact, 60), confirmed),
    lineage: report.splits.map((s) => (s.scheme === "lineage" ? lineageEntry("lineage", n, tp, fp, tn, n, s.folds) : lineageEntry(s.scheme as "random", n, tp, fp, tn, 7, s.folds))),
  };
};

const shaOf = async (text: string): Promise<string> => (await sha256Hex(new TextEncoder().encode(text).buffer as ArrayBuffer))!;

/** Mounts the inspector for `report`, serving its bytes and, when given,
 *  the baselines file; any other path answers 404. */
const mount = async (report: DatasetReport, baselines: unknown | null) => {
  const text = JSON.stringify(report);
  vi.stubGlobal("fetch", async (path: string) => {
    if (path === `/data/datasets/${report.id}.json`) return new Response(text, { headers: { "Content-Type": "application/json" } });
    if (path === `/data/baselines/${report.id}.json` && baselines !== null) return new Response(JSON.stringify(baselines), { headers: { "Content-Type": "application/json" } });
    return new Response("", { status: 404 });
  });
  const row: CatalogRow = { ...manifest.datasets[withMic]!, id: report.id, evaluation: passing };
  const current: Manifest = { ...manifest, evaluation_revision: 14, datasets: [row] };
  const nav: Navigator = { route: { view: "datasets", species: null, dataset: report.id }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
  const host = document.createElement("div");
  document.body.append(host);
  const r = createRoot(host);
  root = r;
  await act(async () => {
    r.render(createElement(ThemeProvider, { children: createElement(DatasetInspector, { id: report.id, manifest: current, nav }) }));
  });
  for (let i = 0; i < 3; i++)
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });
  return text;
};

const marks = () => [...document.querySelectorAll("[data-baseline]")];

describe("the dataset report beside its baselines file", () => {
  it("gives the model, the constant and the difference in one sentence, with the mode in the detail", async () => {
    const report = reportFor("baselines-shown-0");
    const sha = await shaOf(JSON.stringify(report));
    await mount(report, fileFor(report, sha));
    const sentence = document.querySelector('[data-baseline="mic"]')!;
    const n = report.mic!.n_exact;
    expect(sentence.textContent).toContain(`agrees on ${pct(81 / n)} of the same ${n} exact MICs (81 of ${n})`);
    expect(sentence.textContent).toContain(`the model agrees on ${pct(58 / n)} (58 of ${n}), ${((23 / n) * 100).toFixed(1)} points below it.`);
    expect(sentence.textContent).not.toContain("†");
    const text = document.body.textContent ?? "";
    expect(text).toContain(`At the mode of each training fold's exact MICs instead, the constant agrees on ${pct(60 / n)} (60 of ${n})`);
    expect(text).toContain("2 exact-MIC test isolates here have a k-mer-identical copy in their training fold.");
    // The heatmap's band sentence carries the constant.
    expect(text).toContain(`of the measured MIC: ${pct(58 / n)} of ${n} exact MICs for the model, ${pct(81 / n)} for a constant MIC from the training folds`);
    // The MIC folds are their own draw.
    expect(text).toContain("predicted out of fold on its own stratified random folds, which differ from the random split's and hold no lineage out.");
    expect(text).not.toContain("on the random split. ");
  });

  it("sets each learner's random-split score beside lineage identity, with the sign of score minus baseline, and no difference on the lineage split", async () => {
    const report = reportFor("baselines-shown-1");
    const sha = await shaOf(JSON.stringify(report));
    const file = fileFor(report, sha);
    await mount(report, file);
    const base = file.lineage.find((s) => s.scheme === "random")!;
    expect(document.querySelector('[data-baseline="lineage-header"]')?.textContent).toBe(`lineage ${pct(base.balanced_accuracy)}`);
    const random = report.splits.find((s) => s.scheme === "random")!;
    const leads = [...document.querySelectorAll('[data-baseline="lineage-lead"]')];
    const scored = random.learners.filter((l) => l.metrics !== null);
    expect(scored.length).toBeGreaterThan(0);
    expect(leads).toHaveLength(scored.length);
    for (const lead of leads) {
      const tr = lead.closest("tr")!;
      // The lead sits in the Random column, never under lineages held out.
      expect(lead.closest("td")).toBe(tr.querySelectorAll("td")[1]);
      const ba = Number.parseFloat(tr.querySelectorAll("td")[1]!.firstChild!.textContent!) / 100;
      const l = scored.find((x) => pct(x.metrics!.balanced_accuracy) === pct(ba))!;
      const want = signedPoints(leadOver(l.metrics!.balanced_accuracy, base.balanced_accuracy));
      expect(lead.textContent).toBe(`${want} points from lineage identity alone`);
      expect(want.startsWith(l.metrics!.balanced_accuracy > base.balanced_accuracy ? "+" : "−")).toBe(true);
    }
    const note = document.querySelector('[data-baseline="lineage-note"]')!.textContent!;
    expect(note).toContain(`reaches ${pct(base.balanced_accuracy)} balanced accuracy on the random split, where ${pct(7 / base.n)} of test genomes (7 of ${base.n})`);
    expect(note).toContain(`With lineages held out, 100.0% of test genomes (${base.n} of ${base.n}) have a lineage no training genome of their fold shares, so this baseline has nothing of the lineage to read there and no difference is shown.`);
  });

  it("shows nothing of a baseline when the file is absent, and the report renders", async () => {
    const report = reportFor("baselines-absent");
    await mount(report, null);
    expect(marks()).toEqual([]);
    const text = document.body.textContent ?? "";
    expect(text).toContain("Learner comparison");
    expect(text).toContain("Predicted MIC is within one doubling dilution");
    expect(text).not.toContain("constant MIC");
    expect(text).not.toContain("Lineage identity alone");
  });

  it("refuses a file pinned to another report's bytes, and the report still renders", async () => {
    const report = reportFor("baselines-other-hash");
    const other = await shaOf(JSON.stringify({ ...report, id: "another" }));
    await mount(report, fileFor(report, other));
    expect(marks()).toEqual([]);
    expect(document.body.textContent).toContain("Predicted MIC is within one doubling dilution");
    // The same file pinned to the bytes served is shown.
    const again = reportFor("baselines-other-hash-pinned");
    act(() => root?.unmount());
    await mount(again, fileFor(again, await shaOf(JSON.stringify(again))));
    expect(marks().length).toBeGreaterThan(0);
  });

  it("marks a dataset whose held-out predictions do not confirm the folds", async () => {
    const report = reportFor("baselines-unconfirmed");
    const sha = await shaOf(JSON.stringify(report));
    await mount(report, fileFor(report, sha, false));
    expect(document.querySelector('[data-baseline="mic"]')!.textContent).toContain("†");
    expect(document.body.textContent).toContain("† The held-out predictions do not confirm the rebuilt MIC folds for this dataset");
  });
});

// ---- the overview, the datasets table and the methods page --------------------

const constantsFor = (rows: readonly CatalogRow[], pick: (k: number) => { model: number; median: number; confirmed?: boolean } | null) => {
  const files = new Map(
    rows.flatMap((r, k) => {
      const p = pick(k);
      if (r.mic_essential_agreement === null || p === null) return [];
      const f = parseBaselines({ schema: 1, dataset: r.id, report_sha256: "c".repeat(64), evaluation_revision: 14, mic: micPart(100, p.model, p.median, p.median, p.confirmed ?? true), lineage: [] });
      return [[r.id, f] as const];
    }),
  );
  return micConstants(rows, files, 14);
};

describe("the overview's MIC finding beside the constant (synthetic files)", () => {
  // Every model figure a whole share of 100 exact MICs, so a file's model
  // figure is the row's own.
  const rows = manifest.datasets.map((d, k) => (d.mic_essential_agreement === null ? d : { ...d, mic_essential_agreement: (80 + (k % 15)) / 100 }));
  const m: Manifest = { ...manifest, evaluation_revision: 14, datasets: rows };
  const constants = constantsFor(rows, (k) => ({ model: 80 + (k % 15), median: 85 }));

  it("gives the median constant beside the model's median and the count not above it", () => {
    const b = micBaselineSummary(rows, constants);
    expect(b.covered).toBeGreaterThan(0);
    const f = micFinding(m, { page: true, micConstants: constants });
    expect(f.result).toContain(`the median agreement is ${pct(b.medianModel)} for the model and ${pct(b.medianConstant)} for the constant`);
    expect(f.result).toContain(`not above the constant on ${b.notAbove} of ${b.covered}.`);
    expect(f.caption).toContain("by a model fitted without that isolate (out of fold). The MIC folds are random and can put a k-mer-identical copy in training.");
    expect(f.caption).toContain("90% is an adapted research screen");
    expect(f.caption).toContain("the counts at 90% are a description read after the scores");
    const slide = micFinding(m, { page: false, micConstants: constants });
    expect(slide.result).toContain(`reaches a median of ${pct(b.medianConstant)} against the model's ${pct(b.medianModel)}`);
  });

  it("reads as before with no files", () => {
    const f = micFinding(m, { page: true });
    expect(f.result).toBe(micLead(micSummary(rows)));
    expect(f.caption).not.toContain("constant");
    expect(f.caption).toContain("as a research screen");
    expect(micFinding(m, { page: true, micConstants: new Map() }).result).not.toContain("constant");
  });
});

describe("the datasets table's MIC column (synthetic files)", () => {
  const rows = manifest.datasets.filter((d) => d.mic_essential_agreement !== null).map((d) => ({ ...d, mic_essential_agreement: 0.75 }));
  const cellOf = (constants?: ReadonlyMap<string, MicConstant>, k = 0) => {
    const col = datasetColumns(false, 14, constants).find((c) => c.id === "mic")!;
    return renderToStaticMarkup(createElement("div", null, col.cell(rows[k]!)));
  };

  it("puts the constant beneath the model's figure, marked when unconfirmed", () => {
    const confirmed = constantsFor(rows, () => ({ model: 75, median: 81 }));
    const html = cellOf(confirmed);
    expect(html).toContain("75.0%");
    expect(html).toContain('data-baseline="mic-constant"');
    expect(html).toContain("81.0%");
    expect(html).not.toContain("†");
    const unconfirmed = constantsFor(rows, () => ({ model: 75, median: 81, confirmed: false }));
    expect(cellOf(unconfirmed)).toContain("†");
    expect(datasetColumns(false, 14, confirmed).find((c) => c.id === "mic")!.tip!.some((t) => t.term === "Constant MIC (c)")).toBe(true);
  });

  it("shows the model alone with no file", () => {
    const html = cellOf(undefined);
    expect(html).toContain("75.0%");
    expect(html).not.toContain("data-baseline");
    expect(datasetColumns(false, 14).find((c) => c.id === "mic")!.tip).toHaveLength(1);
  });
});

describe("the methods page's 0.90 count (synthetic files)", () => {
  it("labels the counts as a description read after the scores, at an adapted research screen", () => {
    const rows = manifest.datasets.filter((d) => d.mic_essential_agreement !== null).map((d) => ({ ...d, mic_essential_agreement: 0.92 }));
    const b = micBaselineSummary(rows, constantsFor(rows, (k) => ({ model: 92, median: k % 2 === 0 ? 93 : 70 })));
    const s = constantAt90Sentence(b);
    expect(s).toContain("Read after the scores, as a description");
    expect(s).toContain(`the constant reaches 90% on ${b.constantAt90}, the model on ${b.modelAt90}, and both on ${b.bothAt90}`);
    expect(s).toContain("90% is an adapted research screen");
    expect(constantAt90Sentence(null)).toBe("");
  });
});
