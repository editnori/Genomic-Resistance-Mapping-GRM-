// @vitest-environment happy-dom

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseManifest } from "../../data/validate";
import { int } from "../../data/format";
import type { CatalogRow, DatasetReport, Manifest } from "../../data/schema";
import type { Navigator } from "../../app/router";
import { ThemeProvider } from "../../theme/ThemeProvider";
import { DatasetInspector } from "./DatasetInspector";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const manifest = parseManifest(fixtureFile("/data/index.json"));
const row = manifest.datasets[0]!;
const earlier = { revision: 5, valid: true, reasons: [], scope: "Internal cross-validation only." };
const report = { ...(fixtureFile(`/data/datasets/${row.id}.json`) as Record<string, unknown>), evaluation: earlier };
const nav: Navigator = { route: { view: "datasets", species: null, dataset: row.id }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };

let root: Root | null = null;

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

/** Mounts the inspector over a stubbed fetch of `body`. The report cache keys
 *  on the dataset path, so each test opens a different dataset. */
const mount = async (id: string, body: unknown, current: Manifest) => {
  vi.stubGlobal("fetch", async () => new Response(JSON.stringify(body), { headers: { "Content-Type": "application/json" } }));
  const host = document.createElement("div");
  document.body.append(host);
  const r = createRoot(host);
  root = r;
  await act(async () => {
    r.render(createElement(ThemeProvider, { children: createElement(DatasetInspector, { id, manifest: current, nav: { ...nav, route: { ...nav.route, dataset: id } } }) }));
  });
  await act(async () => {
    await new Promise((done) => setTimeout(done, 0));
  });
};

/** The header's status row (the audit chip and the facts; the section pills
 *  sit below it) and, with the Method and coverage fold open, the fold. */
const readInspector = async () => {
  const fold = [...document.querySelectorAll("button")].find((b) => b.textContent === "Method and coverage");
  await act(async () => fold?.click());
  return { header: document.querySelector("header p")?.textContent ?? "", text: document.body.textContent ?? "" };
};

describe("dataset inspector audit", () => {
  it("words a report from an earlier audit revision as provisional in the header and the notice", async () => {
    const current = { ...manifest, datasets: [{ ...row, evaluation: earlier }], evaluation_revision: 7 };
    await mount(row.id, report, current);
    const text = document.body.textContent ?? "";

    expect(text).toContain("does not pass the current evaluation audit");
    expect(text).toContain("current revision is 7");
    expect(text).not.toContain("Internal cross-validation only");
  });
});

describe("dataset inspector method fold", () => {
  /** The inspector's text with the Method and coverage fold open, for the
   *  dataset's fixture report with counts.k set to `k`, or removed. */
  const foldText = async (id: string, k: number | undefined): Promise<string> => {
    const raw = fixtureFile(`/data/datasets/${id}.json`) as Record<string, unknown> & { counts: Record<string, unknown> };
    const { k: _recorded, ...counts } = raw.counts;
    const served = { ...raw, counts: k === undefined ? counts : { ...counts, k } };
    vi.stubGlobal("fetch", async () => new Response(JSON.stringify(served), { headers: { "Content-Type": "application/json" } }));
    act(() => root?.unmount());
    const host = document.createElement("div");
    document.body.replaceChildren(host);
    const r = createRoot(host);
    root = r;
    await act(async () => {
      r.render(createElement(ThemeProvider, { children: createElement(DatasetInspector, { id, manifest, nav }) }));
    });
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });
    const fold = [...document.querySelectorAll("button")].find((b) => b.textContent === "Method and coverage");
    act(() => fold!.click());
    return document.body.textContent ?? "";
  };

  it("gives the k its report records, and says when an older report records none", async () => {
    const recorded = manifest.datasets[4]!;
    const older = manifest.datasets[5]!;
    expect(await foldText(recorded.id, 21)).toContain("(k = 21)");
    expect(await foldText(older.id, undefined)).toContain("(k not recorded)");
  });
});

describe("dataset inspector MIC-inferred labels", () => {
  const open = async (d: CatalogRow, edit: (counts: DatasetReport["counts"]) => DatasetReport["counts"], rowEdit: (r: CatalogRow) => CatalogRow) => {
    const r = fixtureFile(`/data/datasets/${d.id}.json`) as DatasetReport;
    await mount(d.id, { ...r, counts: edit(r.counts) }, { ...manifest, datasets: [rowEdit(d)] });
    return { ...(await readInspector()), counts: edit(r.counts) };
  };

  it("counts the labels left out over the counted genomes with a label, and names the class minimum", async () => {
    const { header, text, counts } = await open(manifest.datasets[1]!, (c) => ({ ...c, labels_from_mic_excluded: 12 }), (r) => ({ ...r, labels_from_mic_excluded: 12 }));

    expect(header).toContain("12 MIC-inferred labels left out");
    expect(text).toContain(`12 of the ${int(counts.samples + 12)} counted genomes with a resistant or susceptible label`);
    expect(text).toContain("25 resistant and 25 susceptible isolates with a reported call");
    expect(text).toContain("Genomes whose label was inferred from MIC at the cohort breakpoint are left out; no learner or MIC regression trains or is scored on them.");
  });

  it("reads a report without the newer counts at the old class minimum of 25 and prints no zero for the left-out labels", async () => {
    const strip = ({ labels_from_mic_excluded: _left, min_class: _min, ...rest }: DatasetReport["counts"]) => rest;
    const stripRow = ({ labels_from_mic_excluded: _left, min_class: _min, ...rest }: CatalogRow) => rest;
    const { header, text } = await open(manifest.datasets[2]!, strip, stripRow);

    expect(header).not.toContain("MIC");
    expect(text).toContain("Class minimum25 resistant and 25 susceptible isolates with a reported call");
    expect(text.split("Not recorded in this report").length - 1).toBe(1);
    expect(text).toContain("MIC-inferred labels left outNot recorded in this report");
    expect(text).not.toContain("Labels inferred from MIC");
    expect(text).toContain("no learner or MIC regression trains or is scored on them");
  });

  it("keeps a pre-audit report's MIC-labeled samples in view", async () => {
    const legacy = ({ labels_from_mic_excluded: _left, min_class: _min, ...rest }: DatasetReport["counts"]) => ({ ...rest, labels_from_mic: 7 });
    const legacyRow = ({ labels_from_mic_excluded: _left, min_class: _min, ...rest }: CatalogRow) => ({ ...rest, labels_from_mic: 7 });
    const { header, text, counts } = await open(manifest.datasets[3]!, legacy, legacyRow);

    expect(header).toContain("7 labels inferred from MIC");
    expect(text).toContain(`Samples labeled from MIC7 of ${int(counts.samples)}`);
    expect(text).not.toContain("MIC-inferred labels left out");
    // Its samples include MIC-inferred labels, so the minimum claims no reported calls.
    expect(text).toContain("Class minimum25 resistant and 25 susceptible isolates");
    expect(text).not.toContain("isolates with a reported call");
  });
});
