// @vitest-environment happy-dom

// Where the clindamycin note stands in the dataset inspector: once, right
// after the learner comparison, on a report the check read, and nowhere on
// any other report. A report must list as many samples as it counts, so the
// synthetic fixture's report cannot carry a real dataset's sample count;
// here the bundled entry for one dataset takes the fixture's count instead.
// clindamycin.test.ts renders the notes against the bundled file itself.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseManifest } from "../../data/validate";
import type { Manifest } from "../../data/schema";
import type { Navigator } from "../../app/router";
import { ThemeProvider } from "../../theme/ThemeProvider";
import { DatasetInspector } from "./DatasetInspector";

const SA = "staphylococcus_aureus__clindamycin";
const CC = "campylobacter_coli__clindamycin";

// Hoisted above the imports, so it names the datasets itself.
vi.mock("../../data/derive/clindamycin", async (importOriginal) => {
  const real = await importOriginal<typeof import("../../data/derive/clindamycin")>();
  const { fixtureFile: fixture } = await import("../../dev/fixtures/files");
  const index = fixture("/data/index.json") as { datasets: { id: string; n: number }[] };
  const file = real.CLINDAMYCIN_LABELS!;
  const datasets = { ...file.datasets };
  for (const id of ["staphylococcus_aureus__clindamycin", "campylobacter_coli__clindamycin"]) datasets[id] = { ...datasets[id]!, samples: index.datasets[0]!.n };
  return { ...real, CLINDAMYCIN_LABELS: { ...file, datasets } };
});

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const manifest = parseManifest(fixtureFile("/data/index.json"));
const row = manifest.datasets[0]!;
const nav: Navigator = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };

let root: Root | null = null;

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

/** Mounts the inspector over a stubbed fetch of `body`; the report cache keys on the path. */
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

/** The fixture's first report under another id and evaluation revision. */
const reportAs = (id: string, revision: number): Record<string, unknown> => ({ ...(fixtureFile(`/data/datasets/${row.id}.json`) as Record<string, unknown>), id, evaluation: { revision, valid: true, reasons: [] } });

describe("the clindamycin note in the dataset inspector", () => {
  it("stands once, right after the learner comparison, on a report the check read", async () => {
    await mount(SA, reportAs(SA, 14), { ...manifest, datasets: [{ ...row, id: SA }], evaluation_revision: 14 });
    const notes = document.querySelectorAll("[data-clindamycin-note]");
    expect(notes.length).toBe(1);
    expect(notes[0]!.getAttribute("data-clindamycin-note")).toBe("measured_dependence");
    expect(notes[0]!.previousElementSibling?.id).toBe("ds-learners");
    expect(notes[0]!.nextElementSibling?.id).toBe("ds-model");
    expect(notes[0]!.textContent).toContain("Clindamycin labels by study: measured dependence.");
  });

  it("is absent when the same dataset's report is of another evaluation revision", async () => {
    // The report cache keys on the dataset path, so this test opens another of the six.
    await mount(CC, reportAs(CC, 13), { ...manifest, datasets: [{ ...row, id: CC }], evaluation_revision: 13 });
    expect(document.getElementById("ds-learners")).not.toBeNull();
    expect(document.querySelectorAll("[data-clindamycin-note]").length).toBe(0);
  });

  it("is absent from a report outside the six, which renders as before", async () => {
    await mount(row.id, fixtureFile(`/data/datasets/${row.id}.json`), manifest);
    expect(document.getElementById("ds-learners")).not.toBeNull();
    expect(document.querySelectorAll("[data-clindamycin-note]").length).toBe(0);
    expect(document.body.textContent).not.toContain("Clindamycin labels by study");
  });
});
