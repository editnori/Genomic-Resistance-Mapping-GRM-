// @vitest-environment happy-dom

// Where the clindamycin note stands in the dataset inspector: once, right
// after the learner comparison, on a report the check read, and nowhere on
// any other report. A report must list as many samples as it counts, so the
// synthetic fixture's report cannot carry a real dataset's sample count;
// here the bundled entries of three datasets take the fixture's count
// instead. The inspector passes the note the SHA-256 of the report's bytes
// as it read them, so one of the three entries also takes the SHA-256 of
// the bytes this test serves; the other two keep the release's, which no
// synthetic report has. clindamycin.test.ts renders the notes against the
// bundled file itself.

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
const EF = "enterococcus_faecium__clindamycin";

// Hoisted above the imports, so it names the datasets itself.
vi.mock("../../data/derive/clindamycin", async (importOriginal) => {
  const real = await importOriginal<typeof import("../../data/derive/clindamycin")>();
  const { fixtureFile: fixture } = await import("../../dev/fixtures/files");
  const index = fixture("/data/index.json") as { datasets: { id: string; n: number }[] };
  const file = real.CLINDAMYCIN_LABELS!;
  const datasets = { ...file.datasets };
  for (const id of ["staphylococcus_aureus__clindamycin", "campylobacter_coli__clindamycin", "enterococcus_faecium__clindamycin"]) datasets[id] = { ...datasets[id]!, samples: index.datasets[0]!.n };
  // The bytes the first test serves, built as reportAs builds them below.
  const sa = "staphylococcus_aureus__clindamycin";
  const served = JSON.stringify({ ...(fixture(`/data/datasets/${index.datasets[0]!.id}.json`) as Record<string, unknown>), id: sa, evaluation: { revision: 14, valid: true, reasons: [] } });
  const hash = await globalThis.crypto.subtle.digest("SHA-256", new TextEncoder().encode(served).buffer as ArrayBuffer);
  datasets[sa] = { ...datasets[sa]!, report_sha256: Array.from(new Uint8Array(hash), (b) => b.toString(16).padStart(2, "0")).join("") };
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
  // The report is read, then its bytes are hashed, then the note is drawn.
  for (let i = 0; i < 3; i++)
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

  it("is absent when the report's bytes are not those the check read", async () => {
    // The id, the revision and the sample count are the entry's; only the bytes differ.
    await mount(EF, reportAs(EF, 14), { ...manifest, datasets: [{ ...row, id: EF }], evaluation_revision: 14 });
    expect(document.getElementById("ds-learners")).not.toBeNull();
    expect(document.querySelectorAll("[data-clindamycin-note]").length).toBe(0);
    expect(document.body.textContent).not.toContain("Clindamycin labels by study");
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
