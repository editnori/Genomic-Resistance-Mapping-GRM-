// @vitest-environment happy-dom

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { capitalize, pct } from "../data/format";
import type { Navigator } from "../app/router";
import { ThemeProvider } from "../theme/ThemeProvider";
import { SpeciesView } from "./SpeciesView";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const manifest = parseManifest(fixtureFile("/data/index.json"));
const nav: Navigator = { route: { view: "species", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };

let root: Root | null = null;

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

describe("species view MIC-inferred labels", () => {
  it("says its profiles, cross-resistance, and call counts include them", async () => {
    vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path)), { headers: { "Content-Type": "application/json" } }));
    const host = document.createElement("div");
    document.body.append(host);
    const r = createRoot(host);
    root = r;
    await act(async () => {
      r.render(createElement(ThemeProvider, { children: createElement(SpeciesView, { manifest, nav }) }));
    });
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });
    const text = document.body.textContent ?? "";

    expect(text).toContain("Corpus facts");
    expect(text).toContain("Profiles, cross-resistance, and the call counts below include labels inferred from MIC; dataset reports leave those out.");
  });
});

describe("species dataset table", () => {
  it("names a dataset with no best learner instead of showing a blank score", async () => {
    vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path)), { headers: { "Content-Type": "application/json" } }));
    const slug = manifest.species.find((s) => s.genomes_analyzed > 0 && manifest.datasets.filter((d) => d.species_slug === s.slug).length > 1)!.slug;
    const [none, kept] = manifest.datasets.filter((d) => d.species_slug === slug);
    const changed = { ...manifest, datasets: manifest.datasets.map((d) => (d.id === none!.id ? { ...d, best_learner: null, best_random: null, best_lineage: null } : d)) };
    const nav: Navigator = { route: { view: "species", species: slug, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
    const host = document.createElement("div");
    document.body.append(host);
    const r = createRoot(host);
    root = r;
    await act(async () => {
      r.render(createElement(ThemeProvider, { children: createElement(SpeciesView, { manifest: changed, nav }) }));
    });
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });
    // Antibiotic, isolates, best random score, lineages held out.
    const cells = (antibiotic: string) =>
      [...document.querySelectorAll("tbody tr")]
        .find((tr) => tr.querySelector("td")?.textContent?.trim() === capitalize(antibiotic))!
        .querySelectorAll("td");
    expect(cells(none!.antibiotic)[2]!.textContent).toBe("No learner calls every isolate");
    expect(cells(kept!.antibiotic)[2]!.textContent).toBe(pct(kept!.best_random));
  });
});
