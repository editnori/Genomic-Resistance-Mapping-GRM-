// @vitest-environment happy-dom

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseCorpus, parseManifest, parseSpecies } from "../data/validate";
import { capitalize, int, pct } from "../data/format";
import type { Manifest } from "../data/schema";
import type { Navigator } from "../app/router";
import { ThemeProvider } from "../theme/ThemeProvider";
import { drugTotals, drugsByCalls, lineageBlocks, rsOf } from "../charts/lineageBlocks";
import { MIN_JOINT, crossPairs, pairSummary } from "../charts/crossPairs";
import { SpeciesView, drugResult } from "./SpeciesView";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const manifest = parseManifest(fixtureFile("/data/index.json"));
const corpus = parseCorpus(fixtureFile("/data/corpus.json"));
const navFor = (species: string | null): Navigator => ({ route: { view: "species", species, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} });

let root: Root | null = null;

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

const render = async (m: Manifest, nav: Navigator) => {
  vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path)), { headers: { "Content-Type": "application/json" } }));
  const host = document.createElement("div");
  document.body.append(host);
  const r = createRoot(host);
  root = r;
  await act(async () => {
    r.render(createElement(ThemeProvider, { children: createElement(SpeciesView, { manifest: m, nav }) }));
  });
  await act(async () => {
    await new Promise((done) => setTimeout(done, 0));
  });
  return host;
};

const first = manifest.species.find((s) => s.genomes_analyzed > 0)!;
const report = parseSpecies(fixtureFile(`/data/species/${first.slug}.json`));

describe("species page hierarchy", () => {
  it("opens with one coverage statement from the manifest, the map drug's result, then a strip of corpus facts", async () => {
    const host = await render(manifest, navFor(null));
    const text = host.textContent ?? "";
    expect(text).toContain(
      `${int(first.genomes_analyzed)} genomes analyzed (${int(first.genomes_labeled)} carry any label, including MIC-inferred calls); ${int(first.datasets_reported)} of ${int(first.datasets_curated)} curated datasets have a report.`,
    );
    expect(text).toContain("Profiles, cross-resistance and call counts on this page include labels inferred from MIC, which dataset reports leave out.");
    const order = drugsByCalls(report);
    expect(text).toContain(drugResult(report, order[0]!));

    const facts = host.querySelector("section[aria-labelledby='sp-facts-title']")!;
    expect(facts.querySelector("h2")?.textContent).toBe("Corpus facts");
    const fact = (label: string) => [...facts.querySelectorAll("dt")].find((dt) => dt.textContent === label)?.nextElementSibling?.textContent;
    const row = corpus.by_species.find((r) => r.species === report.species)!;
    expect(fact("Resistant calls")).toBe(int(row.pairs_resistant));
    expect(fact("Susceptible calls")).toBe(int(row.pairs_susceptible));
    expect(fact("Lineages")).toBe(int(report.n_lineages));
    expect(fact("Catalog gene groups seen")).toBe(int(report.catalog_groups.length));
  });

  it("gives the antibiotic table its own full-width section, apart from the facts", async () => {
    const host = await render(manifest, navFor(null));
    const datasets = host.querySelector("#sp-datasets")!;
    expect(datasets.querySelector("h2")?.textContent).toBe("Antibiotic datasets");
    expect(datasets.querySelector("table")).not.toBeNull();
    expect(datasets.textContent).not.toContain("Corpus facts");
    expect(datasets.querySelector("dl")).toBeNull();
  });

  it("names a dataset with no best learner instead of showing a blank score", async () => {
    const slug = manifest.species.find((s) => s.genomes_analyzed > 0 && manifest.datasets.filter((d) => d.species_slug === s.slug).length > 1)!.slug;
    const [none, kept] = manifest.datasets.filter((d) => d.species_slug === slug);
    const changed = { ...manifest, datasets: manifest.datasets.map((d) => (d.id === none!.id ? { ...d, best_learner: null, best_random: null, best_lineage: null } : d)) };
    const host = await render(changed, navFor(slug));
    // Antibiotic, drug class, isolates, resistant, best random score, lineages held out.
    const cells = (antibiotic: string) =>
      [...host.querySelectorAll("#sp-datasets tbody tr")].find((tr) => tr.querySelector("td")?.textContent?.trim() === capitalize(antibiotic))!.querySelectorAll("td");
    expect(cells(none!.antibiotic)[4]!.textContent).toBe("No learner calls every isolate");
    expect(cells(kept!.antibiotic)[4]!.textContent).toBe(pct(kept!.best_random));
    expect(cells(kept!.antibiotic)[3]!.textContent).toBe(int(kept!.n_resistant));
  });
});

describe("species population structure", () => {
  it("draws lineage blocks by default, with every drug's R or S and I counts beside its row", async () => {
    const host = await render(manifest, navFor(null));
    const blocks = lineageBlocks(report);
    const grid = [...host.querySelectorAll('[role="grid"]')].find((g) => /lineage block/.test(g.getAttribute("aria-label") ?? ""))!;
    expect(grid.getAttribute("aria-label")).toContain(`${int(blocks.length)} lineage block${blocks.length === 1 ? "" : "s"} of ${int(report.genomes.length)} genomes by ${report.drugs.length} drugs`);
    expect(blocks.reduce((a, b) => a + (b.end - b.start), 0)).toBe(report.genomes.length);
    const labels = host.querySelector("#sp-structure")!.textContent ?? "";
    const anyI = report.drugs.some((_, j) => drugTotals(report, j).i > 0);
    report.drugs.forEach((d, j) => {
      const t = drugTotals(report, j);
      expect(labels).toContain(`${capitalize(d)}${int(rsOf(t))}${anyI ? int(t.i) : ""}`);
    });
  });

  it("offers tables whose counts match the chart: the map drug by block, then every genome", async () => {
    const host = await render(manifest, navFor(null));
    const table = [...host.querySelectorAll("button")].find((b) => b.textContent === "Table" && b.closest('[aria-label="Show population structure as"]') !== null)!;
    await act(async () => table.click());
    const order = drugsByCalls(report);
    const drug = order[0]!;
    const blocks = lineageBlocks(report);
    const blockTable = host.querySelector(`[role="region"][aria-label="${capitalize(report.drugs[drug]!)} calls by lineage block"]`)!;
    expect(blockTable.querySelectorAll("tbody tr")).toHaveLength(blocks.length);
    const totals = drugTotals(report, drug);
    const foot = [...blockTable.querySelectorAll("tfoot td")].map((td) => td.textContent);
    expect(foot).toEqual(["All genomes", `1–${int(report.genomes.length)}`, int(report.genomes.length), int(totals.r), int(totals.i), int(totals.s), int(totals.untested)]);
    // Summed body cells give the same totals.
    const sum = (col: number) => [...blockTable.querySelectorAll("tbody tr")].reduce((a, tr) => a + Number(tr.children[col]!.textContent!.replace(/,/g, "")), 0);
    expect([sum(3), sum(4), sum(5), sum(6)]).toEqual([totals.r, totals.i, totals.s, totals.untested]);

    const genomes = host.querySelector('[role="region"][aria-label="Every genome\'s calls"]')!;
    expect(genomes.querySelectorAll("tbody tr")).toHaveLength(report.genomes.length);
    const heads = [...genomes.querySelectorAll("thead th")].map((th) => th.textContent);
    const unit = (j: number) => {
      const t = drugTotals(report, j);
      return `${int(rsOf(t))} R or S${t.i > 0 ? ` · ${int(t.i)} I` : ""}`;
    };
    order.forEach((j) => expect(heads).toContain(`${capitalize(report.drugs[j]!)}${unit(j)}`));
    // The map drug's column holds exactly its counts.
    const col = heads.indexOf(`${capitalize(report.drugs[drug]!)}${unit(drug)}`);
    expect(col).toBeGreaterThan(0);
    const calls = [...genomes.querySelectorAll("tbody tr")].map((tr) => tr.children[col]!.textContent);
    expect([calls.filter((c) => c === "R").length, calls.filter((c) => c === "I").length, calls.filter((c) => c === "S").length, calls.filter((c) => c === "–").length]).toEqual([
      totals.r,
      totals.i,
      totals.s,
      totals.untested,
    ]);
  });
});

describe("species cross-resistance", () => {
  it("states how many pairs are evaluated and lists every pair in its table", async () => {
    const host = await render(manifest, navFor(null));
    const pairs = crossPairs(report);
    const s = pairSummary(pairs);
    const text = host.querySelector("#sp-cross")!.textContent ?? "";
    expect(text).toContain(`${int(s.evaluated)} of ${int(s.total)} pairs have at least ${MIN_JOINT} joint R/S calls`);
    if (s.highest !== null) expect(text).toContain(`${capitalize(s.highest.a)}–${s.highest.b} has highest phi: ${s.highest.phi!.toFixed(2)}`);
    expect(text).not.toContain("goes together most strongly");
    expect(host.querySelectorAll('[role="region"][aria-label="Cross-resistance pairs"] tbody tr')).toHaveLength(s.total);
  });
});
