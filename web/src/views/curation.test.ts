import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { int } from "../data/format";
import type { CatalogRow, Manifest } from "../data/schema";
import type { Navigator } from "../app/router";
import { corpusFinding, CorpusFigure } from "./findings";
import { MethodsView } from "./MethodsView";
import { OverviewView } from "./OverviewView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const nav: Navigator = { route: { view: "overview", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };

/** Rows as a revision-6 manifest writes them, without the class minimum or the left-out count. */
const unrecorded = ({ min_class: _min, labels_from_mic_excluded: _left, ...rest }: CatalogRow): CatalogRow => rest;
const at = (min_class: number): Manifest => ({ ...manifest, datasets: manifest.datasets.map((d) => ({ ...d, min_class })) });
const revision6: Manifest = { ...manifest, datasets: manifest.datasets.map(unrecorded) };

const views = (m: Manifest) => ({
  overview: renderToStaticMarkup(createElement(OverviewView, { manifest: m, nav })),
  methods: renderToStaticMarkup(createElement(MethodsView, { manifest: m })),
  slide: corpusFinding(m, { page: false }).caption,
});

describe("class minimum wording", () => {
  it("states the minimum the reports record in Overview, Methods, and the corpus slide", () => {
    const { overview, methods, slide } = views(at(10));
    const rule = "A report needs 10 resistant and 10 susceptible selected genomes with unconflicted reported R/S calls.";
    expect(overview).toContain(rule);
    expect(methods).toContain(`A dataset is one species and one antibiotic. ${rule}`);
    expect(slide).toContain(rule);
    for (const text of [overview, methods, slide]) expect(text).not.toContain("25 resistant");
  });

  it("reads rows without a minimum as the old learn gate of 25", () => {
    const { overview, methods, slide } = views(revision6);
    const rule = "A report needs 25 resistant and 25 susceptible selected genomes with unconflicted reported R/S calls.";
    expect(overview).toContain(rule);
    expect(methods).toContain(`A dataset is one species and one antibiotic. ${rule}`);
    expect(slide).toContain(rule);
    for (const text of [overview, methods, slide]) expect(text).not.toMatch(/records? the minimum/);
  });

  it("lists each minimum when the reports differ", () => {
    const [first, ...rest] = manifest.datasets;
    for (const others of [rest, rest.map(unrecorded)]) {
      const { overview } = views({ ...manifest, datasets: [{ ...first!, min_class: 10 }, ...others] });
      expect(overview).toContain(`at least 10 of each in 1 report and at least 25 of each in ${rest.length} reports`);
    }
  });
});

describe("MIC-inferred labels", () => {
  const f = manifest.funnel;

  it("leaves them out of the reported pairs in Overview, Methods, and on the corpus slide", () => {
    expect(f.pairs_labeled_from_mic).toBeGreaterThan(0);
    const { overview, methods } = views(manifest);
    expect(overview).toContain(`reported R/S genome–antibiotic pairs, from ${int(f.records)} BV-BRC records`);
    expect(overview).toContain(`>${int(f.pairs_reported)}<`);
    expect(methods).toContain(`${int(f.pairs)} pairs = ${int(f.pairs_reported)} reported + ${int(f.pairs_labeled_from_mic)} from MIC`);
    const slide = renderToStaticMarkup(createElement(CorpusFigure, { manifest }));
    expect(slide).toContain("Reported R/S genome–antibiotic pairs");
    expect(slide).toContain(`>${int(f.pairs_reported)}<`);
    for (const text of [overview, slide]) expect(text).not.toContain(`>${int(f.pairs_reported + f.pairs_labeled_from_mic)}<`);
  });

  it("explains in Methods that no model trains or is scored on them, and which reports count them", () => {
    const { methods } = views(manifest);
    expect(methods).toContain("Kept in the curated records; no model trains or is scored on these labels.");
    expect(methods).toContain("no model trains or is scored on that label. A dataset report counts the labels it left out in its Method and coverage fold; a report written before that count was added shows it as not recorded.");
    expect(methods).toContain("fitted on the MICs of genomes with a reported call; a genome whose label was inferred from its MIC is left out");
  });
});
