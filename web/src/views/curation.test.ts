import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { int } from "../data/format";
import type { CatalogRow, Manifest } from "../data/schema";
import type { Navigator } from "../app/router";
import { corpusFinding, CorpusFigure } from "./findings";
import { curationPath, DataView } from "./DataView";
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
  data: renderToStaticMarkup(createElement(DataView, { manifest: m })),
  methods: renderToStaticMarkup(createElement(MethodsView, { manifest: m })),
  slide: corpusFinding(m, { page: false }).caption,
});

describe("class minimum wording", () => {
  it("states the minimum the reports record in Overview, Methods, and the corpus slide", () => {
    const { overview, data, slide } = views(at(10));
    const rule = "A report needs 10 resistant and 10 susceptible selected genomes with a reported call and no conflicting call.";
    expect(overview).toContain(rule);
    expect(data).toContain(`A dataset is one species and one antibiotic. ${rule}`);
    expect(data).toContain("under 10 resistant or under 10 susceptible genomes");
    expect(slide).toContain(rule);
    for (const text of [overview, data, slide]) expect(text).not.toContain("25 resistant");
  });

  it("reads rows without a minimum as the old learn gate of 25", () => {
    const { overview, data, slide } = views(revision6);
    const rule = "A report needs 25 resistant and 25 susceptible selected genomes with a reported call and no conflicting call.";
    expect(overview).toContain(rule);
    expect(data).toContain(`A dataset is one species and one antibiotic. ${rule}`);
    expect(slide).toContain(rule);
    for (const text of [overview, data, slide]) expect(text).not.toMatch(/records? the minimum/);
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
    const { overview, data } = views(manifest);
    expect(overview).toContain(`genome and antibiotic pairs with a reported R or S call, from ${int(f.records)} BV-BRC records`);
    expect(overview).toContain(`>${int(f.pairs_reported)}<`);
    expect(data).toContain(`${int(f.pairs)} pairs = ${int(f.pairs_reported)} reported + ${int(f.pairs_labeled_from_mic)} from MIC`);
    const steps = curationPath(manifest);
    expect(steps.map((s) => s.key)).toEqual([
      "records",
      "genomes_with_a_record",
      "genomes_passing",
      "pairs",
      "pairs_reported",
      "datasets_curated",
      "genomes_analyzed",
      "datasets_reported",
    ]);
    expect(steps.map((s) => s.count)).toEqual([
      f.records,
      f.genomes_with_records,
      f.genomes_passing_qc,
      f.pairs,
      f.pairs_reported,
      manifest.species.reduce((a, s) => a + s.datasets_curated, 0),
      manifest.species.reduce((a, s) => a + s.genomes_analyzed, 0),
      manifest.species.reduce((a, s) => a + s.datasets_reported, 0),
    ]);
    // Rev14-shaped manifest matching reference/release-summary/rev14-274-20261112/path.json:
    // verifies all 8 steps, the 2,000 cap note, 3 of 277 unreported datasets, and that
    // the 113,547 any-label sum and rebuilt-copy counts (80,967 / 57,206) stay off the path.
    const rev14Species = Array.from({ length: 23 }, (_, i) => ({
      species: `Species ${i + 1}`,
      slug: `species_${i + 1}`,
      datasets_curated: i === 0 ? 255 : 1,
      datasets_reported: i === 0 ? 252 : 1,
      genomes_labeled: i === 0 ? 69_547 : 2_000,
      genomes_analyzed: i < 8 ? 2_000 : i === 8 ? 757 : 500,
      lineages: 10,
      lineage_threshold: 0.005,
    }));
    const rev14Like: Manifest = {
      ...manifest,
      funnel: {
        records: 1_285_111,
        records_class_label: 8_709,
        records_not_antibiotic: 652,
        records_unresolved_name: 395,
        records_no_measure: 42_695,
        genomes_with_records: 123_503,
        genomes_without_metadata: 20,
        genomes_failed_qc: {
          "BV-BRC quality: Poor": 3_164,
          "CheckM completeness low": 633,
          "CheckM contamination high": 599,
          "species check: sequence is another species": 80,
          "genome length off species median": 14,
          "no species": 1,
        },
        genomes_passing_qc: 118_992,
        pairs: 1_063_628,
        pairs_reported: 611_826,
        pairs_labeled_from_mic: 295_219,
        pairs_unlabeled: 137_151,
        pairs_intermediate_only: 15_609,
        pairs_conflict: 3_823,
      },
      species: rev14Species,
    };
    const rev14Path = curationPath(rev14Like);
    expect(rev14Path.map((s) => s.count)).toEqual([1_285_111, 123_503, 118_992, 1_063_628, 611_826, 277, 23_757, 274]);
    expect(rev14Path.map((s) => s.unitMax !== null)).toEqual([false, false, true, false, true, false, true, true]);
    const rev14Html = renderToStaticMarkup(createElement(DataView, { manifest: rev14Like }));
    const pathBlock = rev14Html.slice(rev14Html.indexOf('id="m-path"'), rev14Html.indexOf('id="m-genomes"'));
    expect(pathBlock).toContain("4,491 genomes fail a check and 20 have no metadata");
    expect(pathBlock).toContain("9,756 records on any genome");
    expect(pathBlock).toContain("4,511 genomes that did not pass");
    expect(pathBlock).toContain("42,695 records on passing genomes");
    expect(pathBlock).toContain("295,219 pairs have a label worked out from an MIC");
    expect(pathBlock).toContain("8 of 23 species sit at 2,000");
    expect(pathBlock).toContain("3 of 277 datasets have one class under 25 among the genomes kept");
    expect(pathBlock).not.toContain("113,547");
    expect(rev14Html).toContain("The three name checks (9,756 records) read the records of every genome, passing or not;");
    expect(rev14Html).not.toContain("before QC");
    expect(rev14Html).not.toContain("80,967");
    expect(rev14Html).not.toContain("57,206");
    expect(rev14Html).not.toContain("1,064,823");
    const slide = renderToStaticMarkup(createElement(CorpusFigure, { manifest }));
    expect(slide).toContain("Genome and antibiotic pairs with a reported R or S call");
    expect(slide).toContain(`>${int(f.pairs_reported)}<`);
    for (const text of [overview, slide]) expect(text).not.toContain(`>${int(f.pairs_reported + f.pairs_labeled_from_mic)}<`);
  });

  it("explains in Methods that no model trains or is scored on them, and which reports count them", () => {
    const { data, methods } = views(manifest);
    expect(data).toContain("Kept in the curated records; no model trains or is scored on these labels.");
    expect(data).toContain("no model trains or is scored on that label. A dataset report counts the labels it left out in its Method and coverage fold; a report written before that count was added shows it as not recorded.");
    expect(data).toContain("A pair called susceptible and intermediate is labeled susceptible, and a pair called intermediate and resistant is labeled resistant.");
    const ruleBlock = data.slice(data.indexOf('id="m-label-rule"'), data.indexOf("</table>", data.indexOf('id="m-label-rule"')) + 8);
    expect([...ruleBlock.matchAll(/data-label-rule-row="([^"]+)"/g)].map((m) => m[1])).toEqual(["s", "r", "s-i", "i-r", "i", "conflict", "mic-only"]);
    expect(ruleBlock).toContain("Usable by models");
    expect(ruleBlock).toContain("S and I");
    expect(ruleBlock).toContain("I and R");
    expect(ruleBlock).toContain("R and S (with or without I)");
    expect(ruleBlock).toContain("Inferred R or S from its MIC, else none");
    expect(ruleBlock).toContain("No (curated records only)");
    expect(methods).toContain("fitted on the MICs of genomes with a reported call; a genome whose label was inferred from its MIC is left out");
  });
});
