import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { Manifest } from "../data/schema";
import { LEARNER_LABEL, pct } from "../data/format";
import type { Navigator } from "../app/router";
import { DatasetsView } from "./DatasetsView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const nav: Navigator = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
const render = (datasets: Manifest["datasets"], evaluation_revision?: number) =>
  renderToStaticMarkup(createElement(DatasetsView, { manifest: { ...manifest, datasets, evaluation_revision }, nav }));

/** Cell text of the first table row, in column order: antibiotic, species,
 *  isolates, best learner, random, lineages held out, gap, VME, ME, ... */
const firstRow = (html: string): string[] => {
  const body = html.slice(html.indexOf("<tbody"));
  return [...body.slice(0, body.indexOf("</tr>")).matchAll(/<td[^>]*>([\s\S]*?)<\/td>/g)].map((m) => m[1]!.replace(/<[^>]+>/g, "").trim());
};

describe("dataset catalog", () => {
  it("shows why the lineage split is missing, and the generic words only without a reason", () => {
    const base = manifest.datasets[0]!;
    const reason = "fewer than twice as many lineages as requested folds";
    const withReason = render([{ ...base, best_lineage: null, lineage_unavailable_reason: reason }]);
    const without = render([{ ...base, best_lineage: null, lineage_unavailable_reason: undefined }]);

    expect(withReason).toContain("Too few lineages");
    expect(withReason).toContain("No split");
    expect(without).toContain("No split");
    expect(without).not.toContain("Not available");
  });

  it("names no best learner and shows no score in a row where no learner calls every isolate", () => {
    const base = manifest.datasets.find((d) => d.best_lineage !== null)!;
    const cells = firstRow(render([{ ...base, best_learner: null, best_random: null, best_lineage: null }]));
    expect(cells).toHaveLength(10);
    expect(cells[1]).toContain("No eligible learner: No learner calls every isolate");
    expect([cells[2], cells[3], cells[5], cells[6]]).toEqual(["–", "–", "–", "–"]);
    const scored = firstRow(render([base]));
    expect(scored[1]).toBe(`${pct(base.best_random)}${LEARNER_LABEL[base.best_learner!]}`);
    expect(scored[2]).toBe(pct(base.best_lineage));
  });

  it("marks only the earlier-revision row provisional in a mixed rebuild", () => {
    const [a, b] = manifest.datasets;
    const audit = (revision: number) => ({ revision, valid: true, reasons: [], scope: "Internal cross-validation only." });
    const rows = [
      { ...a!, evaluation: audit(6) },
      { ...b!, evaluation: audit(5) },
    ];
    const mixed = render(rows, 6);
    const older = render(rows);

    expect(mixed.split(">Provisional<").length - 1).toBe(1);
    expect(older).not.toContain(">Provisional<");
  });
});
