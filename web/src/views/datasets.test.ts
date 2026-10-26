import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { Manifest } from "../data/schema";
import { LEARNER_LABEL, pct, speciesShort } from "../data/format";
import { BARS, barWords, verdictOf } from "../data/metrics";
import { CATALOG_CUTOFF, catalogSummary, splitSummary } from "../data/derive/overview";
import type { Navigator } from "../app/router";
import { catalogClose, datasetColumns, DatasetsView, largeGap } from "./DatasetsView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const nav: Navigator = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
const render = (datasets: Manifest["datasets"], evaluation_revision?: number) =>
  renderToStaticMarkup(createElement(DatasetsView, { manifest: { ...manifest, datasets, evaluation_revision }, nav }));

/** Raw markup of each cell of the first table row, in column order:
 *  dataset, random, lineage, gap, isolates, VME, ME, catalog, MIC, rule. */
const firstRowHtml = (html: string): string[] => {
  const body = html.slice(html.indexOf("<tbody"));
  return [...body.slice(0, body.indexOf("</tr>")).matchAll(/<td[^>]*>([\s\S]*?)<\/td>/g)].map((m) => m[1]!);
};
const text = (html: string): string => html.replace(/<[^>]+>/g, "").trim();
/** What a sighted reader sees: the cell without its screen-reader words. */
const visible = (html: string): string => text(html.replace(/<span class="sr-only">[\s\S]*?<\/span>/g, ""));
const firstRow = (html: string): string[] => firstRowHtml(html).map(text);

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

  it("prints a judged value as plain text with its verdict in words, and an absent value quietly", () => {
    const base = manifest.datasets.find((d) => {
      const lite = d.learners.find((l) => l.learner === d.best_learner)?.random;
      return lite?.vme != null && lite.me != null && d.mic_essential_agreement !== null && d.catalog !== null;
    })!;
    const lite = base.learners.find((l) => l.learner === base.best_learner)!.random!;
    const cells = firstRowHtml(render([base]));
    const judged = [
      [5, lite.vme!, BARS.vme],
      [6, lite.me!, BARS.me],
      [8, base.mic_essential_agreement!, BARS.mic],
    ] as const;
    for (const [i, value, bar] of judged) {
      expect(visible(cells[i]!)).toBe(pct(value));
      expect(text(cells[i]!)).toBe(`${pct(value)}, ${verdictOf(value, bar)} ${barWords(bar)}`);
      // No lozenge: no fill, ring, border or rounding on a value.
      expect(cells[i]).not.toMatch(/\b(ring-\S+|border\S*|bg-\S+|rounded\S*)\b/);
    }

    const absent = firstRowHtml(render([{ ...base, mic_essential_agreement: null, catalog: null }]));
    expect(visible(absent[8]!)).toBe("No MIC");
    expect(visible(absent[7]!)).toBe("No rule");
    expect(text(absent[8]!)).not.toMatch(/misses|meets/);
  });

  it("marks the same datasets the Overview's lineage and catalog sentences count, in words", () => {
    const rows = manifest.datasets;
    expect(rows.filter(largeGap)).toHaveLength(splitSummary(rows).nLargeGap);
    expect(rows.filter(catalogClose)).toHaveLength(catalogSummary(rows).within.get(CATALOG_CUTOFF)!);

    const base = rows.find((d) => d.best_learner !== null && d.catalog !== null)!;
    const wide = firstRowHtml(render([{ ...base, best_random: 0.9, best_lineage: 0.8, catalog: { ...base.catalog!, balanced_accuracy: 0.79 } }]));
    const close = firstRowHtml(render([{ ...base, best_random: 0.9, best_lineage: 0.86, catalog: { ...base.catalog!, balanced_accuracy: 0.7 } }]));
    expect(text(wide[3]!)).toBe("10.0, more than 5 points lower with lineages held out");
    expect(text(close[3]!)).toBe("4.0");
    expect(text(wide[7]!)).toBe("79.0%, within 2 points of the lineage score or ahead");
    expect(text(close[7]!)).toBe("70.0%");
    // Random and lineage scores stay neutral.
    expect(wide[1]).not.toContain("data-cue");
    expect(wide[2]).not.toContain("data-cue");
  });

  it("prints what the check, the coral gap and the blue catalog mean, and keeps each hint optional and short", () => {
    const page = text(render(manifest.datasets));
    expect(page).toContain("meets the bar under VME, ME or MIC");
    expect(page).toContain("Gap over 5 points");
    expect(page).toContain("Catalog within 2 points of Lineage, or ahead");

    for (const c of datasetColumns(false, undefined)) {
      const words = (c.tip ?? []).flatMap((l) => [l.term, l.text, l.source]).filter((s) => s !== undefined).join(" ");
      expect(words.split(/\s+/).filter((w) => w !== "").length).toBeLessThanOrEqual(40);
      // A judged column prints its bar in the header, not only in a hint.
      if (c.status !== undefined) expect(c.unit).toMatch(/^[≤≥]\d/);
    }
  });

  it("keeps the species whole with its separator and lets a combination drug break after its slash", () => {
    const base = manifest.datasets[0]!;
    const cell = firstRowHtml(render([{ ...base, antibiotic: "sulfamethoxazole/trimethoprim" }]))[0]!;
    expect(visible(cell)).toContain(`Sulfamethoxazole/trimethoprim${speciesShort(base.species)} · ${base.drug_class}`);
    expect(cell).toContain("/<wbr/>trimethoprim");
    // The separator closes the species' unbreakable run, so no line starts with it.
    expect(cell).toContain(`<i>${speciesShort(base.species)}</i> ·</span>`);
  });
});
