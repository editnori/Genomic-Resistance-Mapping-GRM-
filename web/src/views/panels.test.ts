import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { speciesShort } from "../data/format";
import { withRealProbeRows } from "../dev/fixtures/probeRows";
import type { Lite, Manifest } from "../data/schema";
import { PanelEvidence } from "./PanelEvidence";
import { MethodsView } from "./MethodsView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const render = (m: Manifest) => renderToStaticMarkup(createElement(PanelEvidence, { manifest: m }));

describe("panel evidence and unscored assays", () => {
  it("shows real probe rows with full totals, a separate no-call count, and called-isolate shares", () => {
    const html = render(withRealProbeRows(manifest));
    expect(html).toContain(">794<");
    expect(html).toContain(">897<");
    expect(html).toContain("788 of 794 (99.2%)");
    expect(html).toContain("792 of 889 (89.1%)");
    expect(html).toContain("783 of 897 (87.3%)");
    expect(html).toContain("171 of 1,488");
    expect(html).toContain("6 of 167");
    expect(html.split("8 resistant, 4 susceptible").length - 1).toBe(2);
    expect(html).toContain("93.9%");
    expect(html).not.toContain("count as targets not found");
  });

  it("lists each unsupported assay once with its datasets and reason", () => {
    const reason = "In E. coli, synonymous base changes read as mutations.";
    const row = { panel: "ecoli_qrdr", name: "gyrA and parC QRDR probes", dataset: "escherichia_coli__ciprofloxacin", species: "Escherichia coli", antibiotic: "ciprofloxacin", reason };
    const html = renderToStaticMarkup(createElement(MethodsView, { manifest: { ...manifest, unsupported_panels: [row, { ...row, dataset: "escherichia_coli__levofloxacin", antibiotic: "levofloxacin" }] } }));
    expect(html).toContain("Assays not scored");
    expect(html).toContain("1 assay has no result on the datasets listed");
    expect(html).toContain(`Ciprofloxacin in ${speciesShort("Escherichia coli")}, Levofloxacin in ${speciesShort("Escherichia coli")}`);
    expect(html.split(`Not scored. ${reason}`).length - 1).toBe(1);
  });

  it("keeps each reason of one assay apart", () => {
    const organism = "the panel does not cover this organism: the assay reports none of its targets for it";
    const drug = "no target of this panel is interpreted for this drug in this organism";
    const row = { panel: "gram_negative_bcid", name: "Gram-negative blood-culture panel", dataset: "streptococcus_pneumoniae__meropenem", species: "Streptococcus pneumoniae", antibiotic: "meropenem", reason: organism };
    const rows = [row, { ...row, dataset: "pseudomonas_aeruginosa__colistin", species: "Pseudomonas aeruginosa", antibiotic: "colistin", reason: drug }];
    const html = renderToStaticMarkup(createElement(MethodsView, { manifest: { ...manifest, unsupported_panels: rows } }));
    expect(html).toContain("1 assay has no result on the datasets listed");
    expect(html).toContain(`Not scored. ${organism.charAt(0).toUpperCase()}${organism.slice(1)}.`);
    expect(html).toContain(`Not scored. ${drug.charAt(0).toUpperCase()}${drug.slice(1)}.`);
    expect(html).toContain(`Meropenem in ${speciesShort("Streptococcus pneumoniae")}`);
    expect(html).toContain(`Colistin in ${speciesShort("Pseudomonas aeruginosa")}`);
  });

  it("says how many datasets with a panel have no best learner instead of dropping them silently", () => {
    const base = withRealProbeRows(manifest);
    const none = { ...base.datasets[0]!, id: "x__none", best_learner: null, best_random: null, best_lineage: null, panels: [["carbapenemase", base.panels[0]!.metrics] as [string, Lite]] };
    const html = render({ ...base, datasets: [...base.datasets, none] });
    expect(html).toContain("1 dataset with a simulated panel has no learner scored on every isolate and is left out.");
    expect(html).not.toContain("undefined");
    expect(render(base)).not.toContain("has no learner scored on every isolate");
  });

  it("adds no unsupported section to a manifest without unsupported assays", () => {
    const html = renderToStaticMarkup(createElement(MethodsView, { manifest }));
    expect(html).not.toContain("Assays not scored");
    expect(html).not.toContain("Not scored.");
  });
});
