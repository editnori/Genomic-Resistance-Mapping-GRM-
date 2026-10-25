import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../../data/validate";
import type { DatasetReport, PanelResult } from "../../data/schema";
import { MechanismSection } from "./MechanismSection";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const withPanel = manifest.datasets.find((d) => d.panels.length > 0)!;
const report = parseDataset(fixtureFile(`/data/datasets/${withPanel.id}.json`));
const render = (r: DatasetReport) => renderToStaticMarkup(createElement(MechanismSection, { r }));
const count = (html: string, text: string) => html.split(text).length - 1;

describe("simulated panel wording", () => {
  it("tells a panel without a result apart from a panel whose targets are absent", () => {
    const unavailable = render({ ...report, panels: [] });
    const absent = render({ ...report, panels: [{ ...report.panels[0]!, detected: [] }] });

    expect(unavailable).toContain(`No simulated panel result is available for ${report.antibiotic}`);
    expect(unavailable).not.toContain("No simulated panel targets");
    expect(unavailable).not.toContain("none of its targets occurs here");
    expect(absent).toContain("none of its targets occurs here");
    expect(absent).not.toContain("No simulated panel result is available");
  });
});

describe("probe assay wording", () => {
  const probe: PanelResult = {
    ...report.panels[0]!,
    id: "gc_gyra",
    name: "gyrA codons 91-95 probes",
    mirrors: "gonorrhoea ciprofloxacin genotyping",
    detected: [["mutation detected", 2, 1]],
    metrics: { ...report.panels[0]!.metrics, n: 5, tp: 2, fn: 1, tn: 1, fp: 1 },
  };

  it("states the no-call count and the isolates the metrics cover", () => {
    const html = render({ ...report, panels: [{ ...probe, no_call: { resistant: 2, susceptible: 0 } }] });
    expect(html).toContain("The probes give no call for 2 resistant and 0 susceptible isolates");
    expect(html).toContain(`sensitivity and specificity cover 5 of ${report.counts.samples} isolates`);
    expect(html).toContain("The probes flag 2 of 3 resistant isolates and 1 of 2 susceptible ones");
  });

  it("says every isolate got a call when no isolate lacks one", () => {
    const html = render({ ...report, panels: [{ ...probe, no_call: { resistant: 0, susceptible: 0 } }] });
    expect(html).toContain(`Every one of the ${report.counts.samples} isolates gets a probe call.`);
    expect(html).not.toContain("give no call");
  });

  it("keeps a panel without no-call counts free of no-call wording", () => {
    const html = render(report);
    expect(html).not.toContain("probe call");
    expect(html).not.toContain("give no call");
  });

  it("lists an unsupported assay with its reason in place of metrics", () => {
    const reason = "In E. coli, synonymous base changes in gyrA codons 83-87 and parC codons 80-84 are common.";
    const html = render({
      ...report,
      panels: [],
      unsupported_assays: [{ id: "ecoli_qrdr", name: "gyrA and parC QRDR probes", mirrors: "Quinolone resistance-determining regions", reason }],
    });
    expect(html).toContain("gyrA and parC QRDR probes");
    expect(html).toContain(`Not scored. ${reason} No sensitivity or specificity is reported for this assay.`);
    expect(html).not.toContain("No simulated panel result is available");
    expect(count(html, "Sensitivity")).toBe(count(render({ ...report, panels: [] }), "Sensitivity"));
  });
});

describe("catalog baseline wording", () => {
  it("names GRM's catalog-containment rule and gives the AMRFinderPlus catalog only as its source", () => {
    const version = report.catalog!.version;
    const html = render(report);
    expect(html).toContain(`Catalog-containment rule · AMRFinderPlus catalog ${version}`);
    expect(html).not.toContain(`>AMRFinderPlus ${version}<`);
  });
});
