// @vitest-environment happy-dom

// The clindamycin note on the dataset report, the mark in the datasets
// table and the Methods paragraph, rendered; clindamycin-inspector.test.ts
// places the note in the dataset inspector. The six notes' text is typed
// out here, so a change of wording or of a count shows in review; that each
// count is the page's is the Python suite's check (tests/clindamycin_labels).
// The reports are the synthetic fixture's, renamed to one of the six with the
// bundled file's revision and sample count.

import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../../data/validate";
import type { DatasetReport, Manifest } from "../../data/schema";
import type { Navigator } from "../../app/router";
import { CLINDAMYCIN_LABELS } from "../../data/derive/clindamycin";
import { ClindamycinMark, ClindamycinMethods, ClindamycinNote } from "./ClindamycinNote";
import { DatasetsView } from "../DatasetsView";
import { MethodsView } from "../MethodsView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const row = manifest.datasets[0]!;
const base = parseDataset(fixtureFile(`/data/datasets/${row.id}.json`));
const file = CLINDAMYCIN_LABELS!;
const SA = "staphylococcus_aureus__clindamycin";
const CC = "campylobacter_coli__clindamycin";

const reportAs = (id: string, revision = file.evaluation_revision): DatasetReport => ({
  ...base,
  id,
  evaluation: { revision, valid: true, reasons: [] },
  counts: { ...base.counts, samples: file.datasets[id]!.samples },
});
const text = (html: string): string =>
  html
    .replace(/<\/p>/g, "\n")
    .replace(/<[^>]+>/g, "")
    .replace(/&#x27;/g, "'")
    .trim();

/** Each note as a reader sees it: the first paragraph, then the second. */
const NOTES: Record<string, string> = {
  campylobacter_coli__clindamycin:
    "Clindamycin labels by study: measured dependence. PMID 26519386 holds 32 of the dataset's 89 disagreements between GRM's known-gene rule and the reported phenotype, on 81 of its 389 samples. 294 of the 389 samples carry no PMID, and the study's share of the disagreements stands out only against those samples. Without that study, 9 of the 33 distinct figures tested fall by 0.05 or more with lineage-cluster intervals that exclude 0 (0.05 is an adapted research screen chosen for this check).\n" +
    "Those 9 figures are associated with that one study's labels or genomes. The other 24 do not pass the rule, and neither result shows which reading is wrong. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
  campylobacter_jejuni__clindamycin:
    "Clindamycin labels by study: class balance, dependence not measured. PMID 26519386 holds 17 of the dataset's 24 disagreements between GRM's known-gene rule and the reported phenotype, on 31 of its 508 samples. 475 of the 508 samples carry no PMID, and the study's share of the disagreements stands out only against those samples. Without it 8 reported resistant samples remain, under the class minimum of 25 this release uses, so the dataset cannot be scored without that study.\n" +
    "The dataset depends on that study for its class balance, and whether its figures depend on the study was not measured. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
  enterococcus_faecium__clindamycin:
    "Clindamycin labels by study: one study. Every sample comes from PMID 38219758 (901 of 901), and the known-gene rule and the reported phenotype disagree on 277 of them.\n" +
    "The figures describe that study's genomes and labels alone, and no figure without that study can exist. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
  staphylococcus_aureus__clindamycin:
    "Clindamycin labels by study: measured dependence. PMID 34100643 holds 182 of the dataset's 219 disagreements between GRM's known-gene rule and the reported phenotype, on 659 of its 1,543 samples. Without that study, 19 of the 33 distinct figures tested rise by 0.05 or more with lineage-cluster intervals that exclude 0 (0.05 is an adapted research screen chosen for this check).\n" +
    "Those 19 figures are associated with that one study's labels or genomes. The other 14 do not pass the rule, and neither result shows which reading is wrong. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
  streptococcus_agalactiae__clindamycin:
    "Clindamycin labels by study: class balance, dependence not measured. PMID 26135871 holds the dataset's single disagreement between GRM's known-gene rule and the reported phenotype, on 85 of its 141 samples. Without it 18 reported susceptible samples remain, under the class minimum of 25 this release uses, so the dataset cannot be scored without that study.\n" +
    "The dataset depends on that study for its class balance, and whether its figures depend on the study was not measured. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
  streptococcus_pneumoniae__clindamycin:
    "Clindamycin labels by study: no dependence under the rule. PMID 37823632 holds 18 of the dataset's 22 disagreements between GRM's known-gene rule and the reported phenotype, on 1,461 of its 1,894 samples. Without that study, none of the 33 distinct figures tested changes by 0.05 or more with a lineage-cluster interval that excludes 0 (0.05 is an adapted research screen chosen for this check).\n" +
    "This does not show that the figures do not depend on that study: some intervals are wide, and the fold plans of this dataset were not varied. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
};

describe("the note on a dataset report", () => {
  it("reads as typed here for each of the six datasets", () => {
    expect(Object.keys(NOTES).sort()).toEqual(Object.keys(file.datasets).sort());
    for (const [id, expected] of Object.entries(NOTES)) {
      expect(text(renderToStaticMarkup(createElement(ClindamycinNote, { r: reportAs(id) }))), id).toBe(expected);
    }
  });

  it("with a digest prop, renders only when it is the file's report_sha256; without one, as before", () => {
    const sha = file.datasets[SA]!.report_sha256;
    const note = (props: Record<string, unknown>) => renderToStaticMarkup(createElement(ClindamycinNote, { r: reportAs(SA), ...props }));
    const plain = note({});
    expect(plain).toContain('data-clindamycin-note="measured_dependence"');
    expect(note({ digest: sha })).toBe(plain);
    for (const digest of [undefined, null, file.datasets[CC]!.report_sha256]) expect(note({ digest }), String(digest)).toBe("");
  });

  it("renders nothing for a report outside the six or of another revision", () => {
    expect(renderToStaticMarkup(createElement(ClindamycinNote, { r: base }))).toBe("");
    expect(renderToStaticMarkup(createElement(ClindamycinNote, { r: reportAs(SA, file.evaluation_revision + 1) }))).toBe("");
    expect(renderToStaticMarkup(createElement(ClindamycinNote, { r: reportAs(SA), file: null }))).toBe("");
  });
});

describe("the mark in the datasets table", () => {
  const nav: Navigator = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
  const render = (datasets: Manifest["datasets"]) => renderToStaticMarkup(createElement(DatasetsView, { manifest: { ...manifest, datasets }, nav }));

  it("marks a row of the six with its kind in the tooltip, and no other row", () => {
    const mark = renderToStaticMarkup(createElement(ClindamycinMark, { row: { ...row, id: SA, n: file.datasets[SA]!.samples, evaluation: { revision: 14, valid: true, reasons: [] } } }));
    expect(mark).toContain('title="Clindamycin labels by study: measured dependence. The dataset report gives the counts behind it."');
    expect(text(mark.replace(/<span class="sr-only">[\s\S]*?<\/span>/g, ""))).toBe("Label check");
    const html = render([
      { ...row, id: SA, n: file.datasets[SA]!.samples, evaluation: { revision: 14, valid: true, reasons: [] } },
      { ...row, id: "made_up__drug", evaluation: { revision: 14, valid: true, reasons: [] } },
    ]);
    expect(html.match(/data-clindamycin-mark/g)?.length).toBe(1);
  });

  it("adds nothing to a table of other datasets", () => {
    expect(render(manifest.datasets)).not.toContain("data-clindamycin-mark");
  });
});

describe("the Methods paragraph", () => {
  it("stands in Labels and inclusion for the release the check read, and not for another", () => {
    const html = renderToStaticMarkup(createElement(MethodsView, { manifest: { ...manifest, evaluation_revision: 14 } }));
    const labels = html.slice(html.indexOf('id="m-labels"'), html.indexOf('id="m-models"'));
    expect(labels).toContain("data-clindamycin-methods");
    expect(text(renderToStaticMarkup(createElement(ClindamycinMethods, { manifest: { evaluation_revision: 14 } })))).toBe(
      "A check of release rev14-274-20261112 (evaluation revision 14) asked, for each of its 6 clindamycin datasets, whether one deposited study holds a large share of the samples on which GRM's known-gene rule and the reported phenotype disagree, and whether the released figures depend on that study. A disagreement is a difference between two readings, and neither is called the error. The page sorts the 6 datasets into four kinds: measured dependence (2 datasets); one study (1); class balance, dependence not measured (2); and no dependence under the rule (1). Each of their dataset reports gives its kind with the counts behind it. The check and its limits are in docs/clindamycin-labels.md in the GRM repository.",
    );
    expect(renderToStaticMarkup(createElement(MethodsView, { manifest: { ...manifest, evaluation_revision: 13 } }))).not.toContain("data-clindamycin-methods");
  });
});
