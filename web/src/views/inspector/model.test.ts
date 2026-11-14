// The interpretable-model section's words about a rule term: where its
// k-mers sit, what a list of catalog families is, which named gene groups
// count for the drug's rule, how many features tie it, and whether its
// association holds inside lineages. Every report here is the synthetic
// fixture's, with cards, gene-group rows and verdict tables set by hand.

import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../../data/validate";
import type { DatasetReport, FeatureCard, GroupRow } from "../../data/schema";
import type { LineageControl } from "../../data/derive/lineage";
import { ModelSection } from "./ModelSection";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const several = manifest.datasets.find((d) => d.scm_rules.length >= 2 && d.catalog !== null)!;
const base: DatasetReport = { ...parseDataset(fixtureFile(`/data/datasets/${several.id}.json`)), input_fingerprint: "f".repeat(16), evaluation: { revision: 14, valid: true, reasons: [] } };
const first = base.models.scm.rules[0]!.feature;

const text = (html: string): string =>
  html
    .replace(/<[^>]+>/g, "")
    .replace(/&#x27;/g, "'")
    .replace(/&amp;/g, "&");
const render = (r: DatasetReport, lineageTable?: LineageControl) => text(renderToStaticMarkup(createElement(ModelSection, lineageTable === undefined ? { r } : { r, lineageTable })));
/** The report with its first rule's card changed. */
const withCard = (change: (card: FeatureCard) => FeatureCard, more: Partial<DatasetReport> = {}): DatasetReport => ({ ...base, ...more, features: base.features.map((f) => (f.feature === first ? change(f) : f)) });
/** The same, with every other card's catalog families cleared, so the page's words about families are the first rule's. */
const alone = (r: DatasetReport): DatasetReport => ({ ...r, features: r.features.map((f) => (f.feature === first ? f : { ...f, catalog_families: [] })) });
const group = (name: string, counts: boolean): GroupRow => ({ group: name, class: "", subclass: "", n_resistant: 1, n_susceptible: 1, counts_for_rule: counts });
const NOTE = "Not a list of genes an isolate carries: one k-mer pattern matches all of these families.";
const CAVEAT = "A rule names k-mers that separate the lab calls in these isolates; they can travel with the cause on a plasmid or in a lineage without being the cause.";

describe("a rule term's place", () => {
  it("heads a rule whose k-mers sit outside every gene with their distance from the nearest one", () => {
    const outside = withCard((c) => ({ ...c, catalog_families: [], locus: { ...c.locus!, gene: "mgrB", product: "PhoP/PhoQ regulator MgrB", codons: null, nearest: ["mgrB", 46] } }));
    const page = render(outside);
    expect(page).toContain(`Rule 1${base.models.scm.rules[0]!.presence ? "Presence" : "Absence"} of k-mers 46 bp from the nearest annotated gene, mgrB`);
    const inside = withCard((c) => ({ ...c, catalog_families: [], locus: { ...c.locus!, gene: "mgrB", product: "PhoP/PhoQ regulator MgrB", codons: [3, 12], nearest: null } }));
    expect(render(inside)).toContain(`Rule 1${base.models.scm.rules[0]!.presence ? "Presence" : "Absence"} of mgrB`);
    expect(render(inside)).not.toContain("from the nearest annotated gene");
  });
});

describe("a rule term's catalog families", () => {
  const catalog = { ...base.catalog!, groups: [group("blaNDM", true), group("ble", false), group("blaOXA (OXA-48 family)", true), group("blaOXA (OXA-1 family)", false)] };

  it("says a list of several families is one pattern, and which of their gene groups count for the drug's rule", () => {
    const page = render(withCard((c) => ({ ...c, catalog_families: ["blaNDM", "ble", "blaAFM"] }), { catalog }));
    expect(page).toContain(NOTE);
    expect(page).toContain("Of the gene groups found in these genomes, blaNDM counts for this drug's known-gene rule, but ble does not.");
    // blaAFM has no gene-group row: nothing is said of it.
    expect(page).not.toContain("blaAFM counts");
    expect(page).not.toContain("blaAFM does not");
  });

  it("names the groups of a family that has one counting and one not", () => {
    const page = render(alone(withCard((c) => ({ ...c, catalog_families: ["blaOXA"] }), { catalog })));
    expect(page).toContain("Of the gene groups found in these genomes, blaOXA (OXA-48 family) counts for this drug's known-gene rule, but blaOXA (OXA-1 family) does not.");
    // One family is no list, so the note about lists stays out.
    expect(page).not.toContain(NOTE);
  });

  it("says nothing of the drug's rule where the report holds no row for the family, or no catalog block", () => {
    expect(render(alone(withCard((c) => ({ ...c, catalog_families: ["tet(A)"] }), { catalog })))).not.toContain("known-gene rule");
    expect(render(alone(withCard((c) => ({ ...c, catalog_families: ["blaNDM"] }), { catalog: null })))).not.toContain("known-gene rule");
    expect(render(alone(withCard((c) => ({ ...c, catalog_families: ["blaNDM"] }), { catalog })))).toContain("Of the gene groups found in these genomes, blaNDM counts for this drug's known-gene rule.");
  });
});

describe("the causal caveat and the within-lineage result", () => {
  const table = (verdicts: string[], over: Partial<LineageControl> = {}): LineageControl => ({
    release: "synthetic",
    evaluation_revision: 14,
    min_strata: 5,
    alpha: 0.05,
    terms: verdicts.length,
    counts: {},
    datasets: { [base.id]: { fingerprint: base.input_fingerprint!, terms: base.models.scm.rules.map((r, i) => [r.feature, r.presence, verdicts[i]!]) } },
    ...over,
  });
  const verdicts = base.models.scm.rules.map((_, i) => (i === 0 ? "seen" : "cannot_test"));
  const later = base.models.scm.rules.length === 2 ? "rule 2" : `rules ${base.models.scm.rules.slice(1, -1).map((_, i) => i + 2).join(", ")} and ${base.models.scm.rules.length}`;

  it("keeps the caveat and attaches each rule's result to it", () => {
    const page = render(base, table(verdicts));
    expect(page).toContain(
      `${CAVEAT} For rule 1 the association also holds inside lineages. For ${later} these genomes do not have enough mixed lineages to separate this from lineage. Each rule was tested inside this report's lineage clusters, on the isolates the earlier rules left undecided, and needs 5 clusters that hold both lab calls and both states of the rule (release synthetic). An association that holds inside lineages is still an association, not a cause.`,
    );
  });

  it("keeps the caveat alone on a report the table was not computed on", () => {
    // The bundled table holds the release's datasets; the fixture's are not in it.
    for (const page of [render(base), render({ ...base, input_fingerprint: "0".repeat(16) }, table(verdicts)), render(base, table(verdicts, { evaluation_revision: 15 }))]) {
      expect(page).toContain(CAVEAT);
      expect(page).not.toContain("inside lineages");
      expect(page).not.toContain("mixed lineages");
    }
  });

  it("says nothing of a rule whose verdict has no words", () => {
    const page = render(base, table(base.models.scm.rules.map((_, i) => (i === 0 ? "reversed" : "seen"))));
    expect(page).toContain(`For ${later} the association also holds inside lineages.`);
    expect(page).not.toContain("For rule 1 ");
  });
});
