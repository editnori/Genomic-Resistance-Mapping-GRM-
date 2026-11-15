// The within-lineage verdicts the inspector attaches to its causal caveat.
// No count of the release is typed here: the bundled table's counts are
// checked against its own rows, and that the bundled file is exactly what
// the committed table gives is the Python suite's check
// (tests/lineage_control/test_site_data.py). The reports below are the
// synthetic fixture's, with tables written by hand.

import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../validate";
import type { DatasetReport } from "../schema";
import { int } from "../format";
import { LINEAGE_CONTROL, lineageMethod, lineageSentences, lineageTested, lineageVerdicts, type LineageControl } from "./lineage";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const several = manifest.datasets.find((d) => d.scm_rules.length >= 2)!;
const base = parseDataset(fixtureFile(`/data/datasets/${several.id}.json`));
// A report the table below was computed on: revision 14 and one fingerprint.
const report: DatasetReport = { ...base, input_fingerprint: "f".repeat(16), evaluation: { revision: 14, valid: true, reasons: [] } };
const rules = report.models.scm.rules;

/** The counts a table's verdicts give, as site_data.py writes them. */
const countsOf = (verdicts: readonly string[]): Record<string, number> => {
  const counts: Record<string, number> = { seen: 0, not_shown: 0, cannot_test: 0, reversed: 0 };
  for (const v of verdicts) counts[v] = (counts[v] ?? 0) + 1;
  return counts;
};

const tableFor = (r: DatasetReport, verdicts: string[], over: Partial<LineageControl> = {}): LineageControl => ({
  release: "synthetic",
  evaluation_revision: 14,
  min_strata: 5,
  alpha: 0.05,
  terms: verdicts.length,
  counts: countsOf(verdicts),
  datasets: { [r.id]: { fingerprint: r.input_fingerprint!, terms: r.models.scm.rules.map((rule, i) => [rule.feature, rule.presence, verdicts[i]!]) } },
  ...over,
});
const verdictsFor = rules.map((_, i) => (i === 0 ? "seen" : "cannot_test"));

describe("the bundled within-lineage table", () => {
  it("holds one release's terms, and its counts are the counts of its rows", () => {
    expect([LINEAGE_CONTROL.release, LINEAGE_CONTROL.evaluation_revision, LINEAGE_CONTROL.min_strata, LINEAGE_CONTROL.alpha]).toEqual(["rev14-274-20261112", 14, 5, 0.05]);
    expect(Object.keys(LINEAGE_CONTROL.counts).sort()).toEqual(["cannot_test", "not_shown", "reversed", "seen"]);
    const entries = Object.values(LINEAGE_CONTROL.datasets);
    expect(entries.length).toBeGreaterThan(0);
    const counted: Record<string, number> = { seen: 0, not_shown: 0, cannot_test: 0, reversed: 0 };
    for (const e of entries) {
      expect(e.fingerprint).toMatch(/^[0-9a-f]{16}$/);
      for (const [feature, presence, verdict] of e.terms) {
        expect([typeof feature, typeof presence]).toEqual(["number", "boolean"]);
        counted[verdict as string] = (counted[verdict as string] ?? 0) + 1;
      }
    }
    expect(counted).toEqual(LINEAGE_CONTROL.counts);
    expect(entries.reduce((n, e) => n + e.terms.length, 0)).toBe(LINEAGE_CONTROL.terms);
    expect(Object.values(LINEAGE_CONTROL.counts).reduce((n, c) => n + c, 0)).toBe(LINEAGE_CONTROL.terms);
  });

  it("gives the tested and passed rules from its own counts", () => {
    const { tested, passed } = lineageTested();
    // Every term but the ones that cannot be tested was tested, and `seen` is a pass.
    expect(tested).toBe(LINEAGE_CONTROL.terms - LINEAGE_CONTROL.counts.cannot_test!);
    expect(passed).toBe(LINEAGE_CONTROL.counts.seen);
    expect(tested).toBeGreaterThanOrEqual(passed);
    expect(passed).toBeGreaterThan(0);
    // A synthetic table: 4 seen, 2 not shown, 1 reversed and 5 that cannot be tested.
    const counts = { seen: 4, not_shown: 2, cannot_test: 5, reversed: 1 };
    expect(lineageTested({ ...LINEAGE_CONTROL, counts })).toEqual({ tested: 7, passed: 4 });
    expect(lineageTested({ ...LINEAGE_CONTROL, counts: { cannot_test: 9 } })).toEqual({ tested: 0, passed: 0 });
  });

  it("says nothing of a report it was not computed on", () => {
    // The fixture's datasets are synthetic and are not in the release.
    expect(lineageVerdicts(report)).toBeNull();
    expect(lineageVerdicts(base)).toBeNull();
  });
});

describe("a rule term's verdict", () => {
  it("is given only for the same dataset, revision, fingerprint and rules", () => {
    const table = tableFor(report, verdictsFor);
    expect(lineageVerdicts(report, table)).toEqual(verdictsFor);
    // Another dataset, a later evaluation revision, other inputs.
    expect(lineageVerdicts({ ...report, id: `${report.id}-other` }, table)).toBeNull();
    expect(lineageVerdicts({ ...report, evaluation: { revision: 15, valid: true, reasons: [] } }, table)).toBeNull();
    expect(lineageVerdicts({ ...report, evaluation: undefined }, table)).toBeNull();
    expect(lineageVerdicts({ ...report, input_fingerprint: "0".repeat(16) }, table)).toBeNull();
    expect(lineageVerdicts({ ...report, input_fingerprint: undefined }, table)).toBeNull();
    // A model with another rule at a position, the other sign, or another length.
    const scm = report.models.scm;
    const withRules = (next: typeof rules): DatasetReport => ({ ...report, models: { ...report.models, scm: { ...scm, rules: next } } });
    expect(lineageVerdicts(withRules(rules.map((r, i) => (i === 0 ? { ...r, feature: r.feature + 1 } : r))), table)).toBeNull();
    expect(lineageVerdicts(withRules(rules.map((r, i) => (i === 0 ? { ...r, presence: !r.presence } : r))), table)).toBeNull();
    expect(lineageVerdicts(withRules(rules.slice(1)), table)).toBeNull();
    // A model that is the table's model cut short: every rule it has matches
    // the table's rule at its position, and only the count of rules differs.
    expect(rules.length).toBeGreaterThanOrEqual(2);
    expect(lineageVerdicts(withRules(rules.slice(0, -1)), table)).toBeNull();
    // The other way round: the table holds fewer rules than the model.
    expect(lineageVerdicts(report, tableFor(withRules(rules.slice(0, -1)), verdictsFor.slice(0, -1)))).toBeNull();
    // A dataset named like an object's own method is still no entry.
    expect(lineageVerdicts({ ...report, id: "constructor" }, table)).toBeNull();
  });

  it("has no words for a verdict the site does not word", () => {
    const table = tableFor(report, rules.map((_, i) => (i === 0 ? "reversed" : "seen")));
    const verdicts = lineageVerdicts(report, table)!;
    expect(verdicts[0]).toBeNull();
    expect(lineageSentences(verdicts).join(" ")).not.toContain("rule 1 ");
  });
});

describe("the verdicts in sentences", () => {
  const SEEN = "the association is still seen when isolates are compared only inside their lineage cluster.";
  const NOT_SHOWN = "the test did not show the association inside lineages. That is not a finding of no association.";
  const CANNOT = "too few lineage clusters hold both lab calls and both states of the rule, so the test could not be run.";

  it("name each rule by its number and keep to association", () => {
    expect(lineageSentences(["seen"])).toEqual([`For rule 1 ${SEEN}`]);
    expect(lineageSentences(["cannot_test"])).toEqual([`For rule 1 ${CANNOT}`]);
    expect(lineageSentences(["not_shown"])).toEqual([`For rule 1 ${NOT_SHOWN}`]);
    expect(lineageSentences(["seen", "cannot_test", "seen", "not_shown", "seen", null])).toEqual([`For rules 1, 3 and 5 ${SEEN}`, `For rule 4 ${NOT_SHOWN}`, `For rule 2 ${CANNOT}`]);
    expect(lineageSentences([])).toEqual([]);
    expect(lineageSentences([null, null])).toEqual([]);
  });

  it("do not say the association holds in every lineage, and say what a test that did not pass does not show", () => {
    const all = lineageSentences(["seen", "not_shown", "cannot_test"]).join(" ");
    // One pooled test: nothing claims the association "holds" inside lineages.
    expect(all).not.toContain("holds inside");
    expect(all).not.toContain("also holds");
    expect(all).toContain("That is not a finding of no association.");
    // The sentence for a rule that could not be tested stands by itself.
    expect(all).not.toContain("separate this from lineage");
  });

  it("say how the test was run, how many rules passed by the table's own counts, and that it shows no cause", () => {
    const { tested, passed } = lineageTested();
    expect(lineageMethod()).toBe(
      `Each rule was tested inside this report's lineage clusters, on the isolates the earlier rules left undecided, and needs ${LINEAGE_CONTROL.min_strata} clusters that hold both lab calls and both states of the rule (release ${LINEAGE_CONTROL.release}). One test pools the clusters; it does not show the association in each one. Of the ${int(tested)} rules of this release that could be tested, ${int(passed)} passed (adjusted p under ${LINEAGE_CONTROL.alpha}). The clusters are MinHash clusters, not a phylogeny. The rules were chosen on these same isolates, so the test is not held-out evidence. An association seen inside lineage clusters is still an association, not a cause.`,
    );
    // Another table gives its own counts, level, minimum and release: none is fixed in the words.
    const other = lineageMethod(tableFor(report, verdictsFor, { min_strata: 7, release: "other", alpha: 0.01, counts: { seen: 1234, not_shown: 20, cannot_test: 50, reversed: 3 } }));
    expect(other).toContain("needs 7 clusters that hold both lab calls and both states of the rule (release other).");
    expect(other).toContain("One test pools the clusters; it does not show the association in each one. Of the 1,257 rules of this release that could be tested, 1,234 passed (adjusted p under 0.01). The clusters are MinHash clusters, not a phylogeny.");
    expect(lineageMethod(tableFor(report, ["seen"], { counts: { seen: 1, not_shown: 0, cannot_test: 4, reversed: 0 } }))).toContain("Of the 1 rule of this release that could be tested, 1 passed");
    expect(lineageMethod(tableFor(report, ["cannot_test"]))).toContain("One test pools the clusters; it does not show the association in each one. No rule of this release could be tested. The clusters are MinHash clusters");
  });
});
