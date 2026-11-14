// The within-lineage verdicts the inspector attaches to its causal caveat.
// The bundled table's counts are the ones reference/lineage-control gives
// for its release; that the bundled file is exactly what that table gives is
// the Python suite's check (tests/lineage_control/test_site_data.py). The
// reports below are the synthetic fixture's, with tables written by hand.

import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../validate";
import type { DatasetReport } from "../schema";
import { LINEAGE_CONTROL, lineageMethod, lineageSentences, lineageVerdicts, type LineageControl } from "./lineage";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const several = manifest.datasets.find((d) => d.scm_rules.length >= 2)!;
const base = parseDataset(fixtureFile(`/data/datasets/${several.id}.json`));
// A report the table below was computed on: revision 14 and one fingerprint.
const report: DatasetReport = { ...base, input_fingerprint: "f".repeat(16), evaluation: { revision: 14, valid: true, reasons: [] } };
const rules = report.models.scm.rules;

const tableFor = (r: DatasetReport, verdicts: string[], over: Partial<LineageControl> = {}): LineageControl => ({
  release: "synthetic",
  evaluation_revision: 14,
  min_strata: 5,
  alpha: 0.05,
  terms: verdicts.length,
  counts: {},
  datasets: { [r.id]: { fingerprint: r.input_fingerprint!, terms: r.models.scm.rules.map((rule, i) => [rule.feature, rule.presence, verdicts[i]!]) } },
  ...over,
});
const verdictsFor = rules.map((_, i) => (i === 0 ? "seen" : "cannot_test"));

describe("the bundled within-lineage table", () => {
  it("holds one release's terms, and its counts are the counts of its rows", () => {
    expect([LINEAGE_CONTROL.release, LINEAGE_CONTROL.evaluation_revision, LINEAGE_CONTROL.min_strata, LINEAGE_CONTROL.alpha]).toEqual(["rev14-274-20261112", 14, 5, 0.05]);
    expect(LINEAGE_CONTROL.counts).toEqual({ seen: 327, not_shown: 6, cannot_test: 529, reversed: 0 });
    const entries = Object.values(LINEAGE_CONTROL.datasets);
    expect(entries.length).toBe(274);
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
    expect(LINEAGE_CONTROL.terms).toBe(862);
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
  it("name each rule by its number and keep to association", () => {
    expect(lineageSentences(["seen"])).toEqual(["For rule 1 the association also holds inside lineages."]);
    expect(lineageSentences(["cannot_test"])).toEqual(["For rule 1 these genomes do not have enough mixed lineages to separate this from lineage."]);
    expect(lineageSentences(["not_shown"])).toEqual(["For rule 1 the test did not show the association inside lineages."]);
    expect(lineageSentences(["seen", "cannot_test", "seen", "not_shown", "seen", null])).toEqual([
      "For rules 1, 3 and 5 the association also holds inside lineages.",
      "For rule 4 the test did not show the association inside lineages.",
      "For rule 2 these genomes do not have enough mixed lineages to separate this from lineage.",
    ]);
    expect(lineageSentences([])).toEqual([]);
    expect(lineageSentences([null, null])).toEqual([]);
  });

  it("say how the test was run and that it shows no cause", () => {
    expect(lineageMethod()).toBe(
      "Each rule was tested inside this report's lineage clusters, on the isolates the earlier rules left undecided, and needs 5 clusters that hold both lab calls and both states of the rule (release rev14-274-20261112). The rules were chosen on these same isolates, so the test is not held-out evidence. An association that holds inside lineages is still an association, not a cause.",
    );
    expect(lineageMethod(tableFor(report, verdictsFor, { min_strata: 7, release: "other" }))).toContain("needs 7 clusters that hold both lab calls and both states of the rule (release other).");
  });
});
