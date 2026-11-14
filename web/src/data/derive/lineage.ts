// Whether a rule term's association with the lab calls is still seen when
// isolates are compared only inside their lineage cluster.
// scripts/lineage-control/verdicts.py tests every SCM rule term of one
// release inside its report's lineage clusters, and
// scripts/lineage-control/site_data.py writes data/lineage-control.json from
// its table (reference/lineage-control). A verdict is shown only on the
// report it was computed on; any other report gets none.
//
// A verdict is one test that pools the clusters. It does not show the
// association in each cluster, the clusters are MinHash clusters and no
// phylogeny, and the rules were chosen on the same isolates. The sentences
// below say each of those things beside the verdict.

import bundled from "../lineage-control.json";
import type { DatasetReport } from "../schema";
import { int, listWords, plural } from "../format";

/** `seen`: the pooled test over the lineage clusters passed, with more
 *  resistant isolates than expected where the term holds. `not_shown`: the
 *  term could be tested and the test did not pass. `cannot_test`: too few
 *  clusters hold both lab calls and both states of the term. */
export type LineageVerdict = "seen" | "not_shown" | "cannot_test";

export interface LineageControl {
  release: string;
  evaluation_revision: number;
  /** Fewest informative lineage clusters a term needs to be tested. */
  min_strata: number;
  alpha: number;
  terms: number;
  counts: Record<string, number>;
  /** Per dataset: the report's input fingerprint, and each rule term in
   *  order as [feature, presence, verdict]. */
  datasets: Record<string, { fingerprint: string; terms: (string | number | boolean)[][] }>;
}

export const LINEAGE_CONTROL: LineageControl = bundled;

const WORDED: readonly string[] = ["seen", "not_shown", "cannot_test"];

/** Each rule term's verdict, in the model's order, or null when the table
 *  was not computed on this report: another dataset, evaluation revision or
 *  input fingerprint, or a model whose rules differ. A term whose verdict
 *  has no words here reads as null, so nothing is said of it. */
export const lineageVerdicts = (r: DatasetReport, table: LineageControl = LINEAGE_CONTROL): (LineageVerdict | null)[] | null => {
  const entry = Object.hasOwn(table.datasets, r.id) ? table.datasets[r.id] : undefined;
  const rules = r.models.scm.rules;
  if (entry === undefined || r.evaluation?.revision !== table.evaluation_revision || r.input_fingerprint !== entry.fingerprint) return null;
  if (entry.terms.length !== rules.length || !rules.every((rule, i) => entry.terms[i]![0] === rule.feature && entry.terms[i]![1] === rule.presence)) return null;
  return entry.terms.map(([, , verdict]) => (typeof verdict === "string" && WORDED.includes(verdict) ? (verdict as LineageVerdict) : null));
};

const rulesNamed = (positions: readonly number[]): string => `${positions.length === 1 ? "rule" : "rules"} ${listWords(positions.map((p) => int(p)))}`;

/** The verdicts in sentences, rules named by their number in the model.
 *  Each sentence speaks of an association with the lab calls. "Still seen"
 *  is the claim of one pooled test, not that the association holds in every
 *  lineage. A test that did not pass is followed by what it does not show,
 *  and a rule that could not be tested says why. */
export const lineageSentences = (verdicts: readonly (LineageVerdict | null)[]): string[] => {
  const at = (v: LineageVerdict) => verdicts.flatMap((x, i) => (x === v ? [i + 1] : []));
  const [seen, notShown, cannotTest] = [at("seen"), at("not_shown"), at("cannot_test")];
  return [
    ...(seen.length > 0 ? [`For ${rulesNamed(seen)} the association is still seen when isolates are compared only inside their lineage cluster.`] : []),
    ...(notShown.length > 0 ? [`For ${rulesNamed(notShown)} the test did not show the association inside lineages. That is not a finding of no association.`] : []),
    ...(cannotTest.length > 0 ? [`For ${rulesNamed(cannotTest)} too few lineage clusters hold both lab calls and both states of the rule, so the test could not be run.`] : []),
  ];
};

/** How many rule terms of the table's release could be tested, and how many
 *  of those passed, both read from the table's own counts: every verdict but
 *  `cannot_test` is a term that was tested, and `seen` is one that passed. */
export const lineageTested = (table: LineageControl = LINEAGE_CONTROL): { tested: number; passed: number } => ({
  tested: Object.entries(table.counts).reduce((n, [verdict, terms]) => (verdict === "cannot_test" ? n : n + terms), 0),
  passed: table.counts.seen ?? 0,
});

/** How the verdicts were reached, and what they do not show: one test pools
 *  the clusters, so it does not show the association in each one; how many
 *  of the release's testable rules passed, and at what level; the clusters
 *  are MinHash clusters, not a phylogeny; the test is in-sample; and an
 *  association is not a cause. Every count is the table's. */
export const lineageMethod = (table: LineageControl = LINEAGE_CONTROL): string => {
  const { tested, passed } = lineageTested(table);
  const outcome = tested === 0 ? "No rule of this release could be tested." : `Of the ${plural(tested, "rule")} of this release that could be tested, ${int(passed)} passed (adjusted p under ${table.alpha}).`;
  return `Each rule was tested inside this report's lineage clusters, on the isolates the earlier rules left undecided, and needs ${int(table.min_strata)} clusters that hold both lab calls and both states of the rule (release ${table.release}). One test pools the clusters; it does not show the association in each one. ${outcome} The clusters are MinHash clusters, not a phylogeny. The rules were chosen on these same isolates, so the test is not held-out evidence. An association seen inside lineage clusters is still an association, not a cause.`;
};
