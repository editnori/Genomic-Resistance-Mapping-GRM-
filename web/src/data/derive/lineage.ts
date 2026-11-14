// Whether a rule term's association with the lab calls holds inside
// lineages. scripts/lineage-control/verdicts.py tests every SCM rule term of
// one release inside its report's lineage clusters, and
// scripts/lineage-control/site_data.py writes data/lineage-control.json from
// its table (reference/lineage-control). A verdict is shown only on the
// report it was computed on; any other report gets none.

import bundled from "../lineage-control.json";
import type { DatasetReport } from "../schema";
import { int, listWords } from "../format";

/** `seen`: the association holds inside lineage clusters. `not_shown`: it
 *  could be tested and does not. `cannot_test`: too few clusters hold both
 *  lab calls and both states of the term. */
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
 *  Each sentence speaks of an association with the lab calls. */
export const lineageSentences = (verdicts: readonly (LineageVerdict | null)[]): string[] => {
  const at = (v: LineageVerdict) => verdicts.flatMap((x, i) => (x === v ? [i + 1] : []));
  const [seen, notShown, cannotTest] = [at("seen"), at("not_shown"), at("cannot_test")];
  return [
    ...(seen.length > 0 ? [`For ${rulesNamed(seen)} the association also holds inside lineages.`] : []),
    ...(notShown.length > 0 ? [`For ${rulesNamed(notShown)} the association is not shown inside lineages.`] : []),
    ...(cannotTest.length > 0 ? [`For ${rulesNamed(cannotTest)} these genomes do not have enough mixed lineages to separate this from lineage.`] : []),
  ];
};

/** How the verdicts were reached, and what they do not show. */
export const lineageMethod = (table: LineageControl = LINEAGE_CONTROL): string =>
  `Each rule was tested inside this report's lineage clusters, on the isolates the earlier rules left undecided, and needs ${int(table.min_strata)} clusters that hold both lab calls and both states of the rule (release ${table.release}). An association that holds inside lineages is still an association, not a cause.`;
