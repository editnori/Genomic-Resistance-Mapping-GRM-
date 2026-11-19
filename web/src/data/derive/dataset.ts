// Per-dataset derivations for the inspector. Pure functions of one
// DatasetReport.

import type {
  BinaryMetrics,
  CatalogBaseline,
  DatasetReport,
  FeatureCard,
  Learner,
  LearnerResult,
  Manifest,
  PanelResult,
  Rule,
  Scheme,
  SplitResult,
} from "../schema";
import { decodeBits } from "../bits";
import { capitalize, int, learnerInSentence, listWords, pct } from "../format";
import { CATALOG_CUTOFF } from "./overview";
import { EQUAL_SCORING_REVISION } from "./evaluation";

/** The best learner in words, for every view that names it. report.rs and
 *  bestLearner rank by random-split balanced accuracy, and only learners
 *  whose rates cover every isolate on every split that ran are candidates:
 *  from revision 14 all of them, earlier only those that call every isolate. */
export const BEST_LEARNER_WORDS = "the best random-split learner among those scored on every isolate";

/** A catalog row with no best learner: every learner leaves some isolate
 *  without a call (reports before revision 14). */
export const NO_BEST_LEARNER = "No learner is scored on every isolate";

/** A learner row with at least one call, so it has metrics. */
export type CalledLearner = LearnerResult & { metrics: BinaryMetrics };

/** A split that produced results. A split the pipeline planned but could
 *  not evaluate (no learners, an unavailable_reason) reads as absent here. */
export const splitOf = (r: DatasetReport, scheme: Scheme): SplitResult | undefined =>
  r.splits.find((s) => s.scheme === scheme && s.learners.length > 0);

/** Why a scheme has no results: the pipeline's stated reason, the lineage
 *  floor when the split never ran, or null when it has results. */
export const splitAbsence = (r: DatasetReport, scheme: Scheme): string | null => {
  if (splitOf(r, scheme) !== undefined) return null;
  const planned = r.splits.find((s) => s.scheme === scheme);
  if (planned?.unavailable_reason !== undefined && planned.unavailable_reason !== "") return planned.unavailable_reason;
  return scheme === "lineage" ? "Fewer than 10 lineages, too few to hold whole lineages out in 5 folds." : "Not run.";
};

export const learnerOf = (r: DatasetReport, scheme: Scheme, learner: Learner): LearnerResult | undefined =>
  splitOf(r, scheme)?.learners.find((l) => l.learner === learner);

export { EQUAL_SCORING_REVISION };

/** Whether the report scores a null fit's isolates at their fold's training
 *  prevalence, so every learner's rates cover every isolate. */
export const scoresNullFits = (r: DatasetReport): boolean => (r.evaluation?.revision ?? 0) >= EQUAL_SCORING_REVISION;

/** Isolates a learner row gave no call; null on a catalog row or a row
 *  written before null fits became no-calls. A revision-14 row (it carries
 *  `called_metrics`) calls every isolate, a null fit's at its prevalence. */
export const uncalledOf = (l: LearnerResult): number | null =>
  l.no_call === undefined ? null : l.called_metrics !== undefined ? 0 : l.no_call.resistant + l.no_call.susceptible;

/** Whether the report's learner rows count no-calls, so a null fit left its
 *  fold uncalled instead of scoring it (revisions 7 to 13). */
export const countsNoCalls = (r: DatasetReport): boolean => !scoresNullFits(r) && r.splits.some((s) => s.learners.some((l) => l.no_call !== undefined));

/** Folds where boosting kept no tree or the lasso kept no weight, read from
 *  fold_complexity. Reports without no-call counts scored those folds 0.5;
 *  revision-14 reports at their training prevalence. */
export const nullFolds = (l: LearnerResult): number => (l.learner === "gbt" || l.learner === "lasso" ? l.fold_complexity.filter((c) => c === 0).length : 0);

const callsEveryIsolate = (l: LearnerResult | undefined): l is CalledLearner => l !== undefined && l.metrics !== null && (uncalledOf(l) ?? 0) === 0;

/** The best learner, as report.rs picks it: the highest random-split
 *  balanced accuracy among learners whose rates cover every isolate on every
 *  split that ran (every learner from revision 14). The catalog is a baseline, not a candidate; a tie goes to the
 *  later learner, as Rust's max_by does. */
export const bestLearner = (r: DatasetReport): CalledLearner | undefined => {
  const lineage = splitOf(r, "lineage");
  return splitOf(r, "random")
    ?.learners.filter((l): l is CalledLearner => l.learner !== "catalog" && callsEveryIsolate(l) && (lineage === undefined || callsEveryIsolate(learnerOf(r, "lineage", l.learner))))
    .reduce<CalledLearner | undefined>((best, l) => (best === undefined || l.metrics.balanced_accuracy >= best.metrics.balanced_accuracy ? l : best), undefined);
};

/** Learners that score above the best learner on random splits but leave
 *  isolates without a call, with the uncalled count on each split. */
export const outscoringPartial = (r: DatasetReport): { learner: Learner; ba: number; uncalled: { scheme: Scheme; n: number }[] }[] => {
  const floor = bestLearner(r)?.metrics.balanced_accuracy ?? -Infinity;
  return (splitOf(r, "random")?.learners ?? []).flatMap((l) => {
    if (l.learner === "catalog" || l.metrics === null || l.metrics.balanced_accuracy <= floor) return [];
    const uncalled = (["random", "lineage"] as const).flatMap((scheme) => {
      const row = learnerOf(r, scheme, l.learner);
      const n = row === undefined ? 0 : (uncalledOf(row) ?? 0);
      return n > 0 ? [{ scheme, n }] : [];
    });
    return uncalled.length > 0 ? [{ learner: l.learner, ba: l.metrics.balanced_accuracy, uncalled }] : [];
  });
};

/** A catalog row's lineage-split reason in two or three words for a table
 *  cell; the report's own sentence stays in the inspector. Null when the row
 *  records no reason. */
export const splitReasonShort = (reason: string | undefined): string | null => {
  if (reason === undefined || reason === "") return null;
  if (/lineages as requested folds|fewer than \d+ lineages/i.test(reason)) return "Too few lineages";
  if (/lacks one or both classes/i.test(reason)) return "Fold lacks a class";
  return "Unavailable";
};

/** The interpretable model's rules as a table cell shows them: the first
 *  rule and how many more. Null when the model has none. */
export const ruleSummary = (rules: readonly string[]): { first: string; more: number } | null =>
  rules.length === 0 ? null : { first: rules[0]!, more: rules.length - 1 };

export const featureIndex = (r: DatasetReport): Map<number, FeatureCard> =>
  new Map(r.features.map((f) => [f.feature, f]));

/** The words a feature's catalog families give it, as `grm report` writes
 *  them into a rule term (`annotate::family_label`): the family when it
 *  matches one, or every family it matches, in code-unit order after their
 *  count, since its k-mers may lie in any of them. Null when none. */
export const familyLabel = (families: readonly string[]): string | null => {
  const names = [...new Set(families)].sort();
  if (names.length === 0) return null;
  return names.length === 1 ? names[0]! : `k-mers matching ${names.length} catalog families (${names.join(", ")})`;
};

/** The words a feature goes by: its catalog families, else the gene or
 *  product on the reference assembly, else its id. Views print where the
 *  feature sits beside this name (featurePlace); a rule term, which stands
 *  alone, takes termLabel. */
export const featureLabel = (card: FeatureCard | undefined, id: number): string =>
  (card === undefined ? null : familyLabel(card.catalog_families)) ?? card?.locus?.gene ?? card?.locus?.product ?? `k-mer feature ${id}`;

/** Where the feature sits: codon range inside a gene, or the distance to
 *  the nearest gene when it lies in none. The report holds protein-coding
 *  genes only, so the words say that and not "intergenic": the k-mers may
 *  still lie in an rRNA, tRNA or other RNA gene. Null without a placed
 *  locus. */
export const featurePlace = (card: FeatureCard | undefined): string | null => {
  const locus = card?.locus;
  if (locus === null || locus === undefined) return null;
  const gene = locus.gene ?? locus.product;
  if (locus.codons !== null) {
    const [a, b] = locus.codons;
    const span = a === b ? `codon ${int(a)}` : `codons ${int(a)}–${int(b)}`;
    return gene === null ? span : `${gene} ${span}`;
  }
  if (locus.nearest !== null) return `outside protein-coding genes, ${int(Math.abs(locus.nearest[1]))} bp from ${locus.nearest[0]}`;
  return null;
};

/** The words for k-mers that lie in no annotated protein-coding gene of
 *  the reference contig: the locus then names the nearest one, and the words
 *  give the distance, so a term does not read like one inside that gene.
 *  `grm` reads coding genes only, so the words name that kind and claim
 *  nothing of RNA genes. Null when the k-mers lie inside a gene, or no gene
 *  is named. */
const outsideWords = (card: FeatureCard | undefined): string | null => {
  const locus = card?.locus;
  if (locus === null || locus === undefined) return null;
  const name = locus.gene ?? locus.product;
  if (name === null || locus.codons !== null || locus.nearest === null) return null;
  return `k-mers ${int(Math.abs(locus.nearest[1]))} bp from the nearest annotated protein-coding gene, ${name}`;
};

/** The words a rule term's feature goes by, as `grm report` writes them
 *  into the index (`report::term_label`): its catalog families, else the
 *  gene or product its k-mers lie in, or their distance from the nearest
 *  one when they lie in none, else `unnamed`. */
export const termLabel = (card: FeatureCard | undefined, id: number, unnamed = `k-mer feature ${id}`): string =>
  (card === undefined ? null : familyLabel(card.catalog_families)) ?? outsideWords(card) ?? card?.locus?.gene ?? card?.locus?.product ?? unnamed;

export const ruleWords = (rule: Rule, features: Map<number, FeatureCard>): string =>
  `${rule.presence ? "Presence" : "Absence"} of ${termLabel(features.get(rule.feature), rule.feature)}`;

/** What a rule term's ties add to its words (`report::tie_words`), or null
 *  when nothing ties it. `nEquivalent` counts the other rules with the same
 *  covered and lost counts at the term's greedy step, so the term is one of
 *  that many and one more. A later term's step counts only the isolates the
 *  earlier terms left undecided. */
export const tieWords = (nEquivalent: number, position: number): string | null =>
  nEquivalent > 0
    ? `one of ${int(nEquivalent + 1)} features that separate ${position === 0 ? "these isolates" : "the isolates the earlier rules left undecided"} equally well`
    : null;

/** The family of a catalog gene group (`catalog::family_of_group`): the
 *  group without its " (X family)" suffix. */
export const familyOfGroup = (group: string): string => /^(.*) \(\S+ family\)$/.exec(group)?.[1] ?? group;

/** What GRM's known-gene rule for this drug does with a term's catalog
 *  families, as a predicate: "uses blaNDM but not ble", "does not use armA".
 *  Read from the report's own gene-group rows; null when they say nothing. A
 *  row exists for a gene group GRM's catalog reader called in at least one
 *  genome of the dataset. A named family with no row is left out: the report
 *  does not hold its class, so nothing is said of it. A row that does not
 *  say whether the rule uses it is read as no row. The groups are named,
 *  since one family can have a group the rule uses and one it does not. When
 *  a family has a group with a subfamily and one without, the second carries
 *  the family's bare name and is written "the group named X alone", so that
 *  one name does not stand on both sides. The names are sorted, so the words
 *  do not depend on the order of the report's rows. */
const ruleUseWords = (families: readonly string[], catalog: CatalogBaseline | null): string | null => {
  if (families.length === 0 || catalog === null) return null;
  const named = new Set(families);
  const rows = catalog.groups.filter((g) => typeof g.counts_for_rule === "boolean" && named.has(familyOfGroup(g.group)));
  // The families that have a group with a subfamily among these rows.
  const withSubfamily = new Set(rows.flatMap((g) => (familyOfGroup(g.group) === g.group ? [] : [familyOfGroup(g.group)])));
  const groups = (used: boolean) =>
    rows
      .flatMap((g) => (g.counts_for_rule === used ? [g.group] : []))
      .sort()
      .map((group) => (withSubfamily.has(group) ? `the group named ${group} alone` : group));
  const [uses, not] = [groups(true), groups(false)];
  if (uses.length === 0) return not.length === 0 ? null : `does not use ${listWords(not, "or")}`;
  return not.length === 0 ? `uses ${listWords(uses)}` : `uses ${listWords(uses)} but not ${listWords(not, "or")}`;
};

/** The clause the index gives a rule term that names catalog families
 *  (`report::family_rule_words`), or null when the report's gene-group rows
 *  say nothing of them. GRM's rule is the subject and the group its object,
 *  so the clause describes GRM's rule, which reads the catalog's class and
 *  subclass; it is no statement about what the gene does. */
export const familyRuleWords = (families: readonly string[], catalog: CatalogBaseline | null): string | null => {
  const use = ruleUseWords(families, catalog);
  return use === null ? null : `GRM's known-gene rule for this drug ${use}`;
};

/** The wording `indexRuleWords` writes: `report::RULE_TEXT_VERSION`. */
export const RULE_TEXT_VERSION = 3;

/** A rule term as the index's `scm_rules` words it (`report::rule_text`):
 *  presence or absence of the term's label, then which gene groups of its
 *  catalog families GRM's known-gene rule for this drug uses, then its ties.
 *  `position` is the term's place in the model, from 0. */
export const indexRuleWords = (r: DatasetReport, position: number): string => {
  const rule = r.models.scm.rules[position]!;
  const card = r.features.find((f) => f.feature === rule.feature);
  const clauses = [
    `${rule.presence ? "presence of" : "absence of"} ${termLabel(card, rule.feature, `feature ${rule.feature}`)}`,
    familyRuleWords(card?.catalog_families ?? [], r.catalog),
    tieWords(r.models.scm.steps[position]?.n_equivalent ?? 0, position),
  ];
  return clauses.filter((c) => c !== null).join("; ");
};

/** How many catalog families an index rule term names, when it names
 *  several (familyLabel's words); null for any other term. */
export const familiesNamed = (ruleText: string): number | null => {
  const n = /^(?:presence|absence) of k-mers matching (\d+) catalog families \(/.exec(ruleText)?.[1];
  return n === undefined ? null : Number(n);
};

/** Beside a term that names several catalog families. A feature matches a
 *  family when one of its sampled k-mers lies in any allele of it, so the
 *  k-mers may lie in any one of the families or span them
 *  (`annotate::family_label`): the term is one presence pattern, and the
 *  list does not say which gene an isolate carries. */
export const SEVERAL_FAMILIES_NOTE = "Not a list of genes an isolate carries: one k-mer pattern matches all of these families.";

/** Which wording the index's rule words use, for the Methods page. The
 *  index writes the words when a release is built (`report::rule_text`) and
 *  numbers the wording from version 3 on (`rule_text_version`). An older
 *  index has no number, and the note says so: two releases of one evaluation
 *  revision can differ in wording, so the revision alone does not tell them
 *  apart. A count of rules by the kind of gene they name is taken from these
 *  words, and it moves when the wording does. */
export const ruleWordingNote = (m: Pick<Manifest, "evaluation_revision" | "rule_text_version">): string => {
  const wording = m.rule_text_version === undefined ? "rule wording 2 or earlier: the index gives no number" : `rule wording ${m.rule_text_version}`;
  const named = [...(m.evaluation_revision === undefined ? [] : [`evaluation revision ${m.evaluation_revision}`]), wording];
  return `The rule words in the datasets table come from this release's index (${named.join(", ")}). A later release can word the same rule differently, so a count of rules by what they name holds for one wording only.`;
};

/** The sentences beside a rule term whose feature matches catalog families,
 *  where the view has the report: the note when it matches several, then
 *  which of their gene groups GRM's known-gene rule for this drug uses.
 *  Empty when neither applies. */
export const familyNotes = (families: readonly string[], catalog: CatalogBaseline | null): string[] => {
  const use = ruleUseWords(families, catalog);
  return [...(new Set(families).size > 1 ? [SEVERAL_FAMILIES_NOTE] : []), ...(use === null ? [] : [`Of the gene groups GRM found in these genomes, its known-gene rule for this drug ${use}.`])];
};

const TOKEN_KIND: Record<string, string> = { "C:": "class", "S:": "subclass" };

/** "S:CARBAPENEM" -> "a gene of subclass carbapenem". */
export const catalogTokenWords = (token: string): string => {
  const kind = TOKEN_KIND[token.slice(0, 2)] ?? "class";
  return `a gene of ${kind} ${token.slice(2).toLowerCase()}`;
};

/** The catalog rule (AND of groups, each an OR of tokens) in one sentence. */
export const catalogRuleWords = (rule: readonly (readonly string[])[]): string => {
  if (rule.length === 0) return "No acquired-gene rule applies to this drug.";
  const groups = rule.map((g) => g.map(catalogTokenWords).join(" or "));
  return `Resistant when the genome carries ${groups.join(", and also ")}.`;
};

/** A gene panel names catalog groups in `detected`, so its per-sample calls
 *  can be rebuilt; a probe assay (wild-type region probes) reports only its
 *  totals, under a label such as "mutation detected". */
export const isGenePanel = (r: DatasetReport, panel: PanelResult): boolean =>
  panel.detected.every(([group]) => r.catalog_groups.includes(group));

/** Sample-level panel result: a sample is positive when any group the panel
 *  detected anywhere in the dataset is among the sample's catalog groups.
 *  Null for a probe assay, whose per-sample calls are not in the report. */
export const panelCalls = (r: DatasetReport, panel: PanelResult): boolean[] | null => {
  if (!isGenePanel(r, panel)) return null;
  const detected = new Set(panel.detected.map(([group]) => group));
  const positive = new Set(r.catalog_groups.flatMap((g, i) => (detected.has(g) ? [i] : [])));
  return r.samples.map((s) => s.groups.some((g) => positive.has(g)));
};

export const truthOf = (r: DatasetReport): boolean[] => r.samples.map((s) => s.resistant);

/** Feature ids the final models lean on most, in this order: SCM rules,
 *  then CART splits, then the forest and boosting importances. */
export const modelFeatureOrder = (r: DatasetReport): number[] => {
  const seen = new Set<number>();
  r.models.scm.rules.forEach((rule) => seen.add(rule.feature));
  r.models.cart.nodes.forEach((n) => {
    if (n.split !== null) seen.add(n.split.feature);
  });
  r.models.rf_importance.forEach((w) => seen.add(w.feature));
  r.models.gbt_importance.forEach((w) => seen.add(w.feature));
  return [...seen];
};

export interface UnexplainedFeature {
  feature: number;
  card: FeatureCard;
  /** Missed resistant isolates carrying the feature. */
  inMissed: number;
  /** Susceptible isolates carrying the feature. */
  inSusceptible: number;
  /** Other ranked features with the same label (k-mers elsewhere in the
   *  same gene), merged into this row. */
  alike: number;
}

export interface Unexplained {
  panel: PanelResult;
  /** Resistant isolates the panel calls susceptible. */
  missed: number;
  susceptible: number;
  /** False for a probe assay: the report has no per-sample calls, so the
   *  missed isolates cannot be examined and `features` is empty. */
  perSample: boolean;
  features: UnexplainedFeature[];
}

/** What the panel misses, and which model features set the missed isolates
 *  apart: carried by a larger share of them than of susceptible isolates,
 *  ranked by that difference, so genes every isolate carries drop out. */
export const unexplained = (r: DatasetReport, panel: PanelResult, limit = 5): Unexplained => {
  const truth = truthOf(r);
  const calls = panelCalls(r, panel);
  const nSusceptible = truth.filter((t) => !t).length;
  if (calls === null) return { panel, missed: panel.metrics.fn, susceptible: nSusceptible, perSample: false, features: [] };
  const missed = truth.map((res, i) => res && !calls[i]);
  const nMissed = missed.filter(Boolean).length;
  const lift = (f: { inMissed: number; inSusceptible: number }) => f.inMissed / Math.max(nMissed, 1) - f.inSusceptible / Math.max(nSusceptible, 1);
  const cards = featureIndex(r);
  const features = modelFeatureOrder(r)
    .flatMap((f) => {
      const card = cards.get(f);
      if (card === undefined) return [];
      const present = decodeBits(card.bits, truth.length);
      return [
        {
          feature: f,
          card,
          inMissed: present.filter((p, i) => p && missed[i]).length,
          inSusceptible: present.filter((p, i) => p && !truth[i]).length,
          alike: 0,
        },
      ];
    })
    .filter((f) => f.inMissed > 0 && lift(f) > 0)
    .sort((a, b) => lift(b) - lift(a) || b.inMissed - a.inMissed)
    .reduce<UnexplainedFeature[]>((rows, f) => {
      const same = rows.find((x) => featureLabel(x.card, x.feature) === featureLabel(f.card, f.feature));
      if (same !== undefined) same.alike += 1;
      else rows.push(f);
      return rows;
    }, [])
    .slice(0, limit);
  return { panel, missed: nMissed, susceptible: nSusceptible, perSample: true, features };
};

/** Out-of-fold calls at the 0.5 threshold, and which samples they get
 *  wrong. An isolate without a call is not an error. */
export const oofErrors = (r: DatasetReport, learner: LearnerResult): boolean[] =>
  r.samples.map((s, i) => {
    const score = learner.scores[i];
    return score !== null && score !== undefined && (score >= 0.5) !== s.resistant;
  });

/** What predicts resistance here, and how well, in one sentence: the
 *  interpretable model's rule (named, or counted past two rules), the best
 *  learner on both splits, and the catalog when it does as well. */
export const answerSentence = (r: DatasetReport): string => {
  const best = bestLearner(r);
  const features = featureIndex(r);
  const scm = learnerOf(r, "random", "scm");
  const rules = r.models.scm.rules;
  const lower = (s: string) => s.charAt(0).toLowerCase() + s.slice(1);
  const join = r.models.scm.kind === "disjunction" ? " or " : " and ";
  const ruleText =
    rules.length === 0
      ? null
      : rules.length <= 2
        ? rules.map((rule, i) => (i === 0 ? ruleWords(rule, features) : lower(ruleWords(rule, features)))).join(join)
        : `A ${rules.length}-rule set covering machine led by ${lower(ruleWords(rules[0]!, features))}`;
  // The sentence names the first rule alone, so it says when other features tie it.
  const ties = rules.length === 0 ? null : tieWords(r.models.scm.steps[0]?.n_equivalent ?? 0, 0);
  const tieSentence = ties === null ? "" : ` ${rules.length === 1 ? "The rule" : "The first rule"} is ${ties}.`;
  const lineageBa = (l: Learner) => learnerOf(r, "lineage", l)?.metrics?.balanced_accuracy;
  const bothSplits = (l: CalledLearner) => {
    const lin = lineageBa(l.learner);
    return `${pct(l.metrics.balanced_accuracy)} balanced accuracy on random splits${lin === undefined ? "" : ` and ${pct(lin)} with lineages held out`}`;
  };
  if (best === undefined) {
    const hasLearners = splitOf(r, "random")?.learners.some((l) => l.learner !== "catalog") ?? false;
    const reason = hasLearners ? `${NO_BEST_LEARNER} on every split it ran` : "No learner result is recorded for this dataset";
    return ruleText === null ? `${reason}.` : `${ruleText} is the interpretable model; ${lower(reason)}.${tieSentence}`;
  }
  const catalogBa = learnerOf(r, "random", "catalog")?.metrics?.balanced_accuracy;
  const scmBa = scm?.metrics?.balanced_accuracy;
  const bar = lineageBa(best.learner) ?? best.metrics.balanced_accuracy;
  const catalogClause = catalogBa !== undefined && catalogBa >= bar - CATALOG_CUTOFF ? `; the known-gene catalog alone reaches ${pct(catalogBa)}` : "";
  const text =
    ruleText !== null && best.learner === "scm"
      ? `${ruleText} calls ${r.antibiotic} resistance at ${bothSplits(best)}, the best of the learners scored on every isolate`
      : ruleText !== null && scmBa !== undefined
        ? `${ruleText} calls ${r.antibiotic} resistance at ${pct(scmBa)} balanced accuracy; ${BEST_LEARNER_WORDS}, ${learnerInSentence(best.learner)}, reaches ${bothSplits(best)}`
        : `${capitalize(BEST_LEARNER_WORDS)}, ${learnerInSentence(best.learner)}, reaches ${bothSplits(best)}`;
  return `${text}${catalogClause}.${tieSentence}`;
};

export interface SeqRun {
  text: string;
  /** Part of the k-mer span (uppercase in the source). */
  span: boolean;
  /** Bases left out of a long span, shown as a count. */
  elided?: number;
}

/** Context sequence as runs, the span marked. A span longer than `max`
 *  bases keeps its first and last `keep` bases and states how many it left
 *  out, so a whole accessory gene does not fill the card. */
export const contextRuns = (context: string, max = 120, keep = 48): SeqRun[] =>
  (context.match(/[A-Z]+|[^A-Z]+/g) ?? []).flatMap((run): SeqRun[] => {
    const span = /^[A-Z]/.test(run);
    if (!span || run.length <= max) return [{ text: run, span }];
    return [
      { text: run.slice(0, keep), span: true },
      { text: "", span: true, elided: run.length - 2 * keep },
      { text: run.slice(-keep), span: true },
    ];
  });

/** featurePlace without the label it would repeat: "blaKPC codons 1–290"
 *  under a heading that already says blaKPC reads "Codons 1–290". */
export const featureWhere = (card: FeatureCard | undefined, id: number): string | null => {
  const place = featurePlace(card);
  const label = featureLabel(card, id);
  if (place === null || !place.startsWith(`${label} `)) return place;
  const rest = place.slice(label.length + 1);
  return rest.charAt(0).toUpperCase() + rest.slice(1);
};
