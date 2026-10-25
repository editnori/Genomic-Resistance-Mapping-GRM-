// Per-dataset derivations for the inspector. Pure functions of one
// DatasetReport.

import type {
  BinaryMetrics,
  DatasetReport,
  FeatureCard,
  Learner,
  LearnerResult,
  PanelResult,
  Rule,
  Scheme,
  SplitResult,
} from "../schema";
import { decodeBits } from "../bits";
import { capitalize, int, learnerInSentence, pct } from "../format";
import { CATALOG_CUTOFF } from "./overview";

/** The best learner in words, for every view that names it. report.rs and
 *  bestLearner rank by random-split balanced accuracy, and only learners
 *  that call every isolate on every split that ran are candidates. */
export const BEST_LEARNER_WORDS = "the best random-split learner among those that call every isolate";

/** A catalog row with no best learner: every learner leaves some isolate
 *  without a call. */
export const NO_BEST_LEARNER = "No learner calls every isolate";

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

/** Isolates a learner row gave no call; null on a catalog row or a row
 *  written before null fits became no-calls. */
export const uncalledOf = (l: LearnerResult): number | null => (l.no_call === undefined ? null : l.no_call.resistant + l.no_call.susceptible);

/** Whether the report's learner rows count no-calls, so a null fit left its
 *  fold uncalled instead of scoring it 0.5. */
export const countsNoCalls = (r: DatasetReport): boolean => r.splits.some((s) => s.learners.some((l) => l.no_call !== undefined));

/** Folds where boosting kept no tree or the lasso kept no weight, read from
 *  fold_complexity. Reports without no-call counts scored those folds 0.5. */
export const nullFolds = (l: LearnerResult): number => (l.learner === "gbt" || l.learner === "lasso" ? l.fold_complexity.filter((c) => c === 0).length : 0);

const callsEveryIsolate = (l: LearnerResult | undefined): l is CalledLearner => l !== undefined && l.metrics !== null && (uncalledOf(l) ?? 0) === 0;

/** The best learner, as report.rs picks it: the highest random-split
 *  balanced accuracy among learners that call every isolate on every split
 *  that ran. The catalog is a baseline, not a candidate; a tie goes to the
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

/** The words a feature goes by: its catalog family, else the gene or product
 *  on the reference assembly, else its id. */
export const featureLabel = (card: FeatureCard | undefined, id: number): string =>
  card?.catalog_families[0] ?? card?.locus?.gene ?? card?.locus?.product ?? `k-mer feature ${id}`;

/** Where the feature sits: codon range inside a gene, or the distance to
 *  the nearest gene when intergenic. Null without a placed locus. */
export const featurePlace = (card: FeatureCard | undefined): string | null => {
  const locus = card?.locus;
  if (locus === null || locus === undefined) return null;
  const gene = locus.gene ?? locus.product;
  if (locus.codons !== null) {
    const [a, b] = locus.codons;
    const span = a === b ? `codon ${int(a)}` : `codons ${int(a)}–${int(b)}`;
    return gene === null ? span : `${gene} ${span}`;
  }
  if (locus.nearest !== null) return `intergenic, ${int(Math.abs(locus.nearest[1]))} bp from ${locus.nearest[0]}`;
  return null;
};

export const ruleWords = (rule: Rule, features: Map<number, FeatureCard>): string =>
  `${rule.presence ? "Presence" : "Absence"} of ${featureLabel(features.get(rule.feature), rule.feature)}`;

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
  const lineageBa = (l: Learner) => learnerOf(r, "lineage", l)?.metrics?.balanced_accuracy;
  const bothSplits = (l: CalledLearner) => {
    const lin = lineageBa(l.learner);
    return `${pct(l.metrics.balanced_accuracy)} balanced accuracy on random splits${lin === undefined ? "" : ` and ${pct(lin)} with lineages held out`}`;
  };
  if (best === undefined) {
    const hasLearners = splitOf(r, "random")?.learners.some((l) => l.learner !== "catalog") ?? false;
    const reason = hasLearners ? `${NO_BEST_LEARNER} on every split it ran` : "No learner result is recorded for this dataset";
    return ruleText === null ? `${reason}.` : `${ruleText} is the interpretable model; ${lower(reason)}.`;
  }
  const catalogBa = learnerOf(r, "random", "catalog")?.metrics?.balanced_accuracy;
  const scmBa = scm?.metrics?.balanced_accuracy;
  const bar = lineageBa(best.learner) ?? best.metrics.balanced_accuracy;
  const catalogClause = catalogBa !== undefined && catalogBa >= bar - CATALOG_CUTOFF ? `; the known-gene catalog alone reaches ${pct(catalogBa)}` : "";
  const text =
    ruleText !== null && best.learner === "scm"
      ? `${ruleText} calls ${r.antibiotic} resistance at ${bothSplits(best)}, the best of the learners that call every isolate`
      : ruleText !== null && scmBa !== undefined
        ? `${ruleText} calls ${r.antibiotic} resistance at ${pct(scmBa)} balanced accuracy; ${BEST_LEARNER_WORDS}, ${learnerInSentence(best.learner)}, reaches ${bothSplits(best)}`
        : `${capitalize(BEST_LEARNER_WORDS)}, ${learnerInSentence(best.learner)}, reaches ${bothSplits(best)}`;
  return `${text}${catalogClause}.`;
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
