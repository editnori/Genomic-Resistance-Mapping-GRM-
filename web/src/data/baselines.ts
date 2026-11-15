// No-genome baselines beside a dataset report's scores, from
// baselines/<id>.json beside it (`grm baselines`, src/baselines.rs). Two
// predictors that never read a genome, scored on the report's own folds:
// for an MIC model, constants at each training fold's median and mode MIC,
// scored on the model's exact MICs; for each split, lineage identity, which
// scores a test genome by the resistant share of its lineage among the
// fold's training genomes, or the fold's training prevalence when the
// lineage is unseen. The file is optional, and a file that does not
// describe the report as fetched shows nothing: the dataset report checks
// the report's SHA-256, and the overview and the datasets table, which read
// the index alone, check the dataset, the revision and the model's own
// figure against the index row.

import { median } from "./format";
import type { CatalogRow, DatasetReport, Scheme } from "./schema";
import { ShapeError } from "./validate";

export const BASELINES_SCHEMA = 1;

/** `agree` of `n` exact MICs within one doubling dilution, and the share. */
export interface Agreement {
  agree: number;
  n: number;
  essential_agreement: number;
}

/** One MIC fold: its exact training and test MICs and the two constants. */
export interface FoldConstant {
  train_exact: number;
  test_exact: number;
  /** Log2 dilution; null when the fold's training isolates hold no exact MIC. */
  median: number | null;
  upper_median_differs: boolean;
  mode: number | null;
  mode_tied: boolean;
}

/** Isolates with bit-identical held-out predictions lie in one rebuilt
 *  fold: `groups` groups of two or more, holding `isolates` isolates, leaving
 *  out folds that predict one value for every isolate (`constant_folds`, as
 *  when boosting keeps no tree). The folds are confirmed when a group exists. */
export interface PredictionCheck {
  confirmed: boolean;
  groups: number;
  isolates: number;
  constant_folds: number[];
}

/** The stored training-fold thresholds the rebuilt folds reproduce. Written
 *  beside the prediction check; on about half the datasets a wrong draw
 *  matches them too, so the site reads the prediction check alone. */
export interface ThresholdCheck {
  confirmed: boolean;
  checked: number;
}

export interface MicBaseline {
  folds: number;
  seed: number;
  /** A dataset counts as confirmed by this check. */
  prediction_check: PredictionCheck;
  threshold_check: ThresholdCheck;
  /** Exact-MIC test isolates with a k-mer-identical copy in their training
   *  fold: the MIC folds are random and can split such a pair. */
  exact_with_training_copy: number;
  n: number;
  n_exact: number;
  /** The model's count, recounted from the report's predictions. */
  model: Agreement;
  /** Null when a fold has no exact training MIC. */
  median: Agreement | null;
  mode: Agreement | null;
  fold_constants: FoldConstant[];
}

export interface LineageBaseline {
  scheme: Scheme;
  folds: number;
  n: number;
  tp: number;
  fp: number;
  tn: number;
  fn: number;
  balanced_accuracy: number;
  sensitivity: number;
  specificity: number;
  /** Test genomes whose lineage no training genome of their fold shares. */
  unseen: number;
  unseen_share: number;
}

export interface Baselines {
  schema: number;
  dataset: string;
  report_sha256: string;
  evaluation_revision: number;
  mic: MicBaseline | null;
  lineage: LineageBaseline[];
}

type Obj = Record<string, unknown>;
const isObj = (v: unknown): v is Obj => typeof v === "object" && v !== null && !Array.isArray(v);
const isCount = (v: unknown): v is number => typeof v === "number" && Number.isInteger(v) && v >= 0;
const isDilution = (v: unknown): v is number => typeof v === "number" && Number.isInteger(v);

function need(cond: boolean, message: string): asserts cond {
  if (!cond) throw new ShapeError(message);
}

/** The share as the command divides it, so a file's share must be its own
 *  counts' to the last bit. */
const share = (x: number, n: number): number => x / Math.max(n, 1);
const ratio = (x: number, n: number): number => (n === 0 ? 0 : x / n);

const checkAgreement = (v: unknown, where: string, nExact: number): void => {
  need(isObj(v) && isCount(v.agree) && isCount(v.n) && v.agree <= v.n && v.n === nExact, `${where} should be {agree, n, essential_agreement} over the ${nExact} exact MICs`);
  need(v.essential_agreement === share(v.agree, v.n), `${where}: essential_agreement should be agree / n`);
};

const checkMic = (m: unknown): void => {
  need(isObj(m), "baselines mic should be null or an object");
  need(isCount(m.folds) && m.folds >= 2 && isCount(m.seed), "baselines mic: expected folds and seed");
  need(isCount(m.n) && isCount(m.n_exact) && m.n_exact <= m.n, "baselines mic: expected n and n_exact");
  const folds = m.folds;
  const p = m.prediction_check;
  need(
    isObj(p) && typeof p.confirmed === "boolean" && isCount(p.groups) && isCount(p.isolates) && Array.isArray(p.constant_folds),
    "baselines mic.prediction_check should be {confirmed, groups, isolates, constant_folds}",
  );
  need(p.confirmed === p.groups > 0 && p.isolates >= 2 * p.groups && p.isolates <= m.n, "baselines mic.prediction_check: confirmed should mean a group of two or more isolates, all of them isolates with an MIC");
  need(
    (p.constant_folds as unknown[]).every((f, i, a) => isCount(f) && f < folds && (i === 0 || f > (a[i - 1] as number))),
    "baselines mic.prediction_check: constant_folds should be ascending fold numbers",
  );
  const t = m.threshold_check;
  need(isObj(t) && typeof t.confirmed === "boolean" && isCount(t.checked) && t.confirmed === t.checked > 0 && t.checked <= folds, "baselines mic.threshold_check should be {confirmed, checked}, confirmed when a stored threshold was checked");
  need(isCount(m.exact_with_training_copy) && m.exact_with_training_copy <= m.n_exact, "baselines mic: exact_with_training_copy should count exact-MIC isolates");
  const nExact = m.n_exact;
  checkAgreement(m.model, "baselines mic.model", nExact);
  for (const key of ["median", "mode"] as const) if (m[key] !== null) checkAgreement(m[key], `baselines mic.${key}`, nExact);
  need(Array.isArray(m.fold_constants) && m.fold_constants.length === m.folds, "baselines mic: one fold_constants entry per fold");
  let test = 0;
  (m.fold_constants as unknown[]).forEach((f, i) => {
    const where = `baselines mic.fold_constants[${i}]`;
    need(isObj(f) && isCount(f.train_exact) && isCount(f.test_exact), `${where} should count its exact training and test MICs`);
    need((f.median === null || isDilution(f.median)) && (f.mode === null || isDilution(f.mode)), `${where}: median and mode should be whole dilutions or null`);
    need(typeof f.upper_median_differs === "boolean" && typeof f.mode_tied === "boolean", `${where}: expected upper_median_differs and mode_tied`);
    need((f.median === null) === (f.train_exact === 0) && (f.mode === null) === (f.train_exact === 0), `${where}: a constant is null exactly when the fold has no exact training MIC`);
    test += f.test_exact;
  });
  need(test === nExact, "baselines mic: the folds' exact test MICs should add up to n_exact");
};

const checkLineage = (s: unknown, i: number): void => {
  const where = `baselines lineage[${i}]`;
  need(isObj(s) && (s.scheme === "random" || s.scheme === "lineage") && isCount(s.folds), `${where} should be {scheme, folds, ...} for the random or lineage split`);
  need(isCount(s.n) && isCount(s.tp) && isCount(s.fp) && isCount(s.tn) && isCount(s.fn) && s.tp + s.fp + s.tn + s.fn === s.n, `${where}: tp, fp, tn and fn should be whole numbers that sum to n`);
  const sens = ratio(s.tp, s.tp + s.fn);
  const spec = ratio(s.tn, s.tn + s.fp);
  need(s.sensitivity === sens && s.specificity === spec && s.balanced_accuracy === (sens + spec) / 2, `${where}: the rates should be its counts' own`);
  need(isCount(s.unseen) && s.unseen <= s.n && s.unseen_share === share(s.unseen, s.n), `${where}: unseen should count test genomes and unseen_share be unseen / n`);
};

/** baselines/<id>.json: the schema this view reads, and figures that are
 *  their own counts' to the last bit. */
export function parseBaselines(raw: unknown): Baselines {
  need(isObj(raw), "baselines: expected a JSON object");
  if (raw.schema !== BASELINES_SCHEMA) throw new ShapeError(`baselines file is schema ${String(raw.schema)}; this dashboard reads schema ${BASELINES_SCHEMA}`, raw.schema);
  need(typeof raw.dataset === "string" && typeof raw.report_sha256 === "string" && /^[0-9a-f]{64}$/.test(raw.report_sha256), "baselines: expected dataset and a SHA-256 report_sha256");
  need(isCount(raw.evaluation_revision), "baselines: evaluation_revision should be a whole number");
  if (raw.mic !== null) checkMic(raw.mic);
  need(Array.isArray(raw.lineage), "baselines: lineage should be a list of splits");
  (raw.lineage as unknown[]).forEach(checkLineage);
  need(new Set((raw.lineage as Obj[]).map((s) => s.scheme)).size === raw.lineage.length, "baselines: one lineage entry per split");
  return raw as unknown as Baselines;
}

/** The file as it describes `report`: the same dataset and evaluation
 *  revision, and `digest`, the SHA-256 of the report's bytes as fetched,
 *  equal to the file's pin. A file pinned to other bytes, or a report whose
 *  digest is not known (undefined while it loads, null where the browser
 *  computes none), gives null, and the report shows what it shows without
 *  the file. The MIC part shows only where the report's own MIC model has
 *  its exact-MIC count and figure, and a split's entry only where the
 *  report has that split with the same folds and genomes. */
export const baselinesFor = (file: Baselines | null, report: DatasetReport, digest: string | null | undefined): Baselines | null => {
  if (file === null || digest === null || digest === undefined) return null;
  if (file.dataset !== report.id || file.evaluation_revision !== report.evaluation?.revision || file.report_sha256 !== digest) return null;
  const mic = file.mic !== null && report.mic !== null && file.mic.n_exact === report.mic.n_exact && file.mic.model.essential_agreement === report.mic.essential_agreement ? file.mic : null;
  const lineage = file.lineage.filter((s) => s.n === report.samples.length && report.splits.some((x) => x.scheme === s.scheme && x.folds === s.folds));
  return { ...file, mic, lineage };
};

/** A score's lead over a baseline, as a share: positive when the score is
 *  higher. */
export const leadOver = (score: number, baseline: number): number => score - baseline;

/** A lead in signed points: 0.183 -> "+18.3", -0.021 -> "−2.1". */
export const signedPoints = (lead: number, digits = 1): string => {
  const v = (Math.abs(lead) * 100).toFixed(digits);
  if (Number(v) === 0) return (0).toFixed(digits);
  return `${lead > 0 ? "+" : "−"}${v}`;
};

/** A lead in words beside the baseline it is measured from: "12.6 points
 *  above it", "3.1 points below it", "level with it". */
export const leadWords = (lead: number): string => (lead === 0 ? "level with it" : `${(Math.abs(lead) * 100).toFixed(1)} points ${lead > 0 ? "above" : "below"} it`);

/** The lineage baseline of one split, or undefined. */
export const lineageOf = (b: Baselines | null, scheme: Scheme): LineageBaseline | undefined => b?.lineage.find((s) => s.scheme === scheme);

// ---- the overview and the datasets table ------------------------------------------

/** One dataset's MIC constant, as the index row and the file agree. */
export interface MicConstant {
  /** The model's agreement on exact MICs, the index row's figure. */
  model: number;
  /** A constant at each training fold's lower median MIC. */
  constant: number;
  /** A constant at each training fold's most frequent MIC. */
  mode: number | null;
  /** The held-out predictions confirm the rebuilt folds. */
  confirmed: boolean;
}

/** The constant for an index row, when `file` names that dataset and the
 *  index's revision, and its model figure is the row's own. The index has
 *  no report hash, so this is the check the index allows; the dataset
 *  report checks the hash. */
export const micConstantFor = (row: CatalogRow, file: Baselines | undefined, revision: number | undefined): MicConstant | null => {
  if (file === undefined || file.dataset !== row.id || file.evaluation_revision !== revision) return null;
  const m = file.mic;
  if (m === null || m.median === null || row.mic_essential_agreement === null || m.model.essential_agreement !== row.mic_essential_agreement) return null;
  return { model: m.model.essential_agreement, constant: m.median.essential_agreement, mode: m.mode?.essential_agreement ?? null, confirmed: m.prediction_check.confirmed };
};

/** The constants of every index row a file covers, by dataset id. */
export const micConstants = (rows: readonly CatalogRow[], files: ReadonlyMap<string, Baselines>, revision: number | undefined): ReadonlyMap<string, MicConstant> => {
  const out = new Map<string, MicConstant>();
  for (const row of rows) {
    const c = micConstantFor(row, files.get(row.id), revision);
    if (c !== null) out.set(row.id, c);
  }
  return out;
};

export interface MicBaselineSummary {
  /** Datasets with an MIC model in the index. */
  n: number;
  /** Of those, the ones with a constant. */
  covered: number;
  /** Medians over the covered datasets. */
  medianModel: number | null;
  medianConstant: number | null;
  /** The median over the covered datasets of the model minus its constant. */
  medianLead: number | null;
  /** Covered datasets where the model's figure is at or below the constant's. */
  notAbove: number;
  /** Covered datasets where the constant reaches 0.90, where the model
   *  does, and where both do. */
  constantAt90: number;
  modelAt90: number;
  bothAt90: number;
  /** Covered datasets whose MIC folds the held-out predictions do not confirm. */
  unconfirmed: number;
}

export const micBaselineSummary = (rows: readonly CatalogRow[], constants: ReadonlyMap<string, MicConstant>): MicBaselineSummary => {
  const withModel = rows.filter((r) => r.mic_essential_agreement !== null);
  const covered = withModel.flatMap((r) => {
    const c = constants.get(r.id);
    return c === undefined ? [] : [c];
  });
  return {
    n: withModel.length,
    covered: covered.length,
    medianModel: median(covered.map((c) => c.model)),
    medianConstant: median(covered.map((c) => c.constant)),
    medianLead: median(covered.map((c) => leadOver(c.model, c.constant))),
    notAbove: covered.filter((c) => leadOver(c.model, c.constant) <= 0).length,
    constantAt90: covered.filter((c) => c.constant >= 0.9).length,
    modelAt90: covered.filter((c) => c.model >= 0.9).length,
    bothAt90: covered.filter((c) => c.model >= 0.9 && c.constant >= 0.9).length,
    unconfirmed: covered.filter((c) => !c.confirmed).length,
  };
};
