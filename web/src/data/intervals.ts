// Offline confidence intervals for a dataset report's error rates, from
// intervals/<id>.json beside it (`grm intervals`, src/intervals.rs): a 95%
// Clopper-Pearson interval for each learner row's VME and ME, computed from
// the report's own counts. The file is optional. An interval shows only for
// a row whose counts in the file equal the report's, so a file left from an
// earlier report shows nothing rather than another report's numbers.

import type { DatasetReport, Learner, Scheme } from "./schema";
import { ShapeError } from "./validate";

export const INTERVALS_SCHEMA = 1;

/** x events of n isolates, and the interval around x / n. */
export interface RateInterval {
  x: number;
  n: number;
  lo: number;
  hi: number;
}

export interface IntervalCounts {
  n: number;
  tp: number;
  fp: number;
  tn: number;
  fn: number;
}

export interface IntervalRow {
  learner: Learner;
  counts: IntervalCounts;
  /** Null when no resistant isolate was called (no denominator). */
  vme: RateInterval | null;
  /** Null when no susceptible isolate was called. */
  me: RateInterval | null;
}

export interface Intervals {
  schema: number;
  method: "clopper-pearson";
  confidence: number;
  dataset: string;
  report_sha256: string;
  evaluation_revision: number | null;
  splits: { scheme: Scheme; learners: IntervalRow[] }[];
}

type Obj = Record<string, unknown>;
const isObj = (v: unknown): v is Obj => typeof v === "object" && v !== null && !Array.isArray(v);
const isCount = (v: unknown): v is number => typeof v === "number" && Number.isInteger(v) && v >= 0;

function need(cond: boolean, message: string): asserts cond {
  if (!cond) throw new ShapeError(message);
}

const checkRate = (v: unknown, where: string): void => {
  if (v === null) return;
  need(
    isObj(v) && isCount(v.x) && isCount(v.n) && v.n > 0 && v.x <= v.n && typeof v.lo === "number" && typeof v.hi === "number" && v.lo >= 0 && v.lo <= v.x / v.n && v.x / v.n <= v.hi && v.hi <= 1,
    `${where} should be null or {x, n, lo, hi} with 0 ≤ lo ≤ x/n ≤ hi ≤ 1`,
  );
};

/** intervals/<id>.json: the schema, the method and level this view names,
 *  and each row's counts and intervals, which must be the counts' own. */
export function parseIntervals(raw: unknown): Intervals {
  need(isObj(raw), "intervals: expected a JSON object");
  if (raw.schema !== INTERVALS_SCHEMA) throw new ShapeError(`intervals file is schema ${String(raw.schema)}; this dashboard reads schema ${INTERVALS_SCHEMA}`, raw.schema);
  need(raw.method === "clopper-pearson" && raw.confidence === 0.95, "intervals: expected 95% Clopper-Pearson intervals");
  need(typeof raw.dataset === "string" && typeof raw.report_sha256 === "string" && Array.isArray(raw.splits), "intervals: expected dataset, report_sha256 and splits");
  need(raw.evaluation_revision === null || typeof raw.evaluation_revision === "number", "intervals: evaluation_revision should be a number or null");
  (raw.splits as unknown[]).forEach((s, i) => {
    need(isObj(s) && typeof s.scheme === "string" && Array.isArray(s.learners), `intervals splits[${i}] should be {scheme, learners}`);
    (s.learners as unknown[]).forEach((l, k) => {
      const where = `intervals splits[${i}].learners[${k}]`;
      need(isObj(l) && typeof l.learner === "string" && isObj(l.counts), `${where} should be {learner, counts, vme, me}`);
      const c = l.counts;
      need(isCount(c.n) && isCount(c.tp) && isCount(c.fp) && isCount(c.tn) && isCount(c.fn) && c.tp + c.fp + c.tn + c.fn === c.n, `${where}: counts should be whole numbers that sum to n`);
      checkRate(l.vme, `${where} vme`);
      checkRate(l.me, `${where} me`);
      // Each rate must be the one its counts give.
      need(l.vme === null ? c.tp + c.fn === 0 : (l.vme as Obj).x === c.fn && (l.vme as Obj).n === c.tp + c.fn, `${where}: vme should count fn over tp + fn`);
      need(l.me === null ? c.tn + c.fp === 0 : (l.me as Obj).x === c.fp && (l.me as Obj).n === c.tn + c.fp, `${where}: me should count fp over tn + fp`);
    });
  });
  return raw as unknown as Intervals;
}

/** The intervals of one learner row, keyed `${scheme}/${learner}`. */
export type RowIntervals = ReadonlyMap<string, { vme: RateInterval | null; me: RateInterval | null }>;

export const intervalKey = (scheme: Scheme, learner: Learner): string => `${scheme}/${learner}`;

/** The rows of `file` that describe `report` as it is: the same dataset, the
 *  same evaluation revision, and for each row the same confusion counts. Any
 *  other row is left out, so the view shows that row as it did before. */
export const intervalsFor = (file: Intervals | null, report: DatasetReport): RowIntervals => {
  const out = new Map<string, { vme: RateInterval | null; me: RateInterval | null }>();
  if (file === null || file.dataset !== report.id || file.evaluation_revision !== (report.evaluation?.revision ?? null)) return out;
  for (const split of file.splits) {
    const saved = report.splits.find((s) => s.scheme === split.scheme);
    for (const row of split.learners) {
      const m = saved?.learners.find((l) => l.learner === row.learner)?.metrics;
      if (m === null || m === undefined) continue;
      const c = row.counts;
      if (m.n === c.n && m.tp === c.tp && m.fp === c.fp && m.tn === c.tn && m.fn === c.fn) out.set(intervalKey(split.scheme, row.learner), { vme: row.vme, me: row.me });
    }
  }
  return out;
};
