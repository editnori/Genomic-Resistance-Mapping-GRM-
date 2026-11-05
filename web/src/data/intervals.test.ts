import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseDataset, parseManifest } from "./validate";
import type { DatasetReport } from "./schema";
import { intervalKey, intervalsFor, parseIntervals, type Intervals } from "./intervals";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const report = parseDataset(fixtureFile(`/data/datasets/${manifest.datasets[0]!.id}.json`));

/** The file `grm intervals` writes for `r`, with placeholder limits around
 *  each rate (the arithmetic is src/intervals.rs's, tested there). */
const intervalsOf = (r: DatasetReport): Intervals => {
  const rate = (x: number, n: number) => (n === 0 ? null : { x, n, lo: x / n / 2, hi: Math.min(1, x / n + 0.01) });
  return {
    schema: 1,
    method: "clopper-pearson",
    confidence: 0.95,
    dataset: r.id,
    report_sha256: "0".repeat(64),
    evaluation_revision: r.evaluation?.revision ?? null,
    splits: r.splits.map((s) => ({
      scheme: s.scheme,
      learners: s.learners.flatMap((l) =>
        l.metrics === null ? [] : [{ learner: l.learner, counts: { n: l.metrics.n, tp: l.metrics.tp, fp: l.metrics.fp, tn: l.metrics.tn, fn: l.metrics.fn }, vme: rate(l.metrics.fn, l.metrics.tp + l.metrics.fn), me: rate(l.metrics.fp, l.metrics.tn + l.metrics.fp) }],
      ),
    })),
  };
};

const copy = <T>(v: T): T => JSON.parse(JSON.stringify(v)) as T;

describe("offline intervals file", () => {
  it("parses the file grm intervals writes and keys every row of the report it describes", () => {
    const file = parseIntervals(copy(intervalsOf(report)));
    const rows = intervalsFor(file, report);
    const expected = report.splits.flatMap((s) => s.learners.filter((l) => l.metrics !== null).map((l) => intervalKey(s.scheme, l.learner)));
    expect([...rows.keys()].sort()).toEqual(expected.sort());
    const first = report.splits[0]!.learners.find((l) => l.metrics !== null)!;
    expect(rows.get(intervalKey(report.splits[0]!.scheme, first.learner))!.vme!.n).toBe(first.metrics!.tp + first.metrics!.fn);
  });

  it("refuses another schema, method or level, counts that do not sum, and a rate that is not its counts'", () => {
    const good = intervalsOf(report);
    expect(() => parseIntervals({ ...copy(good), schema: 2 })).toThrow(/schema 2/);
    expect(() => parseIntervals({ ...copy(good), method: "wilson" })).toThrow(/Clopper-Pearson/);
    expect(() => parseIntervals({ ...copy(good), confidence: 0.9 })).toThrow(/95%/);
    const bad = copy(good);
    bad.splits[0]!.learners[0]!.counts.n += 1;
    expect(() => parseIntervals(bad)).toThrow(/sum to n/);
    const swapped = copy(good);
    const row = swapped.splits[0]!.learners[0]!;
    [row.vme, row.me] = [row.me, row.vme];
    expect(() => parseIntervals(swapped)).toThrow(/count/);
    const outside = copy(good);
    outside.splits[0]!.learners[0]!.vme!.hi = 1.5;
    expect(() => parseIntervals(outside)).toThrow(/hi ≤ 1/);
  });

  it("gives no interval for a row whose counts, dataset or revision differ from the report's", () => {
    const changed = copy(report);
    const scheme = changed.splits[0]!.scheme;
    const row = changed.splits[0]!.learners.find((l) => l.metrics !== null)!;
    row.metrics!.tp -= 1;
    row.metrics!.fn += 1;
    const rows = intervalsFor(intervalsOf(report), changed);
    expect(rows.has(intervalKey(scheme, row.learner))).toBe(false);
    expect(rows.size).toBe(intervalsFor(intervalsOf(report), report).size - 1);

    expect(intervalsFor({ ...intervalsOf(report), dataset: "other__drug" }, report).size).toBe(0);
    expect(intervalsFor({ ...intervalsOf(report), evaluation_revision: 99 }, report).size).toBe(0);
    expect(intervalsFor(null, report).size).toBe(0);
  });
});
