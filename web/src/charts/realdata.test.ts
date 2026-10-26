// The published report's species files, when GRM_REPORT_DIR names a report
// directory (the variable vite.config.ts reads): lineage blocks, profile
// totals and pair coverage checked against the files themselves, plus the
// revision-7 M. tuberculosis acceptance figures. Skipped without it, so the
// fixture-only suite stays self-contained.

import { describe, expect, it } from "vitest";
import { parseManifest, parseSpecies } from "../data/validate";
import { blockCalls, blockSize, drugTotals, lineageBlocks, testedOf, treeBreak } from "./lineageBlocks";
import { MIN_JOINT, crossPairs, jointCalls, pairSummary } from "./crossPairs";

const env = (globalThis as { process?: { env: Record<string, string | undefined> } }).process?.env ?? {};
const dir = env.GRM_REPORT_DIR;

const read = async (path: string): Promise<unknown> => {
  // A computed specifier keeps Node's types out of the app's compile.
  const fs = (await import(/* @vite-ignore */ ["node", "fs/promises"].join(":"))) as { readFile(p: string, encoding: string): Promise<string> };
  return JSON.parse(await fs.readFile(`${dir}/${path}`, "utf8")) as unknown;
};

describe.skipIf(dir === undefined)("published species files", () => {
  it("tile every species into blocks that keep every call, and enumerate every pair", async () => {
    const manifest = parseManifest(await read("index.json"));
    const lines: string[] = [];
    for (const row of manifest.species.filter((s) => s.genomes_analyzed > 0)) {
      const s = parseSpecies(await read(`species/${row.slug}.json`));
      const blocks = lineageBlocks(s);
      expect(blocks.reduce((a, b) => a + blockSize(b), 0)).toBe(s.genomes.length);
      blocks.forEach((b, i) => expect(b.start).toBe(i === 0 ? 0 : blocks[i - 1]!.end));
      expect(blocks.length).toBeGreaterThanOrEqual(7);
      expect(blocks.length).toBeLessThanOrEqual(60);
      const calls = blockCalls(s, blocks);
      s.drugs.forEach((_, j) => {
        const want = drugTotals(s, j);
        const got = calls.reduce((a, c) => ({ r: a.r + c[j]!.r, i: a.i + c[j]!.i, s: a.s + c[j]!.s, untested: a.untested + c[j]!.untested }), { r: 0, i: 0, s: 0, untested: 0 });
        expect(got).toEqual(want);
      });
      const pairs = crossPairs(s);
      const summary = pairSummary(pairs);
      expect(summary.total).toBe((s.drugs.length * (s.drugs.length - 1)) / 2);
      expect(summary.evaluated).toBe(s.cross.length);
      pairs.forEach((p) => {
        // The file's evaluated counts are the profile tallies report.rs made.
        expect(jointCalls(s, p.ia, p.ib)).toEqual({ n: p.n, rr: p.rr, rs: p.rs, sr: p.sr, ss: p.ss });
        expect(p.status === "not-evaluated" ? p.n < MIN_JOINT : p.n >= MIN_JOINT).toBe(true);
      });
      const brk = treeBreak(s);
      lines.push(
        `${row.slug}: ${s.genomes.length} genomes, ${blocks.length} blocks (${blocks.filter((b) => b.named).length} named), ${summary.evaluated}/${summary.total} pairs evaluated, ${summary.undefined} undefined, ${summary.notEvaluated} not evaluated, break ${brk === null ? "none" : `below ${brk.below} above ${brk.outliers.map((o) => o.genomes.join("+")).join(", ")}`}`,
      );
    }
    console.log(lines.join("\n"));
  });

  it("gives M. tuberculosis 38 blocks, 173 of 231 pairs and its two outlier genomes", async () => {
    const manifest = parseManifest(await read("index.json"));
    const row = manifest.species.find((r) => r.slug === "mycobacterium_tuberculosis")!;
    expect([row.genomes_analyzed, row.genomes_labeled, row.datasets_reported, row.datasets_curated]).toEqual([2000, 32060, 22, 22]);
    const s = parseSpecies(await read("species/mycobacterium_tuberculosis.json"));
    const blocks = lineageBlocks(s);
    expect(blocks).toHaveLength(38);
    expect(blocks.filter((b) => b.named)).toHaveLength(9);
    expect(blocks.reduce((a, b) => a + blockSize(b), 0)).toBe(2000);
    expect(drugTotals(s, s.drugs.indexOf("rifampin"))).toEqual({ r: 1367, i: 0, s: 600, untested: 33 });
    expect(testedOf(drugTotals(s, s.drugs.indexOf("amoxicillin")))).toBe(57);
    const summary = pairSummary(crossPairs(s));
    expect([summary.evaluated, summary.total, summary.notEvaluated, summary.undefined]).toEqual([173, 231, 58, 5]);
    const brk = treeBreak(s)!;
    expect(brk.outliers.map((o) => o.genomes)).toEqual([["1773.41701"], ["1773.40949"]]);
  });
});
