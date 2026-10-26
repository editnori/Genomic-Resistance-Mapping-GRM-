import { describe, expect, it } from "vitest";
import type { GenomeRow, Merge, SpeciesReport } from "../data/schema";
import { fixtureFile } from "../dev/fixtures/files";
import { parseSpecies } from "../data/validate";
import {
  BLOCK_MIN_GENOMES,
  blockAt,
  blockCalls,
  blockDendrogram,
  blockSize,
  boundaryHeights,
  callsByState,
  drugTotals,
  drugsByCalls,
  lineageBlocks,
  resistantLineages,
  treeBreak,
} from "./lineageBlocks";

const genome = (id: string, lineage: number, profile: string): GenomeRow => ({
  id,
  lineage,
  x: 0,
  y: 0,
  nearest: 0,
  country: null,
  year: null,
  host: null,
  mlst: null,
  length: 0,
  contigs: 0,
  n50: 0,
  groups: [],
  profile,
});

/** Genomes in leaf order 0..n-1 joined left to right: merge k joins
 *  everything before leaf k + 1 to that leaf at heights[k], so the boundary
 *  before leaf p sits at heights[p - 1]. */
const chain = (lineages: readonly number[], heights: readonly number[], profiles?: readonly string[], threshold = 0.001): SpeciesReport => {
  const n = lineages.length;
  const merges: Merge[] = heights.map((height, k) => ({ a: k === 0 ? 0 : n + k - 1, b: k + 1, height, size: k + 2 }));
  return {
    schema: 1,
    species: "Test species",
    slug: "test",
    drugs: ["amikacin", "colistin"],
    drug_classes: ["aminoglycoside", "polymyxin"],
    catalog_groups: [],
    lineage_threshold: threshold,
    n_lineages: new Set(lineages).size,
    leaf_order: lineages.map((_, i) => i),
    merges,
    genomes: lineages.map((l, i) => genome(`G${i}`, l, profiles?.[i] ?? "RS")),
    cross: [],
  };
};

const runs = (spec: readonly [lineage: number, size: number][]): number[] => spec.flatMap(([l, size]) => Array.from({ length: size }, () => l));

describe("lineageBlocks", () => {
  it("makes each lineage of at least 20 genomes a block and pools small neighbours only within 10x the threshold", () => {
    // Lineage 1 (25 genomes), lineages 2-6 (3 genomes each), lineage 7 (20).
    const lineages = runs([
      [1, 25],
      [2, 3],
      [3, 3],
      [4, 3],
      [5, 3],
      [6, 3],
      [7, BLOCK_MIN_GENOMES],
    ]);
    // Small merges inside lineages; joins between neighbours at 0.005 (within
    // 10 x 0.001), except the boundary into lineage 5 at 0.02.
    const heights = lineages.slice(1).map((l, k) => (l === lineages[k] ? 0.0005 : l === 5 ? 0.02 : 0.005));
    const blocks = lineageBlocks(chain(lineages, heights));
    expect(blocks.map((b) => [b.start, b.end, b.named, b.lineages])).toEqual([
      [0, 25, true, [1]],
      [25, 34, false, [2, 3, 4]],
      [34, 40, false, [5, 6]],
      [40, 60, true, [7]],
    ]);
    expect(blocks.map((b) => b.index)).toEqual([0, 1, 2, 3]);
  });

  it("pools at exactly 10x the threshold and splits just above it", () => {
    const lineages = runs([
      [1, 2],
      [2, 2],
      [3, 2],
    ]);
    const at = lineageBlocks(chain(lineages, [0, 0.01, 0, 0.01, 0]));
    expect(at.map((b) => b.lineages)).toEqual([[1, 2, 3]]);
    const above = lineageBlocks(chain(lineages, [0, 0.0101, 0, 0.01, 0]));
    expect(above.map((b) => b.lineages)).toEqual([[1], [2, 3]]);
  });

  it("never pools a large lineage, however close its neighbours", () => {
    const lineages = runs([
      [1, 2],
      [2, BLOCK_MIN_GENOMES],
      [3, 2],
    ]);
    const blocks = lineageBlocks(chain(lineages, lineages.slice(1).map(() => 0)));
    expect(blocks.map((b) => [blockSize(b), b.named])).toEqual([
      [2, false],
      [BLOCK_MIN_GENOMES, true],
      [2, false],
    ]);
  });

  it("tiles the leaf order of both fixture species and keeps every call", () => {
    for (const slug of ["klebsiella_pneumoniae", "staphylococcus_aureus"]) {
      const s = parseSpecies(fixtureFile(`/data/species/${slug}.json`));
      const blocks = lineageBlocks(s);
      expect(blocks[0]!.start).toBe(0);
      blocks.forEach((b, i) => expect(b.start).toBe(i === 0 ? 0 : blocks[i - 1]!.end));
      expect(blocks[blocks.length - 1]!.end).toBe(s.genomes.length);
      const calls = blockCalls(s, blocks);
      s.drugs.forEach((_, j) => {
        // Independent tally straight from the genome rows.
        const want = { r: 0, i: 0, s: 0, untested: 0 };
        s.genomes.forEach((g) => {
          const c = g.profile[j];
          if (c === "R") want.r += 1;
          else if (c === "I") want.i += 1;
          else if (c === "S") want.s += 1;
          else want.untested += 1;
        });
        const got = calls.reduce((a, row) => ({ r: a.r + row[j]!.r, i: a.i + row[j]!.i, s: a.s + row[j]!.s, untested: a.untested + row[j]!.untested }), { r: 0, i: 0, s: 0, untested: 0 });
        expect(got).toEqual(want);
        expect(drugTotals(s, j)).toEqual(want);
      });
    }
  });
});

describe("boundaryHeights", () => {
  it("gives each leaf-order boundary the height of the merge that closes it", () => {
    // ((0,1) at 0.1, (2,3) at 0.2) at 0.5, leaf order 0 1 2 3.
    const s = chain([1, 1, 2, 2], [0, 0, 0]);
    s.merges = [
      { a: 0, b: 1, height: 0.1, size: 2 },
      { a: 2, b: 3, height: 0.2, size: 2 },
      { a: 4, b: 5, height: 0.5, size: 4 },
    ];
    expect([...boundaryHeights(s)]).toEqual([0, 0.1, 0.5, 0.2]);
  });

  it("follows the leaf order, not genome indices", () => {
    const s = chain([1, 1, 2], [0, 0]);
    s.leaf_order = [2, 0, 1];
    s.merges = [
      { a: 0, b: 1, height: 0.01, size: 2 },
      { a: 2, b: 3, height: 0.05, size: 3 },
    ];
    expect([...boundaryHeights(s)]).toEqual([0, 0.05, 0.01]);
  });
});

describe("treeBreak", () => {
  it("breaks above single genomes that join far above every other merge and names them", () => {
    // Two outliers joined last at 0.26 and 0.20; everything else at most 0.0025.
    const lineages = runs([
      [1, 6],
      [2, 1],
      [3, 1],
    ]);
    const heights = [0.0001, 0.0002, 0.001, 0.002, 0.0025, 0.2, 0.26];
    const b = treeBreak(chain(lineages, heights));
    expect(b).not.toBeNull();
    expect(b!.below).toBe(0.0025);
    expect(b!.outliers).toEqual([
      { height: 0.26, genomes: ["G7"] },
      { height: 0.2, genomes: ["G6"] },
    ]);
  });

  it("has no break without a threefold gap at the top", () => {
    const lineages = runs([[1, 8]]);
    expect(treeBreak(chain(lineages, [0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.05]))).toBeNull();
  });

  it("has no break when the distant merge joins a large clade", () => {
    // The top merge joins leaves 0-6 to leaves 7-13 at 0.3; the rest stay low.
    const s = chain(runs([[1, 14]]), Array.from({ length: 13 }, () => 0.001));
    const n = 14;
    const merges: Merge[] = [];
    for (let k = 0; k < 6; k += 1) merges.push({ a: k === 0 ? 0 : n + k - 1, b: k + 1, height: 0.001 * (k + 1), size: k + 2 });
    for (let k = 0; k < 6; k += 1) merges.push({ a: k === 0 ? 7 : n + 6 + k - 1, b: 8 + k, height: 0.001 * (k + 1), size: k + 2 });
    merges.push({ a: n + 5, b: n + 11, height: 0.3, size: 14 });
    s.merges = merges;
    expect(treeBreak(s)).toBeNull();
  });
});

describe("blockDendrogram", () => {
  it("puts every genome at its block's centre and keeps merge heights", () => {
    const lineages = runs([
      [1, BLOCK_MIN_GENOMES],
      [2, BLOCK_MIN_GENOMES],
    ]);
    const heights = lineages.slice(1).map((l, k) => (l === lineages[k] ? 0.0001 : 0.3));
    const s = chain(lineages, heights);
    const links = blockDendrogram(s, lineageBlocks(s));
    // Merges inside block 0 stand at its centre; the join of the blocks spans both centres.
    expect(links[0]).toMatchObject({ xa: 0.5, xb: 0.5, y: 0.0001 });
    const join = links[BLOCK_MIN_GENOMES - 1]!;
    expect(join.y).toBe(0.3);
    expect(join.xa).toBe(0.5);
    expect(join.xb).toBe(1.5);
  });
});

describe("drug summaries", () => {
  const s = chain([1, 1, 2, 2, 3], [0, 0, 0, 0], ["RS", "RI", "R.", "S.", "RR"]);

  it("orders drugs by R or S calls, the long-standing order, leaving intermediate calls out", () => {
    // Amikacin: 3 R or S. Colistin: 1 R or S and 3 intermediate, so it would
    // lead if intermediate calls counted.
    const t = chain([1, 1, 1, 1, 1], [0, 0, 0, 0], ["RI", "SI", "RI", ".S", ".."]);
    expect(drugsByCalls(t)).toEqual([0, 1]);
    expect(drugTotals(t, 1)).toEqual({ r: 0, i: 3, s: 1, untested: 1 });
    expect(drugTotals(s, 1)).toEqual({ r: 1, i: 1, s: 1, untested: 2 });
  });

  it("ranks lineages by resistant genomes, ties by id", () => {
    expect(resistantLineages(s, 0)).toEqual([
      { lineage: 1, resistant: 2, genomes: 2 },
      { lineage: 2, resistant: 1, genomes: 2 },
      { lineage: 3, resistant: 1, genomes: 1 },
    ]);
  });

  it("groups a genome's drugs by call", () => {
    expect(callsByState("R.SI", [3, 2, 1, 0])).toEqual({ r: [0], i: [3], s: [2], untested: [1] });
  });

  it("finds the block holding a position", () => {
    const blocks = lineageBlocks(chain(runs([[1, 3]]), [0.5, 0.5]));
    expect(blocks.map(blockSize)).toEqual([1, 1, 1]);
    expect(blockAt(blocks, 0)?.index).toBe(0);
    expect(blockAt(blocks, 2)?.index).toBe(2);
    expect(blockAt(blocks, 3)).toBeUndefined();
  });
});
