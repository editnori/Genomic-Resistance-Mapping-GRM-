import { describe, expect, it } from "vitest";
import type { TreeNode } from "../schema";
import { dendrogram, treeLayout } from "./layout";

describe("dendrogram", () => {
  it("places leaves by leaf order and merges midway between their children", () => {
    // Genomes 0..3; leaf order puts genome 2 first. Merge 0 makes node 4 from
    // genomes 2 and 3; merge 1 makes node 5 from genomes 0 and 1; merge 2
    // joins nodes 4 and 5.
    const merges = [
      { a: 2, b: 3, height: -0.0, size: 2 },
      { a: 0, b: 1, height: 0.02, size: 2 },
      { a: 4, b: 5, height: 0.1, size: 4 },
    ];
    const d = dendrogram(merges, [2, 3, 0, 1]);
    expect([...d.leafX]).toEqual([2, 3, 0, 1]);
    expect(d.links.map((l) => [l.xa, l.xb, l.y])).toEqual([
      [0, 1, 0],
      [2, 3, 0.02],
      [0.5, 2.5, 0.1],
    ]);
    expect(d.links[2]!.ya).toBe(0);
    expect(d.links[2]!.yb).toBe(0.02);
    expect(d.maxHeight).toBe(0.1);
  });
});

const node = (depth: number, n: number, nR: number, split: TreeNode["split"]): TreeNode => ({ value: nR / n, n, n_resistant: nR, depth, split });

describe("treeLayout", () => {
  // Root splits into absent (1) and present (2); node 2 splits again.
  const nodes: TreeNode[] = [
    node(0, 100, 40, { feature: 7, gain: 0.2, absent: 1, present: 2 }),
    node(1, 60, 5, null),
    node(1, 40, 35, { feature: 9, gain: 0.05, absent: 3, present: 4 }),
    node(2, 10, 6, null),
    node(2, 30, 29, null),
  ];

  it("gives leaves consecutive columns, absent branch first, parents centered", () => {
    const { nodes: laid, columns } = treeLayout(nodes, 8);
    const x = Object.fromEntries(laid.map((l) => [l.id, l.x]));
    expect(columns).toBe(3);
    expect(x).toEqual({ 0: 0.75, 1: 0, 2: 1.5, 3: 1, 4: 2 });
  });

  it("folds splits below the shown depth into a leaf", () => {
    const { nodes: laid, columns } = treeLayout(nodes, 1);
    expect(columns).toBe(2);
    expect(laid.find((l) => l.id === 2)!.folded).toBe(true);
    expect(laid.some((l) => l.id === 3)).toBe(false);
  });
});
