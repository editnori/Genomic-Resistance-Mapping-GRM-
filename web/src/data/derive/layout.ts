// Geometry for the two tree figures: the species dendrogram (scipy-style
// merges over genomes) and the CART decision tree.

import type { Merge, TreeNode } from "../schema";

export interface DendroLink {
  /** Leaf-order positions of the two children (fractional for subtrees). */
  xa: number;
  xb: number;
  /** Heights of the two children and of the merge itself. */
  ya: number;
  yb: number;
  y: number;
}

export interface Dendrogram {
  /** Position of genome i along the leaf axis (0..n-1). */
  leafX: Float64Array;
  links: DendroLink[];
  maxHeight: number;
}

/** Merge k creates node n + k from nodes a and b (ids below n are leaves).
 *  Leaves sit at their leaf-order index; a merge sits midway between its
 *  children at its height. Heights below zero (float -0.0) clamp to zero. */
export const dendrogram = (merges: readonly Merge[], leafOrder: readonly number[]): Dendrogram => {
  const n = leafOrder.length;
  const x = new Float64Array(n + merges.length);
  const y = new Float64Array(n + merges.length);
  leafOrder.forEach((genome, pos) => {
    x[genome] = pos;
  });
  const leafX = x.slice(0, n);
  let maxHeight = 0;
  const links = merges.map((m, k): DendroLink => {
    const id = n + k;
    const h = Math.max(0, m.height);
    x[id] = (x[m.a]! + x[m.b]!) / 2;
    y[id] = h;
    maxHeight = Math.max(maxHeight, h);
    return { xa: x[m.a]!, xb: x[m.b]!, ya: y[m.a]!, yb: y[m.b]!, y: h };
  });
  return { leafX, links, maxHeight };
};

export interface TreeLayoutNode {
  id: number;
  node: TreeNode;
  /** Column position: leaves take 0, 1, 2, ...; parents center on children. */
  x: number;
  depth: number;
  parent: number | null;
  /** Which branch of the parent this node is. */
  branch: "absent" | "present" | null;
  /** True when the node has a split below the shown depth. */
  folded: boolean;
}

/** Tidy layout of a CART tree down to `maxDepth`, absent branch on the left. */
export const treeLayout = (nodes: readonly TreeNode[], maxDepth: number): { nodes: TreeLayoutNode[]; columns: number } => {
  const out: TreeLayoutNode[] = [];
  let nextLeaf = 0;
  const visit = (id: number, parent: number | null, branch: "absent" | "present" | null): number => {
    const node = nodes[id];
    if (node === undefined) return nextLeaf;
    const at = out.length;
    out.push({ id, node, x: 0, depth: node.depth, parent, branch, folded: false });
    const split = node.split;
    if (split === null || node.depth >= maxDepth) {
      out[at]!.x = nextLeaf;
      out[at]!.folded = split !== null;
      nextLeaf += 1;
      return out[at]!.x;
    }
    const left = visit(split.absent, id, "absent");
    const right = visit(split.present, id, "present");
    out[at]!.x = (left + right) / 2;
    return out[at]!.x;
  };
  if (nodes.length > 0) visit(0, null, null);
  return { nodes: out, columns: nextLeaf };
};
