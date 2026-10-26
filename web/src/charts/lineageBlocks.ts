// Lineage blocks: the population figure's default columns. A species file
// lists its genomes in dendrogram leaf order, and each lineage (a cluster
// cut at the lineage threshold) is a contiguous run of that order. A lineage
// of at least BLOCK_MIN_GENOMES genomes is a block of its own; neighbouring
// smaller lineages pool into one block while the merge that joins each next
// genome to the pool stays within POOL_HEIGHT_FACTOR times the threshold.
// So 2,000 genomes read as tens of columns instead of 2,000 slivers.
// Presentation only: every count here tallies the profile calls the species
// file already carries, the four states the legend names.

import type { SpeciesReport } from "../data/schema";

export const BLOCK_MIN_GENOMES = 20;
export const POOL_HEIGHT_FACTOR = 10;

export interface LineageBlock {
  /** Column position, from 0. */
  index: number;
  /** Leaf-order positions [start, end). */
  start: number;
  end: number;
  /** The block's lineages in leaf order. */
  lineages: number[];
  /** One lineage of at least BLOCK_MIN_GENOMES genomes, rather than a pool. */
  named: boolean;
}

export const blockSize = (b: LineageBlock): number => b.end - b.start;

/** Genomes per lineage. */
export const lineageSizes = (report: SpeciesReport): Map<number, number> => {
  const sizes = new Map<number, number>();
  report.genomes.forEach((g) => sizes.set(g.lineage, (sizes.get(g.lineage) ?? 0) + 1));
  return sizes;
};

/** The merge height at each boundary of the leaf order: entry p is the
 *  height of the merge that joins the subtree ending at position p - 1 to the
 *  subtree starting at p (entry 0 is unused). Merge k creates node n + k, so
 *  each merge's later child starts at the boundary it closes. Heights below
 *  zero (float -0.0) clamp to zero. */
export const boundaryHeights = (report: SpeciesReport): Float64Array => {
  const n = report.genomes.length;
  const first = new Float64Array(n + report.merges.length);
  report.leaf_order.forEach((g, pos) => {
    first[g] = pos;
  });
  const out = new Float64Array(report.leaf_order.length);
  report.merges.forEach((m, k) => {
    const a = first[m.a] ?? 0;
    const b = first[m.b] ?? 0;
    first[n + k] = Math.min(a, b);
    const p = Math.max(a, b);
    if (p > 0 && p < out.length) out[p] = Math.max(0, m.height);
  });
  return out;
};

/** The figure's columns, in leaf order. Sizes sum to the genome count and
 *  the ranges tile the leaf order. */
export const lineageBlocks = (report: SpeciesReport): LineageBlock[] => {
  const sizes = lineageSizes(report);
  const heights = boundaryHeights(report);
  const cut = POOL_HEIGHT_FACTOR * report.lineage_threshold;
  const blocks: LineageBlock[] = [];
  report.leaf_order.forEach((g, pos) => {
    const lineage = report.genomes[g]?.lineage ?? -1;
    const named = (sizes.get(lineage) ?? 0) >= BLOCK_MIN_GENOMES;
    const last = blocks[blocks.length - 1];
    const joins = last !== undefined && (named ? last.named && last.lineages[0] === lineage : !last.named && (heights[pos] ?? 0) <= cut);
    if (last !== undefined && joins) {
      last.end = pos + 1;
      if (!last.lineages.includes(lineage)) last.lineages.push(lineage);
    } else {
      blocks.push({ index: blocks.length, start: pos, end: pos + 1, lineages: [lineage], named });
    }
  });
  return blocks;
};

// ---- Calls ------------------------------------------------------------------

/** The four profile states over some genomes for one drug. */
export interface CallCounts {
  r: number;
  i: number;
  s: number;
  untested: number;
}

/** Genomes with a call of any kind: resistant, intermediate or susceptible. */
export const testedOf = (c: CallCounts): number => c.r + c.i + c.s;

/** Genomes called R or S, the drug order's key and the result's denominator. */
export const rsOf = (c: CallCounts): number => c.r + c.s;

const add = (c: CallCounts, call: string | undefined): void => {
  if (call === "R") c.r += 1;
  else if (call === "I") c.i += 1;
  else if (call === "S") c.s += 1;
  else c.untested += 1;
};

const none = (): CallCounts => ({ r: 0, i: 0, s: 0, untested: 0 });

/** Calls for drug column `drug` over every genome in the file. */
export const drugTotals = (report: SpeciesReport, drug: number): CallCounts => {
  const c = none();
  report.genomes.forEach((g) => add(c, g.profile[drug]));
  return c;
};

/** Calls for drug column `drug` over leaf positions [start, end). */
export const rangeCalls = (report: SpeciesReport, start: number, end: number, drug: number): CallCounts => {
  const c = none();
  for (let pos = start; pos < end; pos += 1) add(c, report.genomes[report.leaf_order[pos]!]?.profile[drug]);
  return c;
};

/** Calls per block and drug column: entry [block][drug]. */
export const blockCalls = (report: SpeciesReport, blocks: readonly LineageBlock[]): CallCounts[][] =>
  blocks.map((b) => report.drugs.map((_, j) => rangeCalls(report, b.start, b.end, j)));

/** Drug columns with the most R or S calls first, ties by name: the order
 *  the species page has always used, so its default map drug stays put.
 *  Intermediate calls are counted beside it, never folded into it. */
export const drugsByCalls = (report: SpeciesReport): number[] => {
  const calls = report.drugs.map((_, j) => rsOf(drugTotals(report, j)));
  return report.drugs.map((_, j) => j).sort((a, b) => calls[b]! - calls[a]! || (report.drugs[a] ?? "").localeCompare(report.drugs[b] ?? ""));
};

/** A genome's drugs grouped by its call, in `order`. */
export const callsByState = (profile: string, order: readonly number[]): Record<keyof CallCounts, number[]> => {
  const out: Record<keyof CallCounts, number[]> = { r: [], i: [], s: [], untested: [] };
  order.forEach((j) => {
    const call = profile[j];
    out[call === "R" ? "r" : call === "I" ? "i" : call === "S" ? "s" : "untested"].push(j);
  });
  return out;
};

// ---- Lineages behind a drug's resistance --------------------------------------

export interface ResistantLineage {
  lineage: number;
  resistant: number;
  genomes: number;
}

/** Lineages by how many genomes resistant to `drug` they hold, most first,
 *  ties by lineage id; lineages with none are left out. */
export const resistantLineages = (report: SpeciesReport, drug: number): ResistantLineage[] => {
  const sizes = lineageSizes(report);
  const by = new Map<number, number>();
  report.genomes.forEach((g) => {
    if (g.profile[drug] === "R") by.set(g.lineage, (by.get(g.lineage) ?? 0) + 1);
  });
  return [...by.entries()]
    .map(([lineage, resistant]) => ({ lineage, resistant, genomes: sizes.get(lineage) ?? 0 }))
    .sort((a, b) => b.resistant - a.resistant || a.lineage - b.lineage);
};

/** The block holding a leaf position. */
export const blockAt = (blocks: readonly LineageBlock[], pos: number): LineageBlock | undefined => {
  let lo = 0;
  let hi = blocks.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    const b = blocks[mid]!;
    if (pos < b.start) hi = mid - 1;
    else if (pos >= b.end) lo = mid + 1;
    else return b;
  }
  return undefined;
};

/** The block holding a lineage's first genome. */
export const blockOfLineage = (report: SpeciesReport, blocks: readonly LineageBlock[], lineage: number): LineageBlock | undefined => {
  const pos = report.leaf_order.findIndex((g) => report.genomes[g]?.lineage === lineage);
  return pos < 0 ? undefined : blockAt(blocks, pos);
};

// ---- What a block records ------------------------------------------------------

export interface FieldCoverage {
  /** Genomes with the field recorded. */
  recorded: number;
  /** Distinct recorded values. */
  distinct: number;
}

export interface BlockRecords {
  mlst: FieldCoverage;
  country: FieldCoverage;
  year: FieldCoverage & { min: string | null; max: string | null };
}

/** How many of a block's genomes record MLST, country and year. */
export const blockRecords = (report: SpeciesReport, b: LineageBlock): BlockRecords => {
  const mlst = new Set<string>();
  const country = new Set<string>();
  const years: string[] = [];
  let m = 0;
  let c = 0;
  for (let pos = b.start; pos < b.end; pos += 1) {
    const g = report.genomes[report.leaf_order[pos]!];
    if (g === undefined) continue;
    if (g.mlst !== null) {
      m += 1;
      mlst.add(g.mlst);
    }
    if (g.country !== null) {
      c += 1;
      country.add(g.country);
    }
    if (g.year !== null) years.push(g.year);
  }
  const sorted = [...years].sort();
  return {
    mlst: { recorded: m, distinct: mlst.size },
    country: { recorded: c, distinct: country.size },
    year: { recorded: years.length, distinct: new Set(years).size, min: sorted[0] ?? null, max: sorted[sorted.length - 1] ?? null },
  };
};

// ---- The dendrogram over blocks and its axis break ---------------------------

export interface BlockLink {
  /** Block coordinates: a block's centre is its index + 0.5. */
  xa: number;
  xb: number;
  ya: number;
  yb: number;
  y: number;
}

/** Dendrogram links with every genome at its block's centre, so merges
 *  inside a block draw as one stem and merges between blocks keep their
 *  heights. A merge sits midway between its children. */
export const blockDendrogram = (report: SpeciesReport, blocks: readonly LineageBlock[]): BlockLink[] => {
  const n = report.genomes.length;
  const x = new Float64Array(n + report.merges.length);
  const y = new Float64Array(n + report.merges.length);
  blocks.forEach((b) => {
    for (let pos = b.start; pos < b.end; pos += 1) x[report.leaf_order[pos]!] = b.index + 0.5;
  });
  return report.merges.map((m, k) => {
    const h = Math.max(0, m.height);
    x[n + k] = ((x[m.a] ?? 0) + (x[m.b] ?? 0)) / 2;
    y[n + k] = h;
    return { xa: x[m.a] ?? 0, xb: x[m.b] ?? 0, ya: y[m.a] ?? 0, yb: y[m.b] ?? 0, y: h };
  });
};

/** A top merge is an outlier when it sits at least OUTLIER_GAP times higher
 *  than the next merge down and joins at most OUTLIER_MAX_GENOMES genomes to
 *  the rest; the scan looks at the OUTLIER_SCAN highest merges. */
export const OUTLIER_GAP = 3;
export const OUTLIER_MAX_GENOMES = 5;
export const OUTLIER_SCAN = 10;

export interface OutlierMerge {
  height: number;
  /** Ids of the genomes on the merge's smaller side. */
  genomes: string[];
}

export interface TreeBreak {
  /** The highest merge drawn below the break. */
  below: number;
  /** Merges drawn above the break, highest first. */
  outliers: OutlierMerge[];
}

/** Where the tree's height axis breaks, so a few genomes that join far above
 *  everything else do not flatten every other merge into the baseline. Null
 *  when no top merge is an outlier, or when one that clears the gap joins a
 *  large clade: that is structure, not an outlier. */
export const treeBreak = (report: SpeciesReport): TreeBreak | null => {
  const n = report.genomes.length;
  const order = report.merges.map((m, k) => ({ k, h: Math.max(0, m.height) })).sort((a, b) => b.h - a.h || a.k - b.k);
  let cut = -1;
  for (let i = 0; i < Math.min(OUTLIER_SCAN, order.length - 1); i += 1) {
    const next = order[i + 1]!.h;
    if (next > 0 && order[i]!.h >= OUTLIER_GAP * next) {
      cut = i;
      break;
    }
  }
  if (cut < 0) return null;
  const sizeOf = (id: number): number => (id < n ? 1 : (report.merges[id - n]?.size ?? 1));
  const leaves = (id: number): number[] => {
    const out: number[] = [];
    const stack = [id];
    while (stack.length > 0 && out.length <= OUTLIER_MAX_GENOMES) {
      const cur = stack.pop()!;
      const m = cur < n ? undefined : report.merges[cur - n];
      if (m === undefined) out.push(cur);
      else stack.push(m.b, m.a);
    }
    return out;
  };
  const outliers: OutlierMerge[] = [];
  for (const { k, h } of order.slice(0, cut + 1)) {
    const m = report.merges[k]!;
    const small = sizeOf(m.a) <= sizeOf(m.b) ? m.a : m.b;
    if (sizeOf(small) > OUTLIER_MAX_GENOMES) return null;
    outliers.push({ height: h, genomes: leaves(small).map((g) => report.genomes[g]?.id ?? String(g)) });
  }
  return { below: order[cut + 1]!.h, outliers };
};
