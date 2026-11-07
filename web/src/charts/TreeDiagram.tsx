// The final CART tree as a node-link diagram. Each split names its feature;
// the absent branch goes left, the present branch right. A node's bar splits
// its isolates into resistant (coral) and susceptible (blue) with a 2px gap.
// MarkListbox gives the keyboard path: arrows walk the nodes top down, and
// click, tap or Enter pins a node's counts and split.

import { useEffect, useLayoutEffect, useRef, useState } from "react";
import type { FeatureCard, TreeNode } from "../data/schema";
import { int, pct } from "../data/format";
import { featureLabel, featurePlace } from "../data/derive/dataset";
import { treeLayout } from "../data/derive/layout";
import { MarkListbox, useWidth, type ListMark } from "./kit";
import { ScrollFrame } from "../ui/primitives";

const NODE_H = 50;
const ROW = 92;
// A column is 150 px where the frame has room, and narrows to 140 so the
// three-level view's eight leaves (1,128 px) fit the inspector's 1,136 px
// frame. Below that the tree scrolls inside its frame, as do deeper views.
const COL_MAX = 150;
const COL_MIN = 140;
const NODE_W_MAX = 136;
// At least 8 px between neighbouring nodes.
const NODE_GAP = 8;
// Text starts 10 px into a node and stops 10 px short of its edge, as the bar does.
const TEXT_INSET = 10;
// Before the first measurement (and in a server render) a split's label is
// cut by characters: 19 fit a 136 px node, fewer a narrower one.
const CHARS_PER_PX = 19 / NODE_W_MAX;
// The two lines of a node; the hidden probe measures with the same classes.
const SPLIT_CLASS = "text-[13px] fill-[rgb(var(--ink))] font-semibold";
const LEAF_CLASS = "text-[13px] fill-[rgb(var(--muted-text))]";
const COUNT_CLASS = "fill-[rgb(var(--muted-text))] text-[11.5px] tabular-nums";

const clip = (s: string, n: number) => (s.length > n ? `${s.slice(0, n - 1)}…` : s);

/** One line of a node's text: what it says in full, and a shorter form to
 *  use when the full one is wider than the node. */
interface Line {
  cls: string;
  full: string;
  short?: string;
}
const lineKey = (l: Line) => `${l.cls}\u0000${l.full}`;

/** The text that fits `room` px: the line in full, else its short form,
 *  else the longest prefix that fits with an ellipsis. */
const fitLine = (l: Line, room: number, width: (cls: string, s: string) => number): string => {
  if (width(l.cls, l.full) <= room) return l.full;
  if (l.short !== undefined) return l.short;
  let lo = 0;
  let hi = l.full.length - 1;
  while (lo < hi) {
    const mid = Math.ceil((lo + hi) / 2);
    if (width(l.cls, `${l.full.slice(0, mid).trimEnd()}…`) <= room) lo = mid;
    else hi = mid - 1;
  }
  return `${l.full.slice(0, lo).trimEnd()}…`;
};

export const TreeDiagram = ({ nodes, depth, features }: { nodes: readonly TreeNode[]; depth: number; features: Map<number, FeatureCard> }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const probe = useRef<SVGTextElement | null>(null);
  const avail = useWidth(frame, 0);
  const laid = treeLayout(nodes, depth);
  const maxDepth = Math.max(0, ...laid.nodes.map((n) => n.depth));
  const columns = Math.max(laid.columns, 1);
  const col = Math.max(COL_MIN, Math.min(COL_MAX, Math.floor((avail - 8) / columns)));
  const nodeW = Math.min(NODE_W_MAX, col - NODE_GAP);
  const chars = Math.floor(nodeW * CHARS_PER_PX);
  const width = columns * col + 8;
  const height = (maxDepth + 1) * ROW + 8;
  const pos = new Map(laid.nodes.map((n) => [n.id, { x: 4 + n.x * col + (col - nodeW) / 2, y: 4 + n.depth * ROW }]));

  const room = nodeW - 2 * TEXT_INSET;
  const share = (n: TreeNode) => (n.n === 0 ? 0 : n.n_resistant / n.n);
  const labelLine = (n: TreeNode): Line =>
    n.split === null || n.depth >= depth
      ? { cls: LEAF_CLASS, full: n.value >= 0.5 ? "Calls resistant" : "Calls susceptible" }
      : { cls: SPLIT_CLASS, full: featureLabel(features.get(n.split.feature), n.split.feature) };
  const countLine = (n: TreeNode): Line => ({ cls: COUNT_CLASS, full: `${int(n.n)} isolates · ${pct(share(n), 0)} R`, short: `${int(n.n)} · ${pct(share(n), 0)} R` });

  // Each line trimmed to its rendered width, measured before paint; the
  // node's tooltip and its listbox option keep the label in full.
  const [fonts, setFonts] = useState(0);
  const [fitted, setFitted] = useState<ReadonlyMap<string, string>>(new Map());
  useEffect(() => {
    let live = true;
    void document.fonts?.ready.then(() => live && setFonts((f) => f + 1));
    return () => {
      live = false;
    };
  }, []);
  const lines = laid.nodes.flatMap((h) => [labelLine(h.node), countLine(h.node)]);
  const linesKey = `${room}|${fonts}|${lines.map(lineKey).join("\u0001")}`;
  useLayoutEffect(() => {
    const t = probe.current;
    if (t === null) return;
    const width = (cls: string, s: string) => {
      t.setAttribute("class", cls);
      t.textContent = s;
      return t.getComputedTextLength();
    };
    setFitted(new Map(lines.map((l) => [lineKey(l), fitLine(l, room, width)])));
    // linesKey stands for lines, room and the fonts.
  }, [linesKey]);
  const shown = (l: Line, fallback: string) => fitted.get(lineKey(l)) ?? fallback;
  const label = (n: TreeNode) => {
    const l = labelLine(n);
    return shown(l, l.cls === SPLIT_CLASS ? clip(l.full, chars) : l.full);
  };
  const count = (n: TreeNode) => {
    const l = countLine(n);
    return shown(l, l.full);
  };

  // One mark per node, top down, at the node's top right corner.
  const marks: ListMark[] = laid.nodes.map((h) => {
    const p = pos.get(h.id)!;
    return {
      key: String(h.id),
      x: p.x + nodeW - 14,
      y: p.y,
      tip: () => [
        {
          term: h.node.split !== null && !h.folded ? `Split on ${featureLabel(features.get(h.node.split.feature), h.node.split.feature)}` : h.folded ? "Deeper splits hidden" : "Leaf",
          text: `${int(h.node.n)} training isolates, ${int(h.node.n_resistant)} resistant. Class-balanced resistant share ${pct(h.node.value)}.`,
        },
        ...(h.node.split !== null && !h.folded ? [{ source: featurePlace(features.get(h.node.split.feature)) ?? "Not placed on a reference assembly." }] : []),
      ],
    };
  });
  // A node is its whole box.
  const hit = (px: number, py: number) =>
    laid.nodes.findIndex((h) => {
      const p = pos.get(h.id)!;
      return px >= p.x && px <= p.x + nodeW && py >= p.y && py <= p.y + NODE_H;
    });

  return (
    <div ref={frame} className="relative">
    <svg aria-hidden="true" width={0} height={0} className="pointer-events-none absolute overflow-hidden" style={{ visibility: "hidden" }}>
      <text ref={probe} />
    </svg>
    <ScrollFrame label="Decision tree">
      <MarkListbox label={`Decision tree with ${laid.nodes.length} nodes shown. Arrow keys move between nodes; click or Enter to pin details.`} marks={marks} width={width} hit={hit}>
        <svg width={width} height={height} aria-hidden="true" className="block">
          {laid.nodes.map((n) => {
            if (n.parent === null) return null;
            const p = pos.get(n.parent)!;
            const c = pos.get(n.id)!;
            const x1 = p.x + nodeW / 2;
            const y1 = p.y + NODE_H;
            const x2 = c.x + nodeW / 2;
            const y2 = c.y;
            const my = (y1 + y2) / 2;
            return (
              <g key={`e${n.id}`}>
                <path d={`M${x1},${y1} C${x1},${my} ${x2},${my} ${x2},${y2}`} fill="none" stroke="var(--chart-axis)" strokeWidth={1.5} />
                <text x={(x1 + x2) / 2 + (n.branch === "absent" ? -6 : 6)} y={my + 4} textAnchor={n.branch === "absent" ? "end" : "start"} className="fill-[rgb(var(--muted-text))] text-[11.5px]">
                  {n.branch}
                </text>
              </g>
            );
          })}
          {laid.nodes.map((n) => {
            const { x, y } = pos.get(n.id)!;
            const rShare = share(n.node);
            const barW = nodeW - 20;
            const rW = Math.max(0, barW * rShare - (rShare > 0 && rShare < 1 ? 1 : 0));
            const sW = Math.max(0, barW - rW - (rShare > 0 && rShare < 1 ? 2 : 0));
            return (
              <g key={n.id} className="cursor-pointer">
                <rect x={x} y={y} width={nodeW} height={NODE_H} rx={8} fill="rgb(var(--paper-pure))" stroke="var(--hair)" strokeWidth={1} />
                <text x={x + TEXT_INSET} y={y + 17} className={labelLine(n.node).cls}>
                  {label(n.node)}
                </text>
                <text x={x + TEXT_INSET} y={y + 32} className={COUNT_CLASS}>
                  {count(n.node)}
                </text>
                <rect x={x + 10} y={y + 39} width={rW} height={4} rx={2} fill="var(--series-r)" />
                <rect x={x + 10 + rW + (rShare > 0 && rShare < 1 ? 2 : 0)} y={y + 39} width={sW} height={4} rx={2} fill="var(--series-s)" />
              </g>
            );
          })}
        </svg>
      </MarkListbox>
    </ScrollFrame>
    </div>
  );
};
