// The final CART tree as a node-link diagram. Each split names its feature;
// the absent branch goes left, the present branch right. A node's bar splits
// its isolates into resistant (coral) and susceptible (blue) with a 2px gap.

import { useState } from "react";
import type { FeatureCard, TreeNode } from "../data/schema";
import { int, pct } from "../data/format";
import { featureLabel, featurePlace } from "../data/derive/dataset";
import { treeLayout } from "../data/derive/layout";
import { ChartFrame, ChartTip } from "./kit";
import { ScrollFrame } from "../ui/primitives";

const NODE_W = 136;
const NODE_H = 50;
const COL = 150;
const ROW = 92;

const clip = (s: string, n = 19) => (s.length > n ? `${s.slice(0, n - 1)}…` : s);

export const TreeDiagram = ({ nodes, depth, features }: { nodes: readonly TreeNode[]; depth: number; features: Map<number, FeatureCard> }) => {
  const [hover, setHover] = useState<number | null>(null);
  const laid = treeLayout(nodes, depth);
  const maxDepth = Math.max(0, ...laid.nodes.map((n) => n.depth));
  const width = Math.max(laid.columns, 1) * COL + 8;
  const height = (maxDepth + 1) * ROW + 8;
  const pos = new Map(laid.nodes.map((n) => [n.id, { x: 4 + n.x * COL + (COL - NODE_W) / 2, y: 4 + n.depth * ROW }]));
  const h = hover === null ? undefined : laid.nodes.find((n) => n.id === hover);
  const hp = h === undefined ? undefined : pos.get(h.id);

  const label = (n: TreeNode): string => {
    if (n.split === null || n.depth >= depth) return n.value >= 0.5 ? "Calls resistant" : "Calls susceptible";
    return clip(featureLabel(features.get(n.split.feature), n.split.feature));
  };

  return (
    <ScrollFrame label="Decision tree">
      <ChartFrame>
        <svg width={width} height={height} role="img" aria-label={`Decision tree with ${laid.nodes.length} nodes shown.`} onPointerLeave={() => setHover(null)}>
          {laid.nodes.map((n) => {
            if (n.parent === null) return null;
            const p = pos.get(n.parent)!;
            const c = pos.get(n.id)!;
            const x1 = p.x + NODE_W / 2;
            const y1 = p.y + NODE_H;
            const x2 = c.x + NODE_W / 2;
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
            const rShare = n.node.n === 0 ? 0 : n.node.n_resistant / n.node.n;
            const barW = NODE_W - 20;
            const rW = Math.max(0, barW * rShare - (rShare > 0 && rShare < 1 ? 1 : 0));
            const sW = Math.max(0, barW - rW - (rShare > 0 && rShare < 1 ? 2 : 0));
            const leaf = n.node.split === null || n.folded;
            return (
              <g key={n.id} onPointerEnter={() => setHover(n.id)}>
                <rect x={x} y={y} width={NODE_W} height={NODE_H} rx={8} fill="rgb(var(--paper-pure))" stroke={hover === n.id ? "rgb(var(--blue))" : "var(--hair)"} strokeWidth={hover === n.id ? 1.5 : 1} />
                <text x={x + 10} y={y + 17} className={`text-[13px] ${leaf ? "fill-[rgb(var(--muted-text))]" : "fill-[rgb(var(--ink))] font-semibold"}`}>
                  {label(n.node)}
                </text>
                <text x={x + 10} y={y + 32} className="fill-[rgb(var(--muted-text))] text-[11.5px] tabular-nums">
                  {int(n.node.n)} isolates · {pct(rShare, 0)} R
                </text>
                <rect x={x + 10} y={y + 39} width={rW} height={4} rx={2} fill="var(--series-r)" />
                <rect x={x + 10 + rW + (rShare > 0 && rShare < 1 ? 2 : 0)} y={y + 39} width={sW} height={4} rx={2} fill="var(--series-s)" />
              </g>
            );
          })}
        </svg>
        {h !== undefined && hp !== undefined ? (
          <ChartTip
            x={hp.x + NODE_W - 14}
            y={hp.y}
            width={width}
            lines={[
              {
                term: h.node.split !== null && !h.folded ? `Split on ${featureLabel(features.get(h.node.split.feature), h.node.split.feature)}` : h.folded ? "Deeper splits hidden" : "Leaf",
                text: `${int(h.node.n)} training isolates, ${int(h.node.n_resistant)} resistant. Class-balanced resistant share ${pct(h.node.value)}.`,
              },
              ...(h.node.split !== null && !h.folded ? [{ source: featurePlace(features.get(h.node.split.feature)) ?? "Not placed on a reference assembly." }] : []),
            ]}
          />
        ) : null}
      </ChartFrame>
    </ScrollFrame>
  );
};
