// Cross-resistance: phi between resistance calls for each drug pair, over
// genomes tested for both (pairs with fewer than 30 are blank). The lower
// triangle only, drugs grouped by class with a gap between classes, cells
// never under 24px (the frame scrolls instead of dropping labels), and a
// stepped diverging legend that counts the pairs in each step: coral when
// resistance to one goes with resistance to the other, blue when it goes
// against. A display of the existing phi values; no new analysis.

import { useState } from "react";
import type { SpeciesReport } from "../data/schema";
import { capitalize, int } from "../data/format";
import { ChartFrame, ChartTip, LABEL_CLASS } from "./kit";
import { ScrollFrame } from "../ui/primitives";

type Pair = SpeciesReport["cross"][number];

const STEPS = [
  { key: "strong-together", label: "0.5 or more", fill: "color-mix(in srgb, var(--series-r) 85%, transparent)", test: (phi: number) => phi >= 0.5 },
  { key: "together", label: "0.2 to 0.5", fill: "color-mix(in srgb, var(--series-r) 38%, transparent)", test: (phi: number) => phi >= 0.2 && phi < 0.5 },
  { key: "none", label: "−0.2 to 0.2", fill: "var(--fill-3)", test: (phi: number) => phi > -0.2 && phi < 0.2 },
  { key: "against", label: "−0.5 to −0.2", fill: "color-mix(in srgb, var(--series-s) 38%, transparent)", test: (phi: number) => phi > -0.5 && phi <= -0.2 },
  { key: "strong-against", label: "−0.5 or less", fill: "color-mix(in srgb, var(--series-s) 85%, transparent)", test: (phi: number) => phi <= -0.5 },
] as const;

const stepOf = (phi: number) => STEPS.find((s) => s.test(phi)) ?? STEPS[2];

const CELL = 24;
const GROUP_GAP = 6;
const LABEL = 176;

export const CrossMatrix = ({ report, order }: { report: SpeciesReport; order: readonly number[] }) => {
  const [hover, setHover] = useState<{ pair: Pair; x: number; y: number } | null>(null);
  const byKey = new Map(report.cross.flatMap((p) => [[`${p.a}|${p.b}`, p] as const, [`${p.b}|${p.a}`, p] as const]));
  // Classes in the order of their most-tested drug; drugs keep that order
  // inside a class.
  const groups = new Map<string, number[]>();
  order.forEach((j) => {
    const c = report.drug_classes[j] ?? "";
    groups.set(c, [...(groups.get(c) ?? []), j]);
  });
  const seq = [...groups.values()].flat();
  const groupIndex = new Map(seq.map((j) => [j, [...groups.keys()].indexOf(report.drug_classes[j] ?? "")]));
  const pos = (i: number) => i * CELL + (groupIndex.get(seq[i]!) ?? 0) * GROUP_GAP;
  const names = seq.map((j) => report.drugs[j] ?? "");
  const short = (s: string) => (s.length > 26 ? `${s.slice(0, 25)}…` : s);
  const n = names.length;
  const span = n === 0 ? 0 : pos(n - 1) + CELL;
  const top = 8;
  const bottom = top + span;
  const width = LABEL + span + 24;
  const height = bottom + 150;
  const counts = STEPS.map((s) => report.cross.filter((p) => s.test(p.phi)).length);
  return (
    <div>
      <ScrollFrame label="Cross-resistance matrix">
        <ChartFrame>
          <svg width={width} height={height} role="img" aria-label={`Phi correlation between resistance calls for ${report.cross.length} drug pairs, lower triangle, drugs grouped by class.`} onPointerLeave={() => setHover(null)}>
            {names.map((a, i) => (
              <g key={a}>
                {i > 0 ? (
                  <text x={LABEL - 6} y={top + pos(i) + CELL / 2 + 4} textAnchor="end" className={LABEL_CLASS}>
                    {short(capitalize(a))}
                  </text>
                ) : null}
                {i < n - 1 ? (
                  <text transform={`translate(${LABEL + pos(i) + CELL / 2 - 3} ${bottom + 8}) rotate(55)`} className={LABEL_CLASS}>
                    {short(capitalize(a))}
                  </text>
                ) : null}
                {names.slice(0, i).map((b, k) => {
                  const p = byKey.get(`${a}|${b}`);
                  const x = LABEL + pos(k);
                  const y = top + pos(i);
                  return (
                    <rect
                      key={b}
                      x={x + 1}
                      y={y + 1}
                      width={CELL - 2}
                      height={CELL - 2}
                      fill={p === undefined ? "var(--fill-1)" : stepOf(p.phi).fill}
                      stroke={hover?.pair === p && p !== undefined ? "rgb(var(--ink))" : "none"}
                      onPointerEnter={() => setHover(p === undefined ? null : { pair: p, x: x + CELL, y })}
                    />
                  );
                })}
              </g>
            ))}
          </svg>
          {hover !== null ? (
            <ChartTip
              x={hover.x}
              y={hover.y}
              width={width}
              lines={[
                {
                  term: `${capitalize(hover.pair.a)} and ${hover.pair.b}: phi ${hover.pair.phi.toFixed(2)}`,
                  text: `${int(hover.pair.n)} genomes tested for both: ${int(hover.pair.rr)} resistant to both, ${int(hover.pair.ss)} susceptible to both, ${int(hover.pair.rs + hover.pair.sr)} resistant to one.`,
                },
              ]}
            />
          ) : null}
        </ChartFrame>
      </ScrollFrame>
      <ul className="mt-3 flex flex-wrap gap-x-5 gap-y-1 text-control text-ink" aria-label="Phi steps">
        {STEPS.map((s, i) => (
          <li key={s.key} className="inline-flex items-center gap-1.5">
            <span aria-hidden="true" className="inline-block h-3 w-4 rounded-sm" style={{ background: s.fill }} />
            <span className="text-muted-text">Phi {s.label}</span>
            <span className="font-medium tabular-nums">{int(counts[i] ?? 0)}</span>
          </li>
        ))}
        <li className="inline-flex items-center gap-1.5">
          <span aria-hidden="true" className="inline-block h-3 w-4 rounded-sm bg-surface-subtle ring-1 ring-inset ring-hair" />
          <span className="text-muted-text">Fewer than 30 tested for both</span>
        </li>
      </ul>
    </div>
  );
};
