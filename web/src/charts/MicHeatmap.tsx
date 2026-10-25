// Measured against predicted MIC, one cell per pair of dilution steps,
// counted over isolates with an exact MIC. A one-hue ramp carries the count;
// the dashed diagonals bound essential agreement (within one doubling
// dilution). Counts print inside cells large enough to hold them.

import { useRef } from "react";
import { ScrollFrame } from "../ui/primitives";
import { int, mgL } from "../data/format";
import { MarkListbox, TICK_CLASS, useSvgId, useWidth } from "./kit";

export const MicHeatmap = ({ confusion }: { confusion: readonly [number, number, number][] }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const avail = useWidth(frame);
  const clip = useSvgId("mic-grid");
  if (confusion.length === 0) return null;
  const steps = confusion.flatMap(([o, p]) => [o, p]);
  const lo = Math.min(...steps);
  const hi = Math.max(...steps);
  const n = hi - lo + 1;
  const pad = { l: 56, t: 8, b: 44 };
  const max = Math.max(1, ...confusion.map((c) => c[2]));
  const cell = Math.max(24, String(max).length * 7 + 6, Math.min(34, Math.floor((Math.min(avail, 560) - pad.l - 8) / n)));
  const width = pad.l + cell * n + 8;
  const height = pad.t + cell * n + pad.b;
  const cuts = [...new Set([1, 2, 3, 4].map((i) => Math.ceil(max * i / 4)))];
  const level = (count: number) => count === 0 ? 0 : (cuts.findIndex((cut) => count <= cut) + 1) / cuts.length;
  const fill = (count: number) => {
    const share = level(count);
    // Separate low and high steps enough for their ink and paper labels
    // to stay legible in both themes, including the middle count bins.
    const opacity = count === 0 ? 0.04 : share <= 0.55 ? 0.06 + 0.36 * share : 0.5 + 0.45 * share;
    return `rgb(var(--ink) / ${opacity})`;
  };
  const cx = (p: number) => pad.l + (p - lo) * cell;
  const cy = (o: number) => pad.t + (hi - o) * cell;
  // Band edges |predicted - measured| = 1.5 steps, in pixel space.
  const edge = (offset: number) => {
    const pts: [number, number][] = [];
    for (let o = lo - 1; o <= hi + 1; o += 0.5) pts.push([cx(o + offset) + cell / 2, cy(o) + cell / 2]);
    return pts.map(([x, y], i) => `${i === 0 ? "M" : "L"}${x},${y}`).join(" ");
  };
  const values = Array.from({ length: n }, (_, i) => lo + i);
  return (
    <div ref={frame} className="min-w-0">
      <ScrollFrame label="Measured and predicted MIC">
      <MarkListbox label="Measured against predicted MIC counts. Arrow keys move between cells." width={width} marks={confusion.map(([o, p, count]) => ({ key: `${o},${p}`, x: cx(p) + cell / 2, y: cy(o) + cell / 2, tip: () => [{ term: `${int(count)} isolates`, text: `Measured ${mgL(o)} mg/L, predicted ${mgL(p)} mg/L.` }] }))}>
      <svg width={width} height={height} aria-hidden="true">
        <defs>
          <clipPath id={clip}>
            <rect x={pad.l} y={pad.t} width={cell * n} height={cell * n} />
          </clipPath>
        </defs>
        <rect x={pad.l} y={pad.t} width={cell * n} height={cell * n} fill="var(--fill-1)" />
        {confusion.map(([o, p, c]) => {
          const share = level(c);
          return (
            <g key={`${o},${p}`}>
              <rect x={cx(p) + 1} y={cy(o) + 1} width={cell - 2} height={cell - 2} rx={3} fill={fill(c)} />
              {cell >= 24 ? (
                <text x={cx(p) + cell / 2} y={cy(o) + cell / 2 + 4} textAnchor="middle" className="text-[12px] tabular-nums" fill={share > 0.55 ? "rgb(var(--paper-pure))" : "rgb(var(--ink))"}>
                  {c}
                </text>
              ) : null}
            </g>
          );
        })}
        <g clipPath={`url(#${clip})`} pointerEvents="none">
          <path d={edge(1.5)} fill="none" stroke="rgb(var(--ink) / 0.5)" strokeDasharray="3 3" />
          <path d={edge(-1.5)} fill="none" stroke="rgb(var(--ink) / 0.5)" strokeDasharray="3 3" />
        </g>
        {values.map((v) => (
          <g key={v}>
            <text x={pad.l - 6} y={cy(v) + cell / 2 + 4} textAnchor="end" className={TICK_CLASS}>
              {mgL(v)}
            </text>
            <text x={cx(v) + cell / 2} y={pad.t + cell * n + 15} textAnchor="middle" className={TICK_CLASS}>
              {cell >= 22 || (v - lo) % 2 === 0 ? mgL(v) : ""}
            </text>
          </g>
        ))}
        <text x={pad.l + (cell * n) / 2} y={height - 6} textAnchor="middle" className="fill-[rgb(var(--muted-text))] text-[12px]">
          Predicted MIC (mg/L)
        </text>
        <text transform={`translate(12 ${pad.t + (cell * n) / 2}) rotate(-90)`} textAnchor="middle" className="fill-[rgb(var(--muted-text))] text-[12px]">
          Measured MIC (mg/L)
        </text>
      </svg>
      </MarkListbox>
      </ScrollFrame>
      <div className="mt-3 flex flex-wrap gap-x-4 gap-y-2 text-tiny text-muted-text" aria-label="Cell counts">
        {cuts.map((cut, i) => <span key={cut} className="inline-flex items-center gap-1.5"><span className="size-4" aria-hidden="true" style={{ background: fill(cut) }} />{i === 0 ? 1 : cuts[i - 1]! + 1}–{cut} isolates</span>)}
      </div>
    </div>
  );
};
