// Measured against predicted MIC, one cell per pair of dilution steps,
// counted over isolates with an exact MIC. A one-hue ramp carries the count;
// the dashed diagonals bound essential agreement (within one doubling
// dilution). Counts print inside cells large enough to hold them.

import { useRef } from "react";
import { ScrollFrame } from "../ui/primitives";
import { int, mgL, pct } from "../data/format";
import { MarkListbox, TICK_CLASS, tickWidth, useSvgId, useWidth } from "./kit";

/** The band's share of exact MICs for the model and for a constant MIC from
 *  the training folds, when a baselines file gives the constant. */
export interface BandShares {
  model: number;
  constant: number;
  n: number;
}

export const MicHeatmap = ({ confusion, band }: { confusion: readonly [number, number, number][]; band?: BandShares }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const avail = useWidth(frame);
  const clip = useSvgId("mic-grid");
  if (confusion.length === 0) return null;
  const steps = confusion.flatMap(([o, p]) => [o, p]);
  const lo = Math.min(...steps);
  const hi = Math.max(...steps);
  const n = hi - lo + 1;
  const values = Array.from({ length: n }, (_, i) => lo + i);
  const maxTickW = Math.max(0, ...values.map((v) => tickWidth(mgL(v))));
  const lastTickW = values.length === 0 ? 0 : tickWidth(mgL(values[values.length - 1]!));
  const pad = { l: Math.max(46, maxTickW + 12), r: Math.max(12, Math.ceil(lastTickW / 2) + 4), t: 24, b: 44 };
  const max = Math.max(1, ...confusion.map((c) => c[2]));
  const cell = Math.max(24, String(max).length * 7 + 6, Math.min(34, Math.floor((Math.min(avail, 560) - pad.l - pad.r) / n)));
  const width = pad.l + cell * n + pad.r;
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
  // Thin X-axis dilution labels greedily so multi-character values such as
  // 0.008 or 0.125 never collide when cells are narrow, while narrower labels
  // (1, 2, 4, 8) keep printing every step.
  let prevXRight = -Infinity;
  const xTickShown = values.map((v) => {
    const s = mgL(v);
    const hw = tickWidth(s) / 2;
    const center = cx(v) + cell / 2;
    if (center - hw < prevXRight + 5) return "";
    prevXRight = center + hw;
    return s;
  });
  return (
    <div ref={frame} className="min-w-0">
      <ScrollFrame label="Measured and predicted MIC">
      <MarkListbox label="Measured against predicted MIC counts. Arrow keys move between cells; click or Enter to pin details." width={width} marks={confusion.map(([o, p, count]) => ({ key: `${o},${p}`, x: cx(p) + cell / 2, y: cy(o) + cell / 2, tip: () => [{ term: `${int(count)} isolates`, text: `Measured ${mgL(o)} mg/L, predicted ${mgL(p)} mg/L.` }] }))}>
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
        {values.map((v, i) => (
          <g key={v}>
            <text x={pad.l - 6} y={cy(v) + cell / 2 + 4} textAnchor="end" className={TICK_CLASS}>
              {mgL(v)}
            </text>
            <text x={cx(v) + cell / 2} y={pad.t + cell * n + 16} textAnchor="middle" className={TICK_CLASS}>
              {xTickShown[i]}
            </text>
          </g>
        ))}
        <text x={pad.l + (cell * n) / 2} y={height - 6} textAnchor="middle" className="fill-[rgb(var(--muted-text))] text-[12px]">
          Predicted MIC (mg/L)
        </text>
        <text x={pad.l} y={14} className="fill-[rgb(var(--muted-text))] text-[12px]">
          Measured MIC (mg/L)
        </text>
      </svg>
      </MarkListbox>
      </ScrollFrame>
      {/* The key wraps at the chart's width (w-0 min-w-full): unwrapped, its
          one line set the figure's column at about 925 px and left the text
          beside it 0 px wide at a 1024 px window. */}
      <div className="mt-3 flex w-0 min-w-full flex-wrap items-center gap-x-4 gap-y-2 text-tiny text-muted-text" aria-label="Heatmap key">
        {cuts.map((cut, i) => <span key={cut} className="inline-flex items-center gap-1.5"><span className="size-4" aria-hidden="true" style={{ background: fill(cut) }} />{i === 0 ? 1 : cuts[i - 1]! + 1}–{cut} isolates</span>)}
        <span className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={16} height={12} className="shrink-0">
            <line x1={0} x2={16} y1={6} y2={6} stroke="rgb(var(--ink) / 0.5)" strokeDasharray="3 3" />
          </svg>
          Between the dashed lines: predicted within one doubling dilution of the measured MIC
        </span>
      </div>
      {band === undefined ? null : (
        // Its own line under the key, at zero width to the grid's column
        // sizing (w-0 min-w-full), so the figure's column keeps the width the
        // chart gives it and the text beside it keeps its measure.
        <p className="mt-2 w-0 min-w-full text-tiny text-muted-text" data-baseline="band">
          Between the dashed lines: {pct(band.model)} of {int(band.n)} exact MICs for the model, {pct(band.constant)} for a constant MIC from the training folds.
        </p>
      )}
    </div>
  );
};
