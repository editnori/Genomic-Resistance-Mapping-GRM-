// Labeled genomes per collection year: one series of columns on a
// continuous year axis (missing years are gaps at zero), a 2px surface gap
// between columns, rounded data ends. MarkListbox gives the keyboard path:
// arrows walk the years, and click, tap or Enter pins one year's count.

import { useRef } from "react";
import { int } from "../data/format";
import { AxisY, linear, MarkListbox, TICK_CLASS, ticks, useWidth, type ListMark } from "./kit";

export const YearColumns = ({ years }: { years: readonly [string, number][] }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const numeric = years.flatMap(([y, n]) => (/^\d{4}$/.test(y) ? [[Number(y), n] as const] : []));
  const other = years.filter(([y]) => !/^\d{4}$/.test(y)).reduce((a, [, n]) => a + n, 0);
  if (numeric.length === 0) return null;
  const first = Math.min(...numeric.map(([y]) => y));
  const last = Math.max(...numeric.map(([y]) => y));
  const count = new Map(numeric);
  const max = Math.max(...numeric.map(([, n]) => n));
  const height = 200;
  const pad = { l: 48, r: 8, t: 8, b: 26 };
  const slots = last - first + 1;
  const step = (width - pad.l - pad.r) / slots;
  const barW = Math.max(1, Math.min(24, step - 2));
  const y = linear([0, max], [height - pad.b, pad.t]);
  const decade = slots > 40 ? 20 : slots > 15 ? 10 : 2;
  const labelYears = Array.from({ length: slots }, (_, i) => first + i).filter((yr) => yr % decade === 0);
  // One mark per year with genomes, at its column's top, in year order.
  const marks: ListMark[] = Array.from({ length: slots }, (_, i) => first + i).flatMap((yr) => {
    const n = count.get(yr) ?? 0;
    return n === 0 ? [] : [{ key: String(yr), x: pad.l + (yr - first + 0.5) * step, y: y(n), tip: () => [{ term: String(yr), text: `${int(n)} labeled genomes` }] }];
  });
  // A year is its whole column, not a point at its top.
  const hit = (px: number) => {
    const yr = first + Math.floor((px - pad.l) / step);
    return marks.findIndex((m) => m.key === String(yr));
  };
  return (
    <div ref={frame} className="min-w-0">
      <MarkListbox label={`Labeled genomes per collection year, ${first} to ${last}. Arrow keys move between years; click or Enter to pin details.`} marks={marks} width={width} hit={hit}>
        <svg width={width} height={height} aria-hidden="true" className="block">
          <AxisY y={y} x0={pad.l} x1={width - pad.r} values={ticks(0, max, 4)} format={(v) => (v >= 1000 ? `${Math.round(v / 100) / 10}K` : int(v))} />
          {Array.from({ length: slots }, (_, i) => {
            const yr = first + i;
            const n = count.get(yr) ?? 0;
            if (n === 0) return null;
            const x = pad.l + i * step + (step - barW) / 2;
            const h = Math.max(1, y(0) - y(n));
            const r = Math.min(4, barW / 2, h);
            return (
              <path
                key={yr}
                d={`M${x},${y(0)} V${y(0) - h + r} Q${x},${y(0) - h} ${x + r},${y(0) - h} H${x + barW - r} Q${x + barW},${y(0) - h} ${x + barW},${y(0) - h + r} V${y(0)} Z`}
                fill="var(--chart-mark)"
              />
            );
          })}
          <line x1={pad.l} x2={width - pad.r} y1={y(0)} y2={y(0)} stroke="var(--chart-axis)" />
          {labelYears.map((yr) => (
            <text key={yr} x={pad.l + (yr - first + 0.5) * step} y={height - 8} textAnchor="middle" className={TICK_CLASS}>
              {yr}
            </text>
          ))}
        </svg>
      </MarkListbox>
      {other > 0 ? <p className="text-tiny text-muted-text">{int(other)} genomes carry a year that is not a four-digit year and are left out.</p> : null}
    </div>
  );
};
