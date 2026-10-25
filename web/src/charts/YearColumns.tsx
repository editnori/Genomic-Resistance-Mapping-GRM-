// Labeled genomes per collection year: one series of columns on a
// continuous year axis (missing years are gaps at zero), a 2px surface gap
// between columns, rounded data ends, values on hover.

import { useRef, useState } from "react";
import { int } from "../data/format";
import { AxisY, ChartFrame, ChartTip, linear, TICK_CLASS, ticks, useWidth } from "./kit";

export const YearColumns = ({ years }: { years: readonly [string, number][] }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const [hover, setHover] = useState<number | null>(null);
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
  return (
    <ChartFrame frameRef={frame}>
      <svg
        width={width}
        height={height}
        role="img"
        aria-label={`Labeled genomes per collection year, ${first} to ${last}.`}
        onPointerMove={(e) => {
          const r = e.currentTarget.getBoundingClientRect();
          const i = Math.floor((e.clientX - r.left - pad.l) / step);
          setHover(i >= 0 && i < slots ? first + i : null);
        }}
        onPointerLeave={() => setHover(null)}
      >
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
              opacity={hover === null || hover === yr ? 1 : 0.55}
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
      {other > 0 ? <p className="text-tiny text-muted-text">{int(other)} genomes carry a year that is not a four-digit year and are left out.</p> : null}
      {hover !== null ? (
        <ChartTip x={pad.l + (hover - first + 0.5) * step} y={y(count.get(hover) ?? 0)} width={width} lines={[{ term: String(hover), text: `${int(count.get(hover) ?? 0)} labeled genomes` }]} />
      ) : null}
    </ChartFrame>
  );
};
