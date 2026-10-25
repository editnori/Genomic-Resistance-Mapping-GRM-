// One row of the learner comparison: a balanced-accuracy strip on the axis
// shared by every row. Pooled random split is a filled dot, lineage held
// out is a ring of the same hue, single folds are faint marks of the same
// shape, and the catalog baseline is a dashed line through every row.

import type { Linear } from "./kit";
import { AxisX, ticks } from "./kit";

export const STRIP_HEIGHT = 30;

export const StripAxis = ({ x, width }: { x: Linear; width: number }) => (
  <svg width={width} height={22} aria-hidden="true" className="block overflow-visible">
    <g transform="translate(0,-2)">
      <AxisX x={x} y={6} height={0} values={ticks(x.domain[0], x.domain[1], 5)} format={(v) => `${Math.round(v * 100)}`} />
    </g>
  </svg>
);

export const LearnerStripRow = ({
  x,
  width,
  random,
  lineage,
  foldsRandom,
  foldsLineage,
  catalog,
  label,
}: {
  x: Linear;
  width: number;
  random: number | null;
  lineage: number | null;
  foldsRandom: readonly number[];
  foldsLineage: readonly number[];
  catalog: number | null;
  label: string;
}) => {
  const mid = STRIP_HEIGHT / 2;
  return (
    <svg width={width} height={STRIP_HEIGHT} role="img" aria-label={label} className="block">
      {ticks(x.domain[0], x.domain[1], 5).map((t) => (
        <line key={t} x1={x(t)} x2={x(t)} y1={0} y2={STRIP_HEIGHT} stroke="var(--chart-grid)" />
      ))}
      {catalog !== null ? <line x1={x(catalog)} x2={x(catalog)} y1={0} y2={STRIP_HEIGHT} stroke="rgb(var(--ink) / 0.55)" strokeDasharray="3 3" /> : null}
      {foldsRandom.map((v, i) => (
        <circle key={`r${i}`} cx={x(v)} cy={mid - 6} r={3} fill="var(--chart-mark)" opacity={0.35} />
      ))}
      {foldsLineage.map((v, i) => (
        <circle key={`l${i}`} cx={x(v)} cy={mid + 6} r={2.5} fill="none" stroke="var(--chart-mark)" strokeWidth={1.5} opacity={0.45} />
      ))}
      {random !== null && lineage !== null ? <line x1={x(random)} x2={x(lineage)} y1={mid} y2={mid} stroke="var(--chart-mark)" strokeWidth={2} opacity={0.4} /> : null}
      {lineage !== null ? <circle cx={x(lineage)} cy={mid} r={4.5} fill="rgb(var(--paper-pure))" stroke="var(--chart-mark)" strokeWidth={2} /> : null}
      {random !== null ? <circle cx={x(random)} cy={mid} r={5} fill="var(--chart-mark)" stroke="rgb(var(--paper-pure))" strokeWidth={2} /> : null}
    </svg>
  );
};
