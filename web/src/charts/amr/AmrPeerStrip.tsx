// Compact peer context strip: places one dataset's rate among its peer group
// (for example, all datasets in the same antibiotic class) with an IQR band,
// one tick per peer, a median line, and the active dataset's bold mark.

import { int, pct } from "../../data/format";

export interface AmrPeerStripProps {
  /** Metric name, e.g. "Lineage accuracy". */
  label: string;
  /** This dataset's value (0..1). */
  value: number;
  /** All peer values in the comparison pool (0..1), including this dataset. */
  peers: readonly number[];
  /** Human-readable name of the peer group, e.g. "carbapenem datasets". */
  poolName: string;
  /** Lower-is-better flips the good/bad tone (for VME/ME). Default false. */
  lowerIsBetter?: boolean;
}

const quantile = (sorted: readonly number[], q: number): number => {
  if (sorted.length === 0) return 0;
  const pos = (sorted.length - 1) * q;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  if (lo === hi) return sorted[lo]!;
  return sorted[lo]! + (sorted[hi]! - sorted[lo]!) * (pos - lo);
};

const WIDTH = 240;
const HEIGHT = 48;
const PAD = 10;
const TRACK_Y = 25;

export const AmrPeerStrip = ({ label, value, peers, poolName, lowerIsBetter = false }: AmrPeerStripProps) => {
  if (peers.length < 3) return null;
  const sorted = [...peers].sort((a, b) => a - b);
  const low = sorted[0]!;
  const high = sorted[sorted.length - 1]!;
  const p25 = quantile(sorted, 0.25);
  const med = quantile(sorted, 0.5);
  const p75 = quantile(sorted, 0.75);
  const minVal = Math.min(low, value);
  const maxVal = Math.max(high, value);
  let axisLo = minVal;
  let axisHi = maxVal;
  if (axisHi - axisLo < 0.05) {
    const mid = (axisLo + axisHi) / 2;
    axisLo = Math.max(0, mid - 0.025);
    axisHi = Math.min(1, axisLo + 0.05);
    axisLo = Math.max(0, axisHi - 0.05);
  }
  const span = Math.max(axisHi - axisLo, 0.01);
  const x = (v: number) => PAD + ((v - axisLo) / span) * (WIDTH - 2 * PAD);

  const good = lowerIsBetter ? value <= p25 : value >= p75;
  const bad = lowerIsBetter ? value >= p75 : value <= p25;
  const fill = good ? "rgb(var(--live-text))" : bad ? "rgb(var(--coral-text))" : "rgb(var(--blue))";
  const textClass = good ? "text-live-text" : bad ? "text-coral-text" : "text-ink";

  return (
    <figure data-amr-primitive="peer-strip" className="m-0 inline-block max-w-full">
      <svg
        width={WIDTH}
        height={HEIGHT}
        viewBox={`0 0 ${WIDTH} ${HEIGHT}`}
        role="img"
        aria-label={`${label} ${pct(value)} among ${int(sorted.length)} ${poolName}; median ${pct(med)}, range ${pct(low)} to ${pct(high)}.`}
        className="block max-w-full overflow-visible"
      >
        <rect
          x={x(p25)}
          y={TRACK_Y - 5}
          width={Math.max(1, x(p75) - x(p25))}
          height={10}
          rx={2}
          fill="rgb(var(--ink) / 0.09)"
        />
        <line x1={PAD} x2={WIDTH - PAD} y1={TRACK_Y} y2={TRACK_Y} stroke="var(--chart-axis)" strokeWidth={1} />
        {sorted.map((peer, idx) => (
          <line key={idx} x1={x(peer)} x2={x(peer)} y1={TRACK_Y - 3} y2={TRACK_Y + 3} stroke="rgb(var(--ink) / 0.35)" strokeWidth={1} />
        ))}
        <line x1={x(med)} x2={x(med)} y1={TRACK_Y - 7} y2={TRACK_Y + 7} stroke="rgb(var(--ink) / 0.65)" strokeWidth={1.5} />
        <rect x={x(value) - 2} y={TRACK_Y - 8} width={4} height={16} rx={1} fill={fill} />
        <text
          x={Math.min(Math.max(x(value), 22), WIDTH - 22)}
          y={10}
          textAnchor="middle"
          fontSize={11.5}
          fontWeight={600}
          fill={fill}
        >
          {pct(value)}
        </text>
        <text x={PAD} y={HEIGHT - 2} textAnchor="start" fontSize={11.5} fill="rgb(var(--muted-text))">
          {pct(axisLo)}
        </text>
        <text x={WIDTH - PAD} y={HEIGHT - 2} textAnchor="end" fontSize={11.5} fill="rgb(var(--muted-text))">
          {pct(axisHi)}
        </text>
      </svg>
      <figcaption className="mt-1 text-tiny text-muted-text">
        <span className={`font-semibold tabular-nums ${textClass}`}>
          {label} {pct(value)}
        </span>{" "}
        among {int(sorted.length)} {poolName} · median {pct(med)}
      </figcaption>
    </figure>
  );
};
