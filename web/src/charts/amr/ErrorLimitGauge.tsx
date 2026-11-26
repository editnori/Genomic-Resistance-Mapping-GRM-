// Compact horizontal gauges for VME (vs 1.5% limit) and ME (vs 3.0% limit)
// with exact 95% confidence interval whiskers when available. Clearly separates
// shown-to-meet (green), shown-to-miss (coral), and interval-straddling
// (undecided neutral).

import { pct } from "../../data/format";
import { BARS } from "../../data/metrics";

export interface ErrorRateInput {
  rate: number | null;
  low?: number | null;
  high?: number | null;
}

export interface ErrorLimitGaugeProps {
  vme: ErrorRateInput;
  me: ErrorRateInput;
  splitLabel: string;
}

interface GaugeRowProps {
  name: string;
  limit: number;
  limitText: string;
  input: ErrorRateInput;
}

const GaugeRow = ({ name, limit, limitText, input }: GaugeRowProps) => {
  if (input.rate === null) {
    return (
      <div className="min-w-0">
        <div className="flex items-baseline justify-between gap-2 text-control">
          <span className="font-medium text-ink">{name}</span>
          <span className="tabular-nums text-muted-text">No rate</span>
        </div>
      </div>
    );
  }

  const rate = input.rate;
  const lo = input.low ?? null;
  const hi = input.high ?? null;
  const hasCi = lo !== null && hi !== null;
  const scaleMax = Math.max(0.15, limit * 4, rate * 1.15, (hi ?? rate) * 1.08);
  const pos = (v: number) => `${Math.max(0, Math.min(100, (v / scaleMax) * 100)).toFixed(2)}%`;

  const status: "meets" | "misses" | "open" | "no-interval" = !hasCi
    ? "no-interval"
    : hi <= limit
      ? "meets"
      : lo > limit
        ? "misses"
        : "open";

  // Without an interval, draw neutral ink and say "no interval" (never green on a plain rate).
  const color =
    status === "meets"
      ? "rgb(var(--live-text))"
      : status === "misses"
        ? "rgb(var(--coral-text))"
        : status === "no-interval"
          ? "rgb(var(--ink))"
          : "rgb(var(--muted-text))";
  const statusWords =
    status === "meets"
      ? "interval at or under limit"
      : status === "misses"
        ? "interval over limit"
        : status === "open"
          ? "interval crosses limit"
          : "no interval";

  return (
    <div className="min-w-0" data-gauge-status={status}>
      <div className="flex flex-wrap items-baseline justify-between gap-x-3 gap-y-0.5 text-control">
        <span className="font-medium text-ink">
          {name} <span className="font-normal text-muted-text">({limitText})</span>
        </span>
        <span className="tabular-nums">
          <span className="font-semibold" style={{ color }}>
            {pct(rate)}
          </span>
          {hasCi ? (
            <span className="ml-1.5 text-tiny text-muted-text">
              95% CI {pct(lo)}–{pct(hi)}
            </span>
          ) : null}
        </span>
      </div>

      <div
        role="img"
        aria-label={`${name} ${pct(rate)}${hasCi ? `, 95% interval ${pct(lo)} to ${pct(hi)}` : ""}, against ${limitText}: ${statusWords}.`}
        className="relative mt-1.5 h-4 w-full"
      >
        <span className="absolute inset-x-0 top-1/2 h-1.5 -translate-y-1/2 rounded-full bg-surface-subtle" />
        {/* Reference zone up to limit (neutral tint when no interval is available) */}
        <span
          className="absolute left-0 top-1/2 h-1.5 -translate-y-1/2 rounded-l-full"
          style={{ width: pos(limit), background: hasCi ? "rgb(var(--live) / 0.22)" : "rgb(var(--ink) / 0.08)" }}
        />
        {/* Limit divider */}
        <span
          className="absolute top-0 h-4 w-[1.5px] -translate-x-1/2"
          style={{ left: pos(limit), background: "var(--chart-ref)" }}
        />
        {/* 95% confidence interval whisker */}
        {hasCi ? (
          <span
            className="absolute top-1/2 h-[3px] -translate-y-1/2 rounded-full"
            style={{
              left: pos(lo),
              width: `calc(${pos(hi)} - ${pos(lo)})`,
              background: color,
              opacity: 0.55,
            }}
          />
        ) : null}
        {/* Point estimate mark */}
        <span
          className="absolute top-1/2 size-2.5 -translate-x-1/2 -translate-y-1/2 rounded-full"
          style={{ left: pos(rate), background: color }}
        />
      </div>

      <div className="mt-0.5 flex items-baseline justify-between text-tiny tabular-nums text-muted-text">
        <span>0%</span>
        <span>{statusWords}</span>
        <span>{pct(scaleMax, 0)}</span>
      </div>
    </div>
  );
};

export const ErrorLimitGauge = ({ vme, me, splitLabel }: ErrorLimitGaugeProps) => (
  <div data-amr-primitive="error-limit-gauge" className="grid grid-cols-1 gap-x-8 gap-y-4 border-t border-hair pt-3 sm:grid-cols-2">
    <GaugeRow name={`VME · ${splitLabel}`} limit={BARS.vme.value} limitText="limit ≤1.5%" input={vme} />
    <GaugeRow name={`ME · ${splitLabel}`} limit={BARS.me.value} limitText="limit ≤3%" input={me} />
  </div>
);
