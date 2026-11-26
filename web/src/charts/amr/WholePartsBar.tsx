// Whole-and-parts segmented bar: partitions a single denominator into direct
// named segments, each showing its count and percent beneath the bar.
// Supports diagonal hatch fill for a bucket of a different kind (for example,
// undecided intervals or MIC-inferred labels excluded from model scoring).

import { int, pct } from "../../data/format";

export type PartTone = "live" | "coral" | "blue" | "neutral" | "resistant" | "susceptible";

export interface BarPart {
  key: string;
  label: string;
  count: number;
  tone: PartTone;
  /** Diagonal hatch marks a distinct kind of bucket rather than a verdict. */
  hatched?: boolean;
  sub?: string;
}

export interface WholePartsBarProps {
  /** Heading above the bar naming the whole, e.g. "Datasets tested". */
  wholeLabel: string;
  /** Total denominator; if omitted, uses the sum of part counts. */
  total?: number;
  parts: readonly BarPart[];
  /** Optional accessible label override. */
  ariaLabel?: string;
  /** Compact height when nested inside a table or card. */
  compact?: boolean;
}

const SOLID_FILL: Record<PartTone, string> = {
  live: "rgb(var(--live))",
  coral: "rgb(var(--coral))",
  blue: "rgb(var(--blue))",
  neutral: "var(--series-neutral)",
  resistant: "var(--series-r)",
  susceptible: "var(--series-s)",
};

const TEXT_CLASS: Record<PartTone, string> = {
  live: "text-live-text",
  coral: "text-coral-text",
  blue: "text-blue",
  neutral: "text-muted-text",
  resistant: "text-coral-text",
  susceptible: "text-blue",
};

const hatchFill = (tone: PartTone): string =>
  `repeating-linear-gradient(135deg, ${SOLID_FILL[tone]} 0 3px, rgb(var(--paper-pure)) 3px 6px)`;

export const WholePartsBar = ({ wholeLabel, total, parts, ariaLabel, compact = false }: WholePartsBarProps) => {
  const denom = total ?? parts.reduce((sum, p) => sum + p.count, 0);
  const label =
    ariaLabel ??
    `${wholeLabel}: ${int(denom)} total. ${parts.map((p) => `${p.label}: ${int(p.count)} (${pct(denom === 0 ? 0 : p.count / denom)})`).join("; ")}.`;

  return (
    <div data-amr-primitive="whole-parts-bar" className="w-full max-w-[46rem]">
      <div className="flex flex-wrap items-baseline gap-x-2 gap-y-1">
        <span className="text-control font-semibold text-ink">{wholeLabel}</span>
        <span data-whole-count={denom} className="text-body font-semibold tabular-nums text-ink">{int(denom)}</span>
      </div>

      <div
        role="img"
        aria-label={label}
        className={`mt-2 flex w-full gap-[2px] overflow-hidden rounded-sm bg-surface-subtle ${compact ? "h-4" : "h-7"}`}
      >
        {parts.map((p) => {
          if (p.count <= 0 || denom <= 0) return null;
          const share = (p.count / denom) * 100;
          return (
            <div
              key={p.key}
              className="h-full first:rounded-l-sm last:rounded-r-sm"
              style={{
                width: `${Math.max(share, 0.75).toFixed(2)}%`,
                ...(p.hatched === true
                  ? { backgroundImage: hatchFill(p.tone), border: `1px solid ${SOLID_FILL[p.tone]}` }
                  : { background: SOLID_FILL[p.tone] }),
              }}
              title={`${p.label}: ${int(p.count)} of ${int(denom)} (${pct(p.count / denom)})`}
            />
          );
        })}
      </div>

      <div className={`mt-2.5 flex flex-wrap gap-x-6 gap-y-2 ${compact ? "text-tiny" : "text-control"}`}>
        {parts.map((p) => (
          <div key={p.key} className="min-w-0">
            <div className="flex items-center gap-1.5">
              <span
                aria-hidden="true"
                className="inline-block size-2.5 shrink-0 rounded-[2px]"
                style={
                  p.hatched === true
                    ? { backgroundImage: hatchFill(p.tone), border: `1px solid ${SOLID_FILL[p.tone]}` }
                    : { background: SOLID_FILL[p.tone] }
                }
              />
              <span className="font-medium text-ink">{p.label}</span>
            </div>
            <div className="mt-0.5 flex items-baseline gap-1.5 tabular-nums">
              <span data-part-count={p.count} className={`font-semibold ${TEXT_CLASS[p.tone]}`}>{int(p.count)}</span>
              <span className="text-muted-text">{pct(denom === 0 ? 0 : p.count / denom)}</span>
            </div>
            {p.sub !== undefined ? <div className="text-tiny text-muted-text">{p.sub}</div> : null}
          </div>
        ))}
      </div>
    </div>
  );
};
