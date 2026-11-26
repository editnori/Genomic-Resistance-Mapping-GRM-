// Before-vs-after dumbbell comparison rows on an honest 0%–100% track.
// Used for random-split vs lineage-held-out drops (by antibiotic class or by
// learner) so the move is visible at a glance without toggling.

import { useState } from "react";
import { pct } from "../../data/format";
import { FOCUS_RING } from "../../ui/primitives";

export interface DumbbellRow {
  key: string;
  label: string;
  /** Baseline value (0..1), e.g. random-split balanced accuracy. */
  before: number | null;
  /** Comparison value (0..1), e.g. lineage-held-out balanced accuracy. */
  after: number | null;
  /** Optional denominator text, e.g. "20 of 268 datasets". */
  countText?: string;
  /** Optional explicit paired delta (after - before in 0..1 units). Defaults
   *  to (after - before) when omitted. */
  delta?: number | null;
}

export interface SplitDumbbellProps {
  rows: readonly DumbbellRow[];
  beforeLabel?: string;
  afterLabel?: string;
  /** Drop magnitude (in 0..1 units) at or above which the bar and delta
   *  highlight in coral. Default 0.05 (5 points). */
  largeDrop?: number;
  /** Optional cap on initial rows shown, with a clean expand toggle for long lists. */
  initialVisible?: number;
  /** Noun used in the expand toggle, e.g. "drug classes". */
  itemNoun?: string;
  ariaLabel: string;
}

export const SplitDumbbell = ({
  rows,
  beforeLabel = "Random split",
  afterLabel = "Lineages held out",
  largeDrop = 0.05,
  initialVisible,
  itemNoun = "rows",
  ariaLabel,
}: SplitDumbbellProps) => {
  const [expanded, setExpanded] = useState(false);
  const clamp = (v: number) => Math.max(0, Math.min(1, v));
  const collapsible = initialVisible !== undefined && rows.length > initialVisible;
  const shown = collapsible && !expanded ? rows.slice(0, initialVisible) : rows;

  return (
    <div data-amr-primitive="split-dumbbell" role="region" aria-label={ariaLabel} className="w-full max-w-[54rem]">
      {/* Direct dot key aligned over track, and column headers over numbers */}
      <div className="mb-2 grid grid-cols-1 items-baseline gap-x-4 gap-y-1.5 text-tiny text-muted-text sm:grid-cols-[minmax(9rem,13.5rem)_minmax(8rem,1fr)_8.5rem_5.5rem]">
        <span className="hidden sm:block" />
        <div className="flex flex-wrap items-center gap-x-4 gap-y-1">
          <span className="inline-flex items-center gap-1.5">
            <span aria-hidden="true" className="inline-block size-2.5 rounded-full" style={{ background: "var(--series-s)", opacity: 0.45 }} />
            <span>{beforeLabel}</span>
          </span>
          <span className="inline-flex items-center gap-1.5">
            <span aria-hidden="true" className="inline-block size-2.5 rounded-full bg-ink" />
            <span className="font-medium text-ink">{afterLabel}</span>
          </span>
        </div>
        <span className="hidden text-right font-medium text-ink sm:block">Score</span>
        <span className="hidden text-right font-medium text-ink sm:block">Change</span>
      </div>

      <ul className="divide-y divide-hair border-y border-hair">
        {shown.map((r) => {
          const valid = r.before !== null && r.after !== null;
          const diff = r.delta !== undefined ? r.delta : valid ? r.after! - r.before! : null;
          const pts = diff === null ? null : Math.abs(diff * 100);
          const down = diff !== null && diff < -0.0005;
          const up = diff !== null && diff > 0.0005;
          const big = diff !== null && diff <= -largeDrop;
          const lo = valid ? Math.min(clamp(r.before!), clamp(r.after!)) * 100 : 0;
          const hi = valid ? Math.max(clamp(r.before!), clamp(r.after!)) * 100 : 0;

          return (
            <li key={r.key} className="grid grid-cols-1 items-center gap-x-4 gap-y-1.5 py-2.5 sm:grid-cols-[minmax(9rem,13.5rem)_minmax(8rem,1fr)_8.5rem_5.5rem]">
              <div className="min-w-0">
                <div className="truncate text-control font-medium text-ink" title={r.label}>
                  {r.label}
                </div>
                {r.countText !== undefined ? <div className="text-tiny tabular-nums text-muted-text">{r.countText}</div> : null}
              </div>

              {/* 0%–100% honest track */}
              <div className="relative h-5 w-full" aria-hidden="true">
                <span className="absolute inset-x-0 top-1/2 h-1.5 -translate-y-1/2 rounded-full bg-surface-subtle" />
                <span className="absolute left-1/2 top-1/2 h-3 w-px -translate-y-1/2 bg-hair" />
                {valid ? (
                  <>
                    <span
                      className="absolute top-1/2 h-1.5 -translate-y-1/2 rounded-full"
                      style={{
                        left: `${lo.toFixed(2)}%`,
                        width: `${Math.max(hi - lo, 0.4).toFixed(2)}%`,
                        background: big ? "rgb(var(--coral))" : "var(--chart-axis)",
                      }}
                    />
                    <span
                      className="absolute top-1/2 size-3 -translate-x-1/2 -translate-y-1/2 rounded-full"
                      style={{ left: `${(clamp(r.before!) * 100).toFixed(2)}%`, background: "var(--series-s)", opacity: 0.5 }}
                    />
                    <span
                      className="absolute top-1/2 size-3 -translate-x-1/2 -translate-y-1/2 rounded-full bg-ink"
                      style={{ left: `${(clamp(r.after!) * 100).toFixed(2)}%` }}
                    />
                  </>
                ) : null}
              </div>

              {/* Before -> After values */}
              <div className="flex items-baseline justify-between gap-2 text-control tabular-nums sm:justify-end">
                <span className="text-muted-text sm:hidden">Score</span>
                <span className="whitespace-nowrap">
                  <span className="text-muted-text">{pct(r.before)}</span>
                  <span aria-hidden="true" className="mx-1 text-muted-text">
                    →
                  </span>
                  <span className="font-semibold text-ink">{pct(r.after)}</span>
                </span>
              </div>

              {/* Arrow + unsigned magnitude (arrow owns direction) */}
              <div className="flex items-baseline justify-between gap-2 text-control tabular-nums sm:justify-end">
                <span className="text-muted-text sm:hidden">Change</span>
                <span
                  className={`whitespace-nowrap font-semibold ${
                    big ? "text-coral-text" : up ? "text-live-text" : "text-muted-text"
                  }`}
                >
                  {pts === null ? "–" : `${down ? "▼ " : up ? "▲ " : ""}${pts.toFixed(1)} pts`}
                </span>
              </div>
            </li>
          );
        })}
      </ul>

      <div className="mt-1 hidden grid-cols-[minmax(9rem,13.5rem)_minmax(8rem,1fr)_8.5rem_5.5rem] gap-x-4 text-tiny tabular-nums text-muted-text sm:grid">
        <span />
        <div className="flex justify-between">
          <span>0%</span>
          <span>50%</span>
          <span>100%</span>
        </div>
        <span />
        <span />
      </div>

      {collapsible ? (
        <div className="mt-2">
          <button
            type="button"
            aria-expanded={expanded}
            onClick={() => setExpanded((v) => !v)}
            className={`ui-control-target ui-target-24 inline-flex items-center rounded-sm py-1 text-control font-medium text-blue transition-colors duration-fast ease-out hover:text-ink ${FOCUS_RING}`}
          >
            {expanded ? `Show top ${initialVisible} ${itemNoun}` : `Show all ${rows.length} ${itemNoun}`}
          </button>
        </div>
      ) : null}
    </div>
  );
};
