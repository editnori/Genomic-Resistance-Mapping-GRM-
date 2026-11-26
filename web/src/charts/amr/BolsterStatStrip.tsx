// Boxless hairline-separated headline strip adapted from the bolster grammar:
// bare-noun label above, large tabular figure with optional colored delta
// arrow beside it, and one quiet count/whole line beneath.

import type { MouseEvent } from "react";
import { FOCUS_RING } from "../../ui/primitives";

export interface BolsterStat {
  key: string;
  /** Bare-noun label above the number, e.g. "Lineage drop". */
  label: string;
  /** Primary figure, e.g. "87.7%" or "0 of 268". */
  value: string;
  /** Optional compact delta beside the value, e.g. "▼ 3.8 pts". */
  delta?: string;
  /** Valence of the value or delta. */
  tone?: "ink" | "live" | "coral";
  /** Quiet denominator or context line beneath the figure. */
  sub: string;
  /** Optional in-page anchor target to jump directly to the backing section. */
  href?: string;
  onJump?(e: MouseEvent<HTMLAnchorElement>): void;
}

const TONE_CLASS: Record<NonNullable<BolsterStat["tone"]>, string> = {
  ink: "text-ink",
  live: "text-live-text",
  coral: "text-coral-text",
};

export const BolsterStatStrip = ({ stats, ariaLabel = "Key figures" }: { stats: readonly BolsterStat[]; ariaLabel?: string }) => (
  <ul
    data-amr-primitive="bolster-strip"
    aria-label={ariaLabel}
    className="grid grid-cols-1 gap-x-6 gap-y-5 border-y border-hair py-5 sm:grid-cols-2 lg:grid-cols-4"
  >
    {stats.map((s, i) => {
      const tone = TONE_CLASS[s.tone ?? "ink"];
      const content = (
        <>
          <div className="text-tiny font-medium text-muted-text">{s.label}</div>
          <div className="mt-1 flex flex-wrap items-baseline gap-x-2">
            <span className={`text-stat font-semibold tabular-nums tracking-tight ${s.delta === undefined ? tone : "text-ink"}`}>{s.value}</span>
            {s.delta !== undefined ? <span className={`whitespace-nowrap text-control font-semibold tabular-nums ${tone}`}>{s.delta}</span> : null}
          </div>
          <div className="mt-1 text-note tabular-nums text-muted-text">{s.sub}</div>
        </>
      );
      return (
        <li
          key={s.key}
          className={`min-w-0 ${i > 0 ? "lg:border-l lg:border-dashed lg:border-hair lg:pl-6" : ""} ${i % 2 === 1 ? "sm:border-l sm:border-dashed sm:border-hair sm:pl-6" : ""}`}
        >
          {s.href !== undefined ? (
            <a href={s.href} onClick={s.onJump} className={`block rounded-sm transition-colors duration-fast ease-out hover:text-blue ${FOCUS_RING}`}>
              {content}
            </a>
          ) : (
            content
          )}
        </li>
      );
    })}
  </ul>
);
