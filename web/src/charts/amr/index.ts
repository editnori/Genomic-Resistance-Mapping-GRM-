// Shared barrel and pure data helpers for AMR-specific visual primitives.

import type { CatalogRow } from "../../data/schema";
import { capitalize, int } from "../../data/format";
import type { DumbbellRow } from "./SplitDumbbell";

export { MicPlateStrip, type MicPlateStripProps } from "./MicPlateStrip";
export { WholePartsBar, type BarPart, type PartTone, type WholePartsBarProps } from "./WholePartsBar";
export { SplitDumbbell, type DumbbellRow, type SplitDumbbellProps } from "./SplitDumbbell";
export { BolsterStatStrip, type BolsterStat } from "./BolsterStatStrip";
export { AmrPeerStrip, type AmrPeerStripProps } from "./AmrPeerStrip";
export { ErrorLimitGauge, type ErrorLimitGaugeProps, type ErrorRateInput } from "./ErrorLimitGauge";

const median = (xs: readonly number[]): number | null => {
  if (xs.length === 0) return null;
  const s = [...xs].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 === 1 ? s[m]! : (s[m - 1]! + s[m]!) / 2;
};

/** Antibiotic-class random vs lineage-held-out medians for classes with at
 *  least `minDatasets` paired datasets, sorted from steepest paired median drop
 *  to smallest (matching `splitSummary` in `data/derive/overview.ts`). */
export const classDropRows = (rows: readonly CatalogRow[], minDatasets = 5): DumbbellRow[] => {
  const paired = rows.filter((r) => r.best_random !== null && r.best_lineage !== null);
  if (paired.length === 0) return [];
  const byClass = new Map<string, CatalogRow[]>();
  for (const r of paired) {
    const k = r.drug_class.trim() || "unclassified";
    const list = byClass.get(k);
    if (list === undefined) byClass.set(k, [r]);
    else list.push(r);
  }
  const out: (DumbbellRow & { sortDelta: number })[] = [];
  for (const [cls, list] of byClass.entries()) {
    if (list.length < minDatasets) continue;
    const before = median(list.map((r) => r.best_random!));
    const after = median(list.map((r) => r.best_lineage!));
    // Median of each dataset's own drop (lineage - random, so a drop is negative).
    const pairedDelta = median(list.map((r) => r.best_lineage! - r.best_random!));
    out.push({
      key: cls,
      label: capitalize(cls),
      before,
      after,
      delta: pairedDelta,
      countText: `${int(list.length)} of ${int(paired.length)} datasets`,
      sortDelta: pairedDelta ?? 0,
    });
  }
  return out.sort((a, b) => a.sortDelta - b.sortDelta || a.label.localeCompare(b.label));
};
