// The catalog-containment rule against the best learner with lineages held
// out, one mark per dataset that has both. The diagonal is equal agreement
// and the band around it is within CATALOG_CUTOFF; the three regions carry
// their counts. A mark's shape says whether the rule calls nearly every
// isolate one way (uniformCall), which is why so many sit near 50%.

import { useLayoutEffect, useRef, useState } from "react";
import type { CatalogRow } from "../data/schema";
import { capitalize, int, pct, speciesShort } from "../data/format";
import { CATALOG_CUTOFF, catalogDiff, catalogSummary, uniformCall, type UniformCall } from "../data/derive/overview";
import { AxisX, AxisY, LABEL_CLASS, linear, MarkListbox, RefLineH, ticks, useSvgId, useWidth, type ListMark } from "./kit";

const R = 4;

/** The mark for one dataset: a filled dot, or hollow when the rule calls
 *  nearly every isolate one way (ring: resistant, square: susceptible). */
export const CatalogMark = ({ call, x, y }: { call: UniformCall | null; x: number; y: number }) =>
  call === "susceptible" ? (
    <rect x={x - 3.5} y={y - 3.5} width={7} height={7} fill="rgb(var(--paper-pure))" stroke="var(--chart-mark)" strokeWidth={1.5} />
  ) : call === "resistant" ? (
    <circle cx={x} cy={y} r={R - 0.5} fill="rgb(var(--paper-pure))" stroke="var(--chart-mark)" strokeWidth={1.5} />
  ) : (
    <circle cx={x} cy={y} r={R} fill="var(--chart-mark)" stroke="rgb(var(--paper-pure))" strokeWidth={1.5} />
  );

const Swatch = ({ call }: { call: UniformCall | null }) => (
  <svg aria-hidden="true" width={12} height={12} className="shrink-0">
    <CatalogMark call={call} x={6} y={6} />
  </svg>
);

export const CatalogScatter = ({ rows, height, onOpen }: { rows: readonly CatalogRow[]; height?: number; onOpen?(id: string): void }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const avail = useWidth(frame);
  const footerRef = useRef<HTMLDivElement | null>(null);
  const [footerHeight, setFooterHeight] = useState(156);
  useLayoutEffect(() => {
    const el = footerRef.current;
    if (el === null) return;
    const read = () => setFooterHeight(el.offsetHeight);
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  const summary = catalogSummary(rows);
  const compared = rows.filter((r) => catalogDiff(r) !== null);
  const all = compared.flatMap((r) => [r.catalog!.balanced_accuracy, r.best_lineage!]);
  const lo = Math.max(0, Math.floor(Math.min(0.5, ...all) * 10) / 10);
  // The axis title and wrapped legend reserve their actual height.
  const side = Math.max(160, Math.min(avail, height === undefined ? 520 : height - footerHeight));
  const pad = { l: 46, r: 10, t: 22, b: 34 };
  const width = side;
  const plotH = side - pad.t - pad.b;
  const x = linear([lo, 1], [pad.l, width - pad.r]);
  const y = linear([lo, 1], [pad.t + plotH, pad.t]);
  const clip = useSvgId("catalog-plot");
  const values = ticks(lo, 1, side < 320 ? 3 : 6);
  const band = `M${x(lo)},${y(lo - CATALOG_CUTOFF)} L${x(1)},${y(1 - CATALOG_CUTOFF)} L${x(1)},${y(1 + CATALOG_CUTOFF)} L${x(lo)},${y(lo + CATALOG_CUTOFF)} Z`;
  const marks: (ListMark & { call: UniformCall | null })[] = compared
    .map((r) => {
      const c = r.catalog!;
      const call = uniformCall(c);
      const d = catalogDiff(r)!;
      return {
        key: r.id,
        x: x(r.best_lineage!),
        y: y(c.balanced_accuracy),
        call,
        tip: () => [
          {
            term: `${capitalize(r.antibiotic)}, ${speciesShort(r.species)}`,
            text: `Catalog rule ${pct(c.balanced_accuracy)}, learner ${pct(r.best_lineage)} with lineages held out: ${d >= 0 ? "ahead by" : "behind by"} ${Math.abs(d * 100).toFixed(1)} points.`,
          },
          {
            source:
              call === "resistant"
                ? `The rule calls ${pct(1 - c.specificity)} of susceptible isolates resistant.`
                : call === "susceptible"
                  ? `The rule calls ${pct(1 - c.sensitivity)} of resistant isolates susceptible.`
                  : `Sensitivity ${pct(c.sensitivity)}, specificity ${pct(c.specificity)}.`,
          },
          ...(onOpen === undefined ? [] : [{ source: "Enter or click opens the dataset." }]),
        ],
      };
    })
    .sort((a, b) => a.x - b.x || a.y - b.y);
  const { ahead, within, behind } = summary.regions;
  return (
    <div ref={frame} className="flex min-w-0 flex-col items-center">
      <MarkListbox
        label={`Catalog rule against the best learner with lineages held out, one mark per dataset, ${int(compared.length)} datasets. Arrow keys move between datasets${onOpen === undefined ? "" : "; Enter opens one"}.`}
        marks={marks}
        width={width}
        onOpen={onOpen === undefined ? undefined : (i) => onOpen(marks[i]!.key)}
      >
        <svg aria-hidden="true" width={width} height={side} className="block overflow-visible">
          <defs>
            <clipPath id={clip}>
              <rect x={pad.l} y={pad.t} width={width - pad.l - pad.r} height={plotH} />
            </clipPath>
          </defs>
          <AxisY y={y} x0={pad.l} x1={width - pad.r} values={values} format={(v) => `${Math.round(v * 100)}`} />
          <AxisX x={x} y={pad.t + plotH} height={0} values={values} format={(v) => `${Math.round(v * 100)}`} />
          <path d={band} fill="rgb(var(--ink) / 0.06)" clipPath={`url(#${clip})`} />
          <line x1={x(lo)} y1={y(lo)} x2={x(1)} y2={y(1)} stroke="var(--chart-axis)" clipPath={`url(#${clip})`} />
          <text x={pad.l + 6} y={pad.t + 14} className={`${LABEL_CLASS} font-medium`}>
            {`${int(ahead)} ahead by more than ${Math.round(CATALOG_CUTOFF * 100)}`}
          </text>
          <text x={width - pad.r} y={pad.t - 8} textAnchor="end" className={`${LABEL_CLASS} font-medium`}>
            {`${int(within)} within ${Math.round(CATALOG_CUTOFF * 100)} points (band)`}
          </text>
          <text x={width - pad.r - 6} y={pad.t + plotH - 10} textAnchor="end" className={`${LABEL_CLASS} font-medium`}>
            {`${int(behind)} behind by more than ${Math.round(CATALOG_CUTOFF * 100)}`}
          </text>
          <text transform={`translate(12 ${pad.t + plotH / 2}) rotate(-90)`} textAnchor="middle" className="fill-[rgb(var(--muted-text))] text-[12px]">
            Catalog rule (%)
          </text>
          {marks.map((m) => (
            <g key={m.key} data-dataset={m.key}><CatalogMark call={m.call} x={m.x} y={m.y} /></g>
          ))}
          {lo < 0.5 ? (
            <>
              <RefLineH y={y(0.5)} x0={pad.l} x1={width - pad.r} />
              <text x={pad.l + 4} y={y(0.5) - 10} textAnchor="start" paintOrder="stroke" stroke="rgb(var(--paper))" strokeWidth={3} className={LABEL_CLASS}>
                50%: calls one way
              </text>
            </>
          ) : null}
        </svg>
      </MarkListbox>
      <div ref={footerRef} className="w-full">
      <p className="mx-auto mt-1 text-center text-tiny text-muted-text" style={{ maxWidth: width }}>
        Best learner with lineages held out, balanced accuracy (%)
      </p>
      <ul className="mt-3 space-y-1 text-control text-ink">
        <li className="flex items-center gap-2">
          <Swatch call={null} />
          <span>
            <span className="font-medium tabular-nums">{int(summary.uniform.other)}</span> rules with sensitivity and specificity both above 10%
          </span>
        </li>
        <li className="flex items-center gap-2">
          <Swatch call="resistant" />
          <span>
            <span className="font-medium tabular-nums">{int(summary.uniform.resistant)}</span> call ≥90% of susceptible isolates resistant (specificity ≤10%)
          </span>
        </li>
        <li className="flex items-center gap-2">
          <Swatch call="susceptible" />
          <span>
            <span className="font-medium tabular-nums">{int(summary.uniform.susceptible)}</span> call ≥90% of resistant isolates susceptible (sensitivity ≤10%)
          </span>
        </li>
      </ul>
      </div>
    </div>
  );
};
