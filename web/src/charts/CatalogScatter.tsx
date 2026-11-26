// The catalog-containment rule against the best learner with lineages held
// out, one mark per dataset that has both. The diagonal is equal agreement
// and the band around it is within CATALOG_CUTOFF; the three regions carry
// their counts. A mark's shape says whether the rule calls nearly every
// isolate one way (uniformCall), which is why so many sit near 50%. From
// 760px the legend stands beside the plot; narrower, it goes below.

import { useLayoutEffect, useRef, useState } from "react";
import type { CatalogRow } from "../data/schema";
import { datasetTitle, int, pct } from "../data/format";
import { CATALOG_CUTOFF, catalogDiff, catalogSummary, uniformCall, type UniformCall } from "../data/derive/overview";
import { AxisX, AxisY, LABEL_CLASS, labelWidth, linear, MarkListbox, RefLineH, ticks, useSvgId, useWidth, type ListMark } from "./kit";

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
  <svg aria-hidden="true" width={16} height={12} className="shrink-0">
    <CatalogMark call={call} x={8} y={6} />
  </svg>
);

const HalfSwatch = () => (
  <svg aria-hidden="true" width={16} height={12} className="shrink-0">
    <line x1={0} x2={16} y1={6} y2={6} stroke="var(--chart-ref)" strokeDasharray="4 3" />
  </svg>
);

/** The plot's width beside a legend this wide, from WIDE_FROM. */
const LEGEND_W = 300;
const WIDE_FROM = 760;

export const CatalogScatter = ({ rows, height, onOpen }: { rows: readonly CatalogRow[]; height?: number; onOpen?(id: string): void }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const avail = useWidth(frame);
  const footerRef = useRef<HTMLDivElement | null>(null);
  const [footerHeight, setFooterHeight] = useState(156);
  const wide = avail >= WIDE_FROM;
  // The footer element changes when the legend moves beside the plot, so the
  // measurement follows the layout.
  useLayoutEffect(() => {
    const el = footerRef.current;
    if (el === null) return;
    const read = () => setFooterHeight(el.offsetHeight);
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, [wide]);
  const summary = catalogSummary(rows);
  const compared = rows.filter((r) => catalogDiff(r) !== null);
  const all = compared.flatMap((r) => [r.catalog!.balanced_accuracy, r.best_lineage!]);
  const lo = Math.max(0, Math.floor(Math.min(0.5, ...all) * 10) / 10);
  // Keep the plot box strictly square (plot width === plot height) at every
  // viewport because both axes share one scale and the 45-degree diagonal is
  // the comparison. Header rows sit above the plot box and never squash it.
  const { ahead, within, behind } = summary.regions;
  const cut = Math.round(CATALOG_CUTOFF * 100);
  const aheadText = `${int(ahead)} ahead by more than ${cut}`;
  const withinText = `${int(within)} within ${cut} points`;
  const behindText = `${int(behind)} behind by more than ${cut}`;
  const padL = 40;
  const padR = 18;
  const padB = 30;
  const wAhead = labelWidth(aheadText);
  const wWithin = labelWidth(withinText);
  const wBehind = labelWidth(behindText);
  const maxSvgW = Math.max(160, wide ? Math.min(560, avail - LEGEND_W - 40) : avail);
  const maxPlotW = Math.max(80, Math.min(maxSvgW - padL - padR, height === undefined ? 464 : Infinity));
  let chosen: { plot: number; width: number; headX: number; oneRow: boolean; twoRow: boolean; padT: number } | null = null;
  for (const rowsCount of [1, 2, 3] as const) {
    const padT = 20 + rowsCount * 16 + 8;
    const maxPlotH = height === undefined ? maxPlotW : Math.max(80, height - footerHeight - padT - padB);
    const plot = Math.max(80, Math.min(maxPlotW, maxPlotH));
    const w = plot + padL + padR;
    const headX = w < 300 ? 2 : padL;
    const headSpan = Math.max(0, w - padR - headX);
    const oneRow = wAhead + wWithin + wBehind + 24 <= headSpan;
    const twoRow = !oneRow && wAhead + wWithin + 10 <= headSpan;
    const needed = oneRow ? 1 : twoRow ? 2 : 3;
    if (needed <= rowsCount && (chosen === null || plot > chosen.plot)) {
      chosen = { plot, width: w, headX, oneRow: rowsCount === 1, twoRow: rowsCount === 2, padT };
    }
  }
  const { plot, width, headX, oneRow: oneRowRegions, twoRow: twoRowRegions, padT } = chosen!;
  const pad = { l: padL, r: padR, t: padT, b: padB };
  const plotH = plot;
  const svgHeight = pad.t + plotH + pad.b;
  const yTitle = labelWidth("Known-gene rule, balanced accuracy (%)") <= width - headX - 2 ? "Known-gene rule, balanced accuracy (%)" : "Known-gene rule (%)";
  const x = linear([lo, 1], [pad.l, pad.l + plot]);
  const y = linear([lo, 1], [pad.t + plotH, pad.t]);
  const clip = useSvgId("catalog-plot");
  const values = ticks(lo, 1, plot < 260 ? 3 : 6);
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
            term: datasetTitle(r),
            text: `Known-gene rule ${pct(c.balanced_accuracy)}, learner ${pct(r.best_lineage)} with lineages held out: ${d >= 0 ? "ahead by" : "behind by"} ${Math.abs(d * 100).toFixed(1)} points.`,
          },
          {
            source:
              call === "resistant"
                ? `The rule calls ${pct(1 - c.specificity)} of susceptible isolates resistant.`
                : call === "susceptible"
                  ? `The rule calls ${pct(1 - c.sensitivity)} of resistant isolates susceptible.`
                  : `Sensitivity ${pct(c.sensitivity)}, specificity ${pct(c.specificity)}.`,
          },
        ],
      };
    })
    .sort((a, b) => a.x - b.x || a.y - b.y);
  const axisTitle = (
    <p className="pt-1 text-center text-tiny text-muted-text" style={{ width }}>
      Best learner with lineages held out, balanced accuracy (%). Higher is better on both axes.
    </p>
  );
  const legend = (
    <ul className={`space-y-1.5 text-control text-ink ${wide ? "w-[300px] shrink-0 pt-5" : "mt-3"}`}>
      <li className="flex items-start gap-2">
        <Swatch call={null} />
        <span>
          <span className="font-medium tabular-nums">{int(summary.uniform.other)}</span> rules with sensitivity and specificity both above 10%
        </span>
      </li>
      <li className="flex items-start gap-2">
        <Swatch call="resistant" />
        <span>
          <span className="font-medium tabular-nums">{int(summary.uniform.resistant)}</span> call ≥90% of susceptible isolates resistant (specificity ≤10%)
        </span>
      </li>
      <li className="flex items-start gap-2">
        <Swatch call="susceptible" />
        <span>
          <span className="font-medium tabular-nums">{int(summary.uniform.susceptible)}</span> call ≥90% of resistant isolates susceptible (sensitivity ≤10%)
        </span>
      </li>
      {lo < 0.5 ? (
        <li className="flex items-start gap-2">
          <HalfSwatch />
          <span>50%: a rule that calls every isolate one way</span>
        </li>
      ) : null}
    </ul>
  );
  return (
    <div ref={frame} className={`flex min-w-0 ${wide ? "flex-row items-start gap-10" : "flex-col items-start"}`}>
      <div className="min-w-0">
        <MarkListbox
          label={`Known-gene rule against the best learner with lineages held out, one mark per dataset, ${int(compared.length)} datasets. Arrow keys move between datasets; click or Enter to pin details.`}
          marks={marks}
          width={width}
          onOpen={onOpen === undefined ? undefined : (i) => onOpen(marks[i]!.key)}
          after={wide ? <div ref={footerRef}>{axisTitle}</div> : undefined}
        >
          <svg aria-hidden="true" width={width} height={svgHeight} className="block overflow-visible">
            <defs>
              <clipPath id={clip}>
                <rect x={pad.l} y={pad.t} width={width - pad.l - pad.r} height={plotH} />
              </clipPath>
            </defs>
            <text x={headX} y={13} className="fill-[rgb(var(--muted-text))] text-[12px]">
              {yTitle}
            </text>
            {oneRowRegions ? (
              <>
                <text x={headX} y={31} className={`${LABEL_CLASS} font-medium`}>
                  {aheadText}
                </text>
                <text x={headX + wAhead + 12} y={31} className={`${LABEL_CLASS} font-medium`}>
                  {withinText}
                </text>
                <text x={width - pad.r} y={31} textAnchor="end" className={`${LABEL_CLASS} font-medium`}>
                  {behindText}
                </text>
              </>
            ) : twoRowRegions ? (
              <>
                <text x={headX} y={31} className={`${LABEL_CLASS} font-medium`}>
                  {aheadText}
                </text>
                <text x={width - pad.r} y={31} textAnchor="end" className={`${LABEL_CLASS} font-medium`}>
                  {withinText}
                </text>
                <text x={headX} y={47} className={`${LABEL_CLASS} font-medium`}>
                  {behindText}
                </text>
              </>
            ) : (
              <>
                <text x={headX} y={31} className={`${LABEL_CLASS} font-medium`}>
                  {aheadText}
                </text>
                <text x={headX} y={47} className={`${LABEL_CLASS} font-medium`}>
                  {withinText}
                </text>
                <text x={headX} y={63} className={`${LABEL_CLASS} font-medium`}>
                  {behindText}
                </text>
              </>
            )}
            <AxisY y={y} x0={pad.l} x1={width - pad.r} values={values} format={(v) => `${Math.round(v * 100)}`} />
            <AxisX x={x} y={pad.t + plotH} height={0} values={values} format={(v) => `${Math.round(v * 100)}`} />
            <path d={band} fill="rgb(var(--ink) / 0.06)" clipPath={`url(#${clip})`} />
            <line x1={x(lo)} y1={y(lo)} x2={x(1)} y2={y(1)} stroke="var(--chart-axis)" clipPath={`url(#${clip})`} />
            {marks.map((m) => (
              <g key={m.key} data-dataset={m.key}><CatalogMark call={m.call} x={m.x} y={m.y} /></g>
            ))}
            {lo < 0.5 ? <RefLineH y={y(0.5)} x0={pad.l} x1={width - pad.r} /> : null}
          </svg>
        </MarkListbox>
      </div>
      {wide ? (
        legend
      ) : (
        <div ref={footerRef} className="w-full">
          {axisTitle}
          {legend}
        </div>
      )}
    </div>
  );
};
