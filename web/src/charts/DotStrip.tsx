// One mark per dataset on one value axis: a beeswarm that keeps every
// dataset visible, with labeled reference lines, region counts that sum to
// the plotted datasets, and direct labels on up to three named extremes (two
// below 640px). On a slide the swarm's offsets shrink to the height it is
// given, so a dense strip never scrolls. GapSwarm and MicStrip adapt the
// manifest rows; MarkListbox gives the keyboard path.

import { useLayoutEffect, useRef, useState } from "react";
import type { CatalogRow } from "../data/schema";
import { capitalize, median, pct, speciesShort } from "../data/format";
import { largestGaps } from "../data/derive/overview";
import { REFERENCE_BARS } from "../data/metrics";
import type { TipLine } from "../ui/HoverTip";
import { AxisX, DirectLabel, labelWidth, linear, MarkListbox, PlacedText, placeLabels, RefLine, ticks, useWidth, type LabelSpec, type ListMark } from "./kit";

export interface StripDatum {
  id: string;
  value: number;
  tip: TipLine[];
}

export interface StripRef {
  key: string;
  value: number;
  label: string;
  anchor: "start" | "end";
  /** A baseline such as "no change" is solid; a stated bar is dashed. */
  solid?: boolean;
}

export interface StripRegion {
  key: string;
  label: string;
  /** The value the label sits against, with its side. */
  at: number;
  anchor: "start" | "end";
}

export interface StripNamed {
  id: string;
  lines: [string, string];
}

const R = 4;
const DIAMETER = 9;
const LABEL_ROWS = 3;
const TOP = LABEL_ROWS * 16 + 6;
const AXIS = 22;

/** Greedy swarm: each value takes the slot nearest the midline that does not
 *  overlap an earlier mark within one diameter. Offsets in pixels. */
const swarm = (xs: readonly number[], diameter: number): number[] => {
  const order = xs.map((x, i) => ({ x, i })).sort((a, b) => a.x - b.x);
  const placed: { x: number; y: number }[] = [];
  const ys = new Array<number>(xs.length).fill(0);
  order.forEach(({ x, i }) => {
    for (let k = 0; ; k += 1) {
      const y = (k % 2 === 0 ? 1 : -1) * Math.ceil(k / 2) * diameter;
      if (placed.every((p) => Math.abs(p.x - x) >= diameter || Math.abs(p.y - y) >= diameter)) {
        placed.push({ x, y });
        ys[i] = y;
        break;
      }
    }
  });
  return ys;
};

export const DotStrip = ({
  label,
  data,
  domain,
  format,
  axisTitle,
  refs,
  regions,
  named,
  height,
  onOpen,
}: {
  label: string;
  data: readonly StripDatum[];
  domain: [number, number];
  format(v: number): string;
  axisTitle: string;
  refs: readonly StripRef[];
  regions: readonly StripRegion[];
  named: readonly StripNamed[];
  /** The height to fill; without it the strip takes the height its swarm
   *  needs, up to 220px of band. */
  height?: number;
  onOpen?(id: string): void;
}) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const axisRef = useRef<HTMLParagraphElement | null>(null);
  const [axisHeight, setAxisHeight] = useState(36);
  useLayoutEffect(() => {
    const el = axisRef.current;
    if (el === null) return;
    const read = () => setAxisHeight(el.offsetHeight + 4);
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  const x = linear(domain, [20, width - 20]);
  const xs = data.map((d) => x(d.value));
  const raw = swarm(xs, DIAMETER);
  const maxOffset = Math.max(0, ...raw.map(Math.abs));
  // The axis title is HTML below the plot, 20px.
  const bandHalf = height === undefined ? Math.min(110, maxOffset + R + 4) : Math.max(28, (height - TOP - AXIS - axisHeight) / 2);
  const scale = maxOffset === 0 ? 1 : Math.min(1, (bandHalf - R - 2) / maxOffset);
  const mid = TOP + bandHalf;
  const base = TOP + bandHalf * 2;
  const svgHeight = base + AXIS;
  const marks: ListMark[] = data.map((d, i) => ({ key: d.id, x: xs[i]!, y: mid + raw[i]! * scale, tip: () => d.tip }));
  const ordered = [...marks].sort((a, b) => a.x - b.x || a.y - b.y);

  const specs: LabelSpec[] = [
    ...regions.map((r) => ({ key: r.key, text: r.label, x: x(r.at), anchor: r.anchor })),
    ...refs.map((r) => ({ key: r.key, text: r.label, x: x(r.value), anchor: r.anchor })),
  ];
  const placed = placeLabels(specs, width, LABEL_ROWS);

  // Named marks alternate above and below the midline, nearest first; one
  // that would leave the band or the frame is left to its tooltip.
  const labels = named.slice(0, width < 640 ? 2 : 3).flatMap((n, k) => {
    const m = marks.find((mm) => mm.key === n.id);
    if (m === undefined) return [];
    const above = k % 2 === 0;
    const depth = Math.floor(k / 2) * 34;
    const ly = above ? mid - 14 - depth : mid + 24 + depth;
    if ((above ? mid - ly : ly - mid) + 18 > bandHalf) return [];
    const anchor: "start" | "end" = m.x > width / 2 ? "end" : "start";
    const lx = anchor === "end" ? m.x - 8 : m.x + 8;
    const w = Math.max(...n.lines.map(labelWidth));
    if (anchor === "end" ? lx - w < 0 : lx + w > width) return [];
    return [{ key: n.id, x: m.x, y: m.y, lx, ly: above ? ly - 15 : ly, lines: n.lines, anchor }];
  });

  return (
    <div ref={frame} className="min-w-0">
      <MarkListbox label={label} marks={ordered} width={width} onOpen={onOpen === undefined ? undefined : (i) => onOpen(ordered[i]!.key)}>
        <svg aria-hidden="true" width={width} height={svgHeight} className="block overflow-visible">
          <PlacedText labels={placed} top={0} />
          {refs.map((r) => (
            <RefLine key={r.key} x={x(r.value)} y0={TOP - 4} y1={base} dashed={r.solid !== true} />
          ))}
          <AxisX x={x} y={base} height={0} values={ticks(domain[0], domain[1], width < 480 ? 4 : 8)} format={format} />
          {marks.map((m) => (
            <circle key={m.key} data-dataset={m.key} cx={m.x} cy={m.y} r={R} fill="var(--chart-mark)" stroke="rgb(var(--paper-pure))" strokeWidth={1.5} />
          ))}
          {labels.map((l) => (
            <DirectLabel key={l.key} x={l.x} y={l.y} lx={l.lx} ly={l.ly} lines={l.lines} anchor={l.anchor} />
          ))}
        </svg>
      </MarkListbox>
      <p ref={axisRef} className="mt-1 text-center text-tiny text-muted-text">{axisTitle}</p>
    </div>
  );
};

const datasetName = (r: CatalogRow): string => `${capitalize(r.antibiotic)}, ${speciesShort(r.species)}`;

/** Per paired dataset, random-split minus lineage-held-out balanced
 *  accuracy of the best learner, in points. */
export const GapSwarm = ({ rows, height, onOpen }: { rows: readonly CatalogRow[]; height?: number; onOpen?(id: string): void }) => {
  const paired = rows.filter((r) => r.best_random !== null && r.best_lineage !== null);
  const gaps = paired.map((r) => r.best_random! - r.best_lineage!);
  const mid = median(gaps);
  const above = gaps.filter((g) => g > 0.05).length;
  const lo = Math.floor(Math.min(-0.02, ...gaps) * 100) / 100;
  const hi = Math.ceil(Math.max(0.1, ...gaps) * 100) / 100;
  return (
    <DotStrip
      label={`Balanced accuracy lost when lineages are held out, one mark per dataset, ${paired.length} datasets. Arrow keys move between datasets${onOpen === undefined ? "" : "; Enter opens one"}.`}
      data={paired.map((r) => ({
        id: r.id,
        value: r.best_random! - r.best_lineage!,
        tip: [
          { term: datasetName(r), text: `${((r.best_random! - r.best_lineage!) * 100).toFixed(1)} points: ${pct(r.best_random)} random, ${pct(r.best_lineage)} lineages held out.` },
          ...(onOpen === undefined ? [] : [{ source: "Enter or click opens the dataset." }]),
        ],
      }))}
      domain={[lo, hi]}
      format={(v) => (v * 100).toFixed(0)}
      axisTitle="Points of balanced accuracy lost with lineages held out (random minus lineage)"
      refs={[
        ...(mid === null ? [] : [{ key: "median", value: mid, label: `Median ${(mid * 100).toFixed(1)}`, anchor: "start" as const }]),
        { key: "bar", value: 0.05, label: "5 points", anchor: "end" },
        { key: "zero", value: 0, label: "No change", anchor: "end", solid: true },
      ]}
      regions={[
        { key: "below", label: `${gaps.length - above} at or below 5`, at: lo, anchor: "start" },
        { key: "above", label: `${above} above 5`, at: hi, anchor: "end" },
      ]}
      named={largestGaps(paired, 3).map(({ row }) => ({ id: row.id, lines: [datasetName(row), `${pct(row.best_random)} to ${pct(row.best_lineage)}`] }))}
      height={height}
      onOpen={onOpen}
    />
  );
};

/** Per dataset with an MIC model, the share of exact MICs predicted within
 *  one doubling dilution. Every dataset is drawn. */
export const MicStrip = ({ rows, height, onOpen }: { rows: readonly CatalogRow[]; height?: number; onOpen?(id: string): void }) => {
  const scored = rows.filter((r) => r.mic_essential_agreement !== null);
  const values = scored.map((r) => r.mic_essential_agreement!);
  const mid = median(values);
  const atOrAbove = values.filter((v) => v >= REFERENCE_BARS.mic).length;
  const lo = Math.min(0.8, Math.floor(Math.min(...values) * 10) / 10);
  const lowest = [...scored].sort((a, b) => a.mic_essential_agreement! - b.mic_essential_agreement! || a.id.localeCompare(b.id))[0];
  return (
    <DotStrip
      label={`MIC essential agreement, one mark per dataset, ${scored.length} datasets. Arrow keys move between datasets${onOpen === undefined ? "" : "; Enter opens one"}.`}
      data={scored.map((r) => ({
        id: r.id,
        value: r.mic_essential_agreement!,
        tip: [
          { term: datasetName(r), text: `${pct(r.mic_essential_agreement)} of exact MICs predicted within one doubling dilution.` },
          ...(onOpen === undefined ? [] : [{ source: "Enter or click opens the dataset." }]),
        ],
      }))}
      domain={[lo, 1]}
      format={(v) => `${Math.round(v * 100)}`}
      axisTitle="Essential agreement: exact MICs predicted within one doubling dilution (%)"
      refs={[
        { key: "bar", value: REFERENCE_BARS.mic, label: "90% bar", anchor: "start" },
        ...(mid === null ? [] : [{ key: "median", value: mid, label: `Median ${pct(mid)}`, anchor: "end" as const }]),
      ]}
      regions={[
        { key: "below", label: `${values.length - atOrAbove} below 90%`, at: lo, anchor: "start" },
        { key: "above", label: `${atOrAbove} at or above`, at: 1, anchor: "end" },
      ]}
      named={lowest === undefined ? [] : [{ id: lowest.id, lines: [datasetName(lowest), pct(lowest.mic_essential_agreement)] }]}
      height={height}
      onOpen={onOpen}
    />
  );
};
