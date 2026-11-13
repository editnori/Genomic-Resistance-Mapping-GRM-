// The outside-study test, one mark per dataset: the report's figure with
// whole lineages held out (across) against the same final model's figure on
// newly deposited genomes (up), with the new figure's 95% interval as a
// vertical bar. The diagonal is equality. A hollow mark reads below: its
// whole interval is under the report's figure. Each population has its own
// panel of two plots, sensitivity and specificity, and the panels are never
// combined. Four plots stand in one row from 880px; narrower, each
// population takes a row of its own. On a slide too short for its plots to
// reach MIN_SIDE, each population gets a slide of its own (the deck's paging,
// as for the panel figure); in the one-row layout that page sets the legend
// beside the two plots, so they take its height. A plot never grows past the
// height it is given, down to MIN_FLOOR.

import { useEffect, useLayoutEffect, useRef, useState, type ReactNode } from "react";
import { int, pct } from "../data/format";
import type { Measure, OutsidePoint } from "../data/derive/outside";
import { AxisX, AxisY, linear, MarkListbox, ticks, useWidth, type ListMark } from "./kit";

export interface OutsidePanel {
  key: string;
  /** The population in words, with its count: the panel's heading. */
  heading: string;
  /** Its name inside a plot's label for a screen reader. */
  name: string;
  points: Record<Measure, readonly OutsidePoint[]>;
}

const MEASURES: readonly Measure[] = ["sensitivity", "specificity"];
const MEASURE_NAME: Record<Measure, string> = { sensitivity: "Sensitivity", specificity: "Specificity" };
const PAD = { l: 32, r: 14, t: 8, b: 22 };
/** A panel's heading, and a plot's measure line and its counts line, above
 *  each row of plots, until the laid-out height replaces it. */
const HEAD_H = 66;
const GAP = 24;
const GROUP_GAP = 40;
const FOUR_FROM = 880;
/** The smallest plot a slide draws both populations at; under it, they page. */
const MIN_SIDE = 150;
/** The smallest plot drawn at all, on a slide with no more room. */
const MIN_FLOOR = 72;
/** The legend's width when it stands beside one population's plots. */
const LEGEND_ASIDE_W = 340;

const readingWords = (p: OutsidePoint): string =>
  p.reading === "below"
    ? "the whole interval is under the report's figure, so it reads below"
    : p.reading === "above"
      ? "the whole interval is over the report's figure, so it reads above"
      : "the interval includes the report's figure";

const Mark = ({ below, x, y }: { below: boolean; x: number; y: number }) =>
  below ? (
    <circle cx={x} cy={y} r={3.5} fill="rgb(var(--paper-pure))" stroke="var(--chart-mark)" strokeWidth={1.5} />
  ) : (
    <circle cx={x} cy={y} r={4} fill="var(--chart-mark)" stroke="rgb(var(--paper-pure))" strokeWidth={1} />
  );

const Plot = ({
  panel,
  measure,
  side,
  titleOf,
  onOpen,
  canOpen,
  measuring,
  oneLine = false,
}: {
  panel: OutsidePanel;
  measure: Measure;
  side: number;
  titleOf(id: string): string;
  onOpen?(id: string): void;
  canOpen?(id: string): boolean;
  /** A deck's off-screen measuring copy: the plot alone, with no keyboard layer. */
  measuring: boolean;
  /** Keep the measure line and the counts line on one line each, wider than a small plot if need be. */
  oneLine?: boolean;
}) => {
  const points = panel.points[measure];
  const plotW = side - PAD.l - PAD.r;
  const plotH = side - PAD.t - PAD.b;
  const x = linear([0, 1], [PAD.l, PAD.l + plotW]);
  const y = linear([0, 1], [PAD.t + plotH, PAD.t]);
  const values = ticks(0, 1, side < 220 ? 2 : 4);
  const above = points.filter((p) => p.reading === "above").length;
  const below = points.filter((p) => p.reading === "below").length;
  const marks: ListMark[] = points
    .map((p) => ({
      key: p.dataset,
      x: x(p.report),
      y: y(p.figure.estimate),
      tip: () => [
        {
          term: titleOf(p.dataset),
          text: `${MEASURE_NAME[measure]}: report ${pct(p.report)} with lineages held out; new genomes ${pct(p.figure.estimate)}, 95% interval ${pct(p.figure.lo)} to ${pct(p.figure.hi)}. ${readingWords(p).replace(/^./, (c) => c.toUpperCase())}.`,
        },
      ],
    }))
    .sort((a, b) => a.x - b.x || a.y - b.y);
  const byKey = new Map(points.map((p) => [p.dataset, p]));
  const keyed = (plot: ReactNode) =>
    measuring ? (
      plot
    ) : (
      <MarkListbox
        label={`${MEASURE_NAME[measure]}, ${panel.name}: the report's figure against the figure on the new genomes, ${int(points.length)} datasets. Arrow keys move between datasets; click or Enter to pin details.`}
        marks={marks}
        width={side}
        onOpen={onOpen === undefined ? undefined : (i) => onOpen(marks[i]!.key)}
        canOpen={canOpen === undefined ? undefined : (i) => canOpen(marks[i]!.key)}
      >
        {plot}
      </MarkListbox>
    );
  return (
    <div className={oneLine ? "whitespace-nowrap" : "min-w-0"} style={oneLine ? { minWidth: side } : { width: side }} data-outside-plot={measure}>
      <p className="text-control font-medium text-ink">{MEASURE_NAME[measure]}</p>
      {/* The counts sit above the plot, so no mark or bar is drawn over them. */}
      <p data-outside-counts className="text-control text-muted-text">{`${int(above)} read above, ${int(below)} read below`}</p>
      {keyed(
        <svg aria-hidden="true" width={side} height={side} className="block overflow-visible">
          <AxisY y={y} x0={PAD.l} x1={PAD.l + plotW} values={values} format={(v) => `${Math.round(v * 100)}`} />
          <AxisX x={x} y={PAD.t + plotH} height={0} values={values} format={(v) => `${Math.round(v * 100)}`} />
          <line x1={x(0)} y1={y(0)} x2={x(1)} y2={y(1)} stroke="var(--chart-ref)" strokeDasharray="4 3" />
          {marks.map((m) => {
            const p = byKey.get(m.key)!;
            return (
              <g key={m.key} data-outside={measuring ? undefined : m.key} data-reading={p.reading}>
                <line x1={m.x} x2={m.x} y1={y(p.figure.lo)} y2={y(p.figure.hi)} stroke="var(--chart-mark)" strokeOpacity={0.45} strokeWidth={1.5} />
                <Mark below={p.reading === "below"} x={m.x} y={m.y} />
              </g>
            );
          })}
        </svg>,
      )}
    </div>
  );
};

export const OutsideScatter = ({
  panels,
  height,
  titleOf,
  onOpen,
  canOpen,
  page,
  pages,
  onPages,
}: {
  panels: readonly OutsidePanel[];
  height?: number;
  titleOf(id: string): string;
  onOpen?(id: string): void;
  canOpen?(id: string): boolean;
  page?: number;
  pages?: number;
  onPages?(pages: number): void;
}) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const legendRef = useRef<HTMLDivElement | null>(null);
  const [legendH, setLegendH] = useState(48);
  useLayoutEffect(() => {
    const el = legendRef.current;
    if (el === null) return undefined;
    const read = () => setLegendH(el.offsetHeight);
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  // The laid-out height from a panel's top to its first plot, which replaces
  // HEAD_H once the frame is drawn: text lines are the font's height, not ours.
  const [headH, setHeadH] = useState(HEAD_H);
  useLayoutEffect(() => {
    const panel = frame.current?.querySelector<HTMLElement>("[data-outside-panel]");
    const svg = panel?.querySelector("svg");
    if (panel === undefined || panel === null || svg === undefined || svg === null) return;
    const h = Math.ceil(svg.getBoundingClientRect().top - panel.getBoundingClientRect().top);
    if (h > 0) setHeadH((was) => (was === h ? was : h));
  });
  const four = width >= FOUR_FROM;
  const byWidth = Math.floor(four ? (width - 2 * GAP - GROUP_GAP) / 4 : (width - GAP) / 2);
  const byHeight = (rows: number, legend = legendH) => (height === undefined ? Infinity : Math.floor((height - legend - rows * headH - (rows - 1) * 16) / rows));
  const needed = Math.min(byWidth, byHeight(four ? 1 : panels.length)) >= MIN_SIDE ? 1 : panels.length;
  // An unlaid box measures nothing, so it reports nothing.
  useEffect(() => {
    if (height !== undefined && height > 0) onPages?.(needed);
  }, [onPages, needed, height]);
  const count = Math.max(1, Math.min(pages ?? 1, panels.length));
  const at = Math.min(Math.max(0, page ?? 0), count - 1);
  const shown = count === 1 ? panels : panels.slice(at, at + 1);
  const rows = four ? 1 : shown.length;
  // One population on a wide slide: the legend stands beside its plots.
  const aside = four && count > 1;
  const side = aside
    ? Math.max(MIN_FLOOR, Math.min(300, Math.floor((width - GAP - GROUP_GAP - LEGEND_ASIDE_W) / 2), byHeight(1, 0)))
    : Math.max(MIN_FLOOR, Math.min(300, byWidth, byHeight(rows)));
  const measuring = onPages !== undefined;
  return (
    <div ref={frame} className={aside ? "flex min-w-0 flex-row items-start gap-x-10" : "min-w-0"}>
      <div className={`flex ${four ? "flex-row gap-x-10" : "flex-col gap-y-4"}`}>
        {shown.map((panel) => (
          <div key={panel.key} className="min-w-0" data-outside-panel>
            <h3 className="text-control text-muted-text">{panel.heading}</h3>
            <div className="mt-1 flex gap-x-6">
              {MEASURES.map((m) => (
                <Plot key={m} panel={panel} measure={m} side={side} titleOf={titleOf} onOpen={onOpen} canOpen={canOpen} measuring={measuring} oneLine={aside} />
              ))}
            </div>
          </div>
        ))}
      </div>
      <div
        ref={legendRef}
        className={aside ? "flex flex-col gap-y-1 pt-1 text-control text-ink" : "flex flex-wrap gap-x-5 gap-y-1 pt-2 text-control text-ink"}
        style={aside ? { width: LEGEND_ASIDE_W } : undefined}
      >
        <span className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={12} height={12}>
            <circle cx={6} cy={6} r={4} fill="var(--chart-mark)" />
          </svg>
          Not below on this measure
        </span>
        <span className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={12} height={12}>
            <circle cx={6} cy={6} r={3.5} fill="none" stroke="var(--chart-mark)" strokeWidth={1.5} />
          </svg>
          Reads below: the whole interval is under the line
        </span>
        <span>Reads above: the whole interval is over the line</span>
        <span className="text-muted-text">Across: the report&apos;s figure. Up: the new figure. Both in %.</span>
      </div>
    </div>
  );
};
