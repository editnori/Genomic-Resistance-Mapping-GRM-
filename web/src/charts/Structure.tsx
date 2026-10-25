// Population structure: the average-linkage dendrogram over genomes, a band
// marking each lineage, and a genome-by-drug heatmap in the same leaf order,
// so resistance that follows lineages lines up under its branch. Drawn on
// canvas (a species can carry thousands of genomes); drug names are HTML.
// The figure is one tab stop, a grid of genomes by rows: arrow keys move a
// cursor over it and show the tooltip a hover shows.

import { useEffect, useId, useRef, useState, type KeyboardEvent } from "react";
import type { SpeciesReport } from "../data/schema";
import type { TipLine } from "../ui/HoverTip";
import { FOCUS_RING } from "../ui/primitives";
import { capitalize, int, stLabel } from "../data/format";
import { dendrogram } from "../data/derive/layout";
import { ChartFrame, ChartTip, prepareCanvas, stepCursor, tipText, useCanvasTokens, useWidth } from "./kit";

const LABEL_W = 132;
const DENDRO_H = 120;
const BAND_H = 8;
const ROW_H = 14;
const GAP = 6;

const CALL_WORD: Record<string, string> = { R: "resistant", S: "susceptible", I: "intermediate", ".": "not tested" };

/** The tooltip for the genome at leaf position `leaf`: its call for the drug
 *  in `row`, or how many drugs it resists when `row` is null (the dendrogram
 *  and lineage band). */
const genomeTip = (report: SpeciesReport, drugOrder: readonly number[], leaf: number, row: number | null): TipLine[] | null => {
  const g = report.genomes[report.leaf_order[leaf]!];
  if (g === undefined) return null;
  const drug = row === null ? undefined : drugOrder[row];
  const nR = [...g.profile].filter((c) => c === "R").length;
  const nTested = [...g.profile].filter((c) => c !== ".").length;
  return [
    {
      term: g.id,
      text: drug !== undefined ? `${capitalize(report.drugs[drug] ?? "")}: ${CALL_WORD[g.profile[drug] ?? "."]}.` : `Resistant to ${int(nR)} of ${int(nTested)} drugs tested.`,
    },
    { source: `Lineage ${g.lineage}${g.mlst !== null ? ` · ${stLabel(g.mlst)}` : ""}${g.country !== null ? ` · ${g.country}` : ""}${g.year !== null ? ` ${g.year}` : ""}` },
  ];
};

export const StructureFigure = ({ report, drugOrder, selected }: { report: SpeciesReport; drugOrder: readonly number[]; selected: number | null }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const canvas = useRef<HTMLCanvasElement | null>(null);
  const cellId = useId();
  const width = useWidth(frame);
  const tokens = useCanvasTokens();
  // The active cell, from the pointer or the keyboard; `row` null is the
  // dendrogram and lineage band. `tip` false keeps the keyboard's place while
  // Escape or blur hides the tooltip.
  const [hover, setHover] = useState<{ leaf: number; row: number | null; x: number; y: number; tip: boolean } | null>(null);

  const n = report.genomes.length;
  const plotW = Math.max(width - LABEL_W - 4, 100);
  const top = DENDRO_H + GAP + BAND_H + GAP;
  const height = top + drugOrder.length * ROW_H + 4;
  const d = dendrogram(report.merges, report.leaf_order);
  const leafPx = (pos: number) => LABEL_W + ((pos + 0.5) * plotW) / Math.max(n, 1);
  // Square-root height scale so merges among near-identical genomes, which
  // carry the lineage structure, stay visible beside the deep splits.
  const hPx = (h: number) => DENDRO_H - Math.sqrt(h / Math.max(d.maxHeight, 1e-12)) * (DENDRO_H - 4);

  useEffect(() => {
    const el = canvas.current;
    if (el === null || n === 0) return;
    const ctx = prepareCanvas(el, width, height);
    if (ctx === null) return;
    ctx.clearRect(0, 0, width, height);

    ctx.strokeStyle = tokens.axis;
    ctx.lineWidth = 1;
    ctx.beginPath();
    d.links.forEach((l) => {
      const xa = leafPx(l.xa);
      const xb = leafPx(l.xb);
      const y = hPx(l.y);
      ctx.moveTo(xa, hPx(l.ya));
      ctx.lineTo(xa, y);
      ctx.lineTo(xb, y);
      ctx.lineTo(xb, hPx(l.yb));
    });
    ctx.stroke();

    // Lineage band: alternate two neutral steps per contiguous lineage run.
    const bandY = DENDRO_H + GAP;
    let shade = 0;
    report.leaf_order.forEach((g, pos) => {
      const prev = pos === 0 ? undefined : report.genomes[report.leaf_order[pos - 1]!];
      if (prev !== undefined && prev.lineage !== report.genomes[g]!.lineage) shade = 1 - shade;
      ctx.fillStyle = shade === 0 ? tokens.ink : tokens.muted;
      ctx.globalAlpha = shade === 0 ? 0.55 : 0.2;
      const x0 = LABEL_W + (pos * plotW) / n;
      ctx.fillRect(x0, bandY, plotW / n + 0.35, BAND_H);
    });
    ctx.globalAlpha = 1;

    drugOrder.forEach((j, row) => {
      const y = top + row * ROW_H;
      ctx.fillStyle = tokens.fill;
      ctx.fillRect(LABEL_W, y + 1, plotW, ROW_H - 2);
      report.leaf_order.forEach((g, pos) => {
        const call = report.genomes[g]!.profile[j];
        if (call === "." || call === undefined) return;
        ctx.fillStyle = call === "R" ? tokens.r : call === "S" ? tokens.s : tokens.neutral;
        ctx.fillRect(LABEL_W + (pos * plotW) / n, y + 1, plotW / n + 0.35, ROW_H - 2);
      });
    });
    if (selected !== null) {
      const row = drugOrder.indexOf(selected);
      if (row >= 0) {
        ctx.strokeStyle = tokens.ink;
        ctx.lineWidth = 1.5;
        ctx.strokeRect(LABEL_W - 0.75, top + row * ROW_H + 0.25, plotW + 1.5, ROW_H - 0.5);
      }
    }
  }, [report, drugOrder, selected, width, height, tokens]);

  const onMove = (clientX: number, clientY: number) => {
    const r = canvas.current?.getBoundingClientRect();
    if (r === undefined || n === 0) return;
    const x = clientX - r.left;
    const y = clientY - r.top;
    if (x < LABEL_W || x > LABEL_W + plotW) {
      setHover(null);
      return;
    }
    const leaf = Math.min(n - 1, Math.max(0, Math.floor(((x - LABEL_W) / plotW) * n)));
    const rowIdx = Math.floor((y - top) / ROW_H);
    setHover({ leaf, row: y >= top && rowIdx < drugOrder.length ? rowIdx : null, x, y, tip: true });
  };

  // The keyboard starts on the outlined drug's row.
  const selectedRow = selected === null ? -1 : drugOrder.indexOf(selected);
  const cursor = hover ?? { leaf: 0, row: selectedRow >= 0 ? selectedRow : null };
  const lines = genomeTip(report, drugOrder, cursor.leaf, cursor.row);

  // A keyboard move puts the tooltip on the cell: the genome's column, the
  // row's middle.
  const show = (leaf: number, row: number | null) =>
    setHover({ leaf, row, x: leafPx(leaf), y: row === null ? DENDRO_H + GAP + BAND_H / 2 : top + row * ROW_H + ROW_H / 2, tip: true });

  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key === "Escape" && hover?.tip === true) {
      // Escape closes the tooltip and nothing behind it, as HoverTip does.
      e.preventDefault();
      setHover({ ...hover, tip: false });
      return;
    }
    if (lines === null) return;
    if (e.key === "ArrowUp" || e.key === "ArrowDown") {
      // Rows run from the dendrogram and lineage band (-1 here) down the drugs.
      const row = Math.min(drugOrder.length - 1, Math.max(-1, (cursor.row ?? -1) + (e.key === "ArrowDown" ? 1 : -1)));
      e.preventDefault();
      show(cursor.leaf, row < 0 ? null : row);
      return;
    }
    const leaf = stepCursor(e.key, hover === null ? -1 : hover.leaf, n);
    if (leaf === null) return;
    e.preventDefault();
    show(leaf, cursor.row);
  };

  return (
    <ChartFrame frameRef={frame}>
      <div
        role="grid"
        tabIndex={0}
        aria-label={`Dendrogram of ${int(n)} genomes above their resistance calls for ${drugOrder.length} drugs, in the same order.`}
        aria-rowcount={drugOrder.length + 1}
        aria-colcount={n}
        aria-activedescendant={lines !== null ? cellId : undefined}
        onKeyDown={onKey}
        onFocus={(e) => {
          // Keyboard focus opens the tooltip where the cursor rests; a click
          // leaves it to the pointer.
          if (lines !== null && e.currentTarget.matches(":focus-visible")) show(cursor.leaf, cursor.row);
        }}
        onBlur={() => setHover((v) => (v === null ? null : { ...v, tip: false }))}
        className={FOCUS_RING}
      >
        <canvas ref={canvas} aria-hidden="true" onPointerMove={(e) => onMove(e.clientX, e.clientY)} onPointerLeave={() => setHover(null)} className="block" />
        {lines !== null ? (
          <div role="row" aria-rowindex={cursor.row === null ? 1 : cursor.row + 2} className="sr-only">
            <div role="gridcell" id={cellId} aria-colindex={cursor.leaf + 1}>
              {tipText(lines)}
            </div>
          </div>
        ) : null}
      </div>
      <div aria-hidden="true" className="pointer-events-none absolute left-0 top-0" style={{ width: LABEL_W - 8 }}>
        <div className="flex items-end justify-end text-right text-tiny text-muted-text" style={{ height: DENDRO_H }}>
          Mash distance, square-root scale
        </div>
        <div className="text-right text-tiny text-muted-text" style={{ marginTop: GAP - 3, height: BAND_H + GAP }}>
          Lineages
        </div>
        {drugOrder.map((j) => (
          <div
            key={j}
            className={`truncate text-right text-[11.5px] leading-[14px] ${j === selected ? "font-semibold text-ink" : "text-muted-text"}`}
            style={{ height: ROW_H }}
          >
            {capitalize(report.drugs[j] ?? "")}
          </div>
        ))}
      </div>
      {hover?.tip === true && lines !== null ? (
        <>
          <div aria-hidden="true" className="pointer-events-none absolute w-px bg-ink/60" style={{ left: hover.x, top: 0, height }} />
          <ChartTip x={hover.x} y={hover.y} width={width} lines={lines} />
        </>
      ) : null}
    </ChartFrame>
  );
};
