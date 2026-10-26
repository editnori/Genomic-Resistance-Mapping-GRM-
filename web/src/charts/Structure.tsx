// Population structure for one species, in three views over one selection.
// Blocks, the default, draws tens of equal columns (lineageBlocks.ts) under
// the dendrogram on block positions; each cell splits its block's genomes
// into resistant, intermediate, susceptible and not tested. Genomes opens
// single genomes at no less than MIN_GENOME_PX a column, with a minimap of
// every genome, a draggable, arrow-key brush, Previous/Next and the window's
// position. Table lists the map drug by block, then every genome's calls.
// Click, tap or Enter pins a block or a genome; the inspector beside the
// figure (below it when the width cannot hold both) keeps its details until
// Escape or Clear. Hover only previews an id and a value. The figure is one
// tab stop, a grid whose active cell carries the values in words; marks are
// canvas, labels are HTML.

import { useCallback, useEffect, useId, useMemo, useRef, useState, type KeyboardEvent, type MouseEvent, type PointerEvent, type ReactNode } from "react";
import type { SpeciesReport } from "../data/schema";
import { capitalize, distance, int, plural, stLabel } from "../data/format";
import { dendrogram, type DendroLink } from "../data/derive/layout";
import { FOCUS_RING, GhostButton } from "../ui/primitives";
import { ChartTip, keyboardFocus, labelWidth, prepareCanvas, stepCursor, useCanvasTokens, useElementWidth, type CanvasTokens } from "./kit";
import {
  BLOCK_MIN_GENOMES,
  blockAt,
  blockCalls,
  blockDendrogram,
  blockOfLineage,
  blockRecords,
  blockSize,
  callsByState,
  drugTotals,
  lineageBlocks,
  lineageSizes,
  resistantLineages,
  rsOf,
  treeBreak,
  type BlockLink,
  type CallCounts,
  type LineageBlock,
  type TreeBreak,
} from "./lineageBlocks";
import { LinkedTable, type LinkedColumn } from "./LinkedTable";

export type StructureView = "blocks" | "genomes" | "table";

/** A pinned block (by column) or genome (by leaf-order position). */
export type StructurePick = { kind: "block"; index: number } | { kind: "genome"; pos: number } | null;

/** Everything the views draw from, computed once per species file. */
export interface StructureModel {
  n: number;
  blocks: LineageBlock[];
  /** [block][drug] */
  calls: CallCounts[][];
  links: BlockLink[];
  genomeLinks: DendroLink[];
  maxHeight: number;
  brk: TreeBreak | null;
  singletons: number;
  named: number;
  /** Per drug column: genomes called R or S, and genomes called I. */
  rs: number[];
  inter: number[];
}

export const structureModel = (report: SpeciesReport): StructureModel => {
  const blocks = lineageBlocks(report);
  const d = dendrogram(report.merges, report.leaf_order);
  return {
    n: report.genomes.length,
    blocks,
    calls: blockCalls(report, blocks),
    links: blockDendrogram(report, blocks),
    genomeLinks: d.links,
    maxHeight: d.maxHeight,
    brk: treeBreak(report),
    singletons: [...lineageSizes(report).values()].filter((v) => v === 1).length,
    named: blocks.filter((b) => b.named).length,
    rs: report.drugs.map((_, j) => rsOf(drugTotals(report, j))),
    inter: report.drugs.map((_, j) => drugTotals(report, j).i),
  };
};

/** A drug's calls in words, R or S first: "1,843 R or S · 66 I". */
const callsWords = (model: StructureModel, j: number): string => `${int(model.rs[j] ?? 0)} R or S${(model.inter[j] ?? 0) > 0 ? ` · ${int(model.inter[j] ?? 0)} I` : ""}`;

/** A block's name: its lineage, or how many lineages it pools. */
export const blockTitle = (b: LineageBlock): string => (b.lineages.length === 1 ? `Lineage ${b.lineages[0]}` : `${int(b.lineages.length)} lineages`);

export const STATE_WORD: Record<keyof CallCounts, string> = { r: "resistant", i: "intermediate", s: "susceptible", untested: "not tested" };
const STATES = ["r", "i", "s", "untested"] as const;
const SWATCH: Record<keyof CallCounts, string> = {
  r: "bg-[color:var(--series-r)]",
  i: "bg-[color:var(--series-neutral)]",
  s: "bg-[color:var(--series-s)]",
  untested: "bg-surface-control ring-1 ring-inset ring-hair",
};

/** The legend's swatch for a state; keeps its color in forced colors so it
 *  still matches the canvas, which forced colors leave alone. */
export const StateSwatch = ({ state }: { state: keyof CallCounts }) => (
  <span aria-hidden="true" className={`inline-block h-2.5 w-3.5 shrink-0 rounded-sm [forced-color-adjust:none] ${SWATCH[state]}`} />
);

/** The four-state legend, with counts when given. */
export const StateLegend = ({ counts }: { counts?: CallCounts }) => (
  <ul className="flex flex-wrap gap-x-5 gap-y-1 text-control text-ink" aria-label="Call states">
    {STATES.map((k) => (
      <li key={k} className="inline-flex items-center gap-1.5">
        <StateSwatch state={k} />
        <span className="text-muted-text">{capitalize(STATE_WORD[k])}</span>
        {counts !== undefined ? <span className="font-medium tabular-nums">{int(counts[k])}</span> : null}
      </li>
    ))}
  </ul>
);

const CALL_STATE: Record<string, keyof CallCounts> = { R: "r", I: "i", S: "s" };
const stateOf = (call: string | undefined): keyof CallCounts => (call === undefined ? "untested" : (CALL_STATE[call] ?? "untested"));
const CALL_LETTER: Record<keyof CallCounts, string> = { r: "R", i: "I", s: "S", untested: "–" };

// ---- Geometry ---------------------------------------------------------------

const TREE_H = 112;
const MAIN_TOP = 4;
// With an axis break: outliers in a band at the top, a cleared stripe for the
// break, and the main tree under it.
const BAND_TOP = 3;
const BAND_H = 12;
const STRIPE_TOP = 19;
const STRIPE_H = 5;
const MAIN_TOP_BROKEN = 28;
const HEAD_GAP = 6;
const HEAD_H = 16;
const COLHEAD_H = 18;
const ROW_H = 16;
const MIN_BLOCK_PX = 24;
/** No genome column narrower than this when genomes are open. */
export const MIN_GENOME_PX = 4;
const INSPECTOR_PX = 288;
const SIDE_GAP = 24;
const LABEL_FONT_PX = 6.2;

interface TreeScale {
  y(h: number): number;
  broken: boolean;
  /** Height at the top of the main panel. */
  mainMax: number;
}

/** Square-root height scale, so merges among near-identical genomes, which
 *  carry the lineage structure, stay visible beside the deep splits; with a
 *  break, the main panel ends at the highest merge under it and the
 *  outliers take a thin band above. */
const treeScale = (maxHeight: number, brk: TreeBreak | null): TreeScale => {
  if (brk === null || !(brk.below > 0) || brk.outliers.length === 0) {
    const top = Math.max(maxHeight, 1e-12);
    return { y: (h) => TREE_H - Math.sqrt(Math.min(h, top) / top) * (TREE_H - MAIN_TOP), broken: false, mainMax: maxHeight };
  }
  const hi = brk.outliers[0]!.height;
  return {
    y: (h) =>
      h <= brk.below
        ? TREE_H - Math.sqrt(h / brk.below) * (TREE_H - MAIN_TOP_BROKEN)
        : BAND_TOP + BAND_H - ((h - brk.below) / Math.max(hi - brk.below, 1e-12)) * BAND_H,
    broken: true,
    mainMax: brk.below,
  };
};

const rowsTop = (headH: number) => TREE_H + HEAD_GAP + headH + HEAD_H + COLHEAD_H;

const drawTree = (ctx: CanvasRenderingContext2D, links: readonly { xa: number; xb: number; ya: number; yb: number; y: number }[], xOf: (x: number) => number, scale: TreeScale, width: number, tokens: CanvasTokens) => {
  ctx.save();
  ctx.beginPath();
  ctx.rect(0, 0, width, TREE_H + 1);
  ctx.clip();
  ctx.strokeStyle = tokens.axis;
  ctx.lineWidth = 1;
  ctx.beginPath();
  links.forEach((l) => {
    const xa = xOf(l.xa);
    const xb = xOf(l.xb);
    const y = scale.y(l.y);
    ctx.moveTo(xa, scale.y(l.ya));
    ctx.lineTo(xa, y);
    ctx.lineTo(xb, y);
    ctx.lineTo(xb, scale.y(l.yb));
  });
  ctx.stroke();
  ctx.restore();
  if (!scale.broken) return;
  // The break: a cleared stripe through every stem, marked at the axis
  // with two slants.
  ctx.clearRect(0, STRIPE_TOP, width, STRIPE_H);
  ctx.strokeStyle = tokens.ink;
  ctx.lineWidth = 1.25;
  ctx.beginPath();
  [0, 4].forEach((dx) => {
    ctx.moveTo(dx, STRIPE_TOP + STRIPE_H + 1);
    ctx.lineTo(dx + 6, STRIPE_TOP - 1);
  });
  ctx.stroke();
};

const callColor = (call: string | undefined, tokens: CanvasTokens): string =>
  call === "R" ? tokens.r : call === "S" ? tokens.s : call === "I" ? tokens.neutral : tokens.fill;

// ---- The figure -----------------------------------------------------------------

interface Cursor {
  col: number;
  /** -1 is the tree and the header rows; 0.. are drug rows in `order`. */
  row: number;
}

export const StructureExplorer = ({
  report,
  model,
  order,
  drug,
  view,
  pick,
  onPick,
  onView,
  caption,
}: {
  report: SpeciesReport;
  model: StructureModel;
  /** Drug columns, top row first. */
  order: readonly number[];
  /** The map drug's column. */
  drug: number;
  view: StructureView;
  pick: StructurePick;
  onPick(p: StructurePick): void;
  onView(v: StructureView): void;
  /** Sits right under the figure, above the inspector when that stacks. */
  caption?: ReactNode;
}) => {
  const [outerRef, outerW] = useElementWidth(1200);
  const [scrollerRef, viewW] = useElementWidth(800);
  const scroller = useRef<HTMLDivElement | null>(null);
  const canvas = useRef<HTMLCanvasElement | null>(null);
  const tokens = useCanvasTokens();
  const cellId = useId();
  const { n, blocks } = model;
  const B = blocks.length;
  const D = order.length;
  const drugRow = Math.max(0, order.indexOf(drug));

  // Genome window: its first leaf position.
  const [ws, setWs] = useState(() => (pick?.kind === "block" ? (blocks[pick.index]?.start ?? 0) : pick?.kind === "genome" ? pick.pos : 0));
  const [blockCursor, setBlockCursor] = useState<Cursor>({ col: 0, row: drugRow });
  const [genomeCursor, setGenomeCursor] = useState<Cursor>({ col: 0, row: drugRow });
  const [hover, setHover] = useState<{ col: number; row: number; x: number; y: number } | null>(null);
  const [focused, setFocused] = useState(false);
  // Escape hides the keyboard cursor's preview until the cursor moves.
  const [tipOff, setTipOff] = useState(false);
  const [scrollX, setScrollX] = useState(0);
  const setScroller = useCallback(
    (el: HTMLDivElement | null) => {
      scroller.current = el;
      scrollerRef(el);
    },
    [scrollerRef],
  );

  // Label column: drug names in full, genomes called R or S, and genomes
  // called I where the species has any. Estimated only to decide whether
  // the inspector fits beside the figure.
  const anyI = model.inter.some((v) => v > 0);
  const nameW = Math.max(56, ...order.map((j) => labelWidth(capitalize(report.drugs[j] ?? ""))));
  const side = outerW >= INSPECTOR_PX + SIDE_GAP + nameW + (anyI ? 96 : 64) + B * MIN_BLOCK_PX;

  const genomesView = view === "genomes";
  // The genome window always fits the visible width, whichever view is
  // drawn, so a block opened from the inspector lands inside it.
  const gPlotW = Math.max(viewW, 120);
  const win = Math.max(1, Math.min(n, Math.floor(gPlotW / MIN_GENOME_PX)));
  const gColW = gPlotW / win;
  const plotW = genomesView ? gPlotW : Math.max(viewW, B * MIN_BLOCK_PX);
  const colW = B === 0 ? plotW : plotW / B;
  const start = Math.max(0, Math.min(ws, n - win));

  // Block header labels: horizontal when every one fits its column.
  const longest = Math.max(1, ...blocks.map((b) => Math.max(int(blockSize(b)).length, int(b.lineages.length).length)));
  const flat = genomesView || longest * LABEL_FONT_PX + 4 <= colW;
  const headH = flat ? HEAD_H : Math.ceil(longest * LABEL_FONT_PX + 8);
  const head1 = TREE_H + HEAD_GAP;
  const head2 = head1 + headH;
  const top = rowsTop(headH);
  const height = top + D * ROW_H + 4;
  const scale = useMemo(() => treeScale(model.maxHeight, model.brk), [model]);

  useEffect(() => {
    const el = canvas.current;
    if (el === null || n === 0 || view === "table") return;
    const ctx = prepareCanvas(el, plotW, height);
    if (ctx === null) return;
    ctx.clearRect(0, 0, plotW, height);
    if (!genomesView) {
      drawTree(ctx, model.links, (x) => x * colW, scale, plotW, tokens);
      order.forEach((j, row) => {
        const y = top + row * ROW_H + 1;
        blocks.forEach((b, bi) => {
          const c = model.calls[bi]![j]!;
          const size = blockSize(b);
          const x0 = bi * colW + 1;
          const w = Math.max(1, colW - 2);
          ctx.fillStyle = tokens.fill;
          ctx.fillRect(x0, y, w, ROW_H - 2);
          let x = x0;
          ([
            [c.r, tokens.r],
            [c.i, tokens.neutral],
            [c.s, tokens.s],
          ] as const).forEach(([count, color]) => {
            if (count === 0) return;
            const sw = (w * count) / size;
            ctx.fillStyle = color;
            ctx.fillRect(x, y, sw, ROW_H - 2);
            x += sw;
          });
        });
      });
      return;
    }
    const end = start + win;
    drawTree(
      ctx,
      model.genomeLinks.filter((l) => Math.max(l.xa, l.xb) >= start - 1 && Math.min(l.xa, l.xb) <= end),
      (x) => (x - start + 0.5) * gColW,
      scale,
      plotW,
      tokens,
    );
    const gap = gColW >= 6 ? 1 : 0;
    // Lineage band, then block band: alternate shades per run.
    let shade = 0;
    for (let pos = start; pos < end; pos += 1) {
      const g = report.genomes[report.leaf_order[pos]!]!;
      const prev = pos === start ? undefined : report.genomes[report.leaf_order[pos - 1]!];
      if (prev !== undefined && prev.lineage !== g.lineage) shade = 1 - shade;
      ctx.fillStyle = shade === 0 ? tokens.ink : tokens.muted;
      ctx.globalAlpha = shade === 0 ? 0.55 : 0.2;
      ctx.fillRect((pos - start) * gColW, head1 + 3, gColW + 0.35, headH - 6);
      const b = blockAt(blocks, pos);
      ctx.fillStyle = tokens.blue;
      ctx.globalAlpha = (b?.index ?? 0) % 2 === 0 ? 0.6 : 0.25;
      ctx.fillRect((pos - start) * gColW, head2 + 3, gColW + 0.35, HEAD_H - 6);
      ctx.globalAlpha = 1;
      const profile = g.profile;
      order.forEach((j, row) => {
        ctx.fillStyle = callColor(profile[j], tokens);
        ctx.fillRect((pos - start) * gColW, top + row * ROW_H + 1, gColW - gap + (gap === 0 ? 0.35 : 0), ROW_H - 2);
      });
    }
  }, [report, model, order, view, plotW, height, colW, gColW, start, win, top, head1, head2, headH, scale, tokens]);

  // Keep the keyboard's genome inside the window when the window moves.
  const gCursor: Cursor = { col: Math.max(start, Math.min(genomeCursor.col, start + win - 1)), row: genomeCursor.row };
  const cursor = genomesView ? gCursor : { col: Math.min(blockCursor.col, Math.max(0, B - 1)), row: blockCursor.row };

  const revealBlock = (col: number) => {
    const el = scroller.current;
    if (el === null) return;
    const x0 = col * colW;
    if (x0 < el.scrollLeft) el.scrollLeft = x0;
    else if (x0 + colW > el.scrollLeft + el.clientWidth) el.scrollLeft = x0 + colW - el.clientWidth;
  };

  const moveWindowTo = (pos: number) => setWs(Math.max(0, Math.min(pos, n - win)));

  const pinCol = (col: number) => {
    if (genomesView) onPick(pick?.kind === "genome" && pick.pos === col ? null : { kind: "genome", pos: col });
    else onPick(pick?.kind === "block" && pick.index === col ? null : { kind: "block", index: col });
  };

  const at = (clientX: number, clientY: number, rect: DOMRect): Cursor | null => {
    const x = clientX - rect.left;
    const y = clientY - rect.top;
    if (x < 0 || x > plotW || y < 0 || y > height) return null;
    const row = y >= top ? Math.floor((y - top) / ROW_H) : -1;
    const r = row >= 0 && row < D ? row : -1;
    if (genomesView) return { col: Math.min(start + win - 1, start + Math.floor(x / gColW)), row: r };
    return { col: Math.min(B - 1, Math.floor(x / colW)), row: r };
  };

  const onPointerMove = (e: PointerEvent<HTMLDivElement>) => {
    if (e.pointerType === "touch") return;
    const c = at(e.clientX, e.clientY, e.currentTarget.getBoundingClientRect());
    setHover(c === null ? null : { ...c, x: e.clientX - e.currentTarget.getBoundingClientRect().left, y: e.clientY - e.currentTarget.getBoundingClientRect().top });
  };

  const onClick = (e: MouseEvent<HTMLDivElement>) => {
    const c = at(e.clientX, e.clientY, e.currentTarget.getBoundingClientRect());
    if (c === null) return;
    if (genomesView) setGenomeCursor(c);
    else setBlockCursor(c);
    pinCol(c.col);
  };

  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key === "Escape") {
      // Escape clears the pin first, then a preview; with neither it is not
      // the figure's, so a dialog behind it can close.
      if (pick !== null) {
        e.preventDefault();
        onPick(null);
      } else if (hover !== null || (focused && !tipOff)) {
        e.preventDefault();
        setHover(null);
        setTipOff(true);
      }
      return;
    }
    if (e.key === "Enter" || e.key === " ") {
      e.preventDefault();
      pinCol(cursor.col);
      return;
    }
    let next: Cursor;
    if (e.key === "ArrowUp" || e.key === "ArrowDown") {
      next = { col: cursor.col, row: Math.min(D - 1, Math.max(-1, cursor.row + (e.key === "ArrowDown" ? 1 : -1))) };
    } else {
      const col = stepCursor(e.key, cursor.col, genomesView ? n : B);
      if (col === null) return;
      next = { col, row: cursor.row };
    }
    e.preventDefault();
    setHover(null);
    setTipOff(false);
    if (genomesView) {
      setGenomeCursor(next);
      if (next.col < start) moveWindowTo(next.col);
      else if (next.col >= start + win) moveWindowTo(next.col - win + 1);
    } else {
      setBlockCursor(next);
      revealBlock(next.col);
    }
  };

  // What the active cell says: the values a hover previews, in full.
  const cellText = (c: Cursor): string => {
    const drugName = c.row >= 0 ? capitalize(report.drugs[order[c.row]!] ?? "") : "";
    if (genomesView) {
      const g = report.genomes[report.leaf_order[c.col]!];
      if (g === undefined) return "";
      const pinned = pick?.kind === "genome" && pick.pos === c.col ? " Pinned." : "";
      if (c.row >= 0) return `${g.id}, ${drugName}: ${STATE_WORD[stateOf(g.profile[order[c.row]!])]}.${pinned}`;
      const by = callsByState(g.profile, order);
      return `${g.id}, lineage ${g.lineage}, genome ${int(c.col + 1)} of ${int(n)} in tree order: ${STATES.map((k) => `${int(by[k].length)} ${STATE_WORD[k]}`).join(", ")}.${pinned}`;
    }
    const b = blocks[c.col];
    if (b === undefined) return "";
    const pinned = pick?.kind === "block" && pick.index === c.col ? " Pinned." : "";
    if (c.row >= 0) {
      const k = model.calls[c.col]![order[c.row]!]!;
      return `${blockTitle(b)}, ${drugName}: ${int(k.r)} resistant, ${int(k.i)} intermediate, ${int(k.s)} susceptible, ${int(k.untested)} not tested.${pinned}`;
    }
    return `Block ${c.col + 1} of ${B}, ${blockTitle(b)}: ${plural(blockSize(b), "genome")}, ${int(b.start + 1)} to ${int(b.end)} of ${int(n)} in tree order.${pinned}`;
  };

  // The short preview a hover shows: an id and a value.
  const hoverLines = (c: Cursor) => {
    if (genomesView) {
      const g = report.genomes[report.leaf_order[c.col]!];
      if (g === undefined) return null;
      return [{ term: g.id, text: c.row >= 0 ? `${capitalize(report.drugs[order[c.row]!] ?? "")}: ${STATE_WORD[stateOf(g.profile[order[c.row]!])]}` : `Lineage ${g.lineage}` }];
    }
    const b = blocks[c.col];
    if (b === undefined) return null;
    if (c.row < 0) return [{ term: blockTitle(b), text: plural(blockSize(b), "genome") }];
    const k = model.calls[c.col]![order[c.row]!]!;
    return [{ term: blockTitle(b), text: `${capitalize(report.drugs[order[c.row]!] ?? "")}: ${int(k.r)} R · ${int(k.i)} I · ${int(k.s)} S · ${int(k.untested)} not tested` }];
  };

  const pinnedBlock = pick?.kind === "block" ? blocks[pick.index] : pick?.kind === "genome" ? blockAt(blocks, pick.pos) : undefined;
  const ringFor = (): { left: number; width: number } | null => {
    if (genomesView) {
      if (pick?.kind === "genome" && pick.pos >= start && pick.pos < start + win) return { left: (pick.pos - start) * gColW, width: gColW };
      if (pick?.kind === "block" && pinnedBlock !== undefined) {
        const a = Math.max(pinnedBlock.start, start);
        const z = Math.min(pinnedBlock.end, start + win);
        return z > a ? { left: (a - start) * gColW, width: (z - a) * gColW } : null;
      }
      return null;
    }
    return pinnedBlock === undefined ? null : { left: pinnedBlock.index * colW, width: colW };
  };
  const ring = ringFor();
  const hoverBox = hover === null ? null : genomesView ? { left: (hover.col - start) * gColW, width: gColW } : { left: hover.col * colW, width: colW };
  const cursorBox = !focused ? null : genomesView ? { left: (cursor.col - start) * gColW, width: gColW } : { left: cursor.col * colW, width: colW };
  // The preview follows the pointer, else the keyboard's cell.
  const preview =
    hover !== null
      ? { lines: hoverLines(hover), x: hover.x, y: hover.y }
      : cursorBox !== null && !tipOff
        ? { lines: hoverLines(cursor), x: cursorBox.left + cursorBox.width, y: cursor.row >= 0 ? top + cursor.row * ROW_H : head1 }
        : null;
  const visibleFirst = Math.min(B, Math.floor(scrollX / colW) + 1);
  const visibleLast = Math.min(B, Math.floor((scrollX + viewW - 1) / colW) + 1);

  // The label column: the tree's axis, the header rows' names, then the
  // drug table (full names, R or S counts, I counts where any) in flow, so
  // it sets the width.
  const labels = (
    <div aria-hidden="true" className="relative min-w-[7rem] select-none text-tiny text-muted-text" style={{ height }}>
      <span className="absolute right-0 leading-none" style={{ top: TREE_H - 9 }}>0</span>
      <span className="absolute right-0 whitespace-nowrap leading-none" style={{ top: scale.broken ? MAIN_TOP_BROKEN - 4 : MAIN_TOP - 2 }}>
        {distance(scale.mainMax)}
      </span>
      {scale.broken && model.brk !== null ? (
        <span className="absolute right-0 whitespace-nowrap leading-none" style={{ top: BAND_TOP - 1 }}>
          {distance(model.brk.outliers[0]!.height)}
        </span>
      ) : null}
      <span className="absolute right-0 whitespace-nowrap text-right leading-tight" style={{ top: TREE_H / 2 - 10 }}>
        Mash distance,
        <br />
        square-root scale
      </span>
      <span className="absolute right-0 whitespace-nowrap" style={{ top: head1, lineHeight: `${headH}px` }}>
        {genomesView ? "Lineages" : "Genomes"}
      </span>
      <span className="absolute right-0 whitespace-nowrap" style={{ top: head2, lineHeight: `${HEAD_H}px` }}>
        {genomesView ? "Blocks" : "Lineages"}
      </span>
      <div className={`grid gap-x-3 ${anyI ? "grid-cols-[minmax(0,1fr)_auto_auto]" : "grid-cols-[minmax(0,1fr)_auto]"}`} style={{ paddingTop: top - COLHEAD_H }}>
        <span className="text-right font-medium" style={{ lineHeight: `${COLHEAD_H}px` }}>Drug</span>
        <span className="whitespace-nowrap text-right font-medium" style={{ lineHeight: `${COLHEAD_H}px` }}>R or S</span>
        {anyI ? (
          <span className="text-right font-medium" style={{ lineHeight: `${COLHEAD_H}px` }}>
            I
          </span>
        ) : null}
        {order.map((j) => {
          const on = j === drug;
          const tone = on ? "font-semibold text-ink" : "";
          return (
            <div key={j} className="contents">
              <span className={`whitespace-nowrap text-right text-[12px] leading-[16px] ${tone}`}>{capitalize(report.drugs[j] ?? "")}</span>
              <span className={`text-right text-[12px] leading-[16px] tabular-nums ${tone}`}>{int(model.rs[j] ?? 0)}</span>
              {anyI ? <span className={`text-right text-[12px] leading-[16px] tabular-nums ${tone}`}>{int(model.inter[j] ?? 0)}</span> : null}
            </div>
          );
        })}
      </div>
    </div>
  );

  const figure = (
    <div className="grid grid-cols-[max-content_minmax(0,1fr)] gap-x-3">
      {genomesView ? (
        <>
          <div className="whitespace-nowrap pt-1 text-right text-tiny text-muted-text">All genomes</div>
          <GenomeMinimap report={report} model={model} drug={drug} start={start} win={win} onMove={moveWindowTo} />
        </>
      ) : null}
      {labels}
      <div
        ref={setScroller}
        onScroll={(e) => setScrollX(e.currentTarget.scrollLeft)}
        className="scroll-thin min-w-0 overflow-x-auto overscroll-x-contain"
      >
        <div className="relative" style={{ width: plotW, height }}>
          <div
            role="grid"
            tabIndex={0}
            aria-label={
              genomesView
                ? `Genomes ${int(start + 1)} to ${int(start + win)} of ${int(n)} in tree order, by ${plural(D, "drug")}. Arrow keys move; click or Enter to pin details; Escape clears.`
                : `${plural(B, "lineage block")} of ${int(n)} genomes by ${plural(D, "drug")}. Arrow keys move; click or Enter to pin details; Escape clears.`
            }
            aria-rowcount={D + 1}
            aria-colcount={genomesView ? n : B}
            aria-activedescendant={`${cellId}-${view}-${cursor.row}-${cursor.col}`}
            onKeyDown={onKey}
            onFocus={(e) => {
              setFocused(keyboardFocus(e.currentTarget));
              setTipOff(false);
            }}
            onBlur={() => setFocused(false)}
            onPointerMove={onPointerMove}
            onPointerLeave={() => setHover(null)}
            onClick={onClick}
            className="cursor-pointer rounded-sm focus-visible:outline focus-visible:outline-2 focus-visible:-outline-offset-2 focus-visible:outline-blue"
          >
            <canvas ref={canvas} aria-hidden="true" className="block" />
            <div role="row" aria-rowindex={cursor.row + 2} className="sr-only">
              <div role="gridcell" id={`${cellId}-${view}-${cursor.row}-${cursor.col}`} aria-colindex={cursor.col + 1}>
                {cellText(cursor)}
              </div>
            </div>
          </div>
          {/* The map drug's outline is HTML over the canvas, so a new choice
              glides to its row on the slow clock instead of jumping. */}
          <div
            aria-hidden="true"
            className="pointer-events-none absolute left-0 top-0 rounded-[2px] border-[1.5px] border-ink transition-transform duration-slow ease-out"
            style={{ width: plotW, height: ROW_H, transform: `translateY(${top + drugRow * ROW_H}px)` }}
          />
          {!genomesView
            ? blocks.map((b) => (
                <div key={b.index} aria-hidden="true" className="pointer-events-none absolute text-[11px] tabular-nums" style={{ left: b.index * colW, width: colW, top: head1 }}>
                  <div
                    className={`flex items-center justify-center ${flat ? "" : "[writing-mode:vertical-rl] rotate-180"} ${b.named ? "font-semibold text-ink" : "text-muted-text"}`}
                    style={{ height: headH }}
                  >
                    {int(blockSize(b))}
                  </div>
                  <div className="flex items-center justify-center text-muted-text" style={{ height: HEAD_H }}>
                    {int(b.lineages.length)}
                  </div>
                </div>
              ))
            : null}
          {hoverBox !== null ? <div aria-hidden="true" className="pointer-events-none absolute bg-ink/[0.07]" style={{ left: hoverBox.left, width: hoverBox.width, top: head1, height: height - head1 }} /> : null}
          {ring !== null ? (
            <div aria-hidden="true" className="pointer-events-none absolute rounded-[3px] border-2 border-ink" style={{ left: ring.left - 1, width: ring.width + 2, top: head1 - 2, height: height - head1 }} />
          ) : null}
          {cursorBox !== null ? (
            <div
              aria-hidden="true"
              className="pointer-events-none absolute rounded-[2px] outline outline-2 outline-blue"
              style={cursor.row >= 0 ? { left: cursorBox.left, width: cursorBox.width, top: top + cursor.row * ROW_H, height: ROW_H } : { left: cursorBox.left, width: cursorBox.width, top: head1, height: top - head1 }}
            />
          ) : null}
          {preview !== null && preview.lines !== null ? <ChartTip x={preview.x} y={preview.y} width={scrollX + viewW} lines={preview.lines} /> : null}
        </div>
      </div>
    </div>
  );

  return (
    <div ref={outerRef} className={side ? "grid grid-cols-[minmax(0,1fr)_18rem] items-start gap-6" : "grid grid-cols-1 gap-6"}>
      <div className="min-w-0">
        {view === "table" ? (
          <StructureTables report={report} model={model} order={order} drug={drug} pick={pick} onPick={onPick} />
        ) : (
          <>
            <div className="mb-3 flex flex-wrap items-center justify-between gap-x-6 gap-y-2">
              <StateLegend />
              {genomesView ? (
                <GenomeWindowBar n={n} start={start} win={win} blocks={blocks} onMove={moveWindowTo} />
              ) : (
                <p className="text-control text-muted-text">
                  {`${plural(B, "block")}: ${plural(model.named, "lineage")} of at least ${BLOCK_MIN_GENOMES} genomes, ${plural(B - model.named, "pool")} of smaller ones`}
                </p>
              )}
            </div>
            {figure}
            {!genomesView && plotW > viewW + 1 ? (
              <p className="mt-2 text-control tabular-nums text-muted-text" aria-live="polite">{`Blocks ${visibleFirst}–${visibleLast} of ${B} in view; scroll sideways for the rest.`}</p>
            ) : null}
            {model.brk !== null ? <BreakNote brk={model.brk} /> : null}
          </>
        )}
        {caption !== undefined ? <p className="prose-measure mt-3 text-note text-muted-text">{caption}</p> : null}
      </div>
      <Inspector report={report} model={model} order={order} drug={drug} pick={pick} onPick={onPick} view={view} side={side} onShowGenomes={(pos) => { moveWindowTo(pos); onView("genomes"); }} />
    </div>
  );
};

/** The outliers above the axis break, named for a data-quality check. */
const BreakNote = ({ brk }: { brk: TreeBreak }) => {
  const parts = brk.outliers.map((o) => {
    const who = o.genomes.length === 1 ? `genome ${o.genomes[0]}` : `${o.genomes.length} genomes (${o.genomes.join(", ")})`;
    return `${who} at ${distance(o.height)}`;
  });
  return (
    <p className="prose-measure mt-2 text-note text-muted-text">
      {`The height axis breaks above Mash distance ${distance(brk.below)}, the highest merge under it. Above the break, ${parts.length === 1 ? parts[0] : `${parts.slice(0, -1).join(", ")} and ${parts[parts.length - 1]}`} ${brk.outliers.length === 1 && brk.outliers[0]!.genomes.length === 1 ? "joins" : "join"} the rest of the tree. Joins this distant are worth a data-quality check.`}
    </p>
  );
};

// ---- Genome window: position, Previous/Next, minimap brush -------------------

const GenomeWindowBar = ({ n, start, win, blocks, onMove }: { n: number; start: number; win: number; blocks: readonly LineageBlock[]; onMove(pos: number): void }) => {
  const first = blockAt(blocks, start);
  const last = blockAt(blocks, Math.min(n - 1, start + win - 1));
  return (
    <div className="flex flex-wrap items-center gap-2">
      <GhostButton onClick={() => onMove(start - win)} disabled={start <= 0} ariaLabel={`Previous ${int(win)} genomes`}>
        Previous
      </GhostButton>
      <GhostButton onClick={() => onMove(start + win)} disabled={start + win >= n} ariaLabel={`Next ${int(win)} genomes`}>
        Next
      </GhostButton>
      <p aria-live="polite" className="text-control tabular-nums text-ink">
        {`Genomes ${int(start + 1)}–${int(Math.min(n, start + win))} of ${int(n)}`}
        {first !== undefined && last !== undefined ? <span className="text-muted-text">{` · blocks ${first.index + 1}–${last.index + 1} of ${blocks.length}`}</span> : null}
      </p>
    </div>
  );
};

/** Every genome in tree order, colored by its call for the map drug, with a
 *  brush over the open window: drag it, tap elsewhere to move it there, or
 *  use the arrow keys (Page keys move a whole window, Home and End go to the
 *  ends). */
const GenomeMinimap = ({ report, model, drug, start, win, onMove }: { report: SpeciesReport; model: StructureModel; drug: number; start: number; win: number; onMove(pos: number): void }) => {
  const [ref, width] = useElementWidth(600);
  const canvas = useRef<HTMLCanvasElement | null>(null);
  const tokens = useCanvasTokens();
  const drag = useRef<number | null>(null);
  const { n } = model;
  const H = 20;
  useEffect(() => {
    const el = canvas.current;
    if (el === null || n === 0) return;
    const ctx = prepareCanvas(el, width, H);
    if (ctx === null) return;
    ctx.clearRect(0, 0, width, H);
    for (let pos = 0; pos < n; pos += 1) {
      ctx.fillStyle = callColor(report.genomes[report.leaf_order[pos]!]?.profile[drug], tokens);
      ctx.fillRect((pos * width) / n, 4, width / n + 0.35, H - 8);
    }
    ctx.fillStyle = tokens.muted;
    model.blocks.forEach((b) => {
      if (b.start > 0) ctx.fillRect((b.start * width) / n, 0, 1, 4);
    });
  }, [report, model, drug, width, tokens, n]);

  const posAt = (clientX: number, rect: DOMRect) => Math.floor(((clientX - rect.left) / Math.max(rect.width, 1)) * n);
  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    const step = Math.max(1, Math.round(win / 10));
    const to =
      e.key === "ArrowLeft" || e.key === "ArrowDown"
        ? start - step
        : e.key === "ArrowRight" || e.key === "ArrowUp"
          ? start + step
          : e.key === "PageUp"
            ? start - win
            : e.key === "PageDown"
              ? start + win
              : e.key === "Home"
                ? 0
                : e.key === "End"
                  ? n
                  : null;
    if (to === null) return;
    e.preventDefault();
    onMove(to);
  };
  return (
    <div
      ref={ref}
      className="relative mb-2 touch-pan-y"
      style={{ height: H }}
      onPointerDown={(e) => {
        const rect = e.currentTarget.getBoundingClientRect();
        const pos = posAt(e.clientX, rect);
        drag.current = pos >= start && pos < start + win ? pos - start : Math.floor(win / 2);
        e.currentTarget.setPointerCapture?.(e.pointerId);
        onMove(pos - drag.current);
      }}
      onPointerMove={(e) => {
        if (drag.current === null) return;
        onMove(posAt(e.clientX, e.currentTarget.getBoundingClientRect()) - drag.current);
      }}
      onPointerUp={() => {
        drag.current = null;
      }}
      onPointerCancel={() => {
        drag.current = null;
      }}
    >
      <canvas ref={canvas} aria-hidden="true" className="block" />
      <div
        role="slider"
        tabIndex={0}
        aria-label="Open genomes, as a window over all genomes in tree order"
        aria-valuemin={1}
        aria-valuemax={Math.max(1, n - win + 1)}
        aria-valuenow={start + 1}
        aria-valuetext={`Genomes ${int(start + 1)} to ${int(Math.min(n, start + win))} of ${int(n)}`}
        onKeyDown={onKey}
        className={`absolute top-0 h-full cursor-grab rounded-[3px] border-2 border-ink bg-ink/10 active:cursor-grabbing ${FOCUS_RING}`}
        style={{ left: `${(start / Math.max(n, 1)) * 100}%`, width: `max(6px, ${(win / Math.max(n, 1)) * 100}%)` }}
      />
    </div>
  );
};

// ---- The inspector --------------------------------------------------------------

const Group = ({ title, children }: { title: string; children: ReactNode }) => (
  <div className="min-w-0">
    <h4 className="text-tiny font-semibold uppercase tracking-[0.06em] text-muted-text">{title}</h4>
    <div className="mt-1.5">{children}</div>
  </div>
);

const Inspector = ({
  report,
  model,
  order,
  drug,
  pick,
  onPick,
  view,
  side,
  onShowGenomes,
}: {
  report: SpeciesReport;
  model: StructureModel;
  order: readonly number[];
  drug: number;
  pick: StructurePick;
  onPick(p: StructurePick): void;
  view: StructureView;
  side: boolean;
  onShowGenomes(pos: number): void;
}) => {
  const titleId = useId();
  const drugName = capitalize(report.drugs[drug] ?? "");
  const { n, blocks } = model;
  const frame = side ? "border-l border-hair pl-6" : "border-t border-hair pt-5";
  const block = pick?.kind === "block" ? blocks[pick.index] : undefined;
  const genomePos = pick?.kind === "genome" ? pick.pos : null;
  const genome = genomePos === null ? undefined : report.genomes[report.leaf_order[genomePos]!];

  let kicker = "Map drug, nothing pinned";
  let title = drugName;
  let body: ReactNode;
  if (block !== undefined) {
    kicker = block.named ? `Block ${block.index + 1} of ${blocks.length} · one lineage` : `Block ${block.index + 1} of ${blocks.length} · pooled lineages`;
    title = blockTitle(block);
    const sizes = lineageSizes(report);
    const rec = blockRecords(report, block);
    const size = blockSize(block);
    const singles = block.lineages.filter((l) => sizes.get(l) === 1).length;
    body = (
      <>
        <p className="text-control text-ink">{`${plural(size, "genome")}, ${int(block.start + 1)}–${int(block.end)} of ${int(n)} in tree order.`}</p>
        <div className={`mt-4 grid gap-5 ${side ? "grid-cols-1" : "sm:grid-cols-2"}`}>
          <Group title="Calls by drug">
            <table className="w-full max-w-[24rem] border-separate border-spacing-0 text-tiny tabular-nums">
              <thead>
                <tr className="text-muted-text">
                  <th scope="col" className="pb-1 text-left font-medium">Drug</th>
                  {STATES.map((k) => (
                    <th key={k} scope="col" className="pb-1 pl-2 text-right font-medium">
                      <span aria-hidden="true">{CALL_LETTER[k]}</span>
                      <span className="sr-only">{capitalize(STATE_WORD[k])}</span>
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {order.map((j) => {
                  const c = model.calls[block.index]![j]!;
                  const weight = j === drug ? "font-semibold" : "font-normal";
                  return (
                    <tr key={j} className="text-ink">
                      <th scope="row" className={`py-0.5 pr-2 text-left ${weight}`}>{capitalize(report.drugs[j] ?? "")}</th>
                      {STATES.map((k) => (
                        <td key={k} className={`py-0.5 pl-2 text-right ${weight}`}>{int(c[k])}</td>
                      ))}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </Group>
          <div className="space-y-4">
            <Group title={block.lineages.length === 1 ? "Lineage" : `${int(block.lineages.length)} lineages`}>
              {block.lineages.length === 1 ? (
                <p className="text-control text-ink">{`Lineage ${block.lineages[0]}, ${plural(sizes.get(block.lineages[0]!) ?? 0, "genome")}.`}</p>
              ) : (
                <>
                  <p className="text-control text-ink">{`Each under ${BLOCK_MIN_GENOMES} genomes; ${int(singles)} hold one genome. Lineage (genomes):`}</p>
                  <p className="mt-1 max-h-28 overflow-y-auto text-control tabular-nums text-muted-text scroll-thin">
                    {block.lineages.map((l) => `${l} (${int(sizes.get(l) ?? 0)})`).join(", ")}
                  </p>
                </>
              )}
            </Group>
            <Group title="Recorded metadata">
              <ul className="space-y-0.5 text-control text-ink">
                <li>{`MLST: ${int(rec.mlst.recorded)} of ${int(size)} (${plural(rec.mlst.distinct, "type")})`}</li>
                <li>{`Country: ${int(rec.country.recorded)} of ${int(size)} (${plural(rec.country.distinct, "country", "countries")})`}</li>
                <li>{`Year: ${int(rec.year.recorded)} of ${int(size)}${rec.year.min !== null ? ` (${rec.year.min === rec.year.max ? rec.year.min : `${rec.year.min}–${rec.year.max}`})` : ""}`}</li>
              </ul>
            </Group>
          </div>
        </div>
        {view !== "genomes" ? (
          <div className="mt-4">
            <GhostButton onClick={() => onShowGenomes(block.start)}>Show these genomes</GhostButton>
          </div>
        ) : null}
      </>
    );
  } else if (genome !== undefined && genomePos !== null) {
    const b = blockAt(blocks, genomePos);
    kicker = `Genome ${int(genomePos + 1)} of ${int(n)} in tree order`;
    title = genome.id;
    const by = callsByState(genome.profile, order);
    body = (
      <>
        <p className="text-control text-ink">
          {[
            `Lineage ${genome.lineage}`,
            b !== undefined ? `block ${b.index + 1}` : null,
            genome.mlst !== null ? stLabel(genome.mlst) : "MLST not recorded",
            genome.country ?? "country not recorded",
            genome.year ?? "year not recorded",
          ]
            .filter((x): x is string => x !== null)
            .join(" · ")}
        </p>
        <div className={`mt-4 grid gap-4 ${side ? "grid-cols-1" : "sm:grid-cols-2"}`}>
          {STATES.map((k) => (
            <Group key={k} title={`${capitalize(STATE_WORD[k])} (${by[k].length})`}>
              <p className={`text-control ${by[k].length === 0 ? "text-muted-text" : "text-ink"}`}>
                {by[k].length === 0 ? "None" : by[k].map((j) => capitalize(report.drugs[j] ?? "")).join(", ")}
              </p>
            </Group>
          ))}
        </div>
        {view !== "genomes" ? (
          <div className="mt-4">
            <GhostButton onClick={() => onShowGenomes(Math.max(0, genomePos - 10))}>Show this genome</GhostButton>
          </div>
        ) : null}
      </>
    );
  } else {
    const totals = drugTotals(report, drug);
    const top = resistantLineages(report, drug).slice(0, 3);
    title = `${drugName} across ${int(n)} genomes`;
    body = (
      <div className={`grid gap-x-10 gap-y-4 ${side ? "grid-cols-1" : "sm:grid-cols-2 lg:grid-cols-3"}`}>
        <ul className="space-y-1 text-control text-ink">
          {STATES.map((k) => (
            <li key={k} className="flex items-center gap-2">
              <StateSwatch state={k} />
              <span className="w-24 text-muted-text">{capitalize(STATE_WORD[k])}</span>
              <span className="font-medium tabular-nums">{int(totals[k])}</span>
            </li>
          ))}
        </ul>
        {top.length > 0 ? (
          <div>
            <Group title={top.length === 1 ? "The lineage holding the resistance" : `The ${top.length === 2 ? "two" : "three"} lineages holding the most resistance`}>
              <ul className="space-y-1">
                {top.map((l) => {
                  const b = blockOfLineage(report, blocks, l.lineage);
                  return (
                    <li key={l.lineage}>
                      <button
                        type="button"
                        disabled={b === undefined}
                        onClick={() => b !== undefined && onPick({ kind: "block", index: b.index })}
                        className={`ui-control-target -mx-1 inline-flex min-h-7 items-center rounded-md px-1 text-left text-control text-ink transition-colors duration-fast ease-out hover:bg-surface-subtle hover:text-blue active:scale-[0.98] ${FOCUS_RING}`}
                      >
                        {`Lineage ${l.lineage}: ${int(l.resistant)} of ${int(l.genomes)} resistant${b !== undefined && !b.named ? `, in block ${b.index + 1}` : ""}`}
                      </button>
                    </li>
                  );
                })}
              </ul>
            </Group>
          </div>
        ) : null}
        <p className="text-note text-muted-text">Click, tap or press Enter on a block or genome for its details; Escape clears it.</p>
      </div>
    );
  }

  return (
    <aside aria-labelledby={titleId} className={`min-w-0 ${frame}`}>
      <div className="flex items-start justify-between gap-3">
        <div className="min-w-0">
          <p className="text-tiny font-medium text-muted-text">{kicker}</p>
          <h3 id={titleId} className="mt-0.5 break-words text-body font-semibold text-ink">
            {title}
          </h3>
        </div>
        {pick !== null ? (
          <GhostButton onClick={() => onPick(null)} ariaLabel="Clear the pinned selection">
            Clear
          </GhostButton>
        ) : null}
      </div>
      <div className="mt-3">{body}</div>
      {/* A pin changes the panel, not the focused grid: say what it now holds. */}
      <p className="sr-only" aria-live="polite">
        {pick === null ? "" : `Inspector: ${title}, ${kicker}.`}
      </p>
    </aside>
  );
};

// ---- Tables -------------------------------------------------------------------------

interface GenomeRowAt {
  pos: number;
  index: number;
}

const CALL_RANK: Record<string, number> = { R: 0, I: 1, S: 2 };

/** Generous pixel widths for 13px body text and an 11.5px header with its
 *  sort glyph and padding; used to fix table columns, never to clip. */
const bodyPx = (s: string): number => Math.ceil(s.length * 8) + 24;
const headPx = (s: string): number => Math.ceil(s.length * 7) + 40;

const StructureTables = ({
  report,
  model,
  order,
  drug,
  pick,
  onPick,
}: {
  report: SpeciesReport;
  model: StructureModel;
  order: readonly number[];
  drug: number;
  pick: StructurePick;
  onPick(p: StructurePick): void;
}) => {
  const { n, blocks } = model;
  const drugName = capitalize(report.drugs[drug] ?? "");
  const totals = drugTotals(report, drug);
  const c = (b: LineageBlock) => model.calls[b.index]![drug]!;
  const blockName = (b: LineageBlock) => `${b.index + 1}. ${blockTitle(b)}`;
  const blockColumns: LinkedColumn<LineageBlock>[] = [
    { id: "block", header: "Block", width: `${Math.max(140, ...blocks.map((b) => labelWidth(blockName(b)))) + 28}px`, sortValue: (b) => b.index, firstDir: "asc", cell: blockName },
    { id: "range", header: "Tree order", align: "right", width: "7.5rem", sortValue: (b) => b.start, firstDir: "asc", cell: (b) => `${int(b.start + 1)}–${int(b.end)}` },
    { id: "size", header: "Genomes", align: "right", width: "6rem", sortValue: blockSize, cell: (b) => int(blockSize(b)) },
    ...STATES.map(
      (k): LinkedColumn<LineageBlock> => ({
        id: k,
        header: k === "untested" ? "Not tested" : capitalize(STATE_WORD[k]),
        align: "right",
        width: k === "i" ? "7rem" : "6.5rem",
        sortValue: (b) => c(b)[k],
        cell: (b) => int(c(b)[k]),
      }),
    ),
  ];
  const rows: GenomeRowAt[] = report.leaf_order.map((index, pos) => ({ pos, index }));
  const mlstText = (r: GenomeRowAt) => {
    const m = report.genomes[r.index]?.mlst ?? null;
    return m === null ? "–" : stLabel(m);
  };
  // Every genome column takes a fixed width from the file's own longest
  // value, so filling in, sorting or pinning never moves a column edge.
  const widths = useMemo(() => {
    const fit = (header: string, values: Iterable<string>, unit?: string) => {
      let w = Math.max(headPx(header), unit === undefined ? 0 : headPx(unit));
      for (const v of values) w = Math.max(w, bodyPx(v));
      return `${w}px`;
    };
    const g = report.genomes;
    return {
      id: fit("Genome", g.map((x) => x.id)),
      pos: fit("Tree order", [int(g.length)]),
      lineage: fit("Lineage", [String(Math.max(0, ...g.map((x) => x.lineage)))]),
      block: fit("Block", [String(model.blocks.length)]),
      mlst: fit("MLST", [...new Set(g.map((x) => (x.mlst === null ? "–" : stLabel(x.mlst))))]),
      country: fit("Country", [...new Set(g.map((x) => x.country ?? "–"))]),
      year: fit("Year", [...new Set(g.map((x) => x.year ?? "–"))]),
      drugs: report.drugs.map((d, j) => fit(capitalize(d), ["–"], callsWords(model, j))),
    };
  }, [report, model]);
  const genomeColumns: LinkedColumn<GenomeRowAt>[] = [
    { id: "id", header: "Genome", width: widths.id, sortValue: (r) => report.genomes[r.index]?.id ?? "", firstDir: "asc", cell: (r) => report.genomes[r.index]?.id ?? "" },
    { id: "pos", header: "Tree order", align: "right", width: widths.pos, sortValue: (r) => r.pos, firstDir: "asc", cell: (r) => int(r.pos + 1) },
    { id: "lineage", header: "Lineage", align: "right", width: widths.lineage, sortValue: (r) => report.genomes[r.index]?.lineage ?? null, firstDir: "asc", cell: (r) => String(report.genomes[r.index]?.lineage ?? "") },
    { id: "block", header: "Block", align: "right", width: widths.block, sortValue: (r) => blockAt(blocks, r.pos)?.index ?? null, firstDir: "asc", cell: (r) => String((blockAt(blocks, r.pos)?.index ?? -1) + 1) },
    { id: "mlst", header: "MLST", width: widths.mlst, sortValue: (r) => report.genomes[r.index]?.mlst ?? null, firstDir: "asc", cell: mlstText },
    { id: "country", header: "Country", width: widths.country, sortValue: (r) => report.genomes[r.index]?.country ?? null, firstDir: "asc", cell: (r) => report.genomes[r.index]?.country ?? "–" },
    { id: "year", header: "Year", align: "right", width: widths.year, sortValue: (r) => report.genomes[r.index]?.year ?? null, firstDir: "asc", cell: (r) => report.genomes[r.index]?.year ?? "–" },
    ...order.map(
      (j): LinkedColumn<GenomeRowAt> => ({
        id: `d${j}`,
        header: capitalize(report.drugs[j] ?? ""),
        unit: callsWords(model, j),
        align: "right",
        width: widths.drugs[j],
        sortValue: (r) => CALL_RANK[report.genomes[r.index]?.profile[j] ?? "."] ?? null,
        firstDir: "asc",
        cell: (r) => CALL_LETTER[stateOf(report.genomes[r.index]?.profile[j])],
      }),
    ),
  ];
  return (
    <div className="space-y-8">
      <div>
        <h3 className="text-control font-semibold text-ink">{`${drugName} calls by lineage block`}</h3>
        <p className="mt-0.5 text-note text-muted-text">{`${plural(blocks.length, "block")} in tree order; the totals row matches the map drug's counts. Choose another map drug above for its blocks.`}</p>
        <div className="mt-3">
          <LinkedTable
            label={`${drugName} calls by lineage block`}
            rows={blocks}
            columns={blockColumns}
            rowKey={(b) => `b${b.index}`}
            selectedKey={pick?.kind === "block" ? `b${pick.index}` : null}
            onSelect={(b) => onPick(b === null ? null : { kind: "block", index: b.index })}
            initialSort={{ id: "block", dir: "asc" }}
            maxHeight="22rem"
            foot={{ block: "All genomes", range: `1–${int(n)}`, size: int(n), r: int(totals.r), i: int(totals.i), s: int(totals.s), untested: int(totals.untested) }}
          />
        </div>
      </div>
      <div>
        <h3 className="text-control font-semibold text-ink">Every genome</h3>
        <p className="mt-0.5 text-note text-muted-text">{`${plural(n, "genome")} with lineage, block, MLST, country, year and one call per drug: R resistant, I intermediate, S susceptible, – not tested.`}</p>
        <div className="mt-3">
          <LinkedTable
            label="Every genome's calls"
            rows={rows}
            columns={genomeColumns}
            rowKey={(r) => `g${r.pos}`}
            selectedKey={pick?.kind === "genome" ? `g${pick.pos}` : null}
            onSelect={(r) => onPick(r === null ? null : { kind: "genome", pos: r.pos })}
            initialSort={{ id: "pos", dir: "asc" }}
            maxHeight="28rem"
            stickyFirst
          />
        </div>
      </div>
    </div>
  );
};
