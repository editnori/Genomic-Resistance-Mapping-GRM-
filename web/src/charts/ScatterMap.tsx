// The species map (PCoA of Mash distances) on canvas, since a species can
// carry thousands of genomes. Both axes share one scale, so distance on the
// map approximates Mash distance. Context genomes are small neutral dots;
// marked genomes carry their series color with a surface ring; flagged
// genomes add an ink ring. The nearest genome within 14px is the hover.
// The map is one tab stop, a listbox over the genomes in frame: arrow keys
// walk them left to right and show the tooltip a hover shows. A new frame
// eases in from the last one, so the marks travel to their new places.

import { useEffect, useId, useRef, useState, type KeyboardEvent } from "react";
import type { TipLine } from "../ui/HoverTip";
import { ChartFrame, ChartTip, prepareCanvas, stepCursor, tipText, useCanvasTokens, useWidth, type CanvasTokens } from "./kit";
import { FOCUS_RING, Segmented } from "../ui/primitives";
import { int } from "../data/format";

export type PointTone = "r" | "s" | "neutral" | "blue";

export interface MapPoint {
  x: number;
  y: number;
  /** Context genomes draw first, small and neutral. */
  layer: "context" | "mark";
  tone: PointTone;
  /** Draw an ink ring (an out-of-fold error, a chosen lineage). */
  flag?: boolean;
  tip(): TipLine[];
}

const PAD = 14;

type Bounds = readonly [number, number, number, number];

/** The frame on screen: a new frame (All genomes or Main cluster) eases in
 *  from the last one on the slow clock instead of jumping; reduced motion
 *  snaps to it. */
const useFrameTween = (goal: Bounds): Bounds => {
  const [shown, setShown] = useState<Bounds>(goal);
  const last = useRef<Bounds>(goal);
  const key = goal.join(",");
  useEffect(() => {
    const from = last.current;
    if (from.join(",") === key || window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
      last.current = goal;
      setShown(goal);
      return undefined;
    }
    const start = performance.now();
    let frame = 0;
    const step = (now: number) => {
      const t = Math.min((now - start) / 320, 1);
      const eased = 1 - (1 - t) ** 3;
      const b: Bounds = [0, 1, 2, 3].map((i) => from[i]! + (goal[i]! - from[i]!) * eased) as unknown as Bounds;
      last.current = b;
      setShown(b);
      if (t < 1) frame = requestAnimationFrame(step);
    };
    frame = requestAnimationFrame(step);
    return () => cancelAnimationFrame(frame);
  }, [key]);
  return shown;
};

const medianOf = (v: readonly number[]): number => {
  const s = [...v].sort((a, b) => a - b);
  const k = s.length >> 1;
  return s.length === 0 ? 0 : s.length % 2 === 1 ? s[k]! : (s[k - 1]! + s[k]!) / 2;
};

/** The frame around the main cluster: the median genome plus or minus five
 *  robust standard deviations (1.4826 times the median absolute deviation)
 *  on each axis. For species whose distant clusters (a species complex)
 *  squeeze every lineage into a few dots. Null when the spread is zero. */
export const coreBounds = (xs: readonly number[], ys: readonly number[]): [number, number, number, number] | null => {
  const window = (v: readonly number[]): [number, number] | null => {
    const m = medianOf(v);
    const mad = medianOf(v.map((x) => Math.abs(x - m))) * 1.4826;
    return mad > 0 ? [m - 5 * mad, m + 5 * mad] : null;
  };
  const wx = window(xs);
  const wy = window(ys);
  return wx === null || wy === null ? null : [wx[0], wx[1], wy[0], wy[1]];
};

const smallShare = new Intl.NumberFormat("en-US", { style: "percent", maximumSignificantDigits: 2 });
const sumOfSquares = (v: readonly number[]): number => {
  const m = v.reduce((a, b) => a + b, 0) / v.length;
  return v.reduce((a, b) => a + (b - m) ** 2, 0);
};

/** At the shared scale, a frame where one axis spans under 5% of the other
 *  (the floor the scale uses) draws as a line. The note names the narrow
 *  axis and gives the second eigenvalue as a share of the first: each PCoA
 *  axis is an eigenvector scaled by the square root of its eigenvalue, so
 *  the axes' sums of squares over all genomes are the two eigenvalues.
 *  Null when both axes of the frame are visible. */
export const flatMapNote = (frameX: readonly number[], frameY: readonly number[], allX: readonly number[], allY: readonly number[]): string | null => {
  if (frameX.length < 3) return null;
  const spanX = Math.max(...frameX) - Math.min(...frameX);
  const spanY = Math.max(...frameY) - Math.min(...frameY);
  const wide = Math.max(spanX, spanY);
  if (!(wide > 0) || Math.min(spanX, spanY) >= 0.05 * wide) return null;
  const [narrow, other] = spanY < spanX ? (["vertical", "horizontal"] as const) : (["horizontal", "vertical"] as const);
  return `The coordinates themselves form a line: drawn to one scale, this frame's ${narrow} span is ${smallShare.format(Math.min(spanX, spanY) / wide)} of its ${other} span. Over all genomes, the second principal coordinate's eigenvalue, the sum of its squared coordinates, is ${smallShare.format(sumOfSquares(allY) / sumOfSquares(allX))} of the first's.`;
};

export const ScatterMap = ({ points, height = 380, label }: { points: readonly MapPoint[]; height?: number; label: string }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const canvas = useRef<HTMLCanvasElement | null>(null);
  const optionId = useId();
  const width = useWidth(frame);
  const tokens = useCanvasTokens();
  // The active genome, from the pointer or the keyboard. `tip` false keeps
  // the keyboard's place while Escape or blur hides the tooltip.
  const [hover, setHover] = useState<{ i: number; tip: boolean } | null>(null);
  const [frameKind, setFrameKind] = useState<"all" | "core">("all");

  const xs = points.map((p) => p.x);
  const ys = points.map((p) => p.y);
  const core = coreBounds(xs, ys);
  const goal: Bounds = frameKind === "core" && core !== null ? core : [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
  // The chosen frame decides what is in it; the eased frame decides what is
  // drawn while it travels there.
  const inFrame = (p: MapPoint) => p.x >= goal[0] && p.x <= goal[1] && p.y >= goal[2] && p.y <= goal[3];
  const [x0, x1, y0, y1] = useFrameTween(goal);
  const inView = (p: MapPoint) => p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1;
  const outside = points.filter((p) => !inFrame(p)).length;
  const shown = points.filter(inFrame);
  const flat = flatMapNote(shown.map((p) => p.x), shown.map((p) => p.y), xs, ys);
  const span = Math.max(x1 - x0, y1 - y0, 1e-9);
  const k = Math.min((width - 2 * PAD) / Math.max(x1 - x0, span * 0.05), (height - 2 * PAD) / Math.max(y1 - y0, span * 0.05));
  const ox = (width - (x1 - x0) * k) / 2;
  const oy = (height - (y1 - y0) * k) / 2;
  const px = (x: number) => ox + (x - x0) * k;
  const py = (y: number) => height - (oy + (y - y0) * k);

  useEffect(() => {
    const el = canvas.current;
    if (el === null || points.length === 0) return;
    const ctx = prepareCanvas(el, width, height);
    if (ctx === null) return;
    const color = (t: PointTone, tk: CanvasTokens) => (t === "r" ? tk.r : t === "s" ? tk.s : t === "blue" ? tk.blue : tk.neutral);
    ctx.clearRect(0, 0, width, height);
    ctx.globalAlpha = 0.45;
    points.forEach((p) => {
      if (p.layer !== "context" || !inView(p)) return;
      ctx.fillStyle = tokens.neutral;
      ctx.beginPath();
      ctx.arc(px(p.x), py(p.y), 2, 0, Math.PI * 2);
      ctx.fill();
    });
    ctx.globalAlpha = 1;
    points.forEach((p) => {
      if (p.layer !== "mark" || !inView(p)) return;
      ctx.beginPath();
      ctx.arc(px(p.x), py(p.y), 4, 0, Math.PI * 2);
      ctx.fillStyle = color(p.tone, tokens);
      ctx.strokeStyle = tokens.surface;
      ctx.lineWidth = 1.5;
      ctx.fill();
      ctx.stroke();
    });
    points.forEach((p) => {
      if (p.flag !== true || !inView(p)) return;
      ctx.beginPath();
      ctx.arc(px(p.x), py(p.y), 7, 0, Math.PI * 2);
      ctx.strokeStyle = tokens.ink;
      ctx.lineWidth = 1.5;
      ctx.stroke();
    });
  }, [points, width, height, tokens, x0, x1, y0, y1]);

  const onMove = (clientX: number, clientY: number) => {
    const r = canvas.current?.getBoundingClientRect();
    if (r === undefined) return;
    const mx = clientX - r.left;
    const my = clientY - r.top;
    let best = -1;
    let bestD = 14 * 14;
    points.forEach((p, i) => {
      if (!inFrame(p)) return;
      const d = (px(p.x) - mx) ** 2 + (py(p.y) - my) ** 2;
      // Marked genomes win ties over the context behind them.
      if (d < bestD || (d === bestD && p.layer === "mark")) {
        best = i;
        bestD = d;
      }
    });
    setHover(best >= 0 ? { i: best, tip: true } : null);
  };

  // Keyboard order: genomes in the frame, left to right, then top to bottom.
  const order = points.flatMap((p, i) => (inFrame(p) ? [i] : [])).sort((a, b) => points[a]!.x - points[b]!.x || points[b]!.y - points[a]!.y);
  const at = hover === null ? -1 : order.indexOf(hover.i);
  const cursorIndex = order[Math.max(at, 0)];
  const cursor = cursorIndex === undefined ? undefined : points[cursorIndex];

  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key === "Escape" && hover?.tip === true) {
      // Escape closes the tooltip and nothing behind it, as HoverTip does.
      e.preventDefault();
      setHover({ i: hover.i, tip: false });
      return;
    }
    const next = stepCursor(e.key, at, order.length);
    if (next === null) return;
    e.preventDefault();
    setHover({ i: order[next]!, tip: true });
  };

  const h = hover?.tip === true ? points[hover.i] : undefined;
  // The frame toggle sits above the chart frame, so the ring and tooltip
  // share the canvas's coordinates.
  return (
    <>
      <div className="mb-2 flex flex-wrap items-center justify-between gap-2">
        <Segmented<"all" | "core">
          ariaLabel="Map frame"
          options={[
            { id: "all", label: "All genomes" },
            { id: "core", label: "Main cluster", disabled: core === null },
          ]}
          value={frameKind}
          onSelect={setFrameKind}
        />
        {outside > 0 ? <span className="text-control text-muted-text">{int(outside)} genomes lie outside this frame</span> : null}
      </div>
      {flat !== null ? <p className="mb-2 text-control text-muted-text">{flat}</p> : null}
      <ChartFrame frameRef={frame}>
        <div
          role="listbox"
          tabIndex={0}
          aria-label={label}
          aria-activedescendant={cursor !== undefined ? optionId : undefined}
          onKeyDown={onKey}
          onFocus={(e) => {
            // Keyboard focus opens the tooltip where the cursor rests; a click
            // leaves it to the pointer.
            if (cursorIndex !== undefined && e.currentTarget.matches(":focus-visible")) setHover({ i: cursorIndex, tip: true });
          }}
          onBlur={() => setHover((v) => (v === null ? null : { i: v.i, tip: false }))}
          className={`rounded-lg ${FOCUS_RING}`}
        >
          <canvas
            ref={canvas}
            aria-hidden="true"
            className="block rounded-lg ring-1 ring-hair"
            onPointerMove={(e) => onMove(e.clientX, e.clientY)}
            onPointerLeave={() => setHover(null)}
          />
          {cursor !== undefined ? (
            <div id={optionId} role="option" aria-selected="true" aria-posinset={Math.max(at, 0) + 1} aria-setsize={order.length} className="sr-only">
              {tipText(cursor.tip())}
            </div>
          ) : null}
        </div>
        {h !== undefined ? (
          <>
            <div aria-hidden="true" className="pointer-events-none absolute rounded-full border-2 border-ink" style={{ left: px(h.x) - 8, top: py(h.y) - 8, width: 16, height: 16 }} />
            <ChartTip x={px(h.x)} y={py(h.y)} width={width} lines={h.tip()} />
          </>
        ) : null}
      </ChartFrame>
    </>
  );
};
