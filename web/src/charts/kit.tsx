// Chart kit: scales, ticks, hairline axes, reference lines and direct
// labels, the chart tooltip, the keyboard listbox over a chart's marks, a
// size observer, and theme tokens resolved for canvas. Marks follow the
// dataviz specs: 2px lines, markers of at least 8px with a 2px surface ring,
// solid hairline grids, tick and label text at the 12px rank.

import { useEffect, useId, useLayoutEffect, useRef, useState, type KeyboardEvent, type ReactNode, type RefObject } from "react";
import { TIP_SURFACE, TipBody, type TipLine } from "../ui/HoverTip";
import { useTheme } from "../theme/ThemeProvider";
import { FOCUS_RING } from "../ui/primitives";

export interface Linear {
  (v: number): number;
  domain: [number, number];
  range: [number, number];
  invert(px: number): number;
}

export const linear = (domain: [number, number], range: [number, number]): Linear => {
  const [d0, d1] = domain;
  const [r0, r1] = range;
  const span = d1 - d0 === 0 ? 1 : d1 - d0;
  const f = ((v: number) => r0 + ((v - d0) / span) * (r1 - r0)) as Linear;
  f.domain = domain;
  f.range = range;
  f.invert = (px: number) => d0 + ((px - r0) / (r1 - r0 === 0 ? 1 : r1 - r0)) * span;
  return f;
};

/** Round tick values covering [lo, hi], about `count` of them. */
export const ticks = (lo: number, hi: number, count = 5): number[] => {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / count;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) ?? raw;
  const out: number[] = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.round(v / step) * step);
  return out;
};

export const TICK_CLASS = "fill-[color:var(--chart-tick)] text-[12px]";
/** Reference, region, and direct labels: ink at the 12px rank. */
export const LABEL_CLASS = "fill-[rgb(var(--ink))] text-[12px]";

/** About how wide a 12px Inter label draws; used to keep labels apart and
 *  inside the plot, never to shrink them. */
export const labelWidth = (text: string): number => Math.ceil(text.length * 6.7) + 4;

export interface LabelSpec {
  key: string;
  text: string;
  /** Where the label points: a reference line's x, or a region's edge. */
  x: number;
  /** Preferred side of x; the placer flips it when that side has no room. */
  anchor: "start" | "end";
}

export interface PlacedLabel extends LabelSpec {
  row: number;
  /** The text's x after the 4px offset from the line. */
  tx: number;
}

/** Greedy placement in rows above a plot, in priority order: each label
 *  takes the first row where it overlaps no earlier label and stays inside
 *  [0, width], trying its preferred side first. A label with no room in any
 *  row is left out (and listed in the caller's Data and method). */
export const placeLabels = (specs: readonly LabelSpec[], width: number, rows: number): PlacedLabel[] => {
  const placed: (PlacedLabel & { left: number; right: number })[] = [];
  specs.forEach((s) => {
    const w = labelWidth(s.text);
    for (let row = 0; row < rows; row += 1) {
      for (const anchor of s.anchor === "start" ? (["start", "end"] as const) : (["end", "start"] as const)) {
        const tx = anchor === "start" ? s.x + 4 : s.x - 4;
        const left = anchor === "start" ? tx : tx - w;
        const right = left + w;
        if (left < 0 || right > width) continue;
        if (placed.some((p) => p.row === row && left < p.right + 6 && right > p.left - 6)) continue;
        placed.push({ ...s, anchor, row, tx, left, right });
        return;
      }
    }
  });
  return placed.map(({ left: _l, right: _r, ...p }) => p);
};

/** A reference line with its label beside it; dashed, because it marks a
 *  stated value, not a grid step. */
export const RefLine = ({ x, y0, y1, dashed = true }: { x: number; y0: number; y1: number; dashed?: boolean }) => (
  <line aria-hidden="true" x1={x} x2={x} y1={y0} y2={y1} stroke="var(--chart-ref)" strokeDasharray={dashed ? "4 3" : undefined} />
);

/** A horizontal reference line across the plot. */
export const RefLineH = ({ y, x0, x1 }: { y: number; x0: number; x1: number }) => (
  <line aria-hidden="true" x1={x0} x2={x1} y1={y} y2={y} stroke="var(--chart-ref)" strokeDasharray="4 3" />
);

/** Placed labels as SVG text, one row every 16px from `top`. */
export const PlacedText = ({ labels, top }: { labels: readonly PlacedLabel[]; top: number }) => (
  <g aria-hidden="true">
    {labels.map((l) => (
      <text key={l.key} x={l.tx} y={top + l.row * 16 + 12} textAnchor={l.anchor} className={LABEL_CLASS}>
        {l.text}
      </text>
    ))}
  </g>
);

/** A label for one named mark, with a hairline leader from the mark. */
export const DirectLabel = ({ x, y, lx, ly, lines, anchor }: { x: number; y: number; lx: number; ly: number; lines: readonly string[]; anchor: "start" | "end" }) => (
  <g aria-hidden="true">
    <line x1={x} y1={y} x2={lx + (anchor === "start" ? -2 : 2)} y2={ly - 4} stroke="var(--chart-axis)" />
    <text x={lx} y={ly} textAnchor={anchor} className={LABEL_CLASS}>
      {lines.map((line, i) => (
        <tspan key={i} x={lx} dy={i === 0 ? 0 : 15} className={i === 0 ? "font-medium" : ""}>
          {line}
        </tspan>
      ))}
    </text>
  </g>
);

/** A document-unique id usable inside url(#...) references. */
export const useSvgId = (prefix: string): string => `${prefix}-${useId().replace(/[^a-zA-Z0-9_-]/g, "")}`;

/** Horizontal axis: a baseline hairline, vertical gridlines, labels below. */
export const AxisX = ({ x, y, height, values, format }: { x: Linear; y: number; height: number; values: number[]; format(v: number): string }) => (
  <g aria-hidden="true">
    {values.map((v) => (
      <g key={v}>
        <line x1={x(v)} x2={x(v)} y1={y - height} y2={y} stroke="var(--chart-grid)" />
        <text x={x(v)} y={y + 15} textAnchor="middle" className={TICK_CLASS}>
          {format(v)}
        </text>
      </g>
    ))}
    <line x1={x.range[0]} x2={x.range[1]} y1={y} y2={y} stroke="var(--chart-axis)" />
  </g>
);

/** Vertical axis: horizontal gridlines with labels on the left. */
export const AxisY = ({ y, x0, x1, values, format }: { y: Linear; x0: number; x1: number; values: number[]; format(v: number): string }) => (
  <g aria-hidden="true">
    {values.map((v) => (
      <g key={v}>
        <line x1={x0} x2={x1} y1={y(v)} y2={y(v)} stroke="var(--chart-grid)" />
        <text x={x0 - 6} y={y(v) + 4} textAnchor="end" className={TICK_CLASS}>
          {format(v)}
        </text>
      </g>
    ))}
  </g>
);

/** Width of an element, tracked. Starts at `fallback` before layout. */
export const useWidth = (ref: RefObject<HTMLElement | null>, fallback = 640): number => {
  const [w, setW] = useState(fallback);
  useLayoutEffect(() => {
    const el = ref.current;
    if (el === null) return undefined;
    const read = () => { if (el.clientWidth > 0) setW(el.clientWidth); };
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, [ref]);
  return w;
};

/** A tooltip inside a chart frame, placed beside the hovered mark and kept
 *  inside the frame. Values lead, labels follow. */
export const ChartTip = ({ x, y, width, lines }: { x: number; y: number; width: number; lines: readonly TipLine[] }) => {
  const ref = useRef<HTMLDivElement | null>(null);
  const [w, setW] = useState(220);
  useLayoutEffect(() => {
    if (ref.current !== null) setW(ref.current.offsetWidth);
  });
  const left = x + 14 + w > width ? Math.max(0, x - 14 - w) : x + 14;
  return (
    <div ref={ref} role="tooltip" className={`absolute z-20 max-w-[280px] ${TIP_SURFACE}`} style={{ left, top: Math.max(0, y - 12) }}>
      <TipBody lines={lines} />
    </div>
  );
};

/** A tooltip's lines as one string, for the element assistive technology
 *  reads when the keyboard reaches a mark: the same values a hover shows. */
export const tipText = (lines: readonly TipLine[]): string =>
  lines.map((l) => [l.term, l.text, l.source].filter((s) => s !== undefined).join(", ")).join(" ");

/** The position a key moves a chart's keyboard cursor to among `count` marks
 *  (`at` is -1 before the first move): arrows step one, PageUp and PageDown a
 *  tenth of the marks, since a species carries thousands of genomes, and Home
 *  and End go to either end. Null for any other key, so Tab leaves the chart. */
export const stepCursor = (key: string, at: number, count: number): number | null => {
  const page = Math.max(1, Math.round(count / 10));
  const step =
    key === "ArrowLeft" || key === "ArrowUp"
      ? -1
      : key === "ArrowRight" || key === "ArrowDown"
        ? 1
        : key === "PageUp"
          ? -page
          : key === "PageDown"
            ? page
            : key === "Home"
              ? -count
              : key === "End"
                ? count
                : null;
  return step === null || count === 0 ? null : Math.min(count - 1, Math.max(0, at + step));
};

/** One mark of a keyboard-walkable chart, in plot pixels. */
export interface ListMark {
  key: string;
  x: number;
  y: number;
  tip(): TipLine[];
}

/** The keyboard and pointer layer over an SVG chart's marks, the ScatterMap
 *  pattern: one tab stop, a listbox whose active descendant carries the same
 *  text the tooltip shows. Arrows step one mark, PageUp and PageDown a tenth,
 *  Home and End go to either end, Enter opens the mark, and Escape hides only
 *  the tooltip, so a second Escape reaches the dialog behind it. The pointer
 *  takes the nearest mark within 16px; a click opens it. `marks` come in
 *  keyboard order. */
export const MarkListbox = ({
  label,
  marks,
  width,
  onOpen,
  children,
}: {
  label: string;
  marks: readonly ListMark[];
  width: number;
  onOpen?(index: number): void;
  children: ReactNode;
}) => {
  const optionId = useId();
  // `tip` false keeps the keyboard's place while Escape or blur hides it.
  const [hover, setHover] = useState<{ i: number; tip: boolean } | null>(null);
  const at = hover === null || hover.i >= marks.length ? -1 : hover.i;
  const cursor = marks[Math.max(at, 0)];

  const nearest = (px: number, py: number): number => {
    let best = -1;
    let bestD = 16 * 16;
    marks.forEach((m, i) => {
      const d = (m.x - px) ** 2 + (m.y - py) ** 2;
      if (d < bestD) {
        best = i;
        bestD = d;
      }
    });
    return best;
  };

  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key === "Escape" && hover?.tip === true) {
      e.preventDefault();
      setHover({ i: hover.i, tip: false });
      return;
    }
    if (e.key === "Enter" && onOpen !== undefined && at >= 0) {
      e.preventDefault();
      onOpen(at);
      return;
    }
    const next = stepCursor(e.key, at, marks.length);
    if (next === null) return;
    e.preventDefault();
    setHover({ i: next, tip: true });
  };

  const shown = hover?.tip === true ? marks[hover.i] : undefined;
  return (
    <div className="relative">
      <div
        role="listbox"
        tabIndex={0}
        aria-label={label}
        aria-activedescendant={cursor !== undefined ? optionId : undefined}
        onKeyDown={onKey}
        onFocus={(e) => {
          if (cursor !== undefined && e.currentTarget.matches(":focus-visible")) setHover({ i: Math.max(at, 0), tip: true });
        }}
        onBlur={() => setHover((v) => (v === null ? null : { i: v.i, tip: false }))}
        onPointerMove={(e) => {
          const r = e.currentTarget.getBoundingClientRect();
          const i = nearest(e.clientX - r.left, e.clientY - r.top);
          setHover(i >= 0 ? { i, tip: true } : null);
        }}
        onPointerLeave={() => setHover(null)}
        onClick={() => {
          if (onOpen !== undefined && hover !== null && hover.tip) onOpen(hover.i);
        }}
        className={`rounded-md ${shown !== undefined && onOpen !== undefined ? "cursor-pointer" : ""} ${FOCUS_RING}`}
      >
        {children}
        {cursor !== undefined ? (
          <div id={optionId} role="option" aria-selected="true" aria-posinset={Math.max(at, 0) + 1} aria-setsize={marks.length} className="sr-only">
            {tipText(cursor.tip())}
          </div>
        ) : null}
      </div>
      {shown !== undefined ? (
        <>
          <div aria-hidden="true" className="pointer-events-none absolute rounded-full border-2 border-ink" style={{ left: shown.x - 8, top: shown.y - 8, width: 16, height: 16 }} />
          <ChartTip x={shown.x} y={shown.y} width={width} lines={shown.tip()} />
        </>
      ) : null}
    </div>
  );
};

/** A chart frame: relative container that owns the tooltip layer. */
export const ChartFrame = ({ children, className, frameRef }: { children: ReactNode; className?: string; frameRef?: RefObject<HTMLDivElement | null> }) => (
  <div ref={frameRef} className={`relative ${className ?? ""}`}>
    {children}
  </div>
);

export interface CanvasTokens {
  r: string;
  s: string;
  neutral: string;
  axis: string;
  ink: string;
  muted: string;
  surface: string;
  blue: string;
  fill: string;
}

/** Theme tokens resolved to colors for canvas drawing; re-read when the
 *  theme flips. */
export const useCanvasTokens = (): CanvasTokens => {
  const { theme } = useTheme();
  const read = (): CanvasTokens => {
    const css = getComputedStyle(document.documentElement);
    const v = (name: string) => css.getPropertyValue(name).trim();
    const channel = (name: string) => `rgb(${v(name)})`;
    return {
      r: v("--series-r"),
      s: v("--series-s"),
      neutral: v("--series-neutral"),
      axis: v("--chart-axis"),
      ink: channel("--ink"),
      muted: channel("--muted-text"),
      surface: channel("--paper-pure"),
      blue: channel("--blue"),
      fill: v("--fill-2"),
    };
  };
  const [tokens, setTokens] = useState<CanvasTokens>(read);
  useEffect(() => {
    // The class flips in ThemeProvider's effect; read on the next frame.
    const frame = requestAnimationFrame(() => setTokens(read()));
    return () => cancelAnimationFrame(frame);
  }, [theme]);
  return tokens;
};

/** Sizes a canvas for the device pixel ratio and returns its 2D context. */
export const prepareCanvas = (canvas: HTMLCanvasElement, width: number, height: number): CanvasRenderingContext2D | null => {
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;
  const ctx = canvas.getContext("2d");
  ctx?.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
};

