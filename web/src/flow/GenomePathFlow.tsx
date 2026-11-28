// The genome's path: nine pipeline steps, nine comparisons at the step each
// one changed, and what each led to. Three-column interactive flowchart with
// measured SVG orthogonal connectors and a synced detail panel below.

import { useCallback, useLayoutEffect, useMemo, useRef, useState } from "react";
import { FOCUS_RING } from "../ui/primitives";
import {
  FLOW_DATA,
  buildNodeMap,
  fillTemplate,
  uiText,
  useReducedMotion,
  type FlowData,
  type FlowNode,
} from "./model";
import { FlowGlyph, FlowHeader, WalkControls } from "./shared";

type PlacementKind = "step" | "cmp" | "nxt";

interface Placement {
  kind: PlacementKind;
  row: number;
  col: number;
}

interface SvgLine {
  key: string;
  x1: number;
  y1: number;
  x2: number;
  y2: number;
  variant: "spine" | "live" | "muted";
}

interface SvgPath {
  key: string;
  d: string;
  variant: "live" | "muted";
}

interface SvgArrow {
  key: string;
  points: string;
  live: boolean;
}

const arrowPoints = (x: number, y: number, dir: "r" | "d"): string =>
  dir === "r" ? `${x},${y} ${x - 7},${y - 4} ${x - 7},${y + 4}` : `${x},${y} ${x - 4},${y - 7} ${x + 4},${y - 7}`;

const StatusTag = ({ node, calm }: { node: FlowNode; calm: boolean }) => {
  if (node.tag === undefined || node.tag === "") return null;
  const showDot = (node.state === "running" || node.state === "waiting") && node.kind !== "experiment";
  const isCandidate = node.tag === "candidate";
  return (
    <span
      className={`inline-flex flex-none items-center gap-1.5 whitespace-nowrap rounded-full border px-2 py-1 text-tiny leading-none ${
        isCandidate ? "border-blue font-semibold text-ink" : "border-ink/15 text-muted-text"
      }`}
    >
      {showDot ? (
        <i
          aria-hidden="true"
          className={`inline-block h-[7px] w-[7px] rounded-full ${
            node.state === "running"
              ? `bg-blue ${calm ? "" : "flow-pulse"}`
              : "border-[1.5px] border-muted-text bg-transparent"
          }`}
        />
      ) : null}
      {node.tag}
    </span>
  );
};

export interface GenomePathFlowProps {
  flow?: FlowData;
  startNode?: string;
  calm?: boolean;
}

export const GenomePathFlow = ({ flow = FLOW_DATA, startNode, calm: calmProp }: GenomePathFlowProps) => {
  const reducedMotion = useReducedMotion();
  const calm = calmProp ?? reducedMotion;
  const view = flow.views.genome_path;
  const ui = view.ui;
  const nodes = useMemo(() => buildNodeMap(flow), [flow]);

  const liveIds = useMemo(() => {
    const s = new Set<string>();
    for (const e of flow.edges) {
      if (e.live) {
        s.add(e.from);
        s.add(e.to);
      }
    }
    return s;
  }, [flow.edges]);

  const { placements, stepIds } = useMemo(() => {
    const map: Record<string, Placement> = {};
    const steps: string[] = [];
    view.rows.forEach((r, rowIdx) => {
      if (r[0] !== null) {
        map[r[0]] = { kind: "step", row: rowIdx, col: 1 };
        steps.push(r[0]);
      }
      if (r[1] !== null) {
        map[r[1]] = { kind: "cmp", row: rowIdx, col: 3 };
      }
    });
    for (const [id, rowIdx] of Object.entries(view.next)) {
      map[id] = { kind: "nxt", row: rowIdx, col: 5 };
    }
    return { placements: map, stepIds: steps };
  }, [view.rows, view.next]);

  const initialId = startNode ?? view.start;
  const [selectedId, setSelectedId] = useState<string>(view.walk.includes(initialId) ? initialId : (view.walk[0] ?? ""));
  const [narrow, setNarrow] = useState(false);
  const [edgeGeom, setEdgeGeom] = useState<{ lines: SvgLine[]; paths: SvgPath[]; arrows: SvgArrow[] }>({
    lines: [],
    paths: [],
    arrows: [],
  });

  const gridRef = useRef<HTMLDivElement | null>(null);
  const nodeRefs = useRef<Record<string, HTMLButtonElement | null>>({});

  const computeEdges = useCallback(() => {
    const grid = gridRef.current;
    if (grid === null) return;
    const isNarrow = grid.clientWidth > 0 && grid.clientWidth <= 560;
    setNarrow(isNarrow);
    if (isNarrow || grid.offsetWidth === 0) {
      setEdgeGeom({ lines: [], paths: [], arrows: [] });
      return;
    }

    const base = grid.getBoundingClientRect();
    const z = base.width / grid.offsetWidth || 1;
    const rel = (el: Element) => {
      const r = el.getBoundingClientRect();
      const l = (r.left - base.left) / z;
      const t = (r.top - base.top) / z;
      const w = r.width / z;
      const hh = r.height / z;
      return { l, r: l + w, t, b: t + hh, cx: l + w / 2, cy: t + hh / 2 };
    };

    const lines: SvgLine[] = [];
    const paths: SvgPath[] = [];
    const arrows: SvgArrow[] = [];

    for (let i = 0; i + 1 < stepIds.length; i++) {
      const aEl = nodeRefs.current[stepIds[i]!]?.querySelector("[data-glyph]");
      const bEl = nodeRefs.current[stepIds[i + 1]!]?.querySelector("[data-glyph]");
      if (aEl == null || bEl == null) continue;
      const a = rel(aEl);
      const b = rel(bEl);
      lines.push({
        key: `spine-${stepIds[i]}-${stepIds[i + 1]}`,
        x1: a.cx,
        y1: a.b + 2,
        x2: b.cx,
        y2: b.t - 2,
        variant: "spine",
      });
    }

    flow.edges.forEach((e, idx) => {
      if (e.draw === false) return;
      const fromEl = nodeRefs.current[e.from];
      const toEl = nodeRefs.current[e.to];
      const fromPlace = placements[e.from];
      const toPlace = placements[e.to];
      if (fromEl == null || toEl == null || fromPlace === undefined || toPlace === undefined) return;
      const a = rel(fromEl);
      const b = rel(toEl);
      const variant = e.live ? "live" : "muted";

      if (fromPlace.kind === "step") {
        const y = b.cy;
        const x0 = a.r + 2;
        const x1 = b.l - 1;
        const mid = (x0 + x1) / 2;
        if (fromPlace.row === toPlace.row) {
          lines.push({ key: `e-${idx}`, x1: x0, y1: y, x2: x1, y2: y, variant: "muted" });
        } else {
          paths.push({
            key: `e-${idx}`,
            d: `M${mid} ${a.cy}V${y - 6}Q${mid} ${y} ${mid + 6} ${y}H${x1}`,
            variant: "muted",
          });
        }
      } else if (fromPlace.col === toPlace.col) {
        lines.push({ key: `e-${idx}`, x1: a.cx, y1: a.b + 1, x2: a.cx, y2: b.t - 8, variant });
        arrows.push({ key: `a-${idx}`, points: arrowPoints(a.cx, b.t - 1, "d"), live: e.live === true });
      } else {
        const y = Math.max(a.t + 12, Math.min(a.b - 12, b.cy));
        lines.push({ key: `e-${idx}`, x1: a.r + 1, y1: y, x2: b.l - 8, y2: y, variant });
        arrows.push({ key: `a-${idx}`, points: arrowPoints(b.l - 1, y, "r"), live: e.live === true });
      }
    });

    setEdgeGeom({ lines, paths, arrows });
  }, [flow.edges, placements, stepIds]);

  useLayoutEffect(() => {
    computeEdges();
    const grid = gridRef.current;
    if (grid === null) return undefined;
    const ro = typeof ResizeObserver === "undefined" ? null : new ResizeObserver(computeEdges);
    ro?.observe(grid);
    if (typeof document !== "undefined" && "fonts" in document && document.fonts?.ready !== undefined) {
      void document.fonts.ready.then(computeEdges);
    }
    return () => ro?.disconnect();
  }, [computeEdges]);

  const move = useCallback(
    (delta: number) => {
      const w = view.walk;
      const curIdx = Math.max(0, w.indexOf(selectedId));
      const nextId = w[(curIdx + delta + w.length) % w.length]!;
      setSelectedId(nextId);
      nodeRefs.current[nextId]?.focus({ preventScroll: false });
    },
    [selectedId, view.walk],
  );

  const selectedNode = nodes[selectedId];
  const selectedRecords = selectedNode !== undefined ? [...new Set(selectedNode.src.map((x) => `${x.file}:${x.line}`))] : [];

  return (
    <section data-flow-diagram="genome-path" className="mx-auto w-full max-w-[760px] text-ink">
      <FlowHeader flow={flow} kicker={view.kicker} title={view.title} sub={view.sub} />

      <div className="relative">
        <div
          ref={gridRef}
          onKeyDown={(e) => {
            if (e.key === "ArrowRight" || e.key === "ArrowDown") {
              e.preventDefault();
              move(1);
            } else if (e.key === "ArrowLeft" || e.key === "ArrowUp") {
              e.preventDefault();
              move(-1);
            }
          }}
          className={
            narrow
              ? "relative flex flex-col gap-1.5"
              : "relative grid grid-cols-[minmax(0,1fr)_22px_minmax(0,1.34fr)_34px_minmax(0,1.2fr)] grid-rows-[16px] [grid-auto-rows:minmax(48px,auto)] gap-y-3"
          }
        >
          {!narrow ? (
            <>
              <div className="self-start pl-0.5 text-tiny font-semibold text-muted-text" style={{ gridRow: 1, gridColumn: 1 }}>
                {uiText(ui, "col1")}
              </div>
              <div className="self-start pl-0.5 text-tiny font-semibold text-muted-text" style={{ gridRow: 1, gridColumn: 3 }}>
                {uiText(ui, "col2")}
              </div>
              <div className="self-start pl-0.5 text-tiny font-semibold text-muted-text" style={{ gridRow: 1, gridColumn: 5 }}>
                {uiText(ui, "col3")}
              </div>
            </>
          ) : null}

          {view.walk.map((id) => {
            const node = nodes[id];
            const pl = placements[id];
            if (node === undefined || pl === undefined) return null;
            const isSelected = id === selectedId;
            const isStep = pl.kind === "step";
            const isLive = liveIds.has(id) && !isStep;
            const narrowIndent = narrow ? (pl.kind === "cmp" ? "ml-5" : pl.kind === "nxt" ? "ml-11" : "") : "";

            return (
              <button
                key={id}
                ref={(el) => {
                  nodeRefs.current[id] = el;
                }}
                type="button"
                data-path-node={id}
                aria-pressed={isSelected}
                onClick={() => setSelectedId(id)}
                style={narrow ? undefined : { gridRow: pl.row + 2, gridColumn: pl.col }}
                className={`relative z-[1] grid min-w-0 content-center items-center rounded-[10px] border text-left transition-[border-color,background-color,box-shadow] duration-fast ease-out ${FOCUS_RING} ${narrowIndent} ${
                  isStep
                    ? "grid-cols-[30px_minmax(0,1fr)] gap-x-2.5 gap-y-0.5 py-1.5 pr-2.5 pl-0"
                    : "grid-cols-[minmax(0,1fr)_auto] gap-x-2 gap-y-0.5 px-2.5 py-2"
                } ${
                  isSelected
                    ? "border-ink bg-paper-pure shadow-[0_0_0_1px_rgb(var(--ink))]"
                    : isLive
                      ? "border-blue bg-blue/[0.09] hover:border-ink/60"
                      : isStep
                        ? "border-transparent bg-transparent hover:border-ink/30"
                        : "border-ink/15 bg-paper-pure hover:border-ink/40"
                }`}
              >
                {isStep && node.glyph !== undefined ? (
                  <span
                    data-glyph=""
                    className="row-span-2 grid h-[30px] w-[30px] flex-none place-items-center rounded-full border-[1.5px] border-ink/35 bg-paper-pure text-ink"
                  >
                    <FlowGlyph glyph={node.glyph} />
                  </span>
                ) : null}
                <b className="block text-[13.5px] font-semibold leading-[1.25] text-ink">{fillTemplate(flow, node.title)}</b>
                {!isStep ? <StatusTag node={node} calm={calm} /> : null}
                <small className={`block text-[12.5px] leading-[1.3] text-muted-text ${isStep ? "col-start-2" : "col-span-full"}`}>
                  {fillTemplate(flow, node.line)}
                </small>
              </button>
            );
          })}
        </div>

        {!narrow ? (
          <svg aria-hidden="true" className="pointer-events-none absolute inset-0 h-full w-full overflow-visible">
            {edgeGeom.lines.map((l) => (
              <line
                key={l.key}
                x1={l.x1}
                y1={l.y1}
                x2={l.x2}
                y2={l.y2}
                stroke={l.variant === "live" ? "rgb(var(--blue))" : "rgb(var(--ink) / 0.34)"}
                strokeWidth={l.variant === "muted" ? 1.25 : 2}
              />
            ))}
            {edgeGeom.paths.map((p) => (
              <path
                key={p.key}
                d={p.d}
                fill="none"
                stroke={p.variant === "live" ? "rgb(var(--blue))" : "rgb(var(--ink) / 0.34)"}
                strokeWidth={p.variant === "live" ? 2 : 1.25}
              />
            ))}
            {edgeGeom.arrows.map((a) => (
              <polygon key={a.key} points={a.points} fill={a.live ? "rgb(var(--blue))" : "rgb(var(--ink) / 0.34)"} />
            ))}
          </svg>
        ) : null}
      </div>

      {/* Selected node detail panel */}
      {selectedNode !== undefined ? (
        <section aria-live="polite" className="mt-4 min-h-[190px] border-t border-ink/30 pt-3">
          <div className="flex flex-wrap items-center gap-2.5">
            <h3 className="text-[16px] font-semibold leading-[1.3] text-ink">{fillTemplate(flow, selectedNode.title)}</h3>
            <StatusTag node={selectedNode} calm={calm} />
            <WalkControls
              className="ml-auto"
              backLabel={uiText(ui, "prev")}
              nextLabel={uiText(ui, "next")}
              onBack={() => move(-1)}
              onNext={() => move(1)}
            />
          </div>
          <dl className="mt-2.5 grid grid-cols-1 gap-x-3 gap-y-1.5 text-[13.5px] sm:grid-cols-[5.5rem_minmax(0,1fr)]">
            {(selectedNode.fields ?? []).map(([name, text], idx) => (
              <div key={`${name}-${idx}`} className="contents">
                <dt className="pt-0.5 text-[12.5px] font-semibold text-muted-text">{name}</dt>
                <dd className="m-0 max-w-[64ch] text-ink">{fillTemplate(flow, text)}</dd>
              </div>
            ))}
            <dt className="pt-0.5 text-[12.5px] font-semibold text-muted-text">{uiText(ui, "record")}</dt>
            <dd className="m-0 flex flex-wrap gap-x-3.5 gap-y-0.5 font-mono text-tiny text-muted-text">
              {selectedRecords.map((rec) => (
                <span key={rec} className="[overflow-wrap:anywhere]">
                  {rec}
                </span>
              ))}
            </dd>
          </dl>
        </section>
      ) : null}
    </section>
  );
};
