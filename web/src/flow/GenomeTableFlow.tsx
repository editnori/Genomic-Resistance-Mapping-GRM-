// DNA-to-table animated 5-step walk ("A {k}-letter window turns each genome
// into yes and no"). A drawing of the method: made-up DNA letters from a fixed
// seed and a made-up yes/no grid; all counts and copy come from flow.json.

import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import { ACTION_BUTTON } from "../ui/primitives";
import {
  FLOW_DATA,
  buildNodeMap,
  createSeededRng,
  fillTemplate,
  uiText,
  useReducedMotion,
  type BaseLetter,
  type FlowData,
} from "./model";
import { FlowHeader, StepStrip, WalkControls, baseTileStyle } from "./shared";

const DRAWN = [
  "111111",
  "110100",
  "110100",
  "111111",
  "101011",
  "100000",
  "101011",
  "111110",
  "011010",
  "101011",
  "010110",
  "111111",
] as const;

const CALLS = ["R", "S", "R", "R", "S", "S"] as const;
const NC = DRAWN.length;
const NR = CALLS.length;
const BASE_SEED = 20261127;
const BASES: readonly BaseLetter[] = ["A", "C", "G", "T"];

const countOnes = (col: string): number => [...col].filter((ch) => ch === "1").length;

interface StrandTile {
  id: number;
  base: BaseLetter;
}

export interface GenomeTableFlowProps {
  flow?: FlowData;
  startStep?: string;
  calm?: boolean;
}

export const GenomeTableFlow = ({ flow = FLOW_DATA, startStep, calm: calmProp }: GenomeTableFlowProps) => {
  const reducedMotion = useReducedMotion();
  const calm = calmProp ?? reducedMotion;
  const view = flow.views.genome_table;
  const ui = view.ui;
  const nodes = useMemo(() => buildNodeMap(flow), [flow]);
  const steps = view.steps;

  const kLen = flow.counts.k?.v ?? 31;
  const minIn = flow.counts.min_in?.v ?? 2;
  const minOut = flow.counts.min_out?.v ?? 2;

  const keptMask = useMemo(
    () => DRAWN.map((col) => countOnes(col) >= minIn && NR - countOnes(col) >= minOut),
    [minIn, minOut],
  );
  const firstIdx = useMemo(
    () => DRAWN.map((col, i) => (keptMask[i] ? DRAWN.findIndex((d, j) => keptMask[j] && d === col) : -1)),
    [keptMask],
  );

  const initialIdx = Math.max(0, steps.indexOf(startStep ?? view.start));
  const [step, setStep] = useState(initialIdx);
  const [playing, setPlaying] = useState(!calm);
  const [emitted, setEmitted] = useState<number>(initialIdx === 1 && !calm ? 0 : NC);
  const [stageWidth, setStageWidth] = useState(720);
  const [shiftState, setShiftState] = useState<{ dx: number; ms: number; ease: string }>({ dx: 0, ms: 0, ease: "linear" });

  const stageRef = useRef<HTMLDivElement | null>(null);
  const flyRef = useRef<HTMLDivElement | null>(null);
  const rngRef = useRef<() => number>(createSeededRng(BASE_SEED));
  const nextTileIdRef = useRef(0);

  useLayoutEffect(() => {
    const el = stageRef.current;
    if (el === null) return undefined;
    const read = () => {
      if (el.clientWidth > 0) setStageWidth(el.clientWidth);
    };
    read();
    if (typeof ResizeObserver === "undefined") return undefined;
    const ro = new ResizeObserver(read);
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  const geom = useMemo(() => {
    const W = Math.max(300, stageWidth);
    const pitch = Math.max(8, Math.min(15, Math.floor((W - 24) / (kLen + 4))));
    const winX = Math.round((W - kLen * pitch) / 2);
    const LW = W < 440 ? 92 : 118;
    const cp = Math.max(14, Math.min(28, Math.floor((W - LW - 34) / NC)));
    const cell = Math.round(cp * 0.68);
    const T0 = 92;
    const rp = 27;
    const off = (winX % pitch) - pitch;
    const tileCount = Math.ceil((W - off) / pitch) + 2;
    return { W, pitch, winX, LW, cp, cell, T0, rp, off, tileCount };
  }, [stageWidth, kLen]);

  const [tiles, setTiles] = useState<StrandTile[]>(() => {
    const rng = createSeededRng(BASE_SEED);
    rngRef.current = rng;
    nextTileIdRef.current = 0;
    return Array.from({ length: 64 }, () => ({
      id: nextTileIdRef.current++,
      base: BASES[Math.floor(rng() * 4)]!,
    }));
  });

  // Rebuild strand tiles when geometry requires more tiles than currently allocated.
  useEffect(() => {
    setTiles((prev) => {
      if (prev.length >= geom.tileCount) return prev;
      const extra = Array.from({ length: geom.tileCount - prev.length }, () => ({
        id: nextTileIdRef.current++,
        base: BASES[Math.floor(rngRef.current() * 4)]!,
      }));
      return [...prev, ...extra];
    });
  }, [geom.tileCount]);

  // Precompute deterministic mini-strand gradients for rows 1..NR-1.
  const miniGradients = useMemo(() => {
    const rng = createSeededRng(BASE_SEED + 997);
    const width = NC * geom.cp;
    return Array.from({ length: NR - 1 }, () => {
      const stops: string[] = [];
      let x = 0;
      while (x < width) {
        const ch = ["a", "c", "g", "t"][Math.floor(rng() * 4)]!;
        stops.push(`var(--base-${ch}) ${x}px ${x + 5}px`, `transparent ${x + 5}px ${x + 6}px`);
        x += 6;
      }
      return `linear-gradient(90deg, ${stops.join(",")})`;
    });
  }, [geom.cp]);

  const xLeft = useCallback((n: number) => Math.max(geom.LW, Math.round(geom.W / 2 - (n * geom.cp) / 2)), [geom]);
  const cellX = useCallback((idx: number, n: number) => xLeft(n) + idx * geom.cp + (geom.cp - geom.cell) / 2, [geom, xLeft]);

  const go = useCallback(
    (nextIndex: number, keepPlaying?: boolean) => {
      const wrapped = ((nextIndex % steps.length) + steps.length) % steps.length;
      if (keepPlaying !== undefined) setPlaying(keepPlaying);
      setShiftState({ dx: 0, ms: 0, ease: "linear" });
      setEmitted(wrapped === 1 && !calm ? 0 : NC);
      setStep(wrapped);
    },
    [steps.length, calm],
  );

  // Step choreography (scrolling strand, cell emission on step 1, auto-advance when playing).
  useEffect(() => {
    if (calm) {
      setEmitted(NC);
      setShiftState({ dx: 0, ms: 0, ease: "linear" });
      return undefined;
    }

    let cancelled = false;
    const timers: number[] = [];
    const later = (fn: () => void, ms: number) => {
      const id = window.setTimeout(() => {
        if (!cancelled) fn();
      }, ms);
      timers.push(id);
    };

    const doShift = (ms: number, ease: string, done: () => void) => {
      setShiftState({ dx: -geom.pitch, ms, ease });
      later(() => {
        setShiftState({ dx: 0, ms: 0, ease: "linear" });
        setTiles((prev) => {
          const nextBase = BASES[Math.floor(rngRef.current() * 4)]!;
          return [...prev.slice(1), { id: nextTileIdRef.current++, base: nextBase }];
        });
        later(done, 16);
      }, ms + 16);
    };

    if (step === 0) {
      const loop = () => {
        if (!cancelled) doShift(380, "linear", loop);
      };
      loop();
      if (playing) later(() => go(1), 4200);
    } else if (step === 1) {
      let currentEmitted = 0;
      setEmitted(0);
      const emitNext = () => {
        if (cancelled) return;
        if (currentEmitted >= NC) {
          if (playing) later(() => go(2), 1800);
          return;
        }
        const colIdx = currentEmitted;
        doShift(260, "ease-in-out", () => {
          const finishCol = () => {
            if (cancelled) return;
            currentEmitted = colIdx + 1;
            setEmitted(currentEmitted);
            later(emitNext, 150);
          };
          const colBits = DRAWN[colIdx];
          const flyEl = flyRef.current;
          if (colBits?.[0] !== "1" || flyEl === null || typeof flyEl.animate !== "function") {
            finishCol();
            return;
          }
          const x0 = geom.W / 2 - geom.cell / 2;
          const y0 = 54;
          const x1 = cellX(colIdx, NC);
          const y1 = geom.T0;
          const anim = flyEl.animate(
            [
              { transform: `translate(${x0}px,${y0}px) scaleX(${Math.min(12, (kLen * geom.pitch) / geom.cell)})`, opacity: 0.9 },
              { transform: `translate(${x0}px,${y0 + 8}px) scaleX(1)`, opacity: 1, offset: 0.35 },
              { transform: `translate(${x1}px,${y1}px) scaleX(1)`, opacity: 1 },
            ],
            { duration: 420, easing: "ease-in-out" },
          );
          anim.onfinish = finishCol;
        });
      };
      later(emitNext, 500);
    } else if (playing) {
      later(() => go(step + 1), step === 4 ? 5200 : 3600);
    }

    return () => {
      cancelled = true;
      for (const id of timers) window.clearTimeout(id);
      const flyEl = flyRef.current;
      if (flyEl !== null && typeof flyEl.getAnimations === "function") {
        for (const a of flyEl.getAnimations()) a.cancel();
      }
    };
  }, [step, playing, calm, geom.pitch, geom.W, geom.cell, geom.T0, kLen, cellX, go]);

  const activeId = steps[step] ?? steps[0]!;
  const activeNode = nodes[activeId];
  const toolField = activeNode?.fields?.find((f) => f[0] === uiText(ui, "tool"));

  // Compute column visibility and positions for the active step.
  const kept = DRAWN.map((_c, i) => (step < 2 ? true : keptMask[i]!));
  const shown = DRAWN.map((_c, i) => (step < 3 ? kept[i]! : firstIdx[i] === i));
  const order = shown.map((v, i) => (v ? i : -1)).filter((i) => i >= 0);
  const shownCount = order.length;
  const idxOf = (c: number) => order.indexOf(step >= 3 && keptMask[c] ? firstIdx[c]! : c);
  const callX = xLeft(shownCount) + shownCount * geom.cp + 12;

  const stepItems = steps.map((id) => ({
    id,
    label: fillTemplate(flow, nodes[id]?.title ?? id),
  }));

  return (
    <section
      data-flow-diagram="genome-table"
      tabIndex={-1}
      onKeyDown={(e) => {
        if (e.key === "ArrowRight") go(step + 1, false);
        else if (e.key === "ArrowLeft") go(step - 1, false);
      }}
      className="mx-auto w-full max-w-[760px] text-ink outline-none"
    >
      <FlowHeader flow={flow} kicker={view.kicker} title={view.title} sub={view.sub} />

      <StepStrip
        steps={stepItems}
        activeIndex={step}
        onSelect={(idx) => go(idx, false)}
        ariaLabel={fillTemplate(flow, view.kicker)}
      />

      <p className="mb-2 min-h-[40px] max-w-[60ch] text-body text-ink" aria-live="polite">
        {activeNode !== undefined ? fillTemplate(flow, activeNode.line) : ""}
      </p>

      <div ref={stageRef} className="relative h-[286px] overflow-hidden rounded-lg border border-hair bg-paper-pure px-2">
        <div className="pointer-events-none absolute left-2 top-1.5 whitespace-nowrap text-tiny font-semibold text-ink">
          {uiText(ui, "this")}
        </div>
        <div
          className="pointer-events-none absolute top-1.5 whitespace-nowrap text-center text-tiny font-semibold text-ink transition-opacity duration-med"
          style={{ left: `${geom.winX}px`, width: `${kLen * geom.pitch}px`, opacity: step >= 1 ? 1 : 0 }}
        >
          {uiText(ui, "window")}
        </div>

        {/* Scrolling made-up DNA strand */}
        <div
          className="flow-strand-mask absolute left-0 right-0 top-[24px] h-[26px] overflow-hidden transition-opacity duration-slow"
          style={{ opacity: step <= 1 ? 1 : 0.45 }}
          aria-hidden="true"
        >
          <div
            className="absolute top-0 flex h-[26px] will-change-transform"
            style={{
              left: `${geom.off}px`,
              transform: `translateX(${shiftState.dx}px)`,
              transition: shiftState.ms > 0 ? `transform ${shiftState.ms}ms ${shiftState.ease}` : "none",
            }}
          >
            {tiles.slice(0, geom.tileCount).map((tile) => (
              <div
                key={tile.id}
                className="mr-[1px] flex h-[26px] flex-none items-center justify-center rounded-[3px] text-[12px] font-semibold leading-none"
                style={{ width: `${geom.pitch - 1}px`, ...baseTileStyle(tile.base) }}
              >
                {geom.pitch >= 12 ? uiText(ui, tile.base) : ""}
              </div>
            ))}
          </div>
        </div>

        {/* Fixed 31-letter window border */}
        <div
          aria-hidden="true"
          className="pointer-events-none absolute top-[19px] h-[36px] rounded-[7px] border-2 border-ink transition-opacity duration-med"
          style={{
            left: `${geom.winX - 3}px`,
            width: `${kLen * geom.pitch + 5}px`,
            opacity: step >= 1 ? 1 : 0,
          }}
        />

        {/* Mini DNA strands for other genomes in step 0 */}
        {miniGradients.map((bg, rIdx) => (
          <div
            key={rIdx}
            aria-hidden="true"
            className="absolute h-[10px] rounded-[2px] transition-opacity duration-slow"
            style={{
              top: `${geom.T0 + (rIdx + 1) * geom.rp + (geom.cell - 10) / 2}px`,
              left: `${xLeft(NC)}px`,
              width: `${NC * geom.cp}px`,
              background: bg,
              opacity: step === 0 ? 0.6 : 0,
            }}
          />
        ))}

        {/* Yes/no binary matrix */}
        {DRAWN.map((col, c) => {
          const isTargetPos = shown[c] || (step >= 3 && keptMask[c]);
          const left = isTargetPos ? cellX(idxOf(c), shownCount) : cellX(c, NC);
          return Array.from({ length: NR }, (_unused, r) => {
            const isYes = col[r] === "1";
            const visible = step >= 1 && shown[c] && (r > 0 || step > 1 || c < emitted);
            return (
              <div
                key={`${c}-${r}`}
                aria-hidden="true"
                className="absolute rounded-[4px] border-[1.5px]"
                style={{
                  top: `${geom.T0 + r * geom.rp}px`,
                  left: `${left}px`,
                  width: `${geom.cell}px`,
                  height: `${geom.cell}px`,
                  opacity: visible ? 1 : 0,
                  backgroundColor: isYes
                    ? r === 0
                      ? "rgb(var(--blue))"
                      : "color-mix(in srgb, rgb(var(--ink)) 62%, transparent)"
                    : "transparent",
                  borderColor: isYes ? "transparent" : "rgb(var(--ink) / 0.34)",
                  transition: calm
                    ? "none"
                    : "left 700ms cubic-bezier(0.3, 0, 0.2, 1), opacity 500ms, background-color 250ms, border-color 250ms",
                }}
              />
            );
          });
        })}

        {/* Lab's R/S call column (step 4) */}
        <div
          className="pointer-events-none absolute whitespace-nowrap text-tiny text-muted-text transition-opacity duration-slow"
          style={{
            top: `${geom.T0 - 22}px`,
            left: `${Math.min(callX - 4, geom.W - 62)}px`,
            opacity: step === 4 ? 1 : 0,
          }}
        >
          {uiText(ui, "call")}
        </div>
        {CALLS.map((callVal, r) => (
          <div
            key={r}
            aria-hidden="true"
            className={`absolute flex items-center justify-center rounded-full border-[1.5px] border-ink text-[12px] font-semibold leading-none ${
              callVal === "R" ? "bg-ink text-paper-pure" : "bg-transparent text-ink"
            }`}
            style={{
              top: `${geom.T0 + r * geom.rp - 1}px`,
              left: `${callX}px`,
              width: `${geom.cell + 2}px`,
              height: `${geom.cell + 2}px`,
              opacity: step === 4 ? 1 : 0,
              transition: calm ? "none" : "left 700ms cubic-bezier(0.3, 0, 0.2, 1), opacity 500ms",
            }}
          >
            {uiText(ui, callVal)}
          </div>
        ))}

        {/* Row labels */}
        <div
          className="pointer-events-none absolute whitespace-nowrap text-right text-tiny font-semibold text-ink"
          style={{
            top: `${geom.T0 + (geom.cell - 16) / 2}px`,
            left: `${xLeft(shownCount) - geom.LW}px`,
            width: `${geom.LW - 10}px`,
            opacity: step >= 1 ? 1 : 0,
            transition: calm ? "none" : "left 700ms cubic-bezier(0.3, 0, 0.2, 1), opacity 500ms",
          }}
        >
          {uiText(ui, "this")}
        </div>
        <div
          className="pointer-events-none absolute whitespace-nowrap text-right text-tiny text-muted-text"
          style={{
            top: `${geom.T0 + 3 * geom.rp + (geom.cell - 16) / 2}px`,
            left: `${xLeft(shownCount) - geom.LW}px`,
            width: `${geom.LW - 10}px`,
            transition: calm ? "none" : "left 700ms cubic-bezier(0.3, 0, 0.2, 1), opacity 500ms",
          }}
        >
          {uiText(ui, "others")}
        </div>

        {/* Flying emission dot */}
        <div
          ref={flyRef}
          aria-hidden="true"
          className="pointer-events-none absolute left-0 top-0 rounded-[4px] bg-blue opacity-0"
          style={{ width: `${geom.cell}px`, height: `${geom.cell}px` }}
        />

        {/* Legend */}
        <div
          className="absolute left-0 right-0 flex justify-center gap-4 text-tiny text-muted-text transition-opacity duration-slow"
          style={{ top: `${geom.T0 + NR * geom.rp + 8}px`, opacity: step >= 1 ? 1 : 0 }}
        >
          <span className="inline-flex items-center gap-1.5">
            <i
              aria-hidden="true"
              className="inline-block h-3 w-3 rounded-[3px]"
              style={{ backgroundColor: "color-mix(in srgb, rgb(var(--ink)) 62%, transparent)" }}
            />
            {uiText(ui, "yes")}
          </span>
          <span className="inline-flex items-center gap-1.5">
            <i aria-hidden="true" className="inline-block h-3 w-3 rounded-[3px] border-[1.5px] border-ink/35" />
            {uiText(ui, "no")}
          </span>
        </div>
      </div>

      <WalkControls
        className="mt-3"
        backLabel={uiText(ui, "prev")}
        nextLabel={uiText(ui, "next")}
        onBack={() => go(step - 1, false)}
        onNext={() => go(step + 1, false)}
        middle={
          !calm ? (
            <button type="button" onClick={() => setPlaying((p) => !p)} className={ACTION_BUTTON}>
              {playing ? uiText(ui, "pause") : uiText(ui, "play")}
            </button>
          ) : null
        }
        right={
          toolField !== undefined ? (
            <span className="text-control text-muted-text">
              {uiText(ui, "tool")}{" "}
              <strong className="font-mono font-semibold text-ink">{fillTemplate(flow, toolField[1])}</strong>
            </span>
          ) : null
        }
      />
    </section>
  );
};
