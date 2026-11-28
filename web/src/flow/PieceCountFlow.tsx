// The counting walk: seven interactive steps showing how DNA letters become
// yes-or-no patterns (contigs, sliding 31-letter window, canonical reverse-
// complement mirror, 2-bit integer encoding, sorted deduplicated sketch,
// rare-piece filter, and identical-column merge into unique patterns).

import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import {
  FLOW_DATA,
  buildNodeMap,
  fillTemplate,
  makeLetters,
  uiText,
  useReducedMotion,
  type BaseLetter,
  type FlowData,
  type FlowNode,
} from "./model";
import {
  ChoiceRowsList,
  ContigsStrip,
  DnaTiles,
  FlowBars,
  FlowHeader,
  StepStrip,
  WalkControls,
} from "./shared";

const SHOWN = 46;
const AT = 7;
const SEQ: readonly BaseLetter[] = makeLetters(SHOWN, 20261127);
const PAIR: Record<BaseLetter, BaseLetter> = { A: "T", T: "A", C: "G", G: "C" };
const GRID = [
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

const countOnes = (c: string): number => [...c].filter((x) => x === "1").length;
const baseBits = (c: BaseLetter): [boolean, boolean] => {
  const v = "ACGT".indexOf(c);
  return [v > 1, v % 2 === 1];
};

const StageLabel = ({ bold, text }: { bold?: string; text: string }) => (
  <p className="mb-1 text-tiny text-muted-text">
    {bold !== undefined && bold !== "" ? <b className="mr-1 font-semibold text-ink">{bold}</b> : null}
    {text}
  </p>
);

const WindowStrand = ({
  tileSize,
  kLen,
  winAt,
  windowLabel,
  calm,
}: {
  tileSize: number;
  kLen: number;
  winAt: number | null;
  windowLabel: string;
  calm: boolean;
}) => (
  <div className="relative pt-5">
    <DnaTiles seq={SEQ} size={tileSize} withLetters={tileSize >= 12} />
    {winAt !== null ? (
      <div
        className="pointer-events-none absolute top-[17px] rounded-[4px] border-2 border-ink"
        style={{
          height: `${tileSize + 6}px`,
          width: `${kLen * tileSize + 4}px`,
          left: `${winAt * tileSize - 2}px`,
          transition: calm ? "none" : "left 450ms cubic-bezier(0.2, 0.8, 0.2, 1)",
        }}
      >
        <span className="absolute -top-[19px] left-0 whitespace-nowrap text-tiny font-semibold text-ink">{windowLabel}</span>
      </div>
    ) : null}
  </div>
);

const K2Stage = ({ tileSize, kLen, ui, calm }: { tileSize: number; kLen: number; ui: Record<string, string>; calm: boolean }) => {
  const [idx, setIdx] = useState(calm ? 7 : 0);
  useEffect(() => {
    if (calm) {
      setIdx(7);
      return undefined;
    }
    setIdx(0);
    const timer = window.setInterval(() => {
      setIdx((cur) => (cur + 1) % 8);
    }, 700);
    return () => window.clearInterval(timer);
  }, [calm]);

  return (
    <>
      <WindowStrand tileSize={tileSize} kLen={kLen} winAt={idx} windowLabel={uiText(ui, "window")} calm={calm} />
      <StageLabel text={uiText(ui, "pieces")} />
      <div className="relative h-[78px]" aria-hidden="true">
        {Array.from({ length: 8 }, (_u, i) => {
          const visible = calm ? true : i <= idx;
          const isNew = i === idx;
          return (
            <i
              key={i}
              className={`absolute h-[5px] rounded-[2px] transition-opacity duration-slow ${
                isNew ? "bg-blue" : "bg-ink/35"
              }`}
              style={{
                left: `${i * tileSize}px`,
                top: `${6 + i * 9}px`,
                width: `${kLen * tileSize}px`,
                opacity: visible ? 1 : 0,
              }}
            />
          );
        })}
      </div>
    </>
  );
};

const K3Stage = ({ tileSize, kLen, ui, calm }: { tileSize: number; kLen: number; ui: Record<string, string>; calm: boolean }) => {
  const [settled, setSettled] = useState(calm);
  useEffect(() => {
    if (calm) {
      setSettled(true);
      return undefined;
    }
    setSettled(false);
    const id = window.setTimeout(() => setSettled(true), 1500);
    return () => window.clearTimeout(id);
  }, [calm]);

  const piece = useMemo(() => SEQ.slice(AT, AT + kLen), [kLen]);
  const mirror = useMemo(() => piece.slice().reverse().map((c) => PAIR[c]), [piece]);
  const mineFirst = piece.join("") <= mirror.join("");

  return (
    <>
      <WindowStrand tileSize={tileSize} kLen={kLen} winAt={AT} windowLabel={uiText(ui, "window")} calm={calm} />
      <div className={`transition-opacity duration-slow ${settled && !mineFirst ? "opacity-30" : ""}`}>
        <p className="mb-1 text-tiny text-muted-text">
          <b className="mr-1 font-semibold text-ink">{uiText(ui, "read")}</b>
          {uiText(ui, "piece")}
          {settled && mineFirst ? (
            <span className="ml-2 inline-block rounded-full bg-blue px-2 py-0 text-tiny font-semibold text-paper-pure">
              {uiText(ui, "kept")}
            </span>
          ) : null}
        </p>
        <DnaTiles seq={piece} size={tileSize} withLetters={tileSize >= 12} />
      </div>
      <div
        className={`transition-opacity duration-slow ${calm ? "" : "flow-in"} ${settled && mineFirst ? "opacity-30" : ""}`}
        style={{ animationDelay: calm ? "0ms" : "500ms" }}
      >
        <p className="mb-1 text-tiny text-muted-text">
          <b className="mr-1 font-semibold text-ink">{uiText(ui, "mirror")}</b>
          {uiText(ui, "how")}
          {settled && !mineFirst ? (
            <span className="ml-2 inline-block rounded-full bg-blue px-2 py-0 text-tiny font-semibold text-paper-pure">
              {uiText(ui, "kept")}
            </span>
          ) : null}
        </p>
        <DnaTiles seq={mirror} size={tileSize} withLetters={tileSize >= 12} />
      </div>
    </>
  );
};

const K4Stage = ({ tileSize, kLen, ui, calm }: { tileSize: number; kLen: number; ui: Record<string, string>; calm: boolean }) => {
  const piece = useMemo(() => SEQ.slice(AT, AT + kLen), [kLen]);
  return (
    <>
      <div>
        <StageLabel text={uiText(ui, "bits")} />
        <div className="flex flex-wrap gap-x-5 gap-y-1.5">
          {(["A", "C", "G", "T"] as const).map((base) => {
            const [b0, b1] = baseBits(base);
            return (
              <span key={base} className="inline-flex items-center gap-1.5">
                <DnaTiles seq={[base]} size={18} withLetters />
                <span className="inline-flex" aria-hidden="true">
                  <i className={`h-[10px] w-[11px] rounded-[1px] border ${b0 ? "border-ink bg-ink" : "border-ink/35"}`} />
                  <i className={`h-[10px] w-[11px] rounded-[1px] border ${b1 ? "border-ink bg-ink" : "border-ink/35"}`} />
                </span>
              </span>
            );
          })}
        </div>
      </div>

      <div>
        <StageLabel text={uiText(ui, "piece")} />
        <DnaTiles seq={piece} size={tileSize} withLetters={tileSize >= 12} />
        <div
          aria-hidden="true"
          className={`mt-1 flex ${calm ? "" : "flow-in"}`}
          style={{ animationDelay: calm ? "0ms" : "500ms" }}
        >
          {piece.flatMap((base, idx) => {
            const [b0, b1] = baseBits(base);
            return [
              <i
                key={`${idx}-0`}
                className={`h-[10px] flex-none rounded-[1px] border ${b0 ? "border-ink bg-ink" : "border-ink/35"}`}
                style={{ width: `${tileSize / 2}px` }}
              />,
              <i
                key={`${idx}-1`}
                className={`h-[10px] flex-none rounded-[1px] border ${b1 ? "border-ink bg-ink" : "border-ink/35"}`}
                style={{ width: `${tileSize / 2}px` }}
              />,
            ];
          })}
        </div>
      </div>

      <span
        className={`self-start rounded-full border border-ink/35 px-3 py-0.5 text-control font-semibold text-ink ${calm ? "" : "flow-in"}`}
        style={{ animationDelay: calm ? "0ms" : "1000ms" }}
      >
        &darr;&nbsp;&nbsp;{uiText(ui, "number")}
      </span>
    </>
  );
};

const K5Stage = ({ innerWidth, ui, calm }: { innerWidth: number; ui: Record<string, string>; calm: boolean }) => {
  const [phase, setPhase] = useState<0 | 1 | 2>(calm ? 2 : 0);
  useEffect(() => {
    if (calm) {
      setPhase(2);
      return undefined;
    }
    setPhase(0);
    const t1 = window.setTimeout(() => setPhase(1), 700);
    const t2 = window.setTimeout(() => setPhase(2), 1600);
    return () => {
      window.clearTimeout(t1);
      window.clearTimeout(t2);
    };
  }, [calm]);

  const { vals, dup, keepers } = useMemo(() => {
    const raw = makeLetters(30, 7).map((c, i) => 12 + ((i * 37 + "ACGT".indexOf(c) * 23) % 80));
    raw[11] = raw[4]!;
    raw[23] = raw[17]!;
    const rank = raw.map((_v, i) => i).sort((a, b) => raw[a]! - raw[b]! || a - b);
    const dupSet = new Set(rank.filter((i, r) => r > 0 && raw[rank[r - 1]!] === raw[i]));
    const keepList = rank.filter((i) => !dupSet.has(i));
    return { vals: raw, dup: dupSet, keepers: keepList };
  }, []);

  const gap = Math.max(8, Math.min(18, Math.floor(innerWidth / vals.length)));
  const sketchSet = useMemo(() => new Set(keepers.slice(0, 5)), [keepers]);

  return (
    <>
      <StageLabel text={uiText(ui, "number")} />
      <div className="relative h-[96px]" aria-hidden="true">
        {vals.map((heightPx, i) => {
          const isDup = dup.has(i);
          const isSketch = phase >= 2 && sketchSet.has(i);
          const left = phase >= 1 && !isDup ? keepers.indexOf(i) * gap : i * gap;
          return (
            <i
              key={i}
              className={`absolute bottom-0 w-2 rounded-t-[2px] ${isSketch ? "bg-blue" : "bg-ink/35"}`}
              style={{
                height: `${heightPx}px`,
                left: `${left}px`,
                opacity: phase >= 1 && isDup ? 0 : 1,
                transition: calm ? "none" : "left 700ms cubic-bezier(0.2, 0.8, 0.2, 1), opacity 400ms, background-color 400ms",
              }}
            />
          );
        })}
      </div>
      <div
        className="flex justify-between text-tiny text-muted-text"
        style={{
          width: `${Math.min(innerWidth, vals.length * gap)}px`,
          visibility: phase >= 2 ? "visible" : "hidden",
        }}
      >
        <span>
          <b className="font-semibold text-blue">{uiText(ui, "sketch")} </b>
          {uiText(ui, "small")}
        </span>
        <span>{uiText(ui, "large")}</span>
      </div>
    </>
  );
};

const MatrixGridStage = ({
  merge,
  stays,
  firstIdx,
  calm,
}: {
  merge: boolean;
  stays: boolean[];
  firstIdx: number[];
  calm: boolean;
}) => {
  const [phase, setPhase] = useState<0 | 1 | 2>(calm ? (merge ? 2 : 1) : 0);
  useEffect(() => {
    if (calm) {
      setPhase(merge ? 2 : 1);
      return undefined;
    }
    setPhase(0);
    const t1 = window.setTimeout(() => setPhase(1), 1300);
    const t2 = merge ? window.setTimeout(() => setPhase(2), 2000) : null;
    return () => {
      window.clearTimeout(t1);
      if (t2 !== null) window.clearTimeout(t2);
    };
  }, [merge, calm]);

  const firsts = useMemo(() => [...new Set(firstIdx.filter((x) => x >= 0))], [firstIdx]);
  const cellPitch = 19;

  return (
    <div className="relative" style={{ height: `${GRID[0].length * cellPitch}px` }} aria-hidden="true">
      {GRID.map((col, i) => {
        const keep = stays[i] ?? false;
        const left = phase >= 2 && merge && keep ? firsts.indexOf(firstIdx[i]!) * 24 : i * 24;
        const opacity = phase >= 1 && !keep ? (merge ? 0 : 0.16) : 1;
        return (
          <div
            key={i}
            className={`absolute top-0 flex flex-col gap-[3px] ${calm ? "" : "flow-in"}`}
            style={{
              left: `${left}px`,
              opacity,
              animationDelay: calm ? "0ms" : `${i * 70}ms`,
              transition: calm ? "none" : "left 700ms cubic-bezier(0.2, 0.8, 0.2, 1), opacity 500ms",
            }}
          >
            {[...col].map((bit, r) => (
              <i
                key={r}
                className={`h-4 w-4 rounded-[3px] border-[1.5px] ${
                  bit === "1" ? "border-blue bg-blue" : "border-ink/35 bg-transparent"
                }`}
              />
            ))}
          </div>
        );
      })}
    </div>
  );
};

export interface PieceCountFlowProps {
  flow?: FlowData;
  startStep?: string;
  calm?: boolean;
}

export const PieceCountFlow = ({ flow = FLOW_DATA, startStep, calm: calmProp }: PieceCountFlowProps) => {
  const reducedMotion = useReducedMotion();
  const calm = calmProp ?? reducedMotion;
  const view = flow.views.piece_count;
  const ui = view.ui;
  const nodes = useMemo(() => buildNodeMap(flow), [flow]);
  const steps = view.steps;

  const kLen = flow.counts.k?.v ?? 31;
  const minIn = flow.counts.min_in?.v ?? 2;
  const minOut = flow.counts.min_out?.v ?? 2;

  const stays = useMemo(
    () => GRID.map((c) => countOnes(c) >= minIn && c.length - countOnes(c) >= minOut),
    [minIn, minOut],
  );
  const firstIdx = useMemo(
    () => GRID.map((c, i) => (stays[i] ? GRID.findIndex((d, j) => stays[j] && d === c) : -1)),
    [stays],
  );

  const [stepIdx, setStepIdx] = useState<number>(() => Math.max(0, steps.indexOf(startStep ?? view.start)));
  const [stageWidth, setStageWidth] = useState(720);
  const stageRef = useRef<HTMLDivElement | null>(null);

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

  const go = useCallback(
    (nextIdx: number) => {
      setStepIdx(((nextIdx % steps.length) + steps.length) % steps.length);
    },
    [steps.length],
  );

  const activeId = steps[stepIdx] ?? steps[0]!;
  const node: FlowNode | undefined = nodes[activeId];
  const decisions = node?.decisions ?? [];
  const toolText = (node?.fields ?? [])
    .filter((f) => f[0] === "Tool")
    .map((f) => `${f[0]}: ${fillTemplate(flow, f[1])}`)
    .join("");

  const innerWidth = Math.max(260, stageWidth - 28);
  const tileSize = Math.max(6, Math.min(15, Math.floor(innerWidth / SHOWN)));

  const stepItems = steps.map((id) => ({
    id,
    label: fillTemplate(flow, nodes[id]?.title ?? id),
  }));

  return (
    <section
      data-flow-diagram="piece-count"
      tabIndex={-1}
      onKeyDown={(e) => {
        if (e.key === "ArrowRight") go(stepIdx + 1);
        else if (e.key === "ArrowLeft") go(stepIdx - 1);
      }}
      className="mx-auto w-full max-w-[760px] text-ink outline-none"
    >
      <FlowHeader flow={flow} kicker={view.kicker} title={view.title} sub={view.sub} />

      <StepStrip steps={stepItems} activeIndex={stepIdx} onSelect={go} ariaLabel={fillTemplate(flow, view.kicker)} />

      <WalkControls
        className="mb-3"
        backLabel={uiText(ui, "back")}
        nextLabel={uiText(ui, "next")}
        onBack={() => go(stepIdx - 1)}
        onNext={() => go(stepIdx + 1)}
      />

      {node !== undefined ? (
        <div className="min-h-[340px]">
          <p className="mb-3 max-w-[62ch] text-body text-ink">{fillTemplate(flow, node.line)}</p>

          <div ref={stageRef} className="overflow-hidden rounded-lg border border-hair bg-paper-pure p-3.5">
            <div key={activeId} className="flex min-h-[214px] flex-col gap-2.5">
              {activeId === "k1" ? <ContigsStrip caption={uiText(ui, "drawing")} /> : null}
              {activeId === "k2" ? <K2Stage tileSize={tileSize} kLen={kLen} ui={ui} calm={calm} /> : null}
              {activeId === "k3" ? <K3Stage tileSize={tileSize} kLen={kLen} ui={ui} calm={calm} /> : null}
              {activeId === "k4" ? <K4Stage tileSize={tileSize} kLen={kLen} ui={ui} calm={calm} /> : null}
              {activeId === "k5" ? <K5Stage innerWidth={innerWidth} ui={ui} calm={calm} /> : null}
              {activeId === "k6" ? (
                <>
                  <StageLabel text={`${uiText(ui, "genomes")} \u2193   ${uiText(ui, "pieces")} \u2192`} />
                  <MatrixGridStage merge={false} stays={stays} firstIdx={firstIdx} calm={calm} />
                  <div className="flex flex-wrap gap-x-4 gap-y-1 text-tiny text-muted-text">
                    <span className="inline-flex items-center gap-1.5">
                      <i aria-hidden="true" className="inline-block h-[11px] w-[11px] rounded-[2px] border-[1.5px] border-blue bg-blue" />
                      {uiText(ui, "has")}
                    </span>
                    <span className="inline-flex items-center gap-1.5">
                      <i aria-hidden="true" className="inline-block h-[11px] w-[11px] rounded-[2px] border-[1.5px] border-ink/35" />
                      {uiText(ui, "lacks")}
                    </span>
                    <span className="inline-flex items-center gap-1.5">
                      <i aria-hidden="true" className="inline-block h-[11px] w-[11px] rounded-[2px] border-[1.5px] border-ink/35 opacity-30" />
                      {uiText(ui, "dropped")}
                    </span>
                  </div>
                </>
              ) : null}
              {activeId === "k7" ? (
                <>
                  <StageLabel text={`${uiText(ui, "pattern")} \u2192`} />
                  <MatrixGridStage merge stays={stays} firstIdx={firstIdx} calm={calm} />
                  {node.bars !== undefined ? <FlowBars flow={flow} node={node} calm={calm} /> : null}
                </>
              ) : null}
            </div>
          </div>

          {toolText !== "" ? <p className="mt-2 font-mono text-tiny text-muted-text">{toolText}</p> : null}

          <div className="mt-4">
            <p className="mb-1 text-tiny font-semibold uppercase tracking-[0.06em] text-muted-text">
              {decisions.length > 0 ? uiText(ui, "choices") : uiText(ui, "none")}
            </p>
            {decisions.length > 0 ? (
              <ChoiceRowsList flow={flow} nodes={nodes} ids={decisions} recordLabel={uiText(ui, "record")} />
            ) : null}
          </div>
        </div>
      ) : null}
    </section>
  );
};
