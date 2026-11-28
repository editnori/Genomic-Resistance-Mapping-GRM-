// Shared building blocks for the five flow diagrams, styled with GRM's own
// theme tokens and UI primitives (light and dark follow ThemeProvider).

import { useEffect, useState, type CSSProperties, type ReactNode } from "react";
import { ACTION_BUTTON, FOCUS_RING } from "../ui/primitives";
import {
  fillTemplate,
  formatCount,
  makeLetters,
  uiText,
  type AltState,
  type BaseLetter,
  type CardUse,
  type FlowData,
  type FlowNode,
  type GlyphId,
  type NodeMark,
} from "./model";

// ---- Header & step navigation ----------------------------------------------

export const FlowHeader = ({
  flow,
  kicker,
  title,
  sub,
}: {
  flow: FlowData;
  kicker: string;
  title: string;
  sub: string;
}) => (
  <header className="mb-4">
    <p className="text-tiny font-semibold uppercase tracking-[0.06em] text-muted-text">{fillTemplate(flow, kicker)}</p>
    <h2 className="mt-1 text-headline font-semibold tracking-tight text-ink [text-wrap:balance]">{fillTemplate(flow, title)}</h2>
    <p className="mt-1 text-note text-muted-text">{fillTemplate(flow, sub)}</p>
  </header>
);

export const StepStrip = ({
  steps,
  activeIndex,
  onSelect,
  ariaLabel,
}: {
  steps: { id: string; label: string }[];
  activeIndex: number;
  onSelect(index: number): void;
  ariaLabel: string;
}) => (
  <nav aria-label={ariaLabel} className="mb-3 flex flex-wrap items-center gap-1.5">
    {steps.map((step, i) => {
      const current = i === activeIndex;
      return (
        <span key={step.id} className="inline-flex items-center gap-1.5">
          {i > 0 ? (
            <span aria-hidden="true" className="hidden text-control text-muted-text sm:inline">
              &rarr;
            </span>
          ) : null}
          <button
            type="button"
            data-flow-step={step.id}
            aria-current={current ? "step" : undefined}
            onClick={() => onSelect(i)}
            className={`ui-control-target inline-flex min-h-[28px] items-center whitespace-nowrap rounded-md border px-3 py-1 text-control font-semibold leading-none transition-[background-color,border-color,color,transform] duration-fast ease-out active:scale-[0.98] ${FOCUS_RING} ${
              current
                ? "border-blue bg-blue text-paper-pure"
                : "border-ink/25 bg-paper-pure text-ink hover:border-blue hover:bg-surface-subtle hover:text-blue"
            }`}
          >
            {step.label}
          </button>
        </span>
      );
    })}
  </nav>
);

export const WalkControls = ({
  backLabel,
  nextLabel,
  onBack,
  onNext,
  middle,
  right,
  className = "",
}: {
  backLabel: string;
  nextLabel: string;
  onBack(): void;
  onNext(): void;
  middle?: ReactNode;
  right?: ReactNode;
  className?: string;
}) => (
  <div className={`flex flex-wrap items-center gap-2 ${className}`}>
    <button type="button" onClick={onBack} className={ACTION_BUTTON}>
      {backLabel}
    </button>
    {middle}
    <button type="button" onClick={onNext} className={ACTION_BUTTON}>
      {nextLabel}
    </button>
    {right !== undefined ? <div className="ml-auto">{right}</div> : null}
  </div>
);

// ---- Choice rows & inline detail panel -------------------------------------

const dotClass = (mark: NodeMark | undefined): string => {
  switch (mark) {
    case "kept":
      return "border-ink bg-ink";
    case "cand":
      return "border-blue bg-blue";
    case "run":
      return "border-blue bg-transparent";
    case "none":
      return "border-coral-text bg-transparent";
    case "stated":
    case "wait":
      return "border-dotted border-ink/40 bg-transparent";
    default:
      return "border-ink/40 bg-transparent";
  }
};

const tagClass = (mark: NodeMark | undefined): string => {
  switch (mark) {
    case "cand":
      return "border-blue bg-blue text-paper-pure";
    case "run":
      return "border-blue text-blue";
    case "none":
      return "border-coral-text text-coral-text";
    default:
      return "border-ink/30 text-ink";
  }
};

const altClass = (state: AltState): string => {
  if (state === "pass") return "border-blue text-blue font-semibold";
  if (state === "open") return "border-dashed border-ink/20 text-muted-text";
  return "border-ink/15 text-muted-text";
};

export const ChoiceDetail = ({
  flow,
  nodes,
  id,
  recordLabel,
  className = "",
}: {
  flow: FlowData;
  nodes: Record<string, FlowNode>;
  id: string;
  recordLabel: string;
  className?: string;
}) => {
  const n = nodes[id];
  if (n === undefined) return null;
  const linked = n.test !== undefined ? nodes[n.test] : undefined;
  const entries: [string, string][] = [...(n.fields ?? [])];
  if (linked !== undefined) {
    for (const [k, v] of linked.fields ?? []) {
      if (k === "Tested" || k === "Counted" || k === "Led to") entries.push([k, v]);
    }
  }
  const records = [...new Set([...n.src, ...(linked?.src ?? [])].map((x) => `${x.file.split("/").pop() ?? x.file}:${x.line}`))];

  return (
    <div className={`bg-blue/[0.08] px-3 py-2.5 sm:pl-8 ${className}`}>
      <dl className="m-0 grid grid-cols-1 gap-x-3 gap-y-1.5 text-control sm:grid-cols-[7.25rem_minmax(0,1fr)]">
        {entries.map(([k, v], idx) => (
          <div key={`${k}-${idx}`} className="contents">
            <dt className="text-muted-text">{k}</dt>
            <dd className="m-0 max-w-[64ch] text-ink">{fillTemplate(flow, v)}</dd>
          </div>
        ))}
        <dt className="text-muted-text">{recordLabel}</dt>
        <dd className="m-0 font-mono text-tiny text-muted-text [overflow-wrap:anywhere]">{records.join(", ")}</dd>
      </dl>
    </div>
  );
};

export const ChoiceRowItem = ({
  flow,
  node,
  expanded,
  onSelect,
  branchTick = false,
}: {
  flow: FlowData;
  node: FlowNode;
  expanded: boolean;
  onSelect(): void;
  branchTick?: boolean;
}) => (
  <button
    type="button"
    data-choice-id={node.id}
    aria-expanded={expanded}
    onClick={onSelect}
    className={`relative grid w-full grid-cols-[12px_minmax(0,1fr)] items-start gap-x-2.5 gap-y-1 border-t border-hair px-2.5 py-2 text-left text-control leading-[1.35] transition-colors duration-fast ease-out sm:grid-cols-[12px_minmax(0,1fr)_minmax(0,1.55fr)_5.75rem] ${FOCUS_RING} ${
      expanded ? "bg-blue/[0.09]" : "bg-transparent hover:bg-surface-subtle"
    }`}
    style={node.bare ? { boxShadow: "inset 3px 0 0 rgb(var(--coral-text))" } : undefined}
  >
    {branchTick ? (
      <span aria-hidden="true" className="pointer-events-none absolute -left-[15px] top-[16px] w-[14px] border-t border-ink/35" />
    ) : null}
    <span aria-hidden="true" className={`mt-1 h-2.5 w-2.5 rounded-full border-2 ${dotClass(node.mark)}`} />
    <span className="font-semibold text-ink">{fillTemplate(flow, node.title)}</span>
    {node.alts !== undefined && node.alts.length > 0 ? (
      <span className="col-start-2 flex flex-wrap gap-1 sm:col-start-auto">
        {node.alts.map(([label, state], idx) => (
          <span
            key={`${label}-${idx}`}
            className={`inline-flex items-center rounded-full border px-2 py-0.5 text-tiny leading-tight ${altClass(state)}`}
          >
            {state === "no" ? <span aria-hidden="true">&times;&nbsp;</span> : null}
            {state === "pass" ? <span aria-hidden="true">&#10003;&nbsp;</span> : null}
            {fillTemplate(flow, label)}
          </span>
        ))}
      </span>
    ) : (
      <span className="col-start-2 text-muted-text sm:col-start-auto">{fillTemplate(flow, node.line)}</span>
    )}
    <span className="col-start-2 text-left sm:col-start-auto sm:text-right">
      {node.tag !== undefined && node.tag !== "" ? (
        <span className={`inline-block whitespace-nowrap rounded-full border px-2 py-0.5 text-tiny font-semibold leading-tight ${tagClass(node.mark)}`}>
          {node.tag}
        </span>
      ) : null}
    </span>
  </button>
);

export const ChoiceRowsList = ({
  flow,
  nodes,
  ids,
  recordLabel,
}: {
  flow: FlowData;
  nodes: Record<string, FlowNode>;
  ids: string[];
  recordLabel: string;
}) => {
  const [openId, setOpenId] = useState<string>(ids[0] ?? "");
  useEffect(() => {
    setOpenId(ids[0] ?? "");
  }, [ids.join(",")]);

  return (
    <div>
      {ids.map((id) => {
        const node = nodes[id];
        if (node === undefined) return null;
        const expanded = openId === id;
        return (
          <div key={id}>
            <ChoiceRowItem flow={flow} node={node} expanded={expanded} onSelect={() => setOpenId(id)} />
            {expanded ? <ChoiceDetail flow={flow} nodes={nodes} id={id} recordLabel={recordLabel} /> : null}
          </div>
        );
      })}
    </div>
  );
};

// ---- Shared stage visual primitives (bars, stats, record card, DNA tiles) --

export const FlowBars = ({
  flow,
  node,
  calm = false,
}: {
  flow: FlowData;
  node: FlowNode;
  calm?: boolean;
}) => {
  const bars = node.bars ?? [];
  const [grown, setGrown] = useState(calm);
  useEffect(() => {
    if (calm) {
      setGrown(true);
      return undefined;
    }
    setGrown(false);
    let f2 = 0;
    const f1 = window.requestAnimationFrame(() => {
      f2 = window.requestAnimationFrame(() => setGrown(true));
    });
    return () => {
      window.cancelAnimationFrame(f1);
      window.cancelAnimationFrame(f2);
    };
  }, [calm, node.id]);

  const countVal = node.count !== undefined ? (flow.counts[node.count]?.v ?? 0) : 0;
  const whole = Math.max(countVal, ...bars.map((b) => flow.counts[b[1]]?.v ?? 0), 1);

  return (
    <div className="grid grid-cols-[minmax(0,1fr)_auto] items-center gap-x-3 gap-y-2 text-control sm:grid-cols-[minmax(0,1.1fr)_minmax(0,1.3fr)_auto]">
      {bars.map(([name, key, cls], idx) => {
        const val = flow.counts[key]?.v ?? 0;
        const widthPct = `${Math.max((val / whole) * 100, 0.6)}%`;
        const fillClass =
          cls === "on" ? "bg-blue" : cls === "hatch" ? "flow-bar-hatch" : "bg-[var(--series-neutral)]";
        return (
          <div key={`${key}-${idx}`} className="contents">
            <span className="text-ink">{fillTemplate(flow, name)}</span>
            <span className="order-3 col-span-2 h-3.5 overflow-hidden rounded-[3px] bg-surface-subtle sm:order-none sm:col-span-1">
              <span
                className={`block h-full rounded-[3px] transition-[width] duration-slow ease-out ${fillClass}`}
                style={{ width: grown ? widthPct : "0%" }}
              />
            </span>
            <span className="text-right font-semibold tabular-nums text-ink">{formatCount(flow, key)}</span>
          </div>
        );
      })}
    </div>
  );
};

export const FlowStats = ({ flow, node }: { flow: FlowData; node: FlowNode }) => {
  const stats = node.stats ?? [];
  return (
    <div className="grid grid-cols-1 gap-3 sm:grid-cols-3 sm:gap-0">
      {stats.map(([name, key, note], idx) => (
        <div
          key={`${key}-${idx}`}
          className={`flex flex-col gap-0.5 ${idx > 0 ? "sm:border-l sm:border-dashed sm:border-hair sm:pl-3.5" : ""} sm:pr-3.5`}
        >
          <span className="text-tiny text-muted-text">{fillTemplate(flow, name)}</span>
          <span className="text-stat font-semibold tabular-nums text-ink">{formatCount(flow, key)}</span>
          <span className="text-tiny text-muted-text">{fillTemplate(flow, note)}</span>
        </div>
      ))}
    </div>
  );
};

export const FlowCard = ({
  flow,
  node,
  ui,
}: {
  flow: FlowData;
  node: FlowNode;
  ui: Record<string, string>;
}) => {
  const rows = node.card ?? [];
  const groups: [CardUse, [string, CardUse, string][]][] = [];
  for (const row of rows) {
    const last = groups[groups.length - 1];
    if (last === undefined || last[0] !== row[1]) groups.push([row[1], [row]]);
    else last[1].push(row);
  }

  return (
    <div>
      {node.cardtitle !== undefined ? <p className="mb-2 text-control font-semibold text-ink">{fillTemplate(flow, node.cardtitle)}</p> : null}
      {groups.map(([use, items], gi) => {
        const detailed = items.filter((r) => r[2] !== "");
        const bare = items.filter((r) => r[2] === "");
        const mutedGroup = use === "unread";
        return (
          <div
            key={`${use}-${gi}`}
            className={`grid grid-cols-1 gap-x-3 border-t border-hair py-1.5 sm:grid-cols-[9.5rem_minmax(0,1fr)] ${mutedGroup ? "text-muted-text" : "text-ink"}`}
          >
            <p className={`mb-0.5 text-tiny font-semibold sm:mb-0 ${use === "used" ? "text-blue" : "text-muted-text"}`}>
              {uiText(ui, `k_${use}`)}
            </p>
            <div className="space-y-0.5 text-control">
              {detailed.map((r, ri) => (
                <div key={`${r[0]}-${ri}`} className="grid grid-cols-[7rem_minmax(0,1fr)] gap-x-2.5 sm:grid-cols-[8.25rem_minmax(0,1fr)]">
                  <span className={use === "used" ? "font-semibold" : "font-normal"}>{fillTemplate(flow, r[0])}</span>
                  <span className="text-muted-text">{fillTemplate(flow, r[2])}</span>
                </div>
              ))}
              {bare.length > 0 ? <div className="font-normal">{bare.map((r) => fillTemplate(flow, r[0])).join(", ")}</div> : null}
            </div>
          </div>
        );
      })}
    </div>
  );
};

export const baseTileStyle = (base: BaseLetter): CSSProperties => {
  switch (base) {
    case "A":
      return { background: "var(--base-a)", color: "#ffffff" };
    case "C":
      return { background: "var(--base-c)", color: "#0b0b0b" };
    case "G":
      return { background: "var(--base-g)", color: "#ffffff" };
    case "T":
      return { background: "var(--base-t)", color: "var(--base-t-ink)" };
  }
};

export const DnaTiles = ({
  seq,
  size,
  withLetters,
}: {
  seq: readonly BaseLetter[];
  size: number;
  withLetters: boolean;
}) => (
  <div className="flex" aria-hidden="true">
    {seq.map((base, idx) => (
      <span
        key={idx}
        className="grid flex-none place-items-center rounded-[2px] text-[10px] font-semibold leading-none"
        style={{
          width: `${size}px`,
          height: `${size}px`,
          boxShadow: "inset 0 0 0 0.5px rgb(var(--paper-pure))",
          ...baseTileStyle(base),
        }}
      >
        {withLetters ? base : null}
      </span>
    ))}
  </div>
);

const CONTIG_SEQ = makeLetters(160, 20261127);
const CONTIG_CUTS = [52, 34, 40, 20, 14] as const;

export const ContigsStrip = ({ caption }: { caption: string }) => {
  let offset = 0;
  const stretches = CONTIG_CUTS.map((len) => {
    const slice = CONTIG_SEQ.slice(offset, offset + len);
    offset += len;
    return slice;
  });
  return (
    <div className="flex flex-col gap-[7px]">
      {stretches.map((slice, idx) => (
        <DnaTiles key={idx} seq={slice} size={7} withLetters={false} />
      ))}
      <p className="mt-1 text-tiny text-muted-text">{caption}</p>
    </div>
  );
};

// ---- Step glyphs for GenomePathFlow ----------------------------------------

export const FlowGlyph = ({ glyph }: { glyph: GlyphId }) => {
  const stroke = { fill: "none", stroke: "currentColor", strokeWidth: 1.6, strokeLinecap: "round" as const, strokeLinejoin: "round" as const };
  return (
    <svg viewBox="0 0 20 20" aria-hidden="true" className="block h-[18px] w-[18px]">
      {glyph === "records" ? <path d="M3 5h14M3 10h14M3 15h9" {...stroke} /> : null}
      {glyph === "datasets" ? (
        <>
          <circle cx={5.5} cy={10} r={3.6} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={14.5} cy={10} r={3.6} fill="none" stroke="currentColor" strokeWidth={1.5} />
        </>
      ) : null}
      {glyph === "genome" ? (
        <>
          <path d="M2 5c5 0 5 10 8 10s3-10 8-10" {...stroke} />
          <path d="M2 15c5 0 5-10 8-10s3 10 8 10" {...stroke} />
          <path d="M5.4 8.4v3.2M14.6 8.4v3.2" {...stroke} strokeWidth={1.3} />
        </>
      ) : null}
      {glyph === "pieces" ? (
        <>
          <rect x={1} y={6} width={3.4} height={8} rx={1} fill="var(--base-a)" />
          <rect x={5.4} y={6} width={3.4} height={8} rx={1} fill="var(--base-g)" />
          <rect x={9.8} y={6} width={3.4} height={8} rx={1} fill="var(--base-t)" />
          <rect x={14.2} y={6} width={3.4} height={8} rx={1} fill="var(--base-c)" />
          <path d="M0.6 3.4h12.8v13.2H0.6z" {...stroke} strokeWidth={1.3} />
        </>
      ) : null}
      {glyph === "split" ? (
        <>
          <circle cx={4} cy={6} r={2} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={4} cy={14} r={2} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={10} cy={10} r={2} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <path d="M14.4 3v14" {...stroke} strokeWidth={1.3} strokeDasharray="2 2" />
          <circle cx={18} cy={10} r={2} fill="none" stroke="currentColor" strokeWidth={1.5} />
        </>
      ) : null}
      {glyph === "filter" ? <path d="M2.5 4h15l-5.6 6.6V16l-3.8-1.6v-3.8z" {...stroke} /> : null}
      {glyph === "models" ? (
        <>
          <circle cx={4} cy={5} r={1.9} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={16} cy={5} r={1.9} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={10} cy={10} r={1.9} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={4} cy={15} r={1.9} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={16} cy={15} r={1.9} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
        </>
      ) : null}
      {glyph === "served" ? (
        <>
          <circle cx={4} cy={5} r={1.9} fill="none" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={16} cy={5} r={1.9} fill="none" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={10} cy={10} r={2.6} fill="currentColor" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={4} cy={15} r={1.9} fill="none" stroke="currentColor" strokeWidth={1.5} />
          <circle cx={16} cy={15} r={1.9} fill="none" stroke="currentColor" strokeWidth={1.5} />
        </>
      ) : null}
      {glyph === "score" ? (
        <>
          <path d="M3 14a7 7 0 0 1 14 0" {...stroke} />
          <path d="M10 14l3.6-5" {...stroke} />
        </>
      ) : null}
    </svg>
  );
};
