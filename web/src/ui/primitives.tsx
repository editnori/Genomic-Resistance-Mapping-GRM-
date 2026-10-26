// Shared UI atoms. Token classes only; every control carries the five
// interaction states (rest, hover, 0.98 press, 2px blue focus ring,
// disabled with no hover promise).

import { useId, useLayoutEffect, useRef, useState, type ReactNode } from "react";
import { ChevronDown } from "lucide-react";
import { barWords, verdictOf, type ReferenceBar } from "../data/metrics";
import { pct } from "../data/format";
import { GrmMark } from "../app/Mark";
import { HoverTip } from "./HoverTip";

export const FOCUS_RING =
  "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-blue";

/** The shared content rail: 1560px cap, 16px gutters below sm, 24px above.
 *  The top bar uses the same gutters, so the brand and every view's title
 *  start on one left edge. */
export const VIEW_SHELL = "mx-auto w-full max-w-[1560px] px-4 sm:px-6";

/** An outlined pill action: Present, the inspector's actions, slide steps. */
export const PILL_ACTION = `ui-control-target inline-flex h-8 items-center justify-center gap-1.5 whitespace-nowrap rounded-full border border-hair bg-paper-pure px-3.5 text-control font-medium leading-none text-ink transition-[color,border-color,background-color,transform] duration-fast ease-out active:scale-[0.98] enabled:hover:border-blue enabled:hover:text-blue disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING}`;

// ---- Lozenge: a state or data label, sentence case, four tones -------------

export type LozengeTone = "ink" | "blue" | "coral" | "live";

const LOZENGE_OUTLINE: Record<LozengeTone, string> = {
  ink: "border-hair text-muted-text",
  blue: "border-blue/40 text-blue",
  coral: "border-coral/60 text-coral-text",
  live: "border-live/60 text-live-text",
};

const LOZENGE_SOLID: Record<LozengeTone, string> = {
  ink: "border-transparent bg-ink text-paper-pure",
  blue: "border-transparent bg-blue text-paper-pure",
  coral: "border-transparent bg-coral-text text-paper-pure",
  live: "border-transparent bg-live-text text-paper-pure",
};

/** Outline names a scope or a data label; solid is a status or verdict the
 *  reader should act on. */
export const Lozenge = ({
  children,
  tone = "ink",
  fill = "outline",
}: {
  children: ReactNode;
  tone?: LozengeTone;
  fill?: "outline" | "solid";
}) => (
  <span
    className={`inline-flex shrink-0 items-center whitespace-nowrap rounded-md border px-2 py-0.5 text-tiny font-medium ${
      fill === "solid" ? LOZENGE_SOLID[tone] : LOZENGE_OUTLINE[tone]
    }`}
  >
    {children}
  </span>
);

// ---- Verdict against a reference bar ---------------------------------------

const BAR_SOURCE: Record<string, string> = {
  VME: "FDA-style bars for susceptibility test devices: VME at most 1.5%, ME at most 3%. A reference, not a regulatory claim.",
  ME: "FDA-style bars for susceptibility test devices: VME at most 1.5%, ME at most 3%. A reference, not a regulatory claim.",
  "MIC agreement": "Essential agreement of at least 90% is the bar used for susceptibility test devices; a reference, not a performance claim.",
};

/** A rate judged against a stated bar. The verdict reads by shape before
 *  colour: meeting the bar is a filled pill, missing it an outlined one. */
export const BarPill = ({ value, bar }: { value: number; bar: ReferenceBar }) => {
  const verdict = verdictOf(value, bar);
  const text = pct(value);
  return (
    <HoverTip
      focusable={false}
      lines={[
        { term: `${bar.name} ${text}`, text: `${verdict === "meets" ? "Meets" : "Misses"} the bar: ${barWords(bar)}.` },
        { source: BAR_SOURCE[bar.name] ?? "A reference bar, not a regulatory claim." },
      ]}
    >
      <span
        data-verdict={verdict}
        className={`inline-flex whitespace-nowrap rounded-md px-1.5 py-px text-control font-medium tabular-nums ${
          verdict === "meets" ? "bg-live-text text-paper-pure" : "text-coral-text ring-1 ring-inset ring-coral-text"
        }`}
      >
        {text}
        <span className="sr-only">
          , {verdict === "meets" ? "meets" : "misses"} {barWords(bar)}
        </span>
      </span>
    </HoverTip>
  );
};

/** The same verdict in words, for prose: "VME 4.9%, misses at most 1.5%". */
export const verdictWords = (value: number, bar: ReferenceBar): string =>
  `${bar.name} ${pct(value)}, ${verdictOf(value, bar) === "meets" ? "meets" : "misses"} ${barWords(bar)}`;

// ---- Count pill: label plus a predicting count -----------------------------

export const CountPill = ({
  label,
  count,
  selected,
  onToggle,
  disabled = false,
}: {
  label: string;
  count?: number;
  selected: boolean;
  onToggle(): void;
  disabled?: boolean;
}) => (
  <button
    type="button"
    aria-pressed={selected}
    disabled={disabled}
    onClick={onToggle}
    className={`ui-control-target inline-flex h-7 items-center gap-x-2 rounded-full border pl-3 text-control font-medium tracking-tight transition-[color,background-color,border-color,transform] duration-fast ease-out active:scale-[0.98] disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING} ${
      count !== undefined ? "pr-1" : "pr-3"
    } ${
      selected
        ? "border-blue bg-blue text-paper-pure hover:bg-blue/90"
        : "border-hair bg-paper-pure text-ink enabled:hover:border-blue/50 enabled:hover:bg-blue/5 enabled:hover:text-blue"
    }`}
  >
    <span className="whitespace-nowrap">{label}</span>
    {count !== undefined ? (
      <span
        className={`inline-flex h-5 min-w-5 items-center justify-center rounded-full px-1 text-tiny tabular-nums ${
          selected ? "bg-paper-pure text-ink" : "bg-surface-strong text-muted-text"
        }`}
      >
        {count}
      </span>
    ) : null}
  </button>
);

// ---- Segmented: a mode switch, ink-filled current segment ------------------

export interface SegmentOption<T extends string> {
  id: T;
  label: string;
  disabled?: boolean;
}

export const Segmented = <T extends string>({
  ariaLabel,
  options,
  value,
  onSelect,
  fill = false,
}: {
  ariaLabel: string;
  options: readonly SegmentOption<T>[];
  value: T;
  onSelect(id: T): void;
  /** Inside the rail a segmented row fills its field evenly. */
  fill?: boolean;
}) => (
  <div
    role="group"
    aria-label={ariaLabel}
    className={`inline-flex min-h-7 max-w-full flex-wrap items-center gap-0.5 rounded-[14px] border border-hair bg-paper-pure p-0.5 ${fill ? "w-full" : "w-max"}`}
  >
    {options.map((o) => {
      const selected = o.id === value;
      return (
        <button
          key={o.id}
          type="button"
          aria-pressed={selected}
          disabled={o.disabled === true}
          onClick={() => onSelect(o.id)}
          className={`ui-control-target inline-flex h-[22px] items-center justify-center whitespace-nowrap rounded-full px-3 text-control font-medium tracking-tight transition-[color,background-color,transform] duration-fast ease-out active:scale-[0.98] disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING} ${
            fill ? "flex-1 basis-0" : ""
          } ${selected ? "bg-ink text-paper-pure" : "text-ink enabled:hover:bg-surface-subtle"}`}
        >
          {o.label}
        </button>
      );
    })}
  </div>
);

// ---- Ghost button: quiet secondary action ----------------------------------

export const GhostButton = ({
  children,
  onClick,
  tone = "ink",
  disabled = false,
  ariaLabel,
}: {
  children: ReactNode;
  onClick(): void;
  tone?: "ink" | "blue";
  disabled?: boolean;
  ariaLabel?: string;
}) => (
  <button
    type="button"
    onClick={onClick}
    disabled={disabled}
    aria-label={ariaLabel}
    className={`ui-control-target inline-flex h-7 items-center gap-1.5 whitespace-nowrap rounded-full border border-hair px-3 text-control font-medium tracking-tight transition-[color,background-color,border-color,transform] duration-fast ease-out active:scale-[0.98] disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING} ${
      tone === "blue"
        ? "text-blue enabled:hover:border-blue/60 enabled:hover:bg-blue/5"
        : "text-ink enabled:hover:border-blue/50 enabled:hover:bg-surface-subtle enabled:hover:text-blue"
    }`}
  >
    {children}
  </button>
);

// ---- Panel state: loading, empty, error ------------------------------------

/** A state in place of content, on one quiet panel. Loading runs the mark's
 *  assembly loop; an error is coral and an alert. */
export const PanelState = ({ kind, title, message, action }: { kind: "loading" | "empty" | "error"; title: string; message?: string; action?: ReactNode }) => (
  <div
    role={kind === "error" ? "alert" : "status"}
    aria-busy={kind === "loading" ? true : undefined}
    className={`rounded-xl border px-4 py-4 ${kind === "error" ? "border-coral/50 bg-coral/5" : "border-hair bg-surface-subtle"}`}
  >
    <p className={`flex items-center gap-2.5 text-body font-semibold ${kind === "error" ? "text-coral-text" : "text-ink"}`}>
      {kind === "loading" ? <GrmMark size={20} motion="loop" /> : null}
      {title}
    </p>
    {message !== undefined ? <p className="type-note mt-1">{message}</p> : null}
    {action !== undefined ? <div className="mt-3">{action}</div> : null}
  </div>
);

// ---- Section: one panel answering one question ------------------------------

/** A panel: the title names the artifact, the answer sentence answers it,
 *  then the evidence. A quiet question shows only where the title alone
 *  leaves the reader guessing. Hairline separated, never a card in a card. */
export const Section = ({
  id,
  title,
  question,
  answer,
  right,
  children,
}: {
  id?: string;
  title: string;
  question?: string;
  answer?: ReactNode;
  right?: ReactNode;
  children?: ReactNode;
}) => (
  <section
    id={id}
    aria-labelledby={id === undefined ? undefined : `${id}-title`}
    className="scroll-mt-4 border-t border-hair pt-6 first:border-t-0 first:pt-0"
  >
    <header className="flex flex-wrap items-start justify-between gap-x-6 gap-y-2">
      <div className="min-w-0 flex-1 basis-60">
        <h2 id={id === undefined ? undefined : `${id}-title`} className="text-body font-semibold tracking-tight text-ink">
          {title}
        </h2>
        {question !== undefined ? <p className="mt-0.5 text-note text-muted-text">{question}</p> : null}
      </div>
      {right !== undefined ? <div className="flex min-w-0 flex-wrap items-center gap-2">{right}</div> : null}
    </header>
    {answer !== undefined ? <p className="type-headline mt-3">{answer}</p> : null}
    {children !== undefined ? <div className="mt-5">{children}</div> : null}
  </section>
);

// ---- Disclosure: "Method and coverage" -------------------------------------

export const Disclosure = ({
  label = "Method and coverage",
  subject,
  children,
  defaultOpen = false,
  ruled = true,
}: {
  label?: string;
  /** A second disclosure in one block names its subject after the label. */
  subject?: string;
  children: ReactNode;
  defaultOpen?: boolean;
  /** False where the next block's hairline already closes this one. */
  ruled?: boolean;
}) => {
  const [open, setOpen] = useState(defaultOpen);
  const bodyId = useId();
  return (
    <div className={ruled ? "border-t border-hair pt-3" : ""}>
      <button
        type="button"
        aria-expanded={open}
        aria-controls={bodyId}
        onClick={() => setOpen((v) => !v)}
        className={`ui-control-target -mx-1 inline-flex items-center gap-1.5 rounded-md px-1 py-0.5 text-control font-medium text-ink transition-colors duration-fast ease-out hover:bg-surface-subtle hover:text-blue active:scale-[0.98] ${FOCUS_RING}`}
      >
        {subject === undefined ? label : `${label} · ${subject}`}
        <ChevronDown size={14} strokeWidth={1.75} aria-hidden="true" className={`transition-transform duration-fast ease-out ${open ? "rotate-180" : ""}`} />
      </button>
      {open ? (
        <div id={bodyId} className="motion-disclosure-in mt-3 space-y-3 text-note text-ink">
          {children}
        </div>
      ) : null}
    </div>
  );
};

// ---- Figure block: title, answer, figure, caption, method and coverage ------

/** What a figure is given to fit: the box's measured size on a slide; on a
 *  page the figure takes the width it measures and its own height. */
export interface Fit {
  width?: number;
  height?: number;
}

export type FigureContent = ReactNode | ((fit: Fit) => ReactNode);

const renderFigure = (figure: FigureContent, fit: Fit): ReactNode => (typeof figure === "function" ? figure(fit) : figure);

/** Measures its box and hands the size to the figure, so a slide figure
 *  fills what the slide leaves and nothing scrolls. */
const FitBox = ({ figure }: { figure: (fit: Fit) => ReactNode }) => {
  const ref = useRef<HTMLDivElement | null>(null);
  const [fit, setFit] = useState<Fit | null>(null);
  useLayoutEffect(() => {
    const el = ref.current;
    if (el === null) return undefined;
    const read = () => setFit({ width: el.clientWidth, height: el.clientHeight });
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  return (
    <div ref={ref} data-figure-fit className="relative min-h-0 min-w-0 flex-1">
      {fit !== null ? figure(fit) : null}
    </div>
  );
};

/** A slide's staged parts: Present wraps each beat in its own reveal. */
export type SlidePart = "result" | "figure" | "caption";

/** One analytical block: the title names the measured quantity, the answer
 *  gives the count and its denominator, the figure shows that number, the
 *  caption says what a mark is and how it was made, and Method and coverage
 *  holds the rest. A figure that is its own answer (the corpus counts) leaves
 *  the sentence out. On a slide a measured figure fills the height the text
 *  leaves; any other figure keeps its own height with the caption right
 *  under it, the two centered in the space left. */
export const FigureBlock = ({
  id,
  title,
  question,
  result,
  figure,
  caption,
  details,
  detailsOpen = false,
  actions,
  variant = "page",
  stage = (_part, node) => node,
}: {
  id: string;
  title: string;
  /** Only where the title alone leaves the reader guessing. */
  question?: string;
  result?: ReactNode;
  figure: FigureContent;
  caption: ReactNode;
  details?: ReactNode;
  detailsOpen?: boolean;
  /** A control that changes this figure, beside its title. */
  actions?: ReactNode;
  variant?: "page" | "slide";
  /** Present's staged reveal for each beat of a slide; `grow` is true for
   *  a measured figure, whose beat must fill the height the text leaves. */
  stage?(part: SlidePart, node: ReactNode, grow: boolean): ReactNode;
}) =>
  variant === "slide" ? (
    <section aria-labelledby={`${id}-title`} className="flex h-full min-h-0 flex-col">
      <h2 id={`${id}-title`} className="text-body font-semibold text-muted-text">
        {title}
      </h2>
      {result === undefined ? null : stage("result", <p className="mt-2 max-w-[56rem] text-stat font-semibold text-ink">{result}</p>, false)}
      {typeof figure === "function" ? (
        <>
          {stage(
            "figure",
            <div data-figure className="mt-4 flex h-full min-h-0 flex-1 flex-col">
              <FitBox figure={figure} />
            </div>,
            true,
          )}
          {stage("caption", <p data-caption className="prose-measure mt-3 text-note text-muted-text">{caption}</p>, false)}
        </>
      ) : (
        <div className="flex min-h-0 flex-1 flex-col justify-center py-4">
          {stage("figure", <div data-figure className="min-w-0">{figure}</div>, false)}
          {stage("caption", <p data-caption className="prose-measure mt-6 text-note text-muted-text">{caption}</p>, false)}
        </div>
      )}
    </section>
  ) : (
    <section id={id} aria-labelledby={`${id}-title`} className="scroll-mt-6 border-t border-hair pt-7">
      <div className="flex flex-wrap items-start justify-between gap-x-6 gap-y-2">
        <div className="min-w-0">
          <h2 id={`${id}-title`} className="text-body font-semibold tracking-tight text-ink">
            {title}
          </h2>
          {question !== undefined ? <p className="mt-0.5 text-note text-muted-text">{question}</p> : null}
        </div>
        {actions !== undefined ? <div className="flex flex-wrap items-center gap-2">{actions}</div> : null}
      </div>
      {result === undefined ? null : <p className="type-headline mt-3">{result}</p>}
      <div className="mt-5 min-w-0">{renderFigure(figure, {})}</div>
      <p className="prose-measure mt-3 text-note text-muted-text">{caption}</p>
      {details !== undefined ? (
        <div className="mt-3">
          <Disclosure defaultOpen={detailsOpen} ruled={false}>
            {details}
          </Disclosure>
        </div>
      ) : null}
    </section>
  );

// ---- Fact rows: label column, value column, hairline between ----------------

/** Labels in one column and values left-aligned in the next, so a long value
 *  wraps under itself; on a phone the label sits above its value. */
export const FactRows = ({ rows }: { rows: { label: string; value: ReactNode }[] }) => (
  <dl className="m-0 grid max-w-3xl grid-cols-1 content-start self-start sm:grid-cols-[minmax(0,13rem)_minmax(0,1fr)] sm:gap-x-6">
    {rows.map((r) => (
      <div key={r.label} className="grid grid-cols-subgrid border-b border-hair py-1.5 last:border-b-0 sm:col-span-2">
        <dt className="text-control font-medium text-muted-text">{r.label}</dt>
        <dd className="m-0 min-w-0 text-control tabular-nums text-ink">{r.value}</dd>
      </div>
    ))}
  </dl>
);

/** A keyboard-reachable horizontal scroll owner for wide tables and plots. */
export const ScrollFrame = ({ children, label }: { children: ReactNode; label: string }) => (
  <div role="region" aria-label={label} tabIndex={0} className={`scroll-thin relative overflow-x-auto overscroll-x-contain ${FOCUS_RING}`}>
    {children}
  </div>
);

/** Legend entry: the swatch mirrors the mark, the name and value are text. */
export const LegendItem = ({ swatch, name, value }: { swatch: ReactNode; name: string; value?: string }) => (
  <span className="inline-flex items-center gap-x-1.5 text-control text-ink">
    <span aria-hidden="true" className="inline-flex shrink-0 items-center">
      {swatch}
    </span>
    <span className="text-muted-text">{name}</span>
    {value !== undefined ? <span className="font-medium tabular-nums">{value}</span> : null}
  </span>
);

export const Dot = ({ color, ring = false, size = 8 }: { color: string; ring?: boolean; size?: number }) => (
  <span
    className="inline-block rounded-full"
    style={ring ? { width: size, height: size, border: `2px solid ${color}` } : { width: size, height: size, background: color }}
  />
);
