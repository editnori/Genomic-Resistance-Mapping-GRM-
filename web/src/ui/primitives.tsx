// Shared UI atoms. Token classes only; every control carries the five
// interaction states (rest, hover, 0.98 press, 2px blue focus ring,
// disabled with no hover promise).

import { useId, useLayoutEffect, useRef, useState, type ReactNode } from "react";
import { ChevronDown } from "lucide-react";
import { barWords, verdictOf, type ReferenceBar } from "../data/metrics";
import { pct } from "../data/format";

export const FOCUS_RING =
  "focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-blue";

/** Centered shells: 80rem for the wide tables (Datasets, Methods), 64rem
 *  for the figure column. Gutters 16px, 24px from 640px, 32px from 1280px. */
export const WIDE_SHELL = "mx-auto w-full max-w-[80rem] px-4 sm:px-6 xl:px-8";
export const FIGURE_SHELL = "mx-auto w-full max-w-[64rem] px-4 sm:px-6 xl:px-8";

// ---- Lozenge: the one solid audit status ------------------------------------

export const Lozenge = ({ children, tone = "neutral" }: { children: ReactNode; tone?: "neutral" | "provisional" }) => (
  <span
    className={`inline-flex shrink-0 items-center whitespace-nowrap rounded-md px-2 py-0.5 text-tiny font-medium ${
      tone === "provisional" ? "bg-miss text-miss-fg" : "bg-surface-control text-ink"
    }`}
  >
    {children}
  </span>
);

// ---- Status cell: a value judged against a stated bar -----------------------

/** A table cell filled by its verdict: solid green with ✓ when the value
 *  meets the bar, the quiet magenta tint with × when it misses. The glyph and
 *  the value are visible, the verdict and bar are in the accessible name, so
 *  color is never the only cue. The header states the bar. */
export const StatusTd = ({ value, bar, className = "" }: { value: number; bar: ReferenceBar; className?: string }) => {
  const verdict = verdictOf(value, bar);
  const meets = verdict === "meets";
  return (
    <td data-verdict={verdict} className={`whitespace-nowrap px-2 text-right tabular-nums ${meets ? "bg-meet text-meet-fg" : "bg-miss text-miss-fg"} ${className}`}>
      <span aria-hidden="true">{meets ? "✓" : "×"} </span>
      {pct(value)}
      <span className="sr-only">
        , {meets ? "meets" : "misses"} {bar.name} {barWords(bar)}
      </span>
    </td>
  );
};

/** The same verdict in words, for prose: "VME 4.9%, misses at most 1.5%". */
export const verdictWords = (value: number, bar: ReferenceBar): string =>
  `${bar.name} ${pct(value)}, ${verdictOf(value, bar) === "meets" ? "meets" : "misses"} ${barWords(bar)}`;

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
}: {
  ariaLabel: string;
  options: readonly SegmentOption<T>[];
  value: T;
  onSelect(id: T): void;
}) => (
  <div role="group" aria-label={ariaLabel} className="inline-flex min-h-8 w-fit max-w-full flex-wrap items-center gap-0.5 rounded-md border border-hair bg-paper-pure p-0.5">
    {options.map((o) => {
      const selected = o.id === value;
      return (
        <button
          key={o.id}
          type="button"
          aria-pressed={selected}
          disabled={o.disabled === true}
          onClick={() => onSelect(o.id)}
          className={`ui-control-target inline-flex min-h-[26px] items-center justify-center whitespace-nowrap rounded px-3 text-control font-medium transition-[color,background-color,transform] duration-fast ease-out active:scale-[0.98] disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING} ${
            selected ? "bg-ink text-paper-pure" : "text-ink enabled:hover:bg-surface-subtle"
          }`}
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
  disabled = false,
  ariaLabel,
}: {
  children: ReactNode;
  onClick(): void;
  disabled?: boolean;
  ariaLabel?: string;
}) => (
  <button
    type="button"
    onClick={onClick}
    disabled={disabled}
    aria-label={ariaLabel}
    className={`ui-control-target inline-flex h-8 items-center gap-1.5 whitespace-nowrap rounded-md border border-hair px-3 text-control font-medium text-ink transition-[color,background-color,border-color,transform] duration-fast ease-out active:scale-[0.98] enabled:hover:bg-surface-subtle disabled:cursor-not-allowed disabled:opacity-45 ${FOCUS_RING}`}
  >
    {children}
  </button>
);

/** An inline text link that navigates inside the app. */
export const TextLink = ({ children, onClick }: { children: ReactNode; onClick(): void }) => (
  <button
    type="button"
    onClick={onClick}
    className={`rounded-sm text-left font-medium text-blue underline decoration-blue/30 underline-offset-2 transition-colors duration-fast ease-out hover:decoration-blue active:scale-[0.98] ${FOCUS_RING}`}
  >
    {children}
  </button>
);

// ---- Panel state: loading, empty, error ------------------------------------

export const PanelState = ({ kind, title, message, action }: { kind: "loading" | "empty" | "error"; title: string; message?: string; action?: ReactNode }) => (
  <div
    role={kind === "error" ? "alert" : "status"}
    aria-busy={kind === "loading" ? true : undefined}
    className={`rounded-xl border px-4 py-4 ${kind === "error" ? "border-coral/50 bg-coral/5" : "border-hair bg-surface-subtle"}`}
  >
    <p className={`text-body font-semibold ${kind === "error" ? "text-coral-text" : "text-ink"}`}>{title}</p>
    {message !== undefined ? <p className="type-note mt-1">{message}</p> : null}
    {action !== undefined ? <div className="mt-3">{action}</div> : null}
  </div>
);

// ---- Section: non-figure content ---------------------------------------------

/** A titled block: the answer sentence, then the evidence. Hairline
 *  separated, never a card inside a card. Figures use FigureBlock. */
export const Section = ({
  id,
  title,
  answer,
  right,
  children,
}: {
  id?: string;
  title: string;
  answer?: ReactNode;
  right?: ReactNode;
  children?: ReactNode;
}) => (
  <section id={id} aria-labelledby={id === undefined ? undefined : `${id}-title`} className="border-t border-hair pt-6">
    <header className="flex flex-wrap items-start justify-between gap-x-6 gap-y-2">
      <h2 id={id === undefined ? undefined : `${id}-title`} className="min-w-0 flex-1 basis-60 text-body font-semibold text-ink">
        {title}
      </h2>
      {right !== undefined ? <div className="flex min-w-0 flex-wrap items-center gap-2">{right}</div> : null}
    </header>
    {answer !== undefined ? <p className="type-headline mt-2">{answer}</p> : null}
    {children !== undefined ? <div className="mt-4">{children}</div> : null}
  </section>
);

// ---- Disclosure ----------------------------------------------------------------

export const Disclosure = ({
  label = "Method and coverage",
  subject,
  children,
  defaultOpen = false,
}: {
  label?: string;
  /** A second disclosure in one block names its subject after the label. */
  subject?: string;
  children: ReactNode;
  defaultOpen?: boolean;
}) => {
  const [open, setOpen] = useState(defaultOpen);
  const bodyId = useId();
  return (
    <div className="border-t border-hair pt-3">
      <button
        type="button"
        aria-expanded={open}
        aria-controls={bodyId}
        onClick={() => setOpen((v) => !v)}
        className={`ui-control-target -mx-1 inline-flex items-center gap-1.5 rounded-md px-1 py-0.5 text-control font-medium text-ink transition-colors duration-fast ease-out hover:bg-surface-subtle active:scale-[0.98] ${FOCUS_RING}`}
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

// ---- Figure block: title, result, figure, caption, data and method ----------

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
const FitBox = ({ figure }: { figure: FigureContent }) => {
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
      {fit !== null ? renderFigure(figure, fit) : null}
    </div>
  );
};

/** One analytical block in the figure grammar: the title names the measured
 *  quantity, the result gives the count and its denominator, the figure
 *  shows that number, the caption says what a mark is and the one limit that
 *  changes its reading, and Data and method holds the rest. On a slide the
 *  figure fills the height the text leaves. */
export const FigureBlock = ({
  id,
  title,
  result,
  figure,
  caption,
  details,
  detailsOpen = false,
  actions,
  variant = "page",
}: {
  id: string;
  title: string;
  result: ReactNode;
  figure: FigureContent;
  caption: ReactNode;
  details?: ReactNode;
  detailsOpen?: boolean;
  /** A control that changes this figure, beside its title. */
  actions?: ReactNode;
  variant?: "page" | "slide";
}) =>
  variant === "slide" ? (
    <section aria-labelledby={`${id}-title`} className="flex h-full min-h-0 flex-col">
      <h2 id={`${id}-title`} className="text-body font-semibold text-muted-text">
        {title}
      </h2>
      <p className="mt-2 max-w-[56rem] text-stat font-semibold text-ink">{result}</p>
      <div className="mt-4 flex min-h-0 flex-1 flex-col">
        <FitBox figure={figure} />
      </div>
      <p className="prose-measure mt-3 text-note text-muted-text">{caption}</p>
    </section>
  ) : (
    <section id={id} aria-labelledby={`${id}-title`} className="scroll-mt-6">
      <div className="flex flex-wrap items-center justify-between gap-x-6 gap-y-2">
        <h2 id={`${id}-title`} className="text-body font-semibold text-ink">
          {title}
        </h2>
        {actions !== undefined ? <div className="flex flex-wrap items-center gap-2">{actions}</div> : null}
      </div>
      <p className="type-headline mt-2">{result}</p>
      <div className="mt-4 min-w-0">{renderFigure(figure, {})}</div>
      <p className="prose-measure mt-3 text-note text-muted-text">{caption}</p>
      {details !== undefined ? (
        <div className="mt-4">
          <Disclosure label="Data and method" defaultOpen={detailsOpen}>
            {details}
          </Disclosure>
        </div>
      ) : null}
    </section>
  );

// ---- Fact rows: label left, value right, hairline between -------------------

export const FactRows = ({ rows }: { rows: { label: string; value: ReactNode }[] }) => (
  <dl className="m-0 max-w-xl">
    {rows.map((r) => (
      <div key={r.label} className="flex items-baseline justify-between gap-x-6 border-b border-hair py-1.5 last:border-b-0">
        <dt className="text-control font-medium text-muted-text">{r.label}</dt>
        <dd className="m-0 min-w-0 text-right text-control tabular-nums text-ink">{r.value}</dd>
      </div>
    ))}
  </dl>
);

/** A keyboard-reachable horizontal scroll owner for wide tables and plots. */
export const ScrollFrame = ({ children, label }: { children: ReactNode; label: string }) => (
  <div role="region" aria-label={label} tabIndex={0} className={`relative scroll-thin overflow-x-auto overscroll-x-contain ${FOCUS_RING}`}>
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
