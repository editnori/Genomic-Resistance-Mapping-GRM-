// The one hover recipe: an 11.5px/600 term, a 14px explanation, and a 12px
// source line on a paper sheet. Opens on hover or keyboard focus, closes on
// leave, blur, or Escape. No native title attribute competes with it.

import { useCallback, useEffect, useId, useRef, useState, type ReactNode } from "react";
import { createPortal } from "react-dom";

export interface TipLine {
  term?: string;
  text?: string;
  source?: string;
}

const MAX_WIDTH = 320;
const MARGIN = 8;
const GAP = 6;
const SHOW_DELAY_MS = 200;

interface Placement {
  left: number;
  top: number | null;
  bottom: number | null;
  host: HTMLElement;
}

const place = (anchor: Element): Placement => {
  const rect = anchor.getBoundingClientRect();
  const vw = document.documentElement.clientWidth;
  const vh = document.documentElement.clientHeight;
  const width = Math.min(MAX_WIDTH, vw - MARGIN * 2);
  const left = Math.max(MARGIN, Math.min(rect.left, vw - width - MARGIN));
  const host = anchor.closest<HTMLElement>("dialog[open]") ?? document.body;
  return rect.bottom + GAP + 120 > vh
    ? { left, top: null, bottom: Math.max(MARGIN, vh - rect.top + GAP), host }
    : { left, top: rect.bottom + GAP, bottom: null, host };
};

export const TipBody = ({ lines }: { lines: readonly TipLine[] }) => (
  <div className="space-y-1">
    {lines.map((l, i) => (
      <div key={i}>
        {l.term !== undefined ? <div className="text-tiny font-semibold text-ink">{l.term}</div> : null}
        {l.text !== undefined ? <div className="text-sm leading-snug text-ink">{l.text}</div> : null}
        {l.source !== undefined ? <div className="text-[12px] leading-snug text-muted-text">{l.source}</div> : null}
      </div>
    ))}
  </div>
);

export const TIP_SURFACE =
  "pointer-events-none rounded-xl border border-hair bg-paper-pure px-3 py-2 text-left font-normal normal-case tracking-normal text-ink shadow-[var(--surface-shadow)]";

export const HoverTip = ({
  lines,
  children,
  focusable = true,
  block = false,
}: {
  lines: readonly TipLine[];
  children: ReactNode;
  /** False when the child is already a tab stop; the wrapper then listens to
   *  its focus instead of adding a second stop. */
  focusable?: boolean;
  block?: boolean;
}) => {
  const id = useId();
  const ref = useRef<HTMLElement | null>(null);
  const timer = useRef<number | null>(null);
  const [at, setAt] = useState<Placement | null>(null);

  const hide = useCallback(() => {
    if (timer.current !== null) window.clearTimeout(timer.current);
    timer.current = null;
    setAt(null);
  }, []);

  const show = useCallback((delay: number) => {
    if (timer.current !== null) window.clearTimeout(timer.current);
    timer.current = window.setTimeout(() => {
      if (ref.current !== null) setAt(place(ref.current));
    }, delay);
  }, []);

  useEffect(() => {
    if (at === null) return undefined;
    // An open tip owns Escape: it closes the tip and nothing behind it, so
    // a dialog under the tip stays open.
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "Escape") return;
      e.preventDefault();
      e.stopPropagation();
      hide();
    };
    window.addEventListener("keydown", onKey, true);
    window.addEventListener("scroll", hide, true);
    return () => {
      window.removeEventListener("keydown", onKey, true);
      window.removeEventListener("scroll", hide, true);
    };
  }, [at, hide]);

  useEffect(() => hide, [hide]);

  const Tag = block ? "div" : "span";
  return (
    <Tag
      ref={(el: HTMLElement | null) => {
        ref.current = el;
      }}
      tabIndex={focusable ? 0 : undefined}
      aria-describedby={at !== null ? id : undefined}
      onPointerEnter={() => show(SHOW_DELAY_MS)}
      onPointerLeave={hide}
      onFocus={() => show(0)}
      onBlur={hide}
      className={`${block ? "block" : "inline-flex"} max-w-full rounded-sm focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-blue`}
    >
      {children}
      {at !== null
        ? createPortal(
            <div
              id={id}
              role="tooltip"
              className={`motion-popover-in fixed z-[60] ${TIP_SURFACE}`}
              style={{ left: at.left, top: at.top ?? undefined, bottom: at.bottom ?? undefined, maxWidth: MAX_WIDTH }}
            >
              <TipBody lines={lines} />
            </div>,
            at.host,
          )
        : null}
    </Tag>
  );
};
