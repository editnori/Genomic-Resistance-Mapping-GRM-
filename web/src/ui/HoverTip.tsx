// The one optional hint: an 11.5px/600 term, a 14px explanation, and a 12px
// source line on a paper sheet, at most about 40 words. It never holds
// anything the reader needs: definitions, verdicts and keys a reading
// depends on are printed. A mouse or pen resting on the wrapped content
// opens it, as does keyboard focus on a control inside; it adds no tab stop
// of its own. It is measured into the viewport, takes no pointer events,
// and closes on leave, blur, scroll, or Escape. No native title attribute
// competes with it.

import { useCallback, useEffect, useId, useLayoutEffect, useRef, useState, type ReactNode } from "react";
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

/** Where the tip goes once its own size is known: under the anchor, else
 *  over it, else wherever it fits, always inside the viewport. */
export const placeTip = (
  anchor: { left: number; top: number; bottom: number },
  tip: { width: number; height: number },
  view: { width: number; height: number },
): { left: number; top: number } => {
  const left = Math.max(MARGIN, Math.min(anchor.left, view.width - tip.width - MARGIN));
  const below = anchor.bottom + GAP;
  const above = anchor.top - GAP - tip.height;
  const top =
    below + tip.height <= view.height - MARGIN ? below : above >= MARGIN ? above : Math.max(MARGIN, Math.min(below, view.height - tip.height - MARGIN));
  return { left, top };
};

/** Keyboard focus opens a tip; focus that a click or a dialog's opening
 *  moves does not. */
const keyboardFocus = (el: EventTarget): boolean => {
  try {
    return el instanceof Element && el.matches(":focus-visible");
  } catch {
    return true;
  }
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

interface Opened {
  anchor: { left: number; top: number; bottom: number };
  host: HTMLElement;
}

export const HoverTip = ({
  lines,
  children,
  focusable = false,
  block = false,
}: {
  lines: readonly TipLine[];
  children: ReactNode;
  /** Deprecated: a hint never needs a tab stop of its own. True still adds
   *  one for callers that have not moved their content inline yet. */
  focusable?: boolean;
  block?: boolean;
}) => {
  const id = useId();
  const ref = useRef<HTMLElement | null>(null);
  const tipRef = useRef<HTMLDivElement | null>(null);
  const timer = useRef<number | null>(null);
  const [opened, setOpened] = useState<Opened | null>(null);
  const [at, setAt] = useState<{ left: number; top: number } | null>(null);

  const hide = useCallback(() => {
    if (timer.current !== null) window.clearTimeout(timer.current);
    timer.current = null;
    setOpened(null);
    setAt(null);
  }, []);

  const show = useCallback((delay: number) => {
    if (timer.current !== null) window.clearTimeout(timer.current);
    timer.current = window.setTimeout(() => {
      const el = ref.current;
      if (el === null) return;
      const r = el.getBoundingClientRect();
      setOpened({ anchor: { left: r.left, top: r.top, bottom: r.bottom }, host: el.closest<HTMLElement>("dialog[open]") ?? document.body });
    }, delay);
  }, []);

  // Measured before paint, so the tip never shows outside the viewport.
  useLayoutEffect(() => {
    if (opened === null || tipRef.current === null) return;
    const t = tipRef.current.getBoundingClientRect();
    setAt(placeTip(opened.anchor, { width: t.width, height: t.height }, { width: document.documentElement.clientWidth, height: document.documentElement.clientHeight }));
  }, [opened]);

  useEffect(() => {
    if (opened === null) return undefined;
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
  }, [opened, hide]);

  useEffect(() => hide, [hide]);

  const Tag = block ? "div" : "span";
  return (
    <Tag
      ref={(el: HTMLElement | null) => {
        ref.current = el;
      }}
      tabIndex={focusable ? 0 : undefined}
      aria-describedby={opened !== null ? id : undefined}
      onPointerEnter={(e) => {
        // A touch has no resting pointer; a tap means the control's own
        // action, never a hint.
        if (e.pointerType === "mouse" || e.pointerType === "pen") show(SHOW_DELAY_MS);
      }}
      onPointerLeave={hide}
      onFocus={(e) => {
        if (keyboardFocus(e.target)) show(0);
      }}
      onBlur={hide}
      className={`${block ? "block" : "inline-flex"} max-w-full rounded-sm focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-blue`}
    >
      {children}
      {opened !== null
        ? createPortal(
            <div
              ref={tipRef}
              id={id}
              role="tooltip"
              className={`motion-popover-in fixed z-[60] ${TIP_SURFACE}`}
              style={{ left: at?.left ?? 0, top: at?.top ?? 0, maxWidth: `min(${MAX_WIDTH}px, calc(100vw - ${MARGIN * 2}px))`, visibility: at === null ? "hidden" : undefined }}
            >
              <TipBody lines={lines} />
            </div>,
            opened.host,
          )
        : null}
    </Tag>
  );
};
