// Explicit details: a real button with aria-expanded and, right after it in
// the DOM, a small note of at most three lines. From 640px the note floats
// under the button without moving the text around it; below 640px, or
// inside a box that clips, it opens in the flow. One note is open at a
// time. Escape, or a tap anywhere else, closes it; that tap does nothing
// else. Focus stays on the button.

import { useEffect, useId, useLayoutEffect, useRef, useState, useSyncExternalStore, type ReactNode } from "react";
import { FOCUS_RING } from "./primitives";

// ---- The one open note ------------------------------------------------------

let openId: string | null = null;
const subscribers = new Set<() => void>();

const setOpenId = (id: string | null): void => {
  openId = id;
  subscribers.forEach((s) => s());
};

const subscribe = (s: () => void) => {
  subscribers.add(s);
  return () => {
    subscribers.delete(s);
  };
};

const readOpen = (): string | null => openId;

/** Closes whichever note is open; for a caller that changes the page under
 *  it (a filter, a new view). */
export const closeDetails = (): void => setOpenId(null);

// A tap outside an open note closes it and is spent: the click it would
// make on whatever lies under it is swallowed.
const CLICK_AFTER_TAP_MS = 700;
let swallowUntil = 0;
let swallowing = false;

const swallowClick = (e: MouseEvent): void => {
  if (performance.now() > swallowUntil) return;
  swallowUntil = 0;
  e.preventDefault();
  e.stopPropagation();
};

const listenForSpentClick = (): void => {
  if (swallowing) return;
  swallowing = true;
  document.addEventListener("click", swallowClick, true);
};

// ---- Layout -----------------------------------------------------------------

const NARROW = "(max-width: 639.98px)";

const subscribeNarrow = (cb: () => void) => {
  const q = window.matchMedia(NARROW);
  q.addEventListener("change", cb);
  return () => q.removeEventListener("change", cb);
};

const readNarrow = (): boolean => window.matchMedia(NARROW).matches;

const CLIPS = new Set(["hidden", "clip", "auto", "scroll"]);

/** A box between the button and the view's scroll owner that would cut a
 *  floating note off. */
const insideClip = (el: Element): boolean => {
  for (let n = el.parentElement; n !== null && n !== document.body && !n.hasAttribute("data-view-scroll"); n = n.parentElement) {
    const s = getComputedStyle(n);
    if (CLIPS.has(s.overflowX) || CLIPS.has(s.overflowY)) return true;
  }
  return false;
};

export const DETAILS_BUTTON = `ui-control-target ui-target-24 inline-flex items-center gap-1 rounded-sm text-control font-medium text-blue underline decoration-blue/40 underline-offset-2 transition-colors duration-fast ease-out hover:text-blue-deep hover:decoration-current active:scale-[0.98] ${FOCUS_RING}`;

const NOTE = "block max-w-[min(22rem,calc(100vw-1rem))] rounded-xl border border-hair bg-paper-pure px-3 py-2 text-left text-note font-normal normal-case tracking-normal text-ink shadow-[var(--surface-shadow)]";

export const DetailsButton = ({
  label,
  children,
  className,
}: {
  /** The button's words: what the note explains ("How this is counted"). */
  label: ReactNode;
  /** At most three short lines. */
  children: ReactNode;
  /** Replaces the quiet blue text-button look. */
  className?: string;
}) => {
  const id = useId();
  const noteId = `${id}-note`;
  const open = useSyncExternalStore(subscribe, readOpen, () => null) === id;
  const narrow = useSyncExternalStore(subscribeNarrow, readNarrow, () => false);
  const rootRef = useRef<HTMLSpanElement | null>(null);
  const buttonRef = useRef<HTMLButtonElement | null>(null);
  const noteRef = useRef<HTMLSpanElement | null>(null);
  const [clipped, setClipped] = useState(false);
  const [shift, setShift] = useState(0);
  const inline = narrow || clipped;

  useEffect(() => {
    if (!open) return undefined;
    listenForSpentClick();
    const onDown = (e: PointerEvent) => {
      if (rootRef.current?.contains(e.target as Node) === true) return;
      swallowUntil = performance.now() + CLICK_AFTER_TAP_MS;
      setOpenId(null);
    };
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "Escape") return;
      // Consumed, so a dialog under the note stays open.
      e.preventDefault();
      e.stopPropagation();
      const within = noteRef.current?.contains(document.activeElement) === true;
      setOpenId(null);
      if (within) buttonRef.current?.focus();
    };
    document.addEventListener("pointerdown", onDown, true);
    window.addEventListener("keydown", onKey, true);
    return () => {
      document.removeEventListener("pointerdown", onDown, true);
      window.removeEventListener("keydown", onKey, true);
    };
  }, [open]);

  // A note whose button goes away goes with it.
  useEffect(
    () => () => {
      if (openId === id) setOpenId(null);
    },
    [id],
  );

  // A floating note stays inside the viewport: shifted left from the right
  // edge, never past the left one.
  useLayoutEffect(() => {
    if (!open || inline || noteRef.current === null) {
      setShift(0);
      return;
    }
    const r = noteRef.current.getBoundingClientRect();
    const over = r.right - (document.documentElement.clientWidth - 8);
    setShift(over > 0 ? -Math.min(over, Math.max(0, r.left - 8)) : 0);
  }, [open, inline]);

  return (
    <span ref={rootRef} data-details="" className={inline ? "" : "relative inline-block"}>
      <button
        ref={buttonRef}
        type="button"
        aria-expanded={open}
        aria-controls={open ? noteId : undefined}
        onClick={() => {
          if (!open && buttonRef.current !== null) setClipped(insideClip(buttonRef.current));
          setOpenId(open ? null : id);
        }}
        className={className ?? DETAILS_BUTTON}
      >
        {label}
      </button>
      {open ? (
        <span
          ref={noteRef}
          id={noteId}
          role="note"
          data-placement={inline ? "inline" : "float"}
          className={`${NOTE} ${inline ? "mt-1.5 w-full" : "motion-popover-in absolute left-0 top-full z-[55] mt-1.5 w-max"}`}
          style={shift !== 0 ? { left: shift } : undefined}
        >
          {children}
        </span>
      ) : null}
    </span>
  );
};
