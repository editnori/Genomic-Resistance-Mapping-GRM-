// A press and hold on a touch screen does what a tap does. The hold is a
// touch held about half a second without moving more than a few pixels; a
// move, a lift or a scroll (pointercancel) before then is no hold. The click
// a browser may send when the finger lifts after a hold is swallowed at the
// document, whatever it lands on, so the hold neither undoes itself nor,
// once the rule has shrunk back under the finger, opens the row behind it.
// A click from the keyboard (detail 0) is never swallowed. Mouse and pen
// presses are left alone: for them a click is the way in. A hold target
// should also be HOLD_TARGET, so a long press selects no text.

import { useEffect, useRef, type MouseEvent, type PointerEvent } from "react";

// Under Chromium's own touch long press (500 ms), which would otherwise
// start a text selection and cancel the pointer first.
export const HOLD_MS = 400;
export const HOLD_SLOP_PX = 10;
/** No text selection or callout on a long press. */
export const HOLD_TARGET = "select-none [-webkit-touch-callout:none]";
/** How long after the lift a hold's click is still swallowed. */
export const LIFT_CLICK_MS = 800;

/** Handlers for the element a tap acts on: `act` runs on a click or tap and
 *  on a touch hold, once either way. */
export const useHold = (act: () => void) => {
  const latest = useRef(act);
  latest.current = act;
  const timer = useRef<ReturnType<typeof setTimeout> | undefined>(undefined);
  const start = useRef<{ x: number; y: number } | null>(null);
  const held = useRef(false);
  const unswallow = useRef<(() => void) | null>(null);
  const cancel = () => {
    clearTimeout(timer.current);
    timer.current = undefined;
    start.current = null;
  };
  const stopSwallowing = () => {
    unswallow.current?.();
    unswallow.current = null;
  };
  // Until the click after the lift, the lift's window closing, or the next
  // press, whichever comes first.
  const swallowLiftClick = () => {
    stopSwallowing();
    let lift: ReturnType<typeof setTimeout> | undefined;
    const click = (e: globalThis.MouseEvent) => {
      if (e.detail === 0) return;
      e.preventDefault();
      e.stopPropagation();
      held.current = false;
      stopSwallowing();
    };
    const up = () => {
      clearTimeout(lift);
      lift = setTimeout(stopSwallowing, LIFT_CLICK_MS);
    };
    document.addEventListener("click", click, true);
    document.addEventListener("pointerup", up, true);
    document.addEventListener("pointerdown", stopSwallowing, true);
    unswallow.current = () => {
      clearTimeout(lift);
      document.removeEventListener("click", click, true);
      document.removeEventListener("pointerup", up, true);
      document.removeEventListener("pointerdown", stopSwallowing, true);
    };
  };
  useEffect(
    () => () => {
      cancel();
      stopSwallowing();
    },
    [],
  );
  return {
    onPointerDown: (e: PointerEvent) => {
      cancel();
      held.current = false;
      if (e.pointerType !== "touch") return;
      start.current = { x: e.clientX, y: e.clientY };
      timer.current = setTimeout(() => {
        cancel();
        held.current = true;
        swallowLiftClick();
        latest.current();
      }, HOLD_MS);
    },
    onPointerMove: (e: PointerEvent) => {
      const s = start.current;
      if (s !== null && Math.hypot(e.clientX - s.x, e.clientY - s.y) > HOLD_SLOP_PX) cancel();
    },
    onPointerUp: cancel,
    onPointerCancel: cancel,
    // A long press would otherwise open the browser's own menu over the rule.
    onContextMenu: (e: MouseEvent) => {
      if (held.current || start.current !== null) e.preventDefault();
    },
    onClick: (e: MouseEvent) => {
      if (held.current && e.detail !== 0) {
        held.current = false;
        return;
      }
      held.current = false;
      latest.current();
    },
  };
};
