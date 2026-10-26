// Whether a truncated line is cutting its text off, read after layout and
// again whenever its box changes size or the web font lands. One observer
// serves every watched element.

import { useCallback, useEffect, useRef, useState } from "react";

let observer: ResizeObserver | null = null;
const checks = new Map<Element, () => void>();

const watch = (el: Element, check: () => void): (() => void) => {
  if (typeof ResizeObserver === "undefined") return () => {};
  observer ??= new ResizeObserver((entries) => entries.forEach((e) => checks.get(e.target)?.()));
  checks.set(el, check);
  observer.observe(el);
  return () => {
    checks.delete(el);
    observer?.unobserve(el);
  };
};

export const isClipped = (el: HTMLElement): boolean => el.scrollWidth > el.clientWidth;

/** A ref for the truncated element and whether it is clipped now. While
 *  `paused` (the text is shown whole), readings are ignored and the last
 *  one holds, so the control that opened the text is still there to close
 *  it, with focus on it. */
export const useClipped = <T extends HTMLElement>(paused = false): [(el: T | null) => void, boolean] => {
  const [clipped, setClipped] = useState(false);
  const stop = useRef<(() => void) | null>(null);
  // Read by the observer's callback, so it sees this render's value.
  const hold = useRef(paused);
  hold.current = paused;
  const ref = useCallback((el: T | null) => {
    stop.current?.();
    stop.current = null;
    if (el === null) return;
    const check = () => {
      if (!hold.current) setClipped(isClipped(el));
    };
    check();
    const unwatch = watch(el, check);
    let live = true;
    void document.fonts?.ready.then(() => {
      if (live) check();
    });
    stop.current = () => {
      live = false;
      unwatch();
    };
  }, []);
  useEffect(() => () => stop.current?.(), []);
  return [ref, clipped];
};
