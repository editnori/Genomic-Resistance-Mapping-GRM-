import { useEffect, useRef, useState } from "react";

const reduced = (): boolean => window.matchMedia("(prefers-reduced-motion: reduce)").matches;

/** The value at its value on first paint; a later change (a filter or a
 *  marker choice) counts from the shown number to the new one. */
export const useCountUp = (value: number, durationMs = 320): number => {
  const [shown, setShown] = useState(value);
  const from = useRef(value);
  useEffect(() => {
    if (from.current === value) return undefined;
    if (reduced()) {
      from.current = value;
      setShown(value);
      return undefined;
    }
    const start = performance.now();
    const origin = from.current;
    let frame = 0;
    const tick = (now: number) => {
      const t = Math.min((now - start) / durationMs, 1);
      const eased = 1 - (1 - t) ** 3;
      const v = origin + (value - origin) * eased;
      from.current = v;
      setShown(v);
      if (t < 1) frame = requestAnimationFrame(tick);
      else from.current = value;
    };
    frame = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(frame);
  }, [value, durationMs]);
  return shown;
};
