// The staged-reveal toolkit for deck slides. One useSlideTimeline per slide
// walks its stages on timers; Reveal and AnimatedNumber wait for their stage.
// Advancing mid-animation calls skipToEnd() so nothing makes the viewer wait,
// and reduced motion renders every final state at once.

import { useCallback, useEffect, useRef, useState, type CSSProperties, type ReactNode } from "react";

const REDUCED_QUERY = "(prefers-reduced-motion: reduce)";

const readReduced = (): boolean => typeof window !== "undefined" && window.matchMedia(REDUCED_QUERY).matches;

export const usePrefersReducedMotion = (): boolean => {
  const [reduced, setReduced] = useState(readReduced);
  useEffect(() => {
    const mq = window.matchMedia(REDUCED_QUERY);
    const onChange = (e: MediaQueryListEvent) => setReduced(e.matches);
    mq.addEventListener("change", onChange);
    return () => mq.removeEventListener("change", onChange);
  }, []);
  return reduced;
};

export interface SlideTimeline {
  /** Current stage, 0..steps. A Reveal at={n} shows once stage >= n. */
  stage: number;
  steps: number;
  /** True after a skip: staged elements jump to their final state. */
  skipped: boolean;
  /** Completes every stage; false when already done, so the caller advances. */
  skipToEnd(): boolean;
}

/** Default rhythm: the figure starts at 60ms, later beats follow at 900ms. */
const defaultGap = (stage: number): number => (stage === 1 ? 60 : 900);

export const useSlideTimeline = (steps: number, delaysMs?: readonly number[]): SlideTimeline => {
  const reduced = usePrefersReducedMotion();
  const [stage, setStage] = useState(() => (readReduced() ? steps : 0));
  const [skipped, setSkipped] = useState(false);
  const stageRef = useRef(stage);
  const timers = useRef<number[]>([]);
  const key = JSON.stringify(delaysMs ?? []);

  useEffect(() => {
    if (reduced) {
      stageRef.current = steps;
      setStage(steps);
      return undefined;
    }
    const delays = JSON.parse(key) as number[];
    let at = 0;
    for (let s = 1; s <= steps; s += 1) {
      at += delays[s - 1] ?? defaultGap(s);
      const target = s;
      timers.current.push(
        window.setTimeout(() => {
          if (stageRef.current < target) {
            stageRef.current = target;
            setStage(target);
          }
        }, at),
      );
    }
    return () => {
      timers.current.forEach((t) => window.clearTimeout(t));
      timers.current = [];
    };
  }, [steps, reduced, key]);

  const skipToEnd = useCallback((): boolean => {
    timers.current.forEach((t) => window.clearTimeout(t));
    timers.current = [];
    if (stageRef.current >= steps) return false;
    stageRef.current = steps;
    setStage(steps);
    setSkipped(true);
    return true;
  }, [steps]);

  return { stage, steps, skipped, skipToEnd };
};

export type RevealKind = "fade-up" | "grow";

export const Reveal = ({
  at,
  timeline,
  kind = "fade-up",
  durationMs = 600,
  delayMs = 0,
  className,
  style,
  children,
}: {
  at: number;
  timeline: SlideTimeline;
  kind?: RevealKind;
  durationMs?: number;
  delayMs?: number;
  className?: string;
  style?: CSSProperties;
  children: ReactNode;
}) => {
  const reduced = usePrefersReducedMotion();
  const active = timeline.stage >= at;
  const instant = reduced || timeline.skipped;
  const from = kind === "grow" ? "scaleX(0.001)" : "translateY(var(--motion-shift-surface))";
  return (
    <div
      className={className}
      style={{
        ...style,
        opacity: active ? 1 : 0,
        transform: active ? "none" : from,
        transformOrigin: kind === "grow" ? "left center" : undefined,
        transition: !active || instant ? "none" : `opacity ${durationMs}ms var(--ease-out) ${delayMs}ms, transform ${durationMs}ms var(--ease-out) ${delayMs}ms`,
      }}
    >
      {children}
    </div>
  );
};

/** A count-up that waits for its stage; a skip snaps it to the value. */
export const AnimatedNumber = ({
  value,
  format,
  timeline,
  at,
  durationMs = 900,
}: {
  value: number;
  format(n: number): string;
  timeline: SlideTimeline;
  at: number;
  durationMs?: number;
}) => {
  const reduced = usePrefersReducedMotion();
  const active = timeline.stage >= at;
  const [shown, setShown] = useState(active && (reduced || timeline.skipped) ? value : 0);
  useEffect(() => {
    if (!active) return undefined;
    if (reduced || timeline.skipped) {
      setShown(value);
      return undefined;
    }
    const start = performance.now();
    let frame = 0;
    const tick = (now: number) => {
      const t = Math.min((now - start) / durationMs, 1);
      setShown(value * (1 - (1 - t) ** 3));
      if (t < 1) frame = requestAnimationFrame(tick);
    };
    frame = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(frame);
  }, [active, reduced, timeline.skipped, value, durationMs]);
  return <span>{format(active ? shown : 0)}</span>;
};
