// The GRM mark: an isometric cube read as a genome. Each side face is a
// stack of three k-mer tiles; the two strands swap colour rung by rung, so
// the faces read as a double helix turning, and the top face is the call.
// The motion lives in styles/index.css under "The GRM mark"; this component
// only draws the geometry and names the motion. At rest it is static.

import { useEffect, useState, type CSSProperties } from "react";

export type MarkMotion = "still" | "assemble" | "loop" | "replay" | "present";

const TOP = "12,2.6 20.4,7.3 12,12 3.6,7.3";
const RUNG = 9.4 / 3;
const at = (v: number): string => String(Math.round(v * 100) / 100);
const leftFace = (k: number): string =>
  `3.6,${at(7.3 + RUNG * k)} 12,${at(12 + RUNG * k)} 12,${at(12 + RUNG * (k + 1))} 3.6,${at(7.3 + RUNG * (k + 1))}`;
const rightFace = (k: number): string =>
  `12,${at(12 + RUNG * k)} 20.4,${at(7.3 + RUNG * k)} 20.4,${at(7.3 + RUNG * (k + 1))} 12,${at(12 + RUNG * (k + 1))}`;

/** Rungs bottom first, so the build and the turn both start at the base. */
const RUNGS = [2, 1, 0] as const;

export const GrmMark = ({ size = 24, motion = "still", className = "" }: { size?: number; motion?: MarkMotion; className?: string }) => (
  <svg
    viewBox="0 0 24 24"
    width={size}
    height={size}
    aria-hidden="true"
    focusable="false"
    data-motion={motion}
    className={`grm-mark shrink-0 ${className}`}
  >
    {RUNGS.map((k, i) => (
      <g key={k} className="grm-mark-rung" style={{ "--i": i } as CSSProperties}>
        <polygon className="grm-mark-face" data-strand={k % 2 === 0 ? "a" : "b"} points={leftFace(k)} />
        <polygon className="grm-mark-face" data-strand={k % 2 === 0 ? "b" : "a"} points={rightFace(k)} />
      </g>
    ))}
    <polygon className="grm-mark-top" points={TOP} />
  </svg>
);

/** How long each one-shot motion runs before the mark is at rest. */
const ONE_SHOT_MS: Partial<Record<MarkMotion, number>> = { assemble: 1100, replay: 600, present: 2400 };

/** The brand mark's motion: the loading loop while the report loads, one
 *  assembly when it arrives, then rest; replay() runs the short replay from
 *  rest, and a replay already running is left to finish. */
export const useMarkMotion = (loading: boolean): { motion: MarkMotion; replay(): void } => {
  const [motion, setMotion] = useState<MarkMotion>(loading ? "loop" : "assemble");
  useEffect(() => {
    if (!loading && motion === "loop") setMotion("assemble");
  }, [loading, motion]);
  useEffect(() => {
    const ms = ONE_SHOT_MS[motion];
    if (ms === undefined) return undefined;
    const timer = window.setTimeout(() => setMotion("still"), ms);
    return () => window.clearTimeout(timer);
  }, [motion]);
  return { motion, replay: () => setMotion((m) => (m === "still" ? "replay" : m)) };
};
