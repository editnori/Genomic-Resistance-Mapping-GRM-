// Clinical agreement computed in the browser (the panel designer and the
// unexplained-resistance counts). Definitions match the pipeline's
// BinaryMetrics: VME = resistant called susceptible over all resistant;
// ME = susceptible called resistant over all susceptible.

export interface Agreement {
  n: number;
  tp: number;
  fp: number;
  tn: number;
  fn: number;
  sensitivity: number;
  specificity: number;
  balanced_accuracy: number;
  vme: number;
  me: number;
}

export const agreement = (truth: readonly boolean[], calls: readonly boolean[]): Agreement => {
  let tp = 0;
  let fp = 0;
  let tn = 0;
  let fn = 0;
  truth.forEach((resistant, i) => {
    const called = calls[i] === true;
    if (resistant) {
      if (called) tp += 1;
      else fn += 1;
    } else if (called) fp += 1;
    else tn += 1;
  });
  const sensitivity = tp + fn === 0 ? 0 : tp / (tp + fn);
  const specificity = tn + fp === 0 ? 0 : tn / (tn + fp);
  return {
    n: truth.length,
    tp,
    fp,
    tn,
    fn,
    sensitivity,
    specificity,
    balanced_accuracy: (sensitivity + specificity) / 2,
    vme: 1 - sensitivity,
    me: 1 - specificity,
  };
};

/** A marker-set test: resistant when any (or all) chosen markers are present. */
export const markerCalls = (
  markers: readonly (readonly boolean[])[],
  n: number,
  rule: "any" | "all",
): boolean[] => {
  if (markers.length === 0) return new Array<boolean>(n).fill(false);
  return Array.from({ length: n }, (_, i) =>
    rule === "any" ? markers.some((m) => m[i] === true) : markers.every((m) => m[i] === true),
  );
};

/** FDA-style reference bars for susceptibility test devices. These are the
 *  bars the dashboard judges VME and ME (at most) and MIC essential agreement
 *  (at least) against; they are a reference, not a regulatory claim about
 *  these models. */
export const REFERENCE_BARS = { vme: 0.015, me: 0.03, mic: 0.9 } as const;

export type BarVerdict = "meets" | "misses";

/** An error rate meets its bar at or below it. */
export const barVerdict = (value: number, bar: number): BarVerdict => (value <= bar ? "meets" : "misses");

export interface ReferenceBar {
  /** The measure, as a column header names it. */
  name: string;
  value: number;
  /** Errors are capped from above; agreement is a floor. */
  direction: "at most" | "at least";
}

export const BARS = {
  vme: { name: "VME", value: REFERENCE_BARS.vme, direction: "at most" },
  me: { name: "ME", value: REFERENCE_BARS.me, direction: "at most" },
  mic: { name: "MIC agreement", value: REFERENCE_BARS.mic, direction: "at least" },
} as const satisfies Record<string, ReferenceBar>;

export const verdictOf = (value: number, bar: ReferenceBar): BarVerdict =>
  bar.direction === "at most" ? barVerdict(value, bar.value) : value >= bar.value ? "meets" : "misses";

/** The bar as a header states it: "≤1.5%", "≤3%", "≥90%". */
export const barLabel = (bar: ReferenceBar): string => `${bar.direction === "at most" ? "≤" : "≥"}${Number((bar.value * 100).toFixed(1))}%`;

/** The bar in words for a screen reader or a sentence: "at most 1.5%". */
export const barWords = (bar: ReferenceBar): string => `${bar.direction} ${Number((bar.value * 100).toFixed(1))}%`;
