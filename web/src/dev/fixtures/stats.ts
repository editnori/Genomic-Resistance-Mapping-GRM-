// Ports of the pipeline's metric code (src/learn/mod.rs, src/learn/eval.rs)
// so fixture numbers are computed the same way the real reports are.

import type { BinaryMetrics, ConformalPoint } from "../../data/schema";

export const rng = (seed: number): (() => number) => {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
};

export const normal = (rand: () => number): number => {
  const u = Math.max(rand(), 1e-12);
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rand());
};

const ratio = (a: number, b: number): number => (b === 0 ? 0 : a / b);

const tiedBlocks = (idx: number[], scores: readonly number[], visit: (block: number[]) => void): void => {
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && scores[idx[j + 1]!] === scores[idx[i]!]) j += 1;
    visit(idx.slice(i, j + 1));
    i = j + 1;
  }
};

export const auroc = (truth: readonly boolean[], scores: readonly number[]): number | null => {
  const idx = scores.map((_, i) => i).sort((a, b) => scores[a]! - scores[b]!);
  const nPos = truth.filter(Boolean).length;
  const nNeg = truth.length - nPos;
  if (nPos === 0 || nNeg === 0) return null;
  let rankSum = 0;
  let start = 0;
  tiedBlocks(idx, scores, (block) => {
    const avgRank = (start + start + block.length - 1) / 2 + 1;
    rankSum += block.filter((k) => truth[k]).length * avgRank;
    start += block.length;
  });
  return (rankSum - (nPos * (nPos + 1)) / 2) / (nPos * nNeg);
};

export const auprc = (truth: readonly boolean[], scores: readonly number[]): number | null => {
  const nPos = truth.filter(Boolean).length;
  if (nPos === 0) return null;
  const idx = scores.map((_, i) => i).sort((a, b) => scores[b]! - scores[a]!);
  let tp = 0;
  let seen = 0;
  let ap = 0;
  tiedBlocks(idx, scores, (block) => {
    const blockTp = block.filter((k) => truth[k]).length;
    seen += block.length;
    if (blockTp > 0) {
      tp += blockTp;
      ap += (blockTp / nPos) * (tp / seen);
    }
  });
  return ap;
};

export const binaryMetrics = (truth: readonly boolean[], scores: readonly number[], ranked: boolean): BinaryMetrics => {
  let tp = 0;
  let fp = 0;
  let tn = 0;
  let fn = 0;
  truth.forEach((t, i) => {
    const call = (scores[i] ?? 0) >= 0.5;
    if (t && call) tp += 1;
    else if (!t && call) fp += 1;
    else if (!t) tn += 1;
    else fn += 1;
  });
  const sensitivity = ratio(tp, tp + fn);
  const specificity = ratio(tn, tn + fp);
  const precision = ratio(tp, tp + fp);
  const denom = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return {
    n: truth.length,
    tp,
    fp,
    tn,
    fn,
    accuracy: ratio(tp + tn, truth.length),
    balanced_accuracy: (sensitivity + specificity) / 2,
    sensitivity,
    specificity,
    precision,
    f1: precision + sensitivity === 0 ? 0 : (2 * precision * sensitivity) / (precision + sensitivity),
    mcc: denom === 0 ? 0 : (tp * tn - fp * fn) / denom,
    vme: ratio(fn, tp + fn),
    me: ratio(fp, fp + tn),
    auroc: ranked ? auroc(truth, scores) : null,
    auprc: ranked ? auprc(truth, scores) : null,
  };
};

export const rocCurve = (truth: readonly boolean[], scores: readonly number[], maxPoints = 60): [number, number][] => {
  const nPos = Math.max(truth.filter(Boolean).length, 1);
  const nNeg = Math.max(truth.filter((t) => !t).length, 1);
  const idx = scores.map((_, i) => i).sort((a, b) => scores[b]! - scores[a]!);
  let pts: [number, number][] = [[0, 0]];
  let tp = 0;
  let fp = 0;
  tiedBlocks(idx, scores, (block) => {
    block.forEach((k) => {
      if (truth[k]) tp += 1;
      else fp += 1;
    });
    pts.push([fp / nNeg, tp / nPos]);
  });
  if (pts.length > maxPoints) {
    const step = pts.length / maxPoints;
    const last = pts[pts.length - 1]!;
    pts = Array.from({ length: maxPoints }, (_, i) => pts[Math.floor(i * step)]!);
    pts.push(last);
  }
  return pts;
};

export const ALPHAS = [0.01, 0.025, 0.05, 0.1, 0.2] as const;

/** Mondrian split-conformal, calibrated on one half and measured on the
 *  other, both ways round (src/learn/eval.rs `conformal`). */
export const conformal = (truth: readonly boolean[], scores: readonly number[], halves: readonly number[]): ConformalPoint[] =>
  ALPHAS.map((alpha) => {
    let single = 0;
    let amb = 0;
    let fnrN = 0;
    let fnrD = 0;
    let fprN = 0;
    let fprD = 0;
    for (const half of [0, 1]) {
      const cal = (cls: boolean): number => {
        const nc = truth
          .map((t, i) => ({ t, i }))
          .filter(({ t, i }) => halves[i] === half && t === cls)
          .map(({ i }) => (cls ? 1 - scores[i]! : scores[i]!))
          .sort((a, b) => a - b);
        if (nc.length === 0) return 1;
        const rank = Math.ceil((nc.length + 1) * (1 - alpha));
        return rank > nc.length ? Infinity : nc[rank - 1]!;
      };
      const qR = cal(true);
      const qS = cal(false);
      truth.forEach((t, i) => {
        if (halves[i] === half) return;
        const hasR = 1 - scores[i]! <= qR;
        const hasS = scores[i]! <= qS;
        if (hasR === hasS) amb += 1;
        else {
          single += 1;
          if (hasR && !t) fprN += 1;
          if (hasS && t) fnrN += 1;
        }
        if (t) fnrD += 1;
        else fprD += 1;
      });
    }
    const n = truth.length;
    return { alpha, coverage: single / n, ambiguous: amb / n, vme: fnrN / Math.max(fnrD, 1), me: fprN / Math.max(fprD, 1) };
  });
