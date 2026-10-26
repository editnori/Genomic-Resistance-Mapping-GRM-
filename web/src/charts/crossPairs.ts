// Every drug pair of a species, for the cross-resistance matrix and its
// table. The species file lists the pairs it evaluated: those with at least
// MIN_JOINT genomes called R or S for both drugs (src/report.rs). For the
// pairs it leaves out, the joint counts are tallied here from the same
// profile strings report.rs tallies them from, so the table can show why
// each pair is not evaluated. A pair whose drug has no resistant or no
// susceptible genome among the jointly called ones has no phi: report.rs
// writes 0 there, and this file reads it as undefined, not as no association.

import type { SpeciesReport } from "../data/schema";

/** report.rs evaluates a pair only from this many jointly called genomes. */
export const MIN_JOINT = 30;

export type PairStatus = "evaluated" | "undefined" | "not-evaluated";

export interface PairRow {
  key: string;
  /** Drug columns, a before b in the file's drug order. */
  ia: number;
  ib: number;
  a: string;
  b: string;
  /** Genomes called R or S for both drugs, and the four joint cells. */
  n: number;
  rr: number;
  rs: number;
  sr: number;
  ss: number;
  /** Null unless the pair is evaluated with a defined phi. */
  phi: number | null;
  status: PairStatus;
}

export const pairKey = (a: string, b: string): string => `${a}|${b}`;

/** A margin of the 2x2 table is empty: one drug has no resistant or no
 *  susceptible genome among the jointly called ones. */
export const zeroMargin = (p: { rr: number; rs: number; sr: number; ss: number }): boolean =>
  p.rr + p.rs === 0 || p.sr + p.ss === 0 || p.rr + p.sr === 0 || p.rs + p.ss === 0;

/** Joint R and S calls for drug columns a and b, from the profiles. */
export const jointCalls = (report: SpeciesReport, a: number, b: number): { n: number; rr: number; rs: number; sr: number; ss: number } => {
  let rr = 0;
  let rs = 0;
  let sr = 0;
  let ss = 0;
  report.genomes.forEach((g) => {
    const pa = g.profile[a];
    const pb = g.profile[b];
    if (pa === "R" && pb === "R") rr += 1;
    else if (pa === "R" && pb === "S") rs += 1;
    else if (pa === "S" && pb === "R") sr += 1;
    else if (pa === "S" && pb === "S") ss += 1;
  });
  return { n: rr + rs + sr + ss, rr, rs, sr, ss };
};

/** All pairs in drug order; evaluated pairs carry the file's own counts. */
export const crossPairs = (report: SpeciesReport): PairRow[] => {
  const listed = new Map(report.cross.map((c) => [pairKey(c.a, c.b), c] as const));
  const rows: PairRow[] = [];
  report.drugs.forEach((a, ia) => {
    report.drugs.forEach((b, ib) => {
      if (ib <= ia) return;
      const key = pairKey(a, b);
      const c = listed.get(key);
      const flipped = c === undefined ? listed.get(pairKey(b, a)) : undefined;
      if (c !== undefined || flipped !== undefined) {
        // A pair written in the other order swaps its one-sided cells.
        const counts = c !== undefined ? c : { ...flipped!, rs: flipped!.sr, sr: flipped!.rs };
        const undef = zeroMargin(counts);
        rows.push({ key, ia, ib, a, b, n: counts.n, rr: counts.rr, rs: counts.rs, sr: counts.sr, ss: counts.ss, phi: undef ? null : counts.phi, status: undef ? "undefined" : "evaluated" });
        return;
      }
      rows.push({ key, ia, ib, a, b, ...jointCalls(report, ia, ib), phi: null, status: "not-evaluated" });
    });
  });
  return rows;
};

export interface PairSummary {
  total: number;
  /** At least MIN_JOINT jointly called genomes, phi defined or not. */
  evaluated: number;
  undefined: number;
  notEvaluated: number;
  /** The evaluated pair with the highest defined phi, first in drug order on a tie. */
  highest: PairRow | null;
}

export const pairSummary = (rows: readonly PairRow[]): PairSummary => {
  let highest: PairRow | null = null;
  for (const p of rows) {
    if (p.phi !== null && (highest === null || p.phi > (highest.phi ?? -Infinity))) highest = p;
  }
  return {
    total: rows.length,
    evaluated: rows.filter((p) => p.status !== "not-evaluated").length,
    undefined: rows.filter((p) => p.status === "undefined").length,
    notEvaluated: rows.filter((p) => p.status === "not-evaluated").length,
    highest,
  };
};

/** Why an evaluated pair has no phi, in words. */
export const undefinedWhy = (p: PairRow, name: (drug: string) => string = (d) => d): string => {
  const reasons: string[] = [];
  if (p.rr + p.rs === 0) reasons.push(`none is resistant to ${name(p.a)}`);
  else if (p.sr + p.ss === 0) reasons.push(`all are resistant to ${name(p.a)}`);
  if (p.rr + p.sr === 0) reasons.push(`none is resistant to ${name(p.b)}`);
  else if (p.rs + p.ss === 0) reasons.push(`all are resistant to ${name(p.b)}`);
  return reasons.join(" and ");
};
