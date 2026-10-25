// Revision-5 research for the dataset report on screen. research.json pins
// each entry to the revision-5 cohort its lane read, and a dataset id is not
// that cohort: a later revision can keep the id and change the genomes, their
// order, labels, MICs, or held-out scores. An entry therefore applies only
// when every digest in its guard equals the same digest of the loaded report,
// and an MIC entry only to a report that has an MIC model.
// src/research.rs writes the guards with this encoding, one line per value,
// each ending in "\n", hashed with SHA-256:
//
//   sample_ids     sample id, in report order
//   labels         R or S, then m when the label came from an MIC
//   lineages       lineage cluster, decimal
//   mic_values     - without an MIC, else "<censor> <log2>"
//   mic_predicted  saved predicted log2 MICs: - for null, else the number's
//                  IEEE-754 bits as 16 lowercase hex digits; null without an
//                  MIC model
//   scores         for each split and learner in report order,
//                  "#<scheme> <learner>", then its held-out scores encoded as
//                  mic_predicted is

import type { Resource } from "./load";
import type { CohortGuard, DatasetReport, MicDataset, MicGuard, Research, SelectiveDataset, SelectiveGuard } from "./schema";

const K = Uint32Array.from([
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
  0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
  0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
  0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]);

const rotr = (x: number, n: number): number => (x >>> n) | (x << (32 - n));

/** SHA-256 of a string's UTF-8 bytes, as lowercase hex (FIPS 180-4). It is
 *  synchronous so a render can call it; crypto.subtle is asynchronous. */
export const sha256Hex = (text: string): string => {
  const bytes = new TextEncoder().encode(text);
  const data = new Uint8Array(Math.ceil((bytes.length + 9) / 64) * 64);
  data.set(bytes);
  data[bytes.length] = 0x80;
  const view = new DataView(data.buffer);
  const bits = bytes.length * 8;
  view.setUint32(data.length - 8, Math.floor(bits / 2 ** 32));
  view.setUint32(data.length - 4, bits >>> 0);
  const h = Uint32Array.from([0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]);
  const w = new Uint32Array(64);
  for (let block = 0; block < data.length; block += 64) {
    for (let t = 0; t < 16; t += 1) w[t] = view.getUint32(block + t * 4);
    for (let t = 16; t < 64; t += 1) {
      const a = w[t - 15]!;
      const b = w[t - 2]!;
      w[t] = w[t - 16]! + (rotr(a, 7) ^ rotr(a, 18) ^ (a >>> 3)) + w[t - 7]! + (rotr(b, 17) ^ rotr(b, 19) ^ (b >>> 10));
    }
    let a = h[0]!;
    let b = h[1]!;
    let c = h[2]!;
    let d = h[3]!;
    let e = h[4]!;
    let f = h[5]!;
    let g = h[6]!;
    let k = h[7]!;
    for (let t = 0; t < 64; t += 1) {
      const t1 = (k + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + K[t]! + w[t]!) | 0;
      const t2 = ((rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c))) | 0;
      k = g;
      g = f;
      f = e;
      e = (d + t1) | 0;
      d = c;
      c = b;
      b = a;
      a = (t1 + t2) | 0;
    }
    // Uint32Array stores each sum modulo 2^32.
    h[0]! += a;
    h[1]! += b;
    h[2]! += c;
    h[3]! += d;
    h[4]! += e;
    h[5]! += f;
    h[6]! += g;
    h[7]! += k;
  }
  return Array.from(h, (v) => v.toString(16).padStart(8, "0")).join("");
};

const scratch = new DataView(new ArrayBuffer(8));

const valueLine = (v: number | null): string => {
  if (v === null) return "-";
  scratch.setFloat64(0, v);
  return scratch.getUint32(0).toString(16).padStart(8, "0") + scratch.getUint32(4).toString(16).padStart(8, "0");
};

const digest = (lines: readonly string[]): string => sha256Hex(lines.map((line) => `${line}\n`).join(""));

/** Every guard digest of a report. */
export type ReportGuard = MicGuard & SelectiveGuard;
export type GuardPart = keyof ReportGuard;

const GUARD_PARTS: readonly GuardPart[] = ["samples", "sample_ids", "labels", "lineages", "mic_values", "mic_predicted", "scores"];

export const reportGuard = (report: DatasetReport): ReportGuard => {
  const s = report.samples;
  return {
    samples: s.length,
    // A newline inside an id would let two id lists hash alike, so such a
    // report matches no guard.
    sample_ids: s.some((x) => x.id.includes("\n")) ? "" : digest(s.map((x) => x.id)),
    labels: digest(s.map((x) => (x.resistant ? "R" : "S") + (x.from_mic ? "m" : ""))),
    lineages: digest(s.map((x) => String(x.lineage))),
    mic_values: digest(s.map((x) => (x.mic === null ? "-" : `${x.mic.censor} ${x.mic.log2}`))),
    mic_predicted: report.mic === null ? null : digest(report.mic.predicted.map(valueLine)),
    scores: digest(report.splits.flatMap((split) => split.learners.flatMap((l) => [`#${split.scheme} ${l.learner}`, ...l.scores.map(valueLine)]))),
  };
};

/** The parts of a guard that differ from the report's; empty when the entry
 *  describes this report's cohort. */
const guardFailures = (guard: CohortGuard, actual: ReportGuard): GuardPart[] =>
  GUARD_PARTS.filter((part) => part in guard && (guard as Partial<ReportGuard>)[part] !== actual[part]);

// A parsed report is never mutated, so its guard is hashed once per report
// object, however many sections and renders ask for its research.
const guards = new WeakMap<DatasetReport, ReportGuard>();

const guardOf = (report: DatasetReport): ReportGuard => {
  let guard = guards.get(report);
  if (guard === undefined) {
    guard = reportGuard(report);
    guards.set(report, guard);
  }
  return guard;
};

/** One lane's research for the report on screen. Only "current" carries the
 *  entry. "mismatch" means research exists for this dataset id but describes
 *  another cohort, so none of its values apply here; "loading" means
 *  research.json has not arrived, so a view shows no research yet. */
export type ResearchState<T> =
  | { status: "loading" }
  | { status: "absent" }
  | { status: "mismatch"; failed: GuardPart[] }
  | { status: "current"; entry: T };

/** The MIC lane's research. Its numbers describe GRM's MIC model, so they
 *  apply only to a report with an MIC model whose predictions the guard
 *  checked. "no_model" is a report without one, whether or not the rest of
 *  the guard matches; failed lists the parts that differ all the same. A
 *  mismatch's researchModel says whether the revision-5 report the research
 *  read had an MIC model, so a view can tell predictions the research never
 *  pinned from predictions that differ. "withheld" carries no numbers, only
 *  why, whatever the report's cohort. */
export type MicResearchState =
  | { status: "loading" }
  | { status: "absent" }
  | { status: "mismatch"; failed: GuardPart[]; researchModel: boolean }
  | { status: "no_model"; failed: GuardPart[] }
  | { status: "withheld"; reason: string }
  | { status: "current"; entry: MicDataset };

export interface DatasetResearch {
  mic: MicResearchState;
  selective: ResearchState<SelectiveDataset>;
}

const LOADING = { status: "loading" } as const;
const ABSENT = { status: "absent" } as const;

/** The research that applies to `report`. A missing or invalid
 *  research.json gives "absent" in both lanes, so it never stops the report. */
export const researchFor = (research: Resource<Research>, report: DatasetReport): DatasetResearch => {
  if (research.status === "loading") return { mic: LOADING, selective: LOADING };
  if (research.status === "error") return { mic: ABSENT, selective: ABSENT };
  const { mic, selective } = research.data;
  const check = <T extends { id: string; guard: CohortGuard }>(entries: readonly T[]) => {
    const entry = entries.find((e) => e.id === report.id);
    return entry === undefined ? undefined : { entry, failed: guardFailures(entry.guard, guardOf(report)) };
  };
  const micState = (): MicResearchState => {
    // A withheld model has no numbers for a guard to protect.
    const withheld = mic.withheld.find((w) => w.id === report.id);
    if (withheld !== undefined) return { status: "withheld", reason: withheld.reason };
    const found = check(mic.datasets);
    if (found === undefined) return ABSENT;
    if (report.mic === null) return { status: "no_model", failed: found.failed };
    if (found.failed.length > 0) return { status: "mismatch", failed: found.failed, researchModel: found.entry.guard.mic_predicted !== null };
    return { status: "current", entry: found.entry };
  };
  const sel = check(selective.datasets);
  return {
    mic: micState(),
    selective: sel === undefined ? ABSENT : sel.failed.length === 0 ? { status: "current", entry: sel.entry } : { status: "mismatch", failed: sel.failed },
  };
};
