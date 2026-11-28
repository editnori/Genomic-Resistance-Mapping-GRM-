// Data types, formatting, template filling, and deterministic made-up DNA
// helpers shared by the five flow diagram components. Every count and visible
// string comes from flow.json; none is typed into component code.

import { useEffect, useState } from "react";
import rawFlow from "./flow.json";
import { FLOW_SLUGS, type FlowSlug } from "./slugs";

export { FLOW_SLUGS, type FlowSlug };

export type CountFmt = "dec1" | "dec3" | "dec4" | "pc0" | "pc1" | "pts";

export interface SourceRef {
  root: string;
  file: string;
  line: number;
  text?: string;
  key?: string[];
  col?: string;
  where?: Record<string, string>;
  distinct?: string;
  regex_count?: string;
  word?: string;
}

export interface FlowCount {
  v: number;
  unit?: string;
  fmt?: CountFmt;
  of?: string;
  src: SourceRef;
}

export type CardUse = "used" | "checked" | "kept" | "unread";
export type AltState = "no" | "pass" | "open";
export type BarCls = "on" | "off" | "hatch";
export type NodeMark = "kept" | "cand" | "run" | "none" | "stated" | "wait" | "meas";
export type GlyphId = "records" | "datasets" | "genome" | "pieces" | "split" | "filter" | "models" | "served" | "score";

export interface FlowNode {
  id: string;
  kind: string;
  state: string;
  title: string;
  line?: string;
  note?: string;
  tag?: string;
  mark?: NodeMark;
  bare?: boolean;
  count?: string;
  glyph?: GlyphId;
  draw?: "contigs";
  from?: string;
  step?: string;
  test?: string;
  decisions?: string[];
  fields?: [string, string][];
  alts?: [string, AltState][];
  cardtitle?: string;
  card?: [string, CardUse, string][];
  bars?: [string, string, BarCls][];
  stats?: [string, string, string][];
  tools?: string[];
  src: SourceRef[];
}

export interface FlowEdge {
  from: string;
  to: string;
  label?: string;
  live?: boolean;
  draw?: boolean;
}

export interface BaseViewSpec {
  file: string;
  kicker: string;
  title: string;
  sub: string;
  start: string;
  ui: Record<string, string>;
}

export interface WalkViewSpec extends BaseViewSpec {
  steps: string[];
}

export interface GenomePathViewSpec extends BaseViewSpec {
  rows: [string | null, string | null][];
  next: Record<string, number>;
  walk: string[];
}

export interface DecisionsViewSpec extends BaseViewSpec {
  trunk: [string, string[]][];
}

export interface FlowData {
  about: string;
  pinned: { repo: string; commit: string };
  roots: Record<string, string>;
  states: string[];
  kinds: string[];
  literals: { text: string; src: SourceRef }[];
  counts: Record<string, FlowCount>;
  nodes: FlowNode[];
  edges: FlowEdge[];
  views: {
    genome_table: WalkViewSpec;
    genome_path: GenomePathViewSpec;
    bvbrc: WalkViewSpec;
    piece_count: WalkViewSpec;
    decisions: DecisionsViewSpec;
  };
  sums?: { total: string; parts: string[] }[];
}

export type FlowViewKey = keyof FlowData["views"];

export const FLOW_DATA: FlowData = rawFlow as unknown as FlowData;

/** Formats one count entry from flow.json identically to src/core.js. */
export const formatCount = (flow: FlowData, key: string): string => {
  const c = flow.counts[key];
  if (c === undefined) return `{${key}}`;
  const v = c.v;
  if (c.fmt === "pc1") return `${v.toFixed(1)}%`;
  if (c.fmt === "pts") return `${Math.round(v * 100)} points`;
  if (c.fmt === "dec3") return v.toFixed(3);
  if (c.fmt === "dec4") return v.toFixed(4);
  if (c.fmt === "dec1") return v.toFixed(1);
  if (c.fmt === "pc0") return `${v}%`;
  return v.toLocaleString("en-US");
};

/** Replaces every `{count}` placeholder in a flow.json string with its formatted value. */
export const fillTemplate = (flow: FlowData, template: string | undefined): string =>
  template === undefined ? "" : template.replace(/\{([a-z0-9_]+)\}/g, (_match, key: string) => formatCount(flow, key));

export const buildNodeMap = (flow: FlowData): Record<string, FlowNode> => {
  const map: Record<string, FlowNode> = {};
  for (const n of flow.nodes) map[n.id] = n;
  return map;
};

export const uiText = (ui: Record<string, string>, key: string): string => ui[key] ?? "";

export type BaseLetter = "A" | "C" | "G" | "T";
const BASES: readonly BaseLetter[] = ["A", "C", "G", "T"];

/** Deterministic PRNG (mulberry32) used only for made-up DNA illustrations. */
export const createSeededRng = (initialSeed: number): (() => number) => {
  let a = initialSeed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
};

/** Made-up DNA letters from a fixed seed, never a real genome sequence. */
export const makeLetters = (count: number, seed: number): BaseLetter[] => {
  const next = createSeededRng(seed);
  return Array.from({ length: count }, () => BASES[Math.floor(next() * 4)]!);
};

export const useReducedMotion = (): boolean => {
  const [reduced, setReduced] = useState<boolean>(() => {
    if (typeof window === "undefined" || typeof window.matchMedia !== "function") return false;
    return window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  });
  useEffect(() => {
    if (typeof window === "undefined" || typeof window.matchMedia !== "function") return undefined;
    const mq = window.matchMedia("(prefers-reduced-motion: reduce)");
    const onChange = () => setReduced(mq.matches);
    onChange();
    mq.addEventListener("change", onChange);
    return () => mq.removeEventListener("change", onChange);
  }, []);
  return reduced;
};
