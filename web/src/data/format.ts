import type { Learner } from "./schema";

const grouped = new Intl.NumberFormat("en-US");

/** 1234 -> "1,234". */
export const int = (n: number): string => grouped.format(Math.round(n));

/** A share in [0, 1] as a percent: 0.9412 -> "94.1%". */
export const pct = (share: number | null | undefined, digits = 1): string =>
  share === null || share === undefined || !Number.isFinite(share) ? "–" : `${(share * 100).toFixed(digits)}%`;

/** A difference of two shares in percentage points: 0.031 -> "3.1 points". */
export const points = (diff: number, digits = 1): string => {
  const v = Math.abs(diff * 100).toFixed(digits);
  return `${v} ${v === "1.0" ? "point" : "points"}`;
};

/** Plural noun for a count: plural(1, "dataset") -> "1 dataset". */
export const plural = (n: number, one: string, many = `${one}s`): string => `${int(n)} ${n === 1 ? one : many}`;

/** ["a", "b", "c"] -> "a, b and c". */
export const listWords = (items: readonly string[]): string =>
  items.length <= 1 ? (items[0] ?? "") : `${items.slice(0, -1).join(", ")} and ${items[items.length - 1]}`;

/** "Klebsiella pneumoniae" -> "K. pneumoniae", joined by a non-breaking
 *  space so the abbreviation never ends a line. */
export const speciesShort = (species: string): string => {
  const [genus, ...rest] = species.split(" ");
  return genus !== undefined && rest.length > 0 ? `${genus.charAt(0)}. ${rest.join(" ")}` : species;
};

/** Sentence-case a lowercase drug name for a title or table cell. */
export const capitalize = (s: string): string => (s === "" ? s : s.charAt(0).toUpperCase() + s.slice(1));

export const LEARNER_LABEL: Record<Learner, string> = {
  scm: "SCM",
  cart: "CART",
  rf: "Random forest",
  gbt: "Gradient boosting",
  lasso: "Lasso",
  catalog: "Catalog",
};

export const LEARNER_NAME: Record<Learner, string> = {
  scm: "Set covering machine",
  cart: "Decision tree (CART)",
  rf: "Random forest",
  gbt: "Gradient boosted trees",
  lasso: "L1 logistic regression",
  catalog: "Known-gene catalog",
};

/** A name inside a sentence: "Random forest" becomes "random forest", while
 *  "L1 logistic regression" and "CTX-M PCR" keep their capitals. */
export const inSentence = (name: string): string =>
  /^[A-Z][a-z]/.test(name) ? name.charAt(0).toLowerCase() + name.slice(1) : name;

/** The learner's name inside a sentence. */
export const learnerInSentence = (l: Learner): string => inSentence(LEARNER_NAME[l]);

/** Learners in the order the pipeline runs them; catalog last. */
export const LEARNER_ORDER: readonly Learner[] = ["scm", "cart", "rf", "gbt", "lasso", "catalog"];

/** Lab dilution labels for sub-1 mg/L steps (2^-4 is reported as 0.06). */
const SUB_ONE: Record<number, string> = { [-1]: "0.5", [-2]: "0.25", [-3]: "0.125", [-4]: "0.06", [-5]: "0.03", [-6]: "0.015", [-7]: "0.008", [-8]: "0.004" };

/** log2 dilution step to a concentration label: 3 -> "8", -2 -> "0.25". */
export const mgL = (log2: number): string => {
  if (log2 >= 0) return int(2 ** log2);
  return SUB_ONE[log2] ?? String(2 ** log2);
};

/** ISO timestamp to "Oct 23, 2026". */
export const dateText = (iso: string): string => {
  const d = new Date(iso);
  return Number.isNaN(d.getTime())
    ? iso
    : d.toLocaleDateString("en-US", { month: "short", day: "numeric", year: "numeric", timeZone: "UTC" });
};

/** Median of finite values; null for an empty list. */
export const median = (values: readonly number[]): number | null => {
  const v = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (v.length === 0) return null;
  const mid = v.length >> 1;
  return v.length % 2 === 1 ? v[mid]! : (v[mid - 1]! + v[mid]!) / 2;
};

/** A Mash distance stored as f32 (0.004999999888) to three significant
 *  figures: "0.005". */
export const distance = (v: number): string => String(Number(v.toPrecision(3)));

/** BV-BRC MLST strings ("MLST.klebsiella.512") as sequence types ("ST512");
 *  anything else passes through. */
export const stLabel = (mlst: string): string => {
  const m = /^MLST\.[^.]+\.(\d+)$/.exec(mlst);
  return m === null ? mlst : `ST${m[1]}`;
};

/** A per-dataset difference (random split minus lineages held out) in words,
 *  worded by its sign so a higher lineage score never reads as a loss. */
export const gapWords = (gap: number): string =>
  gap >= 0 ? `${points(gap)} lower with lineages held out` : `${points(-gap)} higher with lineages held out`;
