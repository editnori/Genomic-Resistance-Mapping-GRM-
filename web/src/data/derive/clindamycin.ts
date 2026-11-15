// The clindamycin label check (docs/clindamycin-labels.md). In each of the
// six clindamycin datasets of one release, one deposited study, named by its
// PMID, holds most or a large share of the samples on which GRM's known-gene
// rule and the reported phenotype disagree, and the page sorts the six into
// four kinds. scripts/clindamycin-labels/site_data.py writes
// data/clindamycin-labels.json from reference/clindamycin-labels and the
// release's report checksums: each dataset's kind, the counts behind it and
// the released report's sha256.
//
// A note is given only for a report the check read: one of the six, at the
// file's evaluation revision, with the file's sample count. The app parses
// each report as JSON and keeps no bytes of it, so it cannot compare the
// loaded report's sha256 with the one the file carries. The intervals file
// is pinned more tightly, by revision and each row's confusion counts. Any
// other report, and any report at all when the bundled file does not
// validate, gets nothing.
//
// A disagreement is a difference between two readings, and neither is
// called the error. A change when a study leaves is an association with that
// study's labels or genomes; the words below name no cause and no gene.

import bundled from "../clindamycin-labels.json";
import type { CatalogRow, DatasetReport, Manifest } from "../schema";
import { ShapeError } from "../validate";
import { int, listWords } from "../format";

export type ClindamycinKind = "measured_dependence" | "one_study" | "class_balance" | "no_dependence";

/** The kinds in the page's order. */
export const CLINDAMYCIN_KINDS: readonly ClindamycinKind[] = ["measured_dependence", "one_study", "class_balance", "no_dependence"];

/** The page's name of each kind. */
export const KIND_WORDS: Record<ClindamycinKind, string> = {
  measured_dependence: "Measured dependence",
  one_study: "One study",
  class_balance: "Class balance, dependence not measured",
  no_dependence: "No dependence under the rule",
};

export interface ClindamycinDataset {
  kind: ClindamycinKind;
  /** The study the plan picks. */
  pmid: string;
  samples: number;
  study_samples: number;
  /** Samples where the known-gene rule and the reported phenotype disagree. */
  disagreements: number;
  study_disagreements: number;
  no_pmid_samples: number;
  no_pmid_disagreements: number;
  /** Samples, and their disagreements, under a PMID other than the study's. */
  other_pmid_samples: number;
  other_pmid_disagreements: number;
  /** Reported resistant and susceptible samples outside the study. */
  left_resistant: number;
  left_susceptible: number;
  scored_without: boolean;
  /** Null when the dataset was not scored without its study. */
  rows_passing: number | null;
  rows_tested: number | null;
  distinct_passing: number | null;
  distinct_tested: number | null;
  /** Whether the passing changes are rises or falls; null when none passes. */
  direction: "rises" | "falls" | null;
  /** The released report's sha256, from the release's checksum lines. */
  report_sha256: string;
}

export interface ClindamycinLabels {
  page: string;
  reference: string;
  release: string;
  evaluation_revision: number;
  min_class: number;
  /** The plan's threshold for a change, an adapted research screen. */
  threshold: number;
  sources: Record<string, string>;
  datasets: Record<string, ClindamycinDataset>;
}

const TOP = ["page", "reference", "release", "evaluation_revision", "min_class", "threshold", "sources", "datasets"] as const;
const FIELDS = [
  "kind",
  "pmid",
  "samples",
  "study_samples",
  "disagreements",
  "study_disagreements",
  "no_pmid_samples",
  "no_pmid_disagreements",
  "other_pmid_samples",
  "other_pmid_disagreements",
  "left_resistant",
  "left_susceptible",
  "scored_without",
  "rows_passing",
  "rows_tested",
  "distinct_passing",
  "distinct_tested",
  "direction",
  "report_sha256",
] as const;
const TESTS = ["rows_passing", "rows_tested", "distinct_passing", "distinct_tested"] as const;

type Obj = Record<string, unknown>;
const isObj = (v: unknown): v is Obj => typeof v === "object" && v !== null && !Array.isArray(v);
const isCount = (v: unknown): v is number => typeof v === "number" && Number.isInteger(v) && v >= 0;
const isSha = (v: unknown): v is string => typeof v === "string" && /^[0-9a-f]{64}$/.test(v);
const sameKeys = (o: Obj, keys: readonly string[]): boolean => {
  const own = Object.keys(o);
  return own.length === keys.length && own.every((k, i) => k === keys[i]);
};

function need(cond: boolean, message: string): asserts cond {
  if (!cond) throw new ShapeError(`clindamycin-labels: ${message}`);
}

/** The page's kind of a dataset, by the writer's rule: a dataset not scored
 *  without its study is "one study" when no sample is left without it and
 *  "class balance" otherwise; a scored one is "measured dependence" when a
 *  distinct test passes and "no dependence" when none does. */
export const kindOf = (d: Pick<ClindamycinDataset, "scored_without" | "left_resistant" | "left_susceptible" | "distinct_passing">): ClindamycinKind =>
  !d.scored_without
    ? d.left_resistant === 0 && d.left_susceptible === 0
      ? "one_study"
      : "class_balance"
    : (d.distinct_passing ?? 0) > 0
      ? "measured_dependence"
      : "no_dependence";

const checkDataset = (id: string, d: unknown, file: { min_class: number }): void => {
  const where = `datasets.${id}`;
  need(/^[a-z]+_[a-z_]+__clindamycin$/.test(id), `${id} is not a clindamycin dataset id`);
  need(isObj(d) && sameKeys(d, FIELDS), `${where} should hold exactly ${FIELDS.join(", ")}, in that order, and no field per sample`);
  need(typeof d.kind === "string" && (CLINDAMYCIN_KINDS as readonly string[]).includes(d.kind), `${where}.kind is not one of the four kinds`);
  need(typeof d.pmid === "string" && /^\d+$/.test(d.pmid), `${where}.pmid should be a PMID`);
  for (const k of ["samples", "study_samples", "disagreements", "study_disagreements", "no_pmid_samples", "no_pmid_disagreements", "other_pmid_samples", "other_pmid_disagreements", "left_resistant", "left_susceptible"] as const)
    need(isCount(d[k]), `${where}.${k} should be a count`);
  const [samples, study, dis, studyDis, noPmid, noPmidDis, other, otherDis, leftR, leftS] = [d.samples, d.study_samples, d.disagreements, d.study_disagreements, d.no_pmid_samples, d.no_pmid_disagreements, d.other_pmid_samples, d.other_pmid_disagreements, d.left_resistant, d.left_susceptible] as number[];
  need(samples! > 0 && study! + leftR! + leftS! === samples, `${where}: the study's samples and those left should add up to the dataset's`);
  need(study! + noPmid! + other! === samples && studyDis! + noPmidDis! + otherDis! === dis, `${where}: the study, the samples with no PMID and the other PMIDs should add up to the dataset's samples and disagreements`);
  need(studyDis! <= study! && noPmidDis! <= noPmid! && otherDis! <= other!, `${where}: a part is larger than its whole`);
  need(typeof d.scored_without === "boolean" && d.scored_without === Math.min(leftR!, leftS!) >= file.min_class, `${where}.scored_without should say whether both classes keep the class minimum`);
  if (d.scored_without) {
    for (const k of TESTS) need(isCount(d[k]), `${where}.${k} should be a count`);
    const [rowsPassing, rowsTested, passing, tested] = TESTS.map((k) => d[k] as number);
    need(passing! <= tested! && tested! <= rowsTested! && rowsPassing! <= rowsTested! && passing! <= rowsPassing!, `${where}: more tests pass than were tested`);
    need(passing === 0 ? d.direction === null : d.direction === "rises" || d.direction === "falls", `${where}.direction should be rises or falls when a test passes, and null when none does`);
  } else {
    need(TESTS.every((k) => d[k] === null) && d.direction === null, `${where}: a dataset not scored without its study has no tests`);
  }
  need(isSha(d.report_sha256), `${where}.report_sha256 should be a sha256`);
  need(d.kind === kindOf(d as unknown as ClindamycinDataset), `${where}.kind is ${String(d.kind)}, and its counts give ${kindOf(d as unknown as ClindamycinDataset)}`);
};

/** The bundled file, checked: its fields and only those, every count a
 *  whole number whose parts fit, and each kind the one its counts give. */
export function parseClindamycinLabels(raw: unknown): ClindamycinLabels {
  need(isObj(raw) && sameKeys(raw, TOP), `expected exactly ${TOP.join(", ")}`);
  need(typeof raw.page === "string" && raw.page.endsWith(".md") && typeof raw.reference === "string" && typeof raw.release === "string", "page, reference and release should be names");
  need(isCount(raw.evaluation_revision) && isCount(raw.min_class) && raw.min_class > 0, "evaluation_revision and min_class should be counts");
  need(typeof raw.threshold === "number" && raw.threshold > 0 && raw.threshold < 1, "threshold should lie between 0 and 1");
  need(isObj(raw.sources) && Object.keys(raw.sources).length > 0 && Object.values(raw.sources).every(isSha), "sources should map each file read to its sha256");
  need(isObj(raw.datasets) && Object.keys(raw.datasets).length > 0, "datasets should be an object");
  for (const [id, d] of Object.entries(raw.datasets)) checkDataset(id, d, raw as { min_class: number });
  return raw as unknown as ClindamycinLabels;
}

/** The bundled file, or null when it does not validate, so that no note shows. */
export const CLINDAMYCIN_LABELS: ClindamycinLabels | null = (() => {
  try {
    return parseClindamycinLabels(bundled);
  } catch {
    return null;
  }
})();

/** The file's entry for a dataset the check read: the same dataset id,
 *  evaluation revision and sample count. Null for any other. */
const entryFor = (file: ClindamycinLabels | null, id: string, revision: number | undefined, samples: number): ClindamycinDataset | null => {
  if (file === null || !Object.hasOwn(file.datasets, id)) return null;
  const entry = file.datasets[id]!;
  return revision === file.evaluation_revision && samples === entry.samples ? entry : null;
};

/** The entry for a loaded dataset report, or null. */
export const clindamycinEntry = (r: DatasetReport, file: ClindamycinLabels | null = CLINDAMYCIN_LABELS): ClindamycinDataset | null =>
  entryFor(file, r.id, r.evaluation?.revision, r.counts.samples);

/** The entry for a row of the datasets table, pinned the same way, or null. */
export const clindamycinRowEntry = (row: CatalogRow, file: ClindamycinLabels | null = CLINDAMYCIN_LABELS): ClindamycinDataset | null =>
  entryFor(file, row.id, row.evaluation?.revision, row.n);

export interface ClindamycinNote {
  kind: ClindamycinKind;
  /** "Clindamycin labels by study: measured dependence". */
  title: string;
  /** The counts behind the kind, as sentences. */
  counts: string[];
  /** What the released figure can be read as, by the page. */
  reading: string;
  /** The page, by its path in the repository. */
  page: string;
}

const samplesWord = (n: number, what: string): string => `${int(n)} ${what} ${n === 1 ? "sample" : "samples"}`;

/** The study's share of the disagreements and of the samples. */
const shareSentence = (d: ClindamycinDataset): string => {
  const between = "between GRM's known-gene rule and the reported phenotype";
  if (d.study_samples === d.samples) return `Every sample comes from PMID ${d.pmid} (${int(d.study_samples)} of ${int(d.samples)}), and the known-gene rule and the reported phenotype disagree on ${int(d.disagreements)} of them.`;
  const held = d.disagreements === 1 && d.study_disagreements === 1 ? `the dataset's single disagreement ${between}` : `${int(d.study_disagreements)} of the dataset's ${int(d.disagreements)} disagreements ${between}`;
  return `PMID ${d.pmid} holds ${held}, on ${int(d.study_samples)} of its ${int(d.samples)} samples.`;
};

/** Whether the study's share of the disagreements stands out only against
 *  the samples with no PMID: most samples carry none, the study disagrees
 *  more often than they do, and among the samples with a PMID the study
 *  holds at most one disagreement more than its share of them gives (the
 *  page's "close or reversed"). The rule of site_data.py's docstring, which
 *  tests/clindamycin_labels holds to the two Campylobacter datasets. */
export const noPmidContrast = (d: ClindamycinDataset): boolean =>
  2 * d.no_pmid_samples > d.samples &&
  d.study_disagreements * d.no_pmid_samples > d.no_pmid_disagreements * d.study_samples &&
  (d.study_disagreements - 1) * (d.study_samples + d.other_pmid_samples) <= (d.study_disagreements + d.other_pmid_disagreements) * d.study_samples;

const noPmidSentence = (d: ClindamycinDataset): string[] =>
  noPmidContrast(d) ? [`${int(d.no_pmid_samples)} of the ${int(d.samples)} samples carry no PMID, and the study's share of the disagreements stands out only against those samples.`] : [];

/** The count behind the kind. */
const kindSentence = (d: ClindamycinDataset, file: ClindamycinLabels): string[] => {
  const screen = `(${file.threshold} is an adapted research screen chosen for this check)`;
  switch (d.kind) {
    case "measured_dependence": {
      const passing = d.distinct_passing!;
      const move = d.direction === "rises" ? (passing === 1 ? "rises" : "rise") : passing === 1 ? "falls" : "fall";
      return [`Without that study, ${int(passing)} of the ${int(d.distinct_tested!)} distinct figures tested ${move} by ${file.threshold} or more with lineage-cluster intervals that exclude 0 ${screen}.`];
    }
    case "no_dependence":
      return [`Without that study, none of the ${int(d.distinct_tested!)} distinct figures tested changes by ${file.threshold} or more with a lineage-cluster interval that excludes 0 ${screen}.`];
    case "class_balance": {
      const under = [
        ...(d.left_resistant < file.min_class ? [samplesWord(d.left_resistant, "reported resistant")] : []),
        ...(d.left_susceptible < file.min_class ? [samplesWord(d.left_susceptible, "reported susceptible")] : []),
      ];
      const one = under.length === 1 && (d.left_resistant < file.min_class ? d.left_resistant : d.left_susceptible) === 1;
      return [`Without it ${listWords(under)} ${one ? "remains" : "remain"}, under the class minimum of ${int(file.min_class)} this release uses, so the dataset cannot be scored without that study.`];
    }
    case "one_study":
      return [];
  }
};

/** For measured dependence: the figures that pass, by their count, are
 *  associated with the study, and the others do not pass. */
const passingReading = (d: ClindamycinDataset): string => {
  const passing = d.distinct_passing!;
  const others = d.distinct_tested! - passing;
  const those = passing === 1 ? "That figure is" : `Those ${int(passing)} figures are`;
  if (others === 0) return `${those} associated with that one study's labels or genomes, and the result does not show which reading is wrong.`;
  return `${those} associated with that one study's labels or genomes. The other ${int(others)} ${others === 1 ? "does" : "do"} not pass the rule, and neither result shows which reading is wrong.`;
};

/** What the released figure can be read as, by the page's reading for each
 *  kind; tests/clindamycin_labels and the report list the page's lines. */
const READING: Record<ClindamycinKind, (d: ClindamycinDataset) => string> = {
  measured_dependence: passingReading,
  one_study: () => "The figures describe that study's genomes and labels alone, and no figure without that study can exist.",
  class_balance: () => "The dataset depends on that study for its class balance, and whether its figures depend on the study was not measured.",
  no_dependence: () => "This does not show that the figures do not depend on that study: some intervals are wide, and the fold plans of this dataset were not varied.",
};

/** The note for a loaded dataset report, or null when the check did not
 *  read this report. */
export const clindamycinNote = (r: DatasetReport, file: ClindamycinLabels | null = CLINDAMYCIN_LABELS): ClindamycinNote | null => {
  const d = clindamycinEntry(r, file);
  if (d === null || file === null) return null;
  return {
    kind: d.kind,
    title: `Clindamycin labels by study: ${KIND_WORDS[d.kind].toLowerCase()}`,
    counts: [shareSentence(d), ...noPmidSentence(d), ...kindSentence(d, file)],
    reading: READING[d.kind](d),
    page: file.page,
  };
};

/** The mark's tooltip for a row of the datasets table, or null. */
export const clindamycinMarkTip = (row: CatalogRow, file: ClindamycinLabels | null = CLINDAMYCIN_LABELS): string | null => {
  const d = clindamycinRowEntry(row, file);
  return d === null ? null : `Clindamycin labels by study: ${KIND_WORDS[d.kind].toLowerCase()}. The dataset report gives the counts behind it.`;
};

/** How many datasets the file gives each kind, in the page's order. */
export const kindCounts = (file: ClindamycinLabels): { kind: ClindamycinKind; datasets: number }[] =>
  CLINDAMYCIN_KINDS.map((kind) => ({ kind, datasets: Object.values(file.datasets).filter((d) => d.kind === kind).length }));

/** The Methods paragraph: what the check did, the four kinds with their
 *  counts of datasets, and the page. Null unless the manifest is the
 *  release the check read. */
export const clindamycinMethods = (manifest: Pick<Manifest, "evaluation_revision">, file: ClindamycinLabels | null = CLINDAMYCIN_LABELS): string | null => {
  if (file === null || manifest.evaluation_revision !== file.evaluation_revision) return null;
  const n = Object.keys(file.datasets).length;
  // A kind's name can hold a comma, so the list is joined by semicolons.
  const kinds = kindCounts(file).map(({ kind, datasets }, i) => `${KIND_WORDS[kind].toLowerCase()} (${i === 0 ? `${int(datasets)} ${datasets === 1 ? "dataset" : "datasets"}` : int(datasets)})`);
  const listed = kinds.length <= 1 ? (kinds[0] ?? "") : `${kinds.slice(0, -1).join("; ")}; and ${kinds[kinds.length - 1]}`;
  return `A check of release ${file.release} (evaluation revision ${int(file.evaluation_revision)}) asked, for each of its ${int(n)} clindamycin datasets, whether one deposited study holds a large share of the samples on which GRM's known-gene rule and the reported phenotype disagree, and whether the released figures depend on that study. A disagreement is a difference between two readings, and neither is called the error. The page sorts the ${int(n)} datasets into four kinds: ${listed}. Each of their dataset reports gives its kind with the counts behind it.`;
};
