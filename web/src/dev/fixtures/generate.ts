// A deterministic synthetic report set for the tests. Isolates,
// lineages, and resistance mechanisms are generated from a seed; every metric
// is then computed from those isolates with ports of the pipeline's own
// metric code, so the files agree with themselves the way real reports do.
// Nothing here is a measurement.

import {
  SCHEMA,
  type BinaryMetrics,
  type CatalogRow,
  type Corpus,
  type DatasetReport,
  type FeatureCard,
  type GenomeRow,
  type GroupRow,
  type Learner,
  type LearnerResult,
  type Lite,
  type Manifest,
  type Merge,
  type Mic,
  type MicResult,
  type PanelResult,
  type PanelSummaryRow,
  type Sample,
  type SpeciesReport,
  type SpeciesRow,
  type SplitResult,
  type TreeNode,
} from "../../data/schema";
import { encodeBits } from "../../data/bits";
import { binaryMetrics, conformal, normal, rng, rocCurve } from "./stats";

interface GroupDef {
  family: string;
  class: string;
  subclass: string;
  product: string;
}

const GROUPS: Record<string, GroupDef> = {
  blaKPC: { family: "blaKPC", class: "BETA-LACTAM", subclass: "CARBAPENEM", product: "carbapenem-hydrolyzing class A beta-lactamase KPC" },
  "blaOXA (OXA-48 family)": { family: "blaOXA", class: "BETA-LACTAM", subclass: "CARBAPENEM", product: "OXA-48 family carbapenem-hydrolyzing class D beta-lactamase" },
  blaNDM: { family: "blaNDM", class: "BETA-LACTAM", subclass: "CARBAPENEM", product: "subclass B1 metallo-beta-lactamase NDM" },
  "blaCTX-M": { family: "blaCTX-M", class: "BETA-LACTAM", subclass: "CEPHALOSPORIN", product: "class A extended-spectrum beta-lactamase CTX-M" },
  blaSHV: { family: "blaSHV", class: "BETA-LACTAM", subclass: "BETA-LACTAM", product: "class A beta-lactamase SHV" },
  "aac(6')-Ib-cr": { family: "aac(6')-Ib-cr", class: "AMINOGLYCOSIDE/QUINOLONE", subclass: "AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN", product: "aminoglycoside 6'-N-acetyltransferase AAC(6')-Ib-cr" },
  qnrB: { family: "qnrB", class: "QUINOLONE", subclass: "QUINOLONE", product: "quinolone resistance pentapeptide repeat protein QnrB" },
  armA: { family: "armA", class: "AMINOGLYCOSIDE", subclass: "AMIKACIN/GENTAMICIN/KANAMYCIN/TOBRAMYCIN", product: "16S rRNA (guanine(1405)-N(7))-methyltransferase ArmA" },
  "mcr-1": { family: "mcr-1", class: "COLISTIN", subclass: "COLISTIN", product: "phosphoethanolamine--lipid A transferase MCR-1" },
  "tet(A)": { family: "tet(A)", class: "TETRACYCLINE", subclass: "TETRACYCLINE", product: "tetracycline efflux MFS transporter Tet(A)" },
  fosA: { family: "fosA", class: "FOSFOMYCIN", subclass: "FOSFOMYCIN", product: "fosfomycin resistance glutathione transferase FosA" },
  mecA: { family: "mecA", class: "BETA-LACTAM", subclass: "METHICILLIN", product: "PBP2a family beta-lactam-resistant peptidoglycan transpeptidase MecA" },
  "erm(C)": { family: "erm(C)", class: "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN", subclass: "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN", product: "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(C)" },
  blaZ: { family: "blaZ", class: "BETA-LACTAM", subclass: "PENICILLIN", product: "penicillin-hydrolyzing class A beta-lactamase BlaZ" },
};

interface Mechanism {
  id: string;
  /** Catalog group when the mechanism is an acquired gene. */
  group: string | null;
  gene: string;
  product: string;
  codons: [number, number] | null;
  nearest: [string, number] | null;
  /** Carriage probability for genomes of lineage l. */
  rate(l: number): number;
}

interface DrugSpec {
  antibiotic: string;
  drug_class: string;
  /** Resistant when any of these mechanisms is carried. */
  mechanisms: string[];
  rule: string[][] | null;
  panels: string[];
  mic: { s: number; r: number; lo: number; hi: number; r_above: number } | null;
  kover: { sensitivity: number; specificity: number; n_rules: number; n_examples: number } | null;
}

interface SpeciesSpec {
  species: string;
  slug: string;
  taxon: number;
  genomes: number;
  lineages: number;
  mechanisms: Mechanism[];
  drugs: DrugSpec[];
  seed: number;
}

const gene = (group: string, rate: (l: number) => number): Mechanism => ({
  id: group,
  group,
  gene: GROUPS[group]!.family,
  product: GROUPS[group]!.product,
  codons: [1, 290],
  nearest: null,
  rate,
});

const mutation = (id: string, g: string, product: string, codons: [number, number] | null, nearest: [string, number] | null, rate: (l: number) => number): Mechanism => ({
  id,
  group: null,
  gene: g,
  product,
  codons,
  nearest,
  rate,
});

const KPN: SpeciesSpec = {
  species: "Klebsiella pneumoniae",
  slug: "klebsiella_pneumoniae",
  taxon: 573,
  genomes: 168,
  lineages: 14,
  seed: 7,
  mechanisms: [
    gene("blaKPC", (l) => (l <= 1 ? 0.85 : l === 6 ? 0.3 : 0.02)),
    gene("blaOXA (OXA-48 family)", (l) => (l === 4 ? 0.8 : l === 9 ? 0.4 : 0.01)),
    gene("blaNDM", (l) => (l === 7 ? 0.7 : 0.01)),
    gene("blaCTX-M", (l) => (l <= 5 ? 0.75 : 0.15)),
    gene("blaSHV", () => 0.97),
    gene("aac(6')-Ib-cr", (l) => (l <= 3 ? 0.55 : 0.08)),
    gene("qnrB", (l) => (l === 2 || l === 8 ? 0.5 : 0.04)),
    gene("armA", (l) => (l === 7 || l === 4 ? 0.5 : 0.02)),
    gene("mcr-1", (l) => (l === 10 ? 0.5 : 0.0)),
    gene("tet(A)", (l) => (l % 3 === 0 ? 0.7 : 0.2)),
    gene("fosA", () => 0.98),
    mutation("gyrA-S83", "gyrA", "DNA gyrase subunit A", [83, 87], null, (l) => (l <= 3 || l === 6 ? 0.9 : 0.1)),
    mutation("mgrB", "mgrB", "PhoP/PhoQ regulator MgrB", null, ["mgrB", 46], (l) => (l <= 1 ? 0.35 : 0.03)),
    mutation("ompK36", "ompK36", "outer membrane porin OmpK36", [134, 136], null, (l) => (l <= 1 ? 0.5 : 0.05)),
  ],
  drugs: [
    { antibiotic: "meropenem", drug_class: "carbapenem", mechanisms: ["blaKPC", "blaOXA (OXA-48 family)", "blaNDM"], rule: [["S:CARBAPENEM"]], panels: ["carbapenemase", "gram_negative_bcid"], mic: { s: -4, r: 4, lo: -5, hi: 5, r_above: 1 }, kover: { sensitivity: 0.95, specificity: 0.97, n_rules: 2, n_examples: 1640 } },
    { antibiotic: "ceftriaxone", drug_class: "cephalosporin", mechanisms: ["blaCTX-M", "blaKPC", "blaNDM", "blaOXA (OXA-48 family)"], rule: [["S:CEPHALOSPORIN", "S:CARBAPENEM", "S:METHICILLIN"]], panels: ["gram_negative_bcid", "ctx_m"], mic: null, kover: null },
    { antibiotic: "ciprofloxacin", drug_class: "quinolone", mechanisms: ["gyrA-S83", "qnrB"], rule: [["C:QUINOLONE"]], panels: [], mic: { s: -3, r: 3, lo: -4, hi: 4, r_above: 0 }, kover: { sensitivity: 0.91, specificity: 0.93, n_rules: 3, n_examples: 1408 } },
    { antibiotic: "amikacin", drug_class: "aminoglycoside", mechanisms: ["armA"], rule: [["S:AMIKACIN", "S:AMINOGLYCOSIDE"]], panels: [], mic: { s: 1, r: 6, lo: 0, hi: 6, r_above: 4 }, kover: null },
    { antibiotic: "colistin", drug_class: "polymyxin", mechanisms: ["mgrB", "mcr-1"], rule: [["C:COLISTIN"]], panels: ["mcr"], mic: null, kover: null },
    { antibiotic: "tetracycline", drug_class: "tetracycline", mechanisms: ["tet(A)"], rule: [["C:TETRACYCLINE"]], panels: [], mic: null, kover: null },
  ],
};

const SAU: SpeciesSpec = {
  species: "Staphylococcus aureus",
  slug: "staphylococcus_aureus",
  taxon: 1280,
  genomes: 96,
  lineages: 8,
  seed: 11,
  mechanisms: [
    gene("mecA", (l) => (l <= 2 ? 0.9 : 0.1)),
    gene("erm(C)", (l) => (l === 1 || l === 5 ? 0.7 : 0.1)),
    gene("blaZ", () => 0.8),
    mutation("grlA-S80", "grlA", "DNA topoisomerase IV subunit A", [80, 84], null, (l) => (l <= 2 ? 0.85 : 0.15)),
  ],
  drugs: [
    { antibiotic: "oxacillin", drug_class: "penicillin", mechanisms: ["mecA"], rule: [["S:METHICILLIN"]], panels: ["mec"], mic: null, kover: { sensitivity: 0.99, specificity: 0.98, n_rules: 1, n_examples: 920 } },
    { antibiotic: "ciprofloxacin", drug_class: "quinolone", mechanisms: ["grlA-S80"], rule: [["C:QUINOLONE"]], panels: [], mic: null, kover: null },
    { antibiotic: "erythromycin", drug_class: "macrolide", mechanisms: ["erm(C)"], rule: [["C:MACROLIDE"]], panels: [], mic: null, kover: null },
  ],
};

interface PanelDef {
  id: string;
  name: string;
  mirrors: string;
  targets: string[];
}

const PANEL_DEFS: Record<string, PanelDef> = {
  carbapenemase: { id: "carbapenemase", name: "Carbapenemase PCR", mirrors: "Cepheid Xpert Carba-R targets: KPC, NDM, VIM, IMP, OXA-48-like", targets: ["blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaOXA (OXA-48 family)"] },
  gram_negative_bcid: { id: "gram_negative_bcid", name: "Gram-negative blood-culture panel", mirrors: "BioFire BCID2 resistance targets: CTX-M, KPC, NDM, VIM, IMP, OXA-48-like, mcr-1", targets: ["blaCTX-M", "blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaOXA (OXA-48 family)", "mcr-1"] },
  ctx_m: { id: "ctx_m", name: "CTX-M PCR", mirrors: "ESBL screening by CTX-M detection (BCID2 and similar)", targets: ["blaCTX-M"] },
  mec: { id: "mec", name: "mecA/mecC PCR", mirrors: "MRSA molecular assays (Xpert MRSA, BCID2 mecA/C)", targets: ["mecA", "mecC"] },
  mcr: { id: "mcr", name: "mcr PCR", mirrors: "Plasmid colistin resistance (mcr-1 to mcr-10) assays", targets: ["mcr-1", "mcr-2", "mcr-3", "mcr-4", "mcr-5"] },
};

const COUNTRIES = ["USA", "China", "Italy", "Greece", "India", "Brazil", "Spain", "Thailand"];
const HOSTS = ["Human", "Human", "Human", "Human", "Bovine", "Swine"];
const BASES = "acgt";

interface Genome {
  id: string;
  lineage: number;
  x: number;
  y: number;
  carries: Set<string>;
  calls: Map<string, "R" | "S" | "I" | ".">;
  country: string;
  year: string;
  host: string;
  mlst: string;
  length: number;
  contigs: number;
  n50: number;
}

const pick = <T>(rand: () => number, items: readonly T[]): T => items[Math.floor(rand() * items.length)]!;

const sigmoid = (z: number): number => 1 / (1 + Math.exp(-z));

const makeGenomes = (spec: SpeciesSpec, rand: () => number): Genome[] => {
  const weights = Array.from({ length: spec.lineages }, (_, l) => 1 / (l + 1) ** 0.8);
  const total = weights.reduce((a, b) => a + b, 0);
  const centers = weights.map(() => [normal(rand) * 0.06, normal(rand) * 0.04] as const);
  return Array.from({ length: spec.genomes }, (_, i): Genome => {
    let u = rand() * total;
    let lineage = 0;
    while (lineage < spec.lineages - 1 && u > weights[lineage]!) {
      u -= weights[lineage]!;
      lineage += 1;
    }
    const [cx, cy] = centers[lineage]!;
    const carries = new Set(spec.mechanisms.filter((m) => rand() < m.rate(lineage)).map((m) => m.id));
    const calls = new Map<string, "R" | "S" | "I" | ".">();
    spec.drugs.forEach((d) => {
      const u2 = rand();
      if (u2 < 0.12) calls.set(d.antibiotic, ".");
      else if (u2 < 0.14) calls.set(d.antibiotic, "I");
      else {
        const mech = d.mechanisms.some((m) => carries.has(m));
        calls.set(d.antibiotic, (rand() < 0.04 ? !mech : mech) ? "R" : "S");
      }
    });
    return {
      id: `${spec.taxon}.${10_000 + i * 37}`,
      lineage,
      x: cx + normal(rand) * 0.006,
      y: cy + normal(rand) * 0.006,
      carries,
      calls,
      country: pick(rand, COUNTRIES),
      year: String(2008 + Math.floor(rand() * 15)),
      host: pick(rand, HOSTS),
      mlst: `ST${[258, 11, 147, 307, 15, 101, 14, 231, 37, 392, 45, 20, 23, 17][lineage % 14]}`,
      length: Math.round(5_400_000 + normal(rand) * 120_000),
      contigs: Math.round(40 + rand() * 120),
      n50: Math.round(180_000 + rand() * 220_000),
    };
  });
};

/** Naive average-linkage clustering (O(n^3), fixture sizes only), with
 *  scipy-style ids: merge k creates node n + k. */
const upgma = (dist: (i: number, j: number) => number, n: number): Merge[] => {
  const d = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 0 : dist(i, j))));
  const active = Array.from({ length: n }, () => true);
  const node = Array.from({ length: n }, (_, i) => i);
  const size = Array.from({ length: n }, () => 1);
  const merges: Merge[] = [];
  for (let step = 0; step < n - 1; step += 1) {
    let bi = -1;
    let bj = -1;
    let best = Infinity;
    for (let i = 0; i < n; i += 1) {
      if (!active[i]) continue;
      for (let j = i + 1; j < n; j += 1) {
        if (active[j] && d[i]![j]! < best) {
          best = d[i]![j]!;
          bi = i;
          bj = j;
        }
      }
    }
    merges.push({ a: node[bi]!, b: node[bj]!, height: best, size: size[bi]! + size[bj]! });
    for (let k = 0; k < n; k += 1) {
      if (!active[k] || k === bi || k === bj) continue;
      const v = (d[bi]![k]! * size[bi]! + d[bj]![k]! * size[bj]!) / (size[bi]! + size[bj]!);
      d[bi]![k] = v;
      d[k]![bi] = v;
    }
    active[bj] = false;
    size[bi]! += size[bj]!;
    node[bi] = n + step;
  }
  return merges;
};

/** Port of cluster::leaf_order: roots by descending id, left child first. */
const leafOrder = (merges: readonly Merge[], n: number): number[] => {
  const used = new Array<boolean>(n + merges.length).fill(false);
  merges.forEach((m) => {
    used[m.a] = true;
    used[m.b] = true;
  });
  const out: number[] = [];
  for (let id = n + merges.length - 1; id >= 0; id -= 1) {
    if (used[id]) continue;
    const stack = [id];
    while (stack.length > 0) {
      const cur = stack.pop()!;
      if (cur < n) out.push(cur);
      else {
        const m = merges[cur - n]!;
        stack.push(m.b, m.a);
      }
    }
  }
  return out;
};

const tokenMatch = (token: string, def: GroupDef): boolean => {
  const field = token.startsWith("C:") ? def.class : def.subclass;
  return field.split("/").some((t) => t.trim() === token.slice(2));
};

const ruleCall = (rule: string[][], groups: string[]): boolean =>
  rule.every((alts) => groups.some((g) => alts.some((t) => tokenMatch(t, GROUPS[g]!))));

const context = (rand: () => number, spanLength: number): string => {
  const flank = () => Array.from({ length: 60 }, () => pick(rand, [...BASES])).join("");
  const span = Array.from({ length: spanLength }, () => pick(rand, [...BASES])).join("").toUpperCase();
  return `${flank()}${span}${flank()}`;
};

interface FeatureDraft {
  id: number;
  present: boolean[];
  card: Omit<FeatureCard, "bits" | "n_resistant" | "n_susceptible">;
}

/** Greedy weighted-Gini CART over the fixture features (src/learn/tree.rs). */
const fitTree = (features: FeatureDraft[], truth: boolean[], maxDepth: number): TreeNode[] => {
  const nPos = truth.filter(Boolean).length;
  const nNeg = truth.length - nPos;
  const wPos = truth.length / (2 * Math.max(nPos, 1));
  const wNeg = truth.length / (2 * Math.max(nNeg, 1));
  const gini = (p: number, q: number) => (p + q === 0 ? 0 : 2 * (p / (p + q)) * (1 - p / (p + q)));
  const nodes: TreeNode[] = [];
  const push = (mask: boolean[], depth: number): number => {
    const pos = mask.filter((m, i) => m && truth[i]).length;
    const neg = mask.filter((m, i) => m && !truth[i]).length;
    const wp = pos * wPos;
    const wn = neg * wNeg;
    nodes.push({ value: wp + wn === 0 ? 0.5 : wp / (wp + wn), n: pos + neg, n_resistant: pos, depth, split: null });
    return nodes.length - 1;
  };
  const stack: [number, boolean[]][] = [[push(truth.map(() => true), 0), truth.map(() => true)]];
  while (stack.length > 0) {
    const [id, mask] = stack.pop()!;
    const cur = nodes[id]!;
    if (cur.depth >= maxDepth || cur.n_resistant === 0 || cur.n_resistant === cur.n || cur.n < 2) continue;
    const np = cur.n_resistant * wPos;
    const nn = (cur.n - cur.n_resistant) * wNeg;
    const parent = (np + nn) * gini(np, nn);
    let best: { gain: number; f: FeatureDraft } | null = null;
    for (const f of features) {
      let pp = 0;
      let pn = 0;
      mask.forEach((m, i) => {
        if (m && f.present[i]) {
          if (truth[i]) pp += 1;
          else pn += 1;
        }
      });
      if (pp + pn === 0 || pp + pn === cur.n) continue;
      const lp = pp * wPos;
      const ln = pn * wNeg;
      const gain = parent - (lp + ln) * gini(lp, ln) - (np - lp + nn - ln) * gini(np - lp, nn - ln);
      if (gain > 1e-12 && (best === null || gain > best.gain)) best = { gain, f };
    }
    if (best === null) continue;
    const f = best.f;
    const absentMask = mask.map((m, i) => m && !f.present[i]);
    const presentMask = mask.map((m, i) => m && f.present[i] === true);
    const absent = push(absentMask, cur.depth + 1);
    const present = push(presentMask, cur.depth + 1);
    nodes[id] = { ...cur, split: { feature: f.id, gain: best.gain / truth.length, absent, present } };
    stack.push([absent, absentMask], [present, presentMask]);
  }
  return nodes;
};

const lite = (m: BinaryMetrics): Lite => ({
  balanced_accuracy: m.balanced_accuracy,
  sensitivity: m.sensitivity,
  specificity: m.specificity,
  vme: m.vme,
  me: m.me,
  auroc: m.auroc,
});

const LEARNERS: Learner[] = ["scm", "cart", "rf", "gbt", "lasso"];

/** The class minimum, the CLI default `grm curate` uses. */
const MIN_CLASS = 25;

interface Built {
  report: DatasetReport;
  row: CatalogRow;
  panels: PanelSummaryRow[];
}

const buildDataset = (spec: SpeciesSpec, drug: DrugSpec, genomes: Genome[], rand: () => number): Built | null => {
  const tested = genomes.filter((g) => {
    const c = g.calls.get(drug.antibiotic);
    return c === "R" || c === "S";
  });
  const truth = tested.map((g) => g.calls.get(drug.antibiotic) === "R");
  const n = tested.length;
  const nR = truth.filter(Boolean).length;
  if (nR < MIN_CLASS || n - nR < MIN_CLASS) return null;
  const id = `${spec.slug}__${drug.antibiotic.replace(/[^a-z0-9]+/g, "_")}`;

  // Catalog group dictionary in first-seen order, as the pipeline builds it.
  const dict = new Map<string, number>();
  const groupsOf = tested.map((g) => [...g.carries].filter((m) => GROUPS[m] !== undefined));
  groupsOf.forEach((gs) => gs.forEach((g) => dict.has(g) || dict.set(g, dict.size)));
  const catalogGroups = [...dict.keys()];

  // Features: each mechanism, two lineage markers, three noise features.
  const lineageR = new Map<number, number>();
  tested.forEach((g, i) => lineageR.set(g.lineage, (lineageR.get(g.lineage) ?? 0) + (truth[i] ? 1 : -1)));
  const enriched = new Set([...lineageR.entries()].filter(([, v]) => v > 2).map(([l]) => l));
  const usedIds = new Set<number>();
  const newId = () => {
    let f = Math.floor(rand() * 48_000);
    while (usedIds.has(f)) f += 1;
    usedIds.add(f);
    return f;
  };
  const drafts: FeatureDraft[] = [
    ...spec.mechanisms
      .filter((m) => drug.mechanisms.includes(m.id) || (m.group !== null && rand() < 0.25))
      .map((m): FeatureDraft => {
        const span = m.group !== null ? 31 + Math.floor(rand() * 40) : 31 + 4;
        return {
          id: newId(),
          present: tested.map((g) => g.carries.has(m.id)),
          card: {
            feature: 0,
            n_kmers: m.group !== null ? 600 + Math.floor(rand() * 400) : 31 + Math.floor(rand() * 6),
            kmer: context(rand, 31).slice(60, 91),
            catalog_families: m.group !== null ? [GROUPS[m.group]!.family] : [],
            catalog_products: m.group !== null ? [GROUPS[m.group]!.product] : [],
            locus: {
              genome_id: tested.find((g) => g.carries.has(m.id))?.id ?? tested[0]!.id,
              contig: `NZ_CP0${Math.floor(10_000 + rand() * 80_000)}.1`,
              start: 1_000 + Math.floor(rand() * 5_000_000),
              end: 0,
              kmers_placed: 31,
              gene: m.gene,
              product: m.product,
              patric_id: `fig|${spec.taxon}.${Math.floor(rand() * 90_000)}.peg.${Math.floor(rand() * 5_000)}`,
              gene_strand: rand() < 0.5 ? "+" : "-",
              codons: m.codons,
              nearest: m.nearest,
              context: context(rand, span),
            },
          },
        };
      }),
    ...[0, 1].map(
      (k): FeatureDraft => ({
        id: newId(),
        present: tested.map((g) => (enriched.has(g.lineage) ? rand() < 0.9 : rand() < 0.08 + k * 0.05)),
        card: {
          feature: 0,
          n_kmers: 120 + Math.floor(rand() * 900),
          kmer: context(rand, 31).slice(60, 91),
          catalog_families: [],
          catalog_products: [],
          locus: {
            genome_id: tested[0]!.id,
            contig: `NZ_CP0${Math.floor(10_000 + rand() * 80_000)}.1`,
            start: 1_000 + Math.floor(rand() * 5_000_000),
            end: 0,
            kmers_placed: 31,
            gene: k === 0 ? "wzi" : null,
            product: k === 0 ? "capsule assembly protein Wzi" : "hypothetical protein",
            patric_id: null,
            gene_strand: "+",
            codons: k === 0 ? [40, 52] : null,
            nearest: k === 0 ? null : ["hypothetical protein", 212],
            context: context(rand, 44),
          },
        },
      }),
    ),
    ...[0, 1, 2].map(
      (): FeatureDraft => ({
        id: newId(),
        present: tested.map(() => rand() < 0.3),
        card: { feature: 0, n_kmers: 31, kmer: context(rand, 31).slice(60, 91), catalog_families: [], catalog_products: [], locus: null },
      }),
    ),
  ];
  drafts.forEach((d) => {
    d.card.feature = d.id;
    if (d.card.locus !== null) d.card.locus.end = d.card.locus.start + d.card.locus.context.length - 121;
  });
  const features: FeatureCard[] = drafts.map((d) => ({
    ...d.card,
    n_resistant: d.present.filter((p, i) => p && truth[i]).length,
    n_susceptible: d.present.filter((p, i) => p && !truth[i]).length,
    bits: encodeBits(d.present),
  }));

  // Scores: mechanism signal, lineage leakage, noise. Holding lineages out
  // removes most of the leakage, which is the gap the dashboard reports.
  const signal = tested.map((g) => drug.mechanisms.some((m) => g.carries.has(m)));
  const leak = tested.map((g) => (enriched.has(g.lineage) ? 1 : -1));
  const lineages = [...new Set(tested.map((g) => g.lineage))];
  const nGroups = lineages.length;
  const halves = (() => {
    const out = new Array<number>(n);
    let r = 0;
    let s = 0;
    truth.forEach((t, i) => {
      out[i] = t ? r++ % 2 : s++ % 2;
    });
    return out;
  })();
  const schemeFolds = (scheme: "random" | "lineage"): number[] =>
    scheme === "random" ? tested.map(() => Math.floor(rand() * 5)) : tested.map((g) => lineages.indexOf(g.lineage) % 5);
  const strength: Record<Learner, number> = { scm: 0, cart: 2.4, rf: 3.0, gbt: 3.3, lasso: 2.7, catalog: 0 };
  const splitFor = (scheme: "random" | "lineage"): SplitResult => {
    const folds = schemeFolds(scheme);
    const leakWeight = scheme === "random" ? 0.9 : 0.15;
    const noise = scheme === "random" ? 0.9 : 1.3;
    const learners = LEARNERS.map((learner): LearnerResult => {
      const scores = tested.map((_, i) => {
        if (learner === "scm") return (rand() < (scheme === "random" ? 0.05 : 0.09) ? !signal[i] : signal[i]) ? 1 : 0;
        const z = strength[learner] * (signal[i] ? 1 : -1) + leakWeight * leak[i]! + normal(rand) * noise;
        const p = sigmoid(z);
        return learner === "cart" ? Math.round(p * 4) / 4 : p;
      });
      const ranked = learner !== "scm";
      const foldBa = [0, 1, 2, 3, 4].map((f) => {
        const idx = folds.flatMap((k, i) => (k === f ? [i] : []));
        return binaryMetrics(idx.map((i) => truth[i]!), idx.map((i) => scores[i]!), false).balanced_accuracy;
      });
      return {
        learner,
        metrics: binaryMetrics(truth, scores, ranked),
        fold_balanced_accuracy: foldBa,
        roc: ranked ? rocCurve(truth, scores) : [],
        scores,
        conformal: ranked ? conformal(truth, scores, halves) : [],
        fold_complexity: [0, 1, 2, 3, 4].map(() => (learner === "scm" ? 1 + Math.floor(rand() * 2) : learner === "cart" ? 4 + Math.floor(rand() * 5) : learner === "rf" ? 200 : learner === "gbt" ? 80 + Math.floor(rand() * 60) : 6 + Math.floor(rand() * 20))),
        seconds: rand() * 40,
      };
    });
    return { scheme, folds: 5, n_groups: nGroups, learners };
  };
  const splits: SplitResult[] = [splitFor("random")];
  if (nGroups >= 10) splits.push(splitFor("lineage"));

  // Catalog baseline, appended to every split as the pipeline does.
  const catalogCalls = drug.rule === null ? null : groupsOf.map((gs) => ruleCall(drug.rule!, gs));
  const catalogMetrics = catalogCalls === null ? null : binaryMetrics(truth, catalogCalls.map(Number), false);
  if (catalogCalls !== null && catalogMetrics !== null) {
    splits.forEach((s) =>
      s.learners.push({ learner: "catalog", metrics: catalogMetrics, fold_balanced_accuracy: [], roc: [], scores: catalogCalls.map(Number), conformal: [], fold_complexity: [], seconds: 0 }),
    );
  }
  const groupRows: GroupRow[] =
    drug.rule === null
      ? []
      : catalogGroups
          .map((g): GroupRow => {
            const carriers = groupsOf.map((gs) => gs.includes(g));
            return {
              group: g,
              class: GROUPS[g]!.class,
              subclass: GROUPS[g]!.subclass,
              counts_for_rule: drug.rule!.some((alts) => alts.some((t) => tokenMatch(t, GROUPS[g]!))),
              n_resistant: carriers.filter((c, i) => c && truth[i]).length,
              n_susceptible: carriers.filter((c, i) => c && !truth[i]).length,
            };
          })
          .sort((a, b) => b.n_resistant + b.n_susceptible - (a.n_resistant + a.n_susceptible));

  const panels: PanelResult[] = drug.panels.map((pid) => {
    const def = PANEL_DEFS[pid]!;
    const found = groupsOf.map((gs) => gs.filter((g) => def.targets.includes(g) || def.targets.includes(GROUPS[g]!.family)));
    const detected = new Map<string, [number, number]>();
    found.forEach((f, i) =>
      f.forEach((g) => {
        const e = detected.get(g) ?? [0, 0];
        e[truth[i] ? 0 : 1] += 1;
        detected.set(g, e);
      }),
    );
    return {
      ...def,
      metrics: binaryMetrics(truth, found.map((f) => (f.length > 0 ? 1 : 0)), false),
      detected: [...detected.entries()].sort(([a], [b]) => a.localeCompare(b)).map(([g, [r, s]]) => [g, r, s]),
    };
  });

  // MIC: dilution steps around the susceptible and resistant modes, censored
  // at the tested range; predictions mostly within one step.
  const mics: (Mic | null)[] = tested.map((_, i) => {
    if (drug.mic === null || rand() < 0.1) return null;
    const mode = truth[i] ? drug.mic.r : drug.mic.s;
    const v = Math.round(mode + normal(rand) * 1.1);
    if (v <= drug.mic.lo) return { log2: drug.mic.lo, censor: "AtMost" };
    if (v > drug.mic.hi) return { log2: drug.mic.hi, censor: "Above" };
    return { log2: v, censor: "Exact" };
  });
  // Genomes outside the samples whose label was inferred from MIC, which the
  // pipeline leaves out; one draw per sample with an MIC keeps the seeded
  // sequence the rest of the fixture depends on.
  const leftOut = tested.filter((_, i) => mics[i] !== null && rand() < 0.18).length;
  let mic: MicResult | null = null;
  if (drug.mic !== null) {
    const bp = drug.mic.r_above;
    const predicted = mics.map((m) => (m === null ? null : m.log2 + normal(rand) * 0.8));
    const confusion = new Map<string, number>();
    let ea = 0;
    let eaAll = 0;
    let nExact = 0;
    let has = 0;
    const catTruth: boolean[] = [];
    const catCalls: number[] = [];
    mics.forEach((m, i) => {
      if (m === null) return;
      has += 1;
      const p = Math.round(predicted[i]!);
      let ok: boolean;
      if (m.censor === "Exact") {
        nExact += 1;
        confusion.set(`${m.log2},${p}`, (confusion.get(`${m.log2},${p}`) ?? 0) + 1);
        ok = Math.abs(p - m.log2) <= 1;
        if (ok) ea += 1;
        catTruth.push(m.log2 > bp);
        catCalls.push(p > bp ? 1 : 0);
      } else ok = m.censor === "AtMost" ? p <= m.log2 + 1 : p >= m.log2;
      if (ok) eaAll += 1;
    });
    mic = {
      n: has,
      n_exact: nExact,
      essential_agreement: ea / Math.max(nExact, 1),
      essential_agreement_all: eaAll / Math.max(has, 1),
      categorical: binaryMetrics(catTruth, catCalls, false),
      breakpoint: { r_above: bp, agreement: 0.94 + rand() * 0.05, n: has },
      confusion: [...confusion.entries()].map(([k, c]) => {
        const [o, p] = k.split(",").map(Number);
        return [o!, p!, c] as [number, number, number];
      }).sort((a, b) => a[0] - b[0] || a[1] - b[1]),
      predicted: predicted.map((p) => (p === null ? null : Math.round(p * 1000) / 1000)),
    };
  }

  const samples: Sample[] = tested.map((g, i) => ({
    id: g.id,
    resistant: truth[i]!,
    from_mic: false,
    mic: mics[i] ?? null,
    lineage: g.lineage,
    country: g.country,
    year: g.year,
    host: g.host,
    mlst: g.mlst,
    groups: groupsOf[i]!.map((gr) => dict.get(gr)!),
    catalog_resistant: catalogCalls === null ? null : catalogCalls[i]!,
  }));

  // Final models. SCM: one rule per mechanism feature, as a disjunction.
  const mechDrafts = drafts.filter((d, k) => k < drafts.length - 5 && spec.mechanisms.some((m) => drug.mechanisms.includes(m.id) && d.card.locus?.gene === m.gene));
  const scmRules = mechDrafts.slice(0, 3).map((d) => ({ feature: d.id, presence: true }));
  const covered = scmRules.map((r) => drafts.find((d) => d.id === r.feature)!.present.filter((p, i) => p && truth[i]).length);
  const lost = scmRules.map((r) => drafts.find((d) => d.id === r.feature)!.present.filter((p, i) => p && !truth[i]).length);
  const cart = fitTree(drafts, truth, 3);
  const importance = drafts
    .map((d) => ({ feature: d.id, weight: Math.abs(binaryMetrics(truth, d.present.map(Number), false).balanced_accuracy - 0.5) }))
    .sort((a, b) => b.weight - a.weight);
  const report: DatasetReport = {
    schema: SCHEMA,
    id,
    species: spec.species,
    species_slug: spec.slug,
    antibiotic: drug.antibiotic,
    drug_class: drug.drug_class,
    breakpoint: drug.mic === null ? null : { r_above: drug.mic.r_above, agreement: 0.95, n: n },
    counts: {
      samples: n,
      resistant: nR,
      susceptible: n - nR,
      labels_from_mic: 0,
      labels_from_mic_excluded: leftOut,
      features: 41_000 + Math.floor(rand() * 9_000),
      patterns: 52_000 + Math.floor(rand() * 9_000),
      kmers: 2_400_000 + Math.floor(rand() * 900_000),
      k: 31,
      lineages: nGroups,
      min_support: Math.max(Math.floor(n / 200), 3),
      min_class: MIN_CLASS,
    },
    lineage_threshold: 0.005,
    catalog_groups: catalogGroups,
    samples,
    splits,
    catalog:
      drug.rule === null || catalogMetrics === null
        ? null
        : { version: "2026-08-14.1", rule: drug.rule, metrics: catalogMetrics, groups: groupRows },
    panels,
    mic,
    models: {
      scm: {
        kind: scmRules.length > 1 ? "disjunction" : "conjunction",
        p: 1.0,
        cv_error: 0.04 + rand() * 0.05,
        rules: scmRules,
        steps: scmRules.map((_, k) => ({
          covered: covered[k]!,
          lost: lost[k]!,
          utility: covered[k]! - lost[k]!,
          equivalent: k === 0 ? [{ feature: drafts[drafts.length - 1]!.id, presence: true }] : [],
          n_equivalent: k === 0 ? 3 + Math.floor(rand() * 40) : 0,
        })),
        importance: scmRules.map(() => 0.05 + rand() * 0.3),
      },
      cart: { nodes: cart },
      cart_depth: 3,
      lasso: { intercept: -1.2, coef: importance.slice(0, 4).map((w) => [w.feature, w.weight * 6] as [number, number]), lambda: 0.012 },
      rf_importance: importance,
      gbt_importance: [...importance].reverse().sort((a, b) => b.weight * 0.9 - a.weight * 0.9),
    },
    features,
    seconds: 120 + rand() * 400,
  };

  // Manifest row, mirroring src/report.rs.
  const random = splits[0]!;
  const lineage = splits[1];
  // Every fixture learner calls every isolate, so each has metrics.
  const learnerRows = LEARNERS.map((l) => ({
    learner: l,
    random: lite(random.learners.find((x) => x.learner === l)!.metrics!),
    lineage: lineage === undefined ? null : lite(lineage.learners.find((x) => x.learner === l)!.metrics!),
  }));
  const best = learnerRows.reduce((a, b) => (b.random.balanced_accuracy > a.random.balanced_accuracy ? b : a));
  const byId = new Map(features.map((f) => [f.feature, f]));
  const row: CatalogRow = {
    id,
    species: spec.species,
    species_slug: spec.slug,
    antibiotic: drug.antibiotic,
    drug_class: drug.drug_class,
    n,
    n_resistant: nR,
    labels_from_mic: report.counts.labels_from_mic,
    labels_from_mic_excluded: leftOut,
    min_class: MIN_CLASS,
    features: report.counts.features,
    lineages: nGroups,
    best_learner: best.learner,
    best_random: best.random.balanced_accuracy,
    best_lineage: best.lineage?.balanced_accuracy ?? null,
    learners: learnerRows,
    catalog: catalogMetrics === null ? null : lite(catalogMetrics),
    mic_essential_agreement: mic?.essential_agreement ?? null,
    scm_rules: scmRules.map((r) => {
      const f = byId.get(r.feature);
      return `${r.presence ? "presence of" : "absence of"} ${f?.catalog_families[0] ?? f?.locus?.gene ?? f?.locus?.product ?? `feature ${r.feature}`}`;
    }),
    panels: panels.map((p) => [p.id, lite(p.metrics)]),
    kover:
      drug.kover === null
        ? null
        : { ...drug.kover, balanced_accuracy: (drug.kover.sensitivity + drug.kover.specificity) / 2 },
  };
  const panelRows: PanelSummaryRow[] = panels.map((p) => ({
    panel: p.id,
    name: p.name,
    dataset: id,
    species: spec.species,
    antibiotic: drug.antibiotic,
    metrics: lite(p.metrics),
    best_model: best.random,
    best_learner: best.learner,
  }));
  return { report, row, panels: panelRows };
};

const buildSpecies = (spec: SpeciesSpec, genomes: Genome[], reported: number): { file: SpeciesReport; row: SpeciesRow } => {
  const dist = (i: number, j: number) => Math.hypot(genomes[i]!.x - genomes[j]!.x, genomes[i]!.y - genomes[j]!.y);
  const n = genomes.length;
  const merges = upgma(dist, n);
  const dict = new Map<string, number>();
  genomes.forEach((g) => [...g.carries].forEach((m) => GROUPS[m] !== undefined && !dict.has(m) && dict.set(m, dict.size)));
  const drugs = spec.drugs.map((d) => d.antibiotic);
  const rows: GenomeRow[] = genomes.map((g, i) => ({
    id: g.id,
    lineage: g.lineage,
    x: g.x,
    y: g.y,
    nearest: Math.min(...genomes.map((_, j) => (j === i ? Infinity : dist(i, j)))),
    country: g.country,
    year: g.year,
    host: g.host,
    mlst: g.mlst,
    length: g.length,
    contigs: g.contigs,
    n50: g.n50,
    groups: [...g.carries].filter((m) => dict.has(m)).map((m) => dict.get(m)!),
    profile: drugs.map((d) => g.calls.get(d) ?? ".").join(""),
  }));
  const cross: SpeciesReport["cross"] = [];
  drugs.forEach((a, ai) =>
    drugs.forEach((b, bi) => {
      if (bi <= ai) return;
      let rr = 0;
      let rs = 0;
      let sr = 0;
      let ss = 0;
      rows.forEach((r) => {
        const pa = r.profile[ai];
        const pb = r.profile[bi];
        if (pa === "R" && pb === "R") rr += 1;
        else if (pa === "R" && pb === "S") rs += 1;
        else if (pa === "S" && pb === "R") sr += 1;
        else if (pa === "S" && pb === "S") ss += 1;
      });
      const total = rr + rs + sr + ss;
      if (total < 30) return;
      const denom = Math.sqrt((rr + rs) * (sr + ss) * (rr + sr) * (rs + ss));
      cross.push({ a, b, n: total, rr, rs, sr, ss, phi: denom === 0 ? 0 : (rr * ss - rs * sr) / denom });
    }),
  );
  const labeled = genomes.filter((g) => [...g.calls.values()].some((c) => c !== ".")).length;
  return {
    file: {
      schema: SCHEMA,
      species: spec.species,
      slug: spec.slug,
      drugs,
      drug_classes: spec.drugs.map((d) => d.drug_class),
      catalog_groups: [...dict.keys()],
      lineage_threshold: 0.005,
      n_lineages: spec.lineages,
      leaf_order: leafOrder(merges, n),
      merges,
      genomes: rows,
      cross,
    },
    row: {
      species: spec.species,
      slug: spec.slug,
      datasets_curated: spec.drugs.length,
      datasets_reported: reported,
      genomes_labeled: labeled,
      genomes_analyzed: n,
      lineages: spec.lineages,
      lineage_threshold: 0.005,
    },
  };
};

export interface FixtureSet {
  manifest: Manifest;
  corpus: Corpus;
  species: Map<string, SpeciesReport>;
  datasets: Map<string, DatasetReport>;
}

export const buildFixture = (): FixtureSet => {
  const species = new Map<string, SpeciesReport>();
  const datasets = new Map<string, DatasetReport>();
  const rows: CatalogRow[] = [];
  const panelRows: PanelSummaryRow[] = [];
  const speciesRows: SpeciesRow[] = [];
  const allGenomes: Genome[] = [];
  for (const spec of [KPN, SAU]) {
    const rand = rng(spec.seed);
    const genomes = makeGenomes(spec, rand);
    allGenomes.push(...genomes);
    let reported = 0;
    spec.drugs.forEach((drug) => {
      const built = buildDataset(spec, drug, genomes, rand);
      if (built === null) return;
      reported += 1;
      datasets.set(built.report.id, built.report);
      rows.push(built.row);
      panelRows.push(...built.panels);
    });
    const s = buildSpecies(spec, genomes, reported);
    species.set(spec.slug, s.file);
    speciesRows.push(s.row);
  }
  const pairs = allGenomes.flatMap((g) => [...g.calls.values()]);
  const manifest: Manifest = {
    schema: SCHEMA,
    generated_at: "2026-10-23T12:00:00Z",
    k: 31,
    amr_records: 18_460,
    amrfinder_version: "2026-08-14.1",
    name_vocabulary: "SigDrug reviewed drug vocabulary",
    funnel: {
      records: 18_460,
      records_class_label: 212,
      records_not_antibiotic: 31,
      records_unresolved_name: 18,
      records_no_measure: 140,
      genomes_with_records: 301,
      genomes_without_metadata: 2,
      genomes_failed_qc: { "BV-BRC quality: Poor": 21, "CheckM completeness low": 9, "CheckM contamination high": 8, "genome length off species median": 3 },
      genomes_passing_qc: 258,
      pairs: pairs.length + 210,
      pairs_conflict: 38,
      pairs_intermediate_only: pairs.filter((c) => c === "I").length,
      pairs_reported: pairs.filter((c) => c === "R" || c === "S").length,
      pairs_labeled_from_mic: rows.reduce((a, r) => a + (r.labels_from_mic_excluded ?? 0), 0),
      pairs_unlabeled: pairs.filter((c) => c === ".").length + 172,
    },
    species: speciesRows,
    datasets: rows,
    panels: panelRows,
  };
  const count = <T,>(items: T[]): [T, number][] => {
    const m = new Map<T, number>();
    items.forEach((i) => m.set(i, (m.get(i) ?? 0) + 1));
    return [...m.entries()].sort((a, b) => b[1] - a[1]);
  };
  const corpus: Corpus = {
    schema: SCHEMA,
    by_species: [KPN, SAU].map((spec) => {
      const gs = allGenomes.filter((g) => g.id.startsWith(`${spec.taxon}.`));
      const calls = gs.flatMap((g) => [...g.calls.values()]);
      return { species: spec.species, pairs_resistant: calls.filter((c) => c === "R").length, pairs_susceptible: calls.filter((c) => c === "S").length, genomes: gs.length };
    }),
    years: count(allGenomes.map((g) => g.year)).sort((a, b) => a[0].localeCompare(b[0])),
    countries: count(allGenomes.map((g) => g.country)),
    methods: [["Broth dilution", 9_210], ["Agar dilution", 3_104], ["Disk diffusion", 2_877], ["MIC", 1_960], ["E-test", 640]],
    standards: [["CLSI", 11_820], ["EUCAST", 4_215], ["NARMS", 1_102]],
    breakpoints: [...datasets.values()]
      .filter((d) => d.breakpoint !== null)
      .map((d) => ({ species: d.species, antibiotic: d.antibiotic, r_above_log2: d.breakpoint!.r_above, r_above_mg_l: 2 ** d.breakpoint!.r_above, agreement: d.breakpoint!.agreement, n: d.breakpoint!.n, labels_from_mic: d.counts.labels_from_mic_excluded ?? 0 })),
    drug_names: [
      { records: 4_210, class: "carbapenem", raw: "meropenem", canonical: "meropenem", kind: "drug", components: [{ text: "meropenem", canonical: "meropenem", source: "sigdrug:vocabulary" }] },
      { records: 3_980, class: "quinolone", raw: "Ciprofloxacin", canonical: "ciprofloxacin", kind: "drug", components: [{ text: "ciprofloxacin", canonical: "ciprofloxacin", source: "sigdrug:vocabulary" }] },
      { records: 2_210, class: "folate pathway combination", raw: "trimethoprim/sulfamethoxazole", canonical: "sulfamethoxazole/trimethoprim", kind: "drug", components: [{ text: "trimethoprim", canonical: "trimethoprim", source: "sigdrug:vocabulary" }, { text: "sulfamethoxazole", canonical: "sulfamethoxazole", source: "sigdrug:vocabulary" }] },
      { records: 1_406, class: "penicillin + inhibitor", raw: "amoxicillin_clavulanat", canonical: "amoxicillin/clavulanate", kind: "drug", components: [{ text: "amoxicillin", canonical: "amoxicillin", source: "grm:vocabulary" }, { text: "clavulanat", canonical: "clavulanate", source: "grm:override" }] },
      { records: 380, class: "tetracycline", raw: "Tigecyklin", canonical: "tigecycline", kind: "drug", components: [{ text: "tigecyklin", canonical: "tigecycline", source: "grm:edit-distance-2" }] },
      { records: 211, class: "cephalosporin", raw: "ceftiofur sodium", canonical: "ceftiofur", kind: "drug", components: [{ text: "ceftiofur sodium", canonical: "ceftiofur", source: "sigdrug:salt_base_ingredients.tsv" }] },
      { records: 212, class: "", raw: "beta-lactam", canonical: "beta-lactam", kind: "class", components: [] },
      { records: 31, class: "", raw: "ethidium bromide", canonical: "ethidium bromide", kind: "excluded", components: [] },
      { records: 18, class: "other", raw: "xq-77", canonical: "xq-77", kind: "unresolved", components: [{ text: "xq-77", canonical: "xq-77", source: "unresolved" }] },
    ],
  };
  return { manifest, corpus, species, datasets };
};
