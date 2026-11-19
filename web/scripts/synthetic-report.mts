// Writes a synthetic GRM report tree (index.json, corpus.json, species/,
// datasets/) that `vite preview` or `grm serve` can serve and
// web/scripts/browser-check.cjs can read. Everything comes from the web app's
// committed seeded generator, web/src/dev/fixtures/generate.ts: nothing in it
// is a measurement. probeRows.ts, which holds real probe-panel results, is not
// imported.
//
// Run it from the repository root with Node 22. tsx is not a repo
// dependency: npx downloads it from the npm registry, a tool download, not
// data. tsx keeps its transform cache and an IPC pipe in /tmp/tsx-<uid>/;
// the script itself writes only into <out-dir>. Anything but the two paths
// and one --shape is refused with rc 2.
//
//   npx --yes tsx@4.23.15 web/scripts/synthetic-report.mts web <out-dir> [--shape base|r11|wide|full]
//   (cd <out-dir> && sha256sum -c SOURCE.sha256)
//
// Then serve the tree on 127.0.0.1, either way, and run the check against it
// with playwright-core installed outside the app (GRM_CHROMIUM unset uses
// Playwright's own Chromium). Resolve <out-dir> before the cd into web/,
// where a relative path would miss the tree, and never leave GRM_REPORT_DIR
// out: without it vite serves ../data/report, the real report tree.
//
//   out=$(realpath <out-dir>)
//   (cd web && npm run build && GRM_REPORT_DIR="$out" npm run preview -- --port 4173 --strictPort) &
//   #  or: mkdir -p <root> && ln -s "$out" <root>/report && grm --data <root> serve --app web/dist --port 4173 &
//   node web/scripts/browser-check.cjs <dir>/node_modules/playwright-core <receipt-dir> http://127.0.0.1:4173 "$out"
//
// Identifiers are made generic before writing, so no row can be read as a
// real isolate: genome ids (BV-BRC-shaped "573.10037" in the generator)
// become SYN-<species>-<n>, contigs SYN_CONTIG_<n>.1, and feature ids
// fig|SYN.<n>.peg.<m>. The mapping is consistent across species and dataset
// files.
//
// --shape base (default): the manifest exactly as buildFixture() writes it,
//   which predates the evaluation audit and the panel counts.
// --shape r11: adds what src/report.rs writes at EVALUATION_REVISION 11,
//   derived from the fixture's own files: evaluation_revision, each report's
//   evaluation audit (revision 11, valid, with analysis.rs's scope text),
//   copied to its catalog row, and each panel row's counts (the panel's
//   confusion), best_model_counts (the best learner's random-split confusion)
//   and no_call. The audit is a fixture value, not the result of an audit.
// --shape wide: r11 plus four generically named copies of each species, so
//   the Datasets table has 30 rows, past its 25-row page (see the wide block).
// --shape full: wide plus one dataset for each state the wide tree never
//   reaches, 35 rows in all (see the full block): a lineage split whose fold
//   lacks a class, one with too few lineages, a dataset with no catalog rule,
//   an assay that is not scored, and a single rule long enough to clip at
//   390px. One of them also has a catalog rule below 50%, so the Overview's
//   catalog figure starts under 50% and draws its 50% line. The check's
//   state rows then find all four states, and it asserts the unscored assay,
//   the clipped rule and the 50% line. The full tree also carries no-genome
//   baselines in the format docs/baselines.md gives: baselines/<id>.json for
//   34 of its 35 datasets and baselines.json pinned to its index.json. They
//   are fixture values, not `grm baselines` output (see the baselines block):
//   one model level with its constant, one below it, one dataset whose MIC
//   folds the predictions do not confirm, and one MIC dataset with no file.
//   The other three shapes carry none, so the check also reads a tree
//   without them.
//
// The tree is written to an empty or new directory, with SYNTHETIC.txt and
// SOURCE.sha256 (sha256sum lines for every JSON file, sorted by path).

import { createHash } from "node:crypto";
import { existsSync, mkdirSync, readdirSync, writeFileSync } from "node:fs";
import { join, relative, resolve } from "node:path";
import { pathToFileURL } from "node:url";

const SHAPES = ["base", "r11", "wide", "full"];
const args = process.argv.slice(2);
const shapeAt = args.indexOf("--shape");
const shape = shapeAt === -1 ? "base" : args.splice(shapeAt, 2)[1];
const [webDir, outDir, ...extra] = args;
// Two paths and one known --shape, nothing else: a typo such as `--shap full`
// is refused rather than read as the base shape.
if (webDir === undefined || outDir === undefined || extra.length > 0 || [webDir, outDir].some((a) => a.startsWith("-")) || shape === undefined || !SHAPES.includes(shape)) {
  console.error("usage: tsx web/scripts/synthetic-report.mts <web-dir> <out-dir> [--shape base|r11|wide|full]");
  process.exit(2);
}
const out = resolve(outDir);
if (existsSync(out) && readdirSync(out).length > 0) {
  console.error(`${out} is not empty`);
  process.exit(2);
}

const generator = pathToFileURL(resolve(webDir, "src/dev/fixtures/generate.ts")).href;
const { buildFixture } = (await import(generator)) as typeof import("../src/dev/fixtures/generate");
const stats = pathToFileURL(resolve(webDir, "src/dev/fixtures/stats.ts")).href;
const { binaryMetrics } = (await import(stats)) as typeof import("../src/dev/fixtures/stats");
// The app's copy of the index's rule words (report.rs's rule_text), so a rule
// this script changes is worded as the index of this tree's wording words it.
const derive = pathToFileURL(resolve(webDir, "src/data/derive/dataset.ts")).href;
const { indexRuleWords } = (await import(derive)) as typeof import("../src/data/derive/dataset");
const fixture = buildFixture();

// Generic identifiers, assigned in first-seen order so reruns match.
const genomeIds = new Map<string, string>();
const genome = (id: string): string => {
  let g = genomeIds.get(id);
  if (g === undefined) {
    const taxon = id.split(".")[0] ?? "0";
    g = `SYN-${taxon === "573" ? "A" : taxon === "1280" ? "B" : taxon}-${String(genomeIds.size + 1).padStart(4, "0")}`;
    genomeIds.set(id, g);
  }
  return g;
};
const contigIds = new Map<string, string>();
const contig = (id: string): string => {
  let c = contigIds.get(id);
  if (c === undefined) {
    c = `SYN_CONTIG_${String(contigIds.size + 1).padStart(4, "0")}.1`;
    contigIds.set(id, c);
  }
  return c;
};
let featureN = 0;
const featureId = (): string => `fig|SYN.${++featureN}.peg.${featureN}`;

type Json = null | boolean | number | string | Json[] | { [k: string]: Json };

const species = new Map<string, Json>();
for (const [slug, file] of fixture.species) {
  const f = JSON.parse(JSON.stringify(file)) as { genomes: { id: string }[] };
  f.genomes.forEach((g) => (g.id = genome(g.id)));
  species.set(slug, f as unknown as Json);
}

interface Locus { genome_id: string; contig: string; patric_id: string | null; gene: string | null; product: string | null }
interface Split { scheme: string; folds: number; n_groups: number; learners: { learner: string; metrics: Record<string, number> | null; scores: number[] }[]; unavailable_reason?: string }
interface CatalogGroup { group: string; class: string; subclass: string; counts_for_rule: boolean }
interface DatasetJson {
  id: string;
  samples: { id: string; resistant: boolean; groups: number[]; catalog_resistant: boolean | null }[];
  features: { feature: number; catalog_families: string[]; catalog_products: string[]; locus: Locus | null }[];
  panels: { id: string; metrics: Record<string, number>; no_call?: unknown }[];
  splits: Split[];
  catalog_groups: string[];
  catalog: { rule: string[][]; metrics: Record<string, number>; groups: CatalogGroup[] } | null;
  unsupported_assays?: { id: string; name: string; mirrors: string; reason: string }[];
  models: { scm: { kind: string; rules: { feature: number; presence: boolean }[]; steps: unknown[]; importance: number[] } };
  evaluation?: unknown;
}
const datasets = new Map<string, DatasetJson>();
for (const [id, report] of fixture.datasets) {
  const d = JSON.parse(JSON.stringify(report)) as DatasetJson;
  d.samples.forEach((s) => (s.id = genome(s.id)));
  d.features.forEach((f) => {
    if (f.locus === null) return;
    f.locus.genome_id = genome(f.locus.genome_id);
    f.locus.contig = contig(f.locus.contig);
    if (f.locus.patric_id !== null) f.locus.patric_id = featureId();
  });
  datasets.set(id, d);
}

interface Row {
  id: string;
  species_slug: string;
  antibiotic: string;
  evaluation?: unknown;
  best_learner: string | null;
  best_lineage: number | null;
  lineage_unavailable_reason?: string;
  learners: { lineage: unknown }[];
  catalog: unknown;
  scm_rules: string[];
}
const manifest = JSON.parse(JSON.stringify(fixture.manifest)) as {
  evaluation_revision?: number;
  species: { species: string; slug: string; datasets_reported: number }[];
  datasets: Row[];
  panels: { dataset: string; panel: string; species: string; counts?: unknown; best_model_counts?: unknown; no_call?: unknown }[];
  unsupported_panels?: { panel: string; name: string; dataset: string; species: string; antibiotic: string; reason: string }[];
};

if (shape !== "base") {
  // src/analysis.rs: EVALUATION_REVISION and the audit's scope.
  const REVISION = 11;
  const SCOPE =
    "Internal cross-validation on public BV-BRC genomes scored against reported R/S labels. Feature selection and tuning use only each outer fold's training samples; inner tuning shares that outer-training feature set. MIC thresholds are fitted within training folds. This is not external or clinical validation.";
  const confusion = (m: Record<string, number> | null | undefined) =>
    m == null ? undefined : { n: m.n, tp: m.tp, fp: m.fp, tn: m.tn, fn: m.fn };
  manifest.evaluation_revision = REVISION;
  for (const d of datasets.values()) d.evaluation = { revision: REVISION, valid: true, reasons: [], scope: SCOPE };
  for (const row of manifest.datasets) row.evaluation = datasets.get(row.id)!.evaluation;
  for (const p of manifest.panels) {
    const d = datasets.get(p.dataset)!;
    const panel = d.panels.find((x) => x.id === p.panel);
    const row = manifest.datasets.find((x) => x.id === p.dataset)!;
    const best = d.splits.find((s) => s.scheme === "random")?.learners.find((l) => l.learner === row.best_learner);
    p.counts = confusion(panel?.metrics);
    p.best_model_counts = confusion(best?.metrics);
    if (panel?.no_call !== undefined) p.no_call = panel.no_call;
  }
}

const corpus = JSON.parse(JSON.stringify(fixture.corpus)) as { by_species: { species: string }[]; breakpoints: { species: string; antibiotic: string }[] };

// Generic copies of the generator's two species (letter A for the first, B
// for the second): its species file, dataset reports, catalog and panel rows
// and corpus rows, renamed to another slug, species name and genome prefix.
const originals = manifest.species.map((s) => ({ ...s, letter: s.slug === "klebsiella_pneumoniae" ? "A" : s.slug === "staphylococcus_aureus" ? "B" : s.slug }));
type Original = (typeof originals)[number];
const baseDatasets = [...datasets.values()];
const baseRows = [...manifest.datasets];
const basePanels = [...manifest.panels];
const baseBy = [...corpus.by_species];
const baseBp = [...corpus.breakpoints];
const copySpecies = (from: Original, name: string, keep: (antibiotic: string) => boolean): { slug: string; rows: Row[] } => {
  const slug = `synthetic_species_${name.toLowerCase()}`;
  const copy = <T,>(value: T): T =>
    JSON.parse(
      JSON.stringify(value)
        .replaceAll(`"${from.slug}`, `"${slug}`)
        .replaceAll(`"${from.species}"`, `"Synthetic species ${name}"`)
        .replaceAll(`"SYN-${from.letter}-`, `"SYN-${name}-`),
    ) as T;
  const rows = baseRows.filter((r) => r.species_slug === from.slug && keep(r.antibiotic)).map(copy);
  const ids = new Set(rows.map((r) => r.id));
  species.set(slug, copy(species.get(from.slug)!));
  for (const d of baseDatasets.filter((x) => x.id.startsWith(`${from.slug}__`))) {
    const c = copy(d);
    if (ids.has(c.id)) datasets.set(c.id, c);
  }
  manifest.datasets.push(...rows);
  manifest.panels.push(...basePanels.filter((p) => p.species === from.species).map(copy).filter((p) => ids.has(p.dataset)));
  manifest.species.push({ ...copy(manifest.species.find((s) => s.slug === from.slug)!), datasets_reported: rows.length });
  corpus.by_species.push(...baseBy.filter((b) => b.species === from.species).map(copy));
  corpus.breakpoints.push(...baseBp.filter((b) => b.species === from.species && keep(b.antibiotic)).map(copy));
  return { slug, rows };
};

if (shape === "wide" || shape === "full") {
  // r11 plus copies k = 2..5 of each species under a generic name
  // ("Synthetic species A2", slug synthetic_species_a2, genomes SYN-A2-*).
  // 6 datasets become 30, past the Datasets table's 25-row page, so
  // "Show all" exists. The funnel is not scaled: its totals stay the base set's.
  for (let k = 2; k <= 5; k += 1) for (const from of originals) copySpecies(from, `${from.letter}${k}`, () => true);
}

if (shape === "full") {
  // One more dataset for each state below, on two more copies: "Synthetic
  // species C" (from A) with four of A's datasets and "Synthetic species D"
  // (from B) with one of B's. Each change is the shape src/report.rs and
  // src/analysis.rs write for that state; the numbers stay the generator's.
  const [a, b] = originals as [Original, Original];
  const C = copySpecies(a, "C", (x) => ["meropenem", "ceftriaxone", "ciprofloxacin", "tetracycline"].includes(x));
  const D = copySpecies(b, "D", (x) => x === "oxacillin");
  const at = (rows: Row[], antibiotic: string) => {
    const row = rows.find((r) => r.antibiotic === antibiotic)!;
    return { row, d: datasets.get(row.id)! };
  };
  // A lineage split that could not run: no learners (and so no catalog
  // learner), and the reason on the split and its row (eval.rs, report.rs).
  const noLineageSplit = ({ row, d }: { row: Row; d: DatasetJson }, reason: string) => {
    const random = d.splits.find((s) => s.scheme === "random")!;
    d.splits = [random, { scheme: "lineage", folds: random.folds, n_groups: random.n_groups, learners: [], unavailable_reason: reason }];
    row.best_lineage = null;
    row.learners.forEach((l) => (l.lineage = null));
    row.lineage_unavailable_reason = reason;
  };
  // "Fold lacks a class": learn/mod.rs's fold check, under eval.rs's context.
  noLineageSplit(at(C.rows, "tetracycline"), "invalid lineage outer evaluation: fold 2 training lacks one or both classes");
  // "Too few lineages": B has 8 lineages, under the 10 that 5 lineage folds
  // need (eval.rs's check_split_plan). The generator writes no lineage split
  // for B at all; the pipeline would write this one.
  noLineageSplit(at(D.rows, "oxacillin"), "fewer than twice as many lineages as requested folds");
  // "No rule": the catalog has no rule for the drug, so there is no catalog
  // baseline, catalog learner or catalog call.
  {
    const { row, d } = at(C.rows, "ceftriaxone");
    d.catalog = null;
    d.splits.forEach((s) => (s.learners = s.learners.filter((l) => l.learner !== "catalog")));
    d.samples.forEach((s) => (s.catalog_resistant = null));
    row.catalog = null;
  }
  // An assay that is not scored: listed on its report and in the manifest's
  // unsupported_panels (report.rs) instead of as a panel.
  {
    const { row, d } = at(C.rows, "ciprofloxacin");
    const assay = { id: "syn_qrdr_probes", name: "Synthetic QRDR probe assay", mirrors: "gyrA and parC QRDR probe assays", reason: "synonymous base changes read as mutations" };
    d.unsupported_assays = [assay];
    manifest.unsupported_panels = [{ panel: assay.id, name: assay.name, dataset: row.id, species: "Synthetic species C", antibiotic: row.antibiotic, reason: assay.reason }];
  }
  // A single rule long enough to clip at 390px: the set covering machine's
  // first rule alone, on a feature named by a long product and no gene or
  // catalog family, so report.rs's rule_text names it by the product. The
  // words are indexRuleWords' for the changed report, so they carry what the
  // tree's wording adds (the step's ties) and match the index's
  // rule_text_version.
  {
    const { row, d } = at(C.rows, "meropenem");
    const scm = d.models.scm;
    const rule = scm.rules[0]!;
    const f = d.features.find((x) => x.feature === rule.feature)!;
    const product = "synthetic multidrug efflux RND transporter periplasmic adaptor subunit, a long fixture name for a rule that clips";
    f.catalog_families = [];
    f.catalog_products = [];
    if (f.locus !== null) Object.assign(f.locus, { gene: null, product });
    scm.rules = [rule];
    scm.steps = scm.steps.slice(0, 1);
    scm.importance = scm.importance.slice(0, 1);
    scm.kind = "conjunction";
    row.scm_rules = [indexRuleWords(d as unknown as Parameters<typeof indexRuleWords>[0], 0)];
    if (!row.scm_rules[0]!.includes(product)) throw new Error(`the long rule does not name its product: ${row.scm_rules[0]}`);
  }
  // A catalog rule below 50%, on the same meropenem copy: a too-broad rule,
  // any gene of the BETA-LACTAM subclass, which the intrinsic SHV gives
  // nearly every isolate. It calls almost all of them resistant, so its
  // balanced accuracy is just under 50% and the Overview's catalog figure
  // starts below 50% and draws its 50% line (CatalogScatter.tsx). The calls
  // follow generate.ts's ruleCall over each sample's catalog groups; the
  // metrics, catalog learner, group table and row follow from the calls.
  {
    const { row, d } = at(C.rows, "meropenem");
    const catalog = d.catalog!;
    const rule = [["S:BETA-LACTAM"]];
    const byName = new Map(catalog.groups.map((g) => [g.group, g]));
    const matches = (token: string, g: CatalogGroup) => (token.startsWith("C:") ? g.class : g.subclass).split("/").some((t) => t.trim() === token.slice(2));
    const calls = d.samples.map((s) => rule.every((alts) => s.groups.some((i) => alts.some((t) => matches(t, byName.get(d.catalog_groups[i]!)!)))));
    const metrics = binaryMetrics(d.samples.map((s) => s.resistant), calls.map(Number), false) as unknown as Record<string, number>;
    catalog.rule = rule;
    catalog.metrics = metrics;
    catalog.groups.forEach((g) => (g.counts_for_rule = rule.some((alts) => alts.some((t) => matches(t, g)))));
    d.samples.forEach((s, i) => (s.catalog_resistant = calls[i]!));
    d.splits.forEach((s) => s.learners.filter((l) => l.learner === "catalog").forEach((l) => Object.assign(l, { metrics, scores: calls.map(Number) })));
    const { balanced_accuracy, sensitivity, specificity, vme, me, auroc } = metrics;
    row.catalog = { balanced_accuracy, sensitivity, specificity, vme, me, auroc };
  }
}

const files = new Map<string, Json>([
  ["index.json", manifest as unknown as Json],
  ["corpus.json", corpus as unknown as Json],
  ...[...species].map(([slug, f]) => [`species/${slug}.json`, f] as [string, Json]),
  ...[...datasets].map(([id, d]) => [`datasets/${id}.json`, d as unknown as Json] as [string, Json]),
]);

// Refuse to write if any generator-shaped identifier survived the renaming.
const leftover = /"\d{3,5}\.\d{4,6}"|NZ_CP0\d+|fig\|\d+\./;
for (const [name, value] of files) {
  const text = JSON.stringify(value);
  const hit = leftover.exec(text);
  if (hit !== null) {
    console.error(`${name} still holds a generator identifier: ${hit[0]}`);
    process.exit(1);
  }
}

// The bytes of every file, as written: JSON and a final newline.
const texts = new Map<string, string>([...files].map(([name, value]) => [name, `${JSON.stringify(value)}\n`]));
const sha256 = (text: string): string => createHash("sha256").update(text).digest("hex");

if (shape === "full") {
  // SYNTHETIC no-genome baselines (docs/baselines.md, "The file" and "The
  // summary"): baselines/<id>.json for every dataset but one, and
  // baselines.json beside the index. `grm baselines` cannot write them here:
  // it rebuilds the MIC folds and checks them against the held-out
  // predictions, which the seeded generator does not draw the way a fit
  // does. So this block writes the documented fields itself, as fixture
  // values: the model's count is the report's own, and every other count
  // is chosen here. Each file pins its report's bytes and revision, and the
  // summary pins the bytes of index.json, so the pages take them as they
  // would take the command's. The choices give the cases the browser check
  // reads: a model level with its constant and a learner level with the
  // lineage baseline (a2 meropenem), a model below its constant (a3
  // meropenem), folds the predictions do not confirm (a4 meropenem), folds
  // where the model predicts only its start value (a5 ciprofloxacin), copies
  // in training (the first meropenem dataset), and one MIC dataset with no
  // file, so the summary covers 11 of 12 (C ciprofloxacin).
  const NO_FILE = "synthetic_species_c__ciprofloxacin";
  const LEVEL = "synthetic_species_a2__meropenem";
  const ABOVE = "synthetic_species_a3__meropenem";
  const UNCONFIRMED = "synthetic_species_a4__meropenem";
  const START_VALUE = "synthetic_species_a5__ciprofloxacin";
  const COPIES = "klebsiella_pneumoniae__meropenem";
  // Compact JSON, keys sorted, no final newline; a whole float keeps its point.
  const FLOATS = new Set(["essential_agreement", "sensitivity", "specificity", "balanced_accuracy", "unseen_share", "model", "constant"]);
  const compact = (v: unknown, key?: string): string => {
    if (Array.isArray(v)) return `[${v.map((x) => compact(x)).join(",")}]`;
    if (v !== null && typeof v === "object") {
      const o = v as Record<string, unknown>;
      return `{${Object.keys(o).sort().map((k) => `${JSON.stringify(k)}:${compact(o[k], k)}`).join(",")}}`;
    }
    if (typeof v === "number" && key !== undefined && FLOATS.has(key) && Number.isInteger(v)) return v.toFixed(1);
    return JSON.stringify(v);
  };
  const ratio = (x: number, n: number) => (n === 0 ? 0 : x / n);
  const agreement = (agree: number, n: number) => ({ agree, n, essential_agreement: agree / Math.max(n, 1) });
  const lineageEntry = (scheme: string, folds: number, n: number, [tp, fn, tn, fp]: number[], unseen: number) => {
    const sensitivity = ratio(tp!, tp! + fn!), specificity = ratio(tn!, tn! + fp!);
    return { scheme, folds, n, tp, fn, tn, fp, sensitivity, specificity, balanced_accuracy: (sensitivity + specificity) / 2, unseen, unseen_share: unseen / Math.max(n, 1) };
  };
  type Report = DatasetJson & { evaluation: { revision: number }; mic: { n_exact: number; essential_agreement: number } | null; samples: { resistant: boolean; mic: unknown }[] };
  const entries: { confirmed: boolean; constant: number; dataset: string; model: number; report_sha256: string }[] = [];
  for (const row of manifest.datasets) {
    if (row.id === NO_FILE) continue;
    const r = datasets.get(row.id) as unknown as Report;
    const reportSha = sha256(texts.get(`datasets/${row.id}.json`)!);
    const n = r.samples.length, resistant = r.samples.filter((x) => x.resistant).length, susceptible = n - resistant;
    const lineage = r.splits.filter((x) => typeof x.folds === "number").map((x) => {
      if (x.scheme !== "random") return lineageEntry(x.scheme, x.folds, n, [0, resistant, susceptible, 0], n);
      const scm = x.learners.find((l) => l.learner === "scm")?.metrics;
      const counts = row.id === LEVEL && scm ? [scm.tp!, scm.fn!, scm.tn!, scm.fp!] : [Math.round(resistant * 0.62), resistant - Math.round(resistant * 0.62), Math.round(susceptible * 0.74), susceptible - Math.round(susceptible * 0.74)];
      return lineageEntry("random", x.folds, n, counts, Math.round(n * 0.09));
    });
    let mic = null;
    if (r.mic !== null) {
      const nExact = r.mic.n_exact, modelAgree = Math.round(r.mic.essential_agreement * nExact);
      if (modelAgree / nExact !== r.mic.essential_agreement) throw new Error(`${row.id}: the report's MIC figure is not a count over its exact MICs`);
      const constantAgree = row.id === LEVEL ? modelAgree : row.id === ABOVE ? Math.min(nExact, modelAgree + 4) : modelAgree - (row.antibiotic === "meropenem" ? 15 : 9);
      const folds = 5, withMic = r.samples.filter((x) => x.mic !== null).length;
      const tests = Array.from({ length: folds }, (_, i) => Math.floor(nExact / folds) + (i < nExact % folds ? 1 : 0));
      const confirmed = row.id !== UNCONFIRMED;
      mic = {
        folds, seed: 42,
        prediction_check: { confirmed, groups: confirmed ? 4 : 3, isolates: confirmed ? 9 : 7, spanning: confirmed ? 0 : 1, constant_folds: row.id === START_VALUE ? [2] : [] },
        threshold_check: { confirmed: true, checked: folds },
        exact_with_training_copy: row.id === COPIES ? 2 : 0,
        n: withMic, n_exact: nExact,
        model: agreement(modelAgree, nExact), median: agreement(constantAgree, nExact), mode: agreement(Math.max(0, constantAgree - 3), nExact),
        fold_constants: tests.map((t, i) => ({ train_exact: nExact - t, test_exact: t, median: i % 2, upper_median_differs: false, mode: i % 2, mode_tied: false })),
      };
      entries.push({ confirmed, constant: mic.median.essential_agreement, dataset: row.id, model: mic.model.essential_agreement, report_sha256: reportSha });
    }
    texts.set(`baselines/${row.id}.json`, compact({ schema: 1, dataset: row.id, report_sha256: reportSha, evaluation_revision: r.evaluation.revision, mic, lineage }));
  }
  texts.set("baselines.json", compact({ schema: 1, evaluation_revision: manifest.evaluation_revision, index_sha256: sha256(texts.get("index.json")!), datasets: entries }));
  mkdirSync(join(out, "baselines"), { recursive: true });
}

mkdirSync(join(out, "species"), { recursive: true });
mkdirSync(join(out, "datasets"), { recursive: true });
const sums: string[] = [];
// Byte order, as `LC_ALL=C sort` lists them.
for (const [name, text] of [...texts].sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0))) {
  writeFileSync(join(out, name), text);
  sums.push(`${sha256(text)}  ./${name}`);
}
const baselineFiles = [...texts.keys()].filter((name) => name.startsWith("baselines")).length;
writeFileSync(
  join(out, "SYNTHETIC.txt"),
  `Synthetic GRM report tree, shape ${shape}, from web/src/dev/fixtures/generate.ts (seeded; nothing here is a measurement).\n` +
    `Genome, contig and feature identifiers are generic (SYN-*). Not a release, not real data.\n` +
    (baselineFiles === 0 ? "" : `The ${baselineFiles} baselines files are fixture values written by web/scripts/synthetic-report.mts, not the output of grm baselines.\n`) +
    `Files: ${texts.size}; genomes renamed: ${genomeIds.size}; contigs: ${contigIds.size}; feature ids: ${featureN}.\n`,
);
writeFileSync(join(out, "SOURCE.sha256"), `${sums.join("\n")}\n`);
console.log(`wrote ${texts.size} files to ${relative(process.cwd(), out) || "."} (shape ${shape}; ${genomeIds.size} genomes, ${manifest.datasets.length} datasets, ${manifest.panels.length} panel rows${baselineFiles === 0 ? "" : `, ${baselineFiles} baselines files`})`);
