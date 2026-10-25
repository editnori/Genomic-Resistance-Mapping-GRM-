// Edge validation for the report files and research.json. Fetched JSON arrives as
// unknown; each parser checks the schema version, the top-level shape, and
// the cross-field counts the views index by, then returns the typed value.
// Deeper fields are trusted to the Rust serializers that share schema.ts.

import {
  READ_SCHEMAS,
  RESEARCH_SCHEMA,
  REVISION5_SNAPSHOT_SHA256,
  type Corpus,
  type DatasetReport,
  type Manifest,
  type Research,
  type SelectiveCategory,
  type SpeciesReport,
} from "./schema";

export class ShapeError extends Error {
  constructor(
    message: string,
    /** Set when the file carries a different schema version. */
    readonly foundSchema: unknown = undefined,
  ) {
    super(message);
    this.name = "ShapeError";
  }
}

type Obj = Record<string, unknown>;

const isObj = (v: unknown): v is Obj => typeof v === "object" && v !== null && !Array.isArray(v);

function need(cond: boolean, message: string): asserts cond {
  if (!cond) throw new ShapeError(message);
}

const field = (o: Obj, key: string, kind: "string" | "number" | "object" | "array", where: string): void => {
  const v = o[key];
  const ok =
    kind === "array" ? Array.isArray(v) : kind === "object" ? isObj(v) : typeof v === kind;
  need(ok, `${where}: "${key}" should be ${kind === "array" ? "an array" : `a ${kind}`}`);
};

const fields = (o: Obj, spec: Record<string, "string" | "number" | "object" | "array">, where: string): void => {
  for (const [key, kind] of Object.entries(spec)) field(o, key, kind, where);
};

/** The optional audit block: when present it must have the audit's shape. */
const checkAudit = (v: unknown, where: string): void => {
  if (v === undefined) return;
  need(
    isObj(v) && typeof v.revision === "number" && typeof v.valid === "boolean" && Array.isArray(v.reasons) && (v.scope === undefined || typeof v.scope === "string"),
    `${where}: evaluation should be {revision, valid, reasons, scope?}`,
  );
};

const checkNoCall = (v: unknown, where: string): void => {
  if (v === undefined) return;
  need(isObj(v) && typeof v.resistant === "number" && typeof v.susceptible === "number", `${where}: no_call should be {resistant, susceptible}`);
};

/** Optional counts (min_class, labels_from_mic_excluded): absent, or a whole number. */
const checkCounts = (o: Obj, keys: readonly string[], where: string): void => {
  keys.forEach((key) => {
    const v = o[key];
    need(v === undefined || (typeof v === "number" && Number.isInteger(v) && v >= 0), `${where}: "${key}" should be a whole number when present`);
  });
};

const OPTIONAL_COUNTS = ["min_class", "labels_from_mic_excluded"] as const;

const uncalledTotal = (v: unknown): number => (isObj(v) ? (v.resistant as number) + (v.susceptible as number) : 0);

/** A learner's index rates on one split: all numbers, or all null for a
 *  learner that called no isolate, which then needs no_call to count them. */
const checkLearnerLite = (v: unknown, where: string): void => {
  if (v === null) return;
  need(isObj(v), `${where} should be an object or null`);
  const rates = [v.balanced_accuracy, v.sensitivity, v.specificity, v.vme, v.me];
  const none = rates.every((r) => r === null);
  need(none || rates.every((r) => typeof r === "number"), `${where}: rates should all be numbers, or all null when no isolate is called`);
  checkNoCall(v.no_call, where);
  need(!none || uncalledTotal(v.no_call) > 0, `${where}: null rates need no_call to count the isolates without a call`);
};

const checkList = (o: Obj, key: string, spec: Record<string, "string">, where: string): void => {
  if (o[key] === undefined) return;
  field(o, key, "array", where);
  (o[key] as unknown[]).forEach((item, i) => {
    need(isObj(item), `${where}: ${key}[${i}] should be an object`);
    fields(item, spec, `${where} ${key}[${i}]`);
  });
};

const checkSchema = (raw: unknown, where: string): Obj => {
  need(isObj(raw), `${where}: expected a JSON object`);
  if (typeof raw.schema !== "number" || !READ_SCHEMAS.includes(raw.schema)) {
    throw new ShapeError(`${where} is schema ${String(raw.schema)}; this dashboard reads schemas ${READ_SCHEMAS.join(" and ")}`, raw.schema);
  }
  return raw;
};

export function parseManifest(raw: unknown): Manifest {
  const o = checkSchema(raw, "index.json");
  fields(
    o,
    { generated_at: "string", amr_records: "number", amrfinder_version: "string", funnel: "object", species: "array", datasets: "array", panels: "array" },
    "index.json",
  );
  need(o.k === null || typeof o.k === "number", 'index.json: "k" should be a number or null');
  const winners = new Map<unknown, unknown>();
  (o.datasets as unknown[]).forEach((d, i) => {
    const where = `index.json datasets[${i}]`;
    need(isObj(d), `index.json: datasets[${i}] should be an object`);
    fields(d, { id: "string", species: "string", species_slug: "string", antibiotic: "string", drug_class: "string", n: "number", n_resistant: "number", learners: "array", scm_rules: "array", panels: "array" }, where);
    checkCounts(d, OPTIONAL_COUNTS, where);
    need(d.best_random === null || typeof d.best_random === "number", `${where}: best_random should be a number or null`);
    need(d.best_learner === null || typeof d.best_learner === "string", `${where}: best_learner should be a learner or null`);
    need((d.best_learner === null) === (d.best_random === null), `${where}: best_learner and best_random should both be null only when no learner has full coverage`);
    need(d.best_random !== null || d.best_lineage === null, `${where}: no full-coverage learner cannot have a lineage score`);
    const learners = d.learners as unknown[];
    learners.forEach((l, j) => {
      need(isObj(l), `${where}: learners[${j}] should be an object`);
      checkLearnerLite(l.random, `${where} learners[${j}].random`);
      checkLearnerLite(l.lineage, `${where} learners[${j}].lineage`);
    });
    if (d.best_learner !== null) {
      // The best learner calls every isolate on each split that ran, and
      // best_random is its score. Schema-1 reports may have no lineage split.
      const row = learners.find((l) => isObj(l) && l.learner === d.best_learner) as Obj | undefined;
      const [random, lineage] = [row?.random, row?.lineage];
      const lineageRan = learners.some((l) => isObj(l) && l.lineage !== null && l.lineage !== undefined);
      need(
        isObj(random) && random.balanced_accuracy === d.best_random && uncalledTotal(random.no_call) === 0 && (!lineageRan || (isObj(lineage) && uncalledTotal(lineage.no_call) === 0)),
        `${where}: best learner ${String(d.best_learner)} should call every isolate and score best_random on random splits`,
      );
    }
    winners.set(d.id, d.best_learner);
    checkAudit(d.evaluation, where);
  });
  (o.panels as unknown[]).forEach((p, i) => {
    need(isObj(p), `index.json: panels[${i}] should be an object`);
    fields(p, { panel: "string", dataset: "string", metrics: "object", best_model: "object", best_learner: "string" }, `index.json panels[${i}]`);
    // report.rs writes panel rows only beside a best learner, and names it.
    need(winners.get(p.dataset) === p.best_learner, `index.json panels[${i}]: best_learner ${String(p.best_learner)} is not the best learner of ${String(p.dataset)}`);
    checkNoCall(p.no_call, `index.json panels[${i}]`);
  });
  checkList(o, "unsupported_panels", { panel: "string", name: "string", dataset: "string", reason: "string" }, "index.json");
  return o as unknown as Manifest;
}

export function parseCorpus(raw: unknown): Corpus {
  const o = checkSchema(raw, "corpus.json");
  fields(
    o,
    { by_species: "array", years: "array", countries: "array", methods: "array", standards: "array", breakpoints: "array", drug_names: "array" },
    "corpus.json",
  );
  return o as unknown as Corpus;
}

export function parseSpecies(raw: unknown): SpeciesReport {
  const o = checkSchema(raw, "species file");
  fields(
    o,
    { species: "string", slug: "string", drugs: "array", drug_classes: "array", catalog_groups: "array", lineage_threshold: "number", n_lineages: "number", leaf_order: "array", merges: "array", genomes: "array", cross: "array" },
    "species file",
  );
  const genomes = o.genomes as unknown[];
  const drugs = o.drugs as unknown[];
  need((o.leaf_order as unknown[]).length === genomes.length, "species file: leaf_order and genomes differ in length");
  need((o.merges as unknown[]).length < Math.max(genomes.length, 1), "species file: more merges than a tree over these genomes can have");
  genomes.forEach((g, i) => {
    need(isObj(g) && typeof g.id === "string" && typeof g.profile === "string", `species file: genomes[${i}] lacks an id or profile`);
    need((g.profile as string).length === drugs.length, `species file: genomes[${i}].profile has ${(g.profile as string).length} calls for ${drugs.length} drugs`);
  });
  return o as unknown as SpeciesReport;
}

export function parseDataset(raw: unknown): DatasetReport {
  const o = checkSchema(raw, "dataset report");
  fields(
    o,
    { id: "string", species: "string", species_slug: "string", antibiotic: "string", drug_class: "string", counts: "object", catalog_groups: "array", samples: "array", splits: "array", panels: "array", models: "object", features: "array" },
    "dataset report",
  );
  checkAudit(o.evaluation, "dataset report");
  need(o.eval_config === undefined || isObj(o.eval_config), "dataset report: eval_config should be an object");
  const n = (o.samples as unknown[]).length;
  const counts = o.counts as Obj;
  need(counts.samples === n, `dataset report: counts.samples is ${String(counts.samples)} but ${n} samples are listed`);
  const resistant = (o.samples as unknown[]).map((s) => isObj(s) && s.resistant === true);
  (o.splits as unknown[]).forEach((s, i) => {
    need(isObj(s) && Array.isArray(s.learners), `dataset report: splits[${i}] lacks learners`);
    const folds = s.test_fold;
    need(
      folds === undefined || (Array.isArray(folds) && folds.length === n && folds.every((f) => typeof f === "number" && Number.isInteger(f) && f >= 0 && f < (s.folds as number))),
      `dataset report: splits[${i}].test_fold should give each of the ${n} samples a fold below ${String(s.folds)}`,
    );
    (s.learners as unknown[]).forEach((l, j) => {
      const where = `dataset report: splits[${i}].learners[${j}]`;
      need(isObj(l) && Array.isArray(l.scores) && l.scores.length === n, `${where} scores do not cover the ${n} samples`);
      const scores = l.scores as unknown[];
      need(scores.every((v) => typeof v === "number" || v === null), `${where} scores should be numbers or null`);
      need(l.metrics === null || isObj(l.metrics), `${where}: metrics should be an object, or null when no isolate is called`);
      // A null score is an isolate without a call, so it needs no_call to count it.
      const uncalled = (r: boolean) => scores.filter((v, k) => v === null && resistant[k] === r).length;
      if (l.no_call === undefined) {
        need(uncalled(true) + uncalled(false) === 0, `${where} has scores without a call but no no_call counts`);
        need(l.metrics !== null, `${where} has null metrics but no no_call counts`);
        return;
      }
      checkNoCall(l.no_call, where);
      const x = l.no_call as Record<string, number>;
      if (l.metrics === null) {
        // No call at all: no rate, ROC point, or fold balanced accuracy exists.
        need(x.resistant === counts.resistant && x.susceptible === counts.susceptible, `${where} has null metrics but calls some isolates`);
        need(
          Array.isArray(l.roc) && l.roc.length === 0 && Array.isArray(l.fold_balanced_accuracy) && l.fold_balanced_accuracy.every((v) => v === null),
          `${where} calls no isolate, so it should have no ROC points and no fold balanced accuracy`,
        );
      } else {
        const m = l.metrics as Record<string, number>;
        need(m.n! > 0, `${where} calls no isolate, so its metrics should be null`);
        need(
          m.tp! + m.fn! + x.resistant! === counts.resistant && m.tn! + m.fp! + x.susceptible! === counts.susceptible,
          `${where} calls ${String(m.n)} isolates and gives no call for ${x.resistant! + x.susceptible!}, which does not add up to ${String(counts.resistant)} resistant and ${String(counts.susceptible)} susceptible isolates`,
        );
      }
      need(uncalled(true) === x.resistant && uncalled(false) === x.susceptible, `${where} has null scores that do not match its no_call counts`);
    });
  });
  const models = o.models as Obj;
  need(isObj(models.scm) && isObj(models.cart), "dataset report: models lacks scm or cart");
  checkCounts(counts, OPTIONAL_COUNTS, "dataset report counts");
  const min = counts.min_class as number | undefined;
  need(
    min === undefined || ((counts.resistant as number) >= min && (counts.susceptible as number) >= min),
    `dataset report: a class has fewer isolates than its class minimum of ${String(min)} (${String(counts.resistant)} resistant, ${String(counts.susceptible)} susceptible)`,
  );
  (o.panels as unknown[]).forEach((p, i) => {
    if (!isObj(p) || p.no_call === undefined) return;
    checkNoCall(p.no_call, `dataset report: panels[${i}]`);
    const m = p.metrics as Record<string, number>;
    const x = p.no_call as Record<string, number>;
    need(
      isObj(p.metrics) && m.tp! + m.fn! + x.resistant! === counts.resistant && m.tn! + m.fp! + x.susceptible! === counts.susceptible,
      `dataset report: panels[${i}] scores ${String(m.n)} isolates and gives no call for ${x.resistant! + x.susceptible!}, which does not add up to ${String(counts.resistant)} resistant and ${String(counts.susceptible)} susceptible isolates`,
    );
  });
  checkList(o, "unsupported_assays", { id: "string", name: "string", mirrors: "string", reason: "string" }, "dataset report");
  return o as unknown as DatasetReport;
}

const SELECTIVE_CATEGORIES: readonly SelectiveCategory[] = ["A1", "A2", "B", "C", "D1", "D2", "no_lineage"];

const isDigest = (v: unknown): boolean => typeof v === "string" && /^[0-9a-f]{64}$/.test(v);

const checkGuard = (g: unknown, digests: string[], where: string): void => {
  need(isObj(g) && typeof g.samples === "number", `${where}: guard should be an object with a sample count`);
  for (const key of ["sample_ids", "labels", "lineages", ...digests]) {
    need(isDigest(g[key]) || (key === "mic_predicted" && g[key] === null), `${where}: guard.${key} should be a SHA-256 digest`);
  }
};

const checkRun = (run: unknown, where: string): void => {
  need(
    isObj(run) && typeof run.n_predicted === "number" && typeof run.n_exact === "number" && (run.ea_exact === null || typeof run.ea_exact === "number"),
    `${where} should be {n_predicted, n_exact, ea_exact} with ea_exact a number or null`,
  );
};

/** research.json. Beyond the shape, it checks the pinned snapshot and the
 *  rules the views rely on: a MIC model that scores zone diameters appears
 *  only in mic.withheld, as a reason with no numbers; a rerun MIC model has
 *  both splits and a baseline scored on the same isolates, while a saved
 *  report model has neither; and only selective-calling categories A1 and
 *  A2 carry numbers. */
export function parseResearch(raw: unknown): Research {
  need(isObj(raw), "research.json: expected a JSON object");
  if (raw.schema !== RESEARCH_SCHEMA) {
    throw new ShapeError(`research.json is schema ${String(raw.schema)}; this dashboard reads research schema ${RESEARCH_SCHEMA}`, raw.schema);
  }
  fields(raw, { snapshot: "object", sources: "array", mic: "object", selective: "object" }, "research.json");
  const snapshot = (raw.snapshot as Obj).sha256;
  need(snapshot === REVISION5_SNAPSHOT_SHA256, `research.json was built from snapshot ${String(snapshot)}; this dashboard pins revision-5 snapshot ${REVISION5_SNAPSHOT_SHA256}`);

  const mic = raw.mic as Obj;
  fields(mic, { limits: "object", datasets: "array", withheld: "array" }, "research.json mic");
  fields(mic.limits as Obj, { scope: "string", baseline: "string", exact_minimum: "string", saved_only: "string", zone_fix: "string" }, "research.json mic.limits");
  need(((mic.limits as Obj).zone_fix as string).includes("{removed}"), "research.json mic.limits.zone_fix should hold {removed} for the removed zone diameters");
  const micIds = new Set<string>();
  (mic.datasets as unknown[]).forEach((d, i) => {
    const where = `research.json mic.datasets[${i}]`;
    need(isObj(d) && typeof d.id === "string" && !micIds.has(d.id), `${where}: missing or repeated id`);
    micIds.add(d.id as string);
    checkGuard(d.guard, ["mic_values", "mic_predicted"], where);
    const zone = d.zone_diameters;
    need(isObj(zone) && typeof zone.removed === "number", `${where}: zone_diameters should be {removed, scored}`);
    // A regenerated file that ships a zone-scored model gives no research
    // rather than its numbers.
    need(zone.scored === 0, `${where}: its model scores ${String(zone.scored)} zone diameters as MICs, so it belongs in mic.withheld`);
    const model = d.model;
    need(isObj(model), `${where}: model should be an object`);
    checkRun(model.random, `${where} model.random`);
    if (model.source === "rerun") {
      need(isObj(d.baseline), `${where}: a rerun needs its baseline`);
      checkRun(model.lineage, `${where} model.lineage`);
      checkRun(d.baseline.random, `${where} baseline.random`);
      checkRun(d.baseline.lineage, `${where} baseline.lineage`);
      // The baseline is scored on the model's folds, so a view gives the two
      // one denominator.
      for (const split of ["random", "lineage"] as const) {
        const m = model[split] as Obj;
        const b = (d.baseline as Obj)[split] as Obj;
        need(
          m.n_exact === b.n_exact && m.n_predicted === b.n_predicted,
          `${where}: on the ${split} split the model scores ${String(m.n_exact)} exact MICs of ${String(m.n_predicted)} predictions and the baseline ${String(b.n_exact)} of ${String(b.n_predicted)}, but both score the same isolates`,
        );
      }
    } else {
      need(
        model.source === "saved_report" && model.lineage === null && d.baseline === null && zone.removed === 0,
        `${where}: a model is either a rerun or a saved report model with no lineage run, no baseline, and no removed MICs`,
      );
    }
  });
  (mic.withheld as unknown[]).forEach((w, i) => {
    need(
      isObj(w) && typeof w.id === "string" && !micIds.has(w.id) && typeof w.reason === "string" && w.reason !== "",
      `research.json mic.withheld[${i}]: needs a reason and an id with no other MIC entry`,
    );
    micIds.add(w.id as string);
  });

  const selective = raw.selective as Obj;
  fields(selective, { categories: "array", limits: "object", datasets: "array" }, "research.json selective");
  fields(selective.limits as Obj, { scope: "string", bar: "string", no_call: "string", optimism: "string", procedure: "string", per_100: "string" }, "research.json selective.limits");
  const selectiveIds = new Set<string>();
  (selective.datasets as unknown[]).forEach((d, i) => {
    const where = `research.json selective.datasets[${i}]`;
    need(isObj(d) && typeof d.id === "string" && !selectiveIds.has(d.id), `${where}: missing or repeated id`);
    selectiveIds.add(d.id as string);
    checkGuard(d.guard, ["scores"], where);
    const category = d.category as SelectiveCategory;
    need(SELECTIVE_CATEGORIES.includes(category), `${where}: unknown category ${String(d.category)}`);
    if (category !== "A1" && category !== "A2") {
      need(d.result === undefined, `${where}: category ${category} carries no numbers, but the entry has a result`);
      return;
    }
    const r = d.result;
    need(isObj(r), `${where}: category ${category} needs its result`);
    need(category === "A2" ? r.learner === "gbt" || r.learner === "rf" : r.learner === undefined, `${where}: an A2 result names its refit learner (gbt or rf); an A1 result names none`);
    const n = r as Record<string, number>;
    need(
      [n.n, n.n_r, n.n_s, n.tp, n.fn, n.fp, n.tn].every((v) => Number.isInteger(v)) && n.n === n.n_r! + n.n_s! && n.tp! + n.fn! <= n.n_r! && n.fp! + n.tn! <= n.n_s!,
      `${where}: the result's counts do not partition its isolates`,
    );
  });
  const categories = selective.categories as unknown[];
  need(
    categories.length === SELECTIVE_CATEGORIES.length &&
      categories.every((c, i) => isObj(c) && c.code === SELECTIVE_CATEGORIES[i] && typeof c.label === "string" && c.label.includes("{learner}") === (c.code === "A2")),
    "research.json: selective.categories does not list every category in order with its label, where only A2's holds {learner}",
  );
  return raw as unknown as Research;
}
