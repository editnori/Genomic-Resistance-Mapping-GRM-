import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { ShapeError, parseCorpus, parseDataset, parseManifest, parseSpecies } from "./validate";
import type { Manifest } from "./schema";

const manifest = () => fixtureFile("/data/index.json") as Record<string, unknown>;

describe("edge validation", () => {
  it("accepts every fixture file", () => {
    const m: Manifest = parseManifest(manifest());
    expect(m.datasets.length).toBeGreaterThan(3);
    expect(() => parseCorpus(fixtureFile("/data/corpus.json"))).not.toThrow();
    m.species.forEach((s) => expect(() => parseSpecies(fixtureFile(`/data/species/${s.slug}.json`))).not.toThrow());
    m.datasets.forEach((d) => expect(() => parseDataset(fixtureFile(`/data/datasets/${d.id}.json`))).not.toThrow());
  });

  it("reads schema 1 and schema 2, and rejects another version with the version it found", () => {
    const id = parseManifest(manifest()).datasets[0]!.id;
    for (const schema of [1, 2]) {
      expect(parseManifest({ ...manifest(), schema }).schema).toBe(schema);
      expect(parseDataset({ ...(fixtureFile(`/data/datasets/${id}.json`) as object), schema }).schema).toBe(schema);
    }
    for (const schema of [3, 0, "2"]) {
      try {
        parseManifest({ ...manifest(), schema });
        throw new Error("should have thrown");
      } catch (e) {
        expect(e).toBeInstanceOf(ShapeError);
        expect((e as ShapeError).foundSchema).toBe(schema);
        expect((e as ShapeError).message).toBe(`index.json is schema ${schema}; this dashboard reads schemas 1 and 2`);
      }
    }
  });

  it("rejects a manifest whose datasets field is missing", () => {
    const { datasets: _drop, ...rest } = manifest();
    expect(() => parseManifest(rest)).toThrow(/datasets/);
  });

  it("accepts a manifest whose reports share no k, and rejects a k that is not a number", () => {
    expect(parseManifest({ ...manifest(), k: null }).k).toBeNull();
    expect(() => parseManifest({ ...manifest(), k: "31" })).toThrow(/"k"/);
  });

  describe("best learner of a catalog row", () => {
    type Row = Record<string, unknown> & { learners: { learner: string; random: Record<string, unknown> | null; lineage: Record<string, unknown> | null }[] };
    const edited = (edit: (row: Row, m: Record<string, unknown>) => void) => {
      const m = manifest();
      edit((m.datasets as Row[])[0]!, m);
      return () => parseManifest(m);
    };
    const winnerRow = (row: Row) => row.learners.find((l) => l.learner === row.best_learner)!;
    /** No best learner, and so no panel rows, which report.rs writes only beside one. */
    const noWinner = (row: Row, m: Record<string, unknown>): Row => {
      m.panels = (m.panels as { dataset: unknown }[]).filter((p) => p.dataset !== row.id);
      return Object.assign(row, { best_learner: null, best_random: null, best_lineage: null });
    };

    it("accepts a row with no best learner", () => {
      expect(edited(noWinner)().datasets[0]!.best_learner).toBeNull();
    });

    it("rejects a score or a learner without the other", () => {
      expect(edited((row, m) => (noWinner(row, m).best_random = 0.5))).toThrow(/both be null only when no learner has full coverage/);
      expect(edited((row) => (row.best_random = null))).toThrow(/both be null only when no learner has full coverage/);
      expect(edited((row, m) => (noWinner(row, m).best_lineage = 0.5))).toThrow(/cannot have a lineage score/);
      expect(edited((row) => (row.best_random = "0.9"))).toThrow(/best_random should be a number or null/);
      expect(edited((row) => (row.best_learner = 3))).toThrow(/best_learner should be a learner or null/);
    });

    it("rejects a best learner that leaves isolates uncalled or does not score best_random", () => {
      const message = /should call every isolate and score best_random/;
      expect(edited((row) => (winnerRow(row).random!.no_call = { resistant: 1, susceptible: 0 }))).toThrow(message);
      expect(edited((row) => (winnerRow(row).lineage = { ...winnerRow(row).random!, no_call: { resistant: 0, susceptible: 2 } }))).toThrow(message);
      expect(edited((row) => (winnerRow(row).lineage = null))).toThrow(message);
      expect(edited((row) => (row.best_random = (row.best_random as number) - 0.01))).toThrow(message);
      expect(edited((row) => (winnerRow(row).random!.no_call = { resistant: 0, susceptible: 0 }))).not.toThrow();
    });

    it("rejects a panel row that names a learner other than its dataset's best", () => {
      const panelOf = (m: Record<string, unknown>) => (m.panels as Record<string, unknown>[])[0]!;
      const onRow = (row: Row, m: Record<string, unknown>, best: unknown = row.best_learner) => Object.assign(panelOf(m), { dataset: row.id, best_learner: best });
      expect(edited((row, m) => onRow(row, m))).not.toThrow();
      expect(edited((row, m) => onRow(row, m, row.best_learner === "scm" ? "cart" : "scm"))).toThrow(/is not the best learner of/);
      expect(edited((row, m) => (onRow(row, m), Object.assign(row, { best_learner: null, best_random: null, best_lineage: null })))).toThrow(/is not the best learner of/);
      expect(edited((row, m) => onRow(row, m, null))).toThrow(/"best_learner" should be a string/);
    });

    it("reads null rates only for a learner that called no isolate, with no_call counting them", () => {
      const rates = { balanced_accuracy: null, sensitivity: null, specificity: null, vme: null, me: null, auroc: null };
      const other = (row: Row) => row.learners.find((l) => l.learner !== row.best_learner)!;
      expect(edited((row) => (other(row).lineage = { ...rates, no_call: { resistant: 3, susceptible: 4 } }))).not.toThrow();
      expect(edited((row) => (other(row).lineage = rates))).toThrow(/null rates need no_call/);
      expect(edited((row) => (other(row).lineage = { ...rates, balanced_accuracy: 0, no_call: { resistant: 3, susceptible: 4 } }))).toThrow(/all be numbers, or all null/);
    });
  });

  it("rejects a dataset report whose learner scores do not cover its samples", () => {
    const id = parseManifest(manifest()).datasets[0]!.id;
    const raw = fixtureFile(`/data/datasets/${id}.json`) as { splits: { learners: { scores: number[] }[] }[] };
    raw.splits[0]!.learners[0]!.scores.pop();
    expect(() => parseDataset(raw)).toThrow(/scores do not cover/);
  });

  it("checks that a probe panel's called isolates and no-calls add up to the dataset", () => {
    const id = parseManifest(manifest()).datasets[0]!.id;
    const raw = fixtureFile(`/data/datasets/${id}.json`) as Record<string, unknown> & { counts: { resistant: number; susceptible: number } };
    const { resistant: r, susceptible: s } = raw.counts;
    const probe = (noCall: { resistant: number; susceptible: number }) => ({
      id: "gc_gyra",
      name: "gyrA codons 91-95 probes",
      mirrors: "gonorrhoea genotyping",
      targets: [],
      detected: [["mutation detected", r - 2, 0]],
      metrics: { n: r + s - 3, tp: r - 2, fn: 0, tn: s - 1, fp: 0 },
      no_call: noCall,
    });
    expect(() => parseDataset({ ...raw, panels: [probe({ resistant: 2, susceptible: 1 })] })).not.toThrow();
    expect(() => parseDataset({ ...raw, panels: [probe({ resistant: 1, susceptible: 1 })] })).toThrow(/does not add up/);
  });

  it("checks that a learner's null scores and no-call counts add up to the dataset", () => {
    const id = parseManifest(manifest()).datasets[0]!.id;
    type Row = { learner: string; scores: (number | null)[]; metrics: Record<string, number>; no_call?: { resistant: number; susceptible: number } };
    const raw = fixtureFile(`/data/datasets/${id}.json`) as Record<string, unknown> & {
      samples: { resistant: boolean }[];
      splits: { folds: number; test_fold?: number[]; learners: Row[] }[];
    };
    const split = raw.splits[0]!;
    const folds = raw.samples.map((_, i) => i % split.folds);
    // Fold 0's boosting fit was null, so its isolates have no call.
    const gbt = split.learners.find((l) => l.learner === "gbt")!;
    gbt.scores = gbt.scores.map((s, i) => (i % split.folds === 0 ? null : s));
    const calls = raw.samples.flatMap((s, i) => (gbt.scores[i] === null ? [] : [[s.resistant, gbt.scores[i]! >= 0.5]]));
    const count = (truth: boolean, call: boolean) => calls.filter(([t, c]) => t === truth && c === call).length;
    gbt.metrics = { ...gbt.metrics, n: calls.length, tp: count(true, true), fn: count(true, false), tn: count(false, false), fp: count(false, true) };
    const uncalled = (r: boolean) => raw.samples.filter((s, i) => s.resistant === r && gbt.scores[i] === null).length;
    const noCall = { resistant: uncalled(true), susceptible: uncalled(false) };
    const report = (row: Row, testFold: number[] = folds) => ({ ...raw, splits: [{ ...split, test_fold: testFold, learners: [row] }] });

    expect(() => parseDataset(report({ ...gbt, no_call: noCall }))).not.toThrow();
    expect(() => parseDataset(report({ ...gbt, no_call: undefined }))).toThrow(/without a call but no no_call/);
    expect(() => parseDataset(report({ ...gbt, no_call: { ...noCall, resistant: noCall.resistant + 1 } }))).toThrow(/does not add up/);
    const r0 = gbt.scores.findIndex((s, i) => s === null && raw.samples[i]!.resistant);
    const s0 = gbt.scores.findIndex((s, i) => s !== null && !raw.samples[i]!.resistant);
    const moved = gbt.scores.map((s, i) => (i === r0 ? 0.9 : i === s0 ? null : s));
    expect(() => parseDataset(report({ ...gbt, scores: moved, no_call: noCall }))).toThrow(/do not match its no_call/);
    expect(() => parseDataset(report({ ...gbt, no_call: noCall }, folds.slice(1)))).toThrow(/test_fold/);
  });

  it("reads null metrics only for a learner that called no isolate", () => {
    // Revision 5 scored every lineage fold of E. coli amikacin's boosting 0.5
    // and reported 0% VME; with no isolate called the row has no metrics.
    type Row = Record<string, unknown> & { learner: string; scores: (number | null)[]; fold_balanced_accuracy: (number | null)[]; roc: unknown[]; metrics: object };
    const raw = fixtureFile(`/data/datasets/${parseManifest(manifest()).datasets[0]!.id}.json`) as Record<string, unknown> & {
      counts: { resistant: number; susceptible: number };
      splits: { learners: Row[] }[];
    };
    const split = raw.splits[0]!;
    const gbt = split.learners.find((l) => l.learner === "gbt")!;
    const none = {
      ...gbt,
      metrics: null,
      scores: gbt.scores.map(() => null),
      fold_balanced_accuracy: gbt.fold_balanced_accuracy.map(() => null),
      roc: [],
      no_call: { resistant: raw.counts.resistant, susceptible: raw.counts.susceptible },
    };
    const report = (row: Record<string, unknown>) => ({ ...raw, splits: [{ ...split, learners: [row] }] });

    expect(parseDataset(report(none)).splits[0]!.learners[0]!.metrics).toBeNull();
    expect(() => parseDataset(report({ ...none, roc: gbt.roc }))).toThrow(/no ROC points/);
    expect(() => parseDataset(report({ ...none, fold_balanced_accuracy: gbt.fold_balanced_accuracy }))).toThrow(/no fold balanced accuracy/);
    expect(() => parseDataset(report({ ...none, scores: gbt.scores, no_call: { resistant: 0, susceptible: 0 } }))).toThrow(/null metrics but calls some isolates/);
    expect(() => parseDataset(report({ ...gbt, metrics: null }))).toThrow(/null metrics but no no_call/);
    const zero = { ...gbt.metrics, n: 0, tp: 0, fn: 0, tn: 0, fp: 0, balanced_accuracy: 0, vme: 0, me: 0 };
    expect(() => parseDataset(report({ ...none, metrics: zero }))).toThrow(/its metrics should be null/);
  });

  it("reads a report with or without its evaluation settings", () => {
    const raw = fixtureFile(`/data/datasets/${parseManifest(manifest()).datasets[0]!.id}.json`) as Record<string, unknown>;
    const config = {
      folds: 5,
      inner_folds: 3,
      scm_max_rules: 6,
      cart_max_depth: 8,
      forest: { n_trees: 200, subsample: 0.632, max_depth: 30, min_leaf: 1 },
      gbt: { rounds: 400, learning_rate: 0.1, max_depth: 4, lambda: 1, min_child_weight: 1, colsample: 0.5, patience: 25 },
      screened_features: 4000,
      lasso_lambdas: 30,
      seed: 42,
    };
    expect(parseDataset({ ...raw, eval_config: config }).eval_config?.screened_features).toBe(4000);
    expect(parseDataset(raw).eval_config).toBeUndefined();
    expect(() => parseDataset({ ...raw, eval_config: 42 })).toThrow(/eval_config/);
  });

  it("rejects unsupported-assay entries without a reason", () => {
    const m = parseManifest(manifest());
    const raw = fixtureFile(`/data/datasets/${m.datasets[0]!.id}.json`) as Record<string, unknown>;
    const entry = { id: "ecoli_qrdr", name: "gyrA and parC QRDR probes", mirrors: "QRDR", reason: "synonymous changes" };
    expect(() => parseDataset({ ...raw, unsupported_assays: [entry] })).not.toThrow();
    expect(() => parseDataset({ ...raw, unsupported_assays: [{ ...entry, reason: undefined }] })).toThrow(/unsupported_assays\[0\]/);
    const row = { panel: "ecoli_qrdr", name: "gyrA and parC QRDR probes", dataset: "escherichia_coli__ciprofloxacin", species: "Escherichia coli", antibiotic: "ciprofloxacin", reason: "synonymous changes" };
    expect(() => parseManifest({ ...manifest(), unsupported_panels: [row] })).not.toThrow();
    expect(() => parseManifest({ ...manifest(), unsupported_panels: [{ ...row, reason: 3 }] })).toThrow(/unsupported_panels\[0\]/);
  });

  it("reads manifest rows and dataset reports that record no class minimum or left-out MIC labels", () => {
    const strip = ({ min_class: _min, labels_from_mic_excluded: _left, ...rest }: Record<string, unknown>) => rest;
    const m = manifest() as { datasets: Record<string, unknown>[] };
    expect(() => parseManifest({ ...m, datasets: m.datasets.map(strip) })).not.toThrow();
    const raw = fixtureFile(`/data/datasets/${parseManifest(m).datasets[0]!.id}.json`) as { counts: Record<string, unknown> };
    expect(raw.counts.min_class).toBe(25);
    expect(() => parseDataset({ ...raw, counts: strip(raw.counts) })).not.toThrow();
  });

  it("rejects a class minimum or left-out MIC label count that is not a whole number", () => {
    const m = manifest() as { datasets: Record<string, unknown>[] };
    const row = (key: string, value: unknown) => ({ ...m, datasets: [{ ...m.datasets[0]!, [key]: value }, ...m.datasets.slice(1)] });
    expect(() => parseManifest(row("min_class", "25"))).toThrow(/min_class/);
    expect(() => parseManifest(row("labels_from_mic_excluded", -1))).toThrow(/labels_from_mic_excluded/);
    const raw = fixtureFile(`/data/datasets/${parseManifest(m).datasets[0]!.id}.json`) as { counts: Record<string, unknown> };
    expect(() => parseDataset({ ...raw, counts: { ...raw.counts, labels_from_mic_excluded: 2.5 } })).toThrow(/labels_from_mic_excluded/);
  });

  it("rejects a dataset report with a class below its class minimum", () => {
    const raw = fixtureFile(`/data/datasets/${parseManifest(manifest()).datasets[0]!.id}.json`) as { counts: { resistant: number; susceptible: number } };
    const fewest = Math.min(raw.counts.resistant, raw.counts.susceptible);
    expect(() => parseDataset({ ...raw, counts: { ...raw.counts, min_class: fewest } })).not.toThrow();
    expect(() => parseDataset({ ...raw, counts: { ...raw.counts, min_class: fewest + 1 } })).toThrow(/class minimum/);
  });

  it("rejects a species file whose profiles do not match its drug list", () => {
    const slug = parseManifest(manifest()).species[0]!.slug;
    const raw = fixtureFile(`/data/species/${slug}.json`) as { drugs: string[] };
    raw.drugs.push("an extra drug");
    expect(() => parseSpecies(raw)).toThrow(/profile/);
  });
});
