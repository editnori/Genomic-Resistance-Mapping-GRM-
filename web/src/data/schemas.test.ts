import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { ShapeError, parseDataset, parseManifest } from "./validate";
import { loadErrorText } from "./load";
import { bestLearner } from "./derive/dataset";
import { LearnerSection } from "../views/inspector/LearnerSection";

// A report directory rebuilt after the schema bump holds a schema-2 index
// beside schema-1 reports of datasets not yet recomputed and schema-2 reports
// of datasets that were.

type Raw = Record<string, unknown> & {
  counts: { resistant: number; susceptible: number };
  splits: { scheme: string; learners: (Record<string, unknown> & { learner: string; scores: (number | null)[]; fold_balanced_accuracy: (number | null)[] })[] }[];
};

const manifest = fixtureFile("/data/index.json") as Record<string, unknown>;
const id = parseManifest(manifest).datasets.find((d) => d.best_lineage !== null)!.id;
const raw = () => fixtureFile(`/data/datasets/${id}.json`) as Raw;
const lineageGbt = (edit: (l: Raw["splits"][number]["learners"][number], r: Raw) => Record<string, unknown>, other: (l: Record<string, unknown>) => Record<string, unknown>, schema: number) => {
  const r = raw();
  return {
    ...r,
    schema,
    splits: r.splits.map((s) => ({ ...s, learners: s.learners.map((l) => (s.scheme === "lineage" && l.learner === "gbt" ? edit(l, r) : other(l))) })),
  };
};

/** As revision 5 wrote E. coli amikacin: schema 1, no no-call counts, and
 *  lineage boosting kept no tree in any fold, so every isolate scored 0.5. */
const revision5 = () =>
  lineageGbt(
    (l) => ({ ...l, scores: l.scores.map(() => 0.5), fold_complexity: [0, 0, 0, 0, 0], fold_balanced_accuracy: [0.5, 0.5, 0.5, 0.5, 0.5] }),
    (l) => l,
    1,
  );

/** The same learner in schema 2: no isolate called, so no metrics. */
const revision7 = () =>
  lineageGbt(
    (l, r) => ({
      ...l,
      metrics: null,
      scores: l.scores.map(() => null),
      fold_complexity: [0, 0, 0, 0, 0],
      fold_balanced_accuracy: l.fold_balanced_accuracy.map(() => null),
      roc: [],
      no_call: { resistant: r.counts.resistant, susceptible: r.counts.susceptible },
    }),
    (l) => (l.learner === "catalog" ? l : { ...l, no_call: { resistant: 0, susceptible: 0 } }),
    2,
  );

describe("schema compatibility", () => {
  it("reads a schema-2 index beside a schema-1 and a schema-2 report", () => {
    expect(parseManifest({ ...manifest, schema: 2 }).schema).toBe(2);
    const old = parseDataset(revision5());
    const fresh = parseDataset(revision7());
    expect([old.schema, fresh.schema]).toEqual([1, 2]);

    // The schema-1 report keeps its 0.5 scores, explained as null folds.
    const before = renderToStaticMarkup(createElement(LearnerSection, { r: old }));
    expect(before).toContain("Null folds");
    expect(before).toContain("at 0.5, which the 0.5 cutoff calls resistant");
    expect(before).not.toContain("No call");
    expect(bestLearner(old)).toBeDefined();

    // The schema-2 report shows no rate where nothing was called.
    const after = renderToStaticMarkup(createElement(LearnerSection, { r: fresh }));
    expect(after).toContain(">Called<");
    expect(after).toContain("A learner that calls no isolate on a split has no rate there and shows No call.");
    expect(after).not.toContain("Null folds");
    expect(bestLearner(fresh)?.learner).not.toBe("gbt");
  });

  it("stops at a schema it does not read, with the schema message", () => {
    for (const parse of [() => parseManifest({ ...manifest, schema: 3 }), () => parseDataset({ ...raw(), schema: 3 })]) {
      expect(parse).toThrow(ShapeError);
    }
    let text = "";
    try {
      parseManifest({ ...manifest, schema: 3 });
    } catch (e) {
      text = loadErrorText({ kind: "schema", path: "/data/index.json", message: (e as ShapeError).message, found: (e as ShapeError).foundSchema });
    }
    expect(text).toBe("index.json is schema 3; this dashboard reads schemas 1 and 2. Rebuild the report with the matching `grm` version.");
  });
});
