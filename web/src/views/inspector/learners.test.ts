import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../../data/validate";
import type { DatasetReport, LearnerResult } from "../../data/schema";
import { int, LEARNER_NAME, learnerInSentence, pct } from "../../data/format";
import { bestLearner } from "../../data/derive/dataset";
import { LearnerSection } from "./LearnerSection";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const report = parseDataset(fixtureFile(`/data/datasets/${manifest.datasets[0]!.id}.json`));
const render = (r: DatasetReport) => renderToStaticMarkup(createElement(LearnerSection, { r }));
const n = report.samples.length;

const withRows = (edit: (l: LearnerResult, scheme: string) => LearnerResult, base: DatasetReport = report): DatasetReport => ({
  ...base,
  splits: base.splits.map((s) => ({ ...s, learners: s.learners.map((l) => edit(l, s.scheme)) })),
});

/** Cell text after the strip in the table row headed by `name`: random,
 *  lineages held out, VME, ME, the coverage column if any, AUROC. */
const rowCells = (html: string, name: string): string[] => {
  const row = html.split("<tr").find((r) => r.includes('scope="row"') && r.includes(`>${name}<`))!;
  return [...row.matchAll(/<td[^>]*>([\s\S]*?)<\/td>/g)].map((m) => m[1]!.replace(/<[^>]+>/g, "").trim()).slice(1);
};

describe("learner coverage", () => {
  it("shows the share each row calls beside VME and ME once rows count no-calls", () => {
    // The lasso's fold-0 fit was null; every other trained row called all.
    const nulled = report.samples.filter((_, i) => i % 5 === 0);
    const r = withRows((l, scheme) =>
      l.learner === "catalog"
        ? l
        : scheme === "random" && l.learner === "lasso"
          ? {
              ...l,
              scores: l.scores.map((v, i) => (i % 5 === 0 ? null : v)),
              fold_balanced_accuracy: [null, ...l.fold_balanced_accuracy.slice(1)],
              fold_complexity: [0, ...l.fold_complexity.slice(1)],
              no_call: { resistant: nulled.filter((s) => s.resistant).length, susceptible: nulled.filter((s) => !s.resistant).length },
            }
          : { ...l, no_call: { resistant: 0, susceptible: 0 } },
    );
    const html = render(r);
    expect(html).toContain(">Called<");
    expect(html).toContain(pct((n - nulled.length) / n));
    expect(html).toContain("gets no call: every rate covers the isolates called");
    expect(html).not.toContain("Null folds");
  });

  it("counts null folds from fold_complexity in a report written before no-calls", () => {
    const r = withRows((l, scheme) => (scheme === "random" && l.learner === "gbt" ? { ...l, fold_complexity: [0, 0, ...l.fold_complexity.slice(2)] } : l));
    const html = render(r);
    expect(html).toContain("Null folds");
    expect(html).toContain("2 of 5");
    expect(html).toContain("at 0.5, which the 0.5 cutoff calls resistant");
    expect(html).not.toContain(">Called<");
  });

  it("adds no coverage column to an older report without null folds", () => {
    const html = render(report);
    expect(html).not.toContain(">Called<");
    expect(html).not.toContain("Null folds");
  });
});

describe("best learner and learners without calls", () => {
  const withLineage = parseDataset(fixtureFile(`/data/datasets/${manifest.datasets.find((d) => d.best_lineage !== null)!.id}.json`));
  const total = withLineage.samples.length;
  const nulled = withLineage.samples.filter((_, i) => i % 5 === 0);
  const noCall = { resistant: nulled.filter((s) => s.resistant).length, susceptible: nulled.filter((s) => !s.resistant).length };
  const counted = (l: LearnerResult): LearnerResult => (l.learner === "catalog" ? l : { ...l, no_call: { resistant: 0, susceptible: 0 } });
  const partial = (l: LearnerResult): LearnerResult => ({ ...l, scores: l.scores.map((v, i) => (i % 5 === 0 ? null : v)), no_call: noCall });
  // Every fold null, as E. coli amikacin's lineage boosting in revision 5.
  const uncalled = (l: LearnerResult): LearnerResult => ({
    ...l,
    metrics: null,
    scores: l.scores.map(() => null),
    fold_balanced_accuracy: l.fold_balanced_accuracy.map(() => null),
    fold_complexity: l.fold_complexity.map(() => 0),
    roc: [],
    no_call: { resistant: withLineage.counts.resistant, susceptible: withLineage.counts.susceptible },
  });
  const winner = bestLearner(withLineage)!.learner;
  const other = (["scm", "cart", "rf", "gbt", "lasso"] as const).find((l) => l !== winner)!;

  it("explains in the default view why a higher random score is not marked best", () => {
    // `other` calls every isolate on random splits and scores highest there,
    // but leaves the null-fold isolates uncalled with lineages held out.
    const edit = (l: LearnerResult, scheme: string): LearnerResult => {
      const c = counted(l);
      if (l.learner !== other) return c;
      return scheme === "random" ? { ...c, metrics: { ...c.metrics!, balanced_accuracy: 1 } } : partial(c);
    };
    const html = render(withRows(edit, withLineage));
    expect(html).toContain(
      `${LEARNER_NAME[other]} scores 100.0% on random splits, above ${learnerInSentence(winner)}, but gives no call for ${int(nulled.length)} of ${int(total)} isolates with lineages held out. The best learner has to call every isolate on every split it ran.`,
    );
    expect(rowCells(html, LEARNER_NAME[other])[0]).toBe("100.0%");
    expect(render(withRows(counted, withLineage))).not.toContain("The best learner has to call every isolate");
  });

  it("says learners exist but none calls every isolate, not that no result is recorded", () => {
    const html = render(withRows((l) => (l.learner === "catalog" ? l : partial(l)), withLineage));
    expect(html).toContain("No learner calls every isolate on every split it ran, so none is marked best; Called gives the share each one calls.");
    expect(html).not.toContain("No learner result is recorded");
    expect(html).not.toContain(">Best<");
  });

  it("shows No call, not 0% error, for a learner that called no isolate", () => {
    const amikacin = withRows((l, scheme) => (l.learner === "gbt" && scheme === "lineage" ? uncalled(l) : counted(l)), withLineage);
    const random = withLineage.splits[0]!.learners.find((l) => l.learner === "gbt")!.metrics!.balanced_accuracy;
    expect(rowCells(render(amikacin), LEARNER_NAME.gbt).slice(0, 2)).toEqual([pct(random), "No call"]);
    const neither = withRows((l) => (l.learner === "gbt" ? uncalled(l) : counted(l)), withLineage);
    const html = render(neither);
    const cells = rowCells(html, LEARNER_NAME.gbt);
    // The Called cell shows its share and count in view: none of the isolates.
    expect(cells.slice(0, 4)).toEqual(["No call", "No call", "No call", "No call"]);
    expect(cells[4]).toBe(`0.0%0 of ${int(total)}`);
    expect(cells[5]).toBe("–");
    expect(html).toContain("A learner that calls no isolate on a split has no rate there and shows No call.");
  });
});

describe("no-call details in view", () => {
  it("names each row's null folds and uncalled isolates under the table, with no hover", () => {
    const nulled = report.samples.filter((_, i) => i % 5 === 0);
    const resistant = nulled.filter((s) => s.resistant).length;
    const r = withRows((l, scheme) =>
      l.learner === "catalog"
        ? l
        : scheme === "random" && l.learner === "lasso"
          ? { ...l, scores: l.scores.map((v, i) => (i % 5 === 0 ? null : v)), fold_complexity: [0, ...l.fold_complexity.slice(1)], no_call: { resistant, susceptible: nulled.length - resistant } }
          : { ...l, no_call: { resistant: 0, susceptible: 0 } },
    );
    const html = render(r);
    expect(html).toContain(
      `${LEARNER_NAME.lasso} on random splits: in 1 of 5 folds the lasso kept no weight, so ${int(resistant)} resistant and ${int(nulled.length - resistant)} susceptible isolates get no call; ${int(n - nulled.length)} of ${int(n)} are called.`,
    );
    expect(rowCells(html, LEARNER_NAME.lasso)[4]).toBe(`${pct((n - nulled.length) / n)}${int(n - nulled.length)} of ${int(n)}`);
    expect(html).not.toContain('role="tooltip"');
  });

  it("states a pre-revision-7 report's null folds as text in the cell", () => {
    const r = withRows((l, scheme) => (scheme === "random" && l.learner === "gbt" ? { ...l, fold_complexity: [0, 0, ...l.fold_complexity.slice(2)] } : l));
    expect(rowCells(render(r), LEARNER_NAME.gbt)[4]).toBe("2 of 5");
  });
});
