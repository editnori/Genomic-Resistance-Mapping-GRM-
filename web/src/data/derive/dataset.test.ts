import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../validate";
import type { DatasetReport, FeatureCard, Learner, LearnerResult, Scheme } from "../schema";
import { answerSentence, bestLearner, catalogRuleWords, contextRuns, featureLabel, featurePlace, featureWhere, NO_BEST_LEARNER, oofErrors, outscoringPartial, panelCalls, unexplained } from "./dataset";
import { encodeBits } from "../bits";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const withPanel = manifest.datasets.find((d) => d.panels.length > 0)!;
const report = parseDataset(fixtureFile(`/data/datasets/${withPanel.id}.json`));

describe("panel calls", () => {
  it("reproduce the pipeline's panel metrics from sample catalog groups", () => {
    report.panels.forEach((panel) => {
      const calls = panelCalls(report, panel)!;
      const tp = calls.filter((c, i) => c && report.samples[i]!.resistant).length;
      const fp = calls.filter((c, i) => c && !report.samples[i]!.resistant).length;
      expect([tp, fp]).toEqual([panel.metrics.tp, panel.metrics.fp]);
    });
  });

  it("count every missed resistant isolate as unexplained", () => {
    const panel = report.panels[0]!;
    const u = unexplained(report, panel);
    expect(u.perSample).toBe(true);
    expect(u.missed).toBe(panel.metrics.fn);
    u.features.forEach((f) => expect(f.inMissed).toBeLessThanOrEqual(u.missed));
  });

  it("do not rebuild a probe assay, whose detected label is no catalog group", () => {
    // src/analysis.rs writes probe assays with detected = [["mutation detected", r, s]].
    const base = report.panels[0]!;
    const probe = { ...base, id: "xpert_rif", detected: [["mutation detected", base.metrics.tp, base.metrics.fp]] as [string, number, number][] };
    expect(panelCalls(report, probe)).toBeNull();
    const u = unexplained(report, probe);
    expect([u.perSample, u.missed, u.features.length]).toEqual([false, base.metrics.fn, 0]);
  });
});

describe("words", () => {
  it("states an AND of OR groups as one sentence", () => {
    expect(catalogRuleWords([["C:TRIMETHOPRIM"], ["C:SULFONAMIDE"]])).toBe(
      "Resistant when the genome carries a gene of class trimethoprim, and also a gene of class sulfonamide.",
    );
    expect(catalogRuleWords([["S:CARBAPENEM", "S:CEPHALOSPORIN"]])).toBe(
      "Resistant when the genome carries a gene of subclass carbapenem or a gene of subclass cephalosporin.",
    );
  });

  const card = (locus: Partial<NonNullable<FeatureCard["locus"]>> | null, families: string[] = []): FeatureCard => ({
    feature: 12,
    n_kmers: 31,
    kmer: "",
    catalog_families: families,
    catalog_products: [],
    locus:
      locus === null
        ? null
        : { genome_id: "g", contig: "c", start: 1, end: 31, kmers_placed: 1, gene: null, product: null, patric_id: null, gene_strand: null, codons: null, nearest: null, context: "", ...locus },
    n_resistant: 0,
    n_susceptible: 0,
    bits: "",
  });

  it("labels a feature by catalog family, then gene, then product, then id", () => {
    expect(featureLabel(card({ gene: "gyrA" }, ["blaKPC"]), 12)).toBe("blaKPC");
    expect(featureLabel(card({ gene: "gyrA" }), 12)).toBe("gyrA");
    expect(featureLabel(card({ product: "porin" }), 12)).toBe("porin");
    expect(featureLabel(card(null), 12)).toBe("k-mer feature 12");
  });

  it("places a feature by codon range or by its nearest gene", () => {
    expect(featurePlace(card({ gene: "gyrA", codons: [83, 87] }))).toBe("gyrA codons 83–87");
    expect(featurePlace(card({ gene: "gyrA", codons: [83, 83] }))).toBe("gyrA codon 83");
    expect(featurePlace(card({ nearest: ["mgrB", 46] }))).toBe("intergenic, 46 bp from mgrB");
    expect(featurePlace(card(null))).toBeNull();
  });
});

describe("unexplained resistance", () => {
  it("drops a feature every isolate carries, however many missed isolates have it", () => {
    const panel = report.panels[0]!;
    const everywhere: FeatureCard = { ...report.features[0]!, feature: 999_999, bits: encodeBits(report.samples.map(() => true)) };
    const withHousekeeping = {
      ...report,
      features: [...report.features, everywhere],
      models: { ...report.models, rf_importance: [{ feature: 999_999, weight: 1 }, ...report.models.rf_importance] },
    };
    expect(unexplained(withHousekeeping, panel, 50).features.some((f) => f.feature === 999_999)).toBe(false);
  });
});

describe("context runs", () => {
  it("keeps short spans whole and elides the middle of a long one", () => {
    expect(contextRuns("acgtACGTacgt").map((r) => [r.text, r.span])).toEqual([
      ["acgt", false],
      ["ACGT", true],
      ["acgt", false],
    ]);
    const runs = contextRuns(`aa${"C".repeat(300)}tt`, 120, 48);
    expect(runs.map((r) => r.elided ?? r.text.length)).toEqual([2, 48, 204, 48, 2]);
  });
});

describe("answer sentence", () => {
  const withCatalog = (ba: number) => ({
    ...report,
    splits: report.splits.map((s) => ({
      ...s,
      learners: s.learners.map((l) => (l.learner === "catalog" ? { ...l, metrics: { ...l.metrics!, balanced_accuracy: ba } } : l)),
    })),
  });

  it("names the catalog only when it is within the cutoff of the best learner", () => {
    expect(answerSentence(withCatalog(1))).toContain("the known-gene catalog alone reaches 100.0%");
    expect(answerSentence(withCatalog(0.1))).not.toContain("known-gene catalog");
  });

  it("counts the rules past two instead of listing them", () => {
    const rules = Array.from({ length: 5 }, () => report.models.scm.rules[0]!);
    const text = answerSentence({ ...report, models: { ...report.models, scm: { ...report.models.scm, rules } } });
    expect(text.startsWith("A 5-rule set covering machine led by presence of")).toBe(true);
  });
});

describe("no-calls", () => {
  const withLineage = parseDataset(fixtureFile(`/data/datasets/${manifest.datasets.find((d) => d.best_lineage !== null)!.id}.json`));
  const n = withLineage.samples.length;
  // Every isolate whose index is a multiple of 5 sits in a null fold.
  const nulled = withLineage.samples.filter((_, i) => i % 5 === 0);
  const noCall = { resistant: nulled.filter((s) => s.resistant).length, susceptible: nulled.filter((s) => !s.resistant).length };
  const partial = (l: LearnerResult): LearnerResult => ({ ...l, scores: l.scores.map((v, i) => (i % 5 === 0 ? null : v)), no_call: noCall });
  // Every lineage fold null, as E. coli amikacin's boosting in revision 5.
  const uncalled = (l: LearnerResult): LearnerResult => ({
    ...l,
    metrics: null,
    scores: l.scores.map(() => null),
    fold_balanced_accuracy: l.fold_balanced_accuracy.map(() => null),
    roc: [],
    no_call: { resistant: withLineage.counts.resistant, susceptible: withLineage.counts.susceptible },
  });
  const topRandom = (l: LearnerResult): LearnerResult => ({ ...l, metrics: { ...l.metrics!, balanced_accuracy: 1 } });
  /** Every trained row counts no-calls; `edit` changes rows of one learner. */
  const withRows = (learner: Learner, edit: (l: LearnerResult, scheme: Scheme) => LearnerResult): DatasetReport => ({
    ...withLineage,
    splits: withLineage.splits.map((s) => ({
      ...s,
      learners: s.learners.map((l) => {
        if (l.learner === "catalog") return l;
        const counted = { ...l, no_call: { resistant: 0, susceptible: 0 } };
        return l.learner === learner ? edit(counted, s.scheme) : counted;
      }),
    })),
  });
  const winner = bestLearner(withLineage)!.learner;
  const other = (["scm", "cart", "rf", "gbt", "lasso"] as const).find((l) => l !== winner)!;

  it("keeps the top random score from a learner with random-split no-calls", () => {
    const r = withRows(other, (l, scheme) => (scheme === "random" ? partial(topRandom(l)) : l));
    expect(bestLearner(r)?.learner).toBe(winner);
  });

  it("keeps the top random score from a learner that calls every isolate there but not with lineages held out", () => {
    // The learner calls every isolate on random splits and scores highest there.
    const lineagePartial = withRows(other, (l, scheme) => (scheme === "random" ? topRandom(l) : partial(l)));
    const lineageNone = withRows(other, (l, scheme) => (scheme === "random" ? topRandom(l) : uncalled(l)));
    expect(bestLearner(lineagePartial)?.learner).toBe(winner);
    expect(bestLearner(lineageNone)?.learner).toBe(winner);
    expect(outscoringPartial(lineagePartial)).toEqual([{ learner: other, ba: 1, uncalled: [{ scheme: "lineage", n: noCall.resistant + noCall.susceptible }] }]);
    expect(outscoringPartial(withLineage)).toEqual([]);
  });

  it("does not treat a missing lineage learner row as full coverage", () => {
    const r = withRows(other, (l) => topRandom(l));
    const missing = {
      ...r,
      splits: r.splits.map((s) => (s.scheme === "lineage" ? { ...s, learners: s.learners.filter((l) => l.learner !== other) } : s)),
    };
    expect(bestLearner(missing)?.learner).toBe(winner);
  });

  it("ranks on random splits alone when the lineage split did not run", () => {
    const r = withRows(other, (l) => topRandom(l));
    const randomOnly = { ...r, splits: r.splits.map((s) => (s.scheme === "lineage" ? { ...s, learners: [], unavailable_reason: "fewer than twice as many lineages as requested folds" } : s)) };
    expect(bestLearner(randomOnly)?.learner).toBe(other);
  });

  it("has no best learner when every learner leaves isolates uncalled", () => {
    const none = {
      ...withLineage,
      splits: withLineage.splits.map((s) => ({
        ...s,
        learners: s.learners.map((l) => (l.learner === "catalog" ? l : partial(l))),
      })),
    };
    expect(bestLearner(none)).toBeUndefined();
    // The sentence opens with the reason, or puts it after the SCM rule.
    expect(answerSentence(none).toLowerCase()).toContain(`${NO_BEST_LEARNER.toLowerCase()} on every split it ran.`);
    const untrained = { ...withLineage, splits: withLineage.splits.map((s) => ({ ...s, learners: s.learners.filter((l) => l.learner === "catalog") })) };
    expect(answerSentence(untrained).toLowerCase()).toContain("no learner result is recorded for this dataset.");
  });

  it("counts an isolate without a call as no error", () => {
    const r = withRows("lasso", (l, scheme) => (scheme === "random" ? partial(l) : l));
    const lasso = r.splits[0]!.learners.find((l) => l.learner === "lasso")!;
    const errors = oofErrors(r, lasso);
    expect(errors.length).toBe(n);
    expect(errors.some((e, i) => e && i % 5 === 0)).toBe(false);
    expect(errors.some((e, i) => e && i % 5 !== 0)).toBe(true);
  });
});

describe("featureWhere", () => {
  it("drops the label the place would repeat, and keeps a place that names another gene", () => {
    const base = report.features[0]!;
    const card = (locus: Partial<NonNullable<FeatureCard["locus"]>>, families: string[] = []): FeatureCard => ({
      ...base,
      catalog_families: families,
      locus: { genome_id: "g", contig: "c", start: 1, end: 31, kmers_placed: 1, gene: null, product: null, patric_id: null, gene_strand: null, codons: null, nearest: null, context: "", ...locus },
    });
    expect(featureWhere(card({ gene: "gyrA", codons: [83, 87] }), 1)).toBe("Codons 83–87");
    expect(featureWhere(card({ gene: "ompK36", codons: [1, 9] }, ["blaKPC"]), 1)).toBe("ompK36 codons 1–9");
  });
});
