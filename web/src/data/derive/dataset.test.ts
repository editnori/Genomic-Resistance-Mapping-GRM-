import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../validate";
import type { DatasetReport, FeatureCard, Learner, LearnerResult, Scheme } from "../schema";
import {
  answerSentence,
  bestLearner,
  catalogRuleWords,
  contextRuns,
  familiesNamed,
  familyLabel,
  familyOfGroup,
  familyRuleWords,
  featureLabel,
  featurePlace,
  featureWhere,
  indexRuleWords,
  NO_BEST_LEARNER,
  RULE_TEXT_VERSION,
  ruleWordingNote,
  oofErrors,
  outscoringPartial,
  panelCalls,
  ruleWords,
  SEVERAL_FAMILIES_NOTE,
  termLabel,
  tieWords,
  unexplained,
} from "./dataset";
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

  it("names every catalog family a feature matches, in one order, and never just the first", () => {
    const shared = "k-mers matching 3 catalog families (blaAFM, blaNDM, ble)";
    for (const families of [["blaAFM", "blaNDM", "ble"], ["ble", "blaNDM", "blaAFM"], ["blaNDM", "blaAFM", "ble", "blaAFM"]]) {
      expect(familyLabel(families)).toBe(shared);
      expect(featureLabel(card({ gene: "blaNDM-1" }, families), 12)).toBe(shared);
    }
    expect(familyLabel(["blaKPC", "blaKPC"])).toBe("blaKPC");
    expect(familyLabel([])).toBeNull();
    const features = new Map([[12, card({ gene: "blaNDM-1" }, ["blaNDM", "ble", "blaAFM"])]]);
    expect(ruleWords({ feature: 12, presence: true }, features)).toBe(`Presence of ${shared}`);
    expect(ruleWords({ feature: 12, presence: false }, features)).toBe(`Absence of ${shared}`);
  });

  it("labels a schema-1 report's cards the same way, from the families it already holds", () => {
    const raw = structuredClone(fixtureFile(`/data/datasets/${withPanel.id}.json`)) as { schema: number; features: { catalog_families: string[] }[] };
    raw.schema = 1;
    raw.features[0]!.catalog_families = ["blaNDM", "ble", "blaAFM"];
    const old = parseDataset(raw);
    expect(featureLabel(old.features[0], old.features[0]!.feature)).toBe("k-mers matching 3 catalog families (blaAFM, blaNDM, ble)");
  });

  it("places a feature by codon range or by its nearest gene", () => {
    expect(featurePlace(card({ gene: "gyrA", codons: [83, 87] }))).toBe("gyrA codons 83–87");
    expect(featurePlace(card({ gene: "gyrA", codons: [83, 83] }))).toBe("gyrA codon 83");
    expect(featurePlace(card({ nearest: ["mgrB", 46] }))).toBe("intergenic, 46 bp from mgrB");
    expect(featurePlace(card(null))).toBeNull();
  });

  it("says how far a rule term's k-mers sit from the gene it is named after when they sit outside it", () => {
    // Synthetic cards. Inside a gene the term keeps the gene's name.
    expect(termLabel(card({ gene: "mgrB", codons: [3, 12] }), 12)).toBe("mgrB");
    // Outside every gene the locus names the nearest one, and the term says so.
    const outside = card({ gene: "mgrB", product: "PhoP/PhoQ regulator MgrB", nearest: ["mgrB", 46] });
    expect(termLabel(outside, 12)).toBe("k-mers 46 bp from the nearest annotated gene, mgrB");
    expect(termLabel(card({ product: "hypothetical protein", nearest: ["hypothetical protein", -2502] }), 12)).toBe("k-mers 2,502 bp from the nearest annotated gene, hypothetical protein");
    const features = new Map([[12, outside]]);
    expect(ruleWords({ feature: 12, presence: true }, features)).toBe("Presence of k-mers 46 bp from the nearest annotated gene, mgrB");
    expect(ruleWords({ feature: 12, presence: false }, features)).toBe("Absence of k-mers 46 bp from the nearest annotated gene, mgrB");
    // A catalog family is named wherever the k-mers sit; no name leaves the number.
    expect(termLabel(card({ gene: "tnpA", nearest: ["tnpA", 120] }, ["blaKPC"]), 12)).toBe("blaKPC");
    expect(termLabel(card({ nearest: ["", 9] }), 12)).toBe("k-mer feature 12");
    expect(termLabel(card(null), 12, "feature 12")).toBe("feature 12");
    expect(termLabel(undefined, 12)).toBe("k-mer feature 12");
    // A table that prints the place in its own column keeps the short name.
    expect(featureLabel(outside, 12)).toBe("mgrB");
    expect(featureWhere(outside, 12)).toBe("intergenic, 46 bp from mgrB");
  });

  it("counts a tied rule term among the features that separate the isolates equally well", () => {
    // n_equivalent counts the other rules, so 3 of them make 4 features.
    expect(tieWords(3, 0)).toBe("one of 4 features that separate these isolates equally well");
    expect(tieWords(1, 1)).toBe("one of 2 features that separate the isolates the earlier rules left undecided equally well");
    expect(tieWords(21523, 2)).toBe("one of 21,524 features that separate the isolates the earlier rules left undecided equally well");
    expect(tieWords(0, 0)).toBeNull();
    expect(tieWords(0, 3)).toBeNull();
  });

  it("reads a gene group's family back from its name", () => {
    expect(familyOfGroup("blaOXA (OXA-48 family)")).toBe("blaOXA");
    expect(familyOfGroup("blaKPC")).toBe("blaKPC");
    expect(familyOfGroup("aac(6')-Ib")).toBe("aac(6')-Ib");
    expect(familyOfGroup("tet(A)")).toBe("tet(A)");
    expect(familyOfGroup("x (two words family)")).toBe("x (two words family)");
  });

  it("says which named gene groups count for the drug's rule, and nothing where the report has no row", () => {
    // Synthetic gene-group rows: the groups called in the dataset's genomes.
    const group = (name: string, counts: boolean) => ({ group: name, class: "", subclass: "", n_resistant: 1, n_susceptible: 1, counts_for_rule: counts });
    const groups = [group("ble", false), group("blaKPC", true), group("blaNDM", true), group("blaOXA (OXA-48 family)", true), group("blaOXA (OXA-1 family)", false), group("qacE", false), group("sul1", false)];
    const catalog = { ...report.catalog!, groups };
    const rule = "for this drug's known-gene rule";
    expect(familyRuleWords(["blaKPC"], catalog)).toBe(`blaKPC counts ${rule}`);
    expect(familyRuleWords(["sul1"], catalog)).toBe(`sul1 does not count ${rule}`);
    expect(familyRuleWords(["qacE", "sul1"], catalog)).toBe(`qacE and sul1 do not count ${rule}`);
    expect(familyRuleWords(["blaOXA"], catalog)).toBe(`blaOXA (OXA-48 family) counts ${rule}, but blaOXA (OXA-1 family) does not`);
    expect(familyRuleWords(["sul1", "qacE", "blaKPC", "blaNDM"], catalog)).toBe(`blaKPC and blaNDM count ${rule}, but qacE and sul1 do not`);
    // blaAFM has no row: nothing is said of it, and the rest still stands.
    expect(familyRuleWords(["blaNDM", "ble", "blaAFM"], catalog)).toBe(`blaNDM counts ${rule}, but ble does not`);
    // No row for any named family, no families, or no catalog block: nothing.
    expect(familyRuleWords(["tet(A)"], catalog)).toBeNull();
    expect(familyRuleWords(["blaAFM", "tet(A)"], catalog)).toBeNull();
    expect(familyRuleWords([], catalog)).toBeNull();
    expect(familyRuleWords(["blaKPC"], null)).toBeNull();
    expect(familyRuleWords(["blaKPC"], { ...catalog, groups: [] })).toBeNull();
  });

  it("writes each rule term as the index does, on the cases the index writer is held to", async () => {
    // tests/data/rule-text-cases.json: synthetic reports and the words
    // report::rule_text must give them (src/report.rs reads the same file).
    const fs = (await import(/* @vite-ignore */ `node:${"fs"}`)) as { readFileSync(path: string, encoding: string): string };
    const cwd = (globalThis as unknown as { process: { cwd(): string } }).process.cwd();
    const file = JSON.parse(fs.readFileSync(`${cwd}/../tests/data/rule-text-cases.json`, "utf8")) as { rule_text_version: number; cases: Record<string, { report: DatasetReport; scm_rules: string[] }> };
    expect(file.rule_text_version).toBe(RULE_TEXT_VERSION);
    const names = Object.keys(file.cases);
    expect(names).toEqual(["outside_the_named_gene", "ties", "counts_for_the_drugs_rule", "no_catalog_block"]);
    for (const name of names) {
      const { report: r, scm_rules: words } = file.cases[name]!;
      expect([name, r.models.scm.rules.map((_, i) => indexRuleWords(r, i))]).toEqual([name, words]);
    }
  });

  it("writes the synthetic fixture's index in the same words, and numbers the wording", () => {
    expect(manifest.rule_text_version).toBe(RULE_TEXT_VERSION);
    let terms = 0;
    for (const row of manifest.datasets) {
      const r = parseDataset(fixtureFile(`/data/datasets/${row.id}.json`));
      expect([row.id, row.scm_rules]).toEqual([row.id, r.models.scm.rules.map((_, i) => indexRuleWords(r, i))]);
      terms += row.scm_rules.length;
    }
    expect(terms).toBeGreaterThan(0);
  });

  it("names the release's revision and rule wording beside a warning on counting rules by what they name", () => {
    const tail = " A later release can word the same rule differently, so a count of rules by what they name holds for one wording only.";
    expect(ruleWordingNote({ evaluation_revision: 14, rule_text_version: 3 })).toBe(`The rule words in the datasets table come from this release's index (evaluation revision 14, rule wording 3).${tail}`);
    // The served release's index has no wording number: only its revision is named.
    expect(ruleWordingNote({ evaluation_revision: 14 })).toBe(`The rule words in the datasets table come from this release's index (evaluation revision 14).${tail}`);
    expect(ruleWordingNote({})).toBe(`The rule words in the datasets table come from this release's index.${tail}`);
  });

  it("sees a term that names several catalog families in the index's words, and no other", () => {
    expect(familiesNamed("presence of k-mers matching 3 catalog families (blaAFM, blaNDM, ble)")).toBe(3);
    expect(familiesNamed("absence of k-mers matching 12 catalog families (a, b); blaNDM counts for this drug's known-gene rule, but ble does not")).toBe(12);
    expect(familiesNamed("presence of blaKPC")).toBeNull();
    expect(familiesNamed("presence of blaKPC; blaKPC counts for this drug's known-gene rule")).toBeNull();
    expect(familiesNamed("presence of k-mers 46 bp from the nearest annotated gene, mgrB")).toBeNull();
    expect(familiesNamed("presence of feature 4")).toBeNull();
    expect(SEVERAL_FAMILIES_NOTE).toBe("Not a list of genes an isolate carries: one k-mer pattern matches all of these families.");
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

  it("says when other features tie the rule it names, and nothing when none does", () => {
    // The fixture report with its first step's ties set by hand.
    const tied = (rules: number, nEquivalent: number) => {
      const scm = report.models.scm;
      const first = { ...scm.steps[0]!, n_equivalent: nEquivalent };
      return answerSentence({ ...report, models: { ...report.models, scm: { ...scm, rules: Array.from({ length: rules }, () => scm.rules[0]!), steps: [first, ...scm.steps.slice(1)] } } });
    };
    expect(tied(1, 3).endsWith(". The rule is one of 4 features that separate these isolates equally well.")).toBe(true);
    expect(tied(5, 3).endsWith(". The first rule is one of 4 features that separate these isolates equally well.")).toBe(true);
    expect(tied(2, 21523).endsWith(". The first rule is one of 21,524 features that separate these isolates equally well.")).toBe(true);
    expect(tied(1, 0)).not.toContain("equally well");
    expect(tied(5, 0)).not.toContain("equally well");
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

  it("lets a revision-14 learner with null-fit folds be best, since its rates cover every isolate", () => {
    // Revision 14 scores a null fit's isolates at their fold's training
    // prevalence and keeps the older figure in called_metrics; the same
    // no_call counts on an older row leave it out.
    const scored = (l: LearnerResult): LearnerResult => ({ ...l, no_call: noCall, called_metrics: { ...l.metrics!, n: n - nulled.length } });
    const r14 = { ...withRows(other, (l, scheme) => (scheme === "random" ? topRandom(l) : scored(l))), evaluation: { revision: 14, valid: true, reasons: [] } };
    expect(bestLearner(r14)?.learner).toBe(other);
    expect(outscoringPartial(r14)).toEqual([]);
    const r13 = withRows(other, (l, scheme) => (scheme === "random" ? topRandom(l) : { ...l, no_call: noCall }));
    expect(bestLearner(r13)?.learner).toBe(winner);
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
