// @vitest-environment happy-dom
// The report's copy and its links. The exact sentences come from the
// annotated brief, checked against the revision-7 publication's index.json
// (sha256 eb7064b6…, generated 2026-10-25T13:55:41Z): its summaries and its
// species table are reproduced here, so each lead is tested against the
// numbers it must print rather than against its own template.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { CatalogRow, Manifest, SpeciesRow } from "../data/schema";
import type { Navigator, Route } from "../app/router";
import { catalogWords } from "../app/ReportStatus";
import type { Fitted, LearnerRow, LearnerSummary } from "../data/derive/evidence";
import type { SplitSummary } from "../data/derive/overview";
import { barsLead, catalogLead, cohortLead, coverageNote, learnersLead, lineageLead, micLead } from "./findings";
import { OverviewView } from "./OverviewView";
import { MethodsView } from "./MethodsView";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const base = parseManifest(fixtureFile("/data/index.json"));

/** Revision 7's species table: [species, curated, reported, any label, analyzed]. */
const REV7_SPECIES: [string, number, number, number, number][] = [
  ["Mycobacterium tuberculosis", 22, 22, 32060, 2000],
  ["Streptococcus pneumoniae", 12, 12, 17837, 2000],
  ["Salmonella enterica", 25, 25, 14522, 2000],
  ["Escherichia coli", 37, 36, 12415, 2000],
  ["Neisseria gonorrhoeae", 7, 7, 10223, 2000],
  ["Klebsiella pneumoniae", 39, 39, 8326, 2000],
  ["Staphylococcus aureus", 19, 18, 5374, 2000],
  ["Enterococcus faecium", 14, 14, 3138, 2000],
  ["Acinetobacter baumannii", 18, 18, 2485, 1705],
  ["Shigella sonnei", 20, 20, 1994, 1984],
  ["Pseudomonas aeruginosa", 13, 13, 1462, 1421],
  ["Enterobacter cloacae", 9, 9, 837, 271],
  ["Campylobacter jejuni", 4, 4, 746, 508],
  ["Campylobacter coli", 8, 0, 407, 0],
  ["Shigella flexneri", 4, 0, 366, 0],
  ["Klebsiella aerogenes", 3, 0, 317, 0],
  ["Corynebacterium diphtheriae", 3, 0, 243, 0],
  ["Clostridioides difficile", 4, 0, 204, 0],
  ["Klebsiella quasipneumoniae", 3, 0, 179, 0],
  ["Streptococcus agalactiae", 3, 0, 157, 0],
  ["Staphylococcus haemolyticus", 3, 0, 128, 0],
  ["Klebsiella oxytoca", 4, 0, 116, 0],
  ["Vibrio cholerae", 3, 0, 91, 0],
];

const slug = (s: string) => s.toLowerCase().replace(/ /g, "_");

/** A manifest with revision 7's species table and one report row per
 *  reported dataset; the rows reuse a fixture row, as only their species
 *  and count matter here. */
const rev7Shape = (): Manifest => {
  const species: SpeciesRow[] = REV7_SPECIES.map(([name, curated, reported, labeled, analyzed]) => ({
    species: name,
    slug: slug(name),
    datasets_curated: curated,
    datasets_reported: reported,
    genomes_labeled: labeled,
    genomes_analyzed: analyzed,
    lineages: analyzed === 0 ? null : 10,
    lineage_threshold: analyzed === 0 ? null : 0.005,
  }));
  const datasets: CatalogRow[] = species.flatMap((s) =>
    Array.from({ length: s.datasets_reported }, (_, i) => ({ ...base.datasets[0]!, id: `${s.slug}__drug${i}`, species: s.species, species_slug: s.slug, min_class: 25 })),
  );
  return { ...base, species, datasets, funnel: { ...base.funnel, genomes_passing_qc: 119_072 } };
};

describe("report leads with revision 7's numbers", () => {
  it("leads with the report cohort and one coverage note", () => {
    const m = rev7Shape();
    expect(cohortLead(m)).toBe("237 datasets have a report. They cover 13 species and use 21,889 genomes, chosen from 119,072 genomes that passed quality control (QC).");
    expect(coverageNote(m)).toBe(
      "In this run, selection is capped at 2,000 genomes per species when more are eligible. Some datasets have fewer because genomes failed QC or had no FASTA file, the file that holds the DNA sequence. A report needs 25 resistant and 25 susceptible selected genomes with a reported call and no conflicting call. Two datasets in these species passed curation but have no report. Another 10 species (38 datasets that passed curation) were not analyzed. Their FASTA files are stored apart from the data this run read.",
    );
  });

  it("states the lineage, known-gene, and MIC leads with their denominators", () => {
    expect(lineageLead({ n: 237, nNoBest: 0, nLineage: 231, medianRandom: 0.9253561253561253, medianLineage: 0.859608745684695, medianGap: 0.045254403131115395, nLargeGap: 106 })).toBe(
      "A learner can score well by recognising a genome's close relatives, so each dataset is also tested on whole lineages it never saw. Holding them out lowers median balanced accuracy from 92.5% on a random split to 86.0%. The random-split figure is each dataset's best of five learners, picked on the same folds that score it, so it runs a little above what one learner named in advance would score. 106 of 231 datasets lose more than 5 points; the median loss is 4.5 points.",
    );
    expect(
      catalogLead({
        nWithRule: 217,
        nCompared: 213,
        within: new Map([
          [0.01, 81],
          [0.02, 83],
          [0.05, 92],
        ]),
        medianDiff: -0.04,
        regions: { ahead: 50, within: 33, behind: 130 },
        uniform: { resistant: 43, susceptible: 38, other: 132 },
      }),
    ).toBe(
      "A known-gene rule needs no training: it calls a genome resistant when it carries a known resistance gene for the drug. Against the best learner with lineages held out, these rules are within 2 points or ahead on 83 of 213 datasets: 50 ahead by more than 2 and 33 within 2. They trail by more than 2 on the other 130.",
    );
    expect(micLead({ n: 126, median: 0.8362655650985171, nAbove90: 40 })).toBe(
      "The MIC (minimum inhibitory concentration) is the lowest drug level, in mg/L, that stops the bacterium growing. Labs test levels that double at each step, and a predicted MIC within one step of the measured one counts as essential agreement. For 40 of 126 datasets, at least 90% of exact MICs are predicted within one doubling dilution. The median dataset reaches 83.6%.",
    );
  });

  // Revision 12's index (274 datasets, 268 with lineages held out), as the
  // model of these summaries counts it.
  const alone = (learner: Fitted, medianRandom: number, medianLineage: number, medianGap: number, nNoCall: number, bestHigher: number, bestLower: number): LearnerRow => ({
    learner,
    n: 268,
    medianRandom,
    medianLineage,
    medianGap,
    nNoCall,
    nBest: 268,
    bestHigher,
    bestLower,
    medianBestLead: 0,
  });
  const rev12Learners: LearnerSummary = {
    rows: [
      alone("scm", 0.9156, 0.8708, 0.0242, 0, 106, 79),
      alone("cart", 0.9082, 0.8684, 0.0293, 0, 140, 68),
      alone("rf", 0.8955, 0.8124, 0.0504, 0, 181, 46),
      alone("gbt", 0.9119, 0.8459, 0.0463, 60, 141, 59),
      alone("lasso", 0.9153, 0.8567, 0.0376, 48, 135, 62),
    ],
    nBest: 268,
    bestAlsoTops: 123,
    otherAheadByCutoff: 68,
  };
  const rev12Split: SplitSummary = { n: 274, nNoBest: 0, nLineage: 268, medianRandom: 0.9286458630476951, medianLineage: 0.8790200150467228, medianGap: 0.03711346885361494, nLargeGap: 112 };
  const PICK =
    "Each dataset's best learner is the one of five with the highest random-split score among those that call every isolate. A learner's score with lineages held out takes no part in that pick. With lineages held out, the best learner's median balanced accuracy is 87.9% over 268 datasets.";

  it("sets the best learner beside the strongest learner used alone", () => {
    expect(learnersLead(rev12Learners, rev12Split)).toBe(
      `${PICK} For the set covering machine, used alone for every dataset, that median is 87.1%. The best learner scores higher than the set covering machine on 106 of 268 datasets and lower on 79.`,
    );
  });

  it("names only a learner that calls every isolate as the strongest alone", () => {
    // A higher median over fewer isolates does not stand beside the best learner's.
    const lassoHigher: LearnerSummary = { ...rev12Learners, rows: rev12Learners.rows.map((r) => (r.learner === "lasso" ? { ...r, medianLineage: 0.95 } : r)) };
    expect(learnersLead(lassoHigher, rev12Split)).toContain("For the set covering machine, used alone for every dataset, that median is 87.1%.");
    const nonePlain: LearnerSummary = { ...rev12Learners, rows: rev12Learners.rows.map((r) => ({ ...r, nNoCall: 1 })) };
    expect(learnersLead(nonePlain, rev12Split)).toBe(PICK);
    expect(learnersLead(rev12Learners, { ...rev12Split, nLineage: 0, medianLineage: null })).toBe("No dataset has a lineage-held-out result yet.");
  });

  it("reads both error rates against the bars with the three outcomes and their denominator", () => {
    expect(barsLead({ n: 268, meets: 0, misses: 239, open: 29, openAtBars: 6 })).toBe(
      "The site reads VME against at most 1.5% and ME against at most 3%, limits taken from FDA's 2009 guidance for growth-based susceptibility test systems. The guidance does not cover tests that read resistance from a genome, or tests of anti-mycobacterial drugs, so this site uses its numbers as a research screen adapted from it. They are not a criterion these models were validated against, and meeting them is not clearance of anything. A rate counted on few genomes can sit under its bar by chance, so each rate is read with its 95% interval. With lineages held out, the best learner is shown to meet both bars on 0 of 268 datasets and shown to miss at least one on 239. The other 29 have too few genomes to tell.",
    );
    expect(barsLead({ n: 3, meets: 1, misses: 2, open: 0, openAtBars: 0 })).toBe(
      "The site reads VME against at most 1.5% and ME against at most 3%, limits taken from FDA's 2009 guidance for growth-based susceptibility test systems. The guidance does not cover tests that read resistance from a genome, or tests of anti-mycobacterial drugs, so this site uses its numbers as a research screen adapted from it. They are not a criterion these models were validated against, and meeting them is not clearance of anything. A rate counted on few genomes can sit under its bar by chance, so each rate is read with its 95% interval. With lineages held out, the best learner is shown to meet both bars on 1 of 3 datasets and shown to miss at least one on 2.",
    );
    expect(barsLead({ n: 0, meets: 0, misses: 0, open: 0, openAtBars: 0 })).toBe("No dataset has both error rates with lineages held out yet.");
  });

  it("does not claim a selection cap the data does not show", () => {
    // The fixture's species use every labeled genome they have.
    expect(coverageNote(base)).toContain("In this run, every labeled genome that passed QC and has a FASTA file, the file that holds the DNA sequence, is selected.");
    expect(coverageNote(base)).not.toContain("capped");
  });

  it("names the known-gene catalog with its version, or says the version is not recorded", () => {
    expect(catalogWords({ ...base, amrfinder_version: "2026-08-14.1" })).toBe("Known-gene catalog: AMRFinderPlus 2026-08-14.1");
    expect(catalogWords({ ...base, amrfinder_version: "" })).toBe("Known-gene catalog: AMRFinderPlus version not recorded");
  });
});

// ---- the page --------------------------------------------------------------------

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  window.history.replaceState(null, "", "/");
  vi.unstubAllGlobals();
});

const fakeNav = () => {
  const calls: { go: Partial<Route>[]; open: string[] } = { go: [], open: [] };
  const nav: Navigator = {
    route: { view: "overview", species: null, dataset: null },
    go: (next) => {
      calls.go.push(next);
      window.history.pushState({}, "", next.view === "methods" ? "/methods" : next.view === "datasets" ? "/datasets" : `/species/${next.species ?? ""}`);
    },
    openDataset: (id) => calls.open.push(id),
    closeDataset: () => {},
  };
  return { nav, calls };
};

const mount = async (m: Manifest, nav: Navigator) => {
  const host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  await act(async () => root!.render(createElement(OverviewView, { manifest: m, nav })));
  return host;
};

const click = (el: Element) => el.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, button: 0 }));

/** A section's control by the words on it, and a table by its frame's name. */
const control = (section: Element, text: string) => [...section.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent === text)!;
const framed = (section: Element, label: string) => section.querySelector(`[role="region"][aria-label="${label}"]`)!;

/** The fixture's rows carry no audit, so they read as provisional; this
 *  copy passes the current audit. */
const passing: Manifest = { ...base, evaluation_revision: 7, datasets: base.datasets.map((d) => ({ ...d, evaluation: { revision: 7, valid: true, reasons: [] } })) };

describe("What GRM found", () => {
  it("reads as a report: a plain title, what the site asks, a plain status line, the primer, and eight numbered sections", async () => {
    const host = await mount(passing, fakeNav().nav);
    expect(host.querySelector("h1")?.textContent).toBe("What GRM found");
    // The question comes in words before any chart.
    expect(host.querySelector("header")?.textContent).toContain("asks how well the bacterium's genome, its full DNA sequence, predicts that lab result.");
    expect(host.querySelector("header")?.textContent).toContain("Known-gene catalog: AMRFinderPlus");
    expect([...host.querySelectorAll("h2")].map((h) => h.textContent)).toEqual([
      "How to read this report",
      "1 Genomes and lab calls",
      "2 Scores on unseen lineages",
      "3 The five learners compared",
      "4 Error rates against the bars",
      "5 Known-gene rules against learners",
      "6 Simulated PCR panel",
      "7 Predicting the drug level that stops growth (MIC)",
      "8 Limits",
    ]);
    const scope = host.querySelector<HTMLAnchorElement>('header a[href="#limits"]');
    expect(scope?.textContent).toBe("Tested within one public database only");
    // No chip, no hover wall, no Method and coverage fold.
    expect(host.querySelector('header [tabindex="0"]')).toBeNull();
    expect(host.querySelector('[role="tooltip"]')).toBeNull();
    expect(host.textContent).not.toContain("Method and coverage");
  });

  it("names a provisional audit in bold with its reasons in view, instead of the scope link", async () => {
    const reason = "Computed under evaluation audit revision 6; the current revision is 7, so these metrics are provisional until the dataset is recomputed.";
    const m: Manifest = { ...base, evaluation_revision: 7, datasets: base.datasets.map((d, i) => ({ ...d, evaluation: { revision: i === 0 ? 6 : 7, valid: true, reasons: [] } })) };
    const host = await mount(m, fakeNav().nav);
    expect(host.querySelector("header strong")?.textContent).toBe("Provisional metrics");
    expect(host.querySelector('header a[href="#limits"]')).toBeNull();
    expect(host.querySelector('[role="note"]')?.textContent).toContain(`1 of ${base.datasets.length} dataset evaluations do not pass the evaluation audit`);
    expect(host.querySelector('[role="note"]')?.textContent).toContain(reason);
  });

  it("gives every corpus count a visible action that goes where it says", async () => {
    const { nav, calls } = fakeNav();
    const host = await mount(base, nav);
    const corpus = host.querySelector("#corpus")!;
    const links = [...corpus.querySelectorAll<HTMLAnchorElement>("a")].map((a) => [a.textContent?.replace("→", "").trim(), a.getAttribute("href")]);
    expect(links).toEqual([
      ["See pair curation", "/methods#m-funnel"],
      ["See QC outcomes", "/methods#m-funnel"],
      ["Browse all report datasets", "/datasets"],
    ]);
    await act(async () => click(corpus.querySelector('a[href="/methods#m-funnel"]')!));
    // A router that keeps hashes takes it with the path; this one does not,
    // so the link writes it onto the entry the router pushed.
    expect(calls.go).toEqual([{ view: "methods", species: null, hash: "m-funnel" }]);
    expect(`${window.location.pathname}${window.location.hash}`).toBe("/methods#m-funnel");
    await act(async () => click(corpus.querySelector('a[href="/datasets"]')!));
    expect(calls.go[1]).toEqual({ view: "datasets", species: null });
  });

  it("opens the genomes by species in place, each analyzed species a link, with both denominators named", async () => {
    const { nav, calls } = fakeNav();
    const host = await mount(base, nav);
    const toggle = [...host.querySelectorAll("button")].find((b) => b.textContent?.startsWith("Show genomes by species"))!;
    expect(toggle.getAttribute("aria-expanded")).toBe("false");
    await act(async () => click(toggle));
    expect(toggle.getAttribute("aria-expanded")).toBe("true");
    const table = document.getElementById(toggle.getAttribute("aria-controls")!)!;
    expect([...table.querySelectorAll("thead th")].map((th) => th.textContent)).toEqual(["Species", "Reports", "Genomes analyzed", "Genomes with any label"]);
    const first = base.species.find((s) => s.genomes_analyzed > 0)!;
    const link = table.querySelector<HTMLAnchorElement>(`a[href="/species/${first.slug}"]`)!;
    expect(link.textContent).toBe(first.species);
    await act(async () => click(link));
    expect(calls.go).toEqual([{ view: "species", species: first.slug }]);
    expect(table.textContent).toContain("A report's cohort, the genomes it scores, is smaller: the selected genomes with a reported R or S call and no conflicting call for its drug.");
  });

  it("explains datasets without a lineage-held-out score in words, each name opening its report", async () => {
    const [a, b, c, ...rest] = base.datasets;
    const unpaired = (row: CatalogRow, reason: string): CatalogRow => ({ ...row, best_random: 0.9, best_lineage: null, lineage_unavailable_reason: reason });
    const rows = [
      unpaired(a!, "invalid lineage outer evaluation: fold 2 validation lacks one or both classes"),
      unpaired(b!, "invalid lineage outer evaluation: fold 4 validation lacks one or both classes"),
      unpaired(c!, "fewer than twice as many lineages as requested folds"),
      ...rest.map((r) => ({ ...r, best_random: r.best_random ?? 0.9, best_lineage: r.best_lineage ?? 0.85 })),
    ];
    const { nav, calls } = fakeNav();
    const host = await mount({ ...base, datasets: rows }, nav);
    const lineage = host.querySelector("#lineage")!;
    expect(lineage.textContent).toContain(
      "Three datasets have no score with lineages held out and are left out of these medians: two because a held-out test group (fold) would hold only resistant or only susceptible genomes and one because there are too few lineages.",
    );
    expect(lineage.textContent).not.toContain("lacks one or both classes");
    expect(lineage.textContent).not.toContain("requested folds");
    const name = lineage.querySelector<HTMLAnchorElement>(`a[href="/?dataset=${encodeURIComponent(c!.id)}"]`)!;
    await act(async () => click(name));
    expect(calls.open).toEqual([c!.id]);
  });

  it("sets the five learners under the best learner, with the dataset-by-dataset counts one control away", async () => {
    const host = await mount(base, fakeNav().nav);
    const learners = host.querySelector("#learners")!;
    const medians = framed(learners, "Median balanced accuracy by learner");
    expect([...medians.querySelectorAll("thead th")].map((th) => th.textContent)).toEqual(["Learner", "Datasets", "Random split", "Lineages held out", "Median drop, points"]);
    expect([...medians.querySelectorAll('tbody th[scope="row"]')].map((th) => th.textContent)).toEqual([
      "Best learner for each dataset",
      "Set covering machine",
      "Decision tree (CART)",
      "Random forest",
      "Gradient boosted trees",
      "L1 logistic regression",
    ]);
    expect(learners.querySelectorAll("table")).toHaveLength(1);
    await act(async () => click(control(learners, "Show best learner against each learner")));
    expect(framed(learners, "The best learner against each learner used alone").querySelectorAll('tbody th[scope="row"]')).toHaveLength(5);
    expect(learners.textContent).toContain("Where that learner is itself the best learner, the two score the same.");
  });

  it("counts the best learner's datasets against the bars in three outcomes that add up", async () => {
    const host = await mount(base, fakeNav().nav);
    const bars = host.querySelector("#bars")!;
    const tiles = [...bars.querySelectorAll("ol li")];
    expect(tiles.map((li) => li.lastElementChild?.textContent).slice(1)).toEqual(["Shown to meet both bars", "Shown to miss at least one bar", "Too few genomes to tell"]);
    expect(tiles[0]?.lastElementChild?.textContent).toMatch(/^Datasets? tested with lineages held out$/);
    const counts = tiles.map((li) => Number((li.firstElementChild?.textContent ?? "").replace(/,/g, "")));
    expect(counts).toHaveLength(4);
    expect(counts[0]).toBeGreaterThan(0);
    expect((counts[1] ?? 0) + (counts[2] ?? 0) + (counts[3] ?? 0)).toBe(counts[0]);
    // The reading is stated beside the counts, with its limit.
    expect(bars.textContent).toContain("The interval is exact (Clopper-Pearson)");
    expect(bars.textContent).toContain("so the true interval is likely wider.");
    await act(async () => click(control(bars, "Show each learner against the bars")));
    const each = framed(bars, "Each learner against the bars");
    expect([...each.querySelectorAll('th[scope="rowgroup"]')].map((th) => th.textContent)).toEqual(["Lineages held out", "Random split"]);
    expect(each.querySelectorAll('tbody th[scope="row"]')).toHaveLength(12);
    // The best learner's row with lineages held out repeats the four counts above it.
    expect([...each.querySelectorAll("tbody tr")][1]?.textContent).toBe(`Best learner for each dataset${counts.join("")}`);
  });

  it("groups the known-gene rules by their own score behind one control", async () => {
    const host = await mount(base, fakeNav().nav);
    const catalog = host.querySelector("#catalog")!;
    await act(async () => click(control(catalog, "Show rules grouped by their own score")));
    const groups = framed(catalog, "Known-gene rules grouped by their own balanced accuracy");
    expect([...groups.querySelectorAll('tbody th[scope="row"]')].map((th) => th.textContent)).toEqual(["Under 55%", "55% to under 70%", "70% to under 85%", "85% or more"]);
    // Every rule falls in one group.
    const inGroups = [...groups.querySelectorAll("tbody tr")].map((tr) => Number(tr.children[1]?.textContent));
    expect(inGroups.reduce((a, b) => a + b, 0)).toBe(base.datasets.filter((d) => d.catalog !== null).length);
    expect(catalog.textContent).toContain("A rule that calls every genome one way scores exactly 50%, a coin toss, so the first group ends just above that, at 55%.");
  });

  it("ends with the limits and a link to the audit in Methods", async () => {
    const { nav } = fakeNav();
    const host = await mount(base, nav);
    const limits = host.querySelector("#limits")!;
    expect(limits.textContent).toContain("internal cross-validation on public BV-BRC genomes");
    expect(limits.textContent).toContain("No external, prospective, or clinical validation has been done:");
    expect(limits.querySelector("a")?.getAttribute("href")).toBe("/methods#m-audit");
  });

  it("keeps every panel result behind one explicit control, open on arrival at #panels", async () => {
    const closed = await mount(base, fakeNav().nav);
    // The page says what a PCR panel is before its numbers.
    expect(closed.querySelector("#panels p")?.textContent).toMatch(/^A PCR panel is a lab test that looks for a fixed list of resistance targets in a sample's DNA\. /);
    const toggle = [...closed.querySelectorAll("#panels button")].find((b) => b.hasAttribute("aria-expanded"))!;
    expect(toggle.textContent).toBe("Show all panel results");
    expect(closed.querySelector('#panels [aria-label="Panel results"]')).toBeNull();
    act(() => root?.unmount());
    window.history.replaceState(null, "", "/#panels");
    const open = await mount(base, fakeNav().nav);
    expect([...open.querySelectorAll("#panels button")].find((b) => b.hasAttribute("aria-expanded"))?.getAttribute("aria-expanded")).toBe("true");
  });
});

describe("Data & Methods", () => {
  const html = renderToStaticMarkup(createElement(MethodsView, { manifest: base }));
  const el = document.createElement("div");
  el.innerHTML = html;

  it("is four linear sections, with no pills, nested disclosures, or data browser", () => {
    expect([...el.querySelectorAll("h2")].map((h) => h.textContent)).toEqual(["1 Labels and inclusion", "2 Models and evaluation", "3 How to read scores", "4 Scope, limitations and audit"]);
    expect(el.querySelector("[aria-pressed]")).toBeNull();
    expect(el.querySelector("[aria-expanded]")).toBeNull();
    expect(el.querySelector('input[type="search"]')).toBeNull();
    expect(el.querySelector('[aria-label="Antibiotic names"], [aria-label="Breakpoints"]')).toBeNull();
    expect(html).not.toContain("SigDrug");
    expect(html).not.toContain("Antibiotic name resolution");
  });

  it("keeps the anchors other views link to", () => {
    for (const id of ["m-funnel", "m-audit", "m-limits", "m-labels", "m-models", "m-scores"]) expect(el.querySelector(`#${id}`)).not.toBeNull();
    expect(el.querySelector("#m-funnel h3")?.getAttribute("tabindex")).toBe("-1");
  });

  it("words the scores with their called denominators, both kinds of no-call, and the pre-revision-7 behavior", () => {
    const text = el.textContent ?? "";
    expect(text).toContain("over the resistant isolates called");
    expect(text).toContain("A probe assay gives no call for an isolate whose genome carries none of its probe k-mers.");
    const seven = document.createElement("div");
    seven.innerHTML = renderToStaticMarkup(createElement(MethodsView, { manifest: { ...base, evaluation_revision: 7 } }));
    expect(seven.textContent).toContain("Reports computed under earlier revisions scored those isolates 0.5");
    expect(text).toContain("It does not run the AMRFinderPlus program or test the point mutations");
    expect(text).toContain("not CLSI or EUCAST breakpoints");
  });

  it("leaves learner settings to the dataset reports", () => {
    const text = el.textContent ?? "";
    expect(text).not.toContain("200 randomized trees");
    expect(text).not.toContain("early stopping");
    expect(text).toContain("Each dataset's report gives its learner settings");
  });

  it("lands a link on its section and moves focus to the section's heading", async () => {
    vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path))));
    window.history.replaceState(null, "", "/methods#m-funnel");
    const host = document.createElement("div");
    document.body.append(host);
    root = createRoot(host);
    await act(async () => root!.render(createElement(MethodsView, { manifest: base })));
    expect(document.activeElement?.textContent).toBe("Curation funnel");
  });
});
