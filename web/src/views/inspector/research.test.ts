// @vitest-environment happy-dom

import { act, createElement, type ReactElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, describe, expect, it, vi } from "vitest";
import type { Navigator } from "../../app/router";
import { fixtureFile } from "../../dev/fixtures/files";
import { categoryLabel, MIC_LIMITS, micEntry, micRun, researchWith, SELECTIVE_LIMITS, selectiveEntry, WITHHELD_REASON } from "../../dev/fixtures/research";
import { PATHS, type Resource } from "../../data/load";
import type { DatasetReport, MicDataset, Research, SelectiveCategory } from "../../data/schema";
import { parseDataset, parseManifest } from "../../data/validate";
import { ThemeProvider } from "../../theme/ThemeProvider";
import { SelectiveResearch } from "./SelectiveResearch";
import { DatasetInspector } from "./DatasetInspector";
import { MicSection } from "./MicSection";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const report = (id: string): DatasetReport => parseDataset(fixtureFile(`/data/datasets/${id}.json`));
/** A report with a saved MIC model, and one without (so its guard holds no
 *  saved predictions). */
const modeled = report("klebsiella_pneumoniae__meropenem");
const unmodeled = report("klebsiella_pneumoniae__ceftriaxone");
const copy = <T,>(v: T): T => JSON.parse(JSON.stringify(v)) as T;
const ready = (data: Research): Resource<Research> => ({ status: "ready", data });

/** research.json states that carry no research for any report. */
const NO_RESEARCH: Record<string, Resource<Research>> = {
  loading: { status: "loading" },
  missing: { status: "error", error: { kind: "missing", path: PATHS.research } },
  invalid: { status: "error", error: { kind: "shape", path: PATHS.research, message: 'research.json mic: "withheld" should be an array' } },
};

const render = (element: ReactElement): HTMLElement => {
  const host = document.createElement("div");
  host.innerHTML = renderToStaticMarkup(element);
  return host;
};
const mic = (r: DatasetReport, research: Resource<Research>) => render(createElement(MicSection, { r, research }));
const abstain = (r: DatasetReport, research: Resource<Research>) => render(createElement(SelectiveResearch, { r, research }));

/** The cells of each row of the table in the named scroll region. */
const table = (host: HTMLElement, label: string): string[][] =>
  [...(host.querySelector(`[role="region"][aria-label="${label}"]`)?.querySelectorAll("tr") ?? [])].map((tr) => [...tr.children].map((c) => c.textContent ?? ""));
const notes = (host: HTMLElement): string[] => [...host.querySelectorAll("li")].map((li) => li.textContent ?? "");
/** The research block: everything from its heading on. */
const block = (host: HTMLElement, heading: string): string => {
  const h3 = [...host.querySelectorAll("h3")].find((h) => h.textContent?.startsWith(heading));
  return h3?.parentElement?.textContent ?? "";
};
/** The paragraphs of the research block, in order. */
const paragraphs = (host: HTMLElement, heading: string): string[] => {
  const h3 = [...host.querySelectorAll("h3")].find((h) => h.textContent?.startsWith(heading));
  return [...(h3?.parentElement?.querySelectorAll("p") ?? [])].map((p) => p.textContent ?? "");
};

const HEADING = "No-genome baseline";
const BASELINE = "Essential agreement against a no-genome baseline";
const RESULTS = "Selective calling results";
const RERUN_INTRO =
  "The revision-5 MIC research refit GRM's MIC model on these isolates and scored the no-genome baseline on the same folds. Essential agreement is the share of exact MICs predicted within one doubling dilution.";
const UNUSED = "Revision-5 MIC research exists for this dataset, but";
const NONE_APPLY = "so none of its numbers apply to this report.";

describe("MIC research in the inspector", () => {
  const withEntry = (edit: (d: MicDataset) => void = () => {}, r = modeled) => {
    const entry = micEntry(r);
    edit(entry);
    return ready(researchWith([entry], []));
  };
  /** What the block says after its heading, for a report the research does not describe. */
  const unused = (r: DatasetReport, research: Resource<Research>) => {
    const host = mic(r, research);
    expect(table(host, BASELINE)).toEqual([]);
    expect(notes(host)).toEqual([]);
    return block(host, HEADING).replace(HEADING, "");
  };

  it("shows GRM's essential agreement beside the no-genome baseline on both splits", () => {
    const host = mic(modeled, withEntry());
    expect(paragraphs(host, HEADING)).toEqual([RERUN_INTRO]);
    expect(table(host, BASELINE)).toEqual([
      ["Split", "GRM model", "No-genome baseline"],
      ["Random split", "80.0% of 60", "70.0% of 60"],
      ["Lineages held out", "60.0% of 60", "50.0% of 60"],
    ]);
    expect(notes(host)).toEqual([MIC_LIMITS.scope, MIC_LIMITS.baseline]);
    // The report's own MIC result still leads the section.
    expect(host.textContent).toContain(`for ${(modeled.mic!.essential_agreement * 100).toFixed(1)}% of the ${modeled.mic!.n_exact} isolates with an exact MIC`);
  });

  it("keeps a saved report model's random-split row with no baseline", () => {
    const host = mic(modeled, ready(researchWith([micEntry(modeled, true)], [])));
    expect(paragraphs(host, HEADING)).toEqual([
      "The revision-5 MIC research scored this report's saved MIC predictions. Essential agreement is the share of exact MICs predicted within one doubling dilution.",
    ]);
    expect(table(host, BASELINE)).toEqual([
      ["Split", "GRM model", "No-genome baseline"],
      ["Random split", "80.0% of 60", "Not run"],
    ]);
    expect(notes(host)).toEqual([MIC_LIMITS.scope, MIC_LIMITS.saved_only]);
  });

  it("says where the no-genome baseline beats GRM's model, with both values and their denominator", () => {
    const wins = (edit: (d: MicDataset) => void) => paragraphs(mic(modeled, withEntry(edit)), HEADING).slice(1);
    // GRM's model leads on both splits in the fixture, and a tie is no win.
    expect(wins(() => {})).toEqual([]);
    expect(wins((d) => (d.baseline!.random = micRun(0.8)))).toEqual([]);
    expect(wins((d) => (d.baseline!.lineage = micRun(0.65)))).toEqual([
      "The no-genome baseline has higher essential agreement than GRM's model with lineages held out: 65.0% against 60.0% of the same 60 exact MICs.",
    ]);
    expect(
      wins((d) => {
        d.baseline!.random = micRun(0.9);
        d.baseline!.lineage = micRun(0.65);
      }),
    ).toEqual([
      "The no-genome baseline has higher essential agreement than GRM's model on both splits: 90.0% against 80.0% on the random split and 65.0% against 60.0% with lineages held out, of the same 60 exact MICs.",
    ]);
    // Each split names its own denominator when the two differ.
    expect(
      wins((d) => {
        d.baseline!.random = micRun(0.9);
        d.model.lineage = { ...micRun(0.6, 59), holdout_random_folds: 0 };
        d.baseline!.lineage = micRun(0.65, 59);
      }),
    ).toEqual([
      "The no-genome baseline has higher essential agreement than GRM's model on both splits: 90.0% against 80.0% of 60 exact MICs on the random split and 65.0% against 60.0% of 59 exact MICs with lineages held out.",
    ]);
    // 1,390 and 1,391 of 1,403 both round to 99.1%.
    expect(
      wins((d) => {
        d.model.random = micRun(1390 / 1403, 1403);
        d.baseline!.random = micRun(1391 / 1403, 1403);
      }),
    ).toEqual(["The no-genome baseline has higher essential agreement than GRM's model on the random split: 99.14% against 99.07% of the same 1,403 exact MICs."]);
  });

  it("adds each limit only where its condition holds", () => {
    // Every rerun carries the scope and baseline limits first.
    const extra = (edit: (d: MicDataset) => void) => notes(mic(modeled, withEntry(edit))).slice(2);
    const holdout = (folds: string) => `In ${folds}, the early-stopping holdout (training isolates set aside to decide when boosting stops) fell back to random isolates instead of whole lineages.`;
    const exact = (n: number) => (d: MicDataset) => {
      d.model.lineage!.n_exact = n;
      d.baseline!.lineage.n_exact = n;
    };
    expect(extra(exact(49))).toEqual([MIC_LIMITS.exact_minimum]);
    expect(extra(exact(50))).toEqual([]);
    expect(extra((d) => (d.zone_diameters.removed = 14))).toEqual(["MIC zone-fix limit for 14 disk-diffusion zone diameters."]);
    expect(extra((d) => (d.zone_diameters.removed = 1))).toEqual(["MIC zone-fix limit for 1 disk-diffusion zone diameter."]);
    expect(extra((d) => (d.model.lineage!.holdout_random_folds = 5))).toEqual([holdout("5 lineage-split folds")]);
    expect(extra((d) => (d.model.lineage!.holdout_random_folds = 1))).toEqual([holdout("1 lineage-split fold")]);
  });

  it("words a run without essential agreement instead of printing a number", () => {
    const host = mic(
      modeled,
      withEntry((d) => {
        d.model.random = micRun(null, 0);
        d.baseline!.random = micRun(null, 0);
      }),
    );
    expect(table(host, BASELINE)[1]).toEqual(["Random split", "No exact MIC scored", "No exact MIC scored"]);
    expect(paragraphs(host, HEADING)).toEqual([RERUN_INTRO]);
  });

  it("gives a report without an MIC model no research numbers, even when its guard matches", () => {
    // The entry guards this report exactly, as on K. pneumoniae ciprofloxacin.
    const host = mic(unmodeled, withEntry(() => {}, unmodeled));
    expect(host.textContent).toContain("No MIC regression for this dataset");
    expect(table(host, BASELINE)).toEqual([]);
    expect(notes(host)).toEqual([]);
    expect(block(host, HEADING)).toBe(`${HEADING}${UNUSED} it scored GRM's MIC model and this report has none, ${NONE_APPLY}`);
  });

  it("names the missing MIC model on each kind of revision-6 report without one", () => {
    const research = withEntry();
    const without = (edit: (r: DatasetReport) => void = () => {}) => {
      const r = copy(modeled);
      r.mic = null;
      edit(r);
      return unused(r, research);
    };
    const changeCohort = (r: DatasetReport) => {
      r.samples[0]!.id = "573.900001";
      r.samples.forEach((s) => {
        s.resistant = !s.resistant;
        s.lineage += 1;
      });
      r.samples[r.samples.findIndex((s) => s.mic !== null)]!.mic!.log2 += 1;
    };
    // Only the model is gone (S. enterica ampicillin).
    expect(without()).toBe(`${UNUSED} it scored GRM's MIC model and this report has none, ${NONE_APPLY}`);
    // The isolates changed too, with or without a dropped genome (E. coli
    // aztreonam and piperacillin/tazobactam; S. enterica ceftriaxone and
    // sulfamethoxazole/trimethoprim).
    const cohort = `${UNUSED} it scored GRM's MIC model and this report has none, and this report's isolates, resistant and susceptible labels, lineage clusters, and measured MICs differ from the ones it read, ${NONE_APPLY}`;
    expect(without(changeCohort)).toBe(cohort);
    expect(
      without((r) => {
        changeCohort(r);
        r.samples.pop();
      }),
    ).toBe(cohort);
    // Neither revision has a model and a genome is gone (E. coli
    // ceftriaxone), which changes every part but the predictions.
    const lost = copy(unmodeled);
    lost.samples.pop();
    expect(unused(lost, withEntry(() => {}, unmodeled))).toBe(cohort);
  });

  it("names a model the revision-5 report lacked without describing predictions it never had", () => {
    const research = withEntry((d) => (d.guard.mic_predicted = null));
    expect(unused(modeled, research)).toBe(`${UNUSED} this report has an MIC model and the revision-5 report it read had none, ${NONE_APPLY}`);
    const moved = copy(modeled);
    moved.samples[0]!.id = "573.900001";
    expect(unused(moved, research)).toBe(`${UNUSED} this report has an MIC model and the revision-5 report it read had none, and this report's isolates differ from the ones it read, ${NONE_APPLY}`);
  });

  it("shows a withheld model's reason alone, whatever the report's cohort or model", () => {
    const later = copy(modeled);
    later.samples[0]!.id = "573.900001";
    for (const r of [modeled, later, unmodeled]) {
      const host = mic(r, ready(researchWith([], [], [{ id: r.id, reason: WITHHELD_REASON }])));
      expect(block(host, HEADING)).toBe(`${HEADING}${WITHHELD_REASON}`);
      expect(table(host, BASELINE)).toEqual([]);
      expect(notes(host)).toEqual([]);
    }
  });

  it("gives a report from another cohort no numbers and names what differs", () => {
    const research = withEntry();
    const sentence = (edit: (r: DatasetReport) => void) => {
      const r = copy(modeled);
      edit(r);
      return unused(r, research);
    };
    expect(sentence((r) => (r.mic!.predicted[r.mic!.predicted.findIndex((p) => p !== null)]! += 1))).toBe(`${UNUSED} this report's saved MIC predictions differ from the ones it read, ${NONE_APPLY}`);
    expect(
      sentence((r) => {
        r.samples[0]!.id = "573.900001";
        r.samples[r.samples.findIndex((s) => s.mic !== null)]!.mic!.log2 += 1;
      }),
    ).toBe(`${UNUSED} this report's isolates and measured MICs differ from the ones it read, ${NONE_APPLY}`);
    // A revision that drops a genome changes every value the lane read.
    expect(
      sentence((r) => {
        r.samples.pop();
        r.mic!.predicted.pop();
        r.samples.forEach((s) => {
          s.resistant = !s.resistant;
          s.lineage += 1;
        });
      }),
    ).toBe(`${UNUSED} this report's isolates, resistant and susceptible labels, lineage clusters, measured MICs, and saved MIC predictions differ from the ones it read, ${NONE_APPLY}`);
  });

  it("shows nothing from research while research.json loads or when it is missing, invalid, or has no entry", () => {
    const without = mic(modeled, ready(researchWith([], []))).innerHTML;
    expect(without).not.toMatch(/research|revision|baseline/i);
    for (const state of Object.values(NO_RESEARCH)) expect(mic(modeled, state).innerHTML).toBe(without);
  });
});

describe("selective-calling research in the inspector", () => {
  const selective = (category: SelectiveCategory, r = modeled) => abstain(r, ready(researchWith([], [selectiveEntry(r, category)])));
  /** The fixture report passes the evaluation audit, so the conformal body shows. */
  const audited = { ...copy(modeled), evaluation: { revision: 6, valid: true, reasons: [] } };
  const MISMATCH = "Selective calling, where the genome alone would call a confident subset of isolates and the rest would go to the lab, has not run on this revision's held-out scores.";

  it("says selective calling has not run on a report whose held-out scores differ, with no numbers", () => {
    const research = ready(researchWith([], [selectiveEntry(modeled, "A1")]));
    const rescored = copy(modeled);
    rescored.splits.forEach((s) => s.learners.forEach((l) => (l.scores = l.scores.map((x) => x === null ? null : x / 2))));
    const host = abstain(rescored, research);
    expect(host.textContent).toContain(MISMATCH);
    expect(host.textContent).not.toMatch(/\d/);
    expect(table(host, RESULTS)).toEqual([]);
    // The audited conformal body carries the same sentence after its own table.
    const rescoredAudited = { ...rescored, evaluation: audited.evaluation };
    expect(abstain(rescoredAudited, research).textContent).toContain(MISMATCH);
  });

  it("lays out an A1 result: the bar, coverage and no call, VME and ME with counts and bounds, per-100 rows, and limits", () => {
    for (const r of [modeled, audited]) {
      const host = selective("A1", r);
      expect(paragraphs(host, "Selective calling")).toEqual([
        "Category A1.",
        SELECTIVE_LIMITS.bar,
        "With lineages held out, the genome alone calls 86.0% of the 100 isolates and gives no call for the other 14.0%, which would go to the lab.",
      ]);
      expect(table(host, RESULTS)).toEqual([
        ["Measure", "Resistant", "Susceptible"],
        ["Isolates", "40", "60"],
        ["Given a genome call", "80.0%", "90.0%"],
        ["Wrong among called", "VME 6.3%, 2 of 32", "ME 7.4%, 4 of 54"],
        ["95% confidence interval", "0.8% to 20.8%", "2.1% to 17.9%"],
        ["Per 100 isolates"],
        ["Right call", "75.0 caught", "83.3 cleared"],
        ["Wrong call", "5.0 missed", "6.7 overcalled"],
        ["No call, would go to the lab", "20.0", "10.0"],
      ]);
      expect(notes(host)).toEqual([SELECTIVE_LIMITS.scope, SELECTIVE_LIMITS.no_call, SELECTIVE_LIMITS.optimism, SELECTIVE_LIMITS.procedure, SELECTIVE_LIMITS.per_100]);
    }
  });

  it("names an A2 result's refit learner once and leaves out the A1 procedure limit", () => {
    const host = selective("A2");
    expect(paragraphs(host, "Selective calling").slice(0, 2)).toEqual(["Category A2 with gradient boosted trees.", SELECTIVE_LIMITS.bar]);
    expect(block(host, "Selective calling").match(/gradient boosted trees|refit learner/g)).toEqual(["gradient boosted trees"]);
    expect(table(host, RESULTS)[3]).toEqual(["Wrong among called", "VME 6.3%, 2 of 32", "ME 7.4%, 4 of 54"]);
    expect(notes(host)).toEqual([SELECTIVE_LIMITS.scope, SELECTIVE_LIMITS.no_call, SELECTIVE_LIMITS.optimism, SELECTIVE_LIMITS.per_100]);
  });

  it("gives B, C, D1, and D2 their reason and the bar, and no lineage split its reason alone, none with numbers", () => {
    const shown = (category: SelectiveCategory) => {
      const host = selective(category);
      expect(table(host, RESULTS)).toEqual([]);
      expect(notes(host)).toEqual([]);
      return paragraphs(host, "Selective calling");
    };
    for (const code of ["B", "C", "D1", "D2"] as const) expect(shown(code)).toEqual([categoryLabel(code), SELECTIVE_LIMITS.bar]);
    expect(shown("no_lineage")).toEqual([categoryLabel("no_lineage")]);
  });

  it("words the lab as hypothetical wherever the view mentions it", () => {
    const research = ready(researchWith([], [selectiveEntry(modeled, "A1")]));
    const rescored = copy(modeled);
    rescored.splits.forEach((s) => s.learners.forEach((l) => (l.scores = l.scores.map((x) => x === null ? null : x / 2))));
    for (const host of [selective("A1"), selective("A2"), abstain(rescored, research)]) {
      const lab = host.textContent?.match(/\b\w+ (?:\w+ )?(?:go|goes|sent|send|sends) to the lab/g) ?? [];
      expect(lab.length).toBeGreaterThan(0);
      expect(lab.every((m) => /would (?:\w+ )?go to the lab$/.test(m))).toBe(true);
    }
  });

  it("shows nothing from research while research.json loads or when it is missing, invalid, or has no entry", () => {
    for (const r of [modeled, audited]) {
      const without = abstain(r, ready(researchWith([], []))).innerHTML;
      expect(without).not.toMatch(/selective|research|revision/i);
      for (const state of Object.values(NO_RESEARCH)) expect(abstain(r, state).innerHTML).toBe(without);
    }
  });
});

describe("research in the dataset inspector", () => {
  let root: Root | null = null;

  afterEach(() => {
    act(() => root?.unmount());
    root = null;
    document.body.innerHTML = "";
    vi.unstubAllGlobals();
  });

  it("fetches research.json once beside the report and passes it to both sections", async () => {
    // The MIC entry guards the served report; the selective entry guards
    // other held-out scores, as on every revision-6 report.
    const rescored = copy(modeled);
    rescored.splits.forEach((s) => s.learners.forEach((l) => (l.scores = l.scores.map((x) => x === null ? null : x / 2))));
    const research = researchWith([micEntry(modeled)], [selectiveEntry(rescored, "A1")]);
    const requests: string[] = [];
    vi.stubGlobal("fetch", async (path: string) => {
      requests.push(path);
      return new Response(JSON.stringify(path === PATHS.research ? research : fixtureFile(path)), { headers: { "Content-Type": "application/json" } });
    });
    const manifest = parseManifest(fixtureFile("/data/index.json"));
    const nav: Navigator = { route: { view: "datasets", species: null, dataset: modeled.id }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
    const host = document.createElement("div");
    document.body.append(host);
    const r = createRoot(host);
    root = r;
    await act(async () => {
      r.render(createElement(ThemeProvider, { children: createElement(DatasetInspector, { id: modeled.id, manifest, nav }) }));
    });
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });

    expect(requests.filter((p) => p === PATHS.research)).toHaveLength(1);
    expect(table(document.body, BASELINE)[1]).toEqual(["Random split", "80.0% of 60", "70.0% of 60"]);
    expect(document.getElementById("ds-learners")?.textContent).toContain(
      "Selective calling, where the genome alone would call a confident subset of isolates and the rest would go to the lab, has not run on this revision's held-out scores.",
    );
  });
});
