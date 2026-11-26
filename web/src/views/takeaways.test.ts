// @vitest-environment happy-dom
// The above-the-fold takeaways block, plus the two new study-split and
// other-database sections on Overview.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { Navigator } from "../app/router";
import { classDropRows, MicPlateStrip, WholePartsBar } from "../charts/amr";
import { bottomLineLead, TAKEAWAYS_ANCHOR, takeawayCards } from "./Takeaways";
import { OverviewView } from "./OverviewView";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const base = parseManifest(fixtureFile("/data/index.json"));

let root: Root | null = null;
beforeEach(() => {
  vi.stubGlobal("fetch", async () => new Response("", { status: 404 }));
});
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  window.history.replaceState(null, "", "/");
  vi.unstubAllGlobals();
});

const nav: Navigator = {
  route: { view: "overview", species: null, dataset: null },
  go: () => {},
  openDataset: () => {},
  closeDataset: () => {},
};

const click = (el: Element) => el.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, button: 0 }));

describe("Takeaways and new Overview sections", () => {
  it("renders six takeaway cards, four chapter links, and the scope box above the primer", async () => {
    const host = document.createElement("div");
    document.body.append(host);
    root = createRoot(host);
    await act(async () => root!.render(createElement(OverviewView, { manifest: base, nav })));

    const takeaways = host.querySelector(`#${TAKEAWAYS_ANCHOR}`)!;
    const primer = host.querySelector("#primer")!;
    expect(takeaways.compareDocumentPosition(primer) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
    expect(takeaways.querySelector("h2")?.textContent).toBe("At a glance");
    expect(takeaways.querySelectorAll("[data-takeaway-card]")).toHaveLength(6);
    expect([...takeaways.querySelectorAll('nav[aria-label="Report chapters"] a')].map((a) => a.getAttribute("href"))).toEqual([
      "#corpus",
      "#catalog",
      "#outside",
      "#inputs",
    ]);
    expect(takeaways.querySelector('[data-primer-figure="scope"]')).not.toBeNull();
    expect(takeaways.querySelectorAll('[data-amr-primitive="bolster-strip"] li')).toHaveLength(4);
    const learnerDumbbells = [...host.querySelectorAll('#learners [data-amr-primitive="split-dumbbell"] li[data-dumbbell-row]')];
    expect(learnerDumbbells).toHaveLength(6);
    expect(learnerDumbbells.every((li) => (li.querySelector("[data-dumbbell-delta]")?.getAttribute("data-dumbbell-delta") ?? "").length > 0)).toBe(true);
    for (const sel of ["#bars", "#catalog", "#other-database"]) {
      const bar = host.querySelector(`${sel} [data-amr-primitive="whole-parts-bar"]`)!;
      expect(bar).not.toBeNull();
      const whole = Number(bar.querySelector("[data-whole-count]")?.getAttribute("data-whole-count"));
      const parts = [...bar.querySelectorAll("[data-part-count]")].map((e) => Number(e.getAttribute("data-part-count")));
      expect(parts.length).toBeLessThanOrEqual(4);
      expect(parts.reduce((a, b) => a + b, 0)).toBe(whole);
    }
    const micPlate = host.querySelector('#mic [data-amr-primitive="mic-plate"]')!;
    expect(micPlate.getAttribute("data-plate-mode")).toBe("agreement");
    expect(micPlate.textContent).toContain("±1 doubling well");
    await act(async () => click(primer.querySelector("button[aria-expanded]")!));
    const primerPlate = primer.querySelector('[data-amr-primitive="mic-plate"]')!;
    expect(primerPlate.getAttribute("data-plate-mode")).toBe("cutoff");
    expect(primerPlate.textContent).toContain("drawing, not data");
    expect(primerPlate.textContent).not.toContain("±1 doubling well");

    // Clicking a chapter link moves focus to that section's heading.
    await act(async () => click(takeaways.querySelector('a[href="#outside"]')!));
    expect(document.activeElement).toBe(host.querySelector("#outside-title"));

    // Clicking the guide button opens the 4-part walkthrough.
    const guideBtn = [...takeaways.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent === "Read the guide")!;
    expect(guideBtn.getAttribute("aria-expanded")).toBe("false");
    expect(takeaways.querySelector("[data-start-to-finish]")).toBeNull();
    await act(async () => click(guideBtn));
    expect(guideBtn.getAttribute("aria-expanded")).toBe("true");
    expect(takeaways.querySelectorAll("[data-start-to-finish] h3")).toHaveLength(4);

    // Wording guardrails on takeaways: no em-dashes, no external-validation overclaim.
    const text = takeaways.textContent ?? "";
    expect(text).not.toContain("\u2014");
    expect(text).toContain("No external or clinical validation has been done.");
    expect(text).toContain("growth-based susceptibility test systems");
    expect(text).toContain("12 of 25");
    expect(text).toContain("3 of 4 species where it can be read");
    expect(text).toContain("6 of 8 species (SCM)");
    expect(text).toContain("5 of 5 learners worse");
    expect(text).toContain("the best learner with lineages held out is shown to meet both on");
  });

  it("derives the bottom-line lead and MIC card dynamically with or without baselines", () => {
    expect(bottomLineLead({ ...base, datasets: [] })).toBe("No dataset has a report yet.");
    const noBase = takeawayCards(base);
    expect(noBase).toHaveLength(6);
    expect(noBase[5]!.body).not.toContain("constant training-fold median MIC");

    const micRow = base.datasets.find((d) => d.mic_essential_agreement !== null);
    if (micRow !== undefined) {
      const withBase = takeawayCards(base, new Map([[micRow.id, { constant: 0.75, model: micRow.mic_essential_agreement!, confirmed: true }]]));
      expect(withBase[5]!.body).toContain("constant training-fold median MIC");
    }

    // classDropRows takes the median of each dataset's own drop, not the difference of two medians.
    const synth = [
      { ...base.datasets[0]!, id: "c1", antibiotic: "ciprofloxacin", best_random: 0.9, best_lineage: 0.89 },
      { ...base.datasets[0]!, id: "c2", antibiotic: "ciprofloxacin", best_random: 0.8, best_lineage: 0.78 },
      { ...base.datasets[0]!, id: "c3", antibiotic: "ciprofloxacin", best_random: 0.95, best_lineage: 0.7 },
    ];
    const dropRows = classDropRows(synth, 3);
    expect(dropRows).toHaveLength(1);
    expect(dropRows[0]!.before).toBeCloseTo(0.9, 6);
    expect(dropRows[0]!.after).toBeCloseTo(0.78, 6);
    expect(dropRows[0]!.delta).toBeCloseTo(-0.02, 6);

    // WholePartsBar enforces <= 4 parts and exact sum equality with total.
    expect(() =>
      renderToStaticMarkup(
        createElement(WholePartsBar, {
          wholeLabel: "Bad total",
          total: 10,
          parts: [
            { key: "a", label: "A", count: 4, tone: "ink" },
            { key: "b", label: "B", count: 5, tone: "neutral" },
          ],
        }),
      ),
    ).toThrow(/must equal total/);
    expect(() =>
      renderToStaticMarkup(
        createElement(WholePartsBar, {
          wholeLabel: "Too many",
          total: 5,
          parts: [1, 2, 3, 4, 5].map((n) => ({ key: `p${n}`, label: `P${n}`, count: 1, tone: "neutral" as const })),
        }),
      ),
    ).toThrow(/at most 4 parts/);

    // MicPlateStrip prints unrounded cohort breakpoints and separates cutoff from agreement mode.
    const plateHtml = renderToStaticMarkup(createElement(MicPlateStrip, { mode: "agreement", rAbove: -1, cutoffLabel: "Cohort breakpoint" }));
    expect(plateHtml).toContain("Cohort breakpoint &gt; 0.5 mg/L");
    expect(plateHtml).toContain("±1 doubling well (agrees)");
  });

  it("renders Holding out whole studies and Another public database with expandable detail tables", async () => {
    const host = document.createElement("div");
    document.body.append(host);
    root = createRoot(host);
    await act(async () => root!.render(createElement(OverviewView, { manifest: base, nav })));

    const study = host.querySelector("#study-split")!;
    expect(study.textContent).toContain("retrospective split by study inside one database, BV-BRC (70 datasets, 9 species, 57,914 genome-drug pairs)");
    expect(study.textContent).toContain("never external validation");
    const studyBtn = [...study.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent === "Show repeated draws and matched fold shapes")!;
    expect(studyBtn.getAttribute("aria-expanded")).toBe("false");
    await act(async () => click(studyBtn));
    expect(studyBtn.getAttribute("aria-expanded")).toBe("true");
    expect(study.querySelector('[role="region"][aria-label="Second experiment across repeated fold draws"]')).not.toBeNull();

    const other = host.querySelector("#other-database")!;
    expect(other.textContent).toContain("13,377 NCBI Pathogen Detection genomes");
    expect(other.textContent).toContain("never external validation");
    expect(other.textContent).toContain("3 of the 4 species where it can be read (17 of its 24 datasets across the eight)");
    expect(other.textContent).toContain("decision tree on 12 of 25");
    expect(other.textContent).toContain("1,768 of 2,377");
    expect(other.textContent).toContain("948 of 1,002");
    const otherBtn = [...other.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent === "Show by species and bioproject overlap")!;
    expect(otherBtn.getAttribute("aria-expanded")).toBe("false");
    await act(async () => click(otherBtn));
    expect(otherBtn.getAttribute("aria-expanded")).toBe("true");
    expect(other.querySelector('[role="region"][aria-label="NCBI Pathogen Detection results by species and BioProject overlap"]')).not.toBeNull();

    const comp = host.querySelector("#comparisons")!;
    expect(comp.textContent).toContain("First-public-day split");
    expect(comp.textContent).toContain("Kover 2.0.7 vs GRM SCM");
    expect(comp.textContent).toContain("Dense neural network");
    expect(comp.textContent).toContain("AMRFinderPlus 4.2.7 program");
    expect(comp.textContent).toContain("TBProfiler 6.6.3 (M. tuberculosis)");
    const compBtn1 = [...comp.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent === "Show time split and dense network by learner")!;
    await act(async () => click(compBtn1));
    expect(comp.querySelector('[role="region"][aria-label="First-public-day split and dense neural network by learner"]')).not.toBeNull();
    const compBtn2 = [...comp.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent === "Show amrfinderplus and tbprofiler by learner")!;
    await act(async () => click(compBtn2));
    expect(comp.querySelector('[role="region"][aria-label="AMRFinderPlus 4.2.7 and TBProfiler 6.6.3 against each learner on lineage splits"]')).not.toBeNull();
  });
});
