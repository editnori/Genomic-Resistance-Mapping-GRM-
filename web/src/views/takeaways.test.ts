// @vitest-environment happy-dom
// The above-the-fold takeaways block, plus the two new study-split and
// other-database sections on Overview.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { Navigator } from "../app/router";
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

    // Clicking a chapter link moves focus to that section's heading.
    await act(async () => click(takeaways.querySelector('a[href="#outside"]')!));
    expect(document.activeElement).toBe(host.querySelector("#outside-title"));

    // Clicking the start-to-finish guide button opens the 4-part plain-English summary.
    const guideBtn = [...takeaways.querySelectorAll("button[aria-expanded]")].find((b) => b.textContent?.includes("start-to-finish guide"))!;
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
