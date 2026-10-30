// @vitest-environment happy-dom
// The primer: which real report it reads from, that every number it prints
// is that report's own, that it never claims more than internal
// cross-validation, and that Overview and Present both carry it with
// working navigation, keyboard focus and reduced motion.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { CatalogRow, Manifest } from "../data/schema";
import { int, pct, speciesShort } from "../data/format";
import type { Navigator, Route } from "../app/router";
import { exampleName, exampleNameInline, missedSentence, perHundred, primerExample, primerSteps, PRIMER_ANCHOR } from "./primer";
import { OverviewView } from "./OverviewView";
import { PresentMode } from "../present/PresentMode";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const base = parseManifest(fixtureFile("/data/index.json"));

const withLearner = (id: string, n: number, lineage: { sensitivity: number; specificity: number } | null): CatalogRow => ({
  ...base.datasets[0]!,
  id,
  antibiotic: id,
  n,
  n_resistant: Math.round(n / 2),
  best_learner: "rf",
  learners: [
    {
      learner: "rf",
      random: { balanced_accuracy: 0.95, sensitivity: 0.94, specificity: 0.96, vme: 0.06, me: 0.04, auroc: null },
      lineage:
        lineage === null
          ? null
          : { balanced_accuracy: (lineage.sensitivity + lineage.specificity) / 2, ...lineage, vme: 1 - lineage.sensitivity, me: 1 - lineage.specificity, auroc: null },
    },
  ],
});

/** Every visible text node of a rendered node, as one string. */
const text = (node: unknown): string => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(node as never);
  return el.textContent ?? "";
};

describe("the primer's example report", () => {
  it("prefers E. coli ceftriaxone when it has both splits, and reads every number from that row", () => {
    const ex = primerExample(base);
    expect(ex).not.toBeNull();
    const row = base.datasets.find((d) => d.id === "escherichia_coli__ceftriaxone");
    if (row !== undefined && row.learners.find((l) => l.learner === row.best_learner)?.lineage != null) expect(ex!.row.id).toBe(row.id);
    const steps = primerSteps(base, ex);
    expect(steps.map((s) => s.title)).toEqual(["A dataset", "A genome feature", "Two tests", "Reading one result", "What the numbers support"]);
    const name = exampleName(ex!.row);
    expect(steps[0]!.lead).toBe(`A dataset is one species tested against one antibiotic: ${name}, ${int(ex!.row.n)} genomes, ${int(ex!.row.n_resistant)} resistant.`);
    expect(steps[2]!.lead).toBe(`${name} scores ${pct(ex!.random.balanced_accuracy)} on random splits and ${pct(ex!.lineage.balanced_accuracy)} with whole lineages held out.`);
    const h = perHundred(ex!);
    expect(h.missed).toBe(Math.round((1 - ex!.lineage.sensitivity) * 100));
    expect(steps[3]!.lead).toContain(`called ${h.missed} of every 100 resistant genomes`);
  });

  it("falls back to the largest dataset with a lineage result, and to plain steps when none has one", () => {
    const m: Manifest = { ...base, datasets: [withLearner("small", 100, { sensitivity: 0.8, specificity: 0.9 }), withLearner("big", 900, { sensitivity: 0.7, specificity: 0.95 }), withLearner("unsplit", 5000, null)] };
    const ex = primerExample(m)!;
    expect(ex.row.id).toBe("big");
    expect(perHundred(ex)).toEqual({ missed: 30, caught: 70, flagged: 5, cleared: 95 });
    const none: Manifest = { ...base, datasets: [withLearner("unsplit", 5000, null)] };
    expect(primerExample(none)).toBeNull();
    const steps = primerSteps(none);
    expect(steps).toHaveLength(5);
    expect(steps[0]!.figure).toBeNull();
    expect(steps[3]!.figure).toBeNull();
    // No number stands in for a missing result.
    expect(steps.map((s) => s.lead).join(" ")).not.toMatch(/\d/);
  });

  it("draws the example's errors as filled squares out of 100, counted in words beside them", () => {
    const m: Manifest = { ...base, datasets: [withLearner("big", 900, { sensitivity: 0.83, specificity: 0.95 })] };
    const el = document.createElement("div");
    el.innerHTML = renderToStaticMarkup(primerSteps(m)[3]!.figure as never);
    const grids = el.querySelectorAll('[data-primer-figure="errors"] > div');
    expect(grids).toHaveLength(2);
    expect(grids[0]!.querySelectorAll("[data-error]")).toHaveLength(17);
    expect(grids[1]!.querySelectorAll("[data-error]")).toHaveLength(5);
    expect(el.textContent).toContain("17 called susceptible (VME)");
    expect(el.textContent).toContain("5 called resistant (ME)");
  });

  it("says which drug and species missed which genomes, in one grammatical sentence from the report's own rate", () => {
    const ex = primerExample(base)!;
    const lead = primerSteps(base, ex)[3]!.lead;
    expect(lead).toBe(`With lineages held out, the best learner for ${ex.row.antibiotic} in ${speciesShort(ex.row.species)} called ${perHundred(ex).missed} of every 100 resistant genomes susceptible.`);
    expect(lead).not.toMatch(/genomes in [A-Z]/);
    // The rate, drug and species follow the row; nothing is fixed at 7.
    for (const [antibiotic, species, sensitivity, missed] of [
      ["clindamycin", "Staphylococcus aureus", 0.99, 1],
      ["ciprofloxacin", "Klebsiella pneumoniae", 0.7, 30],
      ["gentamicin", "Escherichia coli", 1, 0],
    ] as const) {
      const row = { ...withLearner("x", 400, { sensitivity, specificity: 0.9 }), antibiotic, species };
      const m: Manifest = { ...base, datasets: [row] };
      const s = primerSteps(m)[3]!;
      expect(s.lead).toBe(`With lineages held out, the best learner for ${antibiotic} in ${speciesShort(species)} called ${missed} of every 100 resistant genomes susceptible.`);
      expect(missedSentence(primerExample(m)!)).toBe(s.lead);
      expect(`${s.lead} ${text(s.body)}`).not.toMatch(/patient|treat|clinical/i);
    }
  });

  it("draws both splits on the same tracks: name and score on one line, then the same five lineage groups", () => {
    const el = document.createElement("div");
    el.innerHTML = renderToStaticMarkup(primerSteps(base)[2]!.figure as never);
    const rows = [...el.querySelectorAll("[data-split-row]")];
    expect(rows).toHaveLength(2);
    for (const row of rows) {
      const head = row.firstElementChild!;
      expect(head.children).toHaveLength(2);
      expect(head.lastElementChild!.textContent).toMatch(/^\d+\.\d%$/);
      const groups = [...row.querySelectorAll("[data-lineage-group]")];
      expect(groups.map((g) => g.children.length)).toEqual([6, 4, 5, 3, 4]);
      // One row of groups with one fixed gap, not a grid that reflows per row.
      expect(groups[0]!.parentElement!.className).toContain("gap-x-2");
    }
    expect(rows[1]!.querySelectorAll("[data-lineage-group]")[2]!.querySelectorAll("[data-hidden]")).toHaveLength(5);
    expect(el.textContent).toContain(`scores are ${exampleNameInline(primerExample(base)!.row)}`);
  });

  it("claims no more than internal cross-validation and labels its drawings", () => {
    const steps = primerSteps(base);
    const all = steps.map((s) => `${s.lead} ${s.caption} ${text(s.body)} ${text(s.figure)}`).join(" ");
    expect(all).toContain("Internal cross-validation only");
    expect(all).toContain("No external, prospective or clinical validation has been done");
    expect(all).toContain("Choosing a treatment for a patient");
    expect(text(steps[1]!.figure)).toContain("Drawing, not data");
    expect(steps[2]!.caption).toContain("the layout is a drawing");
    expect(all).not.toMatch(/clinically validated|proves|diagnos(e|tic accuracy)|let'?s|we'll|in this section/i);
  });
});

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

const fakeNav = () => {
  const calls: { open: string[]; go: Partial<Route>[] } = { open: [], go: [] };
  const nav: Navigator = { route: { view: "overview", species: null, dataset: null }, go: (n) => calls.go.push(n), openDataset: (id) => calls.open.push(id), closeDataset: () => {} };
  return { nav, calls };
};

const motion = (reduced: boolean) =>
  vi.stubGlobal("matchMedia", (query: string) => ({ matches: reduced && query.includes("prefers-reduced-motion"), media: query, addEventListener: () => {}, removeEventListener: () => {} }));

describe("the primer on Overview", () => {
  it("sits above the numbered findings with its five steps, a real report link, and a skip link that moves focus", async () => {
    motion(true);
    const { nav, calls } = fakeNav();
    const host = document.createElement("div");
    document.body.append(host);
    root = createRoot(host);
    await act(async () => root!.render(createElement(OverviewView, { manifest: base, nav })));
    const primer = host.querySelector(`#${PRIMER_ANCHOR}`)!;
    expect(primer.compareDocumentPosition(host.querySelector("#corpus")!) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
    expect([...primer.querySelectorAll("h3")].map((h) => h.textContent)).toEqual(["A dataset", "A genome feature", "Two tests", "Reading one result", "What the numbers support"]);
    // The primer adds no tab stop beyond its two links.
    expect([...primer.querySelectorAll("a, button, [tabindex='0']")].map((e) => e.textContent?.trim())).toEqual(["Skip to the findings↓", exampleName(primerExample(base)!.row)]);
    const example = [...primer.querySelectorAll<HTMLAnchorElement>("a")][1]!;
    await act(async () => example.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, button: 0 })));
    expect(calls.open).toEqual([primerExample(base)!.row.id]);
    const skip = primer.querySelector<HTMLAnchorElement>('a[href="#corpus"]')!;
    await act(async () => skip.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, button: 0 })));
    expect(document.activeElement).toBe(host.querySelector("#corpus-title"));
    // No figure animates: the primer is static on the page.
    expect(primer.querySelector('[class*="motion-"], [class*="animate-"]')).toBeNull();
  });
});

describe("the primer in Present", () => {
  const open = async () => {
    vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path))));
    const host = document.createElement("div");
    document.body.append(host);
    root = createRoot(host);
    await act(async () => root!.render(createElement(PresentMode, { manifest: base, onClose: () => {} })));
    const dialog = document.querySelector<HTMLDialogElement>("dialog")!;
    const status = () => dialog.querySelector('[role="status"]')?.textContent ?? "";
    const key = (k: string) => act(async () => dialog.querySelector("main")!.dispatchEvent(new KeyboardEvent("keydown", { key: k, bubbles: true })));
    return { dialog, status, key };
  };

  it("follows the title with one slide per step, each with its headline, figure and caption in view under reduced motion", async () => {
    motion(true);
    const { dialog, status, key } = await open();
    const steps = primerSteps(base);
    for (const [i, s] of steps.entries()) {
      await key("ArrowRight");
      expect(status()).toContain(`Slide ${i + 2} of `);
      expect(status()).toContain(s.title);
      const section = dialog.querySelector("main section")!;
      expect(section.querySelector("h2")?.textContent).toBe(s.title);
      expect(section.textContent).toContain(s.lead);
      expect(section.querySelector("[data-caption]")?.textContent).toBe(s.caption);
      if (s.figure !== null) expect(section.querySelector("[data-primer-figure]")).not.toBeNull();
      // Reduced motion: every beat is in at once, nothing waits to appear.
      expect(dialog.querySelectorAll("[data-beats] span").length === 0 || [...dialog.querySelectorAll("[data-beats] span")].every((t) => t.className.includes("bg-blue"))).toBe(true);
    }
    await key("ArrowRight");
    expect(status()).toContain("Data");
    await key("Home");
    expect(status()).toContain("Slide 1 of ");
  });

  it("keeps the deck title and every outline label whole, wrapping on a phone instead of cutting off", async () => {
    motion(true);
    const { dialog, key } = await open();
    const title = dialog.querySelector("header [data-deck-title]")!;
    expect(title.textContent?.trim()).toBe("What GRM found · Genomic Resistance Mapping");
    expect(title.className).not.toMatch(/\btruncate\b|text-ellipsis|whitespace-nowrap/);
    expect(title.className).toContain("min-w-0");
    const outline = dialog.querySelector<HTMLButtonElement>('footer [role="combobox"]')!;
    const label = outline.querySelector("[data-dropdown-label]")!;
    expect(label.className).not.toMatch(/\btruncate\b/);
    expect(label.className).toContain("whitespace-normal");
    // A wrapped label grows the control; it never drops under 28px.
    expect(outline.className).toContain("min-h-7");
    expect(outline.className).not.toMatch(/(^|\s)h-7(\s|$)/);
    for (let i = 0; i < 4; i++) await key("ArrowRight");
    expect(label.textContent).toBe("5 · Reading one result");
    expect(outline.getAttribute("aria-label")).toBe("Go to slide: 5 · Reading one result");
    // Phone footer: the outline on its own first row, Previous and Next under
    // it at the two ends; one row from sm up.
    const [prev, next] = [...dialog.querySelectorAll<HTMLButtonElement>("footer button")].filter((b) => b.textContent === "Previous" || b.textContent === "Next");
    expect(prev!.className).toMatch(/\brow-start-2\b.*\bsm:row-start-1\b/);
    expect(next!.className).toMatch(/\bcol-start-3\b.*\brow-start-2\b.*\bsm:row-start-1\b/);
    expect(outline.closest(".col-span-3")?.className).toMatch(/\brow-start-1\b.*\bsm:col-span-1\b/);
  });

  it("sets every fixed primer figure one gap under its answer instead of centring it in the room left", async () => {
    motion(true);
    const { dialog, key } = await open();
    for (let i = 0; i < 5; i++) {
      await key("ArrowRight");
      const box = dialog.querySelector("main section [data-fixed-box]")!;
      expect(box.className).toContain("justify-start");
      expect(box.className).toContain("pt-6");
      expect(box.className).not.toContain("justify-center");
    }
    // The measuring copies use the same box, so page counts match the slide.
    expect([...dialog.querySelectorAll("main [inert] [data-fixed-box]")].every((b) => b.className === dialog.querySelector("main section [data-fixed-box]")!.className)).toBe(true);
  });

  it("stages each primer slide's beats with motion on, and an arrow lands them before moving on", async () => {
    motion(false);
    const { status, key } = await open();
    await key("ArrowRight");
    await key("ArrowRight");
    expect(status()).toContain("Slide 2 of ");
    await key("ArrowRight");
    // The dataset slide was still revealing: the arrow completed it.
    expect(status()).toContain("Slide 2 of ");
    await key("ArrowRight");
    expect(status()).toContain("Slide 3 of ");
  });
});
