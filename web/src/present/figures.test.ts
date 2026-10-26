// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { CatalogRow, Manifest } from "../data/schema";
import { FINDINGS, catalogFinding, lineageFinding, micFinding, panelsFinding } from "../views/findings";
import { PresentMode } from "./PresentMode";
import { ruleDataset } from "./figures";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });
const base = parseManifest(fixtureFile("/data/index.json"));
const row = (id: string, overrides: Partial<CatalogRow>): CatalogRow => ({ ...base.datasets[0]!, id, antibiotic: id, ...overrides });
const figure = (finding: ReturnType<typeof catalogFinding>) => {
  const node = typeof finding.figure === "function" ? finding.figure({}) : finding.figure;
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(node);
  return el;
};

describe("shared distribution figures", () => {
  it("draws every MIC dataset, including zero, and omits missing measurements", () => {
    const m = { ...base, datasets: [row("zero", { mic_essential_agreement: 0 }), ...Array.from({ length: 31 }, (_, k) => row(`d${k}`, { mic_essential_agreement: (60 + k) / 100 })), row("missing", { mic_essential_agreement: null })] };
    const f = micFinding(m, { page: false });
    const html = figure(f);
    expect(html.querySelector('[role="listbox"]')?.getAttribute("aria-label")).toContain("32 datasets");
    expect(html.querySelectorAll("svg circle[data-dataset]")).toHaveLength(32);
    expect(f.result).toContain("1 of 32 datasets");
    expect(f.caption).toContain("out of fold");
  });

  it("draws paired datasets and labels their paired medians", () => {
    const m = { ...base, datasets: [row("a", { best_random: 0.9, best_lineage: 0.8 }), row("b", { best_random: 0.8, best_lineage: 0.78 }), row("unpaired", { best_random: 0.1, best_lineage: null })] };
    const f = lineageFinding(m, { page: false });
    expect(figure(f).querySelectorAll("svg circle[data-dataset]")).toHaveLength(2);
    expect(f.caption).toContain("85.0% random, 79.0% lineage");
  });

  it("shows catalog region totals and explicit near-uniform thresholds", () => {
    const f = catalogFinding(base, { page: false });
    const eligible = base.datasets.filter((d) => d.catalog !== null && d.best_lineage !== null);
    const html = figure(f);
    expect(html.querySelectorAll("svg g[data-dataset]")).toHaveLength(eligible.length);
    expect(f.result).toContain("within 2 points or ahead");
    expect(html.textContent).toContain("≥90% of susceptible isolates resistant (specificity ≤10%)");
    expect(html.textContent).toContain("≥90% of resistant isolates susceptible (sensitivity ≤10%)");
  });

  it("draws every lead-panel dataset without changing its selection", () => {
    const panel = base.panels.find((p) => p.panel === "carbapenemase")!;
    const m: Manifest = { ...base, panels: Array.from({ length: 15 }, (_, k) => ({ ...panel, dataset: `p${k}`, antibiotic: `drug${k}` })) };
    const f = panelsFinding(m, { page: false });
    expect(figure(f).querySelectorAll("svg g[data-dataset]")).toHaveLength(15);
    expect(f.caption).toContain("one of 15 datasets");
  });

  it("keeps the rule showcase's random-score ranking and isolate-count tiebreak", () => {
    const candidates = base.datasets.flatMap((row) => {
      const ba = row.learners.find((l) => l.learner === "scm")?.random?.balanced_accuracy;
      return ba == null || row.scm_rules.length === 0 ? [] : [{ id: row.id, ba, n: row.n }];
    });
    const expected = candidates.sort((a, b) => b.ba - a.ba || b.n - a.n)[0];
    expect(ruleDataset(base)?.row.id).toBe(expected?.id);
    expect(ruleDataset({ ...base, datasets: [] })).toBeNull();
  });
});

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

/** matchMedia that answers the reduced-motion query as given. */
const motion = (reduced: boolean) =>
  vi.stubGlobal("matchMedia", (query: string) => ({
    matches: reduced && query.includes("prefers-reduced-motion"),
    media: query,
    addEventListener: () => {},
    removeEventListener: () => {},
  }));

const open = async () => {
  vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path))));
  const host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  const onClose = vi.fn();
  await act(async () => root!.render(createElement(PresentMode, { manifest: base, onClose })));
  const dialog = document.querySelector<HTMLDialogElement>("dialog")!;
  const status = () => dialog.querySelector('[role="status"]')?.textContent ?? "";
  const next = [...dialog.querySelectorAll("button")].find((b) => b.textContent === "Next")!;
  const right = () => act(async () => next.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowRight", bubbles: true })));
  return { dialog, status, right, onClose };
};

it("opens the deck on its title slide, then the six finding slides, with a named outline and no deck chooser", async () => {
  motion(true);
  const { dialog, status, right, onClose } = await open();
  const titles = ["Genomic Resistance Mapping", ...FINDINGS.map((f) => f(base, { page: false }).title), "One interpretable rule"];
  expect(dialog.textContent).not.toContain("Choose a deck");
  const selector = dialog.querySelector<HTMLButtonElement>('[role="combobox"]')!;
  expect(selector.getAttribute("aria-label")).toBe("Go to slide: 1 · Genomic Resistance Mapping");
  await act(async () => selector.click());
  expect([...document.querySelectorAll('[role="option"]')].map((o) => o.textContent)).toEqual(titles.map((t, i) => `${i + 1} · ${t}`));
  await act(async () => selector.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true })));
  // With reduced motion every beat is already in, so each arrow moves a slide.
  for (let i = 1; i < titles.length; i++) {
    expect(status()).toContain(`Slide ${i} of ${titles.length}`);
    await right();
  }
  expect(status()).toContain(`Slide ${titles.length} of ${titles.length}: One interpretable rule`);
  // The slide control keeps its own arrow keys.
  act(() => selector.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowLeft", bubbles: true })));
  expect(status()).toContain(`Slide ${titles.length} of ${titles.length}`);
  act(() => dialog.dispatchEvent(new Event("cancel", { cancelable: true })));
  expect(onClose).toHaveBeenCalledOnce();
});

it("completes a slide's staged beats before the next arrow moves on", async () => {
  motion(false);
  const { status, right } = await open();
  await right();
  // The title slide was still revealing: the first arrow lands its beats.
  expect(status()).toContain("Slide 1 of 7");
  await right();
  expect(status()).toContain("Slide 2 of 7");
});
