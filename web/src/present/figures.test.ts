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

const open = async (manifest: Manifest = base) => {
  vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path))));
  const host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  const onClose = vi.fn();
  await act(async () => root!.render(createElement(PresentMode, { manifest, onClose })));
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

it("moves only from its own controls: a click on a slide or its caption does nothing", async () => {
  motion(true);
  const { dialog, status, right } = await open();
  await right();
  expect(status()).toContain("Slide 2 of 7");
  const main = dialog.querySelector("main")!;
  await act(async () => main.querySelector("[data-caption]")!.dispatchEvent(new MouseEvent("click", { bubbles: true })));
  await act(async () => main.dispatchEvent(new MouseEvent("click", { bubbles: true })));
  expect(status()).toContain("Slide 2 of 7");
});

it("states the report's scope on the title slide as plain words, with no chip to hover", async () => {
  motion(true);
  // The fixture's rows carry no audit, so they read as provisional.
  const { dialog } = await open();
  const title = dialog.querySelector("main section")!;
  expect(title.querySelector("strong")?.textContent).toBe("Provisional metrics");
  expect(title.querySelector('[tabindex="0"]')).toBeNull();
  act(() => root?.unmount());
  document.body.innerHTML = "";
  const passing: Manifest = { ...base, evaluation_revision: 7, datasets: base.datasets.map((d) => ({ ...d, evaluation: { revision: 7, valid: true, reasons: [] } })) };
  const again = await open(passing);
  expect(again.dialog.querySelector("main section")?.textContent).toContain("Internal cross-validation only");
});

it("shows the rule's sequence inline, never in a second dialog", async () => {
  motion(true);
  const { dialog, status } = await open();
  await act(async () => dialog.querySelector("main")!.dispatchEvent(new KeyboardEvent("keydown", { key: "End", bubbles: true })));
  expect(status()).toContain("Slide 7 of 7: One interpretable rule");
  expect(document.querySelectorAll("dialog")).toHaveLength(1);
  expect([...dialog.querySelectorAll("button")].map((b) => b.textContent)).not.toContain("Show sequence");
  const pick = ruleDataset(base)!;
  const report = fixtureFile(`/data/datasets/${pick.row.id}.json`) as { features: { feature: number; locus: { context: string } | null }[]; models: { scm: { rules: { feature: number }[] } } };
  const context = report.features.find((f) => f.feature === report.models.scm.rules[0]!.feature)?.locus?.context ?? "";
  expect(dialog.querySelector("main section .seq") !== null).toBe(context !== "");
});

/** Gives every element the same laid-out box, as a browser would at one
 *  size: happy-dom lays nothing out. `tall` makes a fixed slide's content
 *  taller than its box. */
const layout = ({ height, tall = false }: { height: number; tall?: boolean }) => {
  const define = (key: string, get: (this: HTMLElement) => number) => Object.defineProperty(HTMLElement.prototype, key, { configurable: true, get });
  define("clientHeight", () => height);
  define("clientWidth", () => 600);
  define("offsetHeight", function () {
    return tall && this.parentElement?.classList.contains("justify-center") === true ? height * 3 : 40;
  });
  return () => ["clientHeight", "clientWidth", "offsetHeight"].forEach((key) => delete (HTMLElement.prototype as unknown as Record<string, unknown>)[key]);
};

describe("continuation slides", () => {
  let restore: (() => void) | null = null;
  afterEach(() => {
    restore?.();
    restore = null;
  });

  it("continues a dense panel figure onto further slides that repeat its caption, with every row reachable", async () => {
    motion(true);
    restore = layout({ height: 300 });
    const panel = base.panels.find((p) => p.panel === "carbapenemase")!;
    const m: Manifest = { ...base, panels: Array.from({ length: 15 }, (_, k) => ({ ...panel, dataset: `p${k}`, antibiotic: `drug${k}` })) };
    vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path))));
    const host = document.createElement("div");
    document.body.append(host);
    root = createRoot(host);
    await act(async () => root!.render(createElement(PresentMode, { manifest: m, onClose: () => {} })));
    const dialog = document.querySelector<HTMLDialogElement>("dialog")!;
    const status = () => dialog.querySelector('[role="status"]')?.textContent ?? "";
    // 300px leaves 202px for rows; 15 rows at 18px need 270, so they split
    // into pages of at most 9 rows: two pages of 8 and 7.
    const caption = panelsFinding(m, { page: false }).caption;
    const seen: string[] = [];
    for (const [n, first, last] of [
      [5, 1, 8],
      [6, 9, 15],
    ] as const) {
      while (!status().startsWith(`Slide ${n} `)) await act(async () => dialog.querySelector("main")!.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowRight", bubbles: true })));
      expect(status()).toBe(`Slide ${n} of 8: PCR panel (${n - 4} of 2)`);
      const section = dialog.querySelector("main section")!;
      const ids = [...section.querySelectorAll("[data-dataset]")].map((g) => g.getAttribute("data-dataset")!);
      expect(ids).toHaveLength(last - first + 1);
      seen.push(...ids);
      expect(section.querySelector("[data-continued]")?.textContent).toBe(`Datasets ${first}–${last} of 15`);
      expect(section.querySelector("[data-caption]")?.textContent).toBe(caption);
    }
    expect(seen.sort()).toEqual(Array.from({ length: 15 }, (_, k) => `p${k}`).sort());
    const selector = dialog.querySelector<HTMLButtonElement>('[role="combobox"]')!;
    await act(async () => selector.click());
    // The outline's options, not the open chart's keyboard option.
    const labels = [...document.querySelectorAll('[role="option"]')].map((o) => o.textContent ?? "").filter((t) => /^\d+ · /.test(t));
    expect(labels).toContain("5 · PCR panel (1 of 2)");
    expect(labels).toContain("6 · PCR panel (2 of 2)");
    expect(labels).toHaveLength(8);
  });

  it("moves the rule's sequence onto its own slide only when the two do not fit together", async () => {
    motion(true);
    restore = layout({ height: 300, tall: true });
    const pick = ruleDataset(base)!;
    const report = fixtureFile(`/data/datasets/${pick.row.id}.json`) as { features: { feature: number; locus: { context: string } | null }[]; models: { scm: { rules: { feature: number }[] } } };
    const context = report.features.find((f) => f.feature === report.models.scm.rules[0]!.feature)?.locus?.context ?? "";
    const { dialog, status } = await open();
    await act(async () => dialog.querySelector("main")!.dispatchEvent(new KeyboardEvent("keydown", { key: "End", bubbles: true })));
    if (context === "") {
      expect(status()).toContain("Slide 7 of 7: One interpretable rule");
      return;
    }
    expect(status()).toBe("Slide 8 of 8: One interpretable rule (2 of 2)");
    expect(dialog.querySelector("main section .seq")).not.toBeNull();
    expect(dialog.querySelector("main section")?.textContent).not.toContain("Resistant isolates carrying it");
    await act(async () => dialog.querySelector("main")!.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowLeft", bubbles: true })));
    expect(status()).toBe("Slide 7 of 8: One interpretable rule (1 of 2)");
    expect(dialog.querySelector("main section .seq")).toBeNull();
    expect(dialog.querySelector("main section")?.textContent).toContain("Resistant isolates carrying it");
    expect(dialog.querySelector("main section")?.textContent).toContain("The rule's sequence follows on the next slide.");
  });
});
