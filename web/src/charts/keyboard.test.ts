// @vitest-environment happy-dom
// Keyboard, tap and pin behaviour of the charts: each is one tab stop whose
// active descendant, one id per mark, carries the mark's words. Click, tap,
// Enter or Space pins a mark until Escape; a hover only previews.

import { act, createElement, useMemo, useState, type ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import type { GenomeRow, Merge, SpeciesReport } from "../data/schema";
import { ThemeProvider } from "../theme/ThemeProvider";
import { DotStrip } from "./DotStrip";
import { stepCursor } from "./kit";
import { ScatterMap, type MapPoint } from "./ScatterMap";
import { StructureExplorer, structureModel, type StructurePick, type StructureView } from "./Structure";
import { CrossPanel } from "./CrossMatrix";
import { crossPairs, pairSummary } from "./crossPairs";
import { YearColumns } from "./YearColumns";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let root: Root | null = null;

const mount = (node: ReactNode): HTMLElement => {
  const host = document.createElement("div");
  document.body.append(host);
  const r = createRoot(host);
  root = r;
  act(() => r.render(createElement(ThemeProvider, { children: node })));
  return host;
};

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
});

const widget = (host: HTMLElement, role: string, name?: RegExp): HTMLElement => {
  const el = [...host.querySelectorAll<HTMLElement>(`[role="${role}"]`)].find((e) => name === undefined || name.test(e.getAttribute("aria-label") ?? ""));
  if (el === undefined) throw new Error(`no ${role}`);
  return el;
};

const press = (el: HTMLElement, key: string): KeyboardEvent => {
  const e = new KeyboardEvent("keydown", { key, bubbles: true, cancelable: true });
  act(() => {
    el.dispatchEvent(e);
  });
  return e;
};

const tap = (el: HTMLElement, x: number, y: number) =>
  act(() => {
    el.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, clientX: x, clientY: y }));
  });

/** What assistive technology reads for the focused widget's active mark. */
const announced = (el: HTMLElement): string | null => {
  const id = el.getAttribute("aria-activedescendant");
  return id === null ? null : (document.getElementById(id)?.textContent ?? null);
};

const tooltips = (host: HTMLElement): string[] => [...host.querySelectorAll('[role="tooltip"]')].map((t) => t.textContent ?? "");

const button = (host: HTMLElement, text: string): HTMLButtonElement | undefined => [...host.querySelectorAll("button")].find((b) => b.textContent === text);

describe("stepCursor", () => {
  it("steps by one, by a tenth of the marks, and to either end", () => {
    expect(stepCursor("ArrowRight", -1, 2000)).toBe(0);
    expect(stepCursor("ArrowRight", 4, 2000)).toBe(5);
    expect(stepCursor("ArrowLeft", 0, 2000)).toBe(0);
    expect(stepCursor("PageDown", 0, 2000)).toBe(200);
    expect(stepCursor("PageUp", 150, 2000)).toBe(0);
    expect(stepCursor("End", -1, 2000)).toBe(1999);
    expect(stepCursor("Home", 1999, 2000)).toBe(0);
  });

  it("leaves Tab and other keys alone", () => {
    expect(stepCursor("Tab", 3, 10)).toBeNull();
    expect(stepCursor("Enter", 3, 10)).toBeNull();
    expect(stepCursor("ArrowRight", -1, 0)).toBeNull();
  });
});

const point = (id: string, x: number, y: number, call: string): MapPoint => ({
  x,
  y,
  layer: "mark",
  tone: "r",
  tip: () => [{ term: id, text: `Amikacin: ${call}.` }, { source: `Lineage ${id.toLowerCase()}` }],
});

describe("ScatterMap keyboard", () => {
  // Keyboard order runs left to right: B, C, A.
  const points = [point("A", 0.3, 0, "resistant"), point("B", 0.1, 0, "susceptible"), point("C", 0.2, 0, "resistant")];

  it("is one named tab stop whose arrows walk the genomes, each under its own id, previewing only the id and value", () => {
    const host = mount(createElement(ScatterMap, { points, label: "Map of 3 genomes" }));
    const map = widget(host, "listbox");
    expect(map.tabIndex).toBe(0);
    expect(map.getAttribute("aria-label")).toBe("Map of 3 genomes");

    act(() => map.focus());
    expect(announced(map)).toBe("B, Amikacin: susceptible. Lineage b");
    expect(tooltips(host).join()).toContain("Amikacin: susceptible.");
    const first = map.getAttribute("aria-activedescendant");

    expect(press(map, "ArrowRight").defaultPrevented).toBe(true);
    expect(announced(map)).toBe("C, Amikacin: resistant. Lineage c");
    expect(map.getAttribute("aria-activedescendant")).not.toBe(first);
    expect(tooltips(host).join()).toContain("Amikacin: resistant.");
    expect(tooltips(host).join()).not.toContain("Lineage c");
    press(map, "End");
    expect(announced(map)).toBe("A, Amikacin: resistant. Lineage a");
    expect(host.querySelector('[role="option"]')?.getAttribute("aria-posinset")).toBe("3");
    press(map, "Home");
    expect(announced(map)).toBe("B, Amikacin: susceptible. Lineage b");
  });

  it("pins with Enter until Escape, then lets Escape hide the preview, then reach a dialog", () => {
    const host = mount(createElement(ScatterMap, { points, label: "Map of 3 genomes" }));
    const map = widget(host, "listbox");
    act(() => map.focus());
    press(map, "ArrowRight");
    expect(press(map, "Enter").defaultPrevented).toBe(true);
    expect(announced(map)).toBe("C, Amikacin: resistant. Lineage c Pinned.");
    // The pinned genome keeps its full details up.
    expect(tooltips(host).join()).toContain("Lineage c");
    expect(press(map, "Tab").defaultPrevented).toBe(false);

    expect(press(map, "Escape").defaultPrevented).toBe(true);
    expect(announced(map)).toBe("C, Amikacin: resistant. Lineage c");
    expect(press(map, "Escape").defaultPrevented).toBe(true);
    expect(tooltips(host)).toEqual([]);
    expect(press(map, "Escape").defaultPrevented).toBe(false);

    press(map, "ArrowRight");
    expect(tooltips(host).join()).toContain("Amikacin: resistant.");
    act(() => map.blur());
    expect(tooltips(host)).toEqual([]);
  });
});

describe("distribution chart pinning", () => {
  const data = [
    { id: "high", value: 0.9, tip: [{ term: "High", text: "90%" }] },
    { id: "low", value: 0.1, tip: [{ term: "Low", text: "10%" }] },
  ];
  const strip = (onOpen: (id: string) => void) =>
    createElement(DotStrip, { label: "Test distribution", data, domain: [0, 1], format: String, axisTitle: "Agreement", refs: [], regions: [], named: [], onOpen });

  it("walks marks by value; Enter pins one, whose strip opens it; Escape clears before a dialog sees it", () => {
    const onOpen = vi.fn();
    const host = mount(strip(onOpen));
    const chart = widget(host, "listbox");
    expect(chart.tabIndex).toBe(0);
    act(() => chart.focus());
    press(chart, "Home");
    expect(announced(chart)).toContain("Low, 10%");
    press(chart, "End");
    expect(announced(chart)).toContain("High, 90%");
    expect(button(host, "Open dataset")).toBeUndefined();
    expect(press(chart, "Enter").defaultPrevented).toBe(true);
    expect(onOpen).not.toHaveBeenCalled();
    expect(announced(chart)).toContain("Pinned.");
    act(() => button(host, "Open dataset")!.click());
    expect(onOpen).toHaveBeenCalledWith("high");

    expect(press(chart, "Escape").defaultPrevented).toBe(true);
    expect(button(host, "Open dataset")).toBeUndefined();
    expect(press(chart, "Escape").defaultPrevented).toBe(true);
    expect(press(chart, "Escape").defaultPrevented).toBe(false);
    expect(press(chart, "Tab").defaultPrevented).toBe(false);
  });

  it("pins a mark from a stationary tap, with no pointer move first", () => {
    const onOpen = vi.fn();
    const host = mount(strip(onOpen));
    const chart = widget(host, "listbox");
    // The 640px fallback frame puts 0.9 at x 560, on the midline at y 62.
    tap(chart, 560, 62);
    expect(host.textContent).toContain("90%");
    act(() => button(host, "Open dataset")!.click());
    expect(onOpen).toHaveBeenCalledWith("high");
    // A second tap on the same mark lets it go.
    tap(chart, 560, 62);
    expect(button(host, "Open dataset")).toBeUndefined();
  });
});

describe("year columns pinning", () => {
  it("walks years by key and pins a tapped column anywhere along its height", () => {
    const host = mount(createElement(YearColumns, { years: [["2010", 5], ["2012", 7], ["unknown", 2]] }));
    const chart = widget(host, "listbox");
    act(() => chart.focus());
    expect(announced(chart)).toBe("2010, 5 labeled genomes");
    press(chart, "End");
    expect(announced(chart)).toBe("2012, 7 labeled genomes");
    act(() => chart.blur());
    // The 640px fallback frame spans 2010-2012 in three columns; 2012's
    // column centres at x 535. A tap low in the column still pins it.
    tap(chart, 535, 160);
    expect(tooltips(host).join()).toContain("7 labeled genomes");
    expect(host.textContent).toContain("2 genomes carry a year that is not a four-digit year");
  });
});

const genome = (id: string, lineage: number, profile: string): GenomeRow => ({
  id,
  lineage,
  x: 0,
  y: 0,
  nearest: 0.001,
  country: null,
  year: null,
  host: null,
  mlst: null,
  length: 0,
  contigs: 0,
  n50: 0,
  groups: [],
  profile,
});

// Three genomes over amikacin and colistin; leaf order puts genome 2 first.
// Rows show colistin, then amikacin; both joins sit within ten times the
// threshold, so the two small lineages pool into one block.
const report: SpeciesReport = {
  schema: 1,
  species: "Test species",
  slug: "test",
  drugs: ["amikacin", "colistin"],
  drug_classes: [],
  catalog_groups: [],
  lineage_threshold: 0.01,
  n_lineages: 2,
  leaf_order: [2, 0, 1],
  merges: [
    { a: 0, b: 1, height: 0.01, size: 2 },
    { a: 2, b: 3, height: 0.05, size: 3 },
  ],
  genomes: [genome("G0", 1, "RS"), genome("G1", 1, "S."), genome("G2", 2, ".R")],
  cross: [],
};

const Harness = ({ s, view: v0 }: { s: SpeciesReport; view: StructureView }) => {
  const [pick, setPick] = useState<StructurePick>(null);
  const [view, setView] = useState(v0);
  const model = useMemo(() => structureModel(s), [s]);
  return createElement(StructureExplorer, { report: s, model, order: [1, 0], drug: 0, view, pick, onPick: setPick, onView: setView });
};

const inspectorTitle = (host: HTMLElement): string | null => host.querySelector("aside h3")?.textContent ?? null;

describe("structure blocks", () => {
  it("is one grid stop whose cells name the block's calls, each under its own id", () => {
    const host = mount(createElement(Harness, { s: report, view: "blocks" }));
    const grid = widget(host, "grid");
    expect(grid.tabIndex).toBe(0);
    expect(grid.getAttribute("aria-label")).toBe("1 lineage block of 3 genomes by 2 drugs. Arrow keys move; click or Enter to pin details; Escape clears.");

    act(() => grid.focus());
    expect(announced(grid)).toBe("2 lineages, Amikacin: 1 resistant, 0 intermediate, 1 susceptible, 1 not tested.");
    const first = grid.getAttribute("aria-activedescendant");
    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("2 lineages, Colistin: 1 resistant, 0 intermediate, 1 susceptible, 1 not tested.");
    expect(grid.getAttribute("aria-activedescendant")).not.toBe(first);
    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("Block 1 of 1, 2 lineages: 3 genomes, 1 to 3 of 3 in tree order.");
    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("Block 1 of 1, 2 lineages: 3 genomes, 1 to 3 of 3 in tree order.");
  });

  it("pins a block by Enter or a tap into the persistent inspector, and Escape clears it", () => {
    const host = mount(createElement(Harness, { s: report, view: "blocks" }));
    const grid = widget(host, "grid");
    expect(inspectorTitle(host)).toBe("Amikacin across 3 genomes");
    act(() => grid.focus());
    expect(press(grid, "Enter").defaultPrevented).toBe(true);
    expect(announced(grid)).toMatch(/Pinned\.$/);
    expect(inspectorTitle(host)).toBe("2 lineages");
    expect(host.querySelector("aside")?.textContent).toContain("3 genomes, 1–3 of 3 in tree order.");
    expect(press(grid, "Tab").defaultPrevented).toBe(false);

    expect(press(grid, "Escape").defaultPrevented).toBe(true);
    expect(inspectorTitle(host)).toBe("Amikacin across 3 genomes");
    expect(press(grid, "Escape").defaultPrevented).toBe(true);
    expect(press(grid, "Escape").defaultPrevented).toBe(false);

    // A stationary tap anywhere in the block's column pins it too.
    tap(grid, 10, 10);
    expect(inspectorTitle(host)).toBe("2 lineages");
    act(() => button(host, "Clear")!.click());
    expect(inspectorTitle(host)).toBe("Amikacin across 3 genomes");
  });
});

describe("structure genomes", () => {
  it("walks single genomes and pins one, listing its calls by state", () => {
    const host = mount(createElement(Harness, { s: report, view: "genomes" }));
    const grid = widget(host, "grid");
    expect(grid.getAttribute("aria-label")).toBe("Genomes 1 to 3 of 3 in tree order, by 2 drugs. Arrow keys move; click or Enter to pin details; Escape clears.");
    act(() => grid.focus());
    expect(announced(grid)).toBe("G2, Amikacin: not tested.");
    press(grid, "ArrowRight");
    expect(announced(grid)).toBe("G0, Amikacin: resistant.");
    press(grid, "ArrowUp");
    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("G0, lineage 1, genome 2 of 3 in tree order: 1 resistant, 0 intermediate, 1 susceptible, 0 not tested.");
    press(grid, "Enter");
    expect(inspectorTitle(host)).toBe("G0");
    const text = host.querySelector("aside")?.textContent ?? "";
    expect(text).toContain("Resistant (1)Amikacin");
    expect(text).toContain("Susceptible (1)Colistin");
    expect(text).toContain("Not tested (0)None");
    expect(press(grid, "Escape").defaultPrevented).toBe(true);
    expect(inspectorTitle(host)).toBe("Amikacin across 3 genomes");
  });

  it("moves its window by brush keys, Previous/Next and the grid, and shows where it is", () => {
    // 500 genomes in one lineage, joined in a chain.
    const n = 500;
    const merges: Merge[] = Array.from({ length: n - 1 }, (_, k) => ({ a: k === 0 ? 0 : n + k - 1, b: k + 1, height: 0.0001, size: k + 2 }));
    const big: SpeciesReport = {
      ...report,
      n_lineages: 1,
      leaf_order: Array.from({ length: n }, (_, i) => i),
      merges,
      genomes: Array.from({ length: n }, (_, i) => genome(`G${i}`, 1, i % 2 === 0 ? "RS" : "S.")),
    };
    const host = mount(createElement(Harness, { s: big, view: "genomes" }));
    const brush = widget(host, "slider");
    const position = () => host.textContent?.match(/Genomes [\d,]+–[\d,]+ of 500/)?.[0];
    // The 800px fallback frame holds 200 genomes at 4px each.
    expect(brush.getAttribute("aria-valuenow")).toBe("1");
    expect(brush.getAttribute("aria-valuemax")).toBe("301");
    expect(position()).toBe("Genomes 1–200 of 500");
    expect(button(host, "Previous")!.disabled).toBe(true);
    press(brush, "ArrowRight");
    expect(brush.getAttribute("aria-valuenow")).toBe("21");
    press(brush, "PageDown");
    expect(brush.getAttribute("aria-valuenow")).toBe("221");
    press(brush, "End");
    expect(brush.getAttribute("aria-valuenow")).toBe("301");
    expect(button(host, "Next")!.disabled).toBe(true);
    act(() => button(host, "Previous")!.click());
    expect(position()).toBe("Genomes 101–300 of 500");
    const grid = widget(host, "grid");
    act(() => grid.focus());
    press(grid, "Home");
    expect(position()).toBe("Genomes 1–200 of 500");
    press(grid, "End");
    expect(position()).toBe("Genomes 301–500 of 500");
    expect(announced(grid)).toBe("G499, Amikacin: susceptible.");
  });
});

// Four drugs: amikacin and colistin evaluated with a defined phi; tetracycline
// evaluated with no resistant genome (phi undefined); meropenem not evaluated.
const pairGenomes = [
  ...Array.from({ length: 10 }, () => "RRSS"),
  ...Array.from({ length: 10 }, () => "RSRS"),
  ...Array.from({ length: 10 }, () => "SR.S"),
  ...Array.from({ length: 10 }, () => "SS.S"),
].map((p, i) => genome(`P${i}`, 0, p));
const pairReport: SpeciesReport = {
  ...report,
  drugs: ["amikacin", "colistin", "meropenem", "tetracycline"],
  drug_classes: ["aminoglycoside", "polymyxin", "carbapenem", "tetracycline"],
  leaf_order: pairGenomes.map((_, i) => i),
  merges: [],
  genomes: pairGenomes,
  cross: [
    { a: "amikacin", b: "colistin", n: 40, rr: 10, rs: 10, sr: 10, ss: 10, phi: 0 },
    { a: "amikacin", b: "tetracycline", n: 40, rr: 0, rs: 20, sr: 0, ss: 20, phi: 0 },
    { a: "colistin", b: "tetracycline", n: 40, rr: 0, rs: 20, sr: 0, ss: 20, phi: 0 },
  ],
};

const CrossHarness = () => {
  const [selected, setSelected] = useState<string | null>(null);
  const pairs = useMemo(() => crossPairs(pairReport), []);
  return createElement(CrossPanel, { report: pairReport, order: [0, 1, 2, 3], pairs, summary: pairSummary(pairs), selected, onSelect: setSelected });
};

const pressed = (host: HTMLElement): string[] => [...host.querySelectorAll('[data-row-open][aria-pressed="true"]')].map((b) => b.textContent ?? "");

describe("cross-resistance matrix and table", () => {
  it("links a pinned cell and its table row in both directions, with exact joint counts", () => {
    const host = mount(createElement(CrossHarness));
    const grid = widget(host, "grid", /^Phi between/);
    act(() => grid.focus());
    expect(announced(grid)).toBe(
      "Amikacin and colistin. Phi 0.00 over 40 genomes called R or S for both: 10 resistant to both, 10 to amikacin only, 10 to colistin only, 10 susceptible to both.",
    );
    expect(press(grid, "Enter").defaultPrevented).toBe(true);
    expect(pressed(host)).toEqual(["Amikacin and colistin"]);
    expect(host.querySelector('[aria-live="polite"]')?.textContent).toContain("Amikacin and colistin. Phi 0.00 over 40 genomes called R or S for both");

    // A table row pins its cell: the matrix's active cell follows it.
    const row = [...host.querySelectorAll<HTMLButtonElement>("[data-row-open]")].find((b) => b.textContent === "Amikacin and tetracycline")!;
    act(() => row.click());
    expect(pressed(host)).toEqual(["Amikacin and tetracycline"]);
    expect(announced(grid)).toContain("Amikacin and tetracycline. Phi undefined: of the 40 genomes called R or S for both, none is resistant to tetracycline.");
    expect(host.querySelector('[aria-live="polite"]')?.textContent).toContain("Phi undefined");

    expect(press(grid, "Escape").defaultPrevented).toBe(true);
    expect(pressed(host)).toEqual([]);
  });

  it("lists every pair, labels undefined and not-evaluated ones, and filters to the missing pairs", () => {
    const host = mount(createElement(CrossHarness));
    const rows = () => [...host.querySelectorAll('[role="region"][aria-label="Cross-resistance pairs"] tbody tr')];
    expect(rows()).toHaveLength(6);
    const phi = (name: string) => rows().find((tr) => tr.querySelector("button")?.textContent === name)?.lastElementChild?.textContent;
    expect(phi("Amikacin and colistin")).toBe("0.00");
    expect(phi("Amikacin and tetracycline")).toBe("undefined");
    expect(phi("Amikacin and meropenem")).toBe("not evaluated");
    act(() => button(host, "Not evaluated 3")!.click());
    expect(rows().map((tr) => tr.querySelector("button")?.textContent)).toEqual(["Amikacin and meropenem", "Colistin and meropenem", "Meropenem and tetracycline"]);
  });
});
