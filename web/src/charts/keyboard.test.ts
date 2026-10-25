// @vitest-environment happy-dom
// Keyboard access to the two canvas figures: each is one tab stop whose
// active descendant carries the values a hover shows.

import { act, createElement, type ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import type { GenomeRow, SpeciesReport } from "../data/schema";
import { ThemeProvider } from "../theme/ThemeProvider";
import { DotStrip } from "./DotStrip";
import { stepCursor } from "./kit";
import { ScatterMap, type MapPoint } from "./ScatterMap";
import { StructureFigure } from "./Structure";

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

const widget = (host: HTMLElement, role: string): HTMLElement => {
  const el = host.querySelector<HTMLElement>(`[role="${role}"]`);
  if (el === null) throw new Error(`no ${role}`);
  return el;
};

const press = (el: HTMLElement, key: string): KeyboardEvent => {
  const e = new KeyboardEvent("keydown", { key, bubbles: true, cancelable: true });
  act(() => {
    el.dispatchEvent(e);
  });
  return e;
};

/** What assistive technology reads for the focused widget's active mark. */
const announced = (el: HTMLElement): string | null => {
  const id = el.getAttribute("aria-activedescendant");
  return id === null ? null : (document.getElementById(id)?.textContent ?? null);
};

const tooltip = (host: HTMLElement): string | null => host.querySelector('[role="tooltip"]')?.textContent ?? null;

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

  it("is one named tab stop whose arrows walk the genomes and expose each tip", () => {
    const host = mount(createElement(ScatterMap, { points, label: "Map of 3 genomes" }));
    const map = widget(host, "listbox");
    expect(map.tabIndex).toBe(0);
    expect(map.getAttribute("aria-label")).toBe("Map of 3 genomes");

    act(() => map.focus());
    expect(announced(map)).toBe("B, Amikacin: susceptible. Lineage b");
    expect(tooltip(host)).toContain("Amikacin: susceptible.");

    expect(press(map, "ArrowRight").defaultPrevented).toBe(true);
    expect(announced(map)).toBe("C, Amikacin: resistant. Lineage c");
    expect(tooltip(host)).toContain("Lineage c");
    press(map, "End");
    expect(announced(map)).toBe("A, Amikacin: resistant. Lineage a");
    expect(host.querySelector('[role="option"]')?.getAttribute("aria-posinset")).toBe("3");
    press(map, "Home");
    expect(announced(map)).toBe("B, Amikacin: susceptible. Lineage b");
  });

  it("lets Tab leave, and Escape hide the tip without losing the place", () => {
    const host = mount(createElement(ScatterMap, { points, label: "Map of 3 genomes" }));
    const map = widget(host, "listbox");
    act(() => map.focus());
    press(map, "ArrowRight");

    expect(press(map, "Tab").defaultPrevented).toBe(false);
    expect(press(map, "Escape").defaultPrevented).toBe(true);
    expect(tooltip(host)).toBeNull();
    // A second Escape is not the chart's, so a dialog behind it can close.
    expect(press(map, "Escape").defaultPrevented).toBe(false);
    press(map, "ArrowRight");
    expect(announced(map)).toBe("A, Amikacin: resistant. Lineage a");
    expect(tooltip(host)).toContain("Lineage a");

    act(() => map.blur());
    expect(tooltip(host)).toBeNull();
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
// Rows show colistin, then amikacin.
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

describe("StructureFigure keyboard", () => {
  it("starts on the outlined drug's row and moves by genome and by row", () => {
    const host = mount(createElement(StructureFigure, { report, drugOrder: [1, 0], selected: 0 }));
    const grid = widget(host, "grid");
    expect(grid.tabIndex).toBe(0);
    expect(grid.getAttribute("aria-label")).toBe("Dendrogram of 3 genomes above their resistance calls for 2 drugs, in the same order.");

    act(() => grid.focus());
    expect(announced(grid)).toBe("G2, Amikacin: not tested. Lineage 2");
    expect(tooltip(host)).toContain("Amikacin: not tested.");

    press(grid, "ArrowRight");
    expect(announced(grid)).toBe("G0, Amikacin: resistant. Lineage 1");
    expect(host.querySelector('[role="gridcell"]')?.getAttribute("aria-colindex")).toBe("2");

    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("G0, Colistin: susceptible. Lineage 1");
    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("G0, Resistant to 1 of 2 drugs tested. Lineage 1");
    expect(host.querySelector('[role="row"]')?.getAttribute("aria-rowindex")).toBe("1");
    press(grid, "ArrowUp");
    expect(announced(grid)).toBe("G0, Resistant to 1 of 2 drugs tested. Lineage 1");

    press(grid, "End");
    expect(announced(grid)).toBe("G1, Resistant to 0 of 1 drugs tested. Lineage 1");
    expect(press(grid, "Tab").defaultPrevented).toBe(false);
    expect(press(grid, "Escape").defaultPrevented).toBe(true);
    expect(tooltip(host)).toBeNull();
  });
});


describe("distribution chart keyboard", () => {
  it("walks marks by value, opens the selected dataset, and dismisses the tooltip before the dialog", () => {
    const onOpen = vi.fn();
    const host = mount(createElement(DotStrip, {
      label: "Test distribution", data: [{ id: "high", value: 0.9, tip: [{ term: "High", text: "90%" }] }, { id: "low", value: 0.1, tip: [{ term: "Low", text: "10%" }] }],
      domain: [0, 1], format: String, axisTitle: "Agreement", refs: [], regions: [], named: [], onOpen,
    }));
    const chart = widget(host, "listbox");
    expect(chart.tabIndex).toBe(0);
    act(() => chart.focus());
    press(chart, "Home");
    expect(announced(chart)).toContain("Low, 10%");
    press(chart, "End");
    expect(announced(chart)).toContain("High, 90%");
    press(chart, "Enter");
    expect(onOpen).toHaveBeenCalledWith("high");
    expect(press(chart, "Escape").defaultPrevented).toBe(true);
    expect(tooltip(host)).toBeNull();
    expect(press(chart, "Escape").defaultPrevented).toBe(false);
    expect(press(chart, "Tab").defaultPrevented).toBe(false);
  });
});
