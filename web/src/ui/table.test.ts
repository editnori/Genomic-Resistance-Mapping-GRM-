// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { DataTable, tableGeometry, type Column, type SortState } from "./DataTable";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

interface Row {
  id: string;
  name: string;
  score: number;
  note: string;
}

const ROWS: Row[] = Array.from({ length: 30 }, (_, i) => ({ id: `r${i}`, name: `Row ${String(i).padStart(2, "0")}`, score: (i * 37) % 100, note: i % 2 === 0 ? "even" : "odd" }));

const COLUMNS: Column<Row>[] = [
  { id: "name", header: "Name", width: "10rem", sortValue: (r) => r.name, firstDir: "asc", cell: (r) => r.name },
  { id: "score", header: "Score", unit: "points", width: "5rem", align: "right", sortValue: (r) => r.score, cell: (r) => String(r.score) },
  { id: "note", header: "Note", minWidth: "8rem", cell: (r) => createElement("button", { type: "button" }, r.note) },
];

let host: HTMLDivElement;
let root: Root;

beforeEach(() => {
  // Reduced motion: a sort lands without its slide.
  window.matchMedia = ((q: string) =>
    ({ matches: q.includes("prefers-reduced-motion"), media: q, onchange: null, addEventListener: () => {}, removeEventListener: () => {}, addListener: () => {}, removeListener: () => {}, dispatchEvent: () => false }) as MediaQueryList) as typeof window.matchMedia;
  host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
});

afterEach(() => {
  act(() => root.unmount());
  host.remove();
  window.getSelection()?.removeAllRanges();
});

type Props = Partial<Parameters<typeof DataTable<Row>>[0]>;

const draw = (props: Props = {}) =>
  act(() => root.render(createElement(DataTable<Row>, { label: "Rows", rows: ROWS, columns: COLUMNS, rowKey: (r) => r.id, onOpen: () => {}, empty: null, ...props })));

const button = (name: string): HTMLButtonElement => {
  const b = [...host.querySelectorAll("button")].find((e) => e.textContent?.trim() === name || e.textContent?.includes(name));
  if (b === undefined) throw new Error(`no button ${name}`);
  return b;
};

const header = (name: string): HTMLElement => [...host.querySelectorAll("th")].find((t) => t.textContent?.includes(name))!;

describe("dense table geometry", () => {
  it("lays out fixed from declared widths, the one column without a width resting on its floor", () => {
    expect(tableGeometry([{}, {}])).toBeNull();
    expect(tableGeometry([{ width: "10rem" }, { width: "5rem" }, { minWidth: "8rem" }])).toEqual({ cols: ["10rem", "5rem", null], minWidth: "calc(10rem + 5rem + 8rem)" });
  });

  it("holds every column width through a sort, Show all and a filtered row set", () => {
    draw();
    const widths = () => {
      const table = host.querySelector("table")!;
      return { cols: [...host.querySelectorAll("col")].map((c) => c.style.width), layout: table.style.tableLayout, floor: table.style.minWidth };
    };
    const arrival = widths();
    expect(arrival).toEqual({ cols: ["10rem", "5rem", ""], layout: "fixed", floor: "calc(10rem + 5rem + 8rem)" });

    act(() => header("Score").querySelector("button")!.click());
    expect(widths()).toEqual(arrival);
    act(() => button("Show all 30").click());
    expect(host.querySelectorAll("tbody tr")).toHaveLength(30);
    expect(widths()).toEqual(arrival);
    draw({ rows: ROWS.slice(0, 3) });
    expect(widths()).toEqual(arrival);
  });

  it("keeps a table without declared widths in its automatic layout", () => {
    draw({ columns: COLUMNS.map((c) => ({ ...c, width: undefined, minWidth: undefined })) });
    expect(host.querySelector("colgroup")).toBeNull();
    expect(host.querySelector("table")!.style.tableLayout).toBe("");
  });

  it("sets the sort arrow in the label's line, before a right-aligned label and after a left-aligned one", () => {
    draw();
    const side = (name: string) => {
      const th = header(name);
      const label = th.querySelector("[data-th-label]")!;
      const glyph = th.querySelector("[data-sort-glyph]")!;
      expect(glyph.parentElement).toBe(label.parentElement);
      return label.compareDocumentPosition(glyph) & Node.DOCUMENT_POSITION_PRECEDING ? "before" : "after";
    };
    expect(side("Score")).toBe("before");
    expect(side("Name")).toBe("after");
    // The unit line sits under the label, outside the arrow's line.
    expect(header("Score").querySelector("[data-th-label]")!.parentElement!.textContent).toBe("Score");
  });
});

describe("row activation", () => {
  it("opens a row from a plain cell, never from a control in it or from a text selection", () => {
    const open = vi.fn();
    draw({ onOpen: open });
    const row = host.querySelector("tbody tr")!;
    const score = row.children[1] as HTMLElement;

    act(() => score.click());
    expect(open).toHaveBeenCalledTimes(1);

    act(() => row.querySelector<HTMLButtonElement>("td:last-child button")!.click());
    expect(open).toHaveBeenCalledTimes(1);

    window.getSelection()!.selectAllChildren(score);
    act(() => score.click());
    expect(open).toHaveBeenCalledTimes(1);
    window.getSelection()!.removeAllRanges();

    act(() => row.querySelector<HTMLButtonElement>("[data-row-open]")!.click());
    expect(open).toHaveBeenCalledTimes(2);
    expect(row.querySelector<HTMLElement>("[data-row-open]")!.dataset.rowKey).toBe(row.getAttribute("data-row-key"));
  });

  it("starts from a kept sort and Show all, and reports each change", () => {
    const sorts: SortState[] = [];
    const alls: boolean[] = [];
    draw({ initialSort: { id: "score", dir: "asc" }, initialShowAll: true, onSortChange: (s) => sorts.push(s), onShowAllChange: (a) => alls.push(a) });
    expect(host.querySelectorAll("tbody tr")).toHaveLength(30);
    expect(header("Score").getAttribute("aria-sort")).toBe("ascending");

    act(() => header("Score").querySelector("button")!.click());
    expect(sorts).toEqual([{ id: "score", dir: "desc" }]);
    act(() => button("Show 25").click());
    expect(alls).toEqual([false]);
    expect(host.querySelectorAll("tbody tr")).toHaveLength(25);
  });
});
