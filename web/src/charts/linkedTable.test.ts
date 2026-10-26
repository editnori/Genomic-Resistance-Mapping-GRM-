// @vitest-environment happy-dom
// The chart table fills a long body in short tasks: a first chunk paints at
// once, every row arrives, and a pinned row is drawn wherever it sits.

import { act, createElement, type ReactNode } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it } from "vitest";
import { LinkedTable, type LinkedColumn } from "./LinkedTable";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let root: Root | null = null;

const mount = (node: ReactNode): HTMLElement => {
  const host = document.createElement("div");
  document.body.append(host);
  const r = createRoot(host);
  root = r;
  act(() => r.render(node));
  return host;
};

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
});

interface Row {
  id: string;
  v: number;
}

// 300 rows by 40 columns: chunks of 100 rows.
const rows: Row[] = Array.from({ length: 300 }, (_, i) => ({ id: `r${i}`, v: i }));
const columns: LinkedColumn<Row>[] = [
  { id: "id", header: "Row", sortValue: (r) => r.v, cell: (r) => r.id },
  ...Array.from({ length: 39 }, (_, c): LinkedColumn<Row> => ({ id: `c${c}`, header: `C${c}`, align: "right", cell: (r) => String(r.v + c) })),
];
const table = (selectedKey: string | null) => createElement(LinkedTable<Row>, { label: "Rows", rows, columns, rowKey: (r) => r.id, selectedKey, onSelect: () => {} });

const region = (host: HTMLElement) => host.querySelector<HTMLElement>('[role="region"]')!;
const count = (host: HTMLElement) => host.querySelectorAll("tbody tr").length;
const tick = () => act(async () => new Promise<void>((done) => setTimeout(done, 0)));

describe("linked table filling", () => {
  it("paints a first chunk, then every row, busy until the last arrives", async () => {
    const host = mount(table(null));
    expect(count(host)).toBe(100);
    expect(region(host).getAttribute("aria-busy")).toBe("true");
    for (let k = 0; k < 10 && count(host) < rows.length; k += 1) await tick();
    expect(count(host)).toBe(300);
    expect(region(host).getAttribute("aria-busy")).toBe("false");
    expect(host.querySelector("tbody tr:last-child button")?.textContent).toBe("r299");
  });

  it("draws a pinned row at once, wherever it sits", () => {
    const host = mount(table("r250"));
    expect(count(host)).toBe(251);
    const pinned = [...host.querySelectorAll<HTMLButtonElement>("[data-row-open]")].filter((b) => b.getAttribute("aria-pressed") === "true");
    expect(pinned.map((b) => b.textContent)).toEqual(["r250"]);
  });
});
