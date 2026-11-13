// @vitest-environment happy-dom
// The input table's paging: for every page count from one to the row count,
// the pages together show every row exactly once, in order, split as evenly
// as they go, with no page empty.
import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { INPUT_COMPARISON } from "../data/derive/inputs";
import { InputRows, pageRows, type InputGroup } from "./InputRows";

const groups: InputGroup[] = [
  { title: "changed", rows: INPUT_COMPARISON.arms.filter((a) => !a.twin).map((a) => ({ key: a.id, name: a.id, ...a.all })) },
  { title: "twins", rows: INPUT_COMPARISON.arms.filter((a) => a.twin).map((a) => ({ key: a.id, name: a.id, ...a.all })) },
];
const keys = groups.flatMap((g) => g.rows.map((r) => r.key));

const shownOn = (page: number, pages: number, slide: boolean) => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(createElement(InputRows, { label: "rows", groups, slide, page, pages }));
  return { rows: [...el.querySelectorAll("tr[data-input]")].map((tr) => tr.getAttribute("data-input")), pager: el.querySelector("[data-continued]")?.textContent ?? null };
};

describe("the input table's pages", () => {
  it("splits n rows into pages no more than one row apart, none empty", () => {
    for (let n = 1; n <= 20; n++)
      for (let count = 1; count <= n; count++) {
        const sizes = Array.from({ length: count }, (_, i) => {
          const [a, b] = pageRows(n, count, i);
          return b - a;
        });
        expect(sizes.reduce((x, y) => x + y, 0)).toBe(n);
        expect(Math.min(...sizes)).toBeGreaterThan(0);
        expect(Math.max(...sizes) - Math.min(...sizes)).toBeLessThanOrEqual(1);
      }
    expect([0, 1, 2, 3].map((i) => pageRows(13, 4, i))).toEqual([
      [0, 3],
      [3, 6],
      [6, 9],
      [9, 13],
    ]);
  });

  it("shows every row exactly once, in order, for every page count from 1 to the row count", () => {
    expect(keys).toHaveLength(13);
    for (const slide of [true, false])
      for (let pages = 1; pages <= keys.length; pages++) {
        const seen = Array.from({ length: pages }, (_, page) => shownOn(page, pages, slide));
        expect(seen.flatMap((s) => s.rows)).toEqual(keys);
        expect(seen.every((s) => s.rows.length > 0)).toBe(true);
        if (pages > 1) expect(seen.at(-1)!.pager).toBe(`Rows ${keys.length - seen.at(-1)!.rows.length + 1}–${keys.length} of ${keys.length}`);
      }
  });
});

describe("the input table's legend", () => {
  // Unlaid, the table is 640 px wide, wide enough for every column.
  const render = (slide: boolean) => {
    const el = document.createElement("div");
    el.innerHTML = renderToStaticMarkup(createElement(InputRows, { label: "rows", groups, slide }));
    return { head: [...el.querySelectorAll("thead th")].map((th) => th.textContent), text: el.textContent ?? "" };
  };

  it("explains the datasets-compared and median columns where it shows them, on a slide in short", () => {
    const slide = render(true);
    expect(slide.head).toEqual(["Change", "Worse", "Worse and better, drawn", "Better", "Of", "Sign p", "Median"]);
    expect(slide.text).toContain("Of: datasets compared. Median: the middle per-dataset change, above zero when better.");
    const page = render(false);
    expect(page.text).toContain("Of: datasets compared; on any not counted worse or better, the score was exactly the same.");
  });
});
