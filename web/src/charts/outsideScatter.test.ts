// @vitest-environment happy-dom
// The outside-study chart against the derived counts: each plot's counts line
// must say what populationSummary counts from the same rows, so a line that
// swapped "above" and "below", or counted another measure, fails here and not
// only in the browser check.
import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { outsidePoints, populationById, populationSummary, type Population } from "../data/derive/outside";
import { OutsideScatter, type OutsidePanel } from "./OutsideScatter";

const panel = (p: Population): OutsidePanel => ({
  key: p.id,
  heading: p.id,
  name: p.id,
  points: { sensitivity: outsidePoints(p, "sensitivity"), specificity: outsidePoints(p, "specificity") },
});

const render = (panels: OutsidePanel[], page?: number, pages?: number, height?: number) => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(createElement(OutsideScatter, { panels, level: "95%", titleOf: (id: string) => id, page, pages, height }));
  return el;
};

describe("the outside-study chart's counts lines", () => {
  const populations = [populationById("primary"), populationById("same-study")];

  it("say, for each population and measure, the reading counts the summary derives", () => {
    const el = render(populations.map(panel));
    const lines = [...el.querySelectorAll("[data-outside-plot]")].map((p) => [p.getAttribute("data-outside-plot"), p.querySelector("[data-outside-counts]")?.textContent]);
    const expected = populations.flatMap((p) => {
      const s = populationSummary(p).scm;
      return [
        ["sensitivity", `${s.sensitivityAbove} read above, ${s.sensitivityBelow} read below`],
        ["specificity", `${s.specificityAbove} read above, ${s.specificityBelow} read below`],
      ];
    });
    expect(lines).toEqual(expected);
    // The primary population's lines as the document gives them.
    expect(lines.slice(0, 2).map((l) => l[1])).toEqual(["36 read above, 27 read below", "25 read above, 35 read below"]);
  });

  it("draws one mark per reading, with its reading", () => {
    const el = render(populations.map(panel));
    const plots = [...el.querySelectorAll("[data-outside-plot]")];
    populations.forEach((p, i) => {
      for (const [k, m] of (["sensitivity", "specificity"] as const).entries()) {
        const marks = [...plots[i * 2 + k]!.querySelectorAll("[data-reading]")].map((g) => g.getAttribute("data-reading"));
        expect(marks.sort()).toEqual(outsidePoints(p, m).map((x) => x.reading).sort());
      }
    });
  });

  it("labels only 0 and 100 on a plot too small for a third label, and keeps every gridline", () => {
    // One population on a page, unlaid: the plot takes the height less the
    // legend's 48 px and the head's 66 px, so 194 px draws an 80 px plot and
    // 250 px a 136 px one, either side of the 100 px that labels 50.
    const plot = (height: number) => render(populations.map(panel), 0, 2, height).querySelector("[data-outside-plot] svg")!;
    const small = plot(194), large = plot(250);
    expect(small.getAttribute("width")).toBe("80");
    expect([...small.querySelectorAll("text")].map((t) => t.textContent)).toEqual(["0", "", "100", "0", "", "100"]);
    expect(small.querySelectorAll("line[stroke='var(--chart-grid)']")).toHaveLength(6);
    expect(large.getAttribute("width")).toBe("136");
    expect([...large.querySelectorAll("text")].map((t) => t.textContent)).toEqual(["0", "50", "100", "0", "50", "100"]);
  });

  it("says in its legend what each axis reads, with lineages held out and the interval's level", () => {
    expect(render(populations.map(panel)).textContent).toContain("Across: the report's figure with lineages held out. Up: the new figure, with its 95% interval as a bar. Both in %.");
  });

  it("shows each population once across two pages", () => {
    const panels = populations.map(panel);
    const seen = [0, 1].flatMap((page) => [...render(panels, page, 2).querySelectorAll("[data-outside-panel] h3")].map((h) => h.textContent));
    expect(seen).toEqual(["primary", "same-study"]);
  });
});
