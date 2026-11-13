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

const render = (panels: OutsidePanel[], page?: number, pages?: number) => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(createElement(OutsideScatter, { panels, titleOf: (id: string) => id, page, pages }));
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

  it("shows each population once across two pages", () => {
    const panels = populations.map(panel);
    const seen = [0, 1].flatMap((page) => [...render(panels, page, 2).querySelectorAll("[data-outside-panel] h3")].map((h) => h.textContent));
    expect(seen).toEqual(["primary", "same-study"]);
  });
});
