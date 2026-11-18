// @vitest-environment happy-dom
// The strip's least height. Given less height than it needs, as a slide
// squeezed by its headline and caption gives it, the strip keeps a band of
// 28 px on each side of its midline and draws past its box; the slide frame
// reads that drawn height to decide whether the caption needs a page of its
// own (web/src/present/paging.ts). A smaller floor would let a squeezed
// strip draw a band too thin for its marks.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it } from "vitest";
import { ThemeProvider } from "../theme/ThemeProvider";
import { DotStrip } from "./DotStrip";
import { MeasuringCopy } from "./kit";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let root: Root | null = null;

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
});

/** The strip given `height`, live or inside a slide's one-page measuring
 *  copy. `read` sees it mounted; it is unmounted after. */
const strip = <T>(height: number, measuring: boolean, read: (host: HTMLElement) => T): T => {
  const host = document.createElement("div");
  document.body.append(host);
  const r = createRoot(host);
  root = r;
  const data = [
    { id: "a", value: 0.2, tip: [{ term: "A", text: "20%" }] },
    { id: "b", value: 0.8, tip: [{ term: "B", text: "80%" }] },
  ];
  const chart = createElement(DotStrip, { label: "Test strip", data, domain: [0, 1], format: String, axisTitle: "Agreement", refs: [], regions: [], named: [], height });
  act(() => r.render(createElement(ThemeProvider, { children: measuring ? createElement(MeasuringCopy.Provider, { value: true }, chart) : chart })));
  const out = read(host);
  act(() => r.unmount());
  root = null;
  host.remove();
  return out;
};

const svgHeight = (host: HTMLElement): number => {
  const svg = host.querySelector("svg[aria-hidden]");
  if (svg === null) throw new Error("no svg");
  return Number(svg.getAttribute("height"));
};

/** The strip's own svg height when the slide gives it `height`. */
const drawn = (height: number): number => strip(height, false, svgHeight);

// Three label rows of 16 px and 6 px above the band, and 22 px for the axis.
const TOP = 3 * 16 + 6;
const AXIS = 22;

describe("the strip's least height", () => {
  it("keeps a band of 28 px each side of its midline however little height it is given", () => {
    expect(drawn(60)).toBe(TOP + 2 * 28 + AXIS);
    expect(drawn(1)).toBe(TOP + 2 * 28 + AXIS);
  });

  it("grows its band with the height it is given above that floor", () => {
    expect(drawn(400)).toBeGreaterThan(drawn(200));
    expect(drawn(200)).toBeGreaterThan(TOP + 2 * 28 + AXIS);
  });
});

describe("the strip in a slide's one-page measuring copy", () => {
  it("draws no keyboard layer, keeps its axis title under the plot, and draws as tall as the live strip", () => {
    const copy = strip(400, true, (host) => {
      const svg = host.querySelector("svg[aria-hidden]");
      const title = [...host.querySelectorAll("p")].find((p) => p.textContent === "Agreement");
      return {
        listbox: host.querySelector('[role="listbox"]') !== null,
        title: title !== undefined,
        under: svg !== null && title !== undefined && (svg.compareDocumentPosition(title) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0,
        height: svgHeight(host),
      };
    });
    // The axis title sets the band the strip draws, so a copy without it
    // would measure a strip of another height than the slide shows.
    expect(copy).toEqual({ listbox: false, title: true, under: true, height: drawn(400) });
  });
});
