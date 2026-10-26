// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi, type Mock } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { Manifest } from "../data/schema";
import { navigate, type Navigator } from "../app/router";
import { DatasetsView } from "./DatasetsView";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const manifest = parseManifest(fixtureFile("/data/index.json"));

/** A wide desktop (the docked rail) with reduced motion, so counts snap and
 *  sorts do not animate. */
const media = (q: string): MediaQueryList =>
  ({
    matches: q.includes("min-width: 1440px") || q.includes("prefers-reduced-motion"),
    media: q,
    onchange: null,
    addEventListener: () => {},
    removeEventListener: () => {},
    addListener: () => {},
    removeListener: () => {},
    dispatchEvent: () => false,
  }) as MediaQueryList;

let host: HTMLDivElement;
let root: Root;
let opened: Mock<(id: string) => void>;
let nav: Navigator;

beforeEach(() => {
  window.matchMedia = media;
  window.history.replaceState(null, "", "/datasets");
  host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  opened = vi.fn<(id: string) => void>();
  nav = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: opened, closeDataset: () => {} };
});

afterEach(() => {
  act(() => root.unmount());
  host.remove();
  vi.restoreAllMocks();
  unclip();
});

const draw = (m: Manifest = manifest) => act(() => root.render(createElement(DatasetsView, { manifest: m, nav })));

const rail = () => host.querySelector<HTMLElement>("aside[aria-label='Filters']")!;
const table = () => host.querySelector<HTMLElement>("[role='region'][aria-label='Datasets']")!;
const live = () => host.querySelector("[aria-live='polite']")!.textContent;

const type = (input: HTMLInputElement, value: string) => {
  Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, "value")!.set!.call(input, value);
  input.dispatchEvent(new Event("input", { bubbles: true }));
};

/** Every truncated line reports itself clipped (or not), as layout would. */
const clipping = (clipped: boolean) => {
  Object.defineProperty(HTMLElement.prototype, "scrollWidth", {
    configurable: true,
    get(this: HTMLElement) {
      return this.classList.contains("truncate") && clipped ? 400 : 50;
    },
  });
  Object.defineProperty(HTMLElement.prototype, "clientWidth", { configurable: true, get: () => 100 });
};

const unclip = () => {
  delete (HTMLElement.prototype as { scrollWidth?: number }).scrollWidth;
  delete (HTMLElement.prototype as { clientWidth?: number }).clientWidth;
};

/** A resize observer the test fires by hand, as a browser does after each
 *  layout: `relayout` reports every observed box again. */
const observed = new Map<Element, ResizeObserverCallback>();
class HandResizeObserver {
  constructor(private readonly callback: ResizeObserverCallback) {}
  observe(el: Element) {
    observed.set(el, this.callback);
  }
  unobserve(el: Element) {
    observed.delete(el);
  }
  disconnect() {}
}
window.ResizeObserver = HandResizeObserver as unknown as typeof ResizeObserver;
const relayout = () =>
  act(() => {
    for (const [el, callback] of observed) callback([{ target: el } as unknown as ResizeObserverEntry], {} as ResizeObserver);
  });

describe("dataset filters", () => {
  it("asks each question once, with native checkboxes and the bars grouped under Meets the bar", () => {
    draw();
    expect(rail().querySelectorAll("[aria-pressed]")).toHaveLength(0);
    const boxes = [...rail().querySelectorAll<HTMLInputElement>("input[type='checkbox']")];
    expect(boxes.map((b) => b.closest("label")!.textContent)).toEqual([
      expect.stringContaining("Has MIC regression"),
      expect.stringContaining("VME ≤1.5%"),
      expect.stringContaining("ME ≤3%"),
    ]);
    const bars = boxes[1]!.closest("fieldset")!;
    expect(bars.querySelector("legend")!.textContent).toBe("Meets the bar");
    expect(boxes[2]!.closest("fieldset")).toBe(bars);

    const labels = [...rail().querySelectorAll(".filter-group-label, .filter-field-label")].map((e) => e.textContent!.toLowerCase());
    expect(new Set(labels).size).toBe(labels.length);
  });

  it("returns the rows a checkbox's count promised", () => {
    draw();
    const vme = rail().querySelectorAll<HTMLInputElement>("input[type='checkbox']")[1]!;
    const promised = Number(vme.closest("label")!.lastElementChild!.textContent!.replace(/,/g, ""));
    act(() => vme.click());
    expect(vme.checked).toBe(true);
    expect(live()).toBe(`${promised} of ${manifest.datasets.length} datasets`);
  });
});

describe("rule cells", () => {
  it("open every rule in place from +N, keep focus on that one toggle, and never open the row", () => {
    draw();
    const more = table().querySelector<HTMLButtonElement>("button[aria-label^='Show all'][aria-label$='rules']")!;
    const cell = more.closest("td")!;
    const row = more.closest("tr")!;
    const n = Number(/Show all (\d+) rules/.exec(more.getAttribute("aria-label")!)![1]);

    more.focus();
    act(() => more.click());
    expect(more.getAttribute("aria-expanded")).toBe("true");
    expect(more.textContent).toBe("Fewer");
    expect(document.activeElement).toBe(more);
    expect(cell.querySelectorAll("li")).toHaveLength(n);

    act(() => more.click());
    expect(more.getAttribute("aria-expanded")).toBe("false");
    expect(more.textContent).toBe(`+${n - 1}`);

    // The first rule's text opens the list too.
    act(() => cell.querySelector<HTMLElement>("[data-row-ignore]")!.click());
    expect(cell.querySelectorAll("li")).toHaveLength(n);
    expect(opened).not.toHaveBeenCalled();

    // Anywhere else on the row still opens the dataset.
    act(() => (row.children[1] as HTMLElement).click());
    expect(opened).toHaveBeenCalledTimes(1);
  });

  it("make a single cut-off rule its own button for the whole text, and leave a rule that fits as text", () => {
    const base = manifest.datasets[0]!;
    const long = "presence of Aminoglycoside N(6')-acetyltransferase (EC 2.3.1.82) => AAC(6')-Ib/AAC(6')-II";
    clipping(true);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: [long] }] });
    const cell = table().querySelector("tbody tr")!.lastElementChild!;
    const whole = cell.querySelector("button")!;
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    expect(whole.textContent).toBe(long);
    whole.focus();
    act(() => whole.click());
    expect(cell.querySelector("button")).toBe(whole);
    expect(whole.getAttribute("aria-expanded")).toBe("true");
    expect(opened).not.toHaveBeenCalled();

    // Shown whole, the rule no longer measures as clipped; after that layout
    // the same button still closes it, and keeps focus while doing so.
    relayout();
    expect(cell.querySelector("button")).toBe(whole);
    act(() => whole.click());
    relayout();
    expect(cell.querySelector("button")).toBe(whole);
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    expect(document.activeElement).toBe(whole);

    clipping(false);
    act(() => root.unmount());
    root = createRoot(host);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: ["presence of blaKPC"] }] });
    const fits = table().querySelector("tbody tr")!.lastElementChild!;
    expect(fits.querySelector("button")).toBeNull();
    act(() => (fits as HTMLElement).click());
    expect(opened).toHaveBeenCalledTimes(1);
  });
});

describe("returning to the catalog", () => {
  it("brings back the search, filters, sort and Show all after a trip away and Back", () => {
    draw();
    act(() => type(host.querySelector<HTMLInputElement>("input[type='search']")!, "a"));
    act(() => rail().querySelectorAll<HTMLInputElement>("input[type='checkbox']")[2]!.click());
    const isolates = [...table().querySelectorAll("th")].find((th) => th.textContent!.startsWith("Isolates"))!;
    act(() => isolates.querySelector("button")!.click());
    const shown = live();

    const writes = vi.spyOn(window.history, "replaceState");
    act(() => {
      navigate("/species");
    });
    const left = writes.mock.calls.at(-1)![0] as { views: Record<string, unknown> };
    expect(left.views["datasets.query"]).toBe("a");

    // Back: the catalog is drawn anew from the entry it left.
    act(() => root.unmount());
    window.history.replaceState(left, "", "/datasets");
    root = createRoot(host);
    draw();
    expect(host.querySelector<HTMLInputElement>("input[type='search']")!.value).toBe("a");
    expect(rail().querySelectorAll<HTMLInputElement>("input[type='checkbox']")[2]!.checked).toBe(true);
    expect([...table().querySelectorAll("th")].find((th) => th.getAttribute("aria-sort") !== null)!.textContent).toContain("Isolates");
    expect(live()).toBe(shown);
  });
});
