// @vitest-environment happy-dom
// The deck's paging rule: when a slide's title, answer, figure at its least
// height and caption do not fit one page, the slide takes a second page with
// the caption. The heights below are the MIC slide's, read in a browser from
// the deck's one-page copy of the slide (the same reading as readSlideFit):
// on the revision 14 release after `grm baselines`, and on the synthetic
// `full` report shape, at 390x844, 1024x768 and 1280x720. The strip's least
// height, 167.03 px, is its half band held at 28 px with its axis title
// wrapped to two lines at that width.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { CatalogRow, Manifest } from "../data/schema";
import { FINDINGS, micFinding } from "../views/findings";
import { primerSteps } from "../views/primer";
import { PresentMode } from "./PresentMode";
import { fitsOnePage, ONE_PAGE, pageShows, readSlideFit, SLACK, slidePlan, type SlideFit } from "./paging";

// The frame's calls to slidePlan go through this spy, so a test can see
// that the deck decides with it.
vi.mock("./paging", async (importOriginal) => {
  const actual = await importOriginal<typeof import("./paging")>();
  return { ...actual, slidePlan: vi.fn(actual.slidePlan) };
});

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

/** The MIC slide at 390x844, revision 14 after `grm baselines`: the figure's
 *  box is 120.5 px and the strip draws 167.03. */
const MIC_390: SlideFit = { room: 663, head: 267.25, figure: 167.03, caption: 275.25 };
/** The same slide on the synthetic `full` shape: a box of 148.09 px. */
const MIC_390_FULL: SlideFit = { room: 663, head: 239.66, figure: 167.03, caption: 275.25 };
/** At 1024x768 the strip fills the 247.88 px it is given. */
const MIC_1024: SlideFit = { room: 623, head: 180.88, figure: 247.52, caption: 194.25 };
/** At 1280x720, 0.35 px to spare without the slack. */
const MIC_1280: SlideFit = { room: 575, head: 180.88, figure: 199.52, caption: 194.25 };

const need = (f: SlideFit) => f.head + f.figure + f.caption;

describe("one page or two", () => {
  it("gives the MIC slide on a phone a second page for its caption, on both trees with baselines files", () => {
    expect(need(MIC_390)).toBeCloseTo(709.53, 2);
    expect(slidePlan(MIC_390, false)).toEqual({ pages: 2, caption: "last" });
    expect(need(MIC_390_FULL)).toBeCloseTo(681.94, 2);
    expect(slidePlan(MIC_390_FULL, false)).toEqual({ pages: 2, caption: "last" });
  });

  it("keeps a slide that fits on one page, at 1024 and 1280 px", () => {
    expect(fitsOnePage(MIC_1024)).toBe(true);
    expect(slidePlan(MIC_1024, false)).toEqual(ONE_PAGE);
    expect(slidePlan(MIC_1280, false)).toEqual(ONE_PAGE);
  });

  it("splits at the first pixel past the slack, and not before", () => {
    // Quarter pixels add up exactly, so the edge is the rule's, not rounding's.
    const at = { room: 708.5, head: 267.25, figure: 167, caption: 275.25 };
    expect(need(at)).toBe(at.room + SLACK);
    expect(slidePlan(at, false)).toEqual(ONE_PAGE);
    // The strip one pixel taller than the slide leaves it.
    expect(slidePlan({ ...at, figure: at.figure + 1 }, false)).toEqual({ pages: 2, caption: "last" });
    expect(slidePlan({ ...at, room: at.room - 1 }, false).pages).toBe(2);
  });

  it("counts every part: a taller title or answer, or a longer caption, moves the caption", () => {
    const at = { ...MIC_1280 };
    expect(slidePlan(at, false).pages).toBe(1);
    expect(slidePlan({ ...at, head: at.head + 2 }, false).pages).toBe(2);
    expect(slidePlan({ ...at, caption: at.caption + 2 }, false).pages).toBe(2);
    expect(slidePlan({ ...at, room: at.room - 2 }, false).pages).toBe(2);
  });

  it("treats an unlaid slide, or one not yet read, as one page", () => {
    expect(slidePlan({ room: 0, head: 0, figure: 400, caption: 0 }, false)).toEqual(ONE_PAGE);
    expect(slidePlan(null, false)).toEqual(ONE_PAGE);
  });

  it("lets a figure that marks itself data-continues keep its caption on both pages", () => {
    expect(slidePlan(MIC_390, true)).toEqual({ pages: 2, caption: "each" });
    expect(slidePlan(MIC_1024, true)).toEqual(ONE_PAGE);
  });

  it("shows the figure without the caption on page 1 of a caption split, and the caption alone on page 2", () => {
    const split = slidePlan(MIC_390, false);
    expect([0, 1].map((p) => pageShows(split, p))).toEqual(["figure", "caption"]);
    expect([0, 1].map((p) => pageShows({ pages: 2, caption: "each" }, p))).toEqual(["all", "all"]);
    expect(pageShows(ONE_PAGE, 0)).toBe("all");
  });
});

/** A laid-out box, as a browser would report it. */
const box = (top: number, height: number) => ({ x: 0, y: top, top, left: 0, right: 300, bottom: top + height, width: 300, height, toJSON: () => ({}) }) as DOMRect;

describe("reading the one-page copy", () => {
  afterEach(() => {
    vi.restoreAllMocks();
    document.body.innerHTML = "";
  });

  it("reads the head to the figure's box, the lowest mark the figure draws, and the caption with its gap and the space under it", () => {
    document.body.innerHTML = `<div id="root"><h2>Title</h2><p>Answer</p><div id="box" style="padding-bottom: 16px"><div data-measure-slot=""><div id="fig"><svg><text id="low">80</text></svg><span class="sr-only" id="sr">every dot</span></div></div><p data-measure-caption="" style="margin-top: 24px">Caption</p></div></div>`;
    const at: Record<string, DOMRect> = { root: box(100, 663), fig: box(367.25, 150), low: box(500, 34.03), sr: box(900, 1) };
    vi.spyOn(Element.prototype, "getBoundingClientRect").mockImplementation(function (this: Element) {
      if (this.hasAttribute("data-measure-slot")) return box(367.25, 120.5);
      if (this.hasAttribute("data-measure-caption")) return box(700, 235.25);
      return at[this.id] ?? box(0, 0);
    });
    const fit = readSlideFit(document.getElementById("root")!)!;
    expect(fit.room).toBe(663);
    expect(fit.head).toBeCloseTo(267.25, 5);
    // The axis label at 534.03 is the lowest mark; the screen-reader line is not drawn.
    expect(fit.figure).toBeCloseTo(166.78, 5);
    expect(fit.caption).toBeCloseTo(235.25 + 24 + 16, 5);
  });

  it("reads nothing before the figure's box and the caption are there", () => {
    document.body.innerHTML = `<div id="root"><h2>Title</h2><div data-measure-slot=""></div></div>`;
    expect(readSlideFit(document.getElementById("root")!)).toBeNull();
  });
});

// ---- the deck ------------------------------------------------------------------

const base = parseManifest(fixtureFile("/data/index.json"));
const row = (id: string, overrides: Partial<CatalogRow>): CatalogRow => ({ ...base.datasets[0]!, id, antibiotic: id, ...overrides });
/** The fixture with 31 more datasets that have an MIC model, so the MIC slide has its strip. */
const withMic: Manifest = { ...base, datasets: [...base.datasets, ...Array.from({ length: 31 }, (_, k) => row(`mic${k}`, { mic_essential_agreement: (60 + k) / 100 }))] };
/** Title, the primer's steps, the findings, the rule: one slide each. */
const DECK = 1 + primerSteps(withMic).length + FINDINGS.length + 1;
const MIC_AT = 1 + primerSteps(withMic).length + FINDINGS.findIndex((f) => f(withMic, { page: false }).id === "mic") + 1;

let root: Root | null = null;
let restore: (() => void) | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  restore?.();
  restore = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
  vi.mocked(slidePlan).mockClear();
});

/** Lays out the deck's one-page copy of the MIC slide with `fit`, as a
 *  browser would: happy-dom lays nothing out, so every other copy reads as
 *  unlaid and stays one page. */
const layoutMic = (fit: SlideFit) => {
  const rect = Element.prototype.getBoundingClientRect;
  Element.prototype.getBoundingClientRect = function (this: Element) {
    const copy = this.closest('[data-slide-measure="mic"]');
    if (copy === null) return rect.call(this);
    if (this === copy.firstElementChild) return box(0, fit.room);
    if (this.hasAttribute("data-measure-slot")) return box(fit.head, fit.room - fit.head - fit.caption);
    if (this.parentElement?.hasAttribute("data-measure-slot") === true) return box(fit.head, fit.figure);
    if (this.hasAttribute("data-measure-caption")) return box(fit.room - fit.caption, fit.caption);
    return box(fit.head, 0);
  };
  return () => {
    Element.prototype.getBoundingClientRect = rect;
  };
};

const open = async (fit: SlideFit) => {
  vi.stubGlobal("matchMedia", (query: string) => ({ matches: query.includes("prefers-reduced-motion"), media: query, addEventListener: () => {}, removeEventListener: () => {} }));
  vi.stubGlobal("fetch", async (path: string) => new Response(JSON.stringify(fixtureFile(path))));
  restore = layoutMic(fit);
  const host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  await act(async () => root!.render(createElement(PresentMode, { manifest: withMic, onClose: () => {} })));
  const dialog = document.querySelector<HTMLDialogElement>("dialog")!;
  const status = () => dialog.querySelector('[role="status"]')?.textContent ?? "";
  const go = async (n: number) => {
    while (!status().startsWith(`Slide ${n} `)) await act(async () => dialog.querySelector("main")!.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowRight", bubbles: true })));
  };
  const section = () => dialog.querySelector("main section")!;
  const beats = () => dialog.querySelectorAll("[data-beats] span").length;
  return { dialog, status, go, section, beats };
};

describe("the deck pages a slide by slidePlan", () => {
  const mic = micFinding(withMic, { page: false });

  it("puts the MIC slide's figure on page 1 under its title and answer, and its caption alone under the title on page 2", async () => {
    const { status, go, section, beats } = await open(MIC_390);
    // The frame read the copy and asked slidePlan.
    expect(vi.mocked(slidePlan)).toHaveBeenCalledWith({ room: 663, head: 267.25, figure: expect.closeTo(167.03, 6), caption: 275.25 }, false);
    // The copy lays out above a copy of the footer row the slide shows as one
    // page, so its room is that page's whichever page is showing.
    const footer = document.querySelector('[data-slide-measure="mic"]')!.parentElement!.nextElementSibling!;
    expect(footer.hasAttribute("data-footer-measure")).toBe(true);
    expect(footer.querySelector("[data-dropdown-label]")?.textContent).toBe(`${MIC_AT} · ${mic.title}`);
    expect(footer.querySelectorAll("span.size-1\\.5")).toHaveLength(3);
    await go(MIC_AT);
    expect(status()).toBe(`Slide ${MIC_AT} of ${DECK + 1}: ${mic.title} (1 of 2)`);
    expect(section().querySelector("h2")?.textContent).toBe(mic.title);
    expect(section().textContent).toContain(mic.result);
    expect(section().querySelector("[data-figure] [data-figure-fit] svg")).not.toBeNull();
    expect(section().querySelectorAll("[data-dataset]")).toHaveLength(withMic.datasets.filter((d) => d.mic_essential_agreement !== null).length);
    expect(section().querySelector("[data-caption]")).toBeNull();
    expect(section().textContent).not.toContain(mic.caption);
    // The answer and the figure are its two beats.
    expect(beats()).toBe(2);
    await go(MIC_AT + 1);
    expect(status()).toBe(`Slide ${MIC_AT + 1} of ${DECK + 1}: ${mic.title} (2 of 2)`);
    expect(section().getAttribute("aria-labelledby")).toBe("mic-title");
    expect(section().querySelector("h2")?.textContent).toBe(mic.title);
    const caption = section().querySelector("h2 + * [data-caption], h2 + [data-caption]");
    expect(caption?.textContent).toBe(mic.caption);
    // A screen reader meets it after the title, as on a one-page slide.
    expect(caption?.closest('[aria-hidden="true"], [inert]')).toBeNull();
    expect(section().querySelector("[data-figure]")).toBeNull();
    expect(section().querySelector("svg")).toBeNull();
    expect(section().textContent).not.toContain(mic.result);
    // One beat, so no ticks.
    expect(beats()).toBe(0);
    // The next slide follows on.
    await go(MIC_AT + 2);
    expect(status()).toContain(`Slide ${MIC_AT + 2} of ${DECK + 1}: Newly deposited genomes`);
  });

  it("keeps the slide whole when the four fit with the slack to spare, and splits it one pixel later", async () => {
    const fits = { room: 708.5, head: 267.25, figure: 167, caption: 275.25 };
    const first = await open(fits);
    expect(vi.mocked(slidePlan)).toHaveBeenCalledWith(fits, false);
    await first.go(MIC_AT);
    expect(first.status()).toBe(`Slide ${MIC_AT} of ${DECK}: ${mic.title}`);
    expect(first.section().querySelector("[data-caption]")?.textContent).toBe(mic.caption);
    expect(first.section().querySelector("[data-figure] svg")).not.toBeNull();
    act(() => root?.unmount());
    restore?.();
    document.body.innerHTML = "";
    const over = { ...fits, figure: fits.figure + 1 };
    const second = await open(over);
    expect(vi.mocked(slidePlan)).toHaveBeenCalledWith(over, false);
    await second.go(MIC_AT);
    expect(second.status()).toBe(`Slide ${MIC_AT} of ${DECK + 1}: ${mic.title} (1 of 2)`);
  });

  it("keeps a one-page copy of every figure slide, and pages none on its own when every copy fits", async () => {
    await open(MIC_1024);
    expect(vi.mocked(slidePlan)).toHaveBeenCalledWith({ room: 623, head: 180.88, figure: expect.closeTo(247.52, 6), caption: 194.25 }, false);
    // Every slide but the title has its copy.
    const ids = new Set([...document.querySelectorAll("[data-slide-measure]")].map((e) => e.getAttribute("data-slide-measure")));
    expect(ids.size).toBe(DECK - 1);
    // A copy adds no keyboard layer: the only listboxes are the live slide's.
    expect(document.querySelectorAll("[data-slide-measure] [role='listbox'], [data-slide-measure] [role='option']")).toHaveLength(0);
    const status = document.querySelector('dialog [role="status"]')?.textContent ?? "";
    expect(status).toContain(`of ${DECK}`);
  });
});
