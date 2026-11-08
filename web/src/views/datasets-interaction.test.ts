// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi, type Mock } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import type { Manifest } from "../data/schema";
import { navigate, type Navigator } from "../app/router";
import { DatasetsView } from "./DatasetsView";
import { HOLD_MS, LIFT_CLICK_MS } from "../ui/useHold";

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
      expect.stringContaining("Has MIC prediction"),
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

/** A pointer event of TYPE from POINTER at (X, Y), as a browser sends it. */
const press = (el: Element, type: string, pointerType: string, x = 10, y = 10) => {
  const e = new MouseEvent(type, { bubbles: true, cancelable: true, clientX: x, clientY: y });
  Object.defineProperty(e, "pointerType", { value: pointerType });
  el.dispatchEvent(e);
  return e;
};
/** A touch held MS milliseconds, then lifted, and the click a browser may
 *  send after it. */
const hold = (el: Element, ms = HOLD_MS + 100, pointerType = "touch") =>
  act(() => {
    press(el, "pointerdown", pointerType);
    vi.advanceTimersByTime(ms);
    press(el, "pointerup", pointerType);
    el.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, detail: 1 }));
  });

describe("rule cells on touch and hover (A34)", () => {
  const base = manifest.datasets[0]!;
  const long = "presence of Aminoglycoside N(6')-acetyltransferase (EC 2.3.1.82) => AAC(6')-Ib/AAC(6')-II";
  beforeEach(() => vi.useFakeTimers({ toFake: ["setTimeout", "clearTimeout"] }));
  afterEach(() => vi.useRealTimers());

  it("open a cut-off rule with a touch hold as a tap does, and close it with another", () => {
    clipping(true);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: [long] }] });
    const whole = table().querySelector("tbody tr")!.lastElementChild!.querySelector("button")!;
    // Held but not long enough: the lift's click is a tap, and opens it.
    hold(whole, HOLD_MS - 100);
    expect(whole.getAttribute("aria-expanded")).toBe("true");
    hold(whole, HOLD_MS - 100);
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    // A hold opens it once: the click after the lift does not undo it.
    hold(whole);
    expect(whole.getAttribute("aria-expanded")).toBe("true");
    hold(whole);
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    // A hold without a click after it opens it as well.
    act(() => {
      press(whole, "pointerdown", "touch");
      vi.advanceTimersByTime(HOLD_MS);
    });
    expect(whole.getAttribute("aria-expanded")).toBe("true");
    act(() => press(whole, "pointerup", "touch"));
    // Then the keyboard still works: Enter's click (detail 0) is never dropped.
    act(() => whole.dispatchEvent(new MouseEvent("click", { bubbles: true, detail: 0 })));
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    expect(opened).not.toHaveBeenCalled();
  });

  it("swallow the click after a hold wherever it lands, so a rule that shrinks under the finger never opens its row", () => {
    clipping(true);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: [long] }] });
    const row = table().querySelector("tbody tr")!;
    const whole = row.lastElementChild!.querySelector("button")!;
    const other = row.children[1] as HTMLElement;
    // Chromium sends the lift's click where the finger now is: on another cell.
    act(() => {
      press(whole, "pointerdown", "touch");
      vi.advanceTimersByTime(HOLD_MS + 100);
      press(whole, "pointerup", "touch");
      other.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, detail: 1 }));
    });
    expect(whole.getAttribute("aria-expanded")).toBe("true");
    expect(opened).not.toHaveBeenCalled();
    // A hold whose click never comes swallows nothing once the lift's window closes.
    act(() => {
      press(whole, "pointerdown", "touch");
      vi.advanceTimersByTime(HOLD_MS + 100);
      press(whole, "pointerup", "touch");
      vi.advanceTimersByTime(LIFT_CLICK_MS + 1);
    });
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    act(() => other.dispatchEvent(new MouseEvent("click", { bubbles: true, cancelable: true, detail: 1 })));
    expect(opened).toHaveBeenCalledTimes(1);
  });

  it("drop a hold that moves, is cancelled or comes from a mouse, and keep the browser's menu off a hold", () => {
    clipping(true);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: [long] }] });
    const whole = table().querySelector("tbody tr")!.lastElementChild!.querySelector("button")!;
    act(() => {
      press(whole, "pointerdown", "touch", 10, 10);
      press(whole, "pointermove", "touch", 30, 10);
      vi.advanceTimersByTime(HOLD_MS * 2);
    });
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    act(() => {
      press(whole, "pointerdown", "touch");
      press(whole, "pointercancel", "touch");
      vi.advanceTimersByTime(HOLD_MS * 2);
    });
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    // A small wobble is still a hold.
    act(() => {
      press(whole, "pointerdown", "touch", 10, 10);
      press(whole, "pointermove", "touch", 14, 13);
      vi.advanceTimersByTime(HOLD_MS);
    });
    expect(whole.getAttribute("aria-expanded")).toBe("true");
    let menu!: Event;
    act(() => {
      menu = press(whole, "contextmenu", "touch");
      press(whole, "pointerup", "touch");
    });
    expect(menu.defaultPrevented).toBe(true);
    // A held mouse button is a click, once, when it lifts.
    hold(whole, HOLD_MS * 2, "mouse");
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    act(() => {
      press(whole, "pointerdown", "mouse");
      vi.advanceTimersByTime(HOLD_MS * 2);
    });
    expect(whole.getAttribute("aria-expanded")).toBe("false");
    act(() => press(whole, "pointerup", "mouse"));
  });

  it("open every rule with a touch hold on the first one, without opening the row", () => {
    draw();
    const more = table().querySelector<HTMLButtonElement>("button[aria-label^='Show all'][aria-label$='rules']")!;
    const cell = more.closest("td")!;
    const n = Number(/Show all (\d+) rules/.exec(more.getAttribute("aria-label")!)![1]);
    hold(cell.querySelector("[data-row-ignore]")!);
    expect(cell.querySelectorAll("li")).toHaveLength(n);
    expect(more.getAttribute("aria-expanded")).toBe("true");
    expect(opened).not.toHaveBeenCalled();
  });

  it("name a cut-off rule's whole text in a native tooltip, only while it is cut off", () => {
    clipping(true);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: [long] }] });
    const cell = table().querySelector("tbody tr")!.lastElementChild!;
    const text = () => cell.querySelector("button > span")!;
    expect(text().getAttribute("title")).toBe(long);
    // Hovering is a resting pointer: nothing opens, and the tooltip is the
    // browser's own, so no element is added.
    const before = cell.innerHTML;
    act(() => {
      press(text(), "pointerover", "mouse");
      press(text(), "pointermove", "mouse");
      vi.advanceTimersByTime(HOLD_MS * 2);
    });
    expect(cell.innerHTML).toBe(before);
    expect(document.querySelector("[role='tooltip']")).toBeNull();
    // Shown whole, it needs no tooltip.
    act(() => cell.querySelector("button")!.click());
    expect(text().hasAttribute("title")).toBe(false);

    // The first of several rules, cut off, names its own text; one that fits has none.
    act(() => root.unmount());
    root = createRoot(host);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: [long, "absence of blaTEM"] }] });
    const first = table().querySelector("tbody tr [data-row-ignore] > span")!;
    expect(first.getAttribute("title")).toBe(long);
    clipping(false);
    act(() => root.unmount());
    root = createRoot(host);
    draw({ ...manifest, datasets: [{ ...base, scm_rules: ["presence of blaKPC"] }] });
    expect(table().querySelector("tbody tr")!.lastElementChild!.querySelector("[title]")).toBeNull();
  });
});

describe("dataset search", () => {
  it("finds a dataset by its species' everyday name, and still by drug then species", () => {
    const base = manifest.datasets[0]!;
    // "gonorrhea" is not part of "Neisseria gonorrhoeae", so only the
    // everyday name can match it.
    const ngo = { ...base, id: "neisseria_gonorrhoeae__ceftriaxone", species: "Neisseria gonorrhoeae", antibiotic: "ceftriaxone" };
    const kpn = { ...base, id: "klebsiella_pneumoniae__ceftriaxone", species: "Klebsiella pneumoniae", antibiotic: "ceftriaxone" };
    draw({ ...manifest, datasets: [ngo, kpn] });
    const search = host.querySelector<HTMLInputElement>("input[type='search']")!;
    act(() => type(search, "gonorrhea"));
    expect(live()).toBe("1 of 2 datasets");
    expect(table().querySelector("tbody tr")!.textContent).toContain("Gonorrhea vs\u00a0ceftriaxone");
    act(() => type(search, "ceftriaxone klebsiella pneumoniae"));
    expect(live()).toBe("1 of 2 datasets");
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
