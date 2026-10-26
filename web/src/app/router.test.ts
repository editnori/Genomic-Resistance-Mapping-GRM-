// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { findFocusToken, focusToken, isAppPath, LIMITS_HREF, navigate, routeHref, useKeptState, useRouter, type Navigator } from "./router";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let host: HTMLDivElement;
let root: Root;

beforeEach(() => {
  host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
});

afterEach(() => {
  act(() => root.unmount());
  host.remove();
  vi.restoreAllMocks();
});

const here = () => `${window.location.pathname}${window.location.search}${window.location.hash}`;

/** Mounts the router and hands back its navigator. */
const withRouter = (): (() => Navigator) => {
  let nav: Navigator | null = null;
  const Probe = () => {
    nav = useRouter();
    return null;
  };
  act(() => root.render(createElement(Probe)));
  return () => nav!;
};

describe("routes", () => {
  it("write a section hash after the place and keep the app's own paths apart from data files", () => {
    expect(routeHref({ view: "methods", species: null, dataset: null }, "m-funnel")).toBe("/methods#m-funnel");
    expect(routeHref({ view: "datasets", species: null, dataset: "a b" }, "#x")).toBe("/datasets?dataset=a%20b#x");
    expect(routeHref({ view: "overview", species: null, dataset: null })).toBe("/");
    expect(LIMITS_HREF).toBe("/#limits");
    for (const p of ["/", "/datasets", "/species", "/species/escherichia_coli", "/methods", "/tests"]) expect(isAppPath(p)).toBe(true);
    for (const p of ["/data/index.json", "/assets/index.js", "/favicon.svg", "/datasets/extra"]) expect(isAppPath(p)).toBe(false);
  });

  it("keep the section hash while the inspector opens over a view and closes again", () => {
    window.history.replaceState(null, "", "/methods#m-labels");
    const nav = withRouter();
    act(() => nav().openDataset("ecoli__amikacin"));
    expect(here()).toBe("/methods?dataset=ecoli__amikacin#m-labels");
    expect(nav().route).toEqual({ view: "methods", species: null, dataset: "ecoli__amikacin" });
    expect((window.history.state as { inspector?: boolean }).inspector).toBe(true);

    // An inspector reached by a link (no entry of its own) closes in place.
    window.history.replaceState(null, "", "/methods?dataset=ecoli__amikacin#m-labels");
    act(() => nav().closeDataset());
    expect(here()).toBe("/methods#m-labels");
  });

  it("go to a place and a section in it", () => {
    window.history.replaceState(null, "", "/");
    const nav = withRouter();
    act(() => nav().go({ view: "methods", hash: "m-funnel" }));
    expect(here()).toBe("/methods#m-funnel");
    expect(nav().route.view).toBe("methods");
  });
});

describe("what an entry keeps", () => {
  it("saves a view's own state, its scroll and the focused control into the entry it leaves", () => {
    window.history.replaceState(null, "", "/datasets");
    let set: ((v: string) => void) | null = null;
    const Keeper = () => {
      const [q, setQ] = useKeptState("test.query", "");
      set = setQ;
      return createElement("button", { type: "button", "data-row-open": "", "data-row-key": "kp__mero" }, q);
    };
    act(() => root.render(createElement(Keeper)));
    act(() => set!("mero"));
    host.querySelector("button")!.focus();

    const writes = vi.spyOn(window.history, "replaceState");
    act(() => {
      navigate("/species/klebsiella_pneumoniae");
    });
    expect(here()).toBe("/species/klebsiella_pneumoniae");
    const left = writes.mock.calls.at(-1)![0] as { views: Record<string, unknown>; focus: unknown; scroll: number };
    expect(left.views["test.query"]).toBe("mero");
    expect(left.focus).toEqual({ row: "kp__mero" });

    // Back: the entry it left starts the view where the reader left it.
    act(() => root.unmount());
    window.history.replaceState(left, "", "/datasets");
    root = createRoot(host);
    act(() => root.render(createElement(Keeper)));
    expect(host.querySelector("button")!.textContent).toBe("mero");
  });

  it("finds the control that left again by row, id, link or name", () => {
    host.innerHTML = `
      <button data-row-open data-row-key="a">Row A</button>
      <button id="present">Present</button>
      <a href="/methods#m-funnel">See QC outcomes</a>
      <button aria-label="Show all 4 rules">+3</button>`;
    const [row, id, link, named] = [...host.querySelectorAll<HTMLElement>("button, a")];
    for (const el of [row!, id!, link!, named!]) {
      const token = focusToken(el)!;
      expect(findFocusToken(token)).toBe(el);
    }
    expect(focusToken(row!)).toEqual({ row: "a" });
    expect(focusToken(link!)).toEqual({ href: "/methods#m-funnel" });
    expect(focusToken(document.body)).toBeNull();
  });

  it("does not follow an href the app does not draw", () => {
    window.history.replaceState(null, "", "/");
    expect(navigate("/data/index.json")).toBe(false);
    expect(navigate("https://example.org/datasets")).toBe(false);
    expect(here()).toBe("/");
  });
});
