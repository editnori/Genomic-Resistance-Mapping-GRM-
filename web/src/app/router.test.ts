// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { findFocusToken, flowPreviewHref, focusToken, isAppPath, LIMITS_HREF, navigate, parseFlowPreview, routeHref, useKeptState, useRouter, type Navigator } from "./router";

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
    expect(routeHref({ view: "data", species: null, dataset: null }, "m-funnel")).toBe("/corpus#m-funnel");
    expect(routeHref({ view: "methods", species: null, dataset: null }, "m-audit")).toBe("/methods#m-audit");
    expect(routeHref({ view: "datasets", species: null, dataset: "a b" }, "#x")).toBe("/datasets?dataset=a%20b#x");
    expect(routeHref({ view: "overview", species: null, dataset: null })).toBe("/");
    expect(LIMITS_HREF).toBe("/#limits");
    for (const p of ["/", "/datasets", "/species", "/species/escherichia_coli", "/corpus", "/methods", "/tests", "/flow", "/flow/genome-table", "/flow/genome-path", "/flow/bvbrc", "/flow/piece-count", "/flow/decision-tree"]) expect(isAppPath(p)).toBe(true);
    for (const p of ["/data/index.json", "/assets/index.js", "/favicon.svg", "/datasets/extra", "/flow/unknown", "/flow/bvbrc/extra"]) expect(isAppPath(p)).toBe(false);
    expect(parseFlowPreview("/flow")).toBe("all");
    expect(parseFlowPreview("/flow/decision-tree")).toBe("decision-tree");
    expect(parseFlowPreview("/flow/unknown")).toBeNull();
    expect(parseFlowPreview("/methods")).toBeNull();
    expect(flowPreviewHref("all")).toBe("/flow");
    expect(flowPreviewHref("genome-path")).toBe("/flow/genome-path");
  });

  it("keep the section hash while the inspector opens over a view and closes again", () => {
    window.history.replaceState(null, "", "/corpus#m-labels");
    const nav = withRouter();
    act(() => nav().openDataset("ecoli__amikacin"));
    expect(here()).toBe("/corpus?dataset=ecoli__amikacin#m-labels");
    expect(nav().route).toEqual({ view: "data", species: null, dataset: "ecoli__amikacin" });
    expect((window.history.state as { inspector?: boolean }).inspector).toBe(true);

    // An inspector reached by a link (no entry of its own) closes in place.
    window.history.replaceState(null, "", "/corpus?dataset=ecoli__amikacin#m-labels");
    act(() => nav().closeDataset());
    expect(here()).toBe("/corpus#m-labels");
  });

  it("go to a place and a section in it, and move legacy /methods data hashes to /corpus", () => {
    window.history.replaceState(null, "", "/");
    const nav = withRouter();
    act(() => nav().go({ view: "data", hash: "m-funnel" }));
    expect(here()).toBe("/corpus#m-funnel");
    expect(nav().route.view).toBe("data");
    act(() => {
      navigate("/methods#m-sources");
    });
    expect(here()).toBe("/corpus#m-sources");
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
