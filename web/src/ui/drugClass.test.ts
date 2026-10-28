// @vitest-environment happy-dom
// Drug-class marks: every curated class maps to one of four families, an
// unknown or empty class falls back to the neutral ring, the mark is never
// the only cue, and the family colours are the validated tokens in both
// themes. Controls: squared, with their states in ARIA rather than shape.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { renderToStaticMarkup } from "react-dom/server";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { byFamily, ClassGlyph, ClassKey, ClassName, classFamily, FAMILIES } from "./drugClass";
import { ACTION_BUTTON, GhostButton, Segmented } from "./primitives";
import { Dropdown } from "./Dropdown";
import { TopBar } from "../app/TopBar";
import { ThemeProvider } from "../theme/ThemeProvider";
import { DatasetsView } from "../views/DatasetsView";
import type { Navigator } from "../app/router";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

/** The 36 class names of the revision-6 publication's index.json. */
const CURATED = [
  "aminoglycoside", "quinolone", "cephalosporin (3rd gen)", "carbapenem", "tetracycline", "cephalosporin (1st-2nd gen)", "penicillin + inhibitor", "penicillin",
  "macrolide", "folate pathway combination", "cephalosporin (4th gen)", "phenicol", "trimethoprim", "cephalosporin + inhibitor", "lincosamide", "polymyxin",
  "glycopeptide", "monobactam", "sulfonamide", "anti-TB (mycolic acid)", "rifamycin", "lipopeptide", "oxazolidinone", "nitrofuran", "anti-TB (pyrazinamide)",
  "streptogramin", "fosfomycin", "glycylcycline", "anti-TB (riminophenazine)", "anti-TB (cycloserine)", "anti-TB (nitroimidazole)", "anti-TB (arabinogalactan)",
  "anti-TB (folate)", "cephalosporin (anti-MRSA)", "fusidane", "mupirocin",
];

const html = (node: unknown) => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(node as never);
  return el;
};

describe("drug-class families", () => {
  it("puts every β-lactam together and leaves rarer or unknown classes neutral", () => {
    const families = new Map(CURATED.map((c) => [c, classFamily(c)]));
    expect(CURATED.filter((c) => families.get(c) === "beta-lactam").sort()).toEqual(
      ["carbapenem", "cephalosporin (1st-2nd gen)", "cephalosporin (3rd gen)", "cephalosporin (4th gen)", "cephalosporin (anti-MRSA)", "cephalosporin + inhibitor", "monobactam", "penicillin", "penicillin + inhibitor"].sort(),
    );
    expect(CURATED.filter((c) => families.get(c) === "aminoglycoside")).toEqual(["aminoglycoside"]);
    expect(CURATED.filter((c) => families.get(c) === "quinolone")).toEqual(["quinolone"]);
    // Colour is spent on three families only, so 25 classes share the ring.
    expect(CURATED.filter((c) => families.get(c) === "other")).toHaveLength(25);
    for (const unknown of ["", "unclassified", "Novel Class X", null, undefined]) expect(classFamily(unknown)).toBe("other");
    expect(classFamily("Carbapenem")).toBe("beta-lactam");
  });

  it("orders a list by family, then count, then name", () => {
    const counts = new Map([["quinolone", 29], ["carbapenem", 15], ["penicillin", 11], ["macrolide", 9], ["aminoglycoside", 32], ["phenicol", 9]]);
    expect(byFamily([...counts.keys()], (c) => counts.get(c) ?? 0)).toEqual(["carbapenem", "penicillin", "aminoglycoside", "quinolone", "macrolide", "phenicol"]);
  });

  it("draws a distinct shape per family, hidden from assistive tech, beside the class name", () => {
    const shapes = FAMILIES.map((f) => html(createElement(ClassGlyph, { family: f.id })).querySelector("svg")!);
    expect(new Set(shapes.map((s) => s.getAttribute("data-shape"))).size).toBe(4);
    for (const s of shapes) {
      expect(s.getAttribute("aria-hidden")).toBe("true");
      expect(s.classList.contains("class-glyph")).toBe(true);
    }
    const named = html(createElement(ClassName, { drugClass: "cephalosporin (3rd gen)" }));
    expect(named.textContent).toBe("cephalosporin (3rd gen)");
    expect(named.querySelector('svg[data-class-family="beta-lactam"]')).not.toBeNull();
    const blank = html(createElement(ClassName, { drugClass: "" }));
    expect(blank.textContent).toBe("Class not recorded");
    expect(blank.querySelector('svg[data-shape="ring"]')).not.toBeNull();
    expect(html(createElement(ClassKey)).textContent).toBe("β-lactamsAminoglycosidesQuinolonesOther classes");
  });

  it("defines each family colour in both themes and draws marks in system ink under forced colors", async () => {
    // Vitest hands CSS imports back empty, so read the stylesheet itself;
    // the app program carries no Node types, hence the untyped import.
    const fs = (await import(/* @vite-ignore */ `node:${"fs"}`)) as { readFileSync(path: string, encoding: string): string };
    const cwd = (globalThis as unknown as { process: { cwd(): string } }).process.cwd();
    const css = fs.readFileSync(`${cwd}/src/styles/index.css`, "utf8");
    const block = (sel: string) => css.slice(css.indexOf(`${sel} {`), css.indexOf("\n}\n", css.indexOf(`${sel} {`)));
    const tokens = (b: string) => Object.fromEntries([...b.matchAll(/--class-([a-z-]+):\s*(#[0-9a-f]{6})/g)].map((m) => [m[1], m[2]]));
    const light = tokens(block(".theme-light"));
    const dark = tokens(block(".theme-dark"));
    // The validated set: dataviz slots blue, yellow and magenta, light and dark steps.
    expect(light).toEqual({ "beta-lactam": "#2a78d6", aminoglycoside: "#eda100", quinolone: "#e87ba4", other: "#8e8e93" });
    expect(dark).toEqual({ "beta-lactam": "#3987e5", aminoglycoside: "#c98500", quinolone: "#d55181", other: "#7a7a80" });
    const forced = css.slice(css.indexOf("@media (forced-colors: active)"));
    expect(forced).toMatch(/\.class-glyph\s*\{\s*forced-color-adjust: none;/);
    expect(forced).toMatch(/fill: CanvasText/);
    expect(forced).toMatch(/stroke: CanvasText/);
  });
});

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
  localStorage.clear();
});

const mount = async (node: unknown) => {
  const host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  await act(async () => root!.render(node as never));
  return host;
};

describe("drug classes in the catalog", () => {
  it("marks each row's class, keys the marks once, and lists the class filter by family with decorative marks", async () => {
    // A wide desktop, so the filters are the docked rail.
    vi.stubGlobal("matchMedia", (query: string) => ({ matches: query.includes("min-width: 1440px") || query.includes("prefers-reduced-motion"), media: query, addEventListener: () => {}, removeEventListener: () => {} }));
    window.history.replaceState(null, "", "/datasets");
    const manifest = parseManifest(fixtureFile("/data/index.json"));
    const nav: Navigator = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };
    const host = await mount(createElement(DatasetsView, { manifest, nav }));
    const rows = [...host.querySelectorAll("tbody tr")];
    expect(rows.length).toBeGreaterThan(0);
    for (const r of rows) {
      const sub = r.querySelector(".dataset-sub")!;
      const glyph = sub.querySelector("svg.class-glyph")!;
      // The class name is printed after the species and its separator.
      const name = (sub.textContent ?? "").split("· ").pop()!;
      expect(manifest.datasets.some((d) => d.drug_class === name)).toBe(true);
      expect(glyph.getAttribute("data-class-family")).toBe(classFamily(name));
    }
    expect(host.querySelectorAll("[data-class-key]")).toHaveLength(1);
    const combo = [...host.querySelectorAll<HTMLButtonElement>('[role="combobox"]')].find((b) => b.getAttribute("aria-label")?.startsWith("Drug class"))!;
    expect(combo.className).toContain("rounded-md");
    await act(async () => combo.click());
    const options = [...document.querySelectorAll('[role="listbox"][aria-label="Drug class"] [role="option"]')].slice(1);
    expect(options).toHaveLength(new Set(manifest.datasets.map((d) => d.drug_class)).size);
    const order = options.map((o) => FAMILIES.findIndex((x) => x.id === o.querySelector("svg")?.getAttribute("data-class-family")));
    expect(order).toEqual([...order].sort((a, b) => a - b));
    for (const o of options) {
      expect(o.querySelector("svg")!.closest('[aria-hidden="true"]')).not.toBeNull();
      expect(o.textContent?.trim().length).toBeGreaterThan(0);
    }
  });
});

describe("squared controls keep their states", () => {
  it("draws no control as a pill", () => {
    expect(ACTION_BUTTON).not.toContain("rounded-full");
    const ghost = html(createElement(GhostButton, { onClick: () => {}, children: "Clear" })).querySelector("button")!;
    expect(ghost.className).not.toContain("rounded-full");
    expect(ghost.className).toContain("ui-control-target");
  });

  it("marks the current segment with aria-pressed and moves it on click", async () => {
    const onSelect = vi.fn();
    const host = await mount(createElement(Segmented<"a" | "b">, { ariaLabel: "Mode", options: [{ id: "a", label: "All" }, { id: "b", label: "Core" }], value: "a", onSelect }));
    const group = host.querySelector('[role="group"]')!;
    expect(group.getAttribute("aria-label")).toBe("Mode");
    expect(group.className).not.toContain("rounded-full");
    expect(group.className).not.toContain("rounded-[14px]");
    const [a, b] = [...group.querySelectorAll("button")];
    expect(a!.getAttribute("aria-pressed")).toBe("true");
    expect(b!.getAttribute("aria-pressed")).toBe("false");
    for (const btn of [a!, b!]) expect(btn.className).not.toContain("rounded-full");
    await act(async () => b!.click());
    expect(onSelect).toHaveBeenCalledWith("b");
  });

  it("keeps the top bar's current view, Present and the theme toggle as squared, labelled, keyboard buttons", async () => {
    vi.stubGlobal("ResizeObserver", class {
      observe() {}
      disconnect() {}
    });
    const onNavigate = vi.fn();
    const host = await mount(createElement(ThemeProvider, null, createElement(TopBar, { active: "datasets", loading: false, onNavigate, onPresent: () => {}, presentDisabled: false })));
    const present = host.querySelector<HTMLButtonElement>('button[aria-label="Present"]')!;
    const theme = host.querySelector<HTMLButtonElement>('button[aria-label^="Switch to"]')!;
    for (const b of [present, theme]) {
      expect(b.className).toContain("rounded-md");
      expect(b.className).not.toContain("rounded-full");
      expect(b.className).toContain("focus-visible:outline");
    }
    await act(async () => theme.click());
    expect(host.querySelector('button[aria-label^="Switch to"]')?.getAttribute("aria-label")).toBe("Switch to light theme");
    // happy-dom lays nothing out, so the compact selector stands in for the row.
    const combo = host.querySelector<HTMLButtonElement>('nav [role="combobox"]');
    if (combo !== null) {
      expect(combo.className).toContain("rounded-md");
      await act(async () => combo.click());
      const species = [...document.querySelectorAll('[role="option"]')].find((o) => o.textContent?.includes("Species"))!;
      await act(async () => species.dispatchEvent(new MouseEvent("click", { bubbles: true })));
      expect(onNavigate).toHaveBeenCalledWith("species");
    }
  });

  it("opens a squared dropdown from the keyboard and closes it with Escape", async () => {
    const onSelect = vi.fn();
    const host = await mount(createElement(Dropdown<"x" | "y">, { ariaLabel: "Pick", options: [{ id: "x", label: "Ex" }, { id: "y", label: "Why" }], mode: { kind: "one", value: "x", onSelect } }));
    const combo = host.querySelector<HTMLButtonElement>('[role="combobox"]')!;
    expect(combo.className).toContain("rounded-md");
    await act(async () => combo.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowDown", bubbles: true })));
    expect(combo.getAttribute("aria-expanded")).toBe("true");
    await act(async () => combo.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true })));
    expect(combo.getAttribute("aria-expanded")).toBe("false");
  });
});
