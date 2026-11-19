// @vitest-environment happy-dom
// Whose error limits these are. The VME, ME and MIC agreement limits come
// from FDA's 2009 guidance for growth-based susceptibility test systems,
// which does not cover tests that read resistance from a genome or tests of
// anti-mycobacterial drugs. Methods must say so in its definition of VME and
// ME, and no view may name FDA without growth-based in the same sentence or
// hint, so a limit never reads as one these models were cleared against.

import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { BARS_SCOPE, BARS_SOURCE } from "../data/metrics";
import type { Navigator } from "../app/router";
import { wordingSummary } from "../data/derive/evidence";
import { barsLead, FINDINGS } from "./findings";
import { datasetColumns } from "./DatasetsView";
import { MethodsView } from "./MethodsView";
import { OverviewView } from "./OverviewView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const nav: Navigator = { route: { view: "overview", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };

const textOf = (html: string): string => {
  const el = document.createElement("div");
  el.innerHTML = html;
  return el.textContent ?? "";
};

/** Sentences of running text: a stop, then space, then a capital. */
const sentences = (text: string): string[] =>
  text
    .replace(/&apos;/g, "'")
    .replace(/\s+/g, " ")
    .split(/(?<=[.!?])\s+(?=[A-Z])/);

/** Each sentence that names FDA and does not say growth-based. */
const unscoped = (text: string): string[] => sentences(text).filter((s) => /FDA/.test(s) && !/growth-based/.test(s));

/** Wording the site dropped: a limit read as a device's or a regulator's. */
const OVERCLAIM = /FDA[- ]style|test devices|lab devices|held to a|regulatory claim|meets? FDA|FDA criteria|FDA[- ]cleared|FDA[- ]approved/i;

describe("the Methods definition of VME and ME", () => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(createElement(MethodsView, { manifest }));
  const dt = [...el.querySelectorAll("dt")].find((d) => d.textContent === "VME and ME");
  const def = dt?.nextElementSibling?.textContent?.replace(/\s+/g, " ") ?? "";

  it("says whose limits they are and what the guidance does not cover", () => {
    expect(dt).toBeDefined();
    expect(def).toContain("limits taken from FDA's 2009 guidance for growth-based susceptibility test systems");
    expect(def).toContain(
      "The guidance does not cover tests that read resistance from a genome, or tests of anti-mycobacterial drugs, so this site uses its numbers as a research screen adapted from it.",
    );
    expect(def).toContain("They are not a criterion these models were validated against, and meeting them is not clearance of anything.");
    expect(def).toContain(BARS_SCOPE);
  });

  it("says the limits are a comparison only for M. tuberculosis, and how the site's rule differs from the guidance's", () => {
    expect(def).toContain("M. tuberculosis datasets test anti-mycobacterial drugs, which the guidance excludes by name, so for them the limits are a comparison only.");
    expect(def).toContain("its lower end under 1.5% and its upper end under 7.5%");
    expect(def).toContain("it has no category for a test that gives no call");
    expect(def).toContain("The site counts each rate over the isolates a learner called");
    expect(def).toContain("it counts ME over every susceptible isolate tested and asks that rate be under 3%, and it has no category");
  });

  it("counts the best learner by the guidance's own wording, and says that reading came after the results", () => {
    const w = wordingSummary(manifest.datasets, "best", "lineage");
    expect(w.n).toBeGreaterThan(0);
    expect(def).toContain(
      `Read by that wording, over the isolates each learner called, the best learner meets on ${w.meets} of ${w.n} datasets with lineages held out. Overview gives the count for each learner. The reading was added after the results were read.`,
    );
  });
});

describe("the Overview's first reading against the limits", () => {
  it("says what the guidance does not cover before any count", () => {
    const lead = barsLead({ n: 268, meets: 0, misses: 239, open: 29, openAtBars: 6 });
    expect(lead.indexOf(BARS_SCOPE)).toBeGreaterThan(0);
    expect(lead.indexOf(BARS_SCOPE)).toBeLessThan(lead.indexOf("shown to meet"));
  });
});

describe("FDA is named only with growth-based", () => {
  it("in the rendered Overview, primer and Methods", () => {
    for (const html of [renderToStaticMarkup(createElement(OverviewView, { manifest, nav })), renderToStaticMarkup(createElement(MethodsView, { manifest }))]) {
      const text = textOf(html);
      expect(text).toMatch(/FDA/);
      expect(unscoped(text)).toEqual([]);
      expect(text).not.toMatch(OVERCLAIM);
    }
  });

  it("in every Present slide's headline and caption", () => {
    for (const build of FINDINGS) {
      const f = build(manifest, { page: false });
      const text = `${f.result} ${typeof f.caption === "string" ? f.caption : ""}`;
      expect(unscoped(text)).toEqual([]);
      expect(text).not.toMatch(OVERCLAIM);
    }
  });

  it("in each Datasets column hint, which names the source wherever a header prints a limit", () => {
    expect(BARS_SOURCE).toMatch(/FDA.*growth-based/);
    for (const c of datasetColumns(false, undefined)) {
      const tip = (c.tip ?? []).map((l) => `${l.term} ${l.text}`).join(" ");
      expect(unscoped(tip)).toEqual([]);
      expect(tip).not.toMatch(OVERCLAIM);
      if (c.status !== undefined) expect(tip).toContain(BARS_SOURCE);
    }
  });

  it("in every source file's text, comments included", () => {
    // Every app source as text; Vite reads them, so no Node types are needed.
    const sources = import.meta.glob<string>(["../**/*.{ts,tsx}", "!../**/*.test.{ts,tsx}"], { query: "?raw", import: "default", eager: true });
    const all = Object.entries(sources);
    expect(all.length).toBeGreaterThan(20);
    expect(all.map(([p]) => p).filter((p) => /MethodsView\.tsx$|metrics\.ts$/.test(p))).toHaveLength(2);
    const offending = all.flatMap(([p, text]) => unscoped(text).map((s) => `${p}: ${s.trim().slice(0, 160)}`));
    expect(offending).toEqual([]);
    expect(all.filter(([, text]) => OVERCLAIM.test(text)).map(([p]) => p)).toEqual([]);
  });
});
