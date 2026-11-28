// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it } from "vitest";
import {
  BvbrcFlow,
  DecisionTreeFlow,
  FLOW_DATA,
  FLOW_SLUGS,
  GenomePathFlow,
  GenomeTableFlow,
  PieceCountFlow,
  fillTemplate,
  formatCount,
  type FlowData,
} from "./index";
import { FlowPreviewView } from "../views/FlowPreviewView";

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
});

const cloneFlow = (): FlowData => JSON.parse(JSON.stringify(FLOW_DATA)) as FlowData;

describe("flow.json data integrity", () => {
  it("resolves every placeholder across all views and nodes without leftovers", () => {
    expect(Object.keys(FLOW_DATA.counts).length).toBe(136);
    expect(FLOW_SLUGS).toHaveLength(5);

    for (const view of Object.values(FLOW_DATA.views)) {
      for (const s of [view.kicker, view.title, view.sub, ...Object.values(view.ui)]) {
        expect(fillTemplate(FLOW_DATA, s)).not.toMatch(/\{[a-z0-9_]+\}/);
      }
    }
    for (const node of FLOW_DATA.nodes) {
      expect(fillTemplate(FLOW_DATA, node.title)).not.toMatch(/\{[a-z0-9_]+\}/);
      expect(fillTemplate(FLOW_DATA, node.line)).not.toMatch(/\{[a-z0-9_]+\}/);
      for (const [, text] of node.fields ?? []) {
        expect(fillTemplate(FLOW_DATA, text)).not.toMatch(/\{[a-z0-9_]+\}/);
      }
    }
  });

  it("formats each count format mode consistently", () => {
    expect(formatCount(FLOW_DATA, "records")).toBe("1,285,111");
    expect(formatCount(FLOW_DATA, "k")).toBe("31");
    expect(formatCount(FLOW_DATA, "median_patterns")).toBe("479,983");
  });
});

describe("five React flow components", () => {
  it("renders GenomeTableFlow and steps through all five steps cleanly", () => {
    act(() => root.render(createElement(GenomeTableFlow, { calm: true })));
    expect(host.textContent).toContain("A 31-letter window turns each genome into yes and no");
    expect(host.textContent).not.toMatch(/\{[a-z0-9_]+\}/);

    const stepButtons = [...host.querySelectorAll<HTMLButtonElement>("button[data-flow-step]")];
    expect(stepButtons).toHaveLength(5);
    expect(stepButtons[1]!.getAttribute("aria-current")).toBe("step");

    act(() => stepButtons[4]!.click());
    expect(stepButtons[4]!.getAttribute("aria-current")).toBe("step");
    expect(host.textContent).toContain("lab's call");
  });

  it("renders GenomePathFlow, switches active node, and reacts to single-file edits", () => {
    act(() => root.render(createElement(GenomePathFlow, { calm: true })));
    expect(host.textContent).toContain("9 comparisons along the genome's path, one still a candidate");
    expect(host.querySelectorAll("button[data-path-node]")).toHaveLength(FLOW_DATA.views.genome_path.walk.length);

    const lastBtn = host.querySelector<HTMLButtonElement>('button[data-path-node="last"]')!;
    act(() => lastBtn.click());
    expect(lastBtn.getAttribute("aria-pressed")).toBe("true");
    expect(host.textContent).toContain("252,148 later lab records");

    // Verify one edit to flow.json updates both GenomePathFlow and DecisionTreeFlow.
    const edited = cloneFlow();
    const lastNode = edited.nodes.find((n) => n.id === "last")!;
    lastNode.line = "updated under ruling {n_comparisons}";
    const mrunNode = edited.nodes.find((n) => n.id === "mrun")!;
    mrunNode.line = "intermediate datasets updated";

    act(() => root.render(createElement(GenomePathFlow, { flow: edited, startNode: "last", calm: true })));
    expect(host.textContent).toContain("updated under ruling 9");
    expect(host.textContent).toContain("intermediate datasets updated");

    act(() => root.render(createElement(DecisionTreeFlow, { flow: edited, startChoice: "last" })));
    expect(host.textContent).toContain("updated under ruling 9");
    expect(host.textContent).toContain("intermediate datasets updated");
  });

  it("renders BvbrcFlow across all five steps and expands step choices", () => {
    act(() => root.render(createElement(BvbrcFlow, { calm: true })));
    expect(host.textContent).toContain("BV-BRC sends 1,285,111 lab records; 611,826 pairs carry a call");
    expect(host.textContent).toContain("One lab record, by what we do with each field");

    const stepButtons = [...host.querySelectorAll<HTMLButtonElement>("button[data-flow-step]")];
    expect(stepButtons).toHaveLength(5);
    for (const btn of stepButtons) {
      act(() => btn.click());
      expect(host.textContent).not.toMatch(/\{[a-z0-9_]+\}/);
    }
  });

  it("renders PieceCountFlow across all seven steps without unresolved placeholders", () => {
    act(() => root.render(createElement(PieceCountFlow, { calm: true })));
    expect(host.textContent).toContain("Counting turns DNA letters into 479,983 yes-or-no patterns");

    const stepButtons = [...host.querySelectorAll<HTMLButtonElement>("button[data-flow-step]")];
    expect(stepButtons).toHaveLength(7);
    for (const btn of stepButtons) {
      act(() => btn.click());
      expect(host.textContent).not.toMatch(/\{[a-z0-9_]+\}/);
    }
    // Final step renders proportional pattern bars.
    expect(host.textContent).toContain("479,983");
    expect(host.textContent).toContain("4,000");
  });

  it("renders DecisionTreeFlow with all trunk choices and steps via Back/Next", () => {
    act(() => root.render(createElement(DecisionTreeFlow)));
    expect(host.textContent).toContain("9 of 22 choices had no reason and no measure");
    const totalRows = FLOW_DATA.views.decisions.trunk.reduce((acc, [, rows]) => acc + rows.length, 0);
    const choiceBtns = [...host.querySelectorAll<HTMLButtonElement>("button[data-choice-id]")];
    expect(choiceBtns).toHaveLength(totalRows);

    const openBefore = host.querySelector<HTMLButtonElement>('button[data-choice-id][aria-expanded="true"]')!;
    expect(openBefore.dataset.choiceId).toBe("as06");

    const nextBtn = [...host.querySelectorAll<HTMLButtonElement>("button")].find((b) => b.textContent === "Next")!;
    act(() => nextBtn.click());
    const openAfter = host.querySelector<HTMLButtonElement>('button[data-choice-id][aria-expanded="true"]')!;
    expect(openAfter.dataset.choiceId).not.toBe("as06");
  });

  it("renders FlowPreviewView for all five and for a single slug", () => {
    let picked = "";
    act(() => root.render(createElement(FlowPreviewView, { target: "all", onSelect: (next) => (picked = next) })));
    expect(host.querySelectorAll("[data-flow-diagram]")).toHaveLength(5);

    const pathTab = [...host.querySelectorAll<HTMLButtonElement>('[role="group"] button')].find((b) => b.textContent === "Path")!;
    act(() => pathTab.click());
    expect(picked).toBe("genome-path");

    act(() => root.render(createElement(FlowPreviewView, { target: "decision-tree", onSelect: () => {} })));
    expect(host.querySelectorAll("[data-flow-diagram]")).toHaveLength(1);
    expect(host.querySelector('[data-flow-diagram="decision-tree"]')).not.toBeNull();
  });
});
