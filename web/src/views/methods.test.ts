// @vitest-environment happy-dom
import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { parseManifest } from "../data/validate";
import { MethodsView } from "./MethodsView";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const render = (k: number | null) => renderToStaticMarkup(createElement(MethodsView, { manifest: { ...manifest, k } }));

describe("methods features", () => {
  it("names the k the reports share, and says when they share none", () => {
    expect(render(21)).toContain("counted into its 21-mers");
    const unshared = render(null);
    expect(unshared).toContain("do not all record the same k");
    expect(unshared).not.toContain("its -mers");
  });
});

describe("methods in plain words", () => {
  const el = document.createElement("div");
  el.innerHTML = renderToStaticMarkup(createElement(MethodsView, { manifest }));

  it("opens with the study's six steps in order, each linking to its detail", () => {
    const steps = [...el.querySelectorAll('ol[aria-label="How the study works"] > li')];
    expect(steps.map((li) => li.querySelector("p")?.textContent)).toEqual([
      "Collect genomes and lab results",
      "Train learners",
      "Test on genomes held out",
      "Score the calls",
      "Count the errors against limits",
      "Read what the numbers support",
    ]);
    expect(steps.map((li) => li.querySelector("a")?.getAttribute("href"))).toEqual(["#m-labels", "#m-models", "#m-models", "#m-scores", "#m-scores", "#m-limits"]);
    for (const step of steps) expect(el.querySelector(step.querySelector("a")!.getAttribute("href")!)).not.toBeNull();
    expect(steps[4]!.textContent).toContain("a VME of at most 1.5% and an ME of at most 3%");
  });

  it("starts every definition with a plain line before its technical detail", () => {
    const defs = [...el.querySelectorAll("dl > div")];
    expect(defs.length).toBeGreaterThanOrEqual(19);
    for (const d of defs) {
      const [plain, detail] = [...d.querySelectorAll("dd > p")];
      expect(plain?.textContent?.length ?? 0).toBeGreaterThan(20);
      expect(detail?.textContent?.length ?? 0).toBeGreaterThan(20);
    }
    const lineage = defs.find((d) => d.querySelector("dt")?.textContent === "Lineage-held-out split");
    expect(lineage?.querySelector("dd > p")?.textContent).toContain("GRM always holds out close relatives (one lineage) together");
  });
});
