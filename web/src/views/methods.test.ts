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

  it("gives the outside-study test an entry that names its documents and links its Overview section", () => {
    const part = el.querySelector("#m-outside")!;
    expect(part.querySelector("h3")?.textContent).toBe("Newly deposited genomes");
    expect(part.textContent).toContain("so it is a retrospective test and not external validation.");
    expect([...part.querySelectorAll("code")].map((c) => c.textContent)).toEqual(["docs/outside-study-test.md", "docs/outside-study-results.md"]);
    expect(part.querySelector("a")?.getAttribute("href")).toBe("/#outside");
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

describe("methods panels", () => {
  const at = (evaluation_revision: number | undefined) => renderToStaticMarkup(createElement(MethodsView, { manifest: { ...manifest, evaluation_revision } }));
  const marker = "it does not prove that the gene causes the isolate&#x27;s resistance";

  it("from revision 13, says what the blood-culture and CTX-M rows read, and that a gene is a marker", () => {
    for (const revision of [13, 14]) {
      const html = at(revision);
      expect(html).toContain("one of its targets that the panel reads for the drug, in a species the assay reports that target for");
      expect(html).toContain("Blood-culture and CTX-M rows");
      expect(html).toContain("also the carbapenemase genes KPC, NDM, VIM and IMP");
      expect(html).toContain("The CTX-M row counts CTX-M alone");
      expect(html).toContain("For aztreonam the blood-culture row adds KPC only");
      expect(html).toContain("OXA-48-like is read for carbapenems only");
      expect(html).toContain(marker);
      expect(html).not.toContain("counts every target of the panel");
    }
  });

  // The site serves revision 12, whose blood-culture rows list all seven
  // targets on every drug, aztreonam included.
  it("before revision 13, says a panel row counts every target on every drug", () => {
    for (const revision of [undefined, 7, 12]) {
      const html = at(revision);
      expect(html).toContain("a panel calls resistant when it finds any gene it looks for");
      expect(html).toContain("counts every target of the panel on every drug the row is shown for, in every species");
      expect(html).toContain("Blood-culture and CTX-M rows");
      expect(html).toContain("CTX-M, the carbapenemase genes KPC, NDM, VIM, IMP and OXA-48-like, and mcr-1");
      expect(html).toContain("The CTX-M row counts CTX-M alone");
      expect(html).toContain("Reports from evaluation revision 13");
      expect(html).toContain(marker);
      expect(html).not.toContain("the panel reads for the drug, in a species");
      expect(html).not.toContain("For aztreonam the blood-culture row adds KPC only");
      expect(html).not.toContain("OXA-48-like is read for carbapenems only");
    }
  });

  it("names the best learner as the highest scorer among those scored on every isolate", () => {
    for (const revision of [12, 14]) expect(at(revision)).toContain("the learner that scored highest on the random split, among those scored on every isolate.");
  });
});
