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

describe("methods panels", () => {
  it("says what the blood-culture and CTX-M rows read, and that a gene is a marker", () => {
    const html = render(21);
    expect(html).toContain("Blood-culture and CTX-M rows");
    expect(html).toContain("also the carbapenemase genes KPC, NDM, VIM and IMP");
    expect(html).toContain("The CTX-M row counts CTX-M alone");
    expect(html).toContain("For aztreonam the blood-culture row adds KPC only");
    expect(html).toContain("it does not prove that the gene causes the isolate&#x27;s resistance");
  });
});
