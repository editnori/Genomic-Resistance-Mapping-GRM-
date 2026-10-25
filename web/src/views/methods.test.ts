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
