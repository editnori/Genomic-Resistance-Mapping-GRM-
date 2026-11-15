// The datasets table mounted over the dev fixture's index, served as bytes,
// with a synthetic summary beside it (baselines.json), for the tests of the
// summary's pin. Records every path the page asks for.

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { vi } from "vitest";
import { fixtureFile } from "./files";
import { sha256Hex, useManifest } from "../../data/load";
import type { Manifest } from "../../data/schema";
import type { Navigator } from "../../app/router";
import { DatasetsView } from "../../views/DatasetsView";

/** The fixture index at evaluation revision 14, as bytes. */
export const indexText = JSON.stringify({ ...(fixtureFile("/data/index.json") as Manifest), evaluation_revision: 14 });
export const micRows = (JSON.parse(indexText) as Manifest).datasets.filter((d) => d.mic_essential_agreement !== null);
export const indexSha = async (): Promise<string> => (await sha256Hex(new TextEncoder().encode(indexText).buffer as ArrayBuffer))!;

const nav: Navigator = { route: { view: "datasets", species: null, dataset: null }, go: () => {}, openDataset: () => {}, closeDataset: () => {} };

const Page = () => {
  const m = useManifest();
  return m.status === "ready" ? createElement(DatasetsView, { manifest: m.data, nav }) : null;
};

/** Mounts the table with `summary` served beside the index; any other path
 *  answers 404. Gives the root and the paths requested. */
export const mountTable = async (summary: unknown): Promise<{ root: Root; requested: string[] }> => {
  const requested: string[] = [];
  vi.stubGlobal("fetch", async (path: string) => {
    requested.push(path);
    if (path === "/data/index.json") return new Response(indexText, { headers: { "Content-Type": "application/json" } });
    if (path === "/data/baselines.json") return new Response(JSON.stringify(summary), { headers: { "Content-Type": "application/json" } });
    return new Response("", { status: 404 });
  });
  const host = document.createElement("div");
  document.body.append(host);
  const root = createRoot(host);
  await act(async () => {
    root.render(createElement(Page));
  });
  for (let i = 0; i < 4; i++)
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });
  return { root, requested };
};
