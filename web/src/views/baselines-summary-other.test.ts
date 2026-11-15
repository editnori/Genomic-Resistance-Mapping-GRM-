// @vitest-environment happy-dom
// A summary pinned to other index bytes is refused whole: the table shows
// the model alone. The index is the dev fixture's; the summary is synthetic.
import { act } from "react";
import { afterEach, describe, expect, it, vi } from "vitest";
import type { Root } from "react-dom/client";
import { summary, summaryEntry } from "../dev/fixtures/baselines";
import { indexSha, micRows, mountTable } from "../dev/fixtures/summaryPage";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

describe("a summary of another index (synthetic)", () => {
  it("is refused whole, and the table renders with the model alone", async () => {
    const entries = micRows.map((r) => summaryEntry(r.id, r.mic_essential_agreement!, 0.5));
    const own = await indexSha();
    const other = own.replace(/^./, (c) => (c === "0" ? "1" : "0"));
    const page = await mountTable(summary(other, entries));
    root = page.root;
    expect(page.requested.filter((p) => p === "/data/baselines.json")).toHaveLength(1);
    expect(document.querySelectorAll('[data-baseline="mic-constant"]')).toHaveLength(0);
    expect(document.querySelector("[role='region'][aria-label='Datasets']")).not.toBeNull();
  });
});
