// @vitest-environment happy-dom
// The datasets table reads every MIC constant from one summary file beside
// the index (baselines.json), pinned to the SHA-256 of the index bytes the
// page read. The index is the dev fixture's; the summary is synthetic. The
// load cache keys on the path, so the refused case has a file of its own
// (baselines-summary-other.test.ts).
import { act } from "react";
import { afterEach, describe, expect, it, vi } from "vitest";
import type { Root } from "react-dom/client";
import { summary, summaryEntry } from "../dev/fixtures/baselines";
import { indexSha, micRows, mountTable } from "../dev/fixtures/summaryPage";
import { pct } from "../data/format";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let root: Root | null = null;
afterEach(() => {
  act(() => root?.unmount());
  root = null;
  document.body.innerHTML = "";
  vi.unstubAllGlobals();
});

describe("the summary beside the index, pinned to its bytes (synthetic)", () => {
  it("is one request, and gives the constant of each dataset it covers", async () => {
    // Half the MIC datasets covered, with a constant of 0.5 and the row's
    // own model figure; one of them unconfirmed.
    const covered = micRows.filter((_, k) => k % 2 === 0);
    const entries = covered.map((r, k) => summaryEntry(r.id, r.mic_essential_agreement!, 0.5, k !== 0));
    const page = await mountTable(summary(await indexSha(), entries));
    root = page.root;
    expect(page.requested.filter((p) => p === "/data/baselines.json")).toHaveLength(1);
    expect(page.requested.filter((p) => p.startsWith("/data/baselines/"))).toEqual([]);
    const cells = [...document.querySelectorAll('[data-baseline="mic-constant"]')];
    expect(cells.length).toBeGreaterThan(0);
    expect(cells.length).toBeLessThanOrEqual(covered.length);
    for (const c of cells) expect(c.textContent).toContain(pct(0.5));
    expect(cells.filter((c) => c.textContent!.includes("†")).length).toBeLessThanOrEqual(1);
  });
});
