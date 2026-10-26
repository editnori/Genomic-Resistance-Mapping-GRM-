// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { HoverTip, placeTip } from "./HoverTip";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

let host: HTMLDivElement;
let root: Root;

beforeEach(() => {
  vi.useFakeTimers();
  host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  act(() => root.render(createElement(HoverTip, { lines: [{ term: "Gap", text: "Random minus lineage-held-out balanced accuracy." }], children: createElement("button", { type: "button" }, "Gap") })));
});

afterEach(() => {
  act(() => root.unmount());
  host.remove();
  vi.useRealTimers();
});

const pointer = (type: string, pointerType: string): Event => {
  const e = new MouseEvent(type, { bubbles: true });
  Object.defineProperty(e, "pointerType", { value: pointerType });
  return e;
};

const tip = () => document.querySelector("[role='tooltip']");

describe("optional hints", () => {
  it("add no tab stop of their own", () => {
    expect(host.querySelector("[tabindex]")).toBeNull();
    expect(host.querySelectorAll("button")).toHaveLength(1);
  });

  it("open for a resting mouse, never for a touch, and close on Escape", () => {
    const button = host.querySelector("button")!;
    act(() => {
      button.dispatchEvent(pointer("pointerover", "touch"));
      vi.advanceTimersByTime(400);
    });
    expect(tip()).toBeNull();
    act(() => {
      button.dispatchEvent(pointer("pointerout", "touch"));
      button.dispatchEvent(pointer("pointerover", "mouse"));
      vi.advanceTimersByTime(400);
    });
    expect(tip()?.textContent).toContain("Random minus lineage-held-out balanced accuracy.");
    act(() => {
      window.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true }));
    });
    expect(tip()).toBeNull();
  });

  it("stay inside the viewport: under the anchor, over it near the bottom, pulled in from the right edge", () => {
    const view = { width: 390, height: 844 };
    expect(placeTip({ left: 20, top: 100, bottom: 120 }, { width: 200, height: 80 }, view)).toEqual({ left: 20, top: 126 });
    expect(placeTip({ left: 20, top: 800, bottom: 820 }, { width: 200, height: 80 }, view)).toEqual({ left: 20, top: 714 });
    expect(placeTip({ left: 300, top: 100, bottom: 120 }, { width: 200, height: 80 }, view).left).toBe(182);
    const tall = placeTip({ left: 8, top: 40, bottom: 60 }, { width: 374, height: 830 }, view);
    expect(tall.top).toBeGreaterThanOrEqual(8);
  });
});
