// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { DetailsButton } from "./Details";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const viewport = (narrow: boolean) => {
  window.matchMedia = ((q: string) =>
    ({
      matches: narrow && q.includes("max-width"),
      media: q,
      onchange: null,
      addEventListener: () => {},
      removeEventListener: () => {},
      addListener: () => {},
      removeListener: () => {},
      dispatchEvent: () => false,
    }) as MediaQueryList) as typeof window.matchMedia;
};

let host: HTMLDivElement;
let root: Root;
let elsewhere: ReturnType<typeof vi.fn>;

const draw = () =>
  act(() =>
    root.render(
      createElement(
        "div",
        null,
        createElement(DetailsButton, { label: "How BA is computed", children: "The mean of sensitivity and specificity." }),
        createElement(DetailsButton, { label: "What a lineage is", children: "Genomes within a Mash distance of 0.0001." }),
        createElement("button", { type: "button", onClick: elsewhere }, "Elsewhere"),
      ),
    ),
  );

const trigger = (name: string) => [...host.querySelectorAll("button")].find((b) => b.textContent === name)!;

beforeEach(() => {
  viewport(false);
  elsewhere = vi.fn();
  host = document.createElement("div");
  document.body.append(host);
  root = createRoot(host);
  draw();
});

afterEach(() => {
  act(() => root.unmount());
  host.remove();
  vi.restoreAllMocks();
});

describe("explicit details", () => {
  it("open one note at a time, right after its button, and leave focus on the button", () => {
    const a = trigger("How BA is computed");
    const b = trigger("What a lineage is");
    a.focus();
    act(() => a.click());
    expect(a.getAttribute("aria-expanded")).toBe("true");
    const note = a.nextElementSibling!;
    expect(note.getAttribute("role")).toBe("note");
    expect(note.id).toBe(a.getAttribute("aria-controls"));
    expect(note.textContent).toBe("The mean of sensitivity and specificity.");
    expect(document.activeElement).toBe(a);

    act(() => b.click());
    expect(a.getAttribute("aria-expanded")).toBe("false");
    expect(b.getAttribute("aria-expanded")).toBe("true");
    expect(host.querySelectorAll("[role='note']")).toHaveLength(1);
  });

  it("close on Escape with focus still on the button", () => {
    const a = trigger("How BA is computed");
    a.focus();
    act(() => a.click());
    act(() => {
      a.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true, cancelable: true }));
    });
    expect(a.getAttribute("aria-expanded")).toBe("false");
    expect(document.activeElement).toBe(a);
  });

  it("close on a tap elsewhere, and that tap does nothing else", () => {
    const a = trigger("How BA is computed");
    const other = trigger("Elsewhere");
    act(() => a.click());
    act(() => {
      other.dispatchEvent(new Event("pointerdown", { bubbles: true }));
      other.click();
    });
    expect(a.getAttribute("aria-expanded")).toBe("false");
    expect(elsewhere).not.toHaveBeenCalled();

    // The next tap is the reader's own again.
    vi.spyOn(performance, "now").mockReturnValue(performance.now() + 5_000);
    act(() => other.click());
    expect(elsewhere).toHaveBeenCalledTimes(1);
  });

  it("float under the button on a wide screen and open in the flow below 640px", () => {
    const a = trigger("How BA is computed");
    act(() => a.click());
    expect(a.nextElementSibling!.getAttribute("data-placement")).toBe("float");
    act(() => a.click());

    act(() => root.unmount());
    viewport(true);
    root = createRoot(host);
    draw();
    const narrow = trigger("How BA is computed");
    act(() => narrow.click());
    expect(narrow.nextElementSibling!.getAttribute("data-placement")).toBe("inline");
  });
});
