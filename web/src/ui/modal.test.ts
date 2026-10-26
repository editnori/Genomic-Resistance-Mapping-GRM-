// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot } from "react-dom/client";
import { expect, it, vi } from "vitest";
import { ModalDialog } from "./ModalDialog";
import { FilterRail } from "./FilterRail";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

it("keeps Escape in a nested dialog from closing its parent", () => {
  const host = document.createElement("div");
  document.body.append(host);
  const root = createRoot(host);
  const parentClose = vi.fn();
  const childClose = vi.fn();
  try {
    act(() => root.render(createElement(ModalDialog, {
      label: "Presentation", onClose: parentClose,
      children: createElement(ModalDialog, { label: "Sequence", onClose: childClose, children: "Sequence context" }),
    })));
    const child = document.querySelector('dialog[aria-label="Sequence"]')!;
    act(() => child.dispatchEvent(new Event("cancel", { bubbles: true, cancelable: true })));
    expect(childClose).toHaveBeenCalledOnce();
    expect(parentClose).not.toHaveBeenCalled();
  } finally {
    act(() => root.unmount());
    host.remove();
  }
});

it("renders the filter sheet as a nonmodal panel and returns focus on Escape", async () => {
  const host = document.createElement("div");
  const toggle = document.createElement("button");
  toggle.dataset.filterToggle = "";
  document.body.append(toggle);
  const root = createRoot(host);
  const setOpen = vi.fn();
  const rail = { inline: false, open: true, setOpen };
  try {
    act(() => root.render(createElement(FilterRail, {
      rail,
      activeCount: 0,
      resultText: "237 datasets",
      onClear: () => {},
      children: createElement("input", { type: "checkbox" }),
    })));
    const panel = host.querySelector("aside[aria-label='Filters']")!;
    expect(panel.getAttribute("role")).not.toBe("dialog");
    expect(host.querySelector("dialog")).toBeNull();

    await act(async () => {
      panel.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape", bubbles: true }));
      await new Promise(requestAnimationFrame);
    });
    expect(setOpen).toHaveBeenCalledWith(false);
    expect(document.activeElement).toBe(toggle);
  } finally {
    act(() => root.unmount());
    toggle.remove();
    host.remove();
  }
});

it("keeps the sheet's Close button focused after clearing its last filter", async () => {
  const host = document.createElement("div");
  document.body.append(host);
  const root = createRoot(host);
  const rail = { inline: false, open: true, setOpen: vi.fn() };
  const render = (activeCount: number) => createElement(FilterRail, {
    rail,
    activeCount,
    resultText: "237 datasets",
    onClear: () => act(() => root.render(render(0))),
    children: createElement("input", { type: "checkbox" }),
  });
  try {
    act(() => root.render(render(1)));
    const clear = [...host.querySelectorAll("button")].find((b) => b.textContent === "Clear")!;
    act(() => clear.click());
    expect(host.querySelector("button")?.textContent).toBe("Close");
    expect(document.activeElement).toBe(host.querySelector("button"));
  } finally {
    act(() => root.unmount());
    host.remove();
  }
});