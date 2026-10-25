// @vitest-environment happy-dom
import { act, createElement } from "react";
import { createRoot } from "react-dom/client";
import { expect, it, vi } from "vitest";
import { ModalDialog } from "./ModalDialog";

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
