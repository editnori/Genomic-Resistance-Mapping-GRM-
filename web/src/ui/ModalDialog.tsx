// Native modal root. showModal() owns focus containment and makes the page
// behind it inert; Escape arrives as a cancel event that the caller owns.

import { useLayoutEffect, useRef, type ReactNode } from "react";
import { createPortal } from "react-dom";

export const ModalDialog = ({
  label,
  labelledBy,
  onClose,
  className,
  children,
}: {
  label?: string;
  labelledBy?: string;
  onClose(): void;
  className?: string;
  children: ReactNode;
}) => {
  const ref = useRef<HTMLDialogElement | null>(null);
  // Layout effect: its cleanup runs while the dialog is still in the page,
  // so close() ends the modal state and focus can return to the opener.
  useLayoutEffect(() => {
    const dialog = ref.current;
    const opener = document.activeElement instanceof HTMLElement ? document.activeElement : null;
    if (dialog !== null && !dialog.open) dialog.showModal();
    return () => {
      dialog?.close();
      if (opener?.isConnected === true) opener.focus({ preventScroll: true });
    };
  }, []);
  return createPortal(
    <dialog
      ref={ref}
      aria-label={label}
      aria-labelledby={labelledBy}
      className={`modal-root ${className ?? ""}`}
      onCancel={(e) => {
        e.preventDefault();
        // Portaled dialogs still share a React event tree. Closing a child
        // must not also cancel its parent presentation or inspector.
        e.stopPropagation();
        onClose();
      }}
    >
      {children}
    </dialog>,
    document.body,
  );
};
