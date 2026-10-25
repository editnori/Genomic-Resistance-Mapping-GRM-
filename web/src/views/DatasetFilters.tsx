// The Datasets filters. From 1440px they are a docked 18rem rail beside the
// table; below that the same fields open in a sheet from the search row's
// Filters button, full width on a phone. The footer holds only Clear; the
// one result count lives in the search row.

import { useEffect, useState, type ReactNode } from "react";
import { SlidersHorizontal, X } from "lucide-react";
import { ModalDialog } from "../ui/ModalDialog";
import { FOCUS_RING, GhostButton } from "../ui/primitives";

/** Must match the docked-rail media query in styles/index.css. */
const RAIL_QUERY = "(min-width: 1440px)";

const readDocked = (): boolean => typeof window !== "undefined" && window.matchMedia(RAIL_QUERY).matches;

/** Whether the rail docks at this width, tracked. */
export const useDocked = (): boolean => {
  const [docked, setDocked] = useState(readDocked);
  useEffect(() => {
    const q = window.matchMedia(RAIL_QUERY);
    const onChange = (e: MediaQueryListEvent) => setDocked(e.matches);
    q.addEventListener("change", onChange);
    return () => q.removeEventListener("change", onChange);
  }, []);
  return docked;
};

export const FiltersButton = ({ activeCount, onOpen }: { activeCount: number; onOpen(): void }) => (
  <button
    type="button"
    onClick={onOpen}
    aria-haspopup="dialog"
    aria-label={activeCount > 0 ? `Filters, ${activeCount} set` : "Filters"}
    className={`ui-control-target inline-flex h-10 items-center gap-1.5 rounded-md border border-hair bg-paper-pure px-3 text-control font-medium text-ink transition-[background-color,transform] duration-fast ease-out hover:bg-surface-subtle active:scale-[0.98] ${FOCUS_RING}`}
  >
    <SlidersHorizontal size={14} strokeWidth={1.75} aria-hidden="true" />
    Filters
    {activeCount > 0 ? (
      <span aria-hidden="true" className="inline-flex h-5 min-w-5 items-center justify-center rounded-md bg-blue px-1 text-tiny tabular-nums text-paper-pure">
        {activeCount}
      </span>
    ) : null}
  </button>
);

/** One filter: a visible label above its control. */
export const FilterField = ({ label, children }: { label: string; children: ReactNode }) => (
  <div className="min-w-0">
    <div className="text-control font-medium text-ink">{label}</div>
    <div className="mt-1.5">{children}</div>
  </div>
);

/** A native checkbox row with the count choosing it would return. */
export const FilterCheck = ({ label, checked, count, onChange }: { label: string; checked: boolean; count: number; onChange(next: boolean): void }) => (
  <label className="flex min-h-11 cursor-pointer items-center gap-2.5 rounded-md px-1 text-control text-ink hover:bg-surface-subtle has-[:focus-visible]:outline has-[:focus-visible]:outline-2 has-[:focus-visible]:outline-blue">
    <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} className="h-4 w-4 shrink-0 accent-[rgb(var(--blue))]" />
    <span className="min-w-0 flex-1">{label}</span>
    <span className="text-tiny tabular-nums text-muted-text">{count}</span>
  </label>
);

export const FilterPanel = ({
  docked,
  open,
  onClose,
  activeCount,
  onClear,
  children,
}: {
  docked: boolean;
  /** The sheet's open state; the docked rail is always shown. */
  open: boolean;
  onClose(): void;
  activeCount: number;
  onClear(): void;
  children: ReactNode;
}) => {
  const footer =
    activeCount > 0 ? (
      <footer className="flex shrink-0 items-center justify-end border-t border-hair px-5 py-3">
        <GhostButton onClick={onClear}>Clear</GhostButton>
      </footer>
    ) : null;
  const body = (
    <div data-modal-scroll-owner="true" className="min-h-0 flex-1 space-y-5 overflow-y-auto overscroll-contain px-5 py-4">
      {children}
    </div>
  );
  if (docked) {
    return (
      <aside aria-label="Filters" data-filter-rail="docked" className="absolute inset-y-0 left-0 z-30 flex w-[var(--rail-width)] flex-col border-r border-hair bg-sidebar-surface text-ink">
        <h2 className="shrink-0 px-5 pt-5 text-body font-semibold text-ink">Filters</h2>
        {body}
        {footer}
      </aside>
    );
  }
  if (!open) return null;
  return (
    <ModalDialog labelledBy="filters-title" onClose={onClose}>
      <div className="absolute inset-0 bg-paper/60 backdrop-blur-[2px]" onClick={onClose} aria-hidden="true" />
      <div className="motion-popover-in absolute inset-0 flex flex-col overflow-hidden bg-sidebar-surface text-ink shadow-[var(--surface-shadow)] sm:inset-y-3 sm:left-auto sm:right-3 sm:w-[24rem] sm:rounded-xl sm:border sm:border-hair">
        <header className="flex shrink-0 items-center justify-between gap-3 border-b border-hair px-5 py-2">
          <h2 id="filters-title" className="text-body font-semibold text-ink">
            Filters
          </h2>
          <button type="button" onClick={onClose} aria-label="Close filters" className={`ui-control-target inline-flex h-10 w-10 items-center justify-center rounded-md hover:bg-surface-control ${FOCUS_RING}`}>
            <X size={18} aria-hidden="true" />
          </button>
        </header>
        {body}
        {footer}
      </div>
    </ModalDialog>
  );
};
