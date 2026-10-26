// The report's one filter rail. From 1440px it is a permanent left rail in
// the view's grid: its title row sits on the search row's line and its
// first group starts level with the table header. Expanded it is the filter
// workspace; collapsed it is a 3.5rem strip carrying the toggle and, when
// filters are set, their count, and the toggle keeps its pixel in both
// states. Below 1440px the same groups open as a sheet from the search
// row's Filters button. The bottom band is always one height: the result
// count until a filter is set, then the filter count and Clear (the search
// row keeps the result count). The collapse preference is shared by every
// view with a rail.

import { useCallback, useEffect, useRef, useState, type ReactNode } from "react";
import { SlidersHorizontal } from "lucide-react";
import { ModalDialog } from "./ModalDialog";
import { FOCUS_RING, GhostButton } from "./primitives";

/** The width where the rail docks beside the table. */
const RAIL_QUERY = "(min-width: 1440px)";
const RAIL_OPEN_KEY = "grm-filter-rail-open";

const readInline = (): boolean => typeof window !== "undefined" && window.matchMedia(RAIL_QUERY).matches;

const readStoredOpen = (): boolean => {
  try {
    return localStorage.getItem(RAIL_OPEN_KEY) !== "0";
  } catch {
    return true;
  }
};

export interface RailState {
  inline: boolean;
  open: boolean;
  setOpen(open: boolean): void;
}

/** Docked: open follows the shared preference. Sheet: starts closed. */
export const useRailState = (): RailState => {
  const [inline, setInline] = useState(readInline);
  const [docked, setDocked] = useState(readStoredOpen);
  const [sheet, setSheet] = useState(false);
  useEffect(() => {
    const q = window.matchMedia(RAIL_QUERY);
    const onChange = (e: MediaQueryListEvent) => setInline(e.matches);
    q.addEventListener("change", onChange);
    return () => q.removeEventListener("change", onChange);
  }, []);
  const setOpen = useCallback(
    (open: boolean) => {
      if (inline) {
        setDocked(open);
        try {
          localStorage.setItem(RAIL_OPEN_KEY, open ? "1" : "0");
        } catch {
          // The rail still toggles for this visit without storage.
        }
      } else setSheet(open);
    },
    [inline],
  );
  return { inline, open: inline ? docked : sheet, setOpen };
};

/** The grid a view lays its rail and its content in, from 1440px. */
export const railGrid = (rail: RailState): string =>
  rail.inline ? `grid gap-x-8 ${rail.open ? "grid-cols-[var(--view-options-width)_minmax(0,1fr)]" : "grid-cols-[var(--view-options-collapsed-width)_minmax(0,1fr)]"}` : "";

/** The trigger that opens the sheet below 1440px; absent where the rail is
 *  permanent. It shares the view's search row. */
export const FiltersButton = ({ rail, activeCount }: { rail: RailState; activeCount: number }) =>
  rail.inline ? null : (
    <button
      type="button"
      onClick={() => rail.setOpen(true)}
      aria-expanded={rail.open}
      aria-haspopup="dialog"
      aria-label={activeCount > 0 ? `Filters, ${activeCount} set` : "Filters"}
      className={`ui-control-target inline-flex h-7 items-center gap-1.5 rounded-full border border-hair bg-paper-pure px-3 text-control font-medium text-ink transition-[color,background-color,border-color,transform] duration-fast ease-out hover:border-blue/35 hover:bg-surface-subtle active:scale-[0.98] ${FOCUS_RING}`}
    >
      <SlidersHorizontal size={13} strokeWidth={1.75} aria-hidden="true" />
      Filters
      {activeCount > 0 ? (
        <span aria-hidden="true" className="inline-flex h-5 min-w-5 items-center justify-center rounded-full bg-blue px-1 text-tiny tabular-nums text-paper-pure">
          {activeCount}
        </span>
      ) : null}
    </button>
  );

/** A standing question. Groups do not fold. */
export const RailGroup = ({ label, children }: { label: string; children: ReactNode }) => (
  <section className="filter-group">
    <h3 className="filter-group-label">{label}</h3>
    <div className="filter-group-body">{children}</div>
  </section>
);

/** One choice: an ink field label, a full-width 28px control, an optional hint. */
export const RailField = ({ label, active = false, hint, children }: { label: string; active?: boolean; hint?: string; children: ReactNode }) => (
  <div className="filter-field" data-active={active ? "true" : "false"}>
    <div className="filter-field-label">{label}</div>
    <div className="filter-field-control">{children}</div>
    {hint !== undefined ? <p className="filter-field-hint">{hint}</p> : null}
  </div>
);

const TITLE = "min-w-0 truncate text-body font-semibold leading-tight tracking-tight text-ink";

export const FilterRail = ({
  rail,
  activeCount,
  resultText,
  onClear,
  children,
}: {
  rail: RailState;
  /** Divergence from the arrival state; a many-of axis counts once. */
  activeCount: number;
  /** "237 datasets": the band's text until a filter is set. */
  resultText: string;
  onClear(): void;
  children: ReactNode;
}) => {
  const toggleRef = useRef<HTMLButtonElement | null>(null);

  // The band is one height in every state: the result count until a filter
  // is set, then the filter count and Clear.
  const band = (
    <footer className="flex h-11 shrink-0 items-center justify-between gap-3 border-t border-hair">
      {activeCount === 0 ? (
        <span className="text-control text-muted-text">{resultText}</span>
      ) : (
        <>
          <span className="text-control text-ink">
            {activeCount} {activeCount === 1 ? "filter" : "filters"}
          </span>
          <GhostButton onClick={onClear}>Clear</GhostButton>
        </>
      )}
    </footer>
  );

  if (rail.inline) {
    return (
      <aside
        aria-label="Filters"
        data-view-options-panel=""
        data-view-options-mode="inline"
        data-state={rail.open ? "open" : "closed"}
        onKeyDown={(e) => {
          if (rail.open && e.key === "Escape") {
            e.stopPropagation();
            rail.setOpen(false);
            toggleRef.current?.focus();
          }
        }}
        className="sticky top-4 flex max-h-[calc(100dvh-7rem)] min-h-0 min-w-0 flex-col self-start border-r border-hair"
      >
        <div className="relative flex h-10 shrink-0 items-center gap-1.5">
          <button
            ref={toggleRef}
            type="button"
            onClick={() => rail.setOpen(!rail.open)}
            aria-expanded={rail.open}
            aria-label={rail.open ? "Collapse filters" : activeCount === 0 ? "Expand filters" : `Expand filters, ${activeCount} set`}
            className={`-ml-2 flex size-9 shrink-0 items-center justify-center rounded-lg text-muted-text transition-colors duration-fast ease-out hover:bg-surface-control hover:text-ink active:scale-[0.98] ${FOCUS_RING}`}
          >
            <SlidersHorizontal size={15} strokeWidth={1.75} aria-hidden="true" />
          </button>
          {rail.open ? (
            <h2 className={`${TITLE} rail-in flex-1`}>Filters</h2>
          ) : activeCount === 0 ? null : (
            <span aria-hidden="true" className="pointer-events-none absolute left-4 top-0 inline-flex h-4 min-w-4 items-center justify-center rounded-full bg-blue px-1 text-tiny tabular-nums text-paper-pure">
              {activeCount}
            </span>
          )}
        </div>
        {rail.open ? (
          <div className="rail-in flex min-h-0 flex-1 flex-col pr-5">
            <div data-modal-scroll-owner="true" className="min-h-0 flex-1 overflow-y-auto overscroll-contain pb-4">
              <div className="filter-rail-body">{children}</div>
            </div>
            {band}
          </div>
        ) : null}
      </aside>
    );
  }

  if (!rail.open) return null;
  return (
    <ModalDialog label="Filters" onClose={() => rail.setOpen(false)}>
      <div className="absolute inset-0 bg-paper/60 backdrop-blur-[2px]" onClick={() => rail.setOpen(false)} aria-hidden="true" />
      <div
        data-view-options-panel=""
        data-view-options-mode="sheet"
        className="motion-sheet-in absolute inset-y-0 right-0 flex w-[26rem] flex-col overflow-hidden bg-sidebar-surface text-ink shadow-[var(--surface-shadow)] sm:inset-y-3 sm:right-3 sm:rounded-2xl sm:border sm:border-hair"
      >
        <header className="flex shrink-0 items-center gap-3 border-b border-hair px-5 py-2">
          <h2 className={`${TITLE} flex-1`}>Filters</h2>
          <GhostButton onClick={() => rail.setOpen(false)}>Close</GhostButton>
        </header>
        <div data-modal-scroll-owner="true" className="min-h-0 flex-1 overflow-y-auto overscroll-contain px-5 pb-4">
          <div className="filter-rail-body">{children}</div>
        </div>
        <div className="px-5">{band}</div>
      </div>
    </ModalDialog>
  );
};
