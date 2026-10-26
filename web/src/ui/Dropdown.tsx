// One custom select for one-of and many-of choices. Closed it reads as a
// 28px control showing the current answer; open it is a paper listbox drawn
// in a layer anchored to the control's screen rectangle, opening upward when
// the room below cannot hold it, capped at twelve rows. Keyboard: arrows,
// Home/End, Enter or Space to choose, Escape or Tab to close, letters to jump.
// The rest option ("All species") is the first row and carries no count.

import { useCallback, useEffect, useId, useLayoutEffect, useRef, useState, type KeyboardEvent } from "react";
import { createPortal } from "react-dom";
import { Check, ChevronDown } from "lucide-react";
import { FOCUS_RING } from "./primitives";

export interface DropdownOption<T extends string> {
  id: T;
  label: string;
  count?: number;
  /** A value the current rows cannot serve: listed, but not choosable. */
  disabled?: boolean;
}

type Mode<T extends string> =
  | { kind: "one"; value: T; onSelect(id: T): void }
  | { kind: "many"; value: ReadonlySet<T>; onChange(next: Set<T>): void };

const ROW_PX = 32;
const MAX_ROWS = 12;
const MARGIN = 8;

interface Layer {
  left: number;
  top: number | null;
  bottom: number | null;
  width: number;
  maxHeight: number;
  host: HTMLElement;
}

const anchorLayer = (button: HTMLElement, rows: number): Layer => {
  const r = button.getBoundingClientRect();
  const vh = document.documentElement.clientHeight;
  const want = Math.min(rows, MAX_ROWS) * ROW_PX + 12;
  const below = vh - r.bottom - MARGIN;
  const above = r.top - MARGIN;
  const host = button.closest<HTMLElement>("dialog[open]") ?? document.body;
  const width = Math.max(r.width, 200);
  const left = Math.max(MARGIN, Math.min(r.left, document.documentElement.clientWidth - width - MARGIN));
  return below >= want || below >= above
    ? { left, top: r.bottom + 4, bottom: null, width, maxHeight: Math.min(want, below), host }
    : { left, top: null, bottom: vh - r.top + 4, width, maxHeight: Math.min(want, above), host };
};

export const Dropdown = <T extends string>({
  ariaLabel,
  options,
  rest,
  mode,
  fullWidth = false,
}: {
  ariaLabel: string;
  options: readonly DropdownOption<T>[];
  /** The rest option's words ("All species"); one-of lists may omit it. */
  rest?: string;
  mode: Mode<T>;
  /** Fill the field, as in the filter rail. */
  fullWidth?: boolean;
}) => {
  const listId = useId();
  const [open, setOpen] = useState(false);
  const [active, setActive] = useState(0);
  const [layer, setLayer] = useState<Layer | null>(null);
  const buttonRef = useRef<HTMLButtonElement | null>(null);
  const listRef = useRef<HTMLUListElement | null>(null);
  const typed = useRef({ text: "", at: 0 });

  // Row 0 is the rest option when present.
  const restRow = rest !== undefined ? 1 : 0;
  const rows = options.length + restRow;
  const isSelected = (id: T) => (mode.kind === "one" ? mode.value === id : mode.value.has(id));
  const restSelected = mode.kind === "many" ? mode.value.size === 0 : false;

  const closedLabel = (() => {
    if (mode.kind === "one") return options.find((o) => o.id === mode.value)?.label ?? rest ?? ariaLabel;
    const chosen = options.filter((o) => mode.value.has(o.id));
    if (chosen.length === 0) return rest ?? ariaLabel;
    return chosen.length === 1 ? chosen[0]!.label : `${chosen[0]!.label} +${chosen.length - 1}`;
  })();

  const close = useCallback(() => setOpen(false), []);

  const openList = () => {
    const first = mode.kind === "one" ? options.findIndex((o) => o.id === mode.value) : -1;
    setActive(first >= 0 ? first + restRow : 0);
    if (buttonRef.current !== null) setLayer(anchorLayer(buttonRef.current, rows));
    setOpen(true);
  };

  useLayoutEffect(() => {
    if (open && buttonRef.current !== null) setLayer(anchorLayer(buttonRef.current, rows));
  }, [open, rows]);

  useEffect(() => {
    if (!open) return undefined;
    const onDown = (e: PointerEvent) => {
      const t = e.target as Node;
      if (!buttonRef.current?.contains(t) && !listRef.current?.contains(t)) setOpen(false);
    };
    const onScroll = (e: Event) => {
      if (!listRef.current?.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener("pointerdown", onDown);
    window.addEventListener("scroll", onScroll, true);
    window.addEventListener("resize", close);
    return () => {
      document.removeEventListener("pointerdown", onDown);
      window.removeEventListener("scroll", onScroll, true);
      window.removeEventListener("resize", close);
    };
  }, [open, close]);

  useEffect(() => {
    if (open) listRef.current?.querySelector<HTMLElement>(`[data-row="${active}"]`)?.scrollIntoView({ block: "nearest" });
  }, [open, active]);

  const choose = (row: number) => {
    if (row < restRow) {
      if (mode.kind === "many") mode.onChange(new Set());
      close();
      buttonRef.current?.focus();
      return;
    }
    const option = options[row - restRow];
    if (option === undefined || option.disabled === true) return;
    if (mode.kind === "one") {
      mode.onSelect(option.id);
      close();
      buttonRef.current?.focus();
    } else {
      const next = new Set(mode.value);
      if (next.has(option.id)) next.delete(option.id);
      else next.add(option.id);
      mode.onChange(next);
    }
  };

  const enabledFrom = (from: number, step: 1 | -1): number => {
    for (let i = from; i >= 0 && i < rows; i += step) {
      if (i < restRow || options[i - restRow]?.disabled !== true) return i;
    }
    return Math.min(Math.max(from - step, 0), rows - 1);
  };

  const onKeyDown = (e: KeyboardEvent<HTMLButtonElement>) => {
    if (!open) {
      if (["ArrowDown", "ArrowUp", "Enter", " "].includes(e.key)) {
        e.preventDefault();
        openList();
      }
      return;
    }
    switch (e.key) {
      case "ArrowDown":
        e.preventDefault();
        setActive((i) => enabledFrom(Math.min(i + 1, rows - 1), 1));
        return;
      case "ArrowUp":
        e.preventDefault();
        setActive((i) => enabledFrom(Math.max(i - 1, 0), -1));
        return;
      case "Home":
        e.preventDefault();
        setActive(enabledFrom(0, 1));
        return;
      case "End":
        e.preventDefault();
        setActive(enabledFrom(rows - 1, -1));
        return;
      case "Enter":
      case " ":
        e.preventDefault();
        choose(active);
        return;
      case "Escape":
        // Consume the key so an enclosing dialog stays open.
        e.preventDefault();
        e.stopPropagation();
        close();
        return;
      case "Tab":
        close();
        return;
      default:
        if (e.key.length === 1 && !e.altKey && !e.ctrlKey && !e.metaKey) {
          const now = Date.now();
          typed.current = { text: now - typed.current.at > 500 ? e.key : typed.current.text + e.key, at: now };
          const q = typed.current.text.toLowerCase();
          const hit = options.findIndex((o) => o.disabled !== true && o.label.toLowerCase().startsWith(q));
          if (hit >= 0) setActive(hit + restRow);
        }
    }
  };

  const rowId = (row: number) => `${listId}-row-${row}`;

  return (
    <div className={`relative inline-block max-w-full ${fullWidth ? "w-full" : ""}`}>
      <button
        ref={buttonRef}
        type="button"
        role="combobox"
        aria-label={`${ariaLabel}: ${closedLabel}`}
        aria-haspopup="listbox"
        aria-expanded={open}
        aria-controls={open ? listId : undefined}
        aria-activedescendant={open ? rowId(active) : undefined}
        onClick={() => (open ? close() : openList())}
        onKeyDown={onKeyDown}
        className={`ui-control-target inline-flex h-7 max-w-full items-center gap-x-1.5 rounded-full border px-3 text-control font-medium tracking-tight transition-[color,background-color,border-color,transform] duration-fast ease-out active:scale-[0.98] ${fullWidth ? "w-full" : ""} ${FOCUS_RING} ${
          open ? "border-blue/60 bg-paper-pure text-blue" : "border-hair bg-paper-pure text-ink hover:border-blue/50 hover:text-blue"
        }`}
      >
        <span className="min-w-0 flex-1 truncate text-left">{closedLabel}</span>
        <ChevronDown size={14} strokeWidth={1.75} aria-hidden="true" className={`shrink-0 transition-transform duration-fast ease-out ${open ? "rotate-180" : ""}`} />
      </button>
      {open && layer !== null
        ? createPortal(
            <ul
              ref={listRef}
              id={listId}
              role="listbox"
              aria-label={ariaLabel}
              aria-multiselectable={mode.kind === "many" ? true : undefined}
              style={{
                position: "fixed",
                left: layer.left,
                top: layer.top ?? undefined,
                bottom: layer.bottom ?? undefined,
                minWidth: layer.width,
                maxHeight: layer.maxHeight,
              }}
              className="motion-popover-in scroll-thin z-[70] max-w-[min(88vw,26rem)] overflow-y-auto overscroll-contain rounded-xl border border-hair bg-paper-pure py-1 shadow-[var(--surface-shadow)]"
            >
              {rest !== undefined ? (
                <li
                  id={rowId(0)}
                  data-row={0}
                  role="option"
                  aria-selected={mode.kind === "one" ? false : restSelected}
                  onMouseDown={(e) => e.preventDefault()}
                  onClick={() => choose(0)}
                  onMouseEnter={() => setActive(0)}
                  className={`flex cursor-pointer items-center gap-x-3 px-3 py-1.5 text-control font-medium text-ink ${active === 0 ? "bg-blue/[0.08]" : ""}`}
                >
                  <span className="min-w-0 flex-1">{rest}</span>
                  <span className="inline-flex w-3.5 shrink-0 justify-end" aria-hidden="true">
                    {restSelected ? <Check size={14} strokeWidth={2} className="text-blue" /> : null}
                  </span>
                </li>
              ) : null}
              {options.map((o, i) => {
                const row = i + restRow;
                const selected = isSelected(o.id);
                return (
                  <li
                    key={o.id}
                    id={rowId(row)}
                    data-row={row}
                    role="option"
                    aria-selected={selected}
                    aria-disabled={o.disabled === true ? true : undefined}
                    onMouseDown={(e) => e.preventDefault()}
                    onClick={() => choose(row)}
                    onMouseEnter={() => {
                      if (o.disabled !== true) setActive(row);
                    }}
                    className={`flex items-center gap-x-3 px-3 py-1.5 text-control font-medium tracking-tight ${
                      o.disabled === true ? "cursor-not-allowed text-muted-text" : `cursor-pointer text-ink ${active === row ? "bg-blue/[0.08]" : ""}`
                    }`}
                  >
                    <span className="min-w-0 flex-1 break-words">{o.label}</span>
                    {o.count !== undefined ? <span className="text-tiny tabular-nums text-muted-text">{o.count}</span> : null}
                    <span className="inline-flex w-3.5 shrink-0 justify-end" aria-hidden="true">
                      {selected ? <Check size={14} strokeWidth={2} className="text-blue" /> : null}
                    </span>
                  </li>
                );
              })}
            </ul>,
            layer.host,
          )
        : null}
    </div>
  );
};
