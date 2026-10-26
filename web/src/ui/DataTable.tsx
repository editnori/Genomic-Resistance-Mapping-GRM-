// The one dense table: blue one- or two-word headers that never wrap, set
// at the top of their cells, each with its sort arrow inline on the side
// away from the values, so a header's edge is its column's edge; an
// optional unit line under a header (the bar a judged column is read
// against); cells that share the first line's baseline, however many lines
// they hold; judged cells as plain text with the verdict in words; 25 rows
// with a Show all toggle, and every match while a search is active. Columns
// that declare widths get a fixed layout: those widths hold through every
// sort, filter and Show all, and the one column without a width takes the
// rest, down to its floor, below which the table scrolls sideways. The
// first column can stay put during that scroll. A click anywhere on a row
// opens it, except on a control or at the end of a text selection.
// ArrowUp and ArrowDown walk the row names; Enter opens the focused row. A
// sort moves each row from its old place to its new one.

import { useLayoutEffect, useMemo, useRef, useState, type KeyboardEvent, type MouseEvent, type ReactNode } from "react";
import { ArrowDown, ArrowUp, ArrowUpDown } from "lucide-react";
import type { ReferenceBar } from "../data/metrics";
import { int } from "../data/format";
import { HoverTip, type TipLine } from "./HoverTip";
import { BarValue, FOCUS_RING, GhostButton, ScrollFrame } from "./primitives";

export interface Column<Row> {
  id: string;
  header: string;
  /** A second header line: the unit or the bar the column is judged against. */
  unit?: string;
  /** An optional hint on the header (at most about 40 words); never the
   *  only place a reading is explained. */
  tip?: TipLine[];
  align?: "left" | "right";
  /** Sort key; omitted columns do not sort. Nulls sort last either way. */
  sortValue?: (row: Row) => number | string | null;
  /** First direction when the column is chosen. */
  firstDir?: "asc" | "desc";
  /** A value judged against a bar shows as plain text with its verdict;
   *  null falls back to cell. */
  status?: (row: Row) => { value: number; bar: ReferenceBar } | null;
  cell(row: Row): ReactNode;
  className?: string;
  /** A fixed column width (any CSS length). Once any column has one, the
   *  table lays out fixed and a column without one takes the rest. */
  width?: string;
  /** The floor of the column that takes the rest; below it the table
   *  scrolls sideways instead of squeezing it. */
  minWidth?: string;
}

export interface SortState {
  id: string;
  dir: "asc" | "desc";
}

export const PAGE_SIZE = 25;

const compare = (a: number | string | null, b: number | string | null, dir: 1 | -1): number => {
  if (a === null && b === null) return 0;
  if (a === null) return 1;
  if (b === null) return -1;
  if (typeof a === "number" && typeof b === "number") return (a - b) * dir;
  return String(a).localeCompare(String(b)) * dir;
};

export const sortRows = <Row,>(rows: readonly Row[], key: (row: Row) => number | string | null, dir: "asc" | "desc"): Row[] =>
  [...rows].sort((x, y) => compare(key(x), key(y), dir === "asc" ? 1 : -1));

/** The fixed layout's widths, one per column (null takes the rest), and
 *  the table's floor; null when no column declares a width. */
export const tableGeometry = (columns: readonly { width?: string; minWidth?: string }[]): { cols: (string | null)[]; minWidth: string } | null =>
  columns.some((c) => c.width !== undefined)
    ? { cols: columns.map((c) => c.width ?? null), minWidth: `calc(${columns.map((c) => c.width ?? c.minWidth ?? "0px").join(" + ")})` }
    : null;

/** A click that lands on a control in the row does that control's job. */
const CONTROL = "button, a[href], input, select, textarea, label, summary, [role='button'], [role='link'], [role='checkbox'], [data-row-ignore]";

/** Whether a click on a row should open it: not on a control, and not the
 *  click that ends a text selection inside the row. */
export const rowClickOpens = (target: EventTarget | null, row: Element, selection: Selection | null): boolean => {
  if (target instanceof Element && target.closest(CONTROL) !== null) return false;
  if (selection !== null && !selection.isCollapsed && selection.toString().trim() !== "") {
    const anchor = selection.anchorNode;
    if (anchor !== null && row.contains(anchor)) return false;
  }
  return true;
};

/** The sticky identity column needs an opaque cell in the surface it sits
 *  on; sticky-first lets the scroll frame shade its edge. */
const STICKY = "sticky-first sticky left-0 z-[1] bg-paper";

/** The sticky cell is opaque, so it mixes the row's hover and focus tints
 *  into its own paper rather than letting them show through. */
const STICKY_TINT =
  "group-focus-within/row:bg-[color:color-mix(in_srgb,rgb(var(--blue))_5%,rgb(var(--paper)))] group-hover/row:bg-[color:color-mix(in_srgb,rgb(var(--blue))_3.5%,rgb(var(--paper)))]";

/** Every cell draws its own rule, so a moving row carries its rule with it.
 *  Cells share the first line's baseline. */
const CELL = "border-b border-hair px-2 py-2 align-baseline text-control";

const reducedMotion = (): boolean => typeof window !== "undefined" && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

export const DataTable = <Row,>({
  label,
  rows,
  columns,
  rowKey,
  onOpen,
  initialSort,
  onSortChange,
  initialShowAll = false,
  onShowAllChange,
  showAllMatches = false,
  stickyFirst = false,
  empty,
}: {
  label: string;
  rows: readonly Row[];
  columns: readonly Column<Row>[];
  rowKey(row: Row): string;
  /** Opens a row; the first column's cell becomes the button that does it. */
  onOpen?(row: Row): void;
  initialSort?: SortState;
  /** Told each new sort, for a view that keeps it across Back. */
  onSortChange?(sort: SortState): void;
  initialShowAll?: boolean;
  onShowAllChange?(all: boolean): void;
  /** An active search shows every match. */
  showAllMatches?: boolean;
  /** Keep the first column in view while the table scrolls sideways; only
   *  on the paper surface. */
  stickyFirst?: boolean;
  empty: ReactNode;
}) => {
  const [sort, setSort] = useState<SortState | null>(initialSort ?? null);
  const [all, setAll] = useState(initialShowAll);
  const bodyRef = useRef<HTMLTableSectionElement | null>(null);
  // Row tops before a sort, read just before the new order renders.
  const before = useRef<Map<string, number> | null>(null);
  const geometry = useMemo(() => tableGeometry(columns), [columns]);

  const sorted = useMemo(() => {
    const col = sort === null ? undefined : columns.find((c) => c.id === sort.id);
    return col?.sortValue === undefined || sort === null ? [...rows] : sortRows(rows, col.sortValue, sort.dir);
  }, [rows, columns, sort]);

  const showAll = all || showAllMatches;
  const visible = showAll ? sorted : sorted.slice(0, PAGE_SIZE);

  const readTops = (): Map<string, number> => {
    const tops = new Map<string, number>();
    bodyRef.current?.querySelectorAll<HTMLTableRowElement>("tr[data-row-key]").forEach((tr) => tops.set(tr.dataset.rowKey!, tr.getBoundingClientRect().top));
    return tops;
  };

  const onSort = (col: Column<Row>) => {
    if (col.sortValue === undefined) return;
    if (!reducedMotion()) before.current = readTops();
    const next: SortState = sort?.id === col.id ? { id: col.id, dir: sort.dir === "asc" ? "desc" : "asc" } : { id: col.id, dir: col.firstDir ?? "desc" };
    setSort(next);
    onSortChange?.(next);
  };

  // FLIP: each row that stayed on screen slides from its old top to its new
  // one on the slow clock; a row that arrives fades in where it lands.
  useLayoutEffect(() => {
    const old = before.current;
    before.current = null;
    if (old === null || bodyRef.current === null) return;
    bodyRef.current.querySelectorAll<HTMLTableRowElement>("tr[data-row-key]").forEach((tr) => {
      const was = old.get(tr.dataset.rowKey!);
      const now = tr.getBoundingClientRect().top;
      const options: KeyframeAnimationOptions = { duration: 320, easing: "cubic-bezier(0.2, 0.8, 0.2, 1)" };
      if (was === undefined) tr.animate([{ opacity: 0 }, { opacity: 1 }], options);
      else if (Math.abs(was - now) > 0.5) tr.animate([{ transform: `translateY(${was - now}px)` }, { transform: "none" }], options);
    });
  }, [sort]);

  const walk = (e: KeyboardEvent<HTMLButtonElement>) => {
    if (e.key !== "ArrowDown" && e.key !== "ArrowUp") return;
    const buttons = [...(bodyRef.current?.querySelectorAll<HTMLButtonElement>("[data-row-open]") ?? [])];
    const at = buttons.indexOf(e.currentTarget);
    const next = buttons[at + (e.key === "ArrowDown" ? 1 : -1)];
    if (next !== undefined) {
      e.preventDefault();
      next.focus();
    }
  };

  const nameHeader = columns[0];
  const rest = columns.slice(1);

  const header = (col: Column<Row>, first: boolean) => {
    const sortedHere = sort?.id === col.id;
    const right = col.align === "right";
    const glyph =
      col.sortValue === undefined ? null : (
        <span data-sort-glyph="" aria-hidden="true" className={`inline-flex w-3 shrink-0 justify-center self-center ${sortedHere ? "" : "opacity-50"}`}>
          {!sortedHere ? <ArrowUpDown size={12} strokeWidth={1.75} /> : sort.dir === "asc" ? <ArrowUp size={12} strokeWidth={2} /> : <ArrowDown size={12} strokeWidth={2} />}
        </span>
      );
    // The arrow is in the label's own line, on the side away from the
    // values: after a left-aligned label, before a right-aligned one. It
    // never leaves the cell and the label keeps the column's edge.
    const words = (
      <span className={`flex flex-col ${right ? "items-end" : "items-start"}`}>
        <span className="inline-flex items-baseline gap-1 whitespace-nowrap">
          {right ? glyph : null}
          <span data-th-label="">{col.header}</span>
          {right ? null : glyph}
        </span>
        {col.unit !== undefined ? <span className="whitespace-nowrap font-normal text-muted-text">{col.unit}</span> : null}
      </span>
    );
    const content =
      col.sortValue !== undefined ? (
        <button
          type="button"
          onClick={() => onSort(col)}
          className={`ui-control-target inline-flex max-w-full items-start rounded-sm text-tiny font-medium text-blue transition-colors duration-fast ease-out hover:text-blue-deep active:scale-[0.98] ${right ? "justify-end" : "justify-start"} ${FOCUS_RING}`}
        >
          {words}
        </button>
      ) : (
        <span className="inline-flex text-tiny font-medium text-blue">{words}</span>
      );
    return (
      <th
        key={col.id}
        scope="col"
        aria-sort={sortedHere ? (sort.dir === "asc" ? "ascending" : "descending") : undefined}
        className={`border-b border-hair px-2 py-2 align-top font-medium ${right ? "text-right" : "text-left"} ${first && stickyFirst ? STICKY : ""}`}
      >
        {col.tip !== undefined ? <HoverTip lines={col.tip}>{content}</HoverTip> : content}
      </th>
    );
  };

  const onRowClick = (e: MouseEvent<HTMLTableRowElement>, row: Row) => {
    if (onOpen === undefined || !rowClickOpens(e.target, e.currentTarget, window.getSelection())) return;
    // The row's own button takes focus first, so whatever opens remembers
    // this row, and focus comes back to it when that closes.
    e.currentTarget.querySelector<HTMLElement>("[data-row-open]")?.focus({ preventScroll: true });
    onOpen(row);
  };

  if (rows.length === 0) return <div className="border-y border-hair py-6">{empty}</div>;

  return (
    <div>
      <ScrollFrame label={label}>
        <table className="w-full border-separate border-spacing-0" style={geometry !== null ? { tableLayout: "fixed", minWidth: geometry.minWidth } : undefined}>
          {geometry !== null ? (
            <colgroup>
              {geometry.cols.map((w, i) => (
                <col key={columns[i]!.id} style={w === null ? undefined : { width: w }} />
              ))}
            </colgroup>
          ) : null}
          <thead>
            <tr>
              {nameHeader !== undefined ? header(nameHeader, true) : null}
              {rest.map((c) => header(c, false))}
            </tr>
          </thead>
          <tbody ref={bodyRef}>
            {visible.map((row) => {
              const key = rowKey(row);
              return (
                <tr
                  key={key}
                  data-row-key={key}
                  onClick={(e) => onRowClick(e, row)}
                  className={`group/row transition-colors duration-fast ease-out hover:bg-blue/[0.035] focus-within:bg-blue/[0.05] ${onOpen !== undefined ? "cursor-pointer" : ""}`}
                >
                  <td className={`${CELL} ${nameHeader?.className ?? ""} ${stickyFirst ? `${STICKY} ${STICKY_TINT}` : ""}`}>
                    {onOpen !== undefined ? (
                      <button
                        type="button"
                        data-row-open=""
                        data-row-key={key}
                        onClick={() => onOpen(row)}
                        onKeyDown={walk}
                        className={`ui-control-target block w-full rounded-sm text-left text-control font-medium text-ink transition-colors duration-fast ease-out hover:text-blue active:scale-[0.98] ${FOCUS_RING}`}
                      >
                        {nameHeader?.cell(row)}
                      </button>
                    ) : (
                      nameHeader?.cell(row)
                    )}
                  </td>
                  {rest.map((col) => {
                    const status = col.status?.(row) ?? null;
                    return (
                      <td key={col.id} className={`${CELL} tabular-nums text-ink ${col.align === "right" ? "text-right" : "text-left"} ${col.className ?? ""}`}>
                        {status !== null ? <BarValue value={status.value} bar={status.bar} /> : col.cell(row)}
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </ScrollFrame>
      {sorted.length > PAGE_SIZE && !showAllMatches ? (
        <div className="mt-3 flex items-center gap-3">
          <GhostButton
            onClick={() => {
              setAll(!showAll);
              onShowAllChange?.(!showAll);
            }}
          >
            {showAll ? `Show ${PAGE_SIZE}` : `Show all ${int(sorted.length)}`}
          </GhostButton>
        </div>
      ) : null}
    </div>
  );
};
