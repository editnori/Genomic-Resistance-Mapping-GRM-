// The one dense table: blue one- or two-word headers that never wrap and
// line up with their values, a sort glyph slot that is always occupied, an
// optional unit line under a header (the bar a judged column is read
// against), middle-aligned rows with a hover and focus tint, judged cells as
// pills (filled meets, outlined misses), 25 rows with a Show all toggle, and
// every match while a search is active. The first column can stay put during
// horizontal scroll. ArrowUp and ArrowDown walk the row names; Enter opens
// the focused row. A sort moves each row from its old place to its new one.

import { useLayoutEffect, useMemo, useRef, useState, type KeyboardEvent, type ReactNode } from "react";
import { ArrowDown, ArrowUp, ArrowUpDown } from "lucide-react";
import type { ReferenceBar } from "../data/metrics";
import { int } from "../data/format";
import { HoverTip, type TipLine } from "./HoverTip";
import { BarPill, FOCUS_RING, GhostButton, ScrollFrame } from "./primitives";

export interface Column<Row> {
  id: string;
  header: string;
  /** A second header line: the unit or the bar the column is judged against. */
  unit?: string;
  /** Definition shown on the header's hover. */
  tip?: TipLine[];
  align?: "left" | "right";
  /** Sort key; omitted columns do not sort. Nulls sort last either way. */
  sortValue?: (row: Row) => number | string | null;
  /** First direction when the column is chosen. */
  firstDir?: "asc" | "desc";
  /** A value judged against a bar shows as a pill; null falls back to cell. */
  status?: (row: Row) => { value: number; bar: ReferenceBar } | null;
  cell(row: Row): ReactNode;
  className?: string;
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

/** The sticky identity column needs an opaque cell in the surface it sits on. */
const STICKY = "sticky left-0 z-[1] bg-paper";

/** The sticky cell is opaque, so it mixes the row's hover and focus tints
 *  into its own paper rather than letting them show through. */
const STICKY_TINT =
  "group-focus-within/row:bg-[color:color-mix(in_srgb,rgb(var(--blue))_5%,rgb(var(--paper)))] group-hover/row:bg-[color:color-mix(in_srgb,rgb(var(--blue))_3.5%,rgb(var(--paper)))]";

/** Every cell draws its own rule, so a moving row carries its rule with it. */
const CELL = "border-b border-hair px-2 py-2 text-control";

const reducedMotion = (): boolean => typeof window !== "undefined" && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

export const DataTable = <Row,>({
  label,
  rows,
  columns,
  rowKey,
  onOpen,
  initialSort,
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
  initialSort?: { id: string; dir: "asc" | "desc" };
  /** An active search shows every match. */
  showAllMatches?: boolean;
  /** Keep the first column in view while the table scrolls sideways; only
   *  on the paper surface. */
  stickyFirst?: boolean;
  empty: ReactNode;
}) => {
  const [sort, setSort] = useState(initialSort ?? null);
  const [all, setAll] = useState(false);
  const bodyRef = useRef<HTMLTableSectionElement | null>(null);
  // Row tops before a sort, read just before the new order renders.
  const before = useRef<Map<string, number> | null>(null);

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
    setSort((s) => (s?.id === col.id ? { id: col.id, dir: s.dir === "asc" ? "desc" : "asc" } : { id: col.id, dir: col.firstDir ?? "desc" }));
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
    const glyph = !sortedHere ? (
      <ArrowUpDown size={12} strokeWidth={1.75} aria-hidden="true" className="opacity-50" />
    ) : sort.dir === "asc" ? (
      <ArrowUp size={12} strokeWidth={2} aria-hidden="true" />
    ) : (
      <ArrowDown size={12} strokeWidth={2} aria-hidden="true" />
    );
    const right = col.align === "right";
    // The glyph sits outside the label's box, on the side away from the
    // column's values, so the label's edge lines up with the values below.
    const words = (
      <span className={`flex flex-col ${right ? "items-end" : "items-start"}`}>
        <span className="relative inline-flex items-center whitespace-nowrap">
          {col.header}
          {col.sortValue !== undefined ? (
            <span className={`absolute top-1/2 inline-flex w-3 -translate-y-1/2 justify-center ${right ? "right-full mr-1" : "left-full ml-1"}`}>{glyph}</span>
          ) : null}
        </span>
        {col.unit !== undefined ? <span className="whitespace-nowrap font-normal text-muted-text">{col.unit}</span> : null}
      </span>
    );
    const content =
      col.sortValue !== undefined ? (
        <button
          type="button"
          onClick={() => onSort(col)}
          className={`ui-control-target inline-flex items-end rounded-sm text-tiny font-medium text-blue transition-colors duration-fast ease-out hover:text-blue-deep active:scale-[0.98] ${FOCUS_RING}`}
        >
          {words}
        </button>
      ) : (
        <span className="text-tiny font-medium text-blue">{words}</span>
      );
    return (
      <th
        key={col.id}
        scope="col"
        aria-sort={sortedHere ? (sort.dir === "asc" ? "ascending" : "descending") : undefined}
        className={`border-b border-hair px-2 py-2 align-bottom font-medium ${right ? "text-right" : "text-left"} ${first && stickyFirst ? STICKY : ""}`}
      >
        {col.tip !== undefined ? (
          <HoverTip lines={col.tip} focusable={col.sortValue === undefined}>
            {content}
          </HoverTip>
        ) : (
          content
        )}
      </th>
    );
  };

  if (rows.length === 0) return <div className="border-y border-hair py-6">{empty}</div>;

  return (
    <div>
      <ScrollFrame label={label}>
        <table className="w-full border-separate border-spacing-0">
          <thead>
            <tr>
              {nameHeader !== undefined ? header(nameHeader, true) : null}
              {rest.map((c) => header(c, false))}
            </tr>
          </thead>
          <tbody ref={bodyRef}>
            {visible.map((row) => (
              <tr
                key={rowKey(row)}
                data-row-key={rowKey(row)}
                onClick={(e) => {
                  if (onOpen === undefined) return;
                  if ((e.target as HTMLElement).closest("button, a, input") !== null) return;
                  onOpen(row);
                }}
                className={`group/row align-middle transition-colors duration-fast ease-out hover:bg-blue/[0.035] focus-within:bg-blue/[0.05] ${onOpen !== undefined ? "cursor-pointer" : ""}`}
              >
                <td className={`${CELL} ${nameHeader?.className ?? ""} ${stickyFirst ? `${STICKY} ${STICKY_TINT}` : ""}`}>
                  {onOpen !== undefined ? (
                    <button
                      type="button"
                      data-row-open=""
                      onClick={() => onOpen(row)}
                      onKeyDown={walk}
                      className={`ui-control-target w-full rounded-sm text-left text-control font-medium text-ink transition-colors duration-fast ease-out hover:text-blue active:scale-[0.98] ${FOCUS_RING}`}
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
                      {status !== null ? <BarPill value={status.value} bar={status.bar} /> : col.cell(row)}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </ScrollFrame>
      {sorted.length > PAGE_SIZE && !showAllMatches ? (
        <div className="mt-3 flex items-center gap-3">
          <GhostButton onClick={() => setAll((v) => !v)}>{showAll ? `Show ${PAGE_SIZE}` : `Show all ${int(sorted.length)}`}</GhostButton>
        </div>
      ) : null}
    </div>
  );
};
