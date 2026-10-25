// The one dense table: neutral one- or two-word headers that never wrap,
// the sort glyph after the label, an optional unit line under a header,
// right-aligned numbers, middle-aligned rows with 8px cell padding, status
// columns filled by their verdict, 25 rows with a Show all toggle, and every
// match while a search is active. The first column can stay put during
// horizontal scroll. ArrowUp and ArrowDown walk the row names; Enter opens
// the focused row.

import { useMemo, useRef, useState, type KeyboardEvent, type ReactNode } from "react";
import { ArrowDown, ArrowUp, ArrowUpDown } from "lucide-react";
import type { ReferenceBar } from "../data/metrics";
import { HoverTip, type TipLine } from "./HoverTip";
import { FOCUS_RING, GhostButton, ScrollFrame, StatusTd } from "./primitives";

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
  /** A value judged against a bar fills its cell; null falls back to cell. */
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

  const sorted = useMemo(() => {
    const col = sort === null ? undefined : columns.find((c) => c.id === sort.id);
    return col?.sortValue === undefined || sort === null ? [...rows] : sortRows(rows, col.sortValue, sort.dir);
  }, [rows, columns, sort]);

  const showAll = all || showAllMatches;
  const visible = showAll ? sorted : sorted.slice(0, PAGE_SIZE);

  const onSort = (col: Column<Row>) => {
    if (col.sortValue === undefined) return;
    setSort((s) => (s?.id === col.id ? { id: col.id, dir: s.dir === "asc" ? "desc" : "asc" } : { id: col.id, dir: col.firstDir ?? "desc" }));
  };

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
    const words = (
      <span className={`flex flex-col ${right ? "items-end" : "items-start"}`}>
        <span className="inline-flex items-center gap-1 whitespace-nowrap">
          {col.header}
          {col.sortValue !== undefined ? <span className="inline-flex w-3 shrink-0 justify-center">{glyph}</span> : null}
        </span>
        {col.unit !== undefined ? <span className="whitespace-nowrap text-tiny font-normal text-muted-text">{col.unit}</span> : null}
      </span>
    );
    const content =
      col.sortValue !== undefined ? (
        <button
          type="button"
          onClick={() => onSort(col)}
          className={`ui-control-target rounded-sm text-control font-medium text-muted-text transition-colors duration-fast ease-out hover:text-ink active:scale-[0.98] ${sortedHere ? "text-ink" : ""} ${FOCUS_RING}`}
        >
          {words}
        </button>
      ) : (
        <span className="text-control font-medium text-muted-text">{words}</span>
      );
    return (
      <th
        key={col.id}
        scope="col"
        aria-sort={sortedHere ? (sort.dir === "asc" ? "ascending" : "descending") : undefined}
        className={`px-2 py-2 align-bottom font-medium ${right ? "text-right" : "text-left"} ${first && stickyFirst ? STICKY : ""}`}
      >
        {col.tip !== undefined ? <HoverTip lines={col.tip} focusable={col.sortValue === undefined}>{content}</HoverTip> : content}
      </th>
    );
  };

  if (rows.length === 0) return <div className="border-y border-hair py-6">{empty}</div>;

  return (
    <div>
      <ScrollFrame label={label}>
        <table className="w-full border-collapse">
          <thead>
            <tr className="border-b border-hair">
              {nameHeader !== undefined ? header(nameHeader, true) : null}
              {rest.map((c) => header(c, false))}
            </tr>
          </thead>
          <tbody ref={bodyRef}>
            {visible.map((row) => (
              <tr
                key={rowKey(row)}
                onClick={(e) => {
                  if (onOpen === undefined) return;
                  if ((e.target as HTMLElement).closest("button, a, input") !== null) return;
                  onOpen(row);
                }}
                className={`border-b border-hair align-middle ${onOpen !== undefined ? "cursor-pointer" : ""}`}
              >
                <td className={`px-2 py-2 text-control ${nameHeader?.className ?? ""} ${stickyFirst ? STICKY : ""}`}>
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
                  return status !== null ? (
                    <StatusTd key={col.id} value={status.value} bar={status.bar} className={`py-2 text-control ${col.className ?? ""}`} />
                  ) : (
                    <td key={col.id} className={`px-2 py-2 text-control tabular-nums text-ink ${col.align === "right" ? "text-right" : "text-left"} ${col.className ?? ""}`}>
                      {col.cell(row)}
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
          <GhostButton onClick={() => setAll((v) => !v)}>{showAll ? `Show ${PAGE_SIZE}` : "Show all"}</GhostButton>
        </div>
      ) : null}
    </div>
  );
};
