// A chart's table alternative: every row of the figure, sortable, with the
// figure's selection shown on its row and set from it. Unlike the shared
// DataTable it never pages, so a mark pinned in the figure always has a row
// to show; instead it scrolls inside its own frame, under a sticky header,
// and brings the selected row into that frame without scrolling the page.
// The first cell is the row's button (Enter or Space pins it, ArrowUp and
// ArrowDown walk the rows); a click anywhere else on the row pins it too.
// A long table (2,000 genomes by 40 drugs) paints its first rows at once and
// adds the rest in short tasks, busy until every row is in, so opening it
// never holds the page for a second.

import { memo, useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState, type KeyboardEvent, type ReactElement, type ReactNode } from "react";
import { ArrowDown, ArrowUp, ArrowUpDown } from "lucide-react";
import { FOCUS_RING } from "../ui/primitives";

export interface LinkedColumn<Row> {
  id: string;
  header: string;
  /** A second header line: a unit, or a count the column carries. */
  unit?: string;
  align?: "left" | "right";
  /** A fixed CSS width. Give every column one, or none: with widths the
   *  table keeps them through sorting and filtering and scrolls sideways
   *  when they outgrow the frame. */
  width?: string;
  /** Sort key; omitted columns do not sort. Nulls sort last either way. */
  sortValue?: (row: Row) => number | string | null;
  firstDir?: "asc" | "desc";
  cell(row: Row): ReactNode;
}

const compare = (a: number | string | null, b: number | string | null, dir: 1 | -1): number => {
  if (a === null && b === null) return 0;
  if (a === null) return 1;
  if (b === null) return -1;
  if (typeof a === "number" && typeof b === "number") return (a - b) * dir;
  return String(a).localeCompare(String(b)) * dir;
};

/** About how many cells one task adds while a long table fills in. */
const CHUNK_CELLS = 4000;

const CELL = "border-b border-hair px-2 py-1.5 text-control";
/** The pinned row's edge: a blue inset bar, and an outline where forced
 *  colors drop both the bar and the row's tint. */
const PINNED = "shadow-[inset_3px_0_0_rgb(var(--blue))] forced-colors:outline forced-colors:outline-2 forced-colors:-outline-offset-2";
const HEAD = "sticky top-0 z-[2] border-b border-hair bg-paper px-2 pb-1.5 pt-2 align-bottom text-tiny font-medium text-blue";

interface BodyRowProps<Row> {
  row: Row;
  rowId: string;
  columns: readonly LinkedColumn<Row>[];
  on: boolean;
  stickyFirst: boolean;
  label?: string;
  onToggle(row: Row, on: boolean): void;
  onWalk(e: KeyboardEvent<HTMLButtonElement>): void;
}

/** One body row; memoized, so filling in or pinning redraws only the rows
 *  that changed. */
const BodyRow = memo(function BodyRow<Row>({ row, rowId, columns, on, stickyFirst, label, onToggle, onWalk }: BodyRowProps<Row>) {
  const [first, ...rest] = columns;
  return (
    <tr
      data-key={rowId}
      onClick={(e) => {
        if ((e.target as HTMLElement).closest("button, a, input") !== null) return;
        if ((window.getSelection()?.toString() ?? "") !== "") return;
        onToggle(row, on);
      }}
      className={`group/row cursor-pointer transition-colors duration-fast ease-out ${on ? "bg-blue/[0.08]" : "hover:bg-blue/[0.035]"}`}
    >
      <td className={`${CELL} ${on ? PINNED : ""} ${stickyFirst ? `sticky left-0 z-[1] ${on ? "bg-[color:color-mix(in_srgb,rgb(var(--blue))_8%,rgb(var(--paper)))]" : "bg-paper"}` : ""}`}>
        <button
          type="button"
          data-row-open=""
          data-row-key={rowId}
          aria-pressed={on}
          aria-label={label}
          onClick={() => onToggle(row, on)}
          onKeyDown={onWalk}
          className={`ui-control-target w-full rounded-sm text-left text-control text-ink transition-colors duration-fast ease-out hover:text-blue active:scale-[0.98] ${on ? "font-semibold" : "font-medium"} ${FOCUS_RING}`}
        >
          {first?.cell(row)}
        </button>
      </td>
      {rest.map((c) => (
        <td key={c.id} className={`${CELL} tabular-nums text-ink ${c.align === "right" ? "text-right" : "text-left"}`}>
          {c.cell(row)}
        </td>
      ))}
    </tr>
  );
}) as <Row>(props: BodyRowProps<Row>) => ReactElement;

export const LinkedTable = <Row,>({
  label,
  rows,
  columns,
  rowKey,
  rowLabel,
  selectedKey,
  onSelect,
  initialSort,
  maxHeight = "26rem",
  stickyFirst = false,
  foot,
}: {
  label: string;
  rows: readonly Row[];
  columns: readonly LinkedColumn<Row>[];
  rowKey(row: Row): string;
  /** The row button's accessible name when the first cell is not enough. */
  rowLabel?(row: Row): string;
  selectedKey: string | null;
  /** Pins a row, or clears the pin when the pinned row is chosen again. */
  onSelect(row: Row | null): void;
  initialSort?: { id: string; dir: "asc" | "desc" };
  maxHeight?: string;
  stickyFirst?: boolean;
  /** A totals row under the body, by column id. */
  foot?: Record<string, ReactNode>;
}) => {
  const [sort, setSort] = useState(initialSort ?? null);
  const frame = useRef<HTMLDivElement | null>(null);
  const sorted = useMemo(() => {
    const col = sort === null ? undefined : columns.find((c) => c.id === sort.id);
    if (col?.sortValue === undefined || sort === null) return [...rows];
    const key = col.sortValue;
    return [...rows].sort((x, y) => compare(key(x), key(y), sort.dir === "asc" ? 1 : -1));
  }, [rows, columns, sort]);
  const fixed = columns.some((c) => c.width !== undefined);

  // Rows drawn so far: a first chunk now, the next in a later task; a pinned
  // row past them comes in at once.
  const step = Math.max(40, Math.floor(CHUNK_CELLS / Math.max(1, columns.length)));
  const [filled, setFilled] = useState(() => Math.min(rows.length, step));
  useEffect(() => {
    if (filled >= sorted.length) return undefined;
    const id = window.setTimeout(() => setFilled((f) => Math.min(sorted.length, f + step)), 0);
    return () => window.clearTimeout(id);
  }, [filled, sorted.length, step]);
  const pinnedAt = selectedKey === null ? -1 : sorted.findIndex((r) => rowKey(r) === selectedKey);
  const visible = sorted.slice(0, Math.max(filled, pinnedAt + 1));
  const busy = visible.length < sorted.length;

  // Bring the pinned row into the frame, wherever the pin came from; a row
  // already in view stays put, and the page itself never scrolls.
  useLayoutEffect(() => {
    const el = frame.current;
    if (el === null || selectedKey === null) return;
    const tr = [...el.querySelectorAll<HTMLTableRowElement>("tbody tr[data-key]")].find((r) => r.dataset.key === selectedKey);
    if (tr === undefined) return;
    const head = el.querySelector("thead")?.getBoundingClientRect().height ?? 0;
    const top = tr.getBoundingClientRect().top - el.getBoundingClientRect().top + el.scrollTop;
    const bottom = top + tr.getBoundingClientRect().height;
    if (top < el.scrollTop + head) el.scrollTop = Math.max(0, top - head);
    else if (bottom > el.scrollTop + el.clientHeight) el.scrollTop = bottom - el.clientHeight;
  }, [selectedKey, sorted]);

  // Stable handlers, so a row the next chunk leaves alone skips rendering.
  const walk = useCallback((e: KeyboardEvent<HTMLButtonElement>) => {
    if (e.key !== "ArrowDown" && e.key !== "ArrowUp") return;
    const buttons = [...(frame.current?.querySelectorAll<HTMLButtonElement>("tbody [data-row-open]") ?? [])];
    const next = buttons[buttons.indexOf(e.currentTarget) + (e.key === "ArrowDown" ? 1 : -1)];
    if (next !== undefined) {
      e.preventDefault();
      next.focus();
    }
  }, []);
  const toggle = useCallback((row: Row, on: boolean) => onSelect(on ? null : row), [onSelect]);

  const onSort = (col: LinkedColumn<Row>) =>
    setSort((s) => (s?.id === col.id ? { id: col.id, dir: s.dir === "asc" ? "desc" : "asc" } : { id: col.id, dir: col.firstDir ?? "desc" }));

  const header = (col: LinkedColumn<Row>, first: boolean) => {
    const here = sort?.id === col.id;
    const right = col.align === "right";
    const glyph = !here ? (
      <ArrowUpDown size={12} strokeWidth={1.75} aria-hidden="true" className="shrink-0 opacity-50" />
    ) : sort.dir === "asc" ? (
      <ArrowUp size={12} strokeWidth={2} aria-hidden="true" className="shrink-0" />
    ) : (
      <ArrowDown size={12} strokeWidth={2} aria-hidden="true" className="shrink-0" />
    );
    // The glyph sits before a right-aligned label and after a left-aligned
    // one, so the label's edge lines up with the values under it.
    const words = (
      <span className={`inline-flex flex-col ${right ? "items-end text-right" : "items-start text-left"}`}>
        <span className={`inline-flex items-center gap-1 whitespace-nowrap ${right ? "flex-row-reverse" : ""}`}>
          {col.header}
          {col.sortValue !== undefined ? glyph : null}
        </span>
        {col.unit !== undefined ? <span className="whitespace-nowrap font-normal text-muted-text">{col.unit}</span> : null}
      </span>
    );
    return (
      <th
        key={col.id}
        scope="col"
        aria-sort={here ? (sort.dir === "asc" ? "ascending" : "descending") : undefined}
        className={`${HEAD} ${right ? "text-right" : "text-left"} ${first && stickyFirst ? "left-0 z-[3]" : ""}`}
      >
        {col.sortValue !== undefined ? (
          <button
            type="button"
            onClick={() => onSort(col)}
            className={`ui-control-target inline-flex items-end rounded-sm transition-colors duration-fast ease-out hover:text-blue-deep active:scale-[0.98] ${FOCUS_RING}`}
          >
            {words}
          </button>
        ) : (
          words
        )}
      </th>
    );
  };

  const [first, ...rest] = columns;
  return (
    <div
      ref={frame}
      role="region"
      aria-label={label}
      aria-busy={busy}
      tabIndex={0}
      className={`scroll-thin relative overflow-auto overscroll-contain rounded-sm ${FOCUS_RING}`}
      style={{ maxHeight }}
    >
      <table className={`w-max min-w-full border-separate border-spacing-0 ${fixed ? "table-fixed" : ""}`}>
        {fixed ? (
          <colgroup>
            {columns.map((c) => (
              <col key={c.id} style={c.width === undefined ? undefined : { width: c.width }} />
            ))}
          </colgroup>
        ) : null}
        <thead>
          <tr>
            {first !== undefined ? header(first, true) : null}
            {rest.map((c) => header(c, false))}
          </tr>
        </thead>
        <tbody>
          {visible.map((row) => {
            const key = rowKey(row);
            return <BodyRow key={key} row={row} rowId={key} columns={columns} on={key === selectedKey} stickyFirst={stickyFirst} label={rowLabel?.(row)} onToggle={toggle} onWalk={walk} />;
          })}
        </tbody>
        {foot !== undefined ? (
          <tfoot>
            <tr>
              {columns.map((c, i) => (
                <td
                  key={c.id}
                  className={`border-t border-hair px-2 py-1.5 text-control font-semibold tabular-nums text-ink ${c.align === "right" ? "text-right" : "text-left"} ${i === 0 && stickyFirst ? "sticky left-0 z-[1] bg-paper" : ""}`}
                >
                  {foot[c.id] ?? null}
                </td>
              ))}
            </tr>
          </tfoot>
        ) : null}
      </table>
    </div>
  );
};
