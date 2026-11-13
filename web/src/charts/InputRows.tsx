// The input comparison, one row per changed input: the datasets where the
// change scored better than the default and worse, over the four species,
// with the sign test's p and the median change. A bar in each row draws the
// two counts out from a centre line, worse to the left and better to the
// right, on one scale for every row; the numbers sit in their own cells, so
// the table reads the same without the bars. The default's repeat runs are a
// group of their own under the changed inputs. On a slide too short for every
// row, the rows page (the deck's paging, as for the panel figure); on a narrow
// slide the row keeps its name, its two counts and the sign test.

import { useEffect, useLayoutEffect, useRef, useState } from "react";
import { int } from "../data/format";
import { ScrollFrame } from "../ui/primitives";
import { useWidth } from "./kit";

export interface InputRow {
  key: string;
  name: string;
  better: number;
  worse: number;
  datasets: number;
  p: number;
  median: number;
}

export interface InputGroup {
  title: string;
  rows: readonly InputRow[];
}

const TH = "px-2 py-2 text-tiny font-medium text-blue";
/** A slide's row, group line and head before the table is laid out; once it
 *  is, the measured heights replace them, since a long name wraps on a
 *  narrow slide. And the continuation line's height. */
const ROW_H = 30;
const HEAD_H = 40;
const GROUP_H = 34;
const PAGER_H = 24;
/** The width under which a slide keeps the counts and the sign test only:
 *  no bars, and the datasets compared and the median stay on the page. */
const BARS_FROM = 560;
const BAR_W = 140;

/** A p value or a median as the document prints them, with a true minus sign. */
export const fixed4 = (v: number): string => (v < 0 ? `−${(-v).toFixed(4)}` : v.toFixed(4));

const Bar = ({ row, scale }: { row: InputRow; scale: number }) => {
  const half = BAR_W / 2;
  const w = (n: number) => (scale === 0 ? 0 : (n / scale) * (half - 2));
  return (
    <svg aria-hidden="true" width={BAR_W} height={14} className="block">
      <rect x={half - w(row.worse)} y={2} width={w(row.worse)} height={10} fill="var(--series-neutral)" rx={1} />
      <rect x={half} y={2} width={w(row.better)} height={10} fill="var(--chart-mark)" rx={1} />
      <line x1={half} x2={half} y1={0} y2={14} stroke="var(--chart-ref)" />
    </svg>
  );
};

export const InputRows = ({
  label,
  groups,
  slide = false,
  height,
  page,
  pages,
  onPages,
}: {
  label: string;
  groups: readonly InputGroup[];
  slide?: boolean;
  height?: number;
  page?: number;
  pages?: number;
  onPages?(pages: number): void;
}) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const legendRef = useRef<HTMLDivElement | null>(null);
  const [legendH, setLegendH] = useState(24);
  useLayoutEffect(() => {
    const el = legendRef.current;
    if (el === null) return undefined;
    const read = () => setLegendH(el.offsetHeight);
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  const tableRef = useRef<HTMLTableElement | null>(null);
  const [sizes, setSizes] = useState({ row: ROW_H, group: GROUP_H, head: HEAD_H });
  // Every row in order, each with its group, so a page can start a group anywhere.
  const all = groups.flatMap((g) => g.rows.map((r) => ({ group: g.title, row: r })));
  const room = (continued: boolean) => (height ?? Infinity) - sizes.head - legendH - groups.length * sizes.group - (continued ? PAGER_H : 0);
  const needed = height === undefined || room(false) >= all.length * sizes.row ? 1 : Math.ceil(all.length / Math.max(1, Math.floor(room(true) / sizes.row)));
  // An unlaid box measures nothing, so it reports nothing.
  useEffect(() => {
    if (height !== undefined && height > 0) onPages?.(needed);
  }, [onPages, needed, height]);
  const count = Math.max(1, Math.min(pages ?? 1, all.length));
  const per = Math.ceil(all.length / count);
  const at = Math.min(Math.max(0, page ?? 0), count - 1);
  const shown = count === 1 ? all : all.slice(at * per, at * per + per);
  const shownGroups = groups.map((g) => ({ title: g.title, rows: shown.filter((s) => s.group === g.title).map((s) => s.row) })).filter((g) => g.rows.length > 0);
  const scale = Math.max(0, ...all.map((s) => Math.max(s.row.better, s.row.worse)));
  const full = !slide || width >= BARS_FROM;
  // The tallest row, group line and head of the whole table at this width;
  // only a copy that shows every row measures, so the pages it reports
  // never depend on the page it shows.
  useLayoutEffect(() => {
    const t = tableRef.current;
    if (t === null || count !== 1) return;
    const tallest = (sel: string, fallback: number) => {
      const h = Math.max(0, ...[...t.querySelectorAll<HTMLElement>(sel)].map((el) => el.offsetHeight));
      return h > 0 ? h : fallback;
    };
    const next = { row: tallest("tbody tr[data-input]", ROW_H), group: tallest("tbody tr:not([data-input])", GROUP_H), head: tallest("thead", HEAD_H) };
    setSizes((was) => (was.row === next.row && was.group === next.group && was.head === next.head ? was : next));
  }, [width, count, full, legendH]);
  const head = full ? ["Change", "Worse", "", "Better", "Of", "Sign p", "Median"] : ["Change", "Worse", "Better", "Sign p"];
  const table = (
    <table ref={tableRef} aria-label={slide ? label : undefined} className={`w-full border-collapse ${slide ? "text-control sm:text-body" : "min-w-[38rem] text-control"}`}>
      <thead>
        <tr className="border-b border-hair">
          {head.map((h, i) => (
            <th key={h === "" ? "bar" : h} scope="col" className={`${TH} ${i === 0 ? "text-left" : h === "" ? "text-center" : "text-right"}`}>
              {h === "" ? <span className="sr-only">Worse and better, drawn</span> : h}
            </th>
          ))}
        </tr>
      </thead>
      {shownGroups.map((g) => (
        <tbody key={g.title} data-input-group={g.title}>
          <tr className="border-b border-hair">
            <th scope="rowgroup" colSpan={head.length} className="px-2 pb-1.5 pt-4 text-left text-tiny font-medium text-muted-text">
              {g.title}
            </th>
          </tr>
          {g.rows.map((r) => (
            <tr key={r.key} data-input={r.key} className="border-b border-hair">
              <th scope="row" className="px-2 py-1 text-left font-normal">
                {r.name}
              </th>
              <td className="px-2 py-1 text-right tabular-nums">{int(r.worse)}</td>
              {full ? (
                <td className="px-2 py-1">
                  <Bar row={r} scale={scale} />
                </td>
              ) : null}
              <td className="px-2 py-1 text-right tabular-nums">{int(r.better)}</td>
              {full ? <td className="px-2 py-1 text-right tabular-nums">{int(r.datasets)}</td> : null}
              <td className="px-2 py-1 text-right tabular-nums">{fixed4(r.p)}</td>
              {full ? <td className="px-2 py-1 text-right tabular-nums">{fixed4(r.median)}</td> : null}
            </tr>
          ))}
        </tbody>
      ))}
    </table>
  );
  const first = count === 1 ? 1 : at * per + 1;
  const last = count === 1 ? all.length : at * per + shown.length;
  return (
    <div ref={frame} className={slide ? "min-w-0" : "min-w-0 max-w-[52rem]"}>
      {count > 1 ? (
        <p data-continued className="text-control text-muted-text" style={{ height: PAGER_H }}>
          {`Rows ${int(first)}–${int(last)} of ${int(all.length)}`}
        </p>
      ) : null}
      {slide ? table : <ScrollFrame label={label}>{table}</ScrollFrame>}
      <div ref={legendRef} className="flex flex-wrap gap-x-5 gap-y-1 pt-2 text-control text-ink">
        {full ? (
          <>
            <span className="inline-flex items-center gap-1.5">
              <svg aria-hidden="true" width={12} height={12}>
                <rect x={1} y={1} width={10} height={10} fill="var(--series-neutral)" rx={1} />
              </svg>
              Datasets where the change scored worse, to the left
            </span>
            <span className="inline-flex items-center gap-1.5">
              <svg aria-hidden="true" width={12} height={12}>
                <rect x={1} y={1} width={10} height={10} fill="var(--chart-mark)" rx={1} />
              </svg>
              Better, to the right
            </span>
          </>
        ) : null}
        {full ? (
          <span className="text-muted-text">
            Of: datasets compared; on any not counted worse or better, the score was exactly the same. Median: the change in log loss on the middle dataset, above zero when the change
            scored better.
          </span>
        ) : null}
      </div>
    </div>
  );
};
