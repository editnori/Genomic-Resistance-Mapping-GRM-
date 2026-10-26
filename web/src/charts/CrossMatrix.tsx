// Cross-resistance: phi between resistance calls for each drug pair, over the
// genomes called R or S for both drugs. The matrix is the pattern overview:
// the lower triangle, drugs grouped by class with a gap between classes,
// full drug names, cells never under 24px (the frame scrolls instead), and a
// stepped diverging legend that counts the pairs in each step: coral when
// resistance to one goes with resistance to the other, blue when it goes
// against. A pair with fewer than MIN_JOINT joint calls is not evaluated
// (an empty cell); an evaluated pair with a zero margin has no phi (a
// hatched cell). Every cell can be pinned by click, tap or Enter: the strip
// under the matrix and the pair's table row then show its exact counts, and
// choosing a table row pins its cell the same way. A display of the existing
// phi values; no new analysis.

import { useEffect, useId, useMemo, useState, type KeyboardEvent, type ReactNode } from "react";
import type { SpeciesReport } from "../data/schema";
import { capitalize, int } from "../data/format";
import { GhostButton, ScrollFrame, Segmented } from "../ui/primitives";
import { ChartFrame, ChartTip, LABEL_CLASS, keyboardFocus, labelWidth, useSvgId } from "./kit";
import { MIN_JOINT, undefinedWhy, type PairRow, type PairSummary } from "./crossPairs";
import { LinkedTable, type LinkedColumn } from "./LinkedTable";

const STEPS = [
  { key: "strong-together", label: "0.5 or more", fill: "color-mix(in srgb, var(--series-r) 85%, transparent)", test: (phi: number) => phi >= 0.5 },
  { key: "together", label: "0.2 to 0.5", fill: "color-mix(in srgb, var(--series-r) 38%, transparent)", test: (phi: number) => phi >= 0.2 && phi < 0.5 },
  { key: "none", label: "−0.2 to 0.2", fill: "var(--fill-3)", test: (phi: number) => phi > -0.2 && phi < 0.2 },
  { key: "against", label: "−0.5 to −0.2", fill: "color-mix(in srgb, var(--series-s) 38%, transparent)", test: (phi: number) => phi > -0.5 && phi <= -0.2 },
  { key: "strong-against", label: "−0.5 or less", fill: "color-mix(in srgb, var(--series-s) 85%, transparent)", test: (phi: number) => phi <= -0.5 },
] as const;

const stepOf = (phi: number) => STEPS.find((s) => s.test(phi)) ?? STEPS[2];

const CELL = 24;
const GROUP_GAP = 6;

const pairName = (p: PairRow): string => `${capitalize(p.a)} and ${p.b}`;

/** The pair in words, with its exact joint counts. */
export const pairWords = (p: PairRow): string => {
  const cells = `${int(p.rr)} resistant to both, ${int(p.rs)} to ${p.a} only, ${int(p.sr)} to ${p.b} only, ${int(p.ss)} susceptible to both`;
  if (p.status === "evaluated") return `Phi ${p.phi!.toFixed(2)} over ${int(p.n)} genomes called R or S for both: ${cells}.`;
  if (p.status === "undefined") return `Phi undefined: of the ${int(p.n)} genomes called R or S for both, ${undefinedWhy(p)}. ${capitalize(cells)}.`;
  return `Not evaluated: ${int(p.n)} ${p.n === 1 ? "genome is" : "genomes are"} called R or S for both, fewer than ${MIN_JOINT}.${p.n > 0 ? ` ${capitalize(cells)}.` : ""}`;
};

const CrossMatrix = ({
  report,
  order,
  pairs,
  selected,
  onSelect,
}: {
  report: SpeciesReport;
  order: readonly number[];
  pairs: readonly PairRow[];
  selected: string | null;
  onSelect(key: string | null): void;
}) => {
  const cellId = useId();
  const hatch = useSvgId("phi-undefined");
  const byCols = useMemo(() => new Map<string, PairRow>(pairs.flatMap((p) => [[`${p.ia}|${p.ib}`, p] as const, [`${p.ib}|${p.ia}`, p] as const])), [pairs]);
  // Classes in the order of their first drug in `order` (most R or S calls
  // first); drugs keep that order inside a class.
  const groups = new Map<string, number[]>();
  order.forEach((j) => {
    const c = report.drug_classes[j] ?? "";
    groups.set(c, [...(groups.get(c) ?? []), j]);
  });
  const seq = [...groups.values()].flat();
  const classes = [...groups.keys()];
  const groupIndex = new Map(seq.map((j) => [j, classes.indexOf(report.drug_classes[j] ?? "")]));
  const pos = (i: number) => i * CELL + (groupIndex.get(seq[i]!) ?? 0) * GROUP_GAP;
  const names = seq.map((j) => capitalize(report.drugs[j] ?? ""));
  const n = names.length;
  const label = Math.max(96, ...names.map(labelWidth)) + 10;
  const span = n === 0 ? 0 : pos(n - 1) + CELL;
  const top = 8;
  const bottom = top + span;
  const width = label + span + 24;
  const height = bottom + Math.ceil(Math.max(0, ...names.map(labelWidth)) * 0.82) + 20;
  const pairAt = (i: number, k: number) => byCols.get(`${seq[i]}|${seq[k]}`);

  const sel = pairs.find((p) => p.key === selected);
  const selI = sel === undefined ? -1 : Math.max(seq.indexOf(sel.ia), seq.indexOf(sel.ib));
  const selK = sel === undefined ? -1 : Math.min(seq.indexOf(sel.ia), seq.indexOf(sel.ib));
  const [cursor, setCursor] = useState({ i: 1, k: 0 });
  // A pin from the table moves the keyboard's cell to it.
  useEffect(() => {
    if (selI > 0) setCursor({ i: selI, k: selK });
  }, [selI, selK]);
  const ci = Math.max(1, Math.min(cursor.i, n - 1));
  const cur = { i: ci, k: Math.max(0, Math.min(cursor.k, ci - 1)) };
  const [hover, setHover] = useState<{ i: number; k: number } | null>(null);
  const [focused, setFocused] = useState(false);
  const [tipOff, setTipOff] = useState(false);

  const toggle = (p: PairRow | undefined) => {
    if (p !== undefined) onSelect(p.key === selected ? null : p.key);
  };

  const onKey = (e: KeyboardEvent<HTMLDivElement>) => {
    if (n < 2) return;
    if (e.key === "Escape") {
      if (selected !== null) {
        e.preventDefault();
        onSelect(null);
      } else if (!tipOff && (focused || hover !== null)) {
        e.preventDefault();
        setHover(null);
        setTipOff(true);
      }
      return;
    }
    if (e.key === "Enter" || e.key === " ") {
      e.preventDefault();
      toggle(pairAt(cur.i, cur.k));
      return;
    }
    const last = n - 1;
    const next =
      e.key === "ArrowRight"
        ? { i: cur.i, k: Math.min(cur.k + 1, cur.i - 1) }
        : e.key === "ArrowLeft"
          ? { i: cur.i, k: Math.max(cur.k - 1, 0) }
          : e.key === "ArrowDown"
            ? { i: Math.min(cur.i + 1, last), k: cur.k }
            : e.key === "ArrowUp"
              ? { i: Math.max(cur.i - 1, 1), k: Math.min(cur.k, Math.max(cur.i - 2, 0)) }
              : e.key === "Home"
                ? { i: cur.i, k: 0 }
                : e.key === "End"
                  ? { i: cur.i, k: cur.i - 1 }
                  : e.key === "PageDown"
                    ? { i: last, k: cur.k }
                    : e.key === "PageUp"
                      ? { i: 1, k: 0 }
                      : null;
    if (next === null) return;
    e.preventDefault();
    setHover(null);
    setTipOff(false);
    setCursor(next);
  };

  // A pair not evaluated is an empty cell with a dashed edge.
  const fill = (p: PairRow | undefined) =>
    p === undefined || p.status === "not-evaluated" ? "transparent" : p.status === "undefined" ? `url(#${hatch})` : stepOf(p.phi!).fill;
  const at = cur;
  const atPair = pairAt(at.i, at.k);
  const shown = hover ?? (focused && !tipOff ? at : null);
  const shownPair = shown === null ? undefined : pairAt(shown.i, shown.k);
  const counts = STEPS.map((s) => pairs.filter((p) => p.phi !== null && s.test(p.phi)).length);
  const undefinedCount = pairs.filter((p) => p.status === "undefined").length;
  const notEvaluated = pairs.filter((p) => p.status === "not-evaluated").length;
  const bold = (i: number) => (sel !== undefined && (i === selI || i === selK) ? "font-semibold" : "");

  return (
    <div>
      <ScrollFrame label="Cross-resistance matrix">
        <ChartFrame>
          <div
            role="grid"
            tabIndex={0}
            aria-label={`Phi between resistance calls for ${int(pairs.length)} drug pairs, lower triangle, drugs grouped by class. Arrow keys move; click or Enter to pin details; Escape clears.`}
            aria-rowcount={n}
            aria-colcount={n}
            aria-activedescendant={atPair === undefined ? undefined : `${cellId}-${at.i}-${at.k}`}
            onKeyDown={onKey}
            onFocus={(e) => {
              setFocused(keyboardFocus(e.currentTarget));
              setTipOff(false);
            }}
            onBlur={() => setFocused(false)}
            onPointerLeave={() => setHover(null)}
            className="w-max rounded-sm focus-visible:outline focus-visible:outline-2 focus-visible:-outline-offset-2 focus-visible:outline-blue"
          >
            <svg width={width} height={height} aria-hidden="true" className="block [forced-color-adjust:none]">
              <defs>
                <pattern id={hatch} patternUnits="userSpaceOnUse" width={5} height={5} patternTransform="rotate(45)">
                  <rect width={5} height={5} fill="var(--fill-1)" />
                  <line x1={0} y1={0} x2={0} y2={5} stroke="var(--chart-axis)" strokeWidth={2} />
                </pattern>
              </defs>
              {names.map((a, i) => (
                <g key={a}>
                  {i > 0 ? (
                    <text x={label - 6} y={top + pos(i) + CELL / 2 + 4} textAnchor="end" className={`${LABEL_CLASS} ${bold(i)}`}>
                      {a}
                    </text>
                  ) : null}
                  {i < n - 1 ? (
                    <text transform={`translate(${label + pos(i) + CELL / 2 - 3} ${bottom + 8}) rotate(55)`} className={`${LABEL_CLASS} ${bold(i)}`}>
                      {a}
                    </text>
                  ) : null}
                  {names.slice(0, i).map((b, k) => {
                    const p = pairAt(i, k);
                    const x = label + pos(k);
                    const y = top + pos(i);
                    return (
                      <rect
                        key={b}
                        x={x + 1}
                        y={y + 1}
                        width={CELL - 2}
                        height={CELL - 2}
                        rx={2}
                        fill={fill(p)}
                        stroke={p === undefined || p.status === "not-evaluated" ? "var(--chart-axis)" : "none"}
                        strokeDasharray={p === undefined || p.status === "not-evaluated" ? "2 2" : undefined}
                        className="cursor-pointer"
                        onPointerEnter={(e) => {
                          if (e.pointerType !== "touch") setHover({ i, k });
                        }}
                        onClick={() => {
                          setCursor({ i, k });
                          toggle(p);
                        }}
                      />
                    );
                  })}
                </g>
              ))}
              {sel !== undefined && selI > 0 ? (
                <rect x={label + pos(selK)} y={top + pos(selI)} width={CELL} height={CELL} rx={3} fill="none" stroke="rgb(var(--ink))" strokeWidth={2} pointerEvents="none" />
              ) : null}
              {focused && n > 1 ? (
                <rect x={label + pos(at.k) - 2} y={top + pos(at.i) - 2} width={CELL + 4} height={CELL + 4} rx={4} fill="none" stroke="rgb(var(--blue))" strokeWidth={2} pointerEvents="none" />
              ) : null}
            </svg>
            {atPair !== undefined ? (
              <div role="row" aria-rowindex={at.i + 1} className="sr-only">
                <div role="gridcell" id={`${cellId}-${at.i}-${at.k}`} aria-colindex={at.k + 1} aria-selected={atPair.key === selected}>
                  {`${pairName(atPair)}. ${pairWords(atPair)}`}
                </div>
              </div>
            ) : null}
          </div>
          {shown !== null && shownPair !== undefined ? (
            <ChartTip
              x={label + pos(shown.k) + CELL}
              y={top + pos(shown.i)}
              width={width}
              lines={[{ term: pairName(shownPair), text: shownPair.status === "evaluated" ? `phi ${shownPair.phi!.toFixed(2)}, ${int(shownPair.n)} joint calls` : shownPair.status === "undefined" ? `phi undefined, ${int(shownPair.n)} joint calls` : `not evaluated, ${int(shownPair.n)} joint calls` }]}
            />
          ) : null}
        </ChartFrame>
      </ScrollFrame>
      <ul className="mt-3 flex flex-wrap gap-x-5 gap-y-1 text-control text-ink" aria-label="Phi steps">
        {STEPS.map((s, i) => (
          <li key={s.key} className="inline-flex items-center gap-1.5">
            <span aria-hidden="true" className="inline-block h-3 w-4 rounded-sm [forced-color-adjust:none]" style={{ background: s.fill }} />
            <span className="text-muted-text">Phi {s.label}</span>
            <span className="font-medium tabular-nums">{int(counts[i] ?? 0)}</span>
          </li>
        ))}
        <li className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={16} height={12} className="shrink-0 rounded-sm [forced-color-adjust:none]">
            <rect width={16} height={12} fill={`url(#${hatch})`} />
          </svg>
          <span className="text-muted-text">Phi undefined (no R or no S for one drug)</span>
          <span className="font-medium tabular-nums">{int(undefinedCount)}</span>
        </li>
        <li className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={16} height={12} className="shrink-0 [forced-color-adjust:none]">
            <rect x={0.5} y={0.5} width={15} height={11} rx={2} fill="transparent" stroke="var(--chart-axis)" strokeDasharray="2 2" />
          </svg>
          <span className="text-muted-text">Not evaluated (fewer than {MIN_JOINT} joint R/S calls)</span>
          <span className="font-medium tabular-nums">{int(notEvaluated)}</span>
        </li>
      </ul>
    </div>
  );
};

type PairFilter = "all" | "evaluated" | "undefined" | "not-evaluated";

/** The matrix, the pinned pair's strip, and every pair as a table, all on
 *  one selection. */
export const CrossPanel = ({
  report,
  order,
  pairs,
  summary,
  selected,
  onSelect,
  caption,
}: {
  report: SpeciesReport;
  order: readonly number[];
  pairs: readonly PairRow[];
  summary: PairSummary;
  selected: string | null;
  onSelect(key: string | null): void;
  /** What phi is and when a pair counts; kept beside the matrix. */
  caption?: ReactNode;
}) => {
  const [filter, setFilter] = useState<PairFilter>("all");
  const sel = pairs.find((p) => p.key === selected);
  // A pin from the matrix shows its row even when the filter hid it.
  const pinFromMatrix = (key: string | null) => {
    const p = pairs.find((q) => q.key === key);
    if (p !== undefined && filter !== "all" && (filter === "evaluated" ? p.status !== "evaluated" : p.status !== filter)) setFilter("all");
    onSelect(key);
  };
  const rows = filter === "all" ? pairs : pairs.filter((p) => p.status === filter);
  // The pair column is set from the species' longest pair name, capped so a
  // long name wraps to a second line on a phone, and sticks while the counts
  // scroll; no filter or sort moves a column edge.
  const pairWidth = `${Math.min(208, Math.max(160, ...pairs.map((p) => labelWidth(pairName(p)))) + 28)}px`;
  const columns: LinkedColumn<PairRow>[] = [
    { id: "pair", header: "Drug pair", width: pairWidth, sortValue: (p) => `${p.a} ${p.b}`, firstDir: "asc", cell: pairName },
    { id: "n", header: "Joint calls", unit: "R or S for both", align: "right", width: "8rem", sortValue: (p) => p.n, cell: (p) => int(p.n) },
    { id: "rr", header: "Both R", align: "right", width: "5rem", sortValue: (p) => p.rr, cell: (p) => int(p.rr) },
    { id: "rs", header: "R, S", unit: "1st, 2nd", align: "right", width: "5rem", sortValue: (p) => p.rs, cell: (p) => int(p.rs) },
    { id: "sr", header: "S, R", unit: "1st, 2nd", align: "right", width: "5rem", sortValue: (p) => p.sr, cell: (p) => int(p.sr) },
    { id: "ss", header: "Both S", align: "right", width: "5rem", sortValue: (p) => p.ss, cell: (p) => int(p.ss) },
    {
      id: "phi",
      header: "Phi",
      align: "right",
      width: "7.5rem",
      sortValue: (p) => p.phi,
      cell: (p) => (p.status === "evaluated" ? p.phi!.toFixed(2) : <span className="text-muted-text">{p.status === "undefined" ? "undefined" : "not evaluated"}</span>),
    },
  ];
  if (summary.total === 0) return <p className="type-note">This species file lists one drug, so it has no pairs.</p>;
  return (
    <div>
      <CrossMatrix report={report} order={order} pairs={pairs} selected={selected} onSelect={pinFromMatrix} />
      {caption !== undefined ? <p className="prose-measure mt-3 text-note text-muted-text">{caption}</p> : null}
      <div aria-live="polite" className="mt-4 flex min-h-[3.5rem] flex-wrap items-start justify-between gap-x-6 gap-y-2 rounded-lg border border-hair px-4 py-3">
        {sel === undefined ? (
          <p className="text-control text-muted-text">Click, tap or press Enter on a cell, or choose a table row, to pin a pair and read its joint counts here.</p>
        ) : (
          <>
            <p className="min-w-0 flex-1 basis-80 text-control text-ink">
              <span className="font-semibold">{pairName(sel)}</span>
              {". "}
              {pairWords(sel)}
            </p>
            <GhostButton onClick={() => onSelect(null)} ariaLabel="Clear the pinned pair">
              Clear
            </GhostButton>
          </>
        )}
      </div>
      <div className="mt-6 flex flex-wrap items-center justify-between gap-x-6 gap-y-2">
        <h3 className="text-control font-semibold text-ink">{`Every pair (${int(summary.total)})`}</h3>
        <Segmented<PairFilter>
          ariaLabel="Show pairs"
          options={[
            { id: "all", label: `All ${int(summary.total)}` },
            { id: "evaluated", label: `Phi defined ${int(summary.evaluated - summary.undefined)}` },
            { id: "undefined", label: `Phi undefined ${int(summary.undefined)}`, disabled: summary.undefined === 0 },
            { id: "not-evaluated", label: `Not evaluated ${int(summary.notEvaluated)}`, disabled: summary.notEvaluated === 0 },
          ]}
          value={filter}
          onSelect={setFilter}
        />
      </div>
      <div className="mt-3">
        <LinkedTable
          label="Cross-resistance pairs"
          rows={rows}
          columns={columns}
          rowKey={(p) => p.key}
          selectedKey={selected}
          onSelect={(p) => onSelect(p === null ? null : p.key)}
          initialSort={{ id: "phi", dir: "desc" }}
          maxHeight="24rem"
          stickyFirst
        />
      </div>
    </div>
  );
};
