import type { DatasetReport, Research, SelectiveResult } from "../../data/schema";
import type { Resource } from "../../data/load";
import { researchFor } from "../../data/research";
import { int, learnerInSentence, pct } from "../../data/format";
import { ScrollFrame } from "../../ui/primitives";

/** The selective-calling research from research.json: on which datasets the
 *  genome alone calls a confident subset with lineages held out. Only an
 *  entry whose guard matches this report shows, and only A1 and A2 carry
 *  numbers; every other category gives its reason. Nothing shows while
 *  research.json loads or when it is missing or invalid. */
export const SelectiveResearch = ({ r, research }: { r: DatasetReport; research: Resource<Research> }) => {
  if (research.status !== "ready") return null;
  const state = researchFor(research, r).selective;
  if (state.status === "loading" || state.status === "absent") return null;
  if (state.status === "mismatch") {
    return (
      <p className="type-note mt-4">
        Selective calling, where the genome alone would call a confident subset of isolates and the rest would go to the lab, has not run on this revision&apos;s held-out scores.
      </p>
    );
  }
  const { categories, limits } = research.data.selective;
  const entry = state.entry;
  const label = categories.find((c) => c.code === entry.category)?.label ?? entry.category;
  return (
    <div className="mt-8">
      <h3 className="flex flex-wrap items-center gap-x-2 gap-y-1 text-control font-semibold text-ink">
        Selective calling{" "}
        {entry.category === "no_lineage" ? null : <span className="font-normal text-muted-text">Category {entry.category}</span>}
      </h3>
      <p className="prose-measure mt-1 text-note text-ink">{entry.category === "A2" ? label.replace("{learner}", learnerInSentence(entry.result.learner)) : label}</p>
      {entry.category === "no_lineage" ? null : <p className="type-note prose-measure mt-2">{limits.bar}</p>}
      {entry.category === "A1" || entry.category === "A2" ? (
        <CalledSubset result={entry.result} notes={[limits.scope, limits.no_call, limits.optimism, ...(entry.category === "A1" ? [limits.procedure] : []), limits.per_100]} />
      ) : null}
    </div>
  );
};

/** A per-100 value, one decimal. */
const per100 = (v: number): string => v.toFixed(1);

/** One measure for resistant and susceptible isolates. A cell given as parts
 *  keeps each part on one line and wraps only between them, so on a narrow
 *  screen "VME 6.3%," and "2 of 32" stack instead of the table overflowing
 *  its scroll frame. */
const PairRow = ({ label, r, s }: { label: string; r: string | readonly string[]; s: string | readonly string[] }) => (
  <tr className="border-b border-hair">
    <th scope="row" className="px-2 py-1.5 text-left font-medium text-ink">
      {label}
    </th>
    {[r, s].map((cell, i) => (
      <td key={i} className="px-2 py-1.5 text-right tabular-nums">
        {(typeof cell === "string" ? [cell] : cell).flatMap((part, j) => [
          ...(j > 0 ? [" "] : []),
          <span key={j} className="whitespace-nowrap">
            {part}
          </span>,
        ])}
      </td>
    ))}
  </tr>
);

/** An A1 or A2 result: VME and ME count called isolates only, so they sit
 *  beside the share that got a call, and the per-100 rows count every
 *  isolate. The lab is hypothetical here, so its rows say would. */
const CalledSubset = ({ result: x, notes }: { result: SelectiveResult; notes: readonly string[] }) => (
  <>
    <p className="prose-measure mt-2 text-note text-ink">
      With lineages held out, the genome alone calls {pct(x.coverage)} of the {int(x.n)} isolates and gives no call for the other {pct(1 - x.coverage)}, which would go to the lab.
    </p>
    <div className="mt-3">
      <ScrollFrame label="Selective calling results">
        <table className="w-full max-w-xl border-collapse text-control">
          <thead>
            <tr className="border-b border-hair">
              <th scope="col" className="px-2 py-2 text-left text-tiny font-medium text-muted-text">
                <span className="sr-only">Measure</span>
              </th>
              {["Resistant", "Susceptible"].map((h) => (
                <th key={h} scope="col" className="px-2 py-2 text-right text-tiny font-medium text-muted-text">
                  {h}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            <PairRow label="Isolates" r={int(x.n_r)} s={int(x.n_s)} />
            <PairRow label="Given a genome call" r={pct(x.coverage_r)} s={pct(x.coverage_s)} />
            <PairRow label="Wrong among called" r={[`VME ${pct(x.vme)},`, `${int(x.fn)} of ${int(x.tp + x.fn)}`]} s={[`ME ${pct(x.me)},`, `${int(x.fp)} of ${int(x.fp + x.tn)}`]} />
            <PairRow label="95% confidence interval" r={[`${pct(x.vme_lo)} to`, pct(x.vme_hi)]} s={[`${pct(x.me_lo)} to`, pct(x.me_hi)]} />
          </tbody>
          <tbody>
            <tr>
              <th scope="rowgroup" colSpan={3} className="px-2 pb-1.5 pt-4 text-left text-tiny font-medium text-muted-text">
                Per 100 isolates
              </th>
            </tr>
            <PairRow label="Right call" r={[per100(x.per100r_caught), "caught"]} s={[per100(x.per100s_cleared), "cleared"]} />
            <PairRow label="Wrong call" r={[per100(x.per100r_missed), "missed"]} s={[per100(x.per100s_overcalled), "overcalled"]} />
            <PairRow label="No call, would go to the lab" r={per100(x.per100r_lab)} s={per100(x.per100s_lab)} />
          </tbody>
        </table>
      </ScrollFrame>
    </div>
    <ul className="type-note mt-3 list-disc space-y-0.5 pl-5">
      {notes.map((n) => (
        <li key={n}>{n}</li>
      ))}
    </ul>
  </>
);
