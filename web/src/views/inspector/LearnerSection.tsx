import { useRef, useState } from "react";
import type { DatasetReport, Learner, LearnerResult, Research, Scheme } from "../../data/schema";
import { int, LEARNER_NAME, LEARNER_ORDER, learnerInSentence, pct } from "../../data/format";
import { BARS } from "../../data/metrics";
import { bestLearner, countsNoCalls, learnerOf, NO_BEST_LEARNER, nullFolds, outscoringPartial, splitAbsence, splitOf, uncalledOf } from "../../data/derive/dataset";
import { linear, useWidth } from "../../charts/kit";
import { LearnerStripRow, StripAxis } from "../../charts/LearnerStrip";
import { HoverTip } from "../../ui/HoverTip";
import { BarPill, Dot, LegendItem, ScrollFrame, Section, Segmented } from "../../ui/primitives";

import type { Resource } from "../../data/load";
import { SelectiveResearch } from "./SelectiveResearch";

const STRIP_MIN = 220;

/** Balanced accuracy, or null for a row that called no isolate. */
const baOf = (l: LearnerResult | undefined): number | null => l?.metrics?.balanced_accuracy ?? null;

/** A balanced accuracy in words: "No call" for a row that called no isolate. */
const baText = (l: LearnerResult | undefined): string => (l?.metrics === null ? "No call" : pct(baOf(l)));

const foldsOf = (l: LearnerResult | undefined): number[] => (l?.fold_balanced_accuracy ?? []).filter((v): v is number => v !== null);

const NULL_FIT: Partial<Record<Learner, string>> = { gbt: "boosting kept no tree", lasso: "the lasso kept no weight" };

/** Share of isolates the row calls, with the null folds behind the rest. */
const CalledCell = ({ r, l }: { r: DatasetReport; l: LearnerResult }) => {
  const uncalled = uncalledOf(l);
  if (uncalled === null || uncalled === 0) return <span className="text-muted-text">All</span>;
  const total = r.counts.samples;
  const x = l.no_call!;
  return (
    <HoverTip
      focusable={false}
      lines={[
        {
          term: `${int(total - uncalled)} of ${int(total)} isolates called`,
          text: `In ${int(nullFolds(l))} of ${int(l.fold_complexity.length)} folds ${NULL_FIT[l.learner] ?? "the fit was null"}, so the ${int(x.resistant)} resistant and ${int(x.susceptible)} susceptible isolates held out there get no call. ${
            l.metrics === null ? "With no isolate called, there is no balanced accuracy, VME, ME, or AUROC." : "Balanced accuracy, VME, ME, and AUROC cover the isolates called."
          }`,
        },
      ]}
    >
      <span className="tabular-nums">{pct((total - uncalled) / total)}</span>
    </HoverTip>
  );
};

/** Null folds in a report written before no-calls, which scored them 0.5. */
const NullFoldCell = ({ l }: { l: LearnerResult }) => {
  const k = nullFolds(l);
  if (k === 0) return <span className="text-muted-text">None</span>;
  return (
    <HoverTip
      focusable={false}
      lines={[
        {
          term: `${int(k)} of ${int(l.fold_complexity.length)} folds`,
          text: `In these folds ${NULL_FIT[l.learner] ?? "the fit was null"}. This report scored their isolates 0.5, which the 0.5 cutoff calls resistant, so its rates count them as resistant calls.`,
        },
      ]}
    >
      <span className="tabular-nums">{`${int(k)} of ${int(l.fold_complexity.length)}`}</span>
    </HoverTip>
  );
};

export const LearnerSection = ({ r, research }: { r: DatasetReport; research?: Resource<Research> }) => {
  const [scheme, setScheme] = useState<Scheme>("random");
  const frame = useRef<HTMLDivElement | null>(null);
  // Wide, the strip takes what the number columns leave, within bounds; on
  // a phone it fits beside the learner names and the numbers scroll.
  const avail = useWidth(frame, 900);
  const stripWidth = avail < 640 ? Math.max(120, avail - 210) : Math.max(STRIP_MIN, Math.min(560, avail - 640));
  const random = splitOf(r, "random");
  const lineage = splitOf(r, "lineage");
  const best = bestLearner(r);
  const catalog = learnerOf(r, "random", "catalog");
  const learners = LEARNER_ORDER.filter((l) => random?.learners.some((x) => x.learner === l));
  const all = learners.flatMap((l) => {
    const a = learnerOf(r, "random", l);
    const b = learnerOf(r, "lineage", l);
    return [baOf(a) ?? 1, baOf(b) ?? 1, ...foldsOf(a), ...foldsOf(b)];
  });
  // Dots, not bars, so the axis may start above zero.
  const lo = Math.max(0, Math.floor((Math.min(...all) - 0.02) * 20) / 20);
  const x = linear([lo, 1], [8, stripWidth - 8]);
  const bestLineage = best === undefined ? undefined : learnerOf(r, "lineage", best.learner);

  // Reports with no-call counts get the share each row calls; older reports
  // get the null folds they scored 0.5, when they have any.
  const noCalls = countsNoCalls(r);
  const rows = learners.flatMap((l) => [learnerOf(r, "random", l), learnerOf(r, "lineage", l)]);
  const anyUncalled = rows.some((l) => l !== undefined && (uncalledOf(l) ?? 0) > 0);
  const anyWithoutCall = rows.some((l) => l?.metrics === null);
  const anyNullFold = !noCalls && rows.some((l) => l !== undefined && nullFolds(l) > 0);
  const extra = noCalls ? "Called" : anyNullFold ? "Null folds" : null;

  const catalogBa = baOf(catalog);
  const answer =
    best === undefined
      ? learners.some((l) => l !== "catalog")
        ? `${NO_BEST_LEARNER} on every split it ran, so none is marked best; Called gives the share each one calls.`
        : "No learner result is recorded for this dataset."
      : `${LEARNER_NAME[best.learner]} scores ${pct(best.metrics.balanced_accuracy)} on random splits${
          bestLineage === undefined ? "" : ` and ${pct(baOf(bestLineage))} with lineages held out`
        }${catalogBa === null ? "." : `; the catalog alone scores ${pct(catalogBa)}.`}`;

  // Higher random-split scores that lose to the best learner because they
  // leave isolates uncalled. The note sits outside the error toggle, so it
  // shows in the default view even when the no-calls are on the lineage split.
  const outscoring = best === undefined ? [] : outscoringPartial(r);
  const uncalledWords = (u: { scheme: Scheme; n: number }) => `${int(u.n)} of ${int(r.counts.samples)} isolates ${u.scheme === "random" ? "on random splits" : "with lineages held out"}`;
  const outscoringNote =
    best === undefined || outscoring.length === 0
      ? null
      : `${outscoring
          .map((o) => `${LEARNER_NAME[o.learner]} scores ${pct(o.ba)} on random splits, above ${learnerInSentence(best.learner)}, but gives no call for ${o.uncalled.map(uncalledWords).join(" and ")}.`)
          .join(" ")} The best learner has to call every isolate on every split it ran.`;

  const errorsOf = (l: Learner) => learnerOf(r, scheme, l);

  return (
    <Section
      id="ds-learners"
      title="Learner comparison"
      answer={answer}
      right={
        lineage !== undefined ? (
          <Segmented<Scheme>
            ariaLabel="Errors measured under"
            options={[
              { id: "random", label: "Errors, random split" },
              { id: "lineage", label: "Errors, lineages held out" },
            ]}
            value={scheme}
            onSelect={setScheme}
          />
        ) : undefined
      }
    >
      <div className="mb-3 flex flex-wrap gap-x-5 gap-y-1">
        <LegendItem swatch={<Dot color="var(--chart-mark)" />} name="Random split, pooled" />
        {lineage !== undefined ? <LegendItem swatch={<Dot color="var(--chart-mark)" ring />} name="Lineages held out, pooled" /> : null}
        <LegendItem swatch={<span className="inline-block h-2 w-2 rounded-full bg-[color:var(--chart-mark)] opacity-40" />} name="Single folds" />
        {catalog !== undefined ? <LegendItem swatch={<span className="inline-block h-0 w-4 border-t-2 border-dashed border-ink/55" />} name="Catalog" /> : null}
      </div>
      {outscoringNote === null ? null : <p className="type-note mb-3">{outscoringNote}</p>}
      <div ref={frame}>
      <ScrollFrame label="Learner comparison">
        <table className="w-full border-collapse text-control">
          <thead>
            <tr className="border-b border-hair text-left">
              <th scope="col" className="px-2 py-2 text-tiny font-medium text-blue">
                Learner
              </th>
              <th scope="col" className="px-2 py-1 text-tiny font-medium text-blue">
                <span className="sr-only">Balanced accuracy, </span>
                <StripAxis x={x} width={stripWidth} />
              </th>
              {["Random", "Lineages held out", "VME", "ME", ...(extra === null ? [] : [extra]), "AUROC"].map((h) => (
                <th key={h} scope="col" className="whitespace-nowrap px-2 py-2 text-right text-tiny font-medium text-blue">
                  {h}
                  {h === "VME" ? <span className="block font-normal text-muted-text">≤1.5%</span> : h === "ME" ? <span className="block font-normal text-muted-text">≤3%</span> : null}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {learners.map((l) => {
              const a = learnerOf(r, "random", l);
              const b = learnerOf(r, "lineage", l);
              const e = errorsOf(l);
              const m = e?.metrics;
              const isBest = best?.learner === l;
              return (
                <tr key={l} className="border-b border-hair">
                  <th scope="row" className="whitespace-nowrap px-2 py-1 text-left font-medium text-ink">
                    {LEARNER_NAME[l]}
                    {isBest ? <span className="ml-2 text-tiny font-medium text-blue">Best</span> : null}
                  </th>
                  <td className="px-2 py-0.5">
                    <div>
                      <LearnerStripRow
                        x={x}
                        width={stripWidth}
                        random={baOf(a)}
                        lineage={l === "catalog" ? null : baOf(b)}
                        foldsRandom={foldsOf(a)}
                        foldsLineage={foldsOf(b)}
                        catalog={catalogBa}
                        label={`${LEARNER_NAME[l]}: ${baText(a)} random${b === undefined || l === "catalog" ? "" : `, ${baText(b)} lineages held out`}`}
                      />
                    </div>
                  </td>
                  <td className={`px-2 py-1 text-right tabular-nums ${a?.metrics === null ? "text-muted-text" : ""}`}>{baText(a)}</td>
                  <td className={`px-2 py-1 text-right tabular-nums ${b?.metrics === null ? "text-muted-text" : ""}`}>
                    {l === "catalog" ? (
                      <span className="text-muted-text">
                        {baText(a)}
                        <span className="sr-only">, the same score: the catalog does not train</span>
                      </span>
                    ) : (
                      baText(b)
                    )}
                  </td>
                  {m === undefined ? (
                    <>
                      <td className="px-2 py-1 text-right">–</td>
                      <td className="px-2 py-1 text-right">–</td>
                    </>
                  ) : m === null ? (
                    <>
                      <td className="px-2 py-1 text-right text-muted-text">No call</td>
                      <td className="px-2 py-1 text-right text-muted-text">No call</td>
                    </>
                  ) : (
                    <>
                      <td className="px-2 py-1 text-right">
                        <BarPill value={m.vme} bar={BARS.vme} />
                      </td>
                      <td className="px-2 py-1 text-right">
                        <BarPill value={m.me} bar={BARS.me} />
                      </td>
                    </>
                  )}
                  {extra === null ? null : (
                    <td className="px-2 py-1 text-right">{e === undefined ? "–" : noCalls ? <CalledCell r={r} l={e} /> : <NullFoldCell l={e} />}</td>
                  )}
                  <td className="px-2 py-1 text-right tabular-nums">{m?.auroc === null || m?.auroc === undefined ? <span aria-label="AUROC unavailable">–</span> : m.auroc.toFixed(3)}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </ScrollFrame>
      </div>
      <p className="type-note mt-3">
        Balanced accuracy over pooled out-of-fold predictions; faint dots are the {random?.folds ?? 5} single folds.
        {lineage !== undefined
          ? ` With lineages held out, each of the ${lineage.n_groups} lineages goes whole to one fold.`
          : ` No lineage-held-out result: ${splitAbsence(r, "lineage") ?? ""}`}
        {anyUncalled
          ? " A fold where boosting kept no tree, or the lasso kept no weight, gets no call: every rate covers the isolates called, and Called gives their share."
          : ""}
        {anyWithoutCall ? " A learner that calls no isolate on a split has no rate there and shows No call." : ""}
        {anyNullFold
          ? " This report scored a fold where boosting kept no tree, or the lasso kept no weight, at 0.5, which the 0.5 cutoff calls resistant; Null folds counts those folds, and every rate includes them."
          : ""}
      </p>
      {r.splits.some((split) => split.learners.some((l) => l.conformal.length > 0)) ? (
        <p className="type-note mt-4">This report's cannot-call rates came from pooled out-of-fold scores with no separate calibration set, so they are not shown as conformal results.</p>
      ) : null}
      {research === undefined ? null : <SelectiveResearch r={r} research={research} />}
    </Section>
  );
};
