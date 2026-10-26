// The five findings Overview and Present share: corpus, lineage cost, the
// known-gene catalog, the lead PCR panel, and MIC prediction. Each is a title,
// one result sentence with its count and denominator (the corpus counts are
// their own result), a figure that shows that number, and one caption that
// says what a mark is and how it was made; Overview adds Method and coverage.
// Every number comes from data/derive at render time.

import { useRef, type ReactNode } from "react";
import type { CatalogRow, Manifest } from "../data/schema";
import { capitalize, int, listWords, pct, plural, speciesShort } from "../data/format";
import { CATALOG_CUTOFF, catalogSummary, classMinimum, classRows, datasetRule, micSummary, noBestWords, speciesCount, splitSummary, type ClassRow } from "../data/derive/overview";
import { leadPanel, panelsWithoutBest, pool, share, type PanelRowJoined } from "../data/derive/tests";
import { BEST_LEARNER_WORDS, splitReasonShort } from "../data/derive/dataset";
import { AxisX, LABEL_CLASS, linear, MarkListbox, PlacedText, placeLabels, RefLine, ticks, useWidth, type ListMark } from "../charts/kit";
import { GapSwarm, MicStrip } from "../charts/DotStrip";
import { CatalogScatter } from "../charts/CatalogScatter";
import type { FigureContent, Fit } from "../ui/primitives";
import { PanelEvidence } from "./PanelEvidence";

export type FindingId = "corpus" | "lineage" | "catalog" | "panels" | "mic";

export interface Finding {
  id: FindingId;
  title: string;
  /** Only where the title alone leaves the reader guessing; page only. */
  question?: string;
  /** Left out when the figure is the result, as the corpus counts are. */
  result?: string;
  figure: FigureContent;
  caption: string;
  /** Overview only. */
  details?: ReactNode;
}

export interface FindingOptions {
  /** Opens a dataset from a figure mark; Present leaves marks unopenable. */
  onOpen?(id: string): void;
  /** Overview builds Method and coverage; Present does not. */
  page: boolean;
}

// ---- corpus -------------------------------------------------------------------

/** The four counts the corpus rests on, each with its unit, side by side
 *  rather than as one funnel because the units differ. */
export const CorpusFigure = ({ manifest, large = false }: { manifest: Manifest; large?: boolean }) => {
  const f = manifest.funnel;
  const stages = [
    { value: manifest.amr_records, unit: "BV-BRC lab records" },
    { value: f.pairs_reported, unit: "Genome–antibiotic pairs with a reported R or S call" },
    { value: f.genomes_passing_qc, unit: "Genomes that passed QC" },
    { value: manifest.datasets.length, unit: `Datasets across ${plural(speciesCount(manifest), "species", "species")}` },
  ];
  return (
    <ol className={`grid grid-cols-2 gap-x-6 md:grid-cols-4 ${large ? "gap-y-8" : "gap-y-5"}`}>
      {stages.map((s) => (
        <li key={s.unit} className="min-w-0">
          <div className={`font-semibold tabular-nums text-ink ${large ? "text-stat sm:text-[2.75rem] sm:leading-[3rem] sm:tracking-tight" : "text-stat"}`}>{int(s.value)}</div>
          <div className={`mt-1 text-muted-text ${large ? "text-note sm:text-body" : "text-note"}`}>{s.unit}</div>
        </li>
      ))}
    </ol>
  );
};

const SpeciesCounts = ({ manifest }: { manifest: Manifest }) => {
  const rows = manifest.species.filter((s) => s.datasets_curated > 0 || s.genomes_analyzed > 0);
  return (
    <table className="w-full max-w-xl border-collapse text-control">
      <thead>
        <tr className="border-b border-hair text-tiny text-blue">
          <th scope="col" className="px-2 py-1.5 text-left font-medium">
            Species
          </th>
          <th scope="col" className="px-2 py-1.5 text-right font-medium">
            Datasets
          </th>
          <th scope="col" className="px-2 py-1.5 text-right font-medium">
            Genomes
          </th>
        </tr>
      </thead>
      <tbody>
        {rows.map((s) => (
          <tr key={s.slug} className="border-b border-hair">
            <th scope="row" className="px-2 py-1.5 text-left font-normal italic text-ink">
              {s.species}
            </th>
            <td className="px-2 py-1.5 text-right tabular-nums">{s.datasets_reported === s.datasets_curated ? int(s.datasets_reported) : `${int(s.datasets_reported)} of ${int(s.datasets_curated)}`}</td>
            <td className="px-2 py-1.5 text-right tabular-nums">{s.genomes_analyzed === 0 ? <span className="text-muted-text">–</span> : int(s.genomes_analyzed)}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
};

export const corpusFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const curated = m.species.reduce((a, s) => a + s.datasets_curated, 0);
  const genomes = m.species.reduce((a, s) => a + s.genomes_analyzed, 0);
  const structured = m.species.filter((s) => s.genomes_analyzed > 0).length;
  return {
    id: "corpus",
    title: "Corpus",
    result: rows.length === 0 ? `No dataset has a report yet; ${plural(curated, "dataset")} passed curation.` : undefined,
    figure: <CorpusFigure manifest={m} large={!opts.page} />,
    caption: `A dataset is one species and one antibiotic with ${datasetRule(classMinimum(rows))}.`,
    details: opts.page ? (
      <>
        <p className="prose-measure">
          {int(genomes)} genomes from {plural(structured, "species", "species")} are counted into k-mers and clustered into lineages.
          {curated > rows.length ? ` ${int(rows.length)} of ${int(curated)} curated datasets ${rows.length === 1 ? "has" : "have"} a report.` : ""}
        </p>
        <SpeciesCounts manifest={m} />
      </>
    ) : undefined,
  };
};

// ---- lineage cost -------------------------------------------------------------

export const lineageFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const s = splitSummary(rows);
  const unpaired = rows.filter((r) => r.best_random !== null && r.best_lineage === null);
  const byReason = new Map<string, CatalogRow[]>();
  unpaired.forEach((r) => {
    const k = splitReasonShort(r.lineage_unavailable_reason) ?? "No reason recorded";
    byReason.set(k, [...(byReason.get(k) ?? []), r]);
  });
  return {
    id: "lineage",
    title: "Lineage cost",
    question: "How much does a score depend on close relatives in training?",
    result:
      s.nLineage === 0 || s.medianGap === null
        ? "No dataset has a lineage-held-out result yet."
        : `${int(s.nLargeGap)} of ${plural(s.nLineage, "paired dataset")} lose more than 5 points with lineages held out.`,
    figure: (fit: Fit) => <GapSwarm rows={rows} height={fit.height} onOpen={opts.onOpen} />,
    caption:
      s.nLineage === 0
        ? "The lineage split holds out whole lineages so no test genome has a near-identical relative in training."
        : `Each dot is one dataset: balanced accuracy on random splits minus balanced accuracy with whole lineages held out, for ${BEST_LEARNER_WORDS}. Medians: ${pct(s.medianRandom)} random, ${pct(s.medianLineage)} lineage, ${((s.medianGap ?? 0) * 100).toFixed(1)} points lost.`,
    details: opts.page ? (
      <>
        <p className="prose-measure">
          Balanced accuracy is the mean of sensitivity and specificity over pooled out-of-fold predictions.
          {unpaired.length > 0 ? ` ${plural(unpaired.length, "dataset")} ${unpaired.length === 1 ? "has" : "have"} no lineage-held-out result and ${unpaired.length === 1 ? "is" : "are"} left out of all three medians:` : ""}
          {noBestWords(s.nNoBest)}
        </p>
        {unpaired.length > 0 ? (
          <ul className="list-disc space-y-1 pl-5">
            {[...byReason.entries()].map(([reason, group]) => (
              <li key={reason}>
                <span className="font-medium">{reason}</span>: {listWords(group.map((r) => `${r.antibiotic} in ${speciesShort(r.species)}`))}
                {group[0]?.lineage_unavailable_reason ? <span className="text-muted-text"> (“{group[0].lineage_unavailable_reason}”{group.length > 1 ? " or the same with another fold" : ""})</span> : null}
              </li>
            ))}
          </ul>
        ) : null}
      </>
    ) : undefined,
  };
};

// ---- known-gene catalog --------------------------------------------------------

/** Per drug class, the median gap between the learner (lineages held out)
 *  and the catalog rule, one mark per class on a points axis with the
 *  cutoff marked. Overview's Method and coverage only. */
const ClassStrip = ({ rows }: { rows: readonly ClassRow[] }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const shown = rows.filter((c) => c.medianGap !== null);
  const labelW = Math.min(200, Math.round(width * 0.42));
  const gaps = shown.map((c) => c.medianGap!);
  const lo = Math.min(-0.02, Math.floor(Math.min(0, ...gaps) * 100) / 100);
  const hi = Math.max(0.05, Math.ceil(Math.max(...gaps) * 100) / 100);
  const x = linear([lo, hi], [labelW + 8, width - 8]);
  const rowH = 26;
  const top = 18;
  const height = top + shown.length * rowH + 22;
  return (
    <div ref={frame}>
      <svg role="img" aria-label={`Median gap between the learner and the catalog rule for ${shown.length} drug classes.`} width={width} height={height} className="block">
        <text x={x(CATALOG_CUTOFF) + 4} y={12} className={LABEL_CLASS}>
          {Math.round(CATALOG_CUTOFF * 100)} points
        </text>
        <RefLine x={x(CATALOG_CUTOFF)} y0={top - 2} y1={top + shown.length * rowH} />
        <line x1={x(0)} x2={x(0)} y1={top} y2={top + shown.length * rowH} stroke="var(--chart-axis)" />
        {shown.map((c, i) => {
          const cy = top + i * rowH + rowH / 2;
          return (
            <g key={c.drugClass}>
              <text x={labelW} y={cy + 4} textAnchor="end" className={LABEL_CLASS}>
                {`${capitalize(c.drugClass)} · ${c.compared}`}
              </text>
              <circle cx={x(c.medianGap!)} cy={cy} r={4} fill="var(--chart-mark)" stroke="rgb(var(--paper-pure))" strokeWidth={1.5} />
            </g>
          );
        })}
        <AxisX x={x} y={top + shown.length * rowH} height={0} values={ticks(lo, hi, 5)} format={(v) => (v * 100).toFixed(0)} />
      </svg>
      <p className="text-center text-tiny text-muted-text">Median points the learner leads the catalog rule, per drug class (count of datasets beside the class)</p>
    </div>
  );
};

export const catalogFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const c = catalogSummary(rows);
  const cut = Math.round(CATALOG_CUTOFF * 100);
  const classes = classRows(rows);
  const matching = classes.filter((k) => k.verdict === "catalog").map((k) => k.drugClass);
  const noRule = classes.filter((k) => k.verdict === "no-rule").map((k) => k.drugClass);
  return {
    id: "catalog",
    title: "Known-gene catalog",
    result:
      c.nCompared === 0
        ? "No dataset has both a catalog rule and a lineage-held-out result yet."
        : `${int(c.within.get(CATALOG_CUTOFF) ?? 0)} of ${plural(c.nCompared, "catalog rule")} are within ${cut} points or ahead of the learner.`,
    figure: (fit: Fit) => <CatalogScatter rows={rows} height={fit.height} onOpen={opts.onOpen} />,
    caption: `Each mark compares a fixed AMRFinderPlus gene rule with ${BEST_LEARNER_WORDS}, scored with lineages held out.`,
    details: opts.page ? (
      <>
        <p className="prose-measure">
          Within 1 point or ahead: {int(c.within.get(0.01) ?? 0)}; within 5 points or ahead: {int(c.within.get(0.05) ?? 0)}.
          {rows.length > c.nWithRule ? ` ${plural(rows.length - c.nWithRule, "dataset")} test drugs with no acquired-gene rule and are not drawn.` : ""}
          {noBestWords(rows.filter((r) => r.catalog !== null && r.best_learner === null).length, " with a catalog rule")} The rule tests only for acquired genes, so mutational resistance counts against it.
        </p>
        <p className="prose-measure">
          Per drug class, the median dataset&apos;s catalog score is within {cut} points of the learner for {matching.length === 0 ? "no class" : listWords(matching)}.
          {noRule.length > 0 ? ` No acquired-gene rule: ${listWords(noRule)}.` : ""}
        </p>
        <ClassStrip rows={classes} />
      </>
    ) : undefined,
  };
};

// ---- lead PCR panel -------------------------------------------------------------

/** The lead panel's datasets, one row each: the filled mark is the share of
 *  resistant isolates with a panel call whose genome carries a target, the
 *  ring the share the best learner calls resistant, with both pooled shares
 *  as reference lines. */
export const PanelFigure = ({ rows, height, onOpen }: { rows: readonly PanelRowJoined[]; height?: number; onOpen?(id: string): void }) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const p = pool(rows);
  const panelShare = share(p.caught, p.calledResistant);
  const modelShare = share(p.modelCaught, p.modelResistant);
  const labelW = Math.min(220, Math.round(width * 0.5));
  const top = 2 * 16 + 4;
  const axis = 22;
  const rowH = height === undefined ? 30 : Math.max(18, Math.min(30, Math.floor((height - top - axis - 52) / Math.max(rows.length, 1))));
  const plotTop = top;
  const plotBottom = top + rows.length * rowH;
  const x = linear([0, 1], [labelW + 10, width - 20]);
  const at = (r: PanelRowJoined) => ({ panel: share(r.caught, r.calledResistant) ?? 0, model: share(r.modelCaught, r.modelResistant) ?? 0 });
  const refs = placeLabels(
    [
      ...(panelShare === null ? [] : [{ key: "panel", text: `Pooled panel ${pct(panelShare)}`, x: x(panelShare), anchor: "end" as const }]),
      ...(modelShare === null ? [] : [{ key: "model", text: `Pooled learner ${pct(modelShare)}`, x: x(modelShare), anchor: "start" as const }]),
    ],
    width,
    2,
  );
  const marks: ListMark[] = rows.map((r, i) => {
    const v = at(r);
    return {
      key: r.dataset,
      x: x(v.panel),
      y: plotTop + i * rowH + rowH / 2,
      tip: () => [
        { term: `${capitalize(r.antibiotic)}, ${speciesShort(r.species)}`, text: `Targets found in ${int(r.caught)} of ${int(r.calledResistant)} resistant isolates with a panel call (${pct(v.panel)}); the learner calls ${int(r.modelCaught)} of ${int(r.modelResistant)} resistant (${pct(v.model)}).` },
        ...(r.noCallResistant + r.noCallSusceptible > 0 ? [{ source: `No probe call: ${int(r.noCallResistant)} resistant, ${int(r.noCallSusceptible)} susceptible.` }] : []),
        ...(onOpen === undefined ? [] : [{ source: "Enter or click opens the dataset." }]),
      ],
    };
  });
  return (
    <div ref={frame} className="min-w-0 max-w-[64rem]">
      <MarkListbox
        label={`${rows[0]?.name ?? "Panel"} by dataset: target-found share against the learner's called-resistant share, ${rows.length} datasets. Arrow keys move between datasets${onOpen === undefined ? "" : "; Enter opens one"}.`}
        marks={marks}
        width={width}
        onOpen={onOpen === undefined ? undefined : (i) => onOpen(rows[i]!.dataset)}
      >
        <svg aria-hidden="true" width={width} height={plotBottom + axis} className="block overflow-visible">
          <PlacedText labels={refs} top={0} />
          {panelShare !== null ? <RefLine x={x(panelShare)} y0={top - 4} y1={plotBottom} /> : null}
          {modelShare !== null ? <RefLine x={x(modelShare)} y0={top - 4} y1={plotBottom} /> : null}
          {rows.map((r, i) => {
            const v = at(r);
            const cy = plotTop + i * rowH + rowH / 2;
            return (
              <g key={r.dataset} data-dataset={r.dataset}>
                <line x1={labelW + 10} x2={width - 20} y1={cy} y2={cy} stroke="var(--chart-grid)" />
                <text x={labelW} y={cy + 4} textAnchor="end" className={LABEL_CLASS}>
                  {`${capitalize(r.antibiotic)}, ${speciesShort(r.species)}`}
                </text>
                <line x1={x(v.panel)} x2={x(v.model)} y1={cy} y2={cy} stroke="var(--chart-axis)" strokeWidth={2} />
                <circle cx={x(v.model)} cy={cy} r={4.5} fill="rgb(var(--paper-pure))" stroke="var(--chart-mark)" strokeWidth={2} />
                <circle cx={x(v.panel)} cy={cy} r={4.5} fill="var(--chart-mark)" stroke="rgb(var(--paper-pure))" strokeWidth={1.5} />
              </g>
            );
          })}
          <AxisX x={x} y={plotBottom} height={0} values={ticks(0, 1, width < 480 ? 4 : 5)} format={(v) => `${Math.round(v * 100)}%`} />
        </svg>
      </MarkListbox>
      <div className="mt-2 flex flex-wrap gap-x-5 gap-y-1 text-control text-ink">
        <span className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={12} height={12}>
            <circle cx={6} cy={6} r={4.5} fill="var(--chart-mark)" />
          </svg>
          Resistant isolates with a panel target
        </span>
        <span className="inline-flex items-center gap-1.5">
          <svg aria-hidden="true" width={12} height={12}>
            <circle cx={6} cy={6} r={4} fill="none" stroke="var(--chart-mark)" strokeWidth={2} />
          </svg>
          Called resistant by the best learner
        </span>
      </div>
    </div>
  );
};

export const panelsFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const lead = leadPanel(m);
  if (lead === null) {
    return {
      id: "panels",
      title: "Simulated PCR panels",
      result: "No simulated panel targets an analyzed drug yet.",
      figure: null,
      caption: "Panels cover carbapenems, extended-spectrum cephalosporins, methicillin, vancomycin, and colistin.",
      details: opts.page ? <PanelEvidence manifest={m} onOpen={opts.onOpen} /> : undefined,
    };
  }
  const p = pool(lead.rows);
  const same = p.calledResistant === p.modelResistant;
  return {
    id: "panels",
    title: `${lead.name} targets`,
    result: `Targets are found in ${pct(share(p.caught, p.calledResistant))} of ${int(p.calledResistant)} resistant isolates${same ? "" : " with a panel call"}; the learner calls ${pct(share(p.modelCaught, p.modelResistant))} of ${same ? "them" : int(p.modelResistant)} resistant.`,
    figure: (fit: Fit) => <PanelFigure rows={lead.rows} height={fit.height} onOpen={opts.onOpen} />,
    caption: `Each row is one of ${int(p.datasets)} datasets. Panel targets are found in assembled genomes by k-mer containment, and the learner is ${BEST_LEARNER_WORDS}.${p.noCallResistant > 0 ? ` ${plural(p.noCallResistant, "resistant isolate")} with no probe call ${p.noCallResistant === 1 ? "is" : "are"} left out of the panel share only.` : ""}${noBestWords(panelsWithoutBest(m), " with a simulated panel")}`,
    details: opts.page ? <PanelEvidence manifest={m} onOpen={opts.onOpen} /> : undefined,
  };
};

// ---- MIC prediction ---------------------------------------------------------------

export const micFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const s = micSummary(rows);
  return {
    id: "mic",
    title: "MIC prediction",
    result:
      s.n === 0
        ? "No dataset had enough MIC values for MIC regression."
        : `${int(s.nAbove90)} of ${plural(s.n, "dataset")} predict ≥90% of exact MICs within one doubling dilution.`,
    figure: s.n === 0 ? null : (fit: Fit) => <MicStrip rows={rows} height={fit.height} onOpen={opts.onOpen} />,
    caption: `Each dot is one dataset: the share of its exact MICs predicted within one doubling dilution, out of fold. Median ${pct(s.median)}; 90% is the essential-agreement bar for susceptibility test devices.`,
    details: opts.page ? (
      <p className="prose-measure">
        Censored gradient-boosted regression on log2 dilution steps, fitted on the MICs of genomes with a reported call and predicted out of fold.
        {rows.length > s.n ? ` ${plural(rows.length - s.n, "dataset")} ${rows.length - s.n === 1 ? "has" : "have"} no MIC regression: it needs at least 100 isolates with an MIC, 50 of them exact.` : ""}
      </p>
    ) : undefined,
  };
};

export const FINDINGS = [corpusFinding, lineageFinding, catalogFinding, panelsFinding, micFinding] as const;
