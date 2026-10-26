// The report's findings, shared by Overview and Present: the data behind the
// reports, what holding out whole lineages costs, GRM's known-gene rules, the
// lead PCR panel, and MIC prediction. Each is a plain title, a lead with its
// count and denominator, a figure that shows that number, and a caption that
// says what a mark is; the page adds short notes on coverage, visible inline.
// Every number comes from data/derive at render time.

import { useEffect, useId, useLayoutEffect, useRef, useState, type MouseEvent, type ReactNode } from "react";
import type { CatalogRow, Manifest } from "../data/schema";
import { capitalize, dateText, int, listWords, pct, plural, points, speciesShort } from "../data/format";
import {
  CATALOG_CUTOFF,
  catalogSummary,
  classMinimum,
  classRows,
  micSummary,
  noBestWords,
  speciesCount,
  splitSummary,
  type CatalogSummary,
  type ClassMinimum,
  type MicSummary,
  type SplitSummary,
} from "../data/derive/overview";
import { leadPanel, panelsWithoutBest, pool, share, type PanelRowJoined } from "../data/derive/tests";
import { BEST_LEARNER_WORDS, splitReasonShort } from "../data/derive/dataset";
import { manifestAudit } from "../data/derive/evaluation";
import { routeHref, type Navigator } from "../app/router";
import { AxisX, LABEL_CLASS, linear, MarkListbox, PlacedText, placeLabels, RefLine, ticks, useWidth, type ListMark } from "../charts/kit";
import { GapSwarm, MicStrip } from "../charts/DotStrip";
import { CatalogScatter } from "../charts/CatalogScatter";
import { FOCUS_RING, ScrollFrame, type Fit } from "../ui/primitives";
import { PanelEvidence } from "./PanelEvidence";

export type FindingId = "corpus" | "lineage" | "catalog" | "panels" | "mic";

/** What a deck gives a figure that may continue onto further slides: the
 *  measured box, the page to draw and the page count, and, while the deck
 *  measures, where to report how many pages the figure needs at that size.
 *  A page passes none of these and gets every row. */
export interface PagedFit extends Fit {
  page?: number;
  pages?: number;
  onPages?(pages: number): void;
}

export type FindingFigure = ReactNode | ((fit: PagedFit) => ReactNode);

export interface Finding {
  id: FindingId;
  /** A plain section name. */
  title: string;
  /** The answer with its count and denominator: the page's lead, or a
   *  slide's one-sentence headline. */
  result: string;
  figure: FindingFigure;
  caption: string;
  /** Page only: short notes on coverage under the caption. */
  note?: ReactNode;
  /** Page only: a centered plot-and-legend unit, or a strip across the
   *  full width. */
  span: "unit" | "wide";
  /** Page only: the evidence behind the figure, under one explicit control. */
  more?: { label: string; body: ReactNode };
  /** Deck only: the figure continues onto further slides when its rows do
   *  not fit at their minimum height. */
  continues?: boolean;
}

export interface FindingOptions {
  /** Overview builds the page's lead and notes; Present a slide's headline. */
  page: boolean;
  /** Page only: marks and names open their reports, and the corpus counts
   *  link to their sources. Present passes none, so the deck stays static. */
  nav?: Navigator;
}

// ---- links and the status line -------------------------------------------------

/** A link's look: blue words, underlined, the focus ring on keyboard focus.
 *  The weight is the caller's. */
export const LINK = `rounded-sm text-blue underline decoration-blue/40 underline-offset-2 transition-colors duration-fast ease-out hover:decoration-blue ${FOCUS_RING}`;

/** A plain click on an in-app link navigates in place; a modified or middle
 *  click keeps the browser's own behavior, such as a new tab. */
export const inApp =
  (go: () => void) =>
  (e: MouseEvent<HTMLAnchorElement>): void => {
    if (e.defaultPrevented || e.button !== 0 || e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
    e.preventDefault();
    go();
  };

/** Data & Methods at one of its sections; MethodsView scrolls to the hash
 *  and moves focus to the section's heading. */
export const methodsHref = (anchor: string): string => `${routeHref({ view: "methods", species: null, dataset: null })}#${anchor}`;

export const openMethodsAt = (nav: Navigator, anchor: string): void => {
  // A router that takes a hash pushes it with the path, keeping this entry's
  // place for Back; one that does not gets it written onto the new entry.
  const next = { view: "methods" as const, species: null, hash: anchor };
  nav.go(next);
  if (window.location.hash !== `#${anchor}`) window.history.replaceState(window.history.state, "", methodsHref(anchor));
};

/** Scrolls a section of the current view into place and moves focus to its
 *  heading, so the keyboard continues from there. */
export const showSection = (id: string, smooth = false): void => {
  const target = document.getElementById(id);
  if (target === null) return;
  const reduced = typeof window.matchMedia === "function" && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  target.scrollIntoView?.({ block: "start", behavior: smooth && !reduced ? "smooth" : "auto" });
  target.querySelector<HTMLElement>("h2[tabindex], h3[tabindex]")?.focus({ preventScroll: true });
};

/** A report section's heading: its number, then its plain name. A link that
 *  jumps to the section moves focus here; it is never a tab stop. */
export const SectionHeading = ({ id, n, children }: { id: string; n: number; children: ReactNode }) => (
  <h2 id={`${id}-title`} tabIndex={-1} className={`flex items-baseline gap-3 rounded-sm text-headline font-semibold tracking-tight text-ink ${FOCUS_RING}`}>
    <span className="tabular-nums text-muted-text">{n}</span> <span>{children}</span>
  </h2>
);

const Sep = () => <span aria-hidden="true">·</span>;

/** The status line's link to where the scope is spelled out: the line's own
 *  muted words, underlined. */
export const STATUS_LINK = `rounded-sm underline decoration-current/40 underline-offset-2 transition-colors duration-fast ease-out hover:text-blue hover:decoration-current ${FOCUS_RING}`;

/** The report's one status line, in plain words: when it was built, the
 *  catalog version, and the scope every number shares. A provisional audit
 *  names itself in bold instead; AuditReasons lists why. */
export const ReportLine = ({ manifest, scope }: { manifest: Manifest; scope?: ReactNode }) => {
  const audit = manifestAudit(manifest);
  return (
    <>
      <span>Report built {dateText(manifest.generated_at)}</span> <Sep />{" "}
      <span>AMRFinderPlus catalog {manifest.amrfinder_version === "" ? "version not recorded" : manifest.amrfinder_version}</span>
      {audit.total === 0 ? null : (
        <>
          {" "}
          <Sep /> {audit.invalid > 0 ? <strong className="font-semibold text-coral-text">Provisional metrics</strong> : (scope ?? <span>Internal cross-validation only</span>)}
        </>
      )}
    </>
  );
};

/** Why the metrics are provisional, under the status line: the count and
 *  each distinct reason, as text. Nothing when every evaluation passes. */
export const AuditReasons = ({ manifest }: { manifest: Manifest }) => {
  const audit = manifestAudit(manifest);
  if (audit.total === 0 || audit.invalid === 0) return null;
  return (
    <div role="note" className="mt-4 max-w-3xl text-note text-ink">
      <p>
        <strong className="font-semibold">Provisional metrics.</strong>{" "}
        {audit.invalid === audit.total
          ? "No dataset evaluation in this report passes the evaluation audit."
          : `${int(audit.invalid)} of ${int(audit.total)} dataset evaluations do not pass the evaluation audit; the rest are internal cross-validation only.`}{" "}
        The numbers describe agreement with lab calls on these sampled isolates; they do not validate a diagnostic test.
      </p>
      <ul className="mt-1.5 list-disc space-y-0.5 pl-5 text-muted-text">
        {audit.reasons.map((r) => (
          <li key={r}>{r}</li>
        ))}
      </ul>
    </div>
  );
};

// ---- the report's words ------------------------------------------------------------

const sum = (values: readonly number[]): number => values.reduce((a, b) => a + b, 0);

const WORDS = ["No", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten"];

/** A count that opens a sentence, spelled out up to ten: "Six datasets". */
export const spelled = (n: number, one: string, many = `${one}s`): string => `${WORDS[n] ?? int(n)} ${n === 1 ? one : many}`;

/** A count inside a sentence, spelled out up to ten: "four because". */
const spelledLower = (n: number): string => (WORDS[n] ?? int(n)).toLowerCase();

/** "237 reports cover 13 species and use 21,889 genomes selected from
 *  119,072 QC-passing genomes." */
export const cohortLead = (m: Manifest): string => {
  const reports = m.datasets.length;
  const one = reports === 1;
  return `${plural(reports, "report")} ${one ? "covers" : "cover"} ${plural(speciesCount(m), "species", "species")} and ${one ? "uses" : "use"} ${int(sum(m.species.map((s) => s.genomes_analyzed)))} genomes selected from ${int(m.funnel.genomes_passing_qc)} QC-passing genomes.`;
};

/** What a report needs, in one sentence: the class minimum every report
 *  records, or each minimum with its report count when they differ. */
export const minimumSentence = (c: ClassMinimum): string => {
  const genomes = "selected genomes with unconflicted reported R/S calls";
  if (c.kind === "value") return `A report needs ${int(c.value)} resistant and ${int(c.value)} susceptible ${genomes}.`;
  if (c.kind === "unavailable") return `A report needs enough resistant and susceptible ${genomes}; no report records the minimum yet.`;
  return `A report needs enough resistant and susceptible ${genomes}: ${listWords(c.reports.map(([minimum, n]) => `at least ${int(minimum)} of each in ${plural(n, "report")}`))}.`;
};

/** The most genomes any species contributes, and whether that reads as a
 *  cap: some species sits exactly there with more labeled genomes than it
 *  used. */
export const selectionCap = (m: Manifest): { cap: number; capped: boolean } => {
  const analyzed = m.species.filter((s) => s.genomes_analyzed > 0);
  const cap = Math.max(0, ...analyzed.map((s) => s.genomes_analyzed));
  return { cap, capped: analyzed.some((s) => s.genomes_analyzed === cap && s.genomes_labeled > cap) };
};

/** What this run covers and why a curated dataset may have no report, in
 *  the four sentences a reader needs before reading any count. */
export const coverageNote = (m: Manifest): string => {
  const analyzed = m.species.filter((s) => s.genomes_analyzed > 0);
  const waiting = m.species.filter((s) => s.genomes_analyzed === 0 && s.datasets_curated > 0);
  const unreported = sum(analyzed.map((s) => Math.max(0, s.datasets_curated - s.datasets_reported)));
  const { cap, capped } = selectionCap(m);
  const selection = capped
    ? `In this run, selection is capped at ${int(cap)} genomes per species when more are eligible, and QC and FASTA availability reduce some cohorts further.`
    : "In this run, QC and FASTA availability decide which labeled genomes are selected.";
  const missing =
    unreported === 0
      ? "Every curated dataset in these species has a report."
      : `${spelled(unreported, "curated dataset")} in these species ${unreported === 1 ? "has" : "have"} no report.`;
  const outside =
    waiting.length === 0
      ? ""
      : ` Another ${plural(waiting.length, "species", "species")} (${plural(sum(waiting.map((s) => s.datasets_curated)), "curated dataset")}) ${waiting.length === 1 ? "was" : "were"} not analyzed; ${waiting.length === 1 ? "its" : "their"} FASTA files are staged separately, not in the pipeline's data root.`;
  return `${selection} ${minimumSentence(classMinimum(m.datasets))} ${missing}${outside}`;
};

/** The lineage lead: both medians, the datasets past 5 points, and the
 *  median loss, over the paired datasets. */
export const lineageLead = (s: SplitSummary): string => {
  if (s.nLineage === 0 || s.medianGap === null || s.medianRandom === null || s.medianLineage === null) return "No dataset has a lineage-held-out result yet.";
  const verb = s.medianLineage <= s.medianRandom ? "lowers" : "raises";
  const lose = s.nLargeGap === 1 ? "loses" : "lose";
  const median = s.medianGap >= 0 ? `the median loss is ${points(s.medianGap)}` : `the median dataset gains ${points(-s.medianGap)}`;
  return `Holding out whole lineages ${verb} median balanced accuracy from ${pct(s.medianRandom)} to ${pct(s.medianLineage)}. ${int(s.nLargeGap)} of ${plural(s.nLineage, "dataset")} ${lose} more than 5 points; ${median}.`;
};

/** The catalog lead: within the cutoff or ahead, split into ahead and
 *  within, then the rest, which trail by more than the cutoff. */
export const catalogLead = (c: CatalogSummary): string => {
  if (c.nCompared === 0) return "No dataset has both a known-gene rule and a lineage-held-out result yet.";
  const cut = Math.round(CATALOG_CUTOFF * 100);
  const { ahead, within, behind } = c.regions;
  return `GRM's known-gene rules are within ${cut} points or ahead on ${int(c.within.get(CATALOG_CUTOFF) ?? 0)} of ${plural(c.nCompared, "dataset")}: ${int(ahead)} ahead by more than ${cut} and ${int(within)} within ${cut}. They trail by more than ${cut} on the other ${int(behind)}.`;
};

/** The MIC lead: datasets at the 90% bar, then the median dataset. */
export const micLead = (s: MicSummary): string =>
  s.n === 0
    ? "No dataset had enough MIC values for MIC regression."
    : `For ${int(s.nAbove90)} of ${plural(s.n, "dataset")}, at least 90% of exact MICs are predicted within one doubling dilution (essential agreement). The median dataset reaches ${pct(s.median)}.`;

// ---- data ------------------------------------------------------------------------

/** The four counts the reports rest on, each with its unit. The deck shows
 *  them static; the page gives each its source in CorpusTiles. */
export const CorpusFigure = ({ manifest, large = false }: { manifest: Manifest; large?: boolean }) => {
  const f = manifest.funnel;
  const stages = [
    { value: f.pairs_reported, unit: "Reported R/S genome–antibiotic pairs" },
    { value: f.genomes_passing_qc, unit: "Genomes passed QC" },
    { value: sum(manifest.species.map((s) => s.genomes_analyzed)), unit: "Genomes analyzed" },
    { value: manifest.datasets.length, unit: manifest.datasets.length === 1 ? "Report" : "Reports" },
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

/** A tile's action: visible words, a 44px target, and an arrow that says
 *  where it goes. The disclosure button reads as a control, not a link. */
const ACTION = `inline-flex min-h-11 items-center gap-1.5 text-control font-medium ${LINK}`;
const DISCLOSE = `inline-flex min-h-11 items-center gap-1.5 rounded-sm text-control font-medium text-blue transition-colors duration-fast ease-out hover:text-ink ${FOCUS_RING}`;

/** Each species' reports and genomes: the selected genomes it contributes
 *  beside its genomes with any label, a different and larger denominator. */
const SpeciesCounts = ({ manifest, nav, id }: { manifest: Manifest; nav: Navigator; id: string }) => {
  const rows = manifest.species.filter((s) => s.datasets_curated > 0 || s.genomes_analyzed > 0);
  const analyzed = rows.filter((s) => s.genomes_analyzed > 0);
  const waiting = rows.filter((s) => s.genomes_analyzed === 0);
  const TH = "px-2 py-2 text-tiny font-medium text-blue";
  const reports = (s: (typeof rows)[number]) =>
    s.datasets_reported === s.datasets_curated ? int(s.datasets_reported) : `${int(s.datasets_reported)} of ${int(s.datasets_curated)}`;
  return (
    <div id={id} className="mt-6">
      <ScrollFrame label="Genomes by species">
        <table className="w-full border-collapse text-control">
          <thead>
            <tr className="border-b border-hair">
              <th scope="col" className={`${TH} text-left`}>
                Species
              </th>
              <th scope="col" className={`${TH} text-right`}>
                Reports
              </th>
              <th scope="col" className={`${TH} text-right`}>
                Genomes analyzed
              </th>
              <th scope="col" className={`${TH} text-right`}>
                Genomes with any label
              </th>
            </tr>
          </thead>
          <tbody>
            {analyzed.map((s) => (
              <tr key={s.slug} className="border-b border-hair">
                <th scope="row" className="px-2 py-1 text-left font-normal">
                  <a
                    href={routeHref({ view: "species", species: s.slug, dataset: null })}
                    onClick={inApp(() => nav.go({ view: "species", species: s.slug }))}
                    className={`ui-control-target inline-flex items-center italic ${LINK}`}
                  >
                    {s.species}
                  </a>
                </th>
                <td className="px-2 py-1 text-right tabular-nums">{reports(s)}</td>
                <td className="px-2 py-1 text-right tabular-nums">{int(s.genomes_analyzed)}</td>
                <td className="px-2 py-1 text-right tabular-nums">{int(s.genomes_labeled)}</td>
              </tr>
            ))}
          </tbody>
          {waiting.length === 0 ? null : (
            <tbody>
              <tr className="border-b border-hair">
                <th scope="rowgroup" colSpan={4} className="px-2 pb-1.5 pt-4 text-left text-tiny font-medium text-muted-text">
                  Not analyzed in this run
                </th>
              </tr>
              {waiting.map((s) => (
                <tr key={s.slug} className="border-b border-hair text-muted-text">
                  <th scope="row" className="px-2 py-1.5 text-left font-normal italic">
                    {s.species}
                  </th>
                  <td className="px-2 py-1.5 text-right tabular-nums">{`0 of ${int(s.datasets_curated)}`}</td>
                  <td className="px-2 py-1.5 text-right">None</td>
                  <td className="px-2 py-1.5 text-right tabular-nums">{int(s.genomes_labeled)}</td>
                </tr>
              ))}
            </tbody>
          )}
        </table>
      </ScrollFrame>
      <p className="type-note mt-2">
        A genome with any label has a curated label for at least one drug, including labels inferred from an MIC. A report&apos;s cohort is smaller: the selected genomes with an
        unconflicted reported R or S call for its drug.
      </p>
    </div>
  );
};

/** The four counts with a way into each: the curation steps behind the
 *  first two, the species behind the genomes, and every report. Native links
 *  and one disclosure button, each with its words in view. */
const CorpusTiles = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const [open, setOpen] = useState(false);
  const tableId = useId();
  const f = manifest.funnel;
  const analyzed = sum(manifest.species.map((s) => s.genomes_analyzed));
  const methods = (label: string, anchor: string) => (
    <a href={methodsHref(anchor)} onClick={inApp(() => openMethodsAt(nav, anchor))} className={ACTION}>
      {label}
      <span aria-hidden="true">→</span>
    </a>
  );
  const tiles: { value: number; unit: string; action: ReactNode }[] = [
    { value: f.pairs_reported, unit: `reported R/S genome–antibiotic pairs, from ${int(f.records)} BV-BRC records`, action: methods("See pair curation", "m-funnel") },
    { value: f.genomes_passing_qc, unit: "genomes passed QC", action: methods("See QC outcomes", "m-funnel") },
    {
      value: analyzed,
      unit: "genomes analyzed",
      action: (
        <button type="button" aria-expanded={open} aria-controls={open ? tableId : undefined} onClick={() => setOpen((v) => !v)} className={DISCLOSE}>
          {open ? "Hide genomes by species" : "Show genomes by species"}
          <span aria-hidden="true">{open ? "↑" : "↓"}</span>
        </button>
      ),
    },
    {
      value: manifest.datasets.length,
      unit: manifest.datasets.length === 1 ? "report" : "reports",
      action: (
        <a href={routeHref({ view: "datasets", species: null, dataset: null })} onClick={inApp(() => nav.go({ view: "datasets", species: null }))} className={ACTION}>
          Browse all report datasets
          <span aria-hidden="true">→</span>
        </a>
      ),
    },
  ];
  return (
    <div>
      <ul className="grid grid-cols-1 gap-x-8 gap-y-6 min-[480px]:grid-cols-2 lg:grid-cols-4">
        {tiles.map((t) => (
          <li key={t.unit} className="flex min-w-0 flex-col border-t border-hair pt-3">
            <span className="text-stat font-semibold tabular-nums text-ink">{int(t.value)}</span> <span className="mt-1 text-note text-muted-text">{t.unit}</span>{" "}
            <span className="mt-auto pt-1">{t.action}</span>
          </li>
        ))}
      </ul>
      {open ? <SpeciesCounts manifest={manifest} nav={nav} id={tableId} /> : null}
    </div>
  );
};

export const corpusFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const curated = sum(m.species.map((s) => s.datasets_curated));
  const { cap, capped } = selectionCap(m);
  return {
    id: "corpus",
    title: "Data",
    result: m.datasets.length === 0 ? `No dataset has a report yet; ${plural(curated, "dataset")} passed curation.` : cohortLead(m),
    figure: opts.page && opts.nav !== undefined ? <CorpusTiles manifest={m} nav={opts.nav} /> : <CorpusFigure manifest={m} large={!opts.page} />,
    caption: opts.page
      ? coverageNote(m)
      : `${capped ? `Selection is capped at ${int(cap)} genomes per species when more are eligible. ` : ""}${minimumSentence(classMinimum(m.datasets))}`,
    span: "unit",
  };
};

// ---- unseen lineages ---------------------------------------------------------------

/** Why a dataset has no lineage-held-out result, in words a reader can use;
 *  the pipeline's own sentence stays in its report. */
const UNPAIRED_WORDS: Record<string, { clause: string; label: string }> = {
  "Fold lacks a class": { clause: "a held-out fold would hold only one class", label: "A fold with one class" },
  "Too few lineages": { clause: "there are too few lineages", label: "Too few lineages" },
  Unavailable: { clause: "the split could not run", label: "Split could not run" },
  none: { clause: "no reason is recorded", label: "No reason recorded" },
};

const datasetName = (r: { antibiotic: string; species: string }): string => `${capitalize(r.antibiotic)} in ${speciesShort(r.species)}`;

/** The datasets with a random-split result but no lineage-held-out one,
 *  counted by why, each name opening its report. */
const UnpairedNote = ({ rows, nav }: { rows: readonly CatalogRow[]; nav?: Navigator }) => {
  const unpaired = rows.filter((r) => r.best_random !== null && r.best_lineage === null);
  if (unpaired.length === 0) return null;
  const groups = new Map<string, CatalogRow[]>();
  unpaired.forEach((r) => {
    const k = splitReasonShort(r.lineage_unavailable_reason) ?? "none";
    groups.set(k, [...(groups.get(k) ?? []), r]);
  });
  const entries = [...groups.entries()].sort((a, b) => b[1].length - a[1].length);
  const words = (k: string) => UNPAIRED_WORDS[k] ?? UNPAIRED_WORDS.Unavailable!;
  const clauses = entries.map(([k, g]) => `${spelledLower(g.length)} because ${words(k).clause}`);
  return (
    <div className="prose-measure mt-3 space-y-1.5 text-note text-ink">
      <p>
        {spelled(unpaired.length, "dataset")} {unpaired.length === 1 ? "has" : "have"} no paired lineage result and {unpaired.length === 1 ? "is" : "are"} left out of these medians
        {entries.length === 1 ? `, because ${words(entries[0]![0]).clause}` : `: ${listWords(clauses)}`}.
      </p>
      {entries.map(([k, g]) => (
        <p key={k}>
          <span className="text-muted-text">{words(k).label}: </span>
          {g.map((r, i) => (
            <span key={r.id}>
              {i > 0 ? ", " : null}
              {nav === undefined ? (
                datasetName(r)
              ) : (
                <a href={routeHref({ ...nav.route, dataset: r.id })} onClick={inApp(() => nav.openDataset(r.id))} className={LINK}>
                  {datasetName(r)}
                </a>
              )}
            </span>
          ))}
        </p>
      ))}
    </div>
  );
};

export const lineageFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const s = splitSummary(rows);
  const none = s.nLineage === 0 || s.medianGap === null;
  return {
    id: "lineage",
    title: "Unseen lineages",
    result: opts.page || none ? lineageLead(s) : `${int(s.nLargeGap)} of ${plural(s.nLineage, "dataset")} lose more than 5 points of balanced accuracy when whole lineages are held out.`,
    figure: (fit: PagedFit) => <GapSwarm rows={rows} height={fit.height} onOpen={opts.nav?.openDataset} />,
    caption: none
      ? "The lineage split holds out whole lineages, so no test genome has a near-identical relative in training."
      : `Each dot is one dataset: balanced accuracy, the mean of sensitivity and specificity, on random splits minus balanced accuracy with whole lineages held out, for ${BEST_LEARNER_WORDS}.${
          opts.page ? "" : ` Medians: ${pct(s.medianRandom)} random, ${pct(s.medianLineage)} lineage, ${((s.medianGap ?? 0) * 100).toFixed(1)} points lost.`
        }`,
    note: opts.page ? (
      <>
        <UnpairedNote rows={rows} nav={opts.nav} />
        {s.nNoBest > 0 ? <p className="prose-measure mt-1.5 text-note text-ink">{noBestWords(s.nNoBest).trim()}</p> : null}
      </>
    ) : undefined,
    span: "wide",
  };
};

// ---- known-gene rules ----------------------------------------------------------------

export const catalogFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const c = catalogSummary(rows);
  const cut = Math.round(CATALOG_CUTOFF * 100);
  const classes = classRows(rows);
  const matching = classes.filter((k) => k.verdict === "catalog").length;
  const noRuleClasses = classes.filter((k) => k.verdict === "no-rule").length;
  const ruled = classes.length - noRuleClasses;
  const noRule = rows.length - c.nWithRule;
  const unpaired = c.nWithRule - c.nCompared;
  const noBest = rows.filter((r) => r.catalog !== null && r.best_learner === null).length;
  return {
    id: "catalog",
    title: "Known-gene rules",
    result:
      opts.page || c.nCompared === 0 ? catalogLead(c) : `Known-gene rules are within ${cut} points or ahead on ${int(c.within.get(CATALOG_CUTOFF) ?? 0)} of ${plural(c.nCompared, "dataset")}.`,
    // On the page the plot may grow past its default size where the unit has room.
    figure: (fit: PagedFit) => <CatalogScatter rows={rows} height={fit.height ?? (opts.page ? 640 : undefined)} onOpen={opts.nav?.openDataset} />,
    // A slide says the same in fewer lines, so a phone's plot keeps the size
    // its 12px labels need.
    caption: opts.page
      ? "GRM detects acquired AMRFinderPlus genes with its own k-mer matching and applies its own drug rules; it does not run the AMRFinderPlus program or test point mutations, so mutational resistance counts against a rule. Each mark is one dataset: the rule's balanced accuracy, which needs no training, against the best learner's with whole lineages held out."
      : "GRM's own k-mer matching finds acquired AMRFinderPlus genes for its own drug rules; it does not run AMRFinderPlus or test point mutations. One mark per dataset.",
    note: opts.page ? (
      <div className="prose-measure mt-3 space-y-1.5 text-note text-ink">
        {noRule + unpaired + noBest > 0 ? (
          <p>
            {[
              ...(noRule > 0 ? [`${plural(noRule, "dataset")} ${noRule === 1 ? "tests a drug" : "test drugs"} with no acquired-gene rule`] : []),
              ...(unpaired > 0 ? [`${plural(unpaired, "dataset")} with a rule ${unpaired === 1 ? "has" : "have"} no lineage-held-out result`] : []),
            ].join(", and ")}
            {noRule + unpaired > 0 ? `; ${noRule + unpaired === 1 ? "it is" : "they are"} not drawn.` : ""}
            {noBestWords(noBest, " with a rule")}
          </p>
        ) : null}
        {ruled > 0 ? (
          <p>
            By drug class, the median dataset&apos;s rule is within {cut} points of the learner in {int(matching)} of {plural(ruled, "class", "classes")} with a rule
            {noRuleClasses > 0 ? `; ${plural(noRuleClasses, "class", "classes")} ${noRuleClasses === 1 ? "has" : "have"} no acquired-gene rule` : ""}.
          </p>
        ) : null}
      </div>
    ) : undefined,
    span: "unit",
  };
};

// ---- lead PCR panel -------------------------------------------------------------------

const PANEL_TOP = 2 * 16 + 4;
const PANEL_AXIS = 22;
/** A slide's rows never get shorter than this; labels keep their 12px. */
const PANEL_ROW_MIN = 18;
const PANEL_ROW_MAX = 30;
/** The continuation line under a figure that spans several slides. */
const PAGER_H = 24;

/** The lead panel's datasets, one row each: the filled mark is the share of
 *  resistant isolates with a panel call whose genome carries a target, the
 *  ring the share the best learner calls resistant, with both pooled shares
 *  as reference lines. On a slide too short for every row at the minimum row
 *  height, the rows continue onto further slides: each page draws its share
 *  of the rows against the same axis and the same pooled lines. */
export const PanelFigure = ({
  rows,
  height,
  onOpen,
  page,
  pages,
  onPages,
}: {
  rows: readonly PanelRowJoined[];
  height?: number;
  onOpen?(id: string): void;
  page?: number;
  pages?: number;
  onPages?(pages: number): void;
}) => {
  const frame = useRef<HTMLDivElement | null>(null);
  const legendRef = useRef<HTMLDivElement | null>(null);
  const width = useWidth(frame);
  const [legendH, setLegendH] = useState(28);
  useLayoutEffect(() => {
    const el = legendRef.current;
    if (el === null) return undefined;
    const read = () => setLegendH(el.offsetHeight);
    read();
    const observer = new ResizeObserver(read);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);
  // The plot's height for rows, with or without the continuation line.
  const room = (continued: boolean) => (height ?? Infinity) - PANEL_TOP - PANEL_AXIS - legendH - (continued ? PAGER_H : 0);
  const needed = height === undefined || room(false) >= rows.length * PANEL_ROW_MIN ? 1 : Math.ceil(rows.length / Math.max(1, Math.floor(room(true) / PANEL_ROW_MIN)));
  // An unlaid box measures nothing, so it reports nothing.
  useEffect(() => {
    if (height !== undefined && height > 0) onPages?.(needed);
  }, [onPages, needed, height]);
  const count = Math.max(1, Math.min(pages ?? 1, rows.length));
  const per = Math.ceil(rows.length / count);
  const at = Math.min(Math.max(0, page ?? 0), count - 1);
  const shown = count === 1 ? rows : rows.slice(at * per, at * per + per);
  const p = pool(rows);
  const panelShare = share(p.caught, p.calledResistant);
  const modelShare = share(p.modelCaught, p.modelResistant);
  const labelW = Math.min(220, Math.round(width * 0.5));
  const rowH = height === undefined ? PANEL_ROW_MAX : Math.max(PANEL_ROW_MIN, Math.min(PANEL_ROW_MAX, Math.floor(room(count > 1) / Math.max(shown.length, 1))));
  const plotTop = PANEL_TOP;
  const plotBottom = PANEL_TOP + shown.length * rowH;
  const x = linear([0, 1], [labelW + 10, width - 20]);
  const at01 = (r: PanelRowJoined) => ({ panel: share(r.caught, r.calledResistant) ?? 0, model: share(r.modelCaught, r.modelResistant) ?? 0 });
  const refs = placeLabels(
    [
      ...(panelShare === null ? [] : [{ key: "panel", text: `Pooled panel ${pct(panelShare)}`, x: x(panelShare), anchor: "end" as const }]),
      ...(modelShare === null ? [] : [{ key: "model", text: `Pooled learner ${pct(modelShare)}`, x: x(modelShare), anchor: "start" as const }]),
    ],
    width,
    2,
  );
  const first = count === 1 ? 1 : at * per + 1;
  const last = count === 1 ? rows.length : at * per + shown.length;
  // A deck's off-screen measuring copy draws the same plot with no keyboard
  // layer, so it adds no listbox or option to the page.
  const measuring = onPages !== undefined;
  const marks: ListMark[] = shown.map((r, i) => {
    const v = at01(r);
    return {
      key: r.dataset,
      x: x(v.panel),
      y: plotTop + i * rowH + rowH / 2,
      tip: () => [
        { term: `${capitalize(r.antibiotic)}, ${speciesShort(r.species)}`, text: `Targets found in ${int(r.caught)} of ${int(r.calledResistant)} resistant isolates with a panel call (${pct(v.panel)}); the learner calls ${int(r.modelCaught)} of ${int(r.modelResistant)} resistant (${pct(v.model)}).` },
        ...(r.noCallResistant + r.noCallSusceptible > 0 ? [{ source: `No probe call: ${int(r.noCallResistant)} resistant, ${int(r.noCallSusceptible)} susceptible.` }] : []),
      ],
    };
  });
  const plot = (
    <svg aria-hidden="true" width={width} height={plotBottom + PANEL_AXIS} className="block overflow-visible">
      <PlacedText labels={refs} top={0} />
      {panelShare !== null ? <RefLine x={x(panelShare)} y0={PANEL_TOP - 4} y1={plotBottom} /> : null}
      {modelShare !== null ? <RefLine x={x(modelShare)} y0={PANEL_TOP - 4} y1={plotBottom} /> : null}
      {shown.map((r, i) => {
        const v = at01(r);
        const cy = plotTop + i * rowH + rowH / 2;
        return (
          <g key={r.dataset} data-dataset={measuring ? undefined : r.dataset}>
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
  );
  return (
    <div ref={frame} className="min-w-0">
      {measuring ? (
        plot
      ) : (
        <MarkListbox
          label={`${rows[0]?.name ?? "Panel"} by dataset: target-found share against the learner's called-resistant share, ${count === 1 ? `${rows.length} datasets` : `datasets ${first} to ${last} of ${rows.length}`}. Arrow keys move between datasets; click or Enter to pin details; Escape clears.`}
          marks={marks}
          width={width}
          onOpen={onOpen === undefined ? undefined : (i) => onOpen(shown[i]!.dataset)}
        >
          {plot}
        </MarkListbox>
      )}
      <div ref={legendRef} className="flex flex-wrap gap-x-5 gap-y-1 pt-2 text-control text-ink">
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
      {count > 1 ? (
        <p data-continued="" className="flex h-6 items-end text-control text-muted-text">
          Datasets {int(first)}–{int(last)} of {int(rows.length)}
        </p>
      ) : null}
    </div>
  );
};

export const panelsFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const lead = leadPanel(m);
  const evidence = opts.page ? { label: "All panel results", body: <PanelEvidence manifest={m} onOpen={opts.nav?.openDataset} /> } : undefined;
  if (lead === null) {
    return {
      id: "panels",
      title: "PCR panel",
      result: "No simulated PCR panel targets an analyzed drug yet.",
      figure: null,
      caption: "Panels cover carbapenems, extended-spectrum cephalosporins, methicillin, vancomycin, and colistin.",
      span: "unit",
      more: evidence,
    };
  }
  const p = pool(lead.rows);
  const same = p.calledResistant === p.modelResistant;
  const found = pct(share(p.caught, p.calledResistant));
  const called = pct(share(p.modelCaught, p.modelResistant));
  return {
    id: "panels",
    title: "PCR panel",
    result: `${lead.name} targets are found in ${found} of ${int(p.calledResistant)} resistant isolates${same ? "" : " with a panel call"}${opts.page ? ` across ${plural(p.datasets, "dataset")}` : ""}; the best learner calls ${called} of ${same ? "them" : int(p.modelResistant)} resistant.`,
    figure: (fit: PagedFit) => <PanelFigure rows={lead.rows} height={fit.height} onOpen={opts.nav?.openDataset} page={fit.page} pages={fit.pages} onPages={fit.onPages} />,
    caption: `Each row is one of ${int(p.datasets)} datasets. Panel targets are found in assembled genomes by k-mer containment, in silico; the learner is ${BEST_LEARNER_WORDS}.${
      p.noCallResistant > 0 ? ` ${plural(p.noCallResistant, "resistant isolate")} with no probe call ${p.noCallResistant === 1 ? "is" : "are"} left out of the panel share only.` : ""
    }${noBestWords(panelsWithoutBest(m), " with a simulated panel")}`,
    span: "unit",
    more: evidence,
    continues: true,
  };
};

// ---- MIC prediction ---------------------------------------------------------------------

export const micFinding = (m: Manifest, opts: FindingOptions): Finding => {
  const rows = m.datasets;
  const s = micSummary(rows);
  const without = rows.length - s.n;
  return {
    id: "mic",
    title: "MIC prediction",
    result: opts.page || s.n === 0 ? micLead(s) : `For ${int(s.nAbove90)} of ${plural(s.n, "dataset")}, at least 90% of exact MICs are predicted within one doubling dilution.`,
    figure: s.n === 0 ? null : (fit: PagedFit) => <MicStrip rows={rows} height={fit.height} onOpen={opts.nav?.openDataset} />,
    caption: `Each dot is one dataset: the share of its exact MICs predicted within one doubling dilution, out of fold. An exact MIC is a single measured value, not a bound such as “≤0.25”. 90% is the essential-agreement bar used for susceptibility test devices.${
      opts.page ? "" : ` Median ${pct(s.median)}.`
    }`,
    note:
      opts.page && without > 0 ? (
        <p className="prose-measure mt-3 text-note text-ink">
          {plural(without, "dataset")} {without === 1 ? "has" : "have"} too few MICs for a model, which needs at least 100 isolates with an MIC, 50 of them exact.
        </p>
      ) : undefined,
    span: "wide",
  };
};

export const FINDINGS = [corpusFinding, lineageFinding, catalogFinding, panelsFinding, micFinding] as const;
