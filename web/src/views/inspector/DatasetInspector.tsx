// The dataset inspector: one near-full-canvas overlay with a fixed header
// and one vertical scroll owner. The header carries the title, one status
// row (the audit chip and the dataset's facts), and the section pills, which
// follow the scroll. The body opens on the learner comparison (how well
// resistance is predicted here, and by which learner), then the model and
// the rest of the evidence in order. Header and body share one column. The
// panel grows from the row that launched it.

import { useEffect, useLayoutEffect, useRef, useState, type RefObject } from "react";
import type { CatalogRow, DatasetReport, Manifest, Research } from "../../data/schema";
import { capitalize, distance, int, mgL, pct, plural, speciesShort } from "../../data/format";
import { loadErrorText, useDataset, useResearch, type Resource } from "../../data/load";
import { splitAbsence, splitOf } from "../../data/derive/dataset";
import { auditOf } from "../../data/derive/evaluation";
import { minClassOf } from "../../data/derive/overview";
import { AuditChip, DatasetNotice } from "../../ui/EvaluationNotice";
import type { Navigator } from "../../app/router";
import { ModalDialog } from "../../ui/ModalDialog";
import { Disclosure, FactRows, FOCUS_RING, GhostButton, PanelState } from "../../ui/primitives";
import { takeLaunch } from "../../ui/launch";
import { LearnerSection } from "./LearnerSection";
import { ModelSection } from "./ModelSection";
import { MechanismSection } from "./MechanismSection";
import { MicSection } from "./MicSection";
import { LineageSection } from "./LineageSection";

/** The section pills, in body order; each names its section's heading. */
const JUMPS = [
  { id: "ds-learners", label: "Learners" },
  { id: "ds-model", label: "Interpretable model" },
  { id: "ds-mechanisms", label: "Genes and panels" },
  { id: "ds-designer", label: "Marker designer" },
  { id: "ds-mic", label: "MIC" },
  { id: "ds-lineages", label: "Lineages" },
] as const;

/** The pill of the section the reader is in: the last whose top has passed
 *  the scroll owner's top edge, or the last section at the very bottom. */
const useCurrentSection = (scroller: RefObject<HTMLDivElement | null>, ready: boolean): string => {
  const [current, setCurrent] = useState<string>(JUMPS[0].id);
  useEffect(() => {
    const el = scroller.current;
    if (el === null || !ready) return undefined;
    const read = () => {
      const edge = el.getBoundingClientRect().top + 24;
      const present = JUMPS.filter((j) => document.getElementById(j.id) !== null);
      let at: string = present[0]?.id ?? JUMPS[0].id;
      for (const j of present) {
        if (document.getElementById(j.id)!.getBoundingClientRect().top <= edge) at = j.id;
      }
      if (el.scrollTop > 0 && el.scrollTop + el.clientHeight >= el.scrollHeight - 2) at = present[present.length - 1]?.id ?? at;
      setCurrent(at);
    };
    // Scroll events already arrive at most once per frame.
    read();
    el.addEventListener("scroll", read, { passive: true });
    return () => el.removeEventListener("scroll", read);
  }, [scroller, ready]);
  return current;
};

/** The header's MIC-label count: the labels left out when the report counts
 *  them, the samples labeled from MIC in a report from before the audit. */
const micLabelWords = (row: CatalogRow): string =>
  row.labels_from_mic_excluded !== undefined
    ? ` · ${plural(row.labels_from_mic_excluded, "MIC-inferred label")} left out`
    : row.labels_from_mic > 0
      ? ` · ${int(row.labels_from_mic)} labels inferred from MIC`
      : "";

const MethodFold = ({ r, current }: { r: DatasetReport; current?: number }) => {
  const audit = auditOf(r, current);
  const c = r.counts;
  const random = splitOf(r, "random");
  const lineage = splitOf(r, "lineage");
  // Reports from before the evaluation audit trained on MIC-inferred labels.
  const legacy = c.labels_from_mic > 0;
  const excluded = c.labels_from_mic_excluded;
  // Counted genomes with an R or S label for the drug: samples plus left out.
  const labeled = c.samples + (excluded ?? 0);
  const minimum = minClassOf(c.min_class);
  return (
    <Disclosure>
      <FactRows
        rows={[
          { label: "Evaluation audit", value: audit.valid ? `Passes, revision ${audit.revision}` : `Provisional. ${audit.reasons.join(" ")}` },
          { label: "Isolates", value: `${int(c.samples)}: ${int(c.resistant)} resistant, ${int(c.susceptible)} susceptible` },
          {
            label: "Class minimum",
            value: `${int(minimum)} resistant and ${int(minimum)} susceptible isolates${legacy ? "" : " with a reported call"}`,
          },
          legacy
            ? { label: "Samples labeled from MIC", value: `${int(c.labels_from_mic)} of ${int(c.samples)} (${pct(c.labels_from_mic / Math.max(c.samples, 1))})` }
            : {
                label: "MIC-inferred labels left out",
                value:
                  excluded === undefined
                    ? "Not recorded in this report"
                    : `${int(excluded)} of the ${int(labeled)} counted genomes with a resistant or susceptible label (${pct(excluded / Math.max(labeled, 1))})`,
              },
          {
            label: "Cohort breakpoint",
            value:
              r.breakpoint === null
                ? "No inferred breakpoint; calls as reported"
                : `Resistant above ${mgL(r.breakpoint.r_above)} mg/L, agreeing with ${pct(r.breakpoint.agreement)} of ${int(r.breakpoint.n)} reported calls; ${legacy ? "used to infer labels for isolates with an MIC but no call" : "it labels genomes with an MIC but no reported call, and this evaluation leaves those labels out"}`,
          },
          { label: "Features", value: `${int(c.features)} presence patterns carried by at least ${int(c.min_support)} isolates` },
          { label: "Patterns", value: int(c.patterns) },
          { label: "K-mers", value: `${int(c.kmers)} (${c.k === undefined ? "k not recorded" : `k = ${c.k}`})` },
          { label: "Lineages", value: `${int(c.lineages)} at Mash distance ${distance(r.lineage_threshold)}` },
          { label: "Random split", value: random === undefined ? `Not available: ${splitAbsence(r, "random")}` : `${random.folds} stratified folds` },
          { label: "Lineage split", value: lineage === undefined ? `Not available: ${splitAbsence(r, "lineage")}` : `${lineage.folds} folds over ${int(lineage.n_groups)} lineages` },
          { label: "Run time", value: `${Math.round(r.seconds / 60)} minutes` },
        ]}
      />
      <p className="type-note">
        {legacy
          ? "Labels inferred from MIC are isolates with a measured MIC but no reported call, called against the cohort breakpoint above."
          : "Every isolate here has a reported resistant or susceptible lab call. Genomes whose label was inferred from MIC at the cohort breakpoint are left out; no learner or MIC regression trains or is scored on them."}
      </p>
    </Disclosure>
  );
};

/** The inspector's one column, shared by the header and the body. */
const COLUMN = "mx-auto w-full max-w-[1200px] px-4 sm:px-8";

const Body = ({ r, row, current, research }: { r: DatasetReport; row?: CatalogRow; current: number | undefined; research: Resource<Research> }) => (
  <div className={`${COLUMN} space-y-8 pb-16 pt-6`}>
    <DatasetNotice audit={auditOf(r, current)} />
    <div className="space-y-10">
      <LearnerSection r={r} research={research} />
      <ModelSection r={r} row={row} />
      <MechanismSection r={r} />
      <MicSection r={r} research={research} />
      <LineageSection r={r} />
      <MethodFold r={r} current={current} />
    </div>
  </div>
);

export const DatasetInspector = ({ id, manifest, nav }: { id: string; manifest: Manifest; nav: Navigator }) => {
  const report = useDataset(id);
  const panel = useRef<HTMLDivElement | null>(null);
  const scroller = useRef<HTMLDivElement | null>(null);
  const section = useCurrentSection(scroller, report.status === "ready");
  // Grow from the launching row's centre, read once as the panel opens.
  useLayoutEffect(() => {
    const at = takeLaunch();
    const el = panel.current;
    if (at === null || el === null) return;
    const box = el.getBoundingClientRect();
    el.style.setProperty("--launch-x", `${Math.round(at.x - box.left)}px`);
    el.style.setProperty("--launch-y", `${Math.round(at.y - box.top)}px`);
  }, []);
  // Fetched beside the report; a missing or invalid file reads as no research.
  const research = useResearch();
  const row = manifest.datasets.find((d) => d.id === id);
  const title = row === undefined ? id : `${capitalize(row.antibiotic)} in ${row.species}`;
  // The report's own audit once loaded; the manifest row's copy before.
  const audit = report.status === "ready" ? auditOf(report.data, manifest.evaluation_revision) : row === undefined ? null : auditOf(row, manifest.evaluation_revision);
  return (
    <ModalDialog labelledBy="inspector-title" onClose={nav.closeDataset}>
      <div className="absolute inset-0 bg-paper/60 backdrop-blur-[2px]" aria-hidden="true" onClick={nav.closeDataset} />
      <div
        ref={panel}
        className="inspector-in absolute inset-0 flex flex-col overflow-hidden border border-hair bg-paper-pure text-ink shadow-[var(--surface-shadow)] sm:inset-3 sm:rounded-xl lg:inset-x-6 lg:inset-y-4"
      >
        <header className="shrink-0 border-b border-hair pb-2.5 pt-3">
          <div className={`${COLUMN} flex flex-wrap items-start justify-between gap-x-6 gap-y-2`}>
          <div className="min-w-0">
            <h2 id="inspector-title" className="text-headline font-semibold tracking-tight text-ink">
              {title}
            </h2>
            {row !== undefined ? (
              <p className="mt-0.5 flex flex-wrap items-center gap-x-2 gap-y-1 text-note text-muted-text">
                {audit !== null ? <AuditChip audit={audit} /> : null}
                <span>
                  {row.drug_class} · {int(row.n)} isolates, {int(row.n_resistant)} resistant · {plural(row.lineages, "lineage")}
                  {micLabelWords(row)}
                </span>
              </p>
            ) : null}
          </div>
          <div className="flex items-center gap-2">
            {row !== undefined ? (
              <GhostButton onClick={() => nav.go({ view: "species", species: row.species_slug })}>
                {speciesShort(row.species)} structure
              </GhostButton>
            ) : null}
            <GhostButton onClick={nav.closeDataset}>Close</GhostButton>
          </div>
          </div>
          {report.status === "ready" ? (
            <nav aria-label="Sections" className={`${COLUMN} mt-2.5`}>
              <div className="scroll-thin -mx-1 flex gap-1.5 overflow-x-auto px-1 py-0.5">
                {JUMPS.map((j) => (
                  <button
                    key={j.id}
                    type="button"
                    aria-current={section === j.id ? "true" : undefined}
                    onClick={() =>
                      document.getElementById(j.id)?.scrollIntoView({ behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth", block: "start" })
                    }
                    className={`ui-control-target inline-flex h-7 shrink-0 items-center whitespace-nowrap rounded-full border px-3 text-control transition-[color,background-color,border-color,transform] duration-fast ease-out active:scale-[0.98] ${FOCUS_RING} ${
                      section === j.id ? "border-transparent bg-surface-control font-semibold text-ink" : "border-hair font-medium text-ink hover:border-blue/50 hover:text-blue"
                    }`}
                  >
                    {j.label}
                  </button>
                ))}
              </div>
            </nav>
          ) : null}
        </header>
        {/* The dialog opens with focus on its body, so arrow keys scroll it
            and no header control opens its tip on arrival. */}
        <div
          ref={(el) => {
            scroller.current = el;
            el?.setAttribute("autofocus", "");
          }}
          tabIndex={-1}
          data-modal-scroll-owner="true"
          className="min-h-0 flex-1 overflow-y-auto overscroll-contain focus:outline-none"
        >
          {report.status === "ready" ? (
            <Body r={report.data} row={row} current={manifest.evaluation_revision} research={research} />
          ) : (
            <div className={`${COLUMN} py-8`}>
              {report.status === "loading" ? (
                <PanelState kind="loading" title="Loading the dataset report" />
              ) : (
                <PanelState kind="error" title="This dataset report could not be read" message={loadErrorText(report.error)} />
              )}
            </div>
          )}
        </div>
      </div>
    </ModalDialog>
  );
};
