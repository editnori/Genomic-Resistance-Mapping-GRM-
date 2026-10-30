// What GRM found: a five-step primer for a reader new to the terms, then the
// report in six numbered sections. The data behind the reports, then the
// four findings Present shares, each with its lead, a figure, a caption that
// says what a mark is, and short notes on coverage in view; last, the limits
// every number shares. The panel section also holds every panel result;
// /tests lands on it through the #panels anchor.

import { useEffect, useId, useState, type ReactNode } from "react";
import { ChevronDown } from "lucide-react";
import type { Manifest } from "../data/schema";
import { int, plural } from "../data/format";
import type { Navigator } from "../app/router";
import { ViewHeader } from "../ui/ViewHeader";
import { FOCUS_RING, PanelState, VIEW_SHELL } from "../ui/primitives";
import { PrimerSection } from "./primer";
import { AuditReasons, FINDINGS, inApp, LINK, methodsHref, openMethodsAt, ReportLine, SectionHeading, selectionCap, showSection, STATUS_LINK, type Finding } from "./findings";

/** The anchor of the panel section, where /tests redirects. */
export const PANELS_ANCHOR = "panels";
/** The anchor of the closing section, where the status line's scope links. */
export const LIMITS_ANCHOR = "limits";

/** Supporting evidence under one explicit control, open on arrival from a
 *  link that names it. */
const More = ({ label, open: initial, children }: { label: string; open: boolean; children: ReactNode }) => {
  const [open, setOpen] = useState(initial);
  const id = useId();
  return (
    <div className="mt-6">
      <button
        type="button"
        aria-expanded={open}
        aria-controls={open ? id : undefined}
        onClick={() => setOpen((v) => !v)}
        className={`inline-flex min-h-11 items-center gap-1.5 rounded-sm text-control font-medium text-blue transition-colors duration-fast ease-out hover:text-ink ${FOCUS_RING}`}
      >
        {open ? `Hide ${label.toLowerCase()}` : `Show ${label.toLowerCase()}`}
        <ChevronDown size={14} strokeWidth={1.75} aria-hidden="true" className={`transition-transform duration-fast ease-out ${open ? "rotate-180" : ""}`} />
      </button>
      {open ? (
        <div id={id} className="motion-disclosure-in mt-4 text-note text-ink">
          {children}
        </div>
      ) : null}
    </div>
  );
};

/** One finding as a report section: heading, lead, then the figure with its
 *  caption and notes as one unit, centered when it is a plot with a legend
 *  and across the width when it is a strip. */
const ReportSection = ({ n, finding, open }: { n: number; finding: Finding; open: boolean }) => {
  const figure = typeof finding.figure === "function" ? finding.figure({}) : finding.figure;
  return (
    <section id={finding.id} aria-labelledby={`${finding.id}-title`} className="scroll-mt-6 border-t border-hair pt-7">
      <SectionHeading id={finding.id} n={n}>
        {finding.title}
      </SectionHeading>
      <p className="prose-measure mt-2 text-body text-ink">{finding.result}</p>
      <div className={finding.span === "unit" ? "mx-auto mt-6 w-full max-w-[60rem]" : "mt-6"}>
        {figure === null || figure === undefined ? null : <div className="min-w-0">{figure}</div>}
        <p className={`prose-measure text-note text-muted-text ${figure === null || figure === undefined ? "" : "mt-3"}`}>{finding.caption}</p>
        {finding.note}
      </div>
      {finding.more === undefined ? null : (
        <More label={finding.more.label} open={open}>
          {finding.more.body}
        </More>
      )}
    </section>
  );
};

/** The limits every number shares, and where the audit lives. */
const Limits = ({ n, manifest, nav }: { n: number; manifest: Manifest; nav: Navigator }) => {
  const { cap, capped } = selectionCap(manifest);
  return (
    <section id={LIMITS_ANCHOR} aria-labelledby={`${LIMITS_ANCHOR}-title`} className="scroll-mt-6 border-t border-hair pt-7">
      <SectionHeading id={LIMITS_ANCHOR} n={n}>
        Limits
      </SectionHeading>
      <ul className="prose-measure mt-3 list-disc space-y-2 pl-5 text-body text-ink">
        <li>Every number here is internal cross-validation on public BV-BRC genomes, scored against the resistant and susceptible calls BV-BRC reports.</li>
        <li>
          {capped
            ? `The genomes are a selection: at most ${int(cap)} per species in this run, and fewer where QC or FASTA availability reduced a cohort.`
            : "The genomes are a selection, reduced where QC or FASTA availability limited a cohort."}
        </li>
        <li>No external, prospective, or clinical validation has been done. The numbers describe agreement with reported lab calls on these isolates, not how a diagnostic test would perform.</li>
      </ul>
      <p className="mt-3 text-control">
        <a href={methodsHref("m-audit")} onClick={inApp(() => openMethodsAt(nav, "m-audit"))} className={`inline-flex min-h-11 items-center gap-1.5 font-medium ${LINK}`}>
          See the evaluation audit in Methods
          <span aria-hidden="true">→</span>
        </a>
      </p>
    </section>
  );
};

export const OverviewView = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const curated = manifest.species.reduce((a, s) => a + s.datasets_curated, 0);
  // Read once per mount: a link to one section, such as /tests landing on
  // the panel evidence, scrolls the view's own scroll owner there.
  const [hash] = useState(() => (typeof window === "undefined" ? "" : decodeURIComponent(window.location.hash.slice(1))));
  useEffect(() => {
    if (hash !== "") document.getElementById(hash)?.scrollIntoView?.({ block: "start" });
  }, [hash]);
  const findings = FINDINGS.map((build) => build(manifest, { page: true, nav }));

  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader
        title="What GRM found"
        status={
          <ReportLine
            manifest={manifest}
            scope={
              <a
                href={`#${LIMITS_ANCHOR}`}
                onClick={(e) => {
                  e.preventDefault();
                  showSection(LIMITS_ANCHOR, true);
                }}
                className={STATUS_LINK}
              >
                Internal cross-validation only
              </a>
            }
          />
        }
      />
      <AuditReasons manifest={manifest} />
      {manifest.datasets.length === 0 ? (
        <div className="mt-8">
          <PanelState kind="empty" title="No dataset has a report yet" message={`${plural(curated, "dataset")} passed curation; reports appear as \`grm learn\` finishes each species.`} />
        </div>
      ) : (
        <div className="mt-10 space-y-12">
          <PrimerSection manifest={manifest} nav={nav} firstFinding={findings[0]!.id} />
          {findings.map((f, i) => (
            <ReportSection key={f.id} n={i + 1} finding={f} open={f.id === PANELS_ANCHOR && hash === PANELS_ANCHOR} />
          ))}
          <Limits n={findings.length + 1} manifest={manifest} nav={nav} />
        </div>
      )}
    </div>
  );
};
