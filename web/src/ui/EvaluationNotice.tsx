// How the evaluation audit shows. Its state is plain words in each view's
// status line and the inspector header: AUDIT_SCOPE when every evaluation
// passes, bold "Provisional metrics" when any does not. In the status line
// the words link to the Overview's limits, where the scope is spelled out;
// no hover holds anything. Only a provisional audit adds a band with its
// reasons where metrics are read (Overview and the dataset inspector). A
// passing audit is never worded as clinical validation.

import type { ReactNode } from "react";
import type { EvaluationAudit } from "../data/schema";
import type { ManifestAudit } from "../data/derive/evaluation";
import { int } from "../data/format";
import type { TipLine } from "./HoverTip";
import { FOCUS_RING } from "./primitives";

/** What every score measures, in plain words: the one reading each audit
 *  line, notice, and status band repeats. */
export const AUDIT_READING =
  "These numbers measure how often the models agreed with the lab's resistant or susceptible calls on these sampled isolates. They do not show that a diagnostic test works.";
const NO_SCOPE = "Internal cross-validation against reported lab calls; not external or clinical validation.";

/** The scope every passing score shares, in plain words: internal
 *  cross-validation, with training and testing inside BV-BRC. The status
 *  line, the inspector header and Present's title slide print it. */
export const AUDIT_SCOPE = "Tested within one public database only";

/** The report audit in lines: the reading, then each scope or reason. */
export const reportTip = (audit: ManifestAudit): TipLine[] =>
  audit.invalid > 0
    ? [
        { term: "Provisional metrics", text: `${int(audit.invalid)} of ${int(audit.total)} dataset evaluations do not pass the evaluation audit. ${AUDIT_READING}` },
        ...audit.reasons.map((r) => ({ source: r })),
      ]
    : [
        { term: AUDIT_SCOPE, text: `Every dataset evaluation passes the internal evaluation audit. ${AUDIT_READING}` },
        ...(audit.scopes.length > 0 ? audit.scopes : [NO_SCOPE]).map((s) => ({ source: s })),
      ];

/** One dataset's audit in lines, for a view that prints its scope. */
export const datasetTip = (audit: EvaluationAudit): TipLine[] =>
  audit.valid
    ? [
        { term: AUDIT_SCOPE, text: `This dataset's evaluation passes the internal audit (revision ${audit.revision}). ${AUDIT_READING}` },
        { source: audit.scope !== undefined && audit.scope !== "" ? audit.scope : NO_SCOPE },
      ]
    : [
        {
          term: audit.revision === 0 ? "Unaudited, provisional" : "Provisional",
          text: audit.revision === 0 ? "This report predates the evaluation audit, so its metrics are unaudited." : `This dataset's evaluation (audit revision ${audit.revision}) does not pass the current evaluation audit.`,
        },
        ...audit.reasons.map((r) => ({ source: r })),
      ];

/** The provisional word: bold, in the warning tone. */
const Provisional = ({ children }: { children: ReactNode }) => <strong className="font-semibold text-coral-text">{children}</strong>;

const STATUS_LINK = `rounded-sm underline decoration-current/40 underline-offset-2 transition-colors duration-fast ease-out hover:text-blue hover:decoration-current ${FOCUS_RING}`;

/** The report's audit state in words; with `href`, the words link to where
 *  the scope and any reasons are printed. */
export const AuditStatus = ({ audit, href }: { audit: ManifestAudit; href?: string }) => {
  if (audit.total === 0) return null;
  const words = audit.invalid > 0 ? <Provisional>Provisional metrics</Provisional> : AUDIT_SCOPE;
  return href === undefined ? (
    <span data-audit-status="">{words}</span>
  ) : (
    <a data-audit-status="" href={href} className={STATUS_LINK}>
      {words}
    </a>
  );
};

/** The report's audit state as plain words, without a link: for Present,
 *  which stays on its slides. */
export const ProvisionalChip = ({ audit }: { audit: ManifestAudit }) => <AuditStatus audit={audit} />;

/** One dataset's audit state in the inspector header. */
export const AuditChip = ({ audit }: { audit: EvaluationAudit }) => (
  <span data-audit-status="">{audit.valid ? AUDIT_SCOPE : <Provisional>{audit.revision === 0 ? "Unaudited, provisional" : "Provisional"}</Provisional>}</span>
);

const Lines = ({ items }: { items: readonly string[] }) =>
  items.length === 0 ? null : (
    <ul className="type-note mt-1.5 list-disc space-y-0.5 pl-5">
      {items.map((r) => (
        <li key={r}>{r}</li>
      ))}
    </ul>
  );

const Band = ({ children }: { children: ReactNode }) => (
  <div role="note" className="border-y border-hair py-3">
    {children}
  </div>
);

export const ReportNotice = ({ audit }: { audit: ManifestAudit }) =>
  audit.total === 0 || audit.invalid === 0 ? null : (
    <Band>
      <p className="prose-measure text-body text-ink">
        <Provisional>Provisional metrics.</Provisional>{" "}
        {audit.invalid === audit.total
          ? "No dataset evaluation in this report passes the evaluation audit, the check that scores were computed under the current rules, so every number is provisional until the datasets are recomputed."
          : `${int(audit.invalid)} of ${int(audit.total)} dataset evaluations do not pass the evaluation audit, the check that scores were computed under the current rules, so their numbers are provisional until those datasets are recomputed. The rest pass it.`}{" "}
        {AUDIT_READING}
      </p>
      <Lines items={audit.reasons} />
    </Band>
  );

export const DatasetNotice = ({ audit }: { audit: EvaluationAudit }) =>
  audit.valid ? null : (
    <Band>
      <p className="prose-measure text-body text-ink">
        <Provisional>Provisional.</Provisional>{" "}
        {audit.revision === 0
          ? "This report predates the evaluation audit, so its metrics are unaudited."
          : `This dataset's evaluation (audit revision ${audit.revision}) does not pass the current evaluation audit.`}{" "}
        {AUDIT_READING}
      </p>
      <Lines items={audit.reasons} />
    </Band>
  );
