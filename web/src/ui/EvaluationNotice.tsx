// How the evaluation audit shows. A caveat that limits the numbers is a
// chip, not a sentence: one chip in every view's status line and in the
// inspector header, with its reading and scope on hover. Only a provisional
// audit adds a band with its reasons where metrics are read (Overview and
// the dataset inspector). A passing audit shows its scope, internal
// cross-validation, and is never worded as clinical validation.

import type { ReactNode } from "react";
import type { EvaluationAudit } from "../data/schema";
import type { ManifestAudit } from "../data/derive/evaluation";
import { int } from "../data/format";
import { HoverTip, type TipLine } from "./HoverTip";
import { Lozenge } from "./primitives";

const READING = "The numbers describe agreement with lab calls on these sampled isolates; they do not validate a diagnostic test.";
const NO_SCOPE = "Internal cross-validation against reported lab calls; not external or clinical validation.";

/** The report chip's hover: the reading, then each scope or reason. */
export const reportTip = (audit: ManifestAudit): TipLine[] =>
  audit.invalid > 0
    ? [
        { term: "Provisional metrics", text: `${int(audit.invalid)} of ${int(audit.total)} dataset evaluations do not pass the evaluation audit. ${READING}` },
        ...audit.reasons.map((r) => ({ source: r })),
      ]
    : [
        { term: "Internal cross-validation only", text: `Every dataset evaluation passes the internal evaluation audit. ${READING}` },
        ...(audit.scopes.length > 0 ? audit.scopes : [NO_SCOPE]).map((s) => ({ source: s })),
      ];

/** One dataset's chip hover. */
export const datasetTip = (audit: EvaluationAudit): TipLine[] =>
  audit.valid
    ? [
        { term: "Internal cross-validation only", text: `This dataset's evaluation passes the internal audit (revision ${audit.revision}). ${READING}` },
        { source: audit.scope !== undefined && audit.scope !== "" ? audit.scope : NO_SCOPE },
      ]
    : [
        {
          term: audit.revision === 0 ? "Unaudited, provisional" : "Provisional",
          text: audit.revision === 0 ? "This report predates the evaluation audit, so its metrics are unaudited." : `This dataset's evaluation (audit revision ${audit.revision}) does not pass the current evaluation audit.`,
        },
        ...audit.reasons.map((r) => ({ source: r })),
      ];

export const ProvisionalChip = ({ audit }: { audit: ManifestAudit }) =>
  audit.total === 0 ? null : (
    <HoverTip lines={reportTip(audit)}>
      {audit.invalid > 0 ? (
        <Lozenge tone="coral" fill="solid">
          Provisional metrics
        </Lozenge>
      ) : (
        <Lozenge tone="ink">Internal cross-validation only</Lozenge>
      )}
    </HoverTip>
  );

/** The audit state of one dataset, as a header chip. */
export const AuditChip = ({ audit }: { audit: EvaluationAudit }) => (
  <HoverTip lines={datasetTip(audit)}>
    {audit.valid ? (
      <Lozenge tone="ink">Internal cross-validation only</Lozenge>
    ) : (
      <Lozenge tone="coral" fill="solid">
        {audit.revision === 0 ? "Unaudited, provisional" : "Provisional"}
      </Lozenge>
    )}
  </HoverTip>
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
        <Lozenge tone="coral" fill="solid">
          Provisional
        </Lozenge>{" "}
        {audit.invalid === audit.total
          ? "No dataset evaluation in this report passes the evaluation audit."
          : `${int(audit.invalid)} of ${int(audit.total)} dataset evaluations do not pass the evaluation audit; the rest are internal cross-validation only.`}{" "}
        {READING}
      </p>
      <Lines items={audit.reasons} />
    </Band>
  );

export const DatasetNotice = ({ audit }: { audit: EvaluationAudit }) =>
  audit.valid ? null : (
    <Band>
      <p className="prose-measure text-body text-ink">
        <Lozenge tone="coral" fill="solid">
          Provisional
        </Lozenge>{" "}
        {audit.revision === 0
          ? "This report predates the evaluation audit, so its metrics are unaudited."
          : `This dataset's evaluation (audit revision ${audit.revision}) does not pass the current evaluation audit.`}{" "}
        {READING}
      </p>
      <Lines items={audit.reasons} />
    </Band>
  );
