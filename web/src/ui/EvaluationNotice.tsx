// How a failed or provisional evaluation audit shows: one band with its
// reasons where metrics are read (Overview and the dataset inspector), and a
// chip in the inspector header. A passing audit shows nothing here; the
// report details dialog in the top bar states its scope once.

import type { ReactNode } from "react";
import type { EvaluationAudit } from "../data/schema";
import type { ManifestAudit } from "../data/derive/evaluation";
import { int } from "../data/format";
import { Lozenge } from "./primitives";

const READING = "The numbers describe agreement with lab calls on these sampled isolates; they do not validate a diagnostic test.";

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
        <Lozenge tone="provisional">Provisional</Lozenge>{" "}
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
        <Lozenge tone="provisional">Provisional</Lozenge>{" "}
        {audit.revision === 0
          ? "This report predates the evaluation audit, so its metrics are unaudited."
          : `This dataset's evaluation (audit revision ${audit.revision}) does not pass the current evaluation audit.`}{" "}
        {READING}
      </p>
      <Lines items={audit.reasons} />
    </Band>
  );

/** The audit state of one dataset as a header chip, only when provisional. */
export const AuditChip = ({ audit }: { audit: EvaluationAudit }) =>
  audit.valid ? null : <Lozenge tone="provisional">{audit.revision === 0 ? "Unaudited, provisional" : "Provisional"}</Lozenge>;
