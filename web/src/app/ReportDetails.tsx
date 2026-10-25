import { useState } from "react";
import { X } from "lucide-react";
import type { Manifest } from "../data/schema";
import { dateText, int } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { ModalDialog } from "../ui/ModalDialog";
import { FOCUS_RING } from "../ui/primitives";

export const ReportDetails = ({ manifest }: { manifest: Manifest }) => {
  const [open, setOpen] = useState(false);
  const audit = manifestAudit(manifest);
  return <>
    <button type="button" onClick={() => setOpen(true)} aria-haspopup="dialog"
      className={`min-h-11 rounded-md px-3 text-control ${FOCUS_RING} ${audit.invalid > 0 ? "bg-[var(--status-miss-bg)] text-[var(--status-miss-fg)]" : "text-muted-text hover:bg-surface-control hover:text-ink"}`}>
      {audit.invalid > 0 ? "Provisional report" : "Report"}
    </button>
    {open ? <ModalDialog label="Report details" onClose={() => setOpen(false)}>
      <div className="absolute inset-0 bg-paper/60" aria-hidden="true" onClick={() => setOpen(false)} />
      <div data-modal-scroll-owner="true" className="absolute inset-x-4 top-1/2 mx-auto max-h-[calc(100dvh-2rem)] max-w-lg -translate-y-1/2 overflow-y-auto rounded-xl border border-hair bg-paper-pure p-6 text-ink shadow-xl">
      <div className="flex items-center justify-between gap-4">
        <h2 className="text-headline font-semibold">Report details</h2>
        <button type="button" onClick={() => setOpen(false)} aria-label="Close report details" className={`inline-flex h-11 w-11 items-center justify-center rounded-md hover:bg-surface-control ${FOCUS_RING}`}><X size={20} aria-hidden="true" /></button>
      </div>
      <dl className="mt-4 grid grid-cols-[auto_1fr] gap-x-6 gap-y-3 text-note">
        <dt className="text-muted-text">Built</dt><dd>{dateText(manifest.generated_at)}</dd>
        <dt className="text-muted-text">Datasets</dt><dd>{int(manifest.datasets.length)}</dd>
        <dt className="text-muted-text">Lab records</dt><dd>{int(manifest.amr_records)} from BV-BRC</dd>
        <dt className="text-muted-text">Catalog</dt><dd>AMRFinderPlus {manifest.amrfinder_version || "version not recorded"}</dd>
        <dt className="text-muted-text">Evaluation</dt><dd>{manifest.evaluation_revision === undefined ? "Revision not recorded" : `Revision ${manifest.evaluation_revision}`}</dd>
      </dl>
      <p className="mt-5 text-note">{audit.invalid > 0 ? `${int(audit.invalid)} of ${int(audit.total)} evaluations are provisional.` : "Internal cross-validation against reported laboratory calls."}</p>
      <p className="mt-2 text-note text-muted-text">These comparisons have not been externally or clinically validated.</p>
      {audit.reasons.length > 0 ? <ul className="mt-3 list-disc space-y-2 pl-4 text-note">{audit.reasons.map(reason => <li key={reason}>{reason}</li>)}</ul> : null}
      </div>
    </ModalDialog> : null}
  </>;
};
