import type { ReactNode } from "react";
import type { Manifest } from "../data/schema";
import { dateText, int } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { ProvisionalChip } from "../ui/EvaluationNotice";

const Sep = () => <span aria-hidden="true">·</span>;

/** The one quiet status line: the evaluation audit's chip, when the report
 *  was built, the catalog version it used, the records it started from,
 *  then view facts. */
export const ReportStatus = ({ manifest, children }: { manifest: Manifest; children?: ReactNode }) => (
  <>
    <ProvisionalChip audit={manifestAudit(manifest)} />
    <span>Report built {dateText(manifest.generated_at)}</span>
    <Sep />
    <span>AMRFinderPlus {manifest.amrfinder_version === "" ? "version not recorded" : manifest.amrfinder_version}</span>
    <Sep />
    <span>{int(manifest.amr_records)} BV-BRC lab records</span>
    {children !== undefined ? (
      <>
        <Sep />
        {children}
      </>
    ) : null}
  </>
);
