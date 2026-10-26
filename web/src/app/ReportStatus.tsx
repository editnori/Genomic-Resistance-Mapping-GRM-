import type { ReactNode } from "react";
import type { Manifest } from "../data/schema";
import { dateText } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { AuditStatus } from "../ui/EvaluationNotice";
import { LIMITS_HREF } from "./router";

/** Facts in a wrapping line, each carrying the separator after it, so a
 *  wrapped line never starts with one. */
export const StatusItems = ({ items }: { items: readonly ReactNode[] }) => (
  <>
    {items.map((item, i) => (
      <span key={i} className="inline-flex min-w-0 items-baseline gap-x-2">
        {item}
        {i < items.length - 1 ? <span aria-hidden="true">·</span> : null}
      </span>
    ))}
  </>
);

/** The one quiet status line, in plain words: when the report was built,
 *  the catalog version it used, and the audit's scope ("Internal
 *  cross-validation only", or bold "Provisional metrics"), which links to
 *  the Overview's limits; then view facts. */
export const ReportStatus = ({ manifest, children }: { manifest: Manifest; children?: ReactNode }) => {
  const audit = manifestAudit(manifest);
  const items: ReactNode[] = [
    `Report built ${dateText(manifest.generated_at)}`,
    `AMRFinderPlus ${manifest.amrfinder_version === "" ? "version not recorded" : manifest.amrfinder_version}`,
  ];
  if (audit.total > 0) items.push(<AuditStatus audit={audit} href={LIMITS_HREF} />);
  if (children !== undefined) items.push(children);
  return <StatusItems items={items} />;
};
