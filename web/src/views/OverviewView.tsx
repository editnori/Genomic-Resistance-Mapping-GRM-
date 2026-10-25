// What GRM found: the five shared findings as figure blocks in one centered
// column, each with its Data and method. The panel block also holds every
// panel result; /tests lands on it through the #panels anchor.

import { useEffect } from "react";
import type { Manifest } from "../data/schema";
import { plural } from "../data/format";
import { manifestAudit } from "../data/derive/evaluation";
import { ReportNotice } from "../ui/EvaluationNotice";
import type { Navigator } from "../app/router";
import { FIGURE_SHELL, FigureBlock, PanelState } from "../ui/primitives";
import { FINDINGS } from "./findings";

/** The anchor of the panel block, where /tests redirects. */
export const PANELS_ANCHOR = "panels";

export const OverviewView = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const audit = manifestAudit(manifest);
  const curated = manifest.species.reduce((a, s) => a + s.datasets_curated, 0);
  // Read once per mount: a /tests link opens the panel evidence and scrolls
  // to it; the view's own scroll owner does not follow the hash by itself.
  const toPanels = typeof window !== "undefined" && window.location.hash === `#${PANELS_ANCHOR}`;
  useEffect(() => {
    if (toPanels) document.getElementById(PANELS_ANCHOR)?.scrollIntoView({ block: "start" });
  }, [toPanels]);

  return (
    <div className={`${FIGURE_SHELL} pb-16 pt-8`}>
      <h1 className="text-display font-semibold text-ink">Overview</h1>
      {audit.invalid > 0 ? (
        <div className="mt-6">
          <ReportNotice audit={audit} />
        </div>
      ) : null}
      {manifest.datasets.length === 0 ? (
        <div className="mt-8">
          <PanelState kind="empty" title="No dataset has a report yet" message={`${plural(curated, "dataset")} passed curation; reports appear as \`grm learn\` finishes each species.`} />
        </div>
      ) : (
        <div className="mt-8 space-y-12 sm:space-y-16">
          {FINDINGS.map((build) => {
            const f = build(manifest, { page: true, onOpen: nav.openDataset });
            return (
              <FigureBlock
                key={f.id}
                id={f.id}
                title={f.title}
                result={f.result}
                figure={f.figure}
                caption={f.caption}
                details={f.details}
                detailsOpen={f.id === PANELS_ANCHOR && toPanels}
              />
            );
          })}
        </div>
      )}
    </div>
  );
};
