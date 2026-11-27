// Preview route (/flow and /flow/<slug>) for the five React flow diagrams.
// Uses the site's own ViewHeader, Segmented control, and theme tokens. Not
// linked from TopBar until browser review and tree sign-off finish.

import { Segmented, VIEW_SHELL } from "../ui/primitives";
import { ViewHeader } from "../ui/ViewHeader";
import {
  BvbrcFlow,
  DecisionTreeFlow,
  FLOW_DATA,
  FLOW_SLUGS,
  GenomePathFlow,
  GenomeTableFlow,
  PieceCountFlow,
  fillTemplate,
  type FlowData,
  type FlowSlug,
} from "../flow";

export type FlowPreviewTarget = "all" | FlowSlug;

const renderDiagram = (slug: FlowSlug, flow: FlowData) => {
  switch (slug) {
    case "genome-table":
      return <GenomeTableFlow flow={flow} />;
    case "genome-path":
      return <GenomePathFlow flow={flow} />;
    case "bvbrc":
      return <BvbrcFlow flow={flow} />;
    case "piece-count":
      return <PieceCountFlow flow={flow} />;
    case "decision-tree":
      return <DecisionTreeFlow flow={flow} />;
  }
};

export const FlowPreviewView = ({
  target,
  onSelect,
  flow = FLOW_DATA,
}: {
  target: FlowPreviewTarget;
  onSelect(next: FlowPreviewTarget): void;
  flow?: FlowData;
}) => {
  const options: { id: FlowPreviewTarget; label: string }[] = [
    { id: "all", label: "All five" },
    ...FLOW_SLUGS.map(({ slug, viewKey }) => ({
      id: slug as FlowPreviewTarget,
      label: fillTemplate(flow, flow.views[viewKey].kicker),
    })),
  ];

  const shownSlugs: readonly FlowSlug[] = target === "all" ? FLOW_SLUGS.map((entry) => entry.slug) : [target];

  return (
    <div className={`${VIEW_SHELL} py-7`}>
      <ViewHeader
        title="Flow diagrams preview"
        purpose="The five genome-path drawings rendered in React from flow.json in the site's light and dark themes."
        actions={<Segmented ariaLabel="Flow diagram" options={options} value={target} onSelect={onSelect} />}
        status={
          <>
            <span>Preview route</span>
            <span aria-hidden="true">&middot;</span>
            <span>Pinned commit {flow.pinned.commit.slice(0, 7)}</span>
            <span aria-hidden="true">&middot;</span>
            <span>{Object.keys(flow.counts).length} sourced counts</span>
          </>
        }
      />

      <div className="mt-7 space-y-10">
        {shownSlugs.map((slug, idx) => (
          <div key={slug} className={idx > 0 ? "border-t border-hair pt-8" : ""}>
            {renderDiagram(slug, flow)}
          </div>
        ))}
      </div>
    </div>
  );
};
