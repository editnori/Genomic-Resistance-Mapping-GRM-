// The decision tree: the plan's choices along the genome's path, one row a
// choice (what we took, what we tried against it, what we found) with a
// single expandable detail panel and keyboard up/down stepping.

import { useCallback, useMemo, useState } from "react";
import { FLOW_DATA, buildNodeMap, fillTemplate, uiText, type FlowData } from "./model";
import { ChoiceDetail, ChoiceRowItem, FlowHeader, WalkControls } from "./shared";

export interface DecisionTreeFlowProps {
  flow?: FlowData;
  startChoice?: string;
}

export const DecisionTreeFlow = ({ flow = FLOW_DATA, startChoice }: DecisionTreeFlowProps) => {
  const view = flow.views.decisions;
  const ui = view.ui;
  const nodes = useMemo(() => buildNodeMap(flow), [flow]);
  const order = useMemo(() => view.trunk.flatMap(([, rows]) => rows), [view.trunk]);

  const initialId = startChoice ?? view.start;
  const [openId, setOpenId] = useState<string>(order.includes(initialId) ? initialId : (order[0] ?? ""));

  const step = useCallback(
    (delta: number) => {
      const idx = Math.max(0, order.indexOf(openId));
      const next = order[(idx + delta + order.length) % order.length];
      if (next !== undefined) setOpenId(next);
    },
    [openId, order],
  );

  return (
    <section
      data-flow-diagram="decision-tree"
      tabIndex={-1}
      onKeyDown={(e) => {
        if (e.key === "ArrowDown") {
          e.preventDefault();
          step(1);
        } else if (e.key === "ArrowUp") {
          e.preventDefault();
          step(-1);
        }
      }}
      className="mx-auto w-full max-w-[760px] text-ink outline-none"
    >
      <FlowHeader flow={flow} kicker={view.kicker} title={view.title} sub={view.sub} />

      <div className="mb-3 flex flex-wrap items-center gap-x-4 gap-y-2 text-tiny text-muted-text">
        <span className="inline-flex items-center">
          <span aria-hidden="true" className="mr-1.5 inline-block h-3.5 w-[3px] bg-coral-text" />
          {uiText(ui, "bare")}
        </span>
        <WalkControls
          className="ml-auto"
          backLabel={uiText(ui, "back")}
          nextLabel={uiText(ui, "next")}
          onBack={() => step(-1)}
          onNext={() => step(1)}
        />
      </div>

      <div>
        {view.trunk.map(([stepId, rowIds], sectionIdx) => {
          const stepNode = nodes[stepId];
          if (stepNode === undefined) return null;
          const isFirst = sectionIdx === 0;
          const isLast = sectionIdx === view.trunk.length - 1;
          return (
            <section key={stepId} className="relative pl-[22px]">
              {/* Vertical trunk spine */}
              <span
                aria-hidden="true"
                className="pointer-events-none absolute left-[6px] w-[2px] bg-ink"
                style={{
                  top: isFirst ? "12px" : 0,
                  bottom: isLast ? "auto" : 0,
                  height: isLast ? "12px" : undefined,
                }}
              />
              <div className="relative flex flex-wrap items-baseline gap-x-2.5 gap-y-0.5 pb-1 pt-1.5">
                <span
                  aria-hidden="true"
                  className="pointer-events-none absolute -left-[22px] top-[9px] h-3.5 w-3.5 rounded-full bg-ink shadow-[0_0_0_3px_rgb(var(--paper))]"
                />
                <span className="text-body font-semibold tracking-tight text-ink">{fillTemplate(flow, stepNode.title)}</span>
                <span className="text-tiny text-muted-text">{fillTemplate(flow, stepNode.line)}</span>
              </div>

              <div className="relative mb-1">
                {rowIds.map((rowId) => {
                  const rowNode = nodes[rowId];
                  if (rowNode === undefined) return null;
                  const expanded = openId === rowId;
                  return (
                    <div key={rowId}>
                      <ChoiceRowItem
                        flow={flow}
                        node={rowNode}
                        expanded={expanded}
                        onSelect={() => setOpenId(rowId)}
                        branchTick
                      />
                      {expanded ? (
                        <ChoiceDetail flow={flow} nodes={nodes} id={rowId} recordLabel={uiText(ui, "record")} />
                      ) : null}
                    </div>
                  );
                })}
              </div>
            </section>
          );
        })}
      </div>
    </section>
  );
};
