// The BV-BRC walk: five interactive steps showing what comes in from BV-BRC,
// what each field is used for, and the choices made on the way to a dataset.

import { useCallback, useMemo, useState } from "react";
import { FLOW_DATA, buildNodeMap, fillTemplate, uiText, useReducedMotion, type FlowData } from "./model";
import {
  ChoiceRowsList,
  ContigsStrip,
  FlowBars,
  FlowCard,
  FlowHeader,
  FlowStats,
  StepStrip,
  WalkControls,
} from "./shared";

export interface BvbrcFlowProps {
  flow?: FlowData;
  startStep?: string;
  calm?: boolean;
}

export const BvbrcFlow = ({ flow = FLOW_DATA, startStep, calm: calmProp }: BvbrcFlowProps) => {
  const reducedMotion = useReducedMotion();
  const calm = calmProp ?? reducedMotion;
  const view = flow.views.bvbrc;
  const ui = view.ui;
  const nodes = useMemo(() => buildNodeMap(flow), [flow]);
  const steps = view.steps;

  const [stepIdx, setStepIdx] = useState<number>(() => Math.max(0, steps.indexOf(startStep ?? view.start)));

  const go = useCallback(
    (nextIdx: number) => {
      setStepIdx(((nextIdx % steps.length) + steps.length) % steps.length);
    },
    [steps.length],
  );

  const activeId = steps[stepIdx] ?? steps[0]!;
  const node = nodes[activeId];
  const decisions = node?.decisions ?? [];
  const toolText = (node?.fields ?? [])
    .filter((f) => f[0] === "Tool")
    .map((f) => `${f[0]}: ${fillTemplate(flow, f[1])}`)
    .join("");

  const stepItems = steps.map((id) => ({
    id,
    label: fillTemplate(flow, nodes[id]?.title ?? id),
  }));

  const hasAboveStats = node !== undefined && (node.card !== undefined || node.draw === "contigs" || node.bars !== undefined);

  return (
    <section
      data-flow-diagram="bvbrc"
      tabIndex={-1}
      onKeyDown={(e) => {
        if (e.key === "ArrowRight") go(stepIdx + 1);
        else if (e.key === "ArrowLeft") go(stepIdx - 1);
      }}
      className="mx-auto w-full max-w-[760px] text-ink outline-none"
    >
      <FlowHeader flow={flow} kicker={view.kicker} title={view.title} sub={view.sub} />

      <StepStrip steps={stepItems} activeIndex={stepIdx} onSelect={go} ariaLabel={fillTemplate(flow, view.kicker)} />

      <WalkControls
        className="mb-3"
        backLabel={uiText(ui, "back")}
        nextLabel={uiText(ui, "next")}
        onBack={() => go(stepIdx - 1)}
        onNext={() => go(stepIdx + 1)}
      />

      {node !== undefined ? (
        <div className="min-h-[340px]">
          <p className="mb-3 max-w-[62ch] text-body text-ink">{fillTemplate(flow, node.line)}</p>

          <div className="overflow-hidden rounded-lg border border-hair bg-paper-pure p-3.5">
            {node.card !== undefined ? <FlowCard flow={flow} node={node} ui={ui} /> : null}
            {node.draw === "contigs" ? <ContigsStrip caption={uiText(ui, "drawing")} /> : null}
            {node.bars !== undefined ? <FlowBars flow={flow} node={node} calm={calm} /> : null}
            {node.stats !== undefined ? (
              <div className={hasAboveStats ? "mt-4" : ""}>
                <FlowStats flow={flow} node={node} />
              </div>
            ) : null}
          </div>

          {toolText !== "" ? <p className="mt-2 font-mono text-tiny text-muted-text">{toolText}</p> : null}

          <div className="mt-4">
            <p className="mb-1 text-tiny font-semibold uppercase tracking-[0.06em] text-muted-text">
              {decisions.length > 0 ? uiText(ui, "choices") : uiText(ui, "none")}
            </p>
            {decisions.length > 0 ? (
              <ChoiceRowsList flow={flow} nodes={nodes} ids={decisions} recordLabel={uiText(ui, "record")} />
            ) : null}
          </div>
        </div>
      ) : null}
    </section>
  );
};
