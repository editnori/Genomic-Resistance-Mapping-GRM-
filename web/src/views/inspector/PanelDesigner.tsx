// Design a marker test in the browser: pick model features (presence from
// their bits) and catalog gene groups (from each isolate's groups), choose
// any or all, and read sensitivity, specificity, VME, and ME at once.

import { useState } from "react";
import type { DatasetReport } from "../../data/schema";
import { decodeBits } from "../../data/bits";
import { int, pct } from "../../data/format";
import { agreement, markerCalls, BARS, type ReferenceBar } from "../../data/metrics";
import { featureIndex, featureLabel, featureWhere, isGenePanel, modelFeatureOrder, truthOf } from "../../data/derive/dataset";
import { FOCUS_RING, GhostButton, Segmented, verdictWords } from "../../ui/primitives";

const SHORT_LIST = 10;

interface Marker {
  key: string;
  label: string;
  detail: string | null;
  r: number;
  s: number;
}

const Stat = ({ label, value, bar }: { label: string; value: number; bar?: ReferenceBar }) => (
  <div className="min-w-0">
    <div className="text-control font-medium text-muted-text">{label}</div>
    <div className="type-stat mt-1">{pct(value)}</div>
    {bar === undefined ? null : <p className="text-note text-muted-text">{verdictWords(value, bar)}</p>}
  </div>
);

const MarkerList = ({ title, markers, chosen, onToggle }: { title: string; markers: Marker[]; chosen: ReadonlySet<string>; onToggle(key: string): void }) => {
  const [all, setAll] = useState(false);
  // Collapsed, the list keeps its first rows plus every chosen marker, in
  // list order, so a choice is never hidden.
  const shown = all ? markers : markers.filter((m, i) => i < SHORT_LIST || chosen.has(m.key));
  return (
    <fieldset className="min-w-0">
      <legend className="text-control font-semibold text-ink">{title}</legend>
      {markers.length === 0 ? <p className="type-note mt-2">None in this dataset.</p> : null}
      <ul className="mt-2 divide-y divide-hair border-y border-hair">
        {shown.map((m) => (
          <li key={m.key}>
            <label className={`flex cursor-pointer items-start gap-3 rounded-md px-1 py-1.5 transition-colors duration-fast ease-out hover:bg-surface-subtle has-[:focus-visible]:outline has-[:focus-visible]:outline-2 has-[:focus-visible]:outline-blue`}>
              <input type="checkbox" checked={chosen.has(m.key)} onChange={() => onToggle(m.key)} className="mt-1 h-4 w-4 shrink-0 accent-[rgb(var(--blue))]" />
              <span className="min-w-0 flex-1">
                <span className="block text-control font-medium text-ink">{m.label}</span>
                {m.detail !== null ? <span className="block text-tiny text-muted-text">{m.detail}</span> : null}
              </span>
              <span className="shrink-0 text-right text-tiny tabular-nums text-muted-text">
                {int(m.r)} R · {int(m.s)} S
              </span>
            </label>
          </li>
        ))}
      </ul>
      {markers.length > shown.length || all ? (
        <div className="mt-2">
          <GhostButton onClick={() => setAll((v) => !v)}>{all ? `Show ${SHORT_LIST}` : `Show all ${markers.length}`}</GhostButton>
        </div>
      ) : null}
    </fieldset>
  );
};

export const PanelDesigner = ({ r }: { r: DatasetReport }) => {
  const n = r.samples.length;
  const truth = truthOf(r);
  const cards = featureIndex(r);
  const featureMarkers: Marker[] = modelFeatureOrder(r).flatMap((f) => {
    const c = cards.get(f);
    return c === undefined ? [] : [{ key: `f${f}`, label: featureLabel(c, f), detail: featureWhere(c, f), r: c.n_resistant, s: c.n_susceptible }];
  });
  const groupMarkers: Marker[] = r.catalog_groups
    .map((g, i) => {
      const carriers = r.samples.map((s) => s.groups.includes(i));
      return { key: `g${i}`, label: g, detail: null, r: carriers.filter((c, k) => c && truth[k]).length, s: carriers.filter((c, k) => c && !truth[k]).length };
    })
    .sort((a, b) => b.r - a.r || a.s - b.s);

  // Arrival: the first gene panel's detected targets, so the designer starts
  // on a real assay's targets (probe assays have no catalog groups to pick);
  // without one, the set covering machine's present rules.
  const genePanel = r.panels.find((p) => isGenePanel(r, p));
  const [chosen, setChosen] = useState<Set<string>>(() => {
    const panel = genePanel;
    if (panel !== undefined) {
      const detected = new Set(panel.detected.map(([g]) => g));
      return new Set(r.catalog_groups.flatMap((g, i) => (detected.has(g) ? [`g${i}`] : [])));
    }
    return new Set(r.models.scm.rules.filter((rule) => rule.presence && cards.has(rule.feature)).map((rule) => `f${rule.feature}`));
  });
  const [rule, setRule] = useState<"any" | "all">("any");

  const toggle = (key: string) =>
    setChosen((cur) => {
      const next = new Set(cur);
      if (next.has(key)) next.delete(key);
      else next.add(key);
      return next;
    });

  const presence = [...chosen].map((key) => {
    if (key.startsWith("f")) {
      const c = cards.get(Number(key.slice(1)));
      return c === undefined ? new Array<boolean>(n).fill(false) : decodeBits(c.bits, n);
    }
    const g = Number(key.slice(1));
    return r.samples.map((s) => s.groups.includes(g));
  });
  const a = agreement(truth, markerCalls(presence, n, rule));
  return (
    <div className="mt-8">
      <div className="flex flex-wrap items-start justify-between gap-3">
        <h3 className="text-body font-semibold text-ink">Marker designer</h3>
        <Segmented<"any" | "all">
          ariaLabel="Call resistant when"
          options={[{ id: "any", label: "Any marker present" }, { id: "all", label: "All markers present" }]}
          value={rule}
          onSelect={setRule}
        />
      </div>
      <p className="prose-measure my-4 text-note text-ink">
        {chosen.size === 0
          ? "Choose one or more markers to see what a test built from them would call."
          : `This ${markerWord(chosen.size)} rule calls ${int(a.tp)} of ${int(a.tp + a.fn)} resistant isolates resistant, and ${int(a.fp)} of ${int(a.fp + a.tn)} susceptible isolates.`}
      </p>
      <div className="grid gap-x-10 gap-y-6 border-b border-hair pb-5 sm:grid-cols-2 lg:grid-cols-4" aria-live="polite">
        <Stat label="Sensitivity" value={a.sensitivity} />
        <Stat label="Specificity" value={a.specificity} />
        <Stat label="VME" value={a.vme} bar={BARS.vme} />
        <Stat label="ME" value={a.me} bar={BARS.me} />
      </div>
      <p className="type-note mt-3">
        Marker calls here use every isolate, with no held-out split, so they describe agreement on these isolates and do not estimate how a test would perform on new ones.
      </p>
      <div className="mt-5 grid gap-x-10 gap-y-6 lg:grid-cols-2">
        <MarkerList title="Model features" markers={featureMarkers} chosen={chosen} onToggle={toggle} />
        <MarkerList title="Catalog gene groups" markers={groupMarkers} chosen={chosen} onToggle={toggle} />
      </div>
      {chosen.size > 0 ? (
        <div className="mt-3">
          <button
            type="button"
            onClick={() => setChosen(new Set())}
            className={`rounded-sm text-control font-medium text-blue underline decoration-blue/30 underline-offset-2 hover:decoration-blue ${FOCUS_RING}`}
          >
            Clear markers
          </button>
        </div>
      ) : null}
    </div>
  );
};

const markerWord = (k: number): string => (k === 1 ? "one-marker" : `${k}-marker`);
