// Molecular-test findings from the manifest's panel rows. Pooled shares are
// weighted by each dataset's resistant isolates, so a large dataset counts
// for its size; panel rows carry exact counts, and older manifests carry only
// rates, scaled by the catalog row's n.

import { capitalize, int } from "../format";
import type { CatalogRow, Manifest, PanelSummaryRow } from "../schema";

export interface PanelRowJoined extends PanelSummaryRow {
  row: CatalogRow | undefined;
  /** Resistant isolates in the dataset, and those the panel calls resistant. */
  resistant: number;
  caught: number;
  calledResistant: number;
  modelResistant: number;
  modelCaught: number;
  susceptible: number;
  calledSusceptible: number;
  falseCalls: number;
  noCallResistant: number;
  noCallSusceptible: number;
}

export const joinPanels = (m: Manifest): PanelRowJoined[] => {
  const byId = new Map(m.datasets.map((d) => [d.id, d]));
  return m.panels.map((p) => {
    const row = byId.get(p.dataset);
    const c = p.counts;
    const model = p.best_model_counts;
    const x = p.no_call ?? { resistant: 0, susceptible: 0 };
    const calledResistant = c !== undefined ? c.tp + c.fn : (row?.n_resistant ?? 0);
    const calledSusceptible = c !== undefined ? c.tn + c.fp : row === undefined ? 0 : row.n - row.n_resistant;
    const resistant = calledResistant + x.resistant;
    return {
      ...p,
      row,
      resistant,
      caught: c !== undefined ? c.tp : Math.round(p.metrics.sensitivity * calledResistant),
      calledResistant,
      modelResistant: model !== undefined ? model.tp + model.fn : resistant,
      modelCaught: model !== undefined ? model.tp : Math.round(p.best_model.sensitivity * resistant),
      susceptible: calledSusceptible + x.susceptible,
      calledSusceptible,
      falseCalls: c !== undefined ? c.fp : Math.round((1 - p.metrics.specificity) * calledSusceptible),
      noCallResistant: x.resistant,
      noCallSusceptible: x.susceptible,
    };
  });
};

export interface Pool {
  datasets: number;
  resistant: number;
  caught: number;
  calledResistant: number;
  modelResistant: number;
  modelCaught: number;
  susceptible: number;
  calledSusceptible: number;
  falseCalls: number;
  noCallResistant: number;
  noCallSusceptible: number;
}

export const pool = (rows: readonly PanelRowJoined[]): Pool =>
  rows.reduce<Pool>(
    (acc, r) => ({
      datasets: acc.datasets + 1,
      resistant: acc.resistant + r.resistant,
      caught: acc.caught + r.caught,
      calledResistant: acc.calledResistant + r.calledResistant,
      modelResistant: acc.modelResistant + r.modelResistant,
      modelCaught: acc.modelCaught + r.modelCaught,
      susceptible: acc.susceptible + r.susceptible,
      calledSusceptible: acc.calledSusceptible + r.calledSusceptible,
      falseCalls: acc.falseCalls + r.falseCalls,
      noCallResistant: acc.noCallResistant + r.noCallResistant,
      noCallSusceptible: acc.noCallSusceptible + r.noCallSusceptible,
    }),
    { datasets: 0, resistant: 0, caught: 0, calledResistant: 0, modelResistant: 0, modelCaught: 0, susceptible: 0, calledSusceptible: 0, falseCalls: 0, noCallResistant: 0, noCallSusceptible: 0 },
  );

/** Datasets with a simulated panel but no best learner. report.rs writes a
 *  panel row only beside a best learner, so the panel rows leave them out. */
export const panelsWithoutBest = (m: Manifest): number => m.datasets.filter((d) => d.best_learner === null && d.panels.length > 0).length;

export const noCallWords = (resistant: number, susceptible: number): string =>
  resistant + susceptible === 0 ? "None" : `${int(resistant)} resistant, ${int(susceptible)} susceptible`;

export const notScored = (reason: string): string => `Not scored. ${capitalize(reason)}${reason.endsWith(".") ? "" : "."}`;

export const share = (part: number, whole: number): number | null => (whole === 0 ? null : part / whole);

/** Distinct datasets behind a set of panel rows, and their resistant
 *  isolates counted once each (two panels on one dataset score the same
 *  isolates twice in `pool`). */
export const distinctResistant = (rows: readonly PanelRowJoined[]): { datasets: number; resistant: number } => {
  const seen = new Map<string, number>();
  rows.forEach((r) => seen.set(r.dataset, r.resistant));
  return { datasets: seen.size, resistant: [...seen.values()].reduce((a, b) => a + b, 0) };
};

/** The panel the panels finding draws: the carbapenemase panel when it ran,
 *  else the panel covering the most resistant isolates; its datasets with the
 *  most resistant isolates first. */
export const leadPanel = (m: Manifest): { id: string; name: string; rows: PanelRowJoined[] } | null => {
  const all = joinPanels(m);
  const ids = [...new Set(all.map((r) => r.panel))];
  if (ids.length === 0) return null;
  const size = (id: string) => all.filter((r) => r.panel === id).reduce((a, r) => a + r.resistant, 0);
  const id = ids.includes("carbapenemase") ? "carbapenemase" : ids.sort((a, b) => size(b) - size(a))[0]!;
  const rows = all.filter((r) => r.panel === id).sort((a, b) => b.resistant - a.resistant);
  return { id, name: rows[0]?.name ?? id, rows };
};

export interface ClassPool extends Pool {
  drugClass: string;
}

/** Per drug class, the share of phenotypic resistance whose genome carries a
 *  panel's targets. Each dataset counts once, with the panel of highest
 *  balanced accuracy, so a broad panel whose extra targets also flag many
 *  susceptible isolates does not win on catches alone. */
export const classPools = (rows: readonly PanelRowJoined[]): ClassPool[] => {
  const ba = (r: PanelRowJoined) => r.metrics.sensitivity + r.metrics.specificity;
  const best = new Map<string, PanelRowJoined>();
  rows.forEach((r) => {
    const cur = best.get(r.dataset);
    if (cur === undefined || ba(r) > ba(cur)) best.set(r.dataset, r);
  });
  const by = new Map<string, PanelRowJoined[]>();
  best.forEach((r) => {
    const cls = r.row?.drug_class ?? "unknown";
    by.set(cls, [...(by.get(cls) ?? []), r]);
  });
  return [...by.entries()]
    .map(([drugClass, group]) => ({ drugClass, ...pool(group) }))
    .sort((a, b) => (share(b.caught, b.calledResistant) ?? 0) - (share(a.caught, a.calledResistant) ?? 0));
};
