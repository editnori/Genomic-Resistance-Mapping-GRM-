// Panel evidence, under the Overview panel figure: every simulated panel
// pooled, every drug class, and every panel on every dataset, filtered by
// species, drug class, and panel. Rows open the dataset inspector.

import { useState } from "react";
import type { CatalogRow, Manifest } from "../data/schema";
import { capitalize, int, LEARNER_LABEL, pct, plural, speciesShort } from "../data/format";
import { applyAxes, countBy, memberTest, type Axis } from "../data/facets";
import { BARS } from "../data/metrics";
import { classPools, distinctResistant, joinPanels, noCallWords, panelsWithoutBest, pool, share, type PanelRowJoined } from "../data/derive/tests";
import { BEST_LEARNER_WORDS } from "../data/derive/dataset";
import { noBestWords } from "../data/derive/overview";
import { Dropdown } from "../ui/Dropdown";
import { DataTable, type Column } from "../ui/DataTable";
import { GhostButton, ScrollFrame } from "../ui/primitives";

interface Filters {
  species: Set<string>;
  classes: Set<string>;
  panels: Set<string>;
}

const ARRIVAL: Filters = { species: new Set(), classes: new Set(), panels: new Set() };

const axesFor = (f: Filters): Axis<PanelRowJoined>[] => [
  { id: "species", test: memberTest(f.species, (r) => r.species) },
  { id: "class", test: memberTest(f.classes, (r) => r.row?.drug_class ?? "") },
  { id: "panel", test: memberTest(f.panels, (r) => r.panel) },
];

const found = (caught: number, of: number): string => `${int(caught)} of ${int(of)} (${pct(share(caught, of))})`;

const rowColumns: Column<PanelRowJoined>[] = [
  {
    id: "dataset",
    header: "Dataset",
    sortValue: (r) => `${r.antibiotic} ${r.species}`,
    firstDir: "asc",
    cell: (r) => (
      <span className="block">
        <span className="block">{capitalize(r.antibiotic)}</span>
        <span className="block text-tiny font-normal italic text-muted-text">{speciesShort(r.species)}</span>
      </span>
    ),
  },
  { id: "panel", header: "Panel", sortValue: (r) => r.name, firstDir: "asc", className: "whitespace-nowrap", cell: (r) => r.name },
  {
    id: "found",
    header: "Found",
    unit: "resistant",
    align: "right",
    sortValue: (r) => share(r.caught, r.calledResistant),
    firstDir: "asc",
    cell: (r) => (
      <span className="block whitespace-nowrap">
        <span className="block">{pct(share(r.caught, r.calledResistant))}</span>
        <span className="block text-tiny text-muted-text">
          {int(r.caught)} of {int(r.calledResistant)}
        </span>
      </span>
    ),
  },
  {
    id: "spec",
    header: "Specificity",
    align: "right",
    sortValue: (r) => r.metrics.specificity,
    cell: (r) => pct(r.metrics.specificity),
  },
  {
    id: "nocall",
    header: "No call",
    align: "right",
    sortValue: (r) => r.noCallResistant + r.noCallSusceptible,
    className: "whitespace-nowrap",
    cell: (r) => (r.noCallResistant + r.noCallSusceptible === 0 ? <span className="text-muted-text">None</span> : noCallWords(r.noCallResistant, r.noCallSusceptible)),
  },
  { id: "vme", header: "VME", unit: "≤1.5%", align: "right", sortValue: (r) => r.metrics.vme, firstDir: "asc", status: (r) => ({ value: r.metrics.vme, bar: BARS.vme }), cell: (r) => pct(r.metrics.vme) },
  { id: "me", header: "ME", unit: "≤3%", align: "right", sortValue: (r) => r.metrics.me, firstDir: "asc", status: (r) => ({ value: r.metrics.me, bar: BARS.me }), cell: (r) => pct(r.metrics.me) },
  {
    id: "learner",
    header: "Learner",
    align: "right",
    sortValue: (r) => r.best_model.sensitivity,
    cell: (r) => (
      <span className="block whitespace-nowrap">
        <span className="block">{pct(r.best_model.sensitivity)}</span>
        <span className="block text-tiny text-muted-text">{LEARNER_LABEL[r.best_learner]}</span>
      </span>
    ),
  },
];

const catalogShare = (rows: readonly CatalogRow[]): { caught: number; resistant: number } =>
  rows.reduce(
    (acc, r) => (r.catalog === null ? acc : { caught: acc.caught + Math.round(r.catalog.sensitivity * r.n_resistant), resistant: acc.resistant + r.n_resistant }),
    { caught: 0, resistant: 0 },
  );

const TH = "whitespace-nowrap px-2 py-2 text-tiny font-medium text-blue";

export const PanelEvidence = ({ manifest, onOpen }: { manifest: Manifest; onOpen?(id: string): void }) => {
  const [f, setF] = useState<Filters>(ARRIVAL);
  const all = joinPanels(manifest);
  const axes = axesFor(f);
  const rows = applyAxes(all, axes);
  const active = f.species.size + f.classes.size + f.panels.size > 0;
  const speciesCounts = countBy(all, axes, "species", (r) => r.species);
  const classCounts = countBy(all, axes, "class", (r) => r.row?.drug_class ?? "");
  const panelCounts = countBy(all, axes, "panel", (r) => r.panel);
  const panelNames = new Map(all.map((r) => [r.panel, r.name]));
  const byPanel = [...new Set(rows.map((r) => r.panel))].map((id) => ({ id, name: panelNames.get(id) ?? id, ...pool(rows.filter((r) => r.panel === id)) })).sort((a, b) => b.resistant - a.resistant);
  const classes = classPools(rows);
  const datasetRows = new Map(manifest.datasets.map((d) => [d.id, d]));
  const pooled = pool(rows);
  const distinct = distinctResistant(rows);

  if (all.length === 0) return <p className="prose-measure">No panel targets an analyzed drug yet.</p>;
  return (
    <div className="space-y-6">
      <div className="flex flex-wrap items-center gap-2" role="group" aria-label="Filter panel evidence">
        <Dropdown
          ariaLabel="Species"
          rest="All species"
          options={[...new Set(all.map((r) => r.species))].sort().map((s) => ({ id: s, label: s, count: speciesCounts.get(s) ?? 0 }))}
          mode={{ kind: "many", value: f.species, onChange: (species) => setF({ ...f, species }) }}
        />
        <Dropdown
          ariaLabel="Drug class"
          rest="All drug classes"
          options={[...new Set(all.map((r) => r.row?.drug_class ?? ""))].sort().map((c) => ({ id: c, label: capitalize(c), count: classCounts.get(c) ?? 0 }))}
          mode={{ kind: "many", value: f.classes, onChange: (classes) => setF({ ...f, classes }) }}
        />
        <Dropdown
          ariaLabel="Panel"
          rest="All panels"
          options={[...panelNames.entries()].map(([id, name]) => ({ id, label: name, count: panelCounts.get(id) ?? 0, disabled: (panelCounts.get(id) ?? 0) === 0 && !f.panels.has(id) }))}
          mode={{ kind: "many", value: f.panels, onChange: (panels) => setF({ ...f, panels }) }}
        />
        <span className="text-control text-muted-text" aria-live="polite">
          {active ? `${int(rows.length)} of ${plural(all.length, "panel result")}` : plural(all.length, "panel result")}
        </span>
        {active ? <GhostButton onClick={() => setF(ARRIVAL)}>Clear</GhostButton> : null}
      </div>

      {rows.length === 0 ? (
        <p className="prose-measure">No panel result matches these filters.</p>
      ) : (
        <>
          <p className="prose-measure">
            Across {plural(byPanel.length, "panel")} and {plural(pooled.datasets, "panel and dataset pair")} over {plural(distinct.datasets, "dataset")} with {int(distinct.resistant)} distinct resistant isolates, targets are found in{" "}
            {pct(share(pooled.caught, pooled.calledResistant))} of {int(pooled.calledResistant)} resistant evaluations with a panel call; {BEST_LEARNER_WORDS} calls {pct(share(pooled.modelCaught, pooled.modelResistant))} of {int(pooled.modelResistant)} resistant.
            {pooled.noCallResistant > 0 ? ` ${plural(pooled.noCallResistant, "resistant evaluation")} with no probe call ${pooled.noCallResistant === 1 ? "is" : "are"} left out of the panel share and kept in the learner's.` : ""} An isolate scored by two panels counts once per panel.
            {noBestWords(panelsWithoutBest(manifest), " with a simulated panel")}
          </p>

          <div>
            <h3 className="mb-2 text-control font-semibold text-ink">Panels</h3>
            <ScrollFrame label="Panels">
              <table className="w-full border-collapse text-control">
                <thead>
                  <tr className="border-b border-hair">
                    {["Panel", "Datasets", "Resistant", "Found", "No call", "Learner", "False positives"].map((h, i) => (
                      <th key={h} scope="col" className={`${TH} ${i === 0 ? "text-left" : "text-right"}`}>
                        {h}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {byPanel.map((p) => (
                    <tr key={p.id} className="border-b border-hair">
                      <th scope="row" className="whitespace-nowrap px-2 py-2 text-left font-medium text-ink">
                        {p.name}
                      </th>
                      <td className="px-2 py-2 text-right tabular-nums">{int(p.datasets)}</td>
                      <td className="px-2 py-2 text-right tabular-nums">{int(p.resistant)}</td>
                      <td className="whitespace-nowrap px-2 py-2 text-right tabular-nums">{found(p.caught, p.calledResistant)}</td>
                      <td className="whitespace-nowrap px-2 py-2 text-right tabular-nums">{noCallWords(p.noCallResistant, p.noCallSusceptible)}</td>
                      <td className="whitespace-nowrap px-2 py-2 text-right tabular-nums">{found(p.modelCaught, p.modelResistant)}</td>
                      <td className="whitespace-nowrap px-2 py-2 text-right tabular-nums">
                        {int(p.falseCalls)} of {int(p.calledSusceptible)}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </ScrollFrame>
            <p className="type-note mt-2">
              Found is resistant isolates whose genome carries a target, over those the panel calls; Learner is resistant isolates the learner calls resistant, over all it scores; False positives is susceptible isolates with a target. Pooled over datasets, weighting each by its resistant isolates; counts are rebuilt from rates and size where a row carries none.
            </p>
          </div>

          <div>
            <h3 className="mb-2 text-control font-semibold text-ink">Drug classes</h3>
            <ScrollFrame label="Drug classes">
              <table className="w-full border-collapse text-control">
                <thead>
                  <tr className="border-b border-hair">
                    {["Drug class", "Datasets", "Resistant", "No call", "Found", "False positives", "Catalog"].map((h, i) => (
                      <th key={h} scope="col" className={`${TH} ${i === 0 ? "text-left" : "text-right"}`}>
                        {h}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {classes.map((c) => {
                    const ids = new Set(rows.filter((r) => (r.row?.drug_class ?? "") === c.drugClass).map((r) => r.dataset));
                    const inClass = [...ids].flatMap((id) => (datasetRows.has(id) ? [datasetRows.get(id)!] : []));
                    const ruled = inClass.filter((d) => d.catalog !== null);
                    const whole = catalogShare(ruled);
                    return (
                      <tr key={c.drugClass} className="border-b border-hair">
                        <th scope="row" className="whitespace-nowrap px-2 py-2 text-left font-medium text-ink">
                          {capitalize(c.drugClass)}
                        </th>
                        <td className="px-2 py-2 text-right tabular-nums">{int(c.datasets)}</td>
                        <td className="px-2 py-2 text-right tabular-nums">{int(c.resistant)}</td>
                        <td className="whitespace-nowrap px-2 py-2 text-right tabular-nums">{noCallWords(c.noCallResistant, c.noCallSusceptible)}</td>
                        <td className="px-2 py-2 text-right tabular-nums">{pct(share(c.caught, c.calledResistant))}</td>
                        <td className="px-2 py-2 text-right tabular-nums">{pct(share(c.falseCalls, c.calledSusceptible))}</td>
                        <td className="px-2 py-2 text-right tabular-nums">
                          {whole.resistant === 0 ? <span className="text-muted-text">No rule</span> : pct(share(whole.caught, whole.resistant))}
                          {ruled.length > 0 && ruled.length < inClass.length ? <span className="block text-tiny text-muted-text">{int(ruled.length)} of {int(inClass.length)} datasets</span> : null}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </ScrollFrame>
            <p className="type-note mt-2">
              Each dataset counts once, with its panel of highest balanced accuracy. Catalog is the catalog rule&apos;s sensitivity on the datasets of that class with a rule.
            </p>
          </div>

          <div>
            <h3 className="mb-2 text-control font-semibold text-ink">Every panel on every dataset</h3>
            <DataTable
              label="Panel results"
              rows={rows}
              columns={rowColumns}
              rowKey={(r) => `${r.panel}|${r.dataset}`}
              onOpen={onOpen === undefined ? undefined : (r) => onOpen(r.dataset)}
              initialSort={{ id: "found", dir: "asc" }}
              empty={<p className="type-note">No panel result matches these filters.</p>}
            />
            <p className="type-note mt-2">
              Found is resistant isolates with a panel call whose assembled genome carries a target, over the resistant isolates the panel calls; Specificity is susceptible isolates with a
              panel call and no target, over the susceptible isolates it calls. No call counts isolates whose genome carries none of a probe assay&apos;s k-mers: they are neither found nor
              missed. Learner is the share of all resistant isolates {BEST_LEARNER_WORDS} calls resistant. VME and ME are read against the FDA-style bars in their headers.
            </p>
          </div>
        </>
      )}
    </div>
  );
};
