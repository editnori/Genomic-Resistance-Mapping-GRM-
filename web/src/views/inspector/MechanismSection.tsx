import type { DatasetReport, GroupRow, PanelResult } from "../../data/schema";
import { int, pct, plural } from "../../data/format";
import { BARS } from "../../data/metrics";
import { catalogRuleWords, featureLabel, featureWhere, unexplained } from "../../data/derive/dataset";
import { notScored } from "../../data/derive/tests";
import { PanelDesigner } from "./PanelDesigner";
import { DataTable, type Column } from "../../ui/DataTable";
import { StatusTd, PanelState, ScrollFrame, Section } from "../../ui/primitives";

const groupColumns: Column<GroupRow>[] = [
  { id: "group", header: "Gene group", sortValue: (g) => g.group, firstDir: "asc", cell: (g) => g.group },
  { id: "class", header: "Class", sortValue: (g) => g.class, firstDir: "asc", cell: (g) => <span className="text-muted-text">{g.class.toLowerCase()}</span> },
  { id: "subclass", header: "Subclass", sortValue: (g) => g.subclass, firstDir: "asc", cell: (g) => <span className="whitespace-normal break-words text-muted-text">{g.subclass.toLowerCase()}</span> },
  { id: "r", header: "Resistant carriers", align: "right", sortValue: (g) => g.n_resistant, cell: (g) => int(g.n_resistant) },
  { id: "s", header: "Susceptible carriers", align: "right", sortValue: (g) => g.n_susceptible, cell: (g) => int(g.n_susceptible) },
  {
    id: "rule",
    header: "Counts for the rule",
    tip: [{ term: "Counts for the rule", text: "Whether the gene's class or subclass satisfies the catalog rule for this drug." }],
    sortValue: (g) => (g.counts_for_rule ? 1 : 0),
    cell: (g) => (g.counts_for_rule ? "Yes" : <span className="text-muted-text">No</span>),
  },
];

const PanelBlock = ({ r, panel }: { r: DatasetReport; panel: PanelResult }) => {
  const u = unexplained(r, panel);
  const m = panel.metrics;
  return (
    <article className="border-t border-hair py-4">
      <h4 className="text-body font-semibold text-ink">{panel.name}</h4>
      <p className="text-note text-muted-text">Mirrors {panel.mirrors}.</p>
      {u.perSample ? (
        <p className="mt-3 text-body text-ink">
          The panel&apos;s targets are absent from {int(u.missed)} of {int(m.tp + m.fn)} resistant isolates
          {panel.detected.length > 0 ? `; found: ${panel.detected.map(([g, res]) => `${g} in ${int(res)}`).join(", ")} resistant isolates` : "; none of its targets occurs here"}.
        </p>
      ) : (
        <p className="mt-3 text-body text-ink">
          The probes flag {int(m.tp)} of {int(m.tp + m.fn)} resistant isolates and {int(m.fp)} of {int(m.fp + m.tn)} susceptible ones. The report carries no per-isolate probe calls, so the {int(u.missed)} resistant isolates the probes miss cannot be examined here.
        </p>
      )}
      {panel.no_call !== undefined ? (
        <p className="mt-2 text-body text-ink">
          {panel.no_call.resistant + panel.no_call.susceptible === 0
            ? `Every one of the ${int(r.counts.samples)} isolates gets a probe call.`
            : `The probes give no call for ${int(panel.no_call.resistant)} resistant and ${int(panel.no_call.susceptible)} susceptible isolates, whose genomes carry none of the probe k-mers, so sensitivity and specificity cover ${int(m.n)} of ${int(r.counts.samples)} isolates.`}
        </p>
      ) : null}
      {u.perSample && u.missed > 0 ? (
        u.features.length === 0 ? (
          <p className="type-note mt-1">No model feature is more common in the missed isolates than in susceptible ones.</p>
        ) : (
          <div className="mt-3">
            <h5 className="text-control font-medium text-ink">Model features the {int(u.missed)} missed isolates carry more often than susceptible isolates do</h5>
            <ScrollFrame label={`Features in isolates ${panel.name} misses`}>
              <table className="mt-1 w-full border-collapse text-control">
                <thead>
                  <tr className="border-b border-hair">
                    {["Feature", "Where", "Missed isolates carrying it", "Susceptible carrying it"].map((h, i) => (
                      <th key={h} scope="col" className={`px-2 py-2 text-tiny font-medium text-muted-text ${i >= 2 ? "text-right" : "text-left"}`}>
                        {h}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {u.features.map((f) => (
                    <tr key={f.feature} className="border-b border-hair">
                      <td className="px-2 py-1.5 font-medium">
                        {featureLabel(f.card, f.feature)}
                        {f.alike > 0 ? <span className="block text-tiny font-normal text-muted-text">and {plural(f.alike, "more feature")} in the same gene</span> : null}
                      </td>
                      <td className="px-2 py-1.5 text-muted-text">{featureWhere(f.card, f.feature) ?? "Not placed"}</td>
                      <td className="px-2 py-1.5 text-right tabular-nums">
                        {int(f.inMissed)} of {int(u.missed)} ({pct(f.inMissed / Math.max(u.missed, 1), 0)})
                      </td>
                      <td className="px-2 py-1.5 text-right tabular-nums">
                        {int(f.inSusceptible)} of {int(u.susceptible)} ({pct(f.inSusceptible / Math.max(u.susceptible, 1), 1)})
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </ScrollFrame>
          </div>
        )
      ) : null}
    </article>
  );
};

export const MechanismSection = ({ r }: { r: DatasetReport }) => {
  const cat = r.catalog;
  const unsupported = r.unsupported_assays ?? [];
  const tests = [...(cat === null ? [] : [{ id: "catalog", name: "Catalog", metrics: cat.metrics }]), ...r.panels];
  return (
    <Section
      id="ds-mechanisms"
      title="Genes and panels"

      answer={
        cat === null
          ? `The catalog has no acquired-gene rule for ${r.antibiotic}.`
          : catalogRuleWords(cat.rule)
      }
    >
      {cat === null ? null : <p className="type-note">Catalog-containment rule · AMRFinderPlus catalog {cat.version}</p>}
      {tests.length === 0 ? null : (
        <ScrollFrame label="Gene and panel results">
          <table className="mt-3 w-full border-collapse text-control">
            <thead>
              <tr className="border-b border-hair">
                {["Test", "Sensitivity", "Specificity", "VME ≤1.5%", "ME ≤3%"].map((h, i) => <th key={h} scope="col" className={`whitespace-nowrap px-2 py-2 font-medium text-muted-text ${i === 0 ? "text-left" : "text-right"}`}>{h}</th>)}
              </tr>
            </thead>
            <tbody>
              {tests.map((t) => <tr key={t.id} className="border-b border-hair">
                <th scope="row" className="px-2 py-2 text-left font-medium">{t.name}</th>
                <td className="px-2 py-2 text-right tabular-nums">{pct(t.metrics.sensitivity)}</td>
                <td className="px-2 py-2 text-right tabular-nums">{pct(t.metrics.specificity)}</td>
                <StatusTd value={t.metrics.vme} bar={BARS.vme} />
                <StatusTd value={t.metrics.me} bar={BARS.me} />
              </tr>)}
            </tbody>
          </table>
        </ScrollFrame>
      )}
      <div className="mt-6">
        <h3 className="text-control font-semibold text-ink">Simulated PCR panels</h3>
        {r.panels.length === 0 && unsupported.length === 0 ? (
          <div className="mt-2">
            <PanelState kind="empty" title={`No simulated panel result is available for ${r.antibiotic}`} message="Panels cover carbapenems, extended-spectrum cephalosporins, methicillin, vancomycin, and colistin." />
          </div>
        ) : (
          <>
            <p className="type-note mt-1">
              A gene panel calls an isolate resistant when its assembled genome carries a target gene, detected by k-mer containment; a probe assay, when a probed region differs from the wild-type sequence. This is target presence in genomes, not an estimate of a PCR assay&apos;s analytical performance.
            </p>
            {r.panels.map((p) => (
              <PanelBlock key={p.id} r={r} panel={p} />
            ))}
            {unsupported.map((a) => (
              <article key={a.id} className="border-t border-hair py-4">
                <h4 className="text-body font-semibold text-ink">{a.name}</h4>
                <p className="text-note text-muted-text">Mirrors {a.mirrors}.</p>
                <p className="mt-2 text-body text-ink">{notScored(a.reason)} No sensitivity or specificity is reported for this assay.</p>
              </article>
            ))}
          </>
        )}
      </div>
      <PanelDesigner r={r} />
      {cat === null ? null : (
        <div className="mt-8">
          <h3 className="mb-3 text-body font-semibold text-ink">Catalog gene groups</h3>
          <DataTable label="Catalog gene groups" rows={cat.groups} columns={groupColumns} rowKey={(g) => g.group} initialSort={{ id: "r", dir: "desc" }} empty={<p className="type-note">No catalog gene occurs in these isolates.</p>} />
        </div>
      )}
    </Section>
  );
};
