// The rule showcase retains its dataset selection and carriage calculations.
import { useState } from "react";
import type { CatalogRow, Manifest } from "../data/schema";
import { capitalize, int, pct, speciesShort } from "../data/format";
import { featureIndex, featureLabel, featurePlace } from "../data/derive/dataset";
import { loadErrorText, useDataset } from "../data/load";
import { GhostButton, PanelState } from "../ui/primitives";
import { ModalDialog } from "../ui/ModalDialog";
import { Seq } from "../ui/Seq";

/** The dataset whose set covering machine scores best on random splits. */
export const ruleDataset = (m: Manifest): { row: CatalogRow; ba: number } | null =>
  m.datasets
    .flatMap((row) => {
      const ba = row.learners.find((l) => l.learner === "scm")?.random?.balanced_accuracy;
      return ba === undefined || ba === null || row.scm_rules.length === 0 ? [] : [{ row, ba }];
    })
    .sort((a, b) => b.ba - a.ba || b.row.n - a.row.n)[0] ?? null;

export const RuleFigure = ({ manifest }: { manifest: Manifest }) => {
  const pick = ruleDataset(manifest);
  const report = useDataset(pick?.row.id ?? null);
  const [sequenceOpen, setSequenceOpen] = useState(false);
  if (pick === null) return <PanelState kind="empty" title="No set covering machine rule in this report" />;
  if (report.status === "loading") return <PanelState kind="loading" title="Loading the dataset report" />;
  if (report.status === "error") return <PanelState kind="error" title="The dataset report could not be read" message={loadErrorText(report.error)} />;
  const r = report.data;
  const rule = r.models.scm.rules[0];
  const cards = featureIndex(r);
  const card = rule === undefined ? undefined : cards.get(rule.feature);
  if (rule === undefined || card === undefined) return <PanelState kind="empty" title="The first rule's feature has no annotation" />;
  const rShare = card.n_resistant / Math.max(r.counts.resistant, 1);
  const sShare = card.n_susceptible / Math.max(r.counts.susceptible, 1);
  const context = card.locus?.context ?? "";
  return (
    <div className="grid h-full grid-cols-1 content-center gap-x-12 gap-y-4 lg:grid-cols-2">
      <div className="min-w-0">
        <p className="text-note text-muted-text">{capitalize(r.antibiotic)} · {speciesShort(r.species)} · rule 1 of {r.models.scm.rules.length}</p>
        <p className="mt-2 break-words text-headline font-semibold text-ink">{rule.presence ? "Presence" : "Absence"} of {featureLabel(card, rule.feature)}</p>
        {card.catalog_families.length === 0 ? null : <p className="mt-2 text-note text-muted-text">Catalog {card.catalog_families.length === 1 ? "family" : "families"} {card.catalog_families.join(", ")}</p>}
        <p className="mt-2 text-note text-ink">{featurePlace(card) ?? "Not placed on a reference assembly."}</p>
        {context === "" ? null : <div className="mt-3"><GhostButton onClick={() => setSequenceOpen(true)}>Show sequence</GhostButton></div>}
      </div>
      <div className="space-y-4 self-center">
        {[
          { label: "Resistant isolates carrying it", n: card.n_resistant, of: r.counts.resistant, share: rShare, color: "var(--series-r)" },
          { label: "Susceptible isolates carrying it", n: card.n_susceptible, of: r.counts.susceptible, share: sShare, color: "var(--series-s)" },
        ].map((b) => (
          <div key={b.label}>
            <div className="flex items-baseline justify-between gap-4">
              <span className="text-note text-ink">{b.label}</span>
              <span className="text-headline font-semibold tabular-nums text-ink">{pct(b.share)}</span>
            </div>
            <div className="mt-2 h-3 overflow-hidden rounded-sm bg-surface-subtle">
              <div className="h-full rounded-r-sm" style={{ background: b.color, width: `${Math.max(0.5, b.share * 100)}%` }} />
            </div>
            <p className="mt-1 text-note text-muted-text">{int(b.n)} of {int(b.of)}</p>
          </div>
        ))}
      </div>
      {sequenceOpen ? <ModalDialog label="Rule sequence" onClose={() => setSequenceOpen(false)}>
        <div className="absolute inset-0 bg-paper/60" onClick={() => setSequenceOpen(false)} />
        <div className="absolute inset-x-4 top-1/2 mx-auto max-h-[80dvh] max-w-3xl -translate-y-1/2 overflow-y-auto rounded-xl border border-hair bg-paper-pure p-6 shadow-xl">
          <div className="flex items-center justify-between gap-4"><h2 className="text-headline font-semibold">Rule sequence</h2><GhostButton onClick={() => setSequenceOpen(false)}>Close sequence</GhostButton></div>
          <p className="type-note mt-3">Contig around the rule's k-mers; the span is in capitals.</p>
          <Seq context={context} className="mt-3" />
        </div>
      </ModalDialog> : null}
    </div>
  );
};
