// The rule showcase: the dataset whose set covering machine scores best on
// random splits, its first rule, where the rule's k-mers sit, the contig
// around them, and how often resistant and susceptible isolates carry them.
// The sequence shows inline beside the rule; on a slide too short for the
// sequence and the carriage bars together, the sequence continues on the
// next slide, so nothing opens a second dialog and nothing scrolls.
import type { CatalogRow, Manifest } from "../data/schema";
import { capitalize, int, pct, speciesShort } from "../data/format";
import { featureIndex, featureLabel, featurePlace } from "../data/derive/dataset";
import { loadErrorText, useDataset } from "../data/load";
import { PanelState } from "../ui/primitives";
import { Seq } from "../ui/Seq";

/** The dataset whose set covering machine scores best on random splits. */
export const ruleDataset = (m: Manifest): { row: CatalogRow; ba: number } | null =>
  m.datasets
    .flatMap((row) => {
      const ba = row.learners.find((l) => l.learner === "scm")?.random?.balanced_accuracy;
      return ba === undefined || ba === null || row.scm_rules.length === 0 ? [] : [{ row, ba }];
    })
    .sort((a, b) => b.ba - a.ba || b.row.n - a.row.n)[0] ?? null;

/** The same set covering machine's balanced accuracy with lineages held
 *  out, beside its random-split score; null when the split did not run. */
export const ruleLineageScore = (row: CatalogRow): number | null => row.learners.find((l) => l.learner === "scm")?.lineage?.balanced_accuracy ?? null;

/** `pages` is the deck's page count for this figure: 2 splits the carriage
 *  bars (page 0) from the sequence (page 1); 1 shows both. Undefined is the
 *  deck's off-screen measuring copy, which shows both and no status panel. */
export const RuleFigure = ({ manifest, page = 0, pages }: { manifest: Manifest; page?: number; pages?: number }) => {
  const pick = ruleDataset(manifest);
  const report = useDataset(pick?.row.id ?? null);
  const measuring = pages === undefined;
  if (pick === null) return measuring ? null : <PanelState kind="empty" title="No set covering machine rule in this report" />;
  if (report.status === "loading") return measuring ? null : <PanelState kind="loading" title="Loading the dataset report" />;
  if (report.status === "error") return measuring ? null : <PanelState kind="error" title="The dataset report could not be read" message={loadErrorText(report.error)} />;
  const r = report.data;
  const rule = r.models.scm.rules[0];
  const cards = featureIndex(r);
  const card = rule === undefined ? undefined : cards.get(rule.feature);
  if (rule === undefined || card === undefined) return measuring ? null : <PanelState kind="empty" title="The first rule's feature has no annotation" />;
  const rShare = card.n_resistant / Math.max(r.counts.resistant, 1);
  const sShare = card.n_susceptible / Math.max(r.counts.susceptible, 1);
  const context = card.locus?.context ?? "";
  const split = pages !== undefined && pages > 1 && context !== "";
  const header = (
    <>
      <p className="text-note text-muted-text">
        {capitalize(r.antibiotic)} · {speciesShort(r.species)} · rule 1 of {r.models.scm.rules.length}
      </p>
      <p className="mt-2 break-words text-headline font-semibold text-ink">
        {rule.presence ? "Presence" : "Absence"} of {featureLabel(card, rule.feature)}
      </p>
      {card.catalog_families.length === 0 ? null : (
        <p className="mt-2 text-note text-muted-text">
          Catalog {card.catalog_families.length === 1 ? "family" : "families"} {card.catalog_families.join(", ")}
        </p>
      )}
      <p className="mt-2 text-note text-ink">{featurePlace(card) ?? "Not placed on a reference assembly."}</p>
    </>
  );
  const sequence =
    context === "" ? null : (
      <div className="mt-3">
        <p className="text-note text-muted-text">The contig around the rule&apos;s k-mers; the k-mers are in capitals.</p>
        <Seq context={context} className="mt-1" />
      </div>
    );
  const bars = [
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
      <p className="mt-1 text-note text-muted-text">
        {int(b.n)} of {int(b.of)}
      </p>
    </div>
  ));
  if (split) {
    return (
      <div className="grid h-full grid-cols-1 content-center gap-y-4">
        <div className="min-w-0">
          {header}
          {page === 1 ? sequence : null}
        </div>
        {page === 0 ? (
          <>
            <div className="space-y-4">{bars}</div>
            <p className="text-note text-muted-text">The rule&apos;s sequence follows on the next slide.</p>
          </>
        ) : null}
      </div>
    );
  }
  // data-continues tells the deck's measure that this figure can move its
  // sequence onto a second slide; without a sequence there is nothing to move.
  return (
    <div data-continues={sequence === null ? undefined : ""} className="grid h-full grid-cols-1 content-center gap-x-12 gap-y-4 lg:grid-cols-2">
      <div className="min-w-0">
        {header}
        {sequence}
      </div>
      <div className="space-y-4 self-center">{bars}</div>
    </div>
  );
};
