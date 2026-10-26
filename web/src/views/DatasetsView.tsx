// The dataset catalog: one row per species and antibiotic, in ten compact
// columns. One search row carries the search, the Filters button below
// 1440px, and the one result count; from 1440px the filters are a rail in
// the same grid, its title row on the search row's line. The inspector
// holds every full value; a rule cell expands in place.

import { useState } from "react";
import { Search } from "lucide-react";
import type { CatalogRow, Learner, LearnerLite, Manifest } from "../data/schema";
import { capitalize, int, LEARNER_LABEL, LEARNER_NAME, pct, speciesShort } from "../data/format";
import { applyAxes, countBy, memberTest, rowsWithout, type Axis } from "../data/facets";
import { BARS, verdictOf } from "../data/metrics";
import { BEST_LEARNER_WORDS, NO_BEST_LEARNER, ruleSummary, splitReasonShort } from "../data/derive/dataset";
import { auditOf, manifestAudit } from "../data/derive/evaluation";
import type { Navigator } from "../app/router";
import { ReportStatus } from "../app/ReportStatus";
import { ViewHeader } from "../ui/ViewHeader";
import { Dropdown } from "../ui/Dropdown";
import { DataTable, type Column } from "../ui/DataTable";
import { FilterRail, FiltersButton, railGrid, RailField, RailGroup, useRailState } from "../ui/FilterRail";
import { CountPill, FOCUS_RING, GhostButton, Lozenge, VIEW_SHELL } from "../ui/primitives";
import { useCountUp } from "../ui/useCountUp";

const MIN_N = [
  { id: "0", label: "Any size" },
  { id: "100", label: "100 or more" },
  { id: "250", label: "250 or more" },
  { id: "500", label: "500 or more" },
  { id: "1000", label: "1,000 or more" },
] as const;

const BEST_LEARNERS: Learner[] = ["scm", "cart", "rf", "gbt", "lasso"];

const bestLite = (r: CatalogRow): LearnerLite | null => r.learners.find((l) => l.learner === r.best_learner)?.random ?? null;
const bestName = (r: CatalogRow): string => (r.best_learner === null ? NO_BEST_LEARNER : LEARNER_LABEL[r.best_learner]);
const gap = (r: CatalogRow): number | null => (r.best_random === null || r.best_lineage === null ? null : r.best_random - r.best_lineage);

const meets = (value: number | null | undefined, bar: (typeof BARS)[keyof typeof BARS]): boolean => value !== null && value !== undefined && verdictOf(value, bar) === "meets";

const searchText = (r: CatalogRow): string => [r.antibiotic, r.species, r.drug_class, bestName(r), ...r.scm_rules].join(" ").toLowerCase();

interface Filters {
  species: Set<string>;
  classes: Set<string>;
  minN: string;
  learners: Set<string>;
  mic: boolean;
  vme: boolean;
  me: boolean;
}

const ARRIVAL: Filters = { species: new Set(), classes: new Set(), minN: "0", learners: new Set(), mic: false, vme: false, me: false };

const axesFor = (f: Filters, q: string): Axis<CatalogRow>[] => [
  // The search lives above the table, but the counts still see it.
  { id: "search", test: q === "" ? null : (r) => searchText(r).includes(q) },
  { id: "species", test: memberTest(f.species, (r) => r.species) },
  { id: "class", test: memberTest(f.classes, (r) => r.drug_class) },
  { id: "minN", test: f.minN === "0" ? null : (r) => r.n >= Number(f.minN) },
  { id: "learner", test: memberTest(f.learners, (r) => r.best_learner ?? "none") },
  { id: "mic", test: f.mic ? (r) => r.mic_essential_agreement !== null : null },
  { id: "vme", test: f.vme ? (r) => meets(bestLite(r)?.vme, BARS.vme) : null },
  { id: "me", test: f.me ? (r) => meets(bestLite(r)?.me, BARS.me) : null },
];

/** A pill counts with only its own axis lifted, so its number is the rows
 *  choosing it returns with every other control as it is. */
const pillCount = (rows: readonly CatalogRow[], axes: Axis<CatalogRow>[], lift: string, test: (r: CatalogRow) => boolean): number =>
  applyAxes(rows, axes.filter((a) => a.id !== lift)).filter(test).length;

const Two = ({ top, below, right = false }: { top: string; below: string; right?: boolean }) => (
  <span className={`block whitespace-nowrap ${right ? "text-right" : ""}`}>
    <span className="block">{top}</span>
    <span className="block text-tiny font-normal text-muted-text">{below}</span>
  </span>
);

/** The model's first rule and how many more, on one line; the count opens
 *  every rule in place, one per line. */
const RuleCell = ({ rules }: { rules: readonly string[] }) => {
  const [open, setOpen] = useState(false);
  const s = ruleSummary(rules);
  if (s === null) {
    return (
      <span className="text-muted-text">
        –<span className="sr-only">The set covering machine found no rule</span>
      </span>
    );
  }
  const toggle = `ui-control-target shrink-0 rounded-full px-1.5 text-tiny font-medium text-blue transition-colors duration-fast ease-out hover:bg-blue/5 hover:text-blue-deep active:scale-[0.98] ${FOCUS_RING}`;
  return open ? (
    <span className="block max-w-[22rem] whitespace-normal">
      <ol className="space-y-0.5">
        {rules.map((rule, i) => (
          <li key={`${i}-${rule}`} className="break-words">
            {rule}
          </li>
        ))}
      </ol>
      <button type="button" aria-expanded={true} onClick={() => setOpen(false)} className={`-ml-1.5 mt-0.5 ${toggle}`}>
        Fewer
      </button>
    </span>
  ) : (
    <span className="flex max-w-[16rem] items-center gap-1">
      <span className="min-w-0 truncate">{s.first}</span>
      {s.more > 0 ? (
        <button type="button" aria-expanded={false} aria-label={`Show all ${rules.length} rules`} onClick={() => setOpen(true)} className={toggle}>
          +{s.more}
        </button>
      ) : null}
    </span>
  );
};

/** Rows wear a provisional marker only when the report mixes audit states;
 *  when every row shares one, the status-line chip says it once. */
const columnsFor = (markProvisional: boolean, current: number | undefined): Column<CatalogRow>[] => [
  {
    id: "dataset",
    header: "Dataset",
    className: "w-40 min-w-40 max-w-40 sm:w-52 sm:min-w-52 sm:max-w-52",
    sortValue: (r) => `${r.antibiotic} ${r.species}`,
    firstDir: "asc",
    cell: (r) => (
      <span className="block whitespace-nowrap">
        <span className="flex items-center gap-x-2">
          <span className="truncate">{capitalize(r.antibiotic)}</span>
          {markProvisional && !auditOf(r, current).valid ? (
            <Lozenge tone="coral" fill="solid">
              Provisional
            </Lozenge>
          ) : null}
        </span>
        <span className="block truncate text-tiny font-normal text-muted-text">
          <i>{speciesShort(r.species)}</i> · {r.drug_class}
        </span>
      </span>
    ),
  },
  {
    id: "random",
    header: "Random",
    tip: [{ term: "Random split", text: `Balanced accuracy of ${BEST_LEARNER_WORDS}, stratified 5-fold, with that learner beneath. A learner that leaves any isolate without a call is not a candidate, and neither is the catalog.` }],
    align: "right",
    sortValue: (r) => r.best_random,
    cell: (r) =>
      r.best_random === null || r.best_learner === null ? (
        <span className="block whitespace-nowrap text-right">
          <span className="block">–</span>
          <span className="block text-tiny text-muted-text">
            No eligible learner<span className="sr-only">: {NO_BEST_LEARNER}</span>
          </span>
        </span>
      ) : (
        <Two top={pct(r.best_random)} below={LEARNER_LABEL[r.best_learner]} right />
      ),
  },
  {
    id: "lineage",
    header: "Lineage",
    tip: [{ term: "Lineages held out", text: "The same learner when whole MinHash lineage clusters are held out together, so no test genome has a near-identical relative in training. A missing split names its reason; the inspector has the report's own words." }],
    align: "right",
    sortValue: (r) => r.best_lineage,
    cell: (r) => {
      if (r.best_lineage !== null) return pct(r.best_lineage);
      if (r.best_learner === null) return "–";
      const reason = splitReasonShort(r.lineage_unavailable_reason);
      return reason === null ? <span className="text-muted-text">No split</span> : <Two top="No split" below={reason} right />;
    },
  },
  {
    id: "gap",
    header: "Gap",
    unit: "points",
    tip: [{ term: "Gap", text: "Random minus lineage-held-out balanced accuracy, in points: an estimate, for this sample, of how much of the random-split score depends on near-identical relatives in training." }],
    align: "right",
    sortValue: gap,
    cell: (r) => {
      const g = gap(r);
      return g === null ? "–" : (g * 100).toFixed(1);
    },
  },
  {
    id: "n",
    header: "Isolates",
    tip: [{ term: "Isolates", text: "Genomes with a reported resistant or susceptible call, and how many are resistant." }],
    align: "right",
    sortValue: (r) => r.n,
    cell: (r) => <Two top={int(r.n)} below={`${int(r.n_resistant)} R`} right />,
  },
  {
    id: "vme",
    header: "VME",
    unit: "≤1.5%",
    tip: [
      { term: "Very major error, best learner, random split", text: "Resistant isolates called susceptible, over all resistant isolates." },
      { source: "Filled: meets the bar. Outlined: misses it. FDA-style bar, not a regulatory claim." },
    ],
    align: "right",
    sortValue: (r) => bestLite(r)?.vme ?? null,
    firstDir: "asc",
    status: (r) => {
      const v = bestLite(r)?.vme;
      return v === null || v === undefined ? null : { value: v, bar: BARS.vme };
    },
    cell: () => "–",
  },
  {
    id: "me",
    header: "ME",
    unit: "≤3%",
    tip: [
      { term: "Major error, best learner, random split", text: "Susceptible isolates called resistant, over all susceptible isolates." },
      { source: "Filled: meets the bar. Outlined: misses it. FDA-style bar, not a regulatory claim." },
    ],
    align: "right",
    sortValue: (r) => bestLite(r)?.me ?? null,
    firstDir: "asc",
    status: (r) => {
      const v = bestLite(r)?.me;
      return v === null || v === undefined ? null : { value: v, bar: BARS.me };
    },
    cell: () => "–",
  },
  {
    id: "catalog",
    header: "Catalog",
    tip: [{ term: "Known-gene catalog", text: "Balanced accuracy of the catalog-containment rule alone: AMRFinderPlus genes found by k-mer containment, with a fixed rule per drug. It does not train, so one number serves both splits." }],
    align: "right",
    sortValue: (r) => r.catalog?.balanced_accuracy ?? null,
    cell: (r) => (r.catalog === null ? <span className="text-muted-text">No rule</span> : pct(r.catalog.balanced_accuracy)),
  },
  {
    id: "mic",
    header: "MIC",
    unit: "≥90%",
    tip: [
      { term: "MIC essential agreement", text: "Share of isolates with an exact MIC whose predicted MIC is within one doubling dilution. MIC regression needs at least 100 isolates with an MIC, 50 of them exact." },
      { source: "Filled: meets the 90% bar. Outlined: misses it. A reference, not a performance claim." },
    ],
    align: "right",
    sortValue: (r) => r.mic_essential_agreement,
    status: (r) => (r.mic_essential_agreement === null ? null : { value: r.mic_essential_agreement, bar: BARS.mic }),
    cell: () => <span className="whitespace-nowrap text-muted-text">No MIC</span>,
  },
  {
    id: "rule",
    header: "Rule",
    tip: [{ term: "Set covering machine", text: "The interpretable model's first rule and how many more; the count shows them all, and the inspector places each rule's k-mers." }],
    cell: (r) => <RuleCell rules={r.scm_rules} />,
  },
];

export const DatasetsView = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const rail = useRailState();
  const [f, setF] = useState<Filters>(ARRIVAL);
  const [query, setQuery] = useState("");
  const rows = manifest.datasets;
  const audit = manifestAudit(manifest);

  const q = query.trim().toLowerCase();
  const axes = axesFor(f, q);
  const shown = applyAxes(rows, axes);
  const shownCount = Math.round(useCountUp(shown.length));

  const active =
    (f.species.size > 0 ? 1 : 0) + (f.classes.size > 0 ? 1 : 0) + (f.minN !== "0" ? 1 : 0) + (f.learners.size > 0 ? 1 : 0) + (f.mic ? 1 : 0) + (f.vme ? 1 : 0) + (f.me ? 1 : 0);

  const speciesCounts = countBy(rows, axes, "species", (r) => r.species);
  const classCounts = countBy(rows, axes, "class", (r) => r.drug_class);
  const learnerCounts = countBy(rows, axes, "learner", (r) => r.best_learner ?? "none");
  const withoutMinN = rowsWithout(rows, axes, "minN");
  const speciesNames = [...new Set(rows.map((r) => r.species))].sort((a, b) => (speciesCounts.get(b) ?? 0) - (speciesCounts.get(a) ?? 0) || a.localeCompare(b));
  const classNames = [...new Set(rows.map((r) => r.drug_class))].sort((a, b) => (classCounts.get(b) ?? 0) - (classCounts.get(a) ?? 0) || a.localeCompare(b));
  const clearAll = () => {
    setF(ARRIVAL);
    setQuery("");
  };

  const filterRail = (
    <FilterRail rail={rail} activeCount={active} resultText={`${int(rows.length)} datasets`} onClear={() => setF(ARRIVAL)}>
      <RailGroup label="Species">
        <RailField label="Species" active={f.species.size > 0}>
          <Dropdown
            ariaLabel="Species"
            rest="All species"
            options={speciesNames.map((s) => ({ id: s, label: s, count: speciesCounts.get(s) ?? 0 }))}
            mode={{ kind: "many", value: f.species, onChange: (species) => setF({ ...f, species }) }}
          />
        </RailField>
      </RailGroup>
      <RailGroup label="Drug class">
        <RailField label="Drug class" active={f.classes.size > 0}>
          <Dropdown
            ariaLabel="Drug class"
            rest="All drug classes"
            options={classNames.map((c) => ({ id: c, label: capitalize(c), count: classCounts.get(c) ?? 0 }))}
            mode={{ kind: "many", value: f.classes, onChange: (classes) => setF({ ...f, classes }) }}
          />
        </RailField>
      </RailGroup>
      <RailGroup label="Data">
        <RailField label="Minimum isolates" active={f.minN !== "0"}>
          <Dropdown
            ariaLabel="Minimum isolates"
            options={MIN_N.map((o) => ({ id: o.id, label: o.label, count: o.id === "0" ? undefined : withoutMinN.filter((r) => r.n >= Number(o.id)).length }))}
            mode={{ kind: "one", value: f.minN, onSelect: (minN) => setF({ ...f, minN }) }}
          />
        </RailField>
        <div className="filter-pill-row">
          <CountPill label="Has MIC regression" selected={f.mic} onToggle={() => setF({ ...f, mic: !f.mic })} count={pillCount(rows, axes, "mic", (r) => r.mic_essential_agreement !== null)} />
        </div>
      </RailGroup>
      <RailGroup label="Model">
        <RailField label="Best learner" active={f.learners.size > 0}>
          <Dropdown
            ariaLabel="Best learner"
            rest="All learners"
            options={BEST_LEARNERS.map((l) => ({ id: l, label: LEARNER_NAME[l], count: learnerCounts.get(l) ?? 0, disabled: (learnerCounts.get(l) ?? 0) === 0 && !f.learners.has(l) }))}
            mode={{ kind: "many", value: f.learners, onChange: (learners) => setF({ ...f, learners }) }}
          />
        </RailField>
        <div className="filter-pill-row">
          <CountPill label="Meets VME ≤1.5%" selected={f.vme} onToggle={() => setF({ ...f, vme: !f.vme })} count={pillCount(rows, axes, "vme", (r) => meets(bestLite(r)?.vme, BARS.vme))} />
          <CountPill label="Meets ME ≤3%" selected={f.me} onToggle={() => setF({ ...f, me: !f.me })} count={pillCount(rows, axes, "me", (r) => meets(bestLite(r)?.me, BARS.me))} />
        </div>
      </RailGroup>
    </FilterRail>
  );

  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader title="Datasets" status={<ReportStatus manifest={manifest} />} />

      <div className={`mt-7 ${railGrid(rail)}`}>
        {rail.inline ? filterRail : null}
        <div className="min-w-0">
          <div className="flex flex-wrap items-center gap-3">
            <label className="relative flex min-w-[12rem] flex-1 items-center">
              <span className="sr-only">Search datasets</span>
              <Search size={16} strokeWidth={1.75} aria-hidden="true" className="pointer-events-none absolute left-3.5 text-muted-text" />
              <input
                type="search"
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                placeholder="Search antibiotic, species, drug class, or gene"
                className={`ui-control-target h-10 w-full rounded-full border border-hair bg-paper-pure pl-10 pr-4 text-body text-ink placeholder:text-muted-text transition-colors duration-fast ease-out hover:border-blue/40 ${FOCUS_RING}`}
              />
            </label>
            <FiltersButton rail={rail} activeCount={active} />
            <span className="text-control tabular-nums text-muted-text" aria-live="polite">
              {shownCount === rows.length ? `${int(rows.length)} datasets` : `${int(shownCount)} of ${int(rows.length)} datasets`}
            </span>
          </div>

          <div className="mt-4">
            <DataTable
              label="Datasets"
              rows={shown}
              columns={columnsFor(audit.invalid > 0 && audit.invalid < audit.total, manifest.evaluation_revision)}
              rowKey={(r) => r.id}
              onOpen={(r) => nav.openDataset(r.id)}
              initialSort={{ id: "dataset", dir: "asc" }}
              showAllMatches={q !== ""}
              stickyFirst
              empty={
                rows.length === 0 ? (
                  <span className="text-body text-ink">No dataset has a report yet.</span>
                ) : (
                  <div className="flex flex-wrap items-center gap-3">
                    <span className="text-body text-ink">No dataset matches {q !== "" ? `“${query.trim()}” and ` : ""}these filters.</span>
                    <GhostButton onClick={clearAll}>Clear filters and search</GhostButton>
                  </div>
                )
              }
            />
          </div>
        </div>
      </div>
      {rail.inline ? null : filterRail}
    </div>
  );
};
