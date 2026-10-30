// The dataset catalog: one row per species and antibiotic, in ten columns.
// Every number column has one fixed width and Rule takes the rest, so no
// sort, filter or Show all moves a column; below the table's floor it
// scrolls sideways under a sticky Dataset column. Values are plain text: a
// check marks a value that meets its bar, coral marks a gap over 5 points
// and blue a catalog rule within 2 points of Lineage or ahead, a shape and
// colour before each drug class name mark its family, and a key over the
// table says so in words. One search row carries the search, the
// Filters button below 1440px, and the one result count; below 640px the
// search takes its own row. From 1440px the filters are a rail in the same
// grid. Search, filters, sort and Show all survive a trip away and Back.
// The inspector holds every full value; a clipped rule opens in place.

import { Fragment, useState, type ReactNode } from "react";
import { Search } from "lucide-react";
import type { CatalogRow, Learner, LearnerLite, Manifest } from "../data/schema";
import { capitalize, int, LEARNER_LABEL, LEARNER_NAME, pct, speciesShort } from "../data/format";
import { applyAxes, countBy, memberTest, rowsWithout, type Axis } from "../data/facets";
import { BARS, verdictOf } from "../data/metrics";
import { BEST_LEARNER_WORDS, NO_BEST_LEARNER, ruleSummary, splitReasonShort } from "../data/derive/dataset";
import { CATALOG_CUTOFF, catalogDiff } from "../data/derive/overview";
import { auditOf, manifestAudit } from "../data/derive/evaluation";
import { useKeptState, type Navigator } from "../app/router";
import { ReportStatus } from "../app/ReportStatus";
import { ViewHeader } from "../ui/ViewHeader";
import { Dropdown } from "../ui/Dropdown";
import { DataTable, type Column, type SortState } from "../ui/DataTable";
import { FilterRail, FiltersButton, railGrid, RailCheck, RailChecks, RailControl, RailField, RailGroup, useRailState } from "../ui/FilterRail";
import { FOCUS_RING, GhostButton, MeetsMark, VIEW_SHELL } from "../ui/primitives";
import { useClipped } from "../ui/useClipped";
import { useCountUp } from "../ui/useCountUp";
import { byFamily, ClassGlyph, ClassKey, ClassName } from "../ui/drugClass";

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

/** The gap the Overview's lineage headline counts: more than 5 points of
 *  balanced accuracy lost when whole lineages are held out. */
export const LARGE_GAP = 0.05;

export const largeGap = (r: CatalogRow): boolean => {
  const g = gap(r);
  return g !== null && g > LARGE_GAP;
};

/** The comparison the Overview's catalog headline counts: the known-gene
 *  rule within 2 points of the learner with lineages held out, or ahead. */
export const catalogClose = (r: CatalogRow): boolean => {
  const d = catalogDiff(r);
  return d !== null && d >= -CATALOG_CUTOFF;
};

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

const ARRIVAL_SORT: SortState = { id: "dataset", dir: "asc" };

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

/** A checkbox counts with only its own axis lifted, so its number is the
 *  rows choosing it returns with every other control as it is. */
const checkCount = (rows: readonly CatalogRow[], axes: Axis<CatalogRow>[], lift: string, test: (r: CatalogRow) => boolean): number =>
  applyAxes(rows, axes.filter((a) => a.id !== lift)).filter(test).length;

/** Column widths. Identity narrows on a phone and wraps; each number column
 *  holds its widest value and its header with the sort arrow; Rule takes
 *  the rest, never less than its floor. */
export const DATASET_WIDTHS = {
  dataset: "clamp(8.5rem, 30vw, 13rem)",
  random: "7.75rem",
  lineage: "5rem",
  gap: "4rem",
  n: "5.25rem",
  vme: "4.75rem",
  me: "4.75rem",
  catalog: "5rem",
  mic: "5rem",
  ruleFloor: "12rem",
} as const;

const Muted = ({ children }: { children: ReactNode }) => <span className="text-muted-text">{children}</span>;

/** A value over a quiet second line, both right-aligned. */
const Two = ({ top, below }: { top: ReactNode; below: ReactNode }) => (
  <span className="block text-right">
    <span className="block whitespace-nowrap">{top}</span>
    <span className="block text-tiny font-normal text-muted-text">{below}</span>
  </span>
);

/** A drug name that may break after its "/" rather than be cut. */
const DrugName = ({ name }: { name: string }) => (
  <>
    {capitalize(name)
      .split("/")
      .map((part, i) =>
        i === 0 ? (
          part
        ) : (
          <Fragment key={i}>
            /<wbr />
            {part}
          </Fragment>
        ),
      )}
  </>
);

/** The toggle beside a rule list: "+N" closed, "Fewer" open, in one slot at
 *  the cell's right edge on the first line's baseline. */
const RULE_TOGGLE = `ui-control-target ui-target-24 inline-flex items-center justify-center justify-self-end rounded-sm px-1 text-control font-medium tabular-nums text-blue transition-colors duration-fast ease-out hover:bg-blue/[0.06] hover:text-blue-deep active:scale-[0.98] ${FOCUS_RING}`;

/** The set covering machine's rules. The first rule shows on one line; when
 *  there are more, "+N" opens every rule in place, one per line, at the same
 *  width, and "Fewer" closes them from the same slot. A single rule that is
 *  cut off is itself the button that shows it whole. None of it opens the
 *  row. */
const RuleCell = ({ rules }: { rules: readonly string[] }) => {
  const [open, setOpen] = useState(false);
  // Read only while the rule is cut to one line: shown whole it no longer
  // measures as clipped, and the button must stay to close it.
  const [textRef, clipped] = useClipped<HTMLSpanElement>(open);
  const s = ruleSummary(rules);
  if (s === null) {
    return (
      <Muted>
        –<span className="sr-only">The set covering machine found no rule</span>
      </Muted>
    );
  }
  const line = (
    <span ref={textRef} className={open ? "block whitespace-normal break-words" : "block truncate"}>
      {s.first}
    </span>
  );
  if (s.more === 0) {
    return open || clipped ? (
      <button
        type="button"
        aria-expanded={open}
        onClick={() => setOpen(!open)}
        className={`ui-control-target block w-full rounded-sm text-left transition-colors duration-fast ease-out hover:text-blue ${FOCUS_RING}`}
      >
        {line}
      </button>
    ) : (
      line
    );
  }
  return (
    <div className="grid grid-cols-[minmax(0,1fr)_2.75rem] items-baseline gap-x-2">
      {open ? (
        <ol className="min-w-0 space-y-0.5">
          {rules.map((rule, i) => (
            <li key={`${i}-${rule}`} className="break-words">
              {rule}
            </li>
          ))}
        </ol>
      ) : (
        // A click on the first rule opens the list as well; the count is
        // the way in by keyboard.
        <span data-row-ignore="" className="min-w-0 cursor-pointer" onClick={() => setOpen(true)}>
          {line}
        </span>
      )}
      <button type="button" aria-expanded={open} aria-label={open ? undefined : `Show all ${rules.length} rules`} onClick={() => setOpen(!open)} className={RULE_TOGGLE}>
        {open ? "Fewer" : `+${s.more}`}
      </button>
    </div>
  );
};

/** Rows carry a provisional marker only when the report mixes audit states;
 *  when every row shares one, the status line says it once. Header hints
 *  are optional; the key over the table prints what each mark means. */
export const datasetColumns = (markProvisional: boolean, current: number | undefined): Column<CatalogRow>[] => [
  {
    id: "dataset",
    header: "Dataset",
    width: DATASET_WIDTHS.dataset,
    sortValue: (r) => `${r.antibiotic} ${r.species}`,
    firstDir: "asc",
    cell: (r) => (
      <span className="block">
        <span className="flex items-baseline gap-x-2">
          <span className="min-w-0">
            <DrugName name={r.antibiotic} />
          </span>
          {markProvisional && !auditOf(r, current).valid ? <strong className="shrink-0 text-tiny font-semibold text-coral-text">Provisional</strong> : null}
        </span>
        {/* Species and class wrap to two lines at most; the class, last, is
            what an ellipsis cuts, and the separator stays with the species. */}
        <span className="dataset-sub text-tiny font-normal text-muted-text">
          <span className="whitespace-nowrap">
            <i>{speciesShort(r.species)}</i> ·
          </span>{" "}
          <ClassName drugClass={r.drug_class} />
        </span>
      </span>
    ),
  },
  {
    id: "random",
    header: "Random",
    width: DATASET_WIDTHS.random,
    tip: [{ term: "Random split", text: `Balanced accuracy, the mean of sensitivity and specificity, of ${BEST_LEARNER_WORDS}, stratified 5-fold. That learner is named beneath.` }],
    align: "right",
    sortValue: (r) => r.best_random,
    cell: (r) =>
      r.best_random === null || r.best_learner === null ? (
        <Two
          top={<Muted>–</Muted>}
          below={
            <>
              No eligible learner<span className="sr-only">: {NO_BEST_LEARNER}</span>
            </>
          }
        />
      ) : (
        <Two top={pct(r.best_random)} below={LEARNER_LABEL[r.best_learner]} />
      ),
  },
  {
    id: "lineage",
    header: "Lineage",
    width: DATASET_WIDTHS.lineage,
    tip: [{ term: "Lineages held out", text: "The same learner with whole lineages held out together, so no test genome has a near-identical relative in training. A missing split names its reason." }],
    align: "right",
    sortValue: (r) => r.best_lineage,
    cell: (r) => {
      if (r.best_lineage !== null) return pct(r.best_lineage);
      if (r.best_learner === null) return <Muted>–</Muted>;
      const reason = splitReasonShort(r.lineage_unavailable_reason);
      return reason === null ? <Muted>No split</Muted> : <Two top={<Muted>No split</Muted>} below={reason} />;
    },
  },
  {
    id: "gap",
    header: "Gap",
    unit: "points",
    width: DATASET_WIDTHS.gap,
    tip: [{ term: "Gap", text: "Random minus lineage-held-out balanced accuracy, in points: how much of the random-split score may lean on near-identical relatives in training." }],
    align: "right",
    sortValue: gap,
    cell: (r) => {
      const g = gap(r);
      if (g === null) return <Muted>–</Muted>;
      const points = (g * 100).toFixed(1).replace(/^-/, "−");
      return largeGap(r) ? (
        <span data-cue="gap" className="font-medium text-coral-text">
          {points}
          <span className="sr-only">, more than 5 points lower with lineages held out</span>
        </span>
      ) : (
        points
      );
    },
  },
  {
    id: "n",
    header: "Isolates",
    width: DATASET_WIDTHS.n,
    tip: [{ term: "Isolates", text: "Genomes with a reported resistant or susceptible call, and how many are resistant." }],
    align: "right",
    sortValue: (r) => r.n,
    cell: (r) => <Two top={int(r.n)} below={`${int(r.n_resistant)} R`} />,
  },
  {
    id: "vme",
    header: "VME",
    unit: "≤1.5%",
    width: DATASET_WIDTHS.vme,
    tip: [{ term: "Very major error", text: "Resistant isolates called susceptible, over all resistant isolates; best learner, random split. The bar is FDA-style, not a regulatory claim." }],
    align: "right",
    sortValue: (r) => bestLite(r)?.vme ?? null,
    firstDir: "asc",
    status: (r) => {
      const v = bestLite(r)?.vme;
      return v === null || v === undefined ? null : { value: v, bar: BARS.vme };
    },
    cell: () => <Muted>–</Muted>,
  },
  {
    id: "me",
    header: "ME",
    unit: "≤3%",
    width: DATASET_WIDTHS.me,
    tip: [{ term: "Major error", text: "Susceptible isolates called resistant, over all susceptible isolates; best learner, random split. The bar is FDA-style, not a regulatory claim." }],
    align: "right",
    sortValue: (r) => bestLite(r)?.me ?? null,
    firstDir: "asc",
    status: (r) => {
      const v = bestLite(r)?.me;
      return v === null || v === undefined ? null : { value: v, bar: BARS.me };
    },
    cell: () => <Muted>–</Muted>,
  },
  {
    id: "catalog",
    header: "Catalog",
    width: DATASET_WIDTHS.catalog,
    tip: [{ term: "Known-gene catalog", text: "Balanced accuracy of the catalog rule alone: AMRFinderPlus genes found by k-mer containment, with a fixed rule per drug. It does not train, so one number serves both splits." }],
    align: "right",
    sortValue: (r) => r.catalog?.balanced_accuracy ?? null,
    cell: (r) => {
      if (r.catalog === null) return <Muted>No rule</Muted>;
      const value = pct(r.catalog.balanced_accuracy);
      return catalogClose(r) ? (
        <span data-cue="catalog" className="font-medium text-blue">
          {value}
          <span className="sr-only">, within 2 points of the lineage score or ahead</span>
        </span>
      ) : (
        value
      );
    },
  },
  {
    id: "mic",
    header: "MIC",
    unit: "≥90%",
    width: DATASET_WIDTHS.mic,
    tip: [{ term: "MIC essential agreement", text: "Share of isolates with an exact MIC whose predicted MIC is within one doubling dilution. A model needs at least 100 isolates with an MIC, 50 of them exact." }],
    align: "right",
    sortValue: (r) => r.mic_essential_agreement,
    status: (r) => (r.mic_essential_agreement === null ? null : { value: r.mic_essential_agreement, bar: BARS.mic }),
    cell: () => <Muted>No MIC</Muted>,
  },
  {
    id: "rule",
    header: "Rule",
    unit: "set covering machine",
    minWidth: DATASET_WIDTHS.ruleFloor,
    tip: [{ term: "Set covering machine", text: "The interpretable model's first rule and how many more. A clipped rule opens in place; the inspector places each rule's k-mers." }],
    cell: (r) => <RuleCell rules={r.scm_rules} />,
  },
];

/** What the table's marks mean, printed once over it. Each mark is shown
 *  as it appears in the cells, then said in words. The drug-class key
 *  always takes the second line, so the key is two lines at every width
 *  the rail docks and the rail's first group stays on the header line. */
const TableKey = () => (
  <p className="mt-3 flex flex-wrap items-baseline gap-x-5 gap-y-1 text-control text-muted-text">
    <span className="inline-flex items-baseline gap-1">
      <span className="inline-flex items-baseline text-live-text">
        <MeetsMark />
      </span>
      meets the bar under VME, ME or MIC
    </span>
    <span>
      <span className="font-medium text-coral-text">Gap</span> over 5 points
    </span>
    <span>
      <span className="font-medium text-blue">Catalog</span> within 2 points of Lineage, or ahead
    </span>
    <ClassKey className="basis-full" />
  </p>
);

export const DatasetsView = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const rail = useRailState();
  const [f, setF] = useKeptState<Filters>("datasets.filters", ARRIVAL);
  const [query, setQuery] = useKeptState("datasets.query", "");
  const [sort, setSort] = useKeptState<SortState>("datasets.sort", ARRIVAL_SORT);
  const [all, setAll] = useKeptState("datasets.all", false);
  const rows = manifest.datasets;
  const audit = manifestAudit(manifest);

  const q = query.trim().toLowerCase();
  const axes = axesFor(f, q);
  const shown = applyAxes(rows, axes);
  const shownCount = Math.round(useCountUp(shown.length));
  const countText = (n: number) => (n === rows.length ? `${int(rows.length)} datasets` : `${int(n)} of ${int(rows.length)} datasets`);

  const active =
    (f.species.size > 0 ? 1 : 0) + (f.classes.size > 0 ? 1 : 0) + (f.minN !== "0" ? 1 : 0) + (f.learners.size > 0 ? 1 : 0) + (f.mic ? 1 : 0) + (f.vme ? 1 : 0) + (f.me ? 1 : 0);

  const speciesCounts = countBy(rows, axes, "species", (r) => r.species);
  const classCounts = countBy(rows, axes, "class", (r) => r.drug_class);
  const learnerCounts = countBy(rows, axes, "learner", (r) => r.best_learner ?? "none");
  const withoutMinN = rowsWithout(rows, axes, "minN");
  const speciesNames = [...new Set(rows.map((r) => r.species))].sort((a, b) => (speciesCounts.get(b) ?? 0) - (speciesCounts.get(a) ?? 0) || a.localeCompare(b));
  // Classes list by family, so each family's mark runs in one block.
  const classNames = byFamily([...new Set(rows.map((r) => r.drug_class))], (c) => classCounts.get(c) ?? 0);
  const clearAll = () => {
    setF(ARRIVAL);
    setQuery("");
  };

  // One label per control: a group holding a single control is named by
  // its group label alone.
  const filterRail = (
    <FilterRail rail={rail} activeCount={active} resultText={`${int(rows.length)} datasets`} onClear={() => setF(ARRIVAL)}>
      <RailGroup label="Species" active={f.species.size > 0}>
        <RailControl>
          <Dropdown
            ariaLabel="Species"
            rest="All species"
            options={speciesNames.map((s) => ({ id: s, label: s, count: speciesCounts.get(s) ?? 0 }))}
            mode={{ kind: "many", value: f.species, onChange: (species) => setF({ ...f, species }) }}
          />
        </RailControl>
      </RailGroup>
      <RailGroup label="Drug class" active={f.classes.size > 0}>
        <RailControl>
          <Dropdown
            ariaLabel="Drug class"
            rest="All drug classes"
            options={classNames.map((c) => ({ id: c, label: capitalize(c), count: classCounts.get(c) ?? 0, icon: <ClassGlyph drugClass={c} /> }))}
            mode={{ kind: "many", value: f.classes, onChange: (classes) => setF({ ...f, classes }) }}
          />
        </RailControl>
      </RailGroup>
      <RailGroup label="Data" active={f.minN !== "0" || f.mic}>
        <RailField label="Minimum isolates" active={f.minN !== "0"}>
          <Dropdown
            ariaLabel="Minimum isolates"
            options={MIN_N.map((o) => ({ id: o.id, label: o.label, count: o.id === "0" ? undefined : withoutMinN.filter((r) => r.n >= Number(o.id)).length }))}
            mode={{ kind: "one", value: f.minN, onSelect: (minN) => setF({ ...f, minN }) }}
          />
        </RailField>
        <RailCheck label="Has MIC regression" checked={f.mic} onChange={(mic) => setF({ ...f, mic })} count={checkCount(rows, axes, "mic", (r) => r.mic_essential_agreement !== null)} />
      </RailGroup>
      <RailGroup label="Best learner" active={f.learners.size > 0}>
        <RailControl>
          <Dropdown
            ariaLabel="Best learner"
            rest="All learners"
            options={BEST_LEARNERS.map((l) => ({ id: l, label: LEARNER_NAME[l], count: learnerCounts.get(l) ?? 0, disabled: (learnerCounts.get(l) ?? 0) === 0 && !f.learners.has(l) }))}
            mode={{ kind: "many", value: f.learners, onChange: (learners) => setF({ ...f, learners }) }}
          />
        </RailControl>
      </RailGroup>
      <RailChecks label="Meets the bar" active={f.vme || f.me}>
        <RailCheck label="VME ≤1.5%" checked={f.vme} onChange={(vme) => setF({ ...f, vme })} count={checkCount(rows, axes, "vme", (r) => meets(bestLite(r)?.vme, BARS.vme))} />
        <RailCheck label="ME ≤3%" checked={f.me} onChange={(me) => setF({ ...f, me })} count={checkCount(rows, axes, "me", (r) => meets(bestLite(r)?.me, BARS.me))} />
      </RailChecks>
    </FilterRail>
  );

  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader title="Datasets" status={<ReportStatus manifest={manifest} />} />

      <div className={`mt-7 ${railGrid(rail)}`}>
        {rail.inline ? filterRail : null}
        <div className="min-w-0">
          {/* Below 640px the search takes the whole row; Filters and the
              count wrap under it. */}
          <div className="flex flex-wrap items-center gap-x-3 gap-y-2">
            <label className="relative flex w-full min-w-0 items-center sm:w-auto sm:min-w-[12rem] sm:flex-1">
              <span className="sr-only">Search datasets</span>
              <Search size={16} strokeWidth={1.75} aria-hidden="true" className="pointer-events-none absolute left-3.5 text-muted-text" />
              <input
                type="search"
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                placeholder="Search antibiotic, species, drug class, or gene"
                className={`ui-control-target h-10 w-full rounded-md border border-hair bg-paper-pure pl-10 pr-4 text-body text-ink placeholder:text-muted-text transition-colors duration-fast ease-out hover:border-blue/40 ${FOCUS_RING}`}
              />
            </label>
            <FiltersButton rail={rail} activeCount={active} />
            <span aria-hidden="true" className="text-control tabular-nums text-muted-text">
              {countText(shownCount)}
            </span>
            <span className="sr-only" aria-live="polite">
              {countText(shown.length)}
            </span>
          </div>

          {shown.length > 0 ? <TableKey /> : null}

          <div className="mt-3">
            <DataTable
              label="Datasets"
              rows={shown}
              columns={datasetColumns(audit.invalid > 0 && audit.invalid < audit.total, manifest.evaluation_revision)}
              rowKey={(r) => r.id}
              onOpen={(r) => nav.openDataset(r.id)}
              initialSort={sort}
              onSortChange={setSort}
              initialShowAll={all}
              onShowAllChange={setAll}
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
