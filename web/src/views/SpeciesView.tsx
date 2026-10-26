// One species: how its analyzed genomes cluster into lineages, and how
// resistance follows them. The title is the species; choosing a species and
// a map drug is the page's primary action, so the selectors take the row
// under the status line, followed by one coverage statement, the map drug's
// result and a strip of corpus facts. Then the population structure (lineage
// blocks by default), cross-resistance, and the species' antibiotic
// datasets as one full-width table. The genome map lives in the dataset
// inspector.

import { useMemo, type ReactNode } from "react";
import type { Manifest, SpeciesReport, SpeciesRow } from "../data/schema";
import { capitalize, distance, int, pct } from "../data/format";
import { loadErrorText, useCorpus, useSpecies } from "../data/load";
import { BEST_LEARNER_WORDS, NO_BEST_LEARNER } from "../data/derive/dataset";
import { useKeptState, type Navigator } from "../app/router";
import { ReportStatus } from "../app/ReportStatus";
import { ViewHeader } from "../ui/ViewHeader";
import { Dropdown } from "../ui/Dropdown";
import { FigureBlock, PanelState, Section, Segmented, VIEW_SHELL } from "../ui/primitives";
import { DataTable } from "../ui/DataTable";
import { MIN_GENOME_PX, StructureExplorer, structureModel, type StructurePick, type StructureView } from "../charts/Structure";
import { CrossPanel } from "../charts/CrossMatrix";
import { MIN_JOINT, crossPairs, pairSummary } from "../charts/crossPairs";
import { BLOCK_MIN_GENOMES, POOL_HEIGHT_FACTOR, drugTotals, drugsByCalls, resistantLineages, rsOf } from "../charts/lineageBlocks";

/** What the reader set on this page for the species shown. */
interface SpeciesUi {
  slug: string;
  drug: number | null;
  view: StructureView;
  pick: StructurePick;
  pair: string | null;
}

/** A selector with its ink label on the same line. */
const Field = ({ label, children }: { label: string; children: ReactNode }) => (
  <div className="flex min-w-0 items-center gap-2">
    <span className="text-tiny font-semibold text-ink">{label}</span>
    {children}
  </div>
);

/** The map drug's result on the page's long-standing footing: resistant of
 *  the genomes called R or S, with intermediate calls counted beside it,
 *  never folded in, and how much of the resistance the three largest
 *  holders carry. */
export const drugResult = (s: SpeciesReport, drug: number): string => {
  const t = drugTotals(s, drug);
  const name = s.drugs[drug] ?? "";
  const rs = rsOf(t);
  const more = (lead: string) => (t.i > 0 ? `${lead}${int(t.i)} more ${t.i === 1 ? "is" : "are"} intermediate` : "");
  if (rs === 0) return `No genome in this species file is called R or S for ${name}${t.i > 0 ? `; ${int(t.i)} ${t.i === 1 ? "is" : "are"} intermediate` : ""}.`;
  if (t.r === 0) return `None of the ${int(rs)} genomes called R or S for ${name} is resistant${more("; ")}.`;
  const lineages = resistantLineages(s, drug);
  const held = lineages.slice(0, 3).reduce((a, l) => a + l.resistant, 0);
  const them = t.i > 0 ? "the resistant ones" : "them";
  const holders =
    lineages.length === 1 ? `one lineage holds all of ${them}` : lineages.length === 2 ? `two lineages hold all of ${them}` : `three lineages hold ${pct(held / t.r)} of ${them}`;
  return `${int(t.r)} of ${int(rs)} genomes called R or S for ${name} are resistant${more(", and ")}; ${holders}.`;
};

/** One coverage statement from the manifest, and the two facts that change
 *  how every count on the page reads. The count with unconflicted reported
 *  R/S calls, the pool the analysis draws from, is not in the manifest, so
 *  the statement never offers the any-label count as that pool. Profiles,
 *  cross-resistance and the corpus call counts come from each phenotype's
 *  consolidated label in src/report.rs, MIC-inferred labels included; drop
 *  that sentence when they use reported calls only. */
const Coverage = ({ row }: { row: SpeciesRow }) => (
  <div className="mt-5 space-y-1">
    <p className="text-body text-ink">
      {`${int(row.genomes_analyzed)} genomes analyzed (${int(row.genomes_labeled)} carry any label, including MIC-inferred calls); ${int(row.datasets_reported)} of ${int(row.datasets_curated)} curated datasets have a report.`}
    </p>
    <p className="type-note max-w-[60rem]">
      Analysis genomes are selected per species, up to a cap, from those with unconflicted reported R/S calls; QC and assembly availability can reduce a cohort. Profiles, cross-resistance and call counts on this page include labels inferred from MIC, which dataset reports leave out.
    </p>
  </div>
);

const SpeciesBody = ({
  s,
  manifest,
  nav,
  order,
  drug,
  ui,
  update,
}: {
  s: SpeciesReport;
  manifest: Manifest;
  nav: Navigator;
  order: readonly number[];
  drug: number;
  ui: SpeciesUi;
  update(patch: Partial<SpeciesUi>): void;
}) => {
  const corpus = useCorpus();
  const model = useMemo(() => structureModel(s), [s]);
  const pairs = useMemo(() => crossPairs(s), [s]);
  const summary = useMemo(() => pairSummary(pairs), [pairs]);
  const drugName = s.drugs[drug] ?? "";
  const corpusRow = corpus.status === "ready" ? corpus.data.by_species.find((r) => r.species === s.species) : undefined;
  const datasets = manifest.datasets.filter((d) => d.species_slug === s.slug);
  const top = summary.highest;
  const facts = [
    { label: "Resistant calls", value: corpusRow === undefined ? "–" : int(corpusRow.pairs_resistant) },
    { label: "Susceptible calls", value: corpusRow === undefined ? "–" : int(corpusRow.pairs_susceptible) },
    { label: "Lineages", value: int(s.n_lineages) },
    { label: "Single-genome lineages", value: int(model.singletons) },
    { label: "Catalog gene groups seen", value: int(s.catalog_groups.length) },
  ];
  const caption =
    ui.view === "genomes"
      ? `Each column is one genome in tree order, at least ${MIN_GENOME_PX}px wide; the minimap shows all ${int(model.n)} genomes by their ${drugName} call, with the open window outlined, and the bands under the tree mark lineages and blocks. Rows are drugs, most R or S calls first, with ${drugName} outlined.`
      : ui.view === "table"
        ? `The tables hold the figure's counts: ${drugName} by lineage block, then every genome's call for every drug. A row pins its block or genome in the inspector.`
        : `Each column is a lineage block: one lineage of at least ${BLOCK_MIN_GENOMES} genomes, or neighbouring smaller lineages pooled while each joining merge stays within ${POOL_HEIGHT_FACTOR} times the lineage threshold (Mash distance ${distance(POOL_HEIGHT_FACTOR * s.lineage_threshold)}). Lineages are clusters at Mash distance ${distance(s.lineage_threshold)} in the average-linkage tree of Mash distances, drawn above the blocks. Rows are drugs, most R or S calls first, with ${drugName} outlined; beside each drug are its genomes called R or S, and those called I where the species has any. A cell splits its block's genomes into resistant, intermediate, susceptible and not tested.`;

  return (
    <div className="mt-6 space-y-10">
      <div>
        <p className="type-headline">{drugResult(s, drug)}</p>
        <section aria-labelledby="sp-facts-title" className="mt-4 flex flex-wrap items-baseline gap-x-8 gap-y-2">
          <h2 id="sp-facts-title" className="text-tiny font-semibold text-muted-text">
            Corpus facts
          </h2>
          <dl className="m-0 flex flex-wrap items-baseline gap-x-8 gap-y-2">
            {facts.map((f) => (
              <div key={f.label} className="flex items-baseline gap-2">
                <dt className="text-control text-muted-text">{f.label}</dt>
                <dd className="m-0 text-control font-semibold tabular-nums text-ink">{f.value}</dd>
              </div>
            ))}
          </dl>
        </section>
      </div>

      {s.genomes.length === 0 ? (
        <PanelState kind="empty" title="This species file lists no genomes" />
      ) : (
        <FigureBlock
          id="sp-structure"
          title="Lineages and resistance profiles"
          actions={
            <Segmented<StructureView>
              ariaLabel="Show population structure as"
              options={[
                { id: "blocks", label: "Blocks" },
                { id: "genomes", label: "Genomes" },
                { id: "table", label: "Table" },
              ]}
              value={ui.view}
              onSelect={(view) => update({ view })}
            />
          }
          figure={
            <StructureExplorer
              key={s.slug}
              report={s}
              model={model}
              order={order}
              drug={drug}
              view={ui.view}
              pick={ui.pick}
              onPick={(pick) => update({ pick })}
              onView={(view) => update({ view })}
              caption={caption}
            />
          }
          caption={null}
        />
      )}

      <FigureBlock
        id="sp-cross"
        title="Cross-resistance"
        result={
          summary.total === 0
            ? undefined
            : `${int(summary.evaluated)} of ${int(summary.total)} pairs have at least ${MIN_JOINT} joint R/S calls${summary.notEvaluated > 0 ? `; ${int(summary.notEvaluated)} are not evaluated.` : "."}${
                top === null ? "" : ` ${capitalize(top.a)}–${top.b} has highest phi: ${top.phi!.toFixed(2)} over ${int(top.n)} genomes called R or S for both.`
              }`
        }
        figure={
          <CrossPanel
            key={s.slug}
            report={s}
            order={order}
            pairs={pairs}
            summary={summary}
            selected={ui.pair}
            onSelect={(pair) => update({ pair })}
            caption={`Phi is the correlation between two drugs' resistant calls over the genomes called R or S for both; intermediate and untested genomes never count, and the calls include labels inferred from MIC. A pair is evaluated only with at least ${MIN_JOINT} such genomes, and its phi is undefined when one drug has no resistant or no susceptible genome among them.`}
          />
        }
        caption={null}
      />

      <Section id="sp-datasets" title="Antibiotic datasets">
        <DataTable
          label={`${s.species} datasets`}
          rows={datasets}
          rowKey={(d) => d.id}
          onOpen={(d) => nav.openDataset(d.id)}
          initialSort={{ id: "n", dir: "desc" }}
          columns={[
            { id: "antibiotic", header: "Antibiotic", sortValue: (d) => d.antibiotic, firstDir: "asc", cell: (d) => capitalize(d.antibiotic) },
            { id: "class", header: "Drug class", sortValue: (d) => d.drug_class, firstDir: "asc", cell: (d) => d.drug_class },
            { id: "n", header: "Isolates", align: "right", sortValue: (d) => d.n, cell: (d) => int(d.n) },
            { id: "resistant", header: "Resistant", align: "right", sortValue: (d) => d.n_resistant, cell: (d) => int(d.n_resistant) },
            {
              id: "random",
              header: "Random",
              tip: [{ term: "Random split", text: `Balanced accuracy on random splits of ${BEST_LEARNER_WORDS}.` }],
              align: "right",
              sortValue: (d) => d.best_random,
              cell: (d) => (d.best_random === null ? NO_BEST_LEARNER : pct(d.best_random)),
            },
            { id: "lineage", header: "Lineage", align: "right", sortValue: (d) => d.best_lineage, cell: (d) => pct(d.best_lineage) },
          ]}
          empty={<p className="type-note">No dataset for this species has a report.</p>}
        />
      </Section>
    </div>
  );
};

export const SpeciesView = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const options = manifest.species.filter((s) => s.genomes_analyzed > 0);
  const slug = nav.route.species ?? options[0]?.slug ?? null;
  const species = useSpecies(slug);
  const [kept, setKept] = useKeptState<SpeciesUi | null>("species.ui", null);
  const report = species.status === "ready" ? species.data : null;
  const order = useMemo(() => (report === null ? [] : drugsByCalls(report)), [report]);
  const row = manifest.species.find((r) => r.slug === slug);
  const ui: SpeciesUi | null =
    report === null ? null : kept !== null && kept.slug === report.slug ? kept : { slug: report.slug, drug: null, view: "blocks", pick: null, pair: null };
  const drug = report === null ? 0 : (ui?.drug ?? order[0] ?? 0);
  const update = (patch: Partial<SpeciesUi>) => {
    if (report === null) return;
    const fresh: SpeciesUi = { slug: report.slug, drug: null, view: "blocks", pick: null, pair: null };
    setKept((k) => ({ ...(k !== null && k.slug === report.slug ? k : fresh), ...patch }));
  };

  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader title={report === null ? "Species" : report.species} purpose="How genomes cluster into lineages, and how resistance follows them." status={<ReportStatus manifest={manifest} />} />
      <div className="mt-6 flex min-h-10 flex-wrap items-center gap-x-6 gap-y-3">
        <Field label="Species">
          <Dropdown
            ariaLabel="Species"
            options={manifest.species.map((s) => ({ id: s.slug, label: s.species, disabled: s.genomes_analyzed === 0 }))}
            mode={{ kind: "one", value: slug ?? "", onSelect: (v) => nav.go({ view: "species", species: v }) }}
          />
        </Field>
        {report !== null ? (
          <Field label="Map drug">
            <Dropdown
              ariaLabel="Map drug"
              options={order.map((j) => ({ id: String(j), label: capitalize(report.drugs[j] ?? "") }))}
              mode={{ kind: "one", value: String(drug), onSelect: (v) => update({ drug: Number(v) }) }}
            />
          </Field>
        ) : null}
      </div>
      {row !== undefined && row.genomes_analyzed > 0 ? <Coverage row={row} /> : null}
      {slug === null ? (
        <div className="mt-8">
          <PanelState kind="empty" title="No species has a structure file yet" message="Species files appear once `grm structure` and `grm report` have run." />
        </div>
      ) : species.status === "loading" ? (
        <div className="mt-8">
          <PanelState kind="loading" title="Loading the species file" />
        </div>
      ) : species.status === "error" ? (
        <div className="mt-8">
          <PanelState kind="error" title="This species file could not be read" message={loadErrorText(species.error)} />
        </div>
      ) : report !== null && ui !== null ? (
        <SpeciesBody s={report} manifest={manifest} nav={nav} order={order} drug={drug} ui={ui} update={update} />
      ) : null}
    </div>
  );
};
