// Population structure for one species: lineages and the resistance profile
// of every genome in dendrogram order, cross-resistance, and the data behind
// them. The title is the species; choosing a species and a drug is the
// page's primary action, so the selectors take the row under the status
// line. The genome map lives in the dataset inspector.

import { useState, type ReactNode } from "react";
import type { Manifest, SpeciesReport } from "../data/schema";
import { capitalize, distance, int, pct, plural, speciesShort } from "../data/format";
import { loadErrorText, useCorpus, useSpecies } from "../data/load";
import { BEST_LEARNER_WORDS, NO_BEST_LEARNER } from "../data/derive/dataset";
import type { Navigator } from "../app/router";
import { ReportStatus } from "../app/ReportStatus";
import { ViewHeader } from "../ui/ViewHeader";
import { Dropdown } from "../ui/Dropdown";
import { FactRows, FigureBlock, LegendItem, PanelState, Section, Segmented, VIEW_SHELL } from "../ui/primitives";
import { DataTable } from "../ui/DataTable";
import { StructureFigure } from "../charts/Structure";
import { CrossMatrix } from "../charts/CrossMatrix";

const tally = (s: SpeciesReport, j: number) => {
  let r = 0;
  let tested = 0;
  s.genomes.forEach((g) => {
    const c = g.profile[j];
    if (c === "R") r += 1;
    if (c === "R" || c === "S") tested += 1;
  });
  return { r, tested };
};

/** Drugs with the most resistant-or-susceptible calls first. */
const drugOrder = (s: SpeciesReport): number[] =>
  s.drugs.map((_, j) => j).sort((a, b) => tally(s, b).tested - tally(s, a).tested || (s.drugs[a] ?? "").localeCompare(s.drugs[b] ?? ""));

const lineageStats = (s: SpeciesReport) => {
  const sizes = new Map<number, number>();
  s.genomes.forEach((g) => sizes.set(g.lineage, (sizes.get(g.lineage) ?? 0) + 1));
  const largest = Math.max(0, ...sizes.values());
  return { count: sizes.size, largestShare: s.genomes.length === 0 ? 0 : largest / s.genomes.length };
};

const concentration = (s: SpeciesReport, j: number) => {
  const by = new Map<number, number>();
  let total = 0;
  s.genomes.forEach((g) => {
    if (g.profile[j] === "R") {
      by.set(g.lineage, (by.get(g.lineage) ?? 0) + 1);
      total += 1;
    }
  });
  const top = [...by.values()].sort((a, b) => b - a).slice(0, 3).reduce((a, b) => a + b, 0);
  return { lineages: by.size, share: total === 0 ? 0 : top / total };
};

/** A selector with its ink label on the same line. */
const Field = ({ label, children }: { label: string; children: ReactNode }) => (
  <div className="flex min-w-0 items-center gap-2">
    <span className="text-tiny font-semibold text-ink">{label}</span>
    {children}
  </div>
);

const Swatch = ({ className }: { className: string }) => <span className={`inline-block h-2.5 w-3.5 rounded-sm ${className}`} />;

// Profiles, cross-resistance, and the species and corpus call counts come from
// each phenotype's label in src/report.rs, MIC-inferred labels included. The
// notes below say so; drop them when those counts use reported calls only.
const SpeciesBody = ({ s, manifest, nav, drug }: { s: SpeciesReport; manifest: Manifest; nav: Navigator; drug: number }) => {
  const corpus = useCorpus();
  const [crossView, setCrossView] = useState<"matrix" | "table">("matrix");
  const order = drugOrder(s);
  const lin = lineageStats(s);
  const t = tally(s, drug);
  const conc = concentration(s, drug);
  const drugName = s.drugs[drug] ?? "";
  const strongest = [...s.cross].sort((a, b) => b.phi - a.phi)[0];
  const row = manifest.species.find((r) => r.slug === s.slug);
  const corpusRow = corpus.status === "ready" ? corpus.data.by_species.find((r) => r.species === s.species) : undefined;
  const datasets = manifest.datasets.filter((d) => d.species_slug === s.slug);

  return (
    <div className="mt-8 space-y-10">
      <FigureBlock
        id="sp-structure"
        title="Lineages and resistance profiles"
        result={
          t.r === 0
            ? `No genome is resistant to ${drugName}.`
            : `${int(t.r)} of ${int(t.tested)} genomes tested for ${drugName} are resistant; three lineages hold ${pct(conc.share)} of them.`
        }
        figure={
          <>
            <div className="mb-3 flex flex-wrap gap-x-5 gap-y-1">
              <LegendItem swatch={<Swatch className="bg-[color:var(--series-r)]" />} name="Resistant" />
              <LegendItem swatch={<Swatch className="bg-[color:var(--series-s)]" />} name="Susceptible" />
              <LegendItem swatch={<Swatch className="bg-[color:var(--series-neutral)]" />} name="Intermediate" />
              <LegendItem swatch={<Swatch className="bg-surface-control" />} name="Not tested" />
            </div>
            <StructureFigure report={s} drugOrder={order} selected={drug} />
          </>
        }
        caption={`Each column is one of ${int(s.genomes.length)} genomes, ordered by average-linkage clustering of Mash distances; each row is one drug, most-tested first, with ${drugName} outlined. The band marks the ${plural(lin.count, "lineage")} at Mash distance ${distance(s.lineage_threshold)}; the largest holds ${pct(lin.largestShare)} of genomes. Profiles, cross-resistance, and the call counts below include labels inferred from MIC; dataset reports leave those out.`}
      />

      <FigureBlock
        id="sp-cross"
        title="Cross-resistance"
        result={
          strongest === undefined
            ? "No drug pair has 30 or more genomes tested for both."
            : `Resistance to ${strongest.a} and ${strongest.b} goes together most strongly: phi ${strongest.phi.toFixed(2)} over ${int(strongest.n)} genomes tested for both.`
        }
        actions={
          s.cross.length > 0 ? (
            <Segmented<"matrix" | "table">
              ariaLabel="Show cross-resistance as"
              options={[
                { id: "matrix", label: "Matrix" },
                { id: "table", label: "Table" },
              ]}
              value={crossView}
              onSelect={setCrossView}
            />
          ) : undefined
        }
        figure={
          s.cross.length === 0 ? null : crossView === "table" ? (
            <DataTable
              label="Cross-resistance pairs"
              rows={s.cross}
              rowKey={(p) => `${p.a}|${p.b}`}
              initialSort={{ id: "phi", dir: "desc" }}
              columns={[
                { id: "pair", header: "Drug pair", sortValue: (p) => `${p.a} ${p.b}`, firstDir: "asc", cell: (p) => `${capitalize(p.a)} and ${p.b}` },
                { id: "n", header: "Tested", align: "right", sortValue: (p) => p.n, cell: (p) => int(p.n) },
                { id: "rr", header: "Both R", align: "right", sortValue: (p) => p.rr, cell: (p) => int(p.rr) },
                { id: "split", header: "One R", align: "right", sortValue: (p) => p.rs + p.sr, cell: (p) => int(p.rs + p.sr) },
                { id: "ss", header: "Both S", align: "right", sortValue: (p) => p.ss, cell: (p) => int(p.ss) },
                { id: "phi", header: "Phi", align: "right", sortValue: (p) => p.phi, cell: (p) => p.phi.toFixed(2) },
              ]}
              empty={null}
            />
          ) : (
            <CrossMatrix report={s} order={order} />
          )
        }
        caption="Phi is the correlation between two drugs' resistant calls over genomes tested for both; pairs with fewer than 30 such genomes are blank."
      />

      <Section id="sp-facts" title="Corpus facts">
        <div className="grid grid-cols-1 gap-x-12 gap-y-6 lg:grid-cols-2">
          <FactRows
            rows={[
              { label: "Genomes with a call", value: row === undefined ? "–" : int(row.genomes_labeled) },
              { label: "Datasets", value: row === undefined ? "–" : `${int(row.datasets_reported)} reported of ${int(row.datasets_curated)} curated` },
              { label: "Resistant calls", value: corpusRow === undefined ? "–" : int(corpusRow.pairs_resistant) },
              { label: "Susceptible calls", value: corpusRow === undefined ? "–" : int(corpusRow.pairs_susceptible) },
              { label: "Catalog gene groups seen", value: int(s.catalog_groups.length) },
            ]}
          />
          <DataTable
            label={`${s.species} datasets`}
            rows={datasets}
            rowKey={(d) => d.id}
            onOpen={(d) => nav.openDataset(d.id)}
            initialSort={{ id: "n", dir: "desc" }}
            columns={[
              { id: "antibiotic", header: "Antibiotic", sortValue: (d) => d.antibiotic, firstDir: "asc", cell: (d) => capitalize(d.antibiotic) },
              { id: "n", header: "Isolates", align: "right", sortValue: (d) => d.n, cell: (d) => int(d.n) },
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
        </div>
      </Section>
    </div>
  );
};

export const SpeciesView = ({ manifest, nav }: { manifest: Manifest; nav: Navigator }) => {
  const options = manifest.species.filter((s) => s.genomes_analyzed > 0);
  const slug = nav.route.species ?? options[0]?.slug ?? null;
  const species = useSpecies(slug);
  const [drugChoice, setDrugChoice] = useState<{ slug: string; drug: number } | null>(null);
  const report = species.status === "ready" ? species.data : null;
  const drug = report === null ? 0 : drugChoice?.slug === report.slug ? drugChoice.drug : (drugOrder(report)[0] ?? 0);

  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader
        title={report === null ? "Species" : report.species}
        purpose="How genomes cluster into lineages, and how resistance follows them."
        status={
          <ReportStatus manifest={manifest}>
            {report === null ? undefined : <span>{`${speciesShort(report.species)}: ${int(report.genomes.length)} genomes, ${plural(report.drugs.length, "drug")}`}</span>}
          </ReportStatus>
        }
      />
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
              options={drugOrder(report).map((j) => ({ id: String(j), label: capitalize(report.drugs[j] ?? "") }))}
              mode={{ kind: "one", value: String(drug), onSelect: (v) => setDrugChoice({ slug: report.slug, drug: Number(v) }) }}
            />
          </Field>
        ) : null}
      </div>
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
      ) : report !== null ? (
        <SpeciesBody s={report} manifest={manifest} nav={nav} drug={drug} />
      ) : null}
    </div>
  );
};
