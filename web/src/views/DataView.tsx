// Data: where GRM's genomes and lab calls come from, what passes curation,
// what is set aside, and which rules decide inclusion. Four linear sections
// put the numbers first (the curation funnel and isolate sources), then
// species coverage, then the label and inclusion rules. Every count comes
// from the manifest or the corpus file at render time.

import { useEffect, useState, type MouseEvent, type ReactNode } from "react";
import type { Corpus, Manifest } from "../data/schema";
import { int, pct } from "../data/format";
import { loadErrorText, useCorpus } from "../data/load";
import { classMinimum } from "../data/derive/overview";
import type { Navigator } from "../app/router";
import { PanelState, ScrollFrame, VIEW_SHELL } from "../ui/primitives";
import { ViewHeader } from "../ui/ViewHeader";
import { ReportStatus } from "../app/ReportStatus";
import { YearColumns } from "../charts/YearColumns";
import { BolsterStatStrip, type BolsterStat } from "../charts/amr";
import { minimumSentence, SectionHeading, selectionCap, showSection, SpeciesCounts, STATUS_LINK } from "./findings";
import { ClindamycinMethods } from "./inspector/ClindamycinNote";

/** The anchors other views link to on Data. */
export const DATA_ANCHORS = { funnel: "m-funnel", sources: "m-sources", species: "m-species", labels: "m-labels" } as const;

const P = ({ children }: { children: ReactNode }) => <p className="prose-measure text-body text-ink">{children}</p>;

const Detail = ({ children }: { children: ReactNode }) => <p className="prose-measure text-note text-ink">{children}</p>;

const DataSection = ({ id, n, title, children }: { id: string; n: number; title: string; children: ReactNode }) => (
  <section id={id} aria-labelledby={`${id}-title`} className="scroll-mt-6 border-t border-hair pt-7">
    <SectionHeading id={id} n={n}>
      {title}
    </SectionHeading>
    <div className="mt-4 space-y-4">{children}</div>
  </section>
);

/** A count's length on a track of the whole: blue for what was kept, hatched
 *  blue for a distinct kind kept in records but excluded from model training,
 *  and grey for what was set aside. */
const Bar = ({ value, max, tone = "kept" }: { value: number; max: number; tone?: "kept" | "dropped" | "excluded" }) => (
  <span className="relative block h-2.5 w-full min-w-[6rem] overflow-hidden rounded-sm bg-surface-subtle" aria-hidden="true">
    <span
      className={`absolute inset-y-0 left-0 rounded-r-sm ${tone === "kept" ? "bg-[color:var(--series-s)]" : tone === "dropped" ? "bg-[color:var(--series-neutral)]" : ""}`}
      style={{
        width: `${max === 0 ? 0 : Math.max(0.5, (value / max) * 100)}%`,
        ...(tone === "excluded"
          ? {
              backgroundImage: "repeating-linear-gradient(135deg, var(--series-s) 0 3px, rgb(var(--paper-pure)) 3px 6px)",
              borderRight: "1px solid var(--series-s)",
            }
          : {}),
      }}
    />
  </span>
);

interface Stage {
  label: string;
  count: number;
  tone: "kept" | "dropped" | "excluded";
  note?: string;
}

const StageTable = ({ label, stages, whole }: { label: string; stages: Stage[]; whole: number }) => (
  <ScrollFrame label={label}>
    <table className="w-full border-collapse text-control">
      <tbody>
        {stages.map((s) => (
          <tr key={s.label} className="border-b border-hair">
            <th scope="row" className={`py-1.5 pl-0 pr-2 text-left align-top font-medium ${s.tone === "kept" ? "text-ink" : "text-muted-text"}`}>
              {s.label}
              {s.note !== undefined ? <span className="block text-tiny font-normal text-muted-text">{s.note}</span> : null}
            </th>
            <td className="w-24 px-2 py-1.5 text-right align-top tabular-nums">{int(s.count)}</td>
            <td className="w-16 py-1.5 pl-2 pr-0 text-right align-top tabular-nums text-muted-text sm:pr-2">{pct(whole === 0 ? 0 : s.count / whole)}</td>
            <td className="hidden w-[34%] py-1.5 pl-2 pr-0 align-top sm:table-cell">
              <span className="flex h-[1.2rem] items-center">
                <Bar value={s.count} max={whole} tone={s.tone} />
              </span>
            </td>
          </tr>
        ))}
      </tbody>
    </table>
  </ScrollFrame>
);

const BarList = ({ rows, label }: { rows: [string, number][]; label: string }) => {
  const max = Math.max(0, ...rows.map((r) => r[1]));
  return (
    <ScrollFrame label={label}>
      <table className="w-full border-collapse text-control">
        <tbody>
          {rows.map(([k, v]) => (
            <tr key={k} className="border-b border-hair">
              <th scope="row" className="max-w-[12rem] py-1 pl-0 pr-2 text-left font-medium text-ink">
                {k}
              </th>
              <td className="w-20 py-1 pl-2 pr-0 text-right tabular-nums sm:pr-2">{int(v)}</td>
              <td className="hidden w-[45%] py-1 pl-2 pr-0 sm:table-cell">
                <span className="flex h-[1.2rem] items-center">
                  <Bar value={v} max={max} />
                </span>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </ScrollFrame>
  );
};

// ---- 1 curation funnel -----------------------------------------------------------------

export interface PathStep {
  key:
    | "records"
    | "genomes_with_a_record"
    | "genomes_passing"
    | "pairs"
    | "pairs_reported"
    | "datasets_curated"
    | "genomes_analyzed"
    | "datasets_reported";
  step: string;
  unit: "Lab records" | "Genomes" | "Pairs" | "Datasets";
  count: number;
  /** Starting count of this step's unit when the step nests inside it. */
  unitMax: number | null;
  leftOut: string;
}

/** The eight curation steps from lab records to datasets with a report, in the
 *  order curation runs, derived solely from the release manifest (the same
 *  fields `scripts/release-summary/path.py` writes to `path.json`). */
export const curationPath = (manifest: Manifest): readonly PathStep[] => {
  const f = manifest.funnel;
  const failedTotal = Object.values(f.genomes_failed_qc).reduce((a, b) => a + b, 0);
  const notPassed = failedTotal + f.genomes_without_metadata;
  const nameDropped = f.records_class_label + f.records_not_antibiotic + f.records_unresolved_name;
  const datasetsCurated = manifest.species.reduce((a, s) => a + s.datasets_curated, 0);
  const datasetsReported = manifest.species.reduce((a, s) => a + s.datasets_reported, 0);
  const genomesAnalyzed = manifest.species.reduce((a, s) => a + s.genomes_analyzed, 0);
  const { cap, capped } = selectionCap(manifest);
  const speciesAtCap = cap > 0 ? manifest.species.filter((s) => s.genomes_analyzed === cap).length : 0;
  const minClass = classMinimum(manifest.datasets);
  const minWord = minClass.kind === "value" ? int(minClass.value) : null;
  const unreported = Math.max(0, datasetsCurated - datasetsReported);
  const waitingDatasets = manifest.species
    .filter((s) => s.genomes_analyzed === 0 && s.datasets_curated > 0)
    .reduce((a, s) => a + s.datasets_curated, 0);
  const nests = (value: number, max: number): number | null => (max > 0 && value <= max ? max : null);

  return [
    {
      key: "records",
      step: "Lab records",
      unit: "Lab records",
      count: f.records,
      unitMax: null,
      leftOut: "",
    },
    {
      key: "genomes_with_a_record",
      step: "Genomes those records cover",
      unit: "Genomes",
      count: f.genomes_with_records,
      unitMax: null,
      leftOut: "",
    },
    {
      key: "genomes_passing",
      step: "Genomes that pass quality checks",
      unit: "Genomes",
      count: f.genomes_passing_qc,
      unitMax: nests(f.genomes_passing_qc, f.genomes_with_records),
      leftOut: `${int(failedTotal)} ${failedTotal === 1 ? "genome fails" : "genomes fail"} a check and ${int(f.genomes_without_metadata)} ${f.genomes_without_metadata === 1 ? "has" : "have"} no metadata`,
    },
    {
      key: "pairs",
      step: "Genome-and-drug pairs on those genomes",
      unit: "Pairs",
      count: f.pairs,
      unitMax: null,
      leftOut: `${int(nameDropped)} records on any genome name a drug class, no antibiotic, or an unresolved name; the other records on the ${int(notPassed)} ${notPassed === 1 ? "genome" : "genomes"} that did not pass are skipped with no count; ${int(f.records_no_measure)} records on passing genomes have no call and no MIC`,
    },
    {
      key: "pairs_reported",
      step: "Pairs with a reported resistant or susceptible call",
      unit: "Pairs",
      count: f.pairs_reported,
      unitMax: nests(f.pairs_reported, f.pairs),
      leftOut: `${int(f.pairs_labeled_from_mic)} pairs have a label worked out from an MIC, ${int(f.pairs_unlabeled)} have no label, ${int(f.pairs_intermediate_only)} are intermediate only, and ${int(f.pairs_conflict)} conflict`,
    },
    {
      key: "datasets_curated",
      step: "Datasets curated",
      unit: "Datasets",
      count: datasetsCurated,
      unitMax: null,
      leftOut:
        minWord !== null
          ? `species and drug combinations with under ${minWord} resistant or under ${minWord} susceptible genomes; the index does not count them`
          : "species and drug combinations under the class minimum of resistant or susceptible genomes; the index does not count them",
    },
    {
      key: "genomes_analyzed",
      step: "Genomes analyzed",
      unit: "Genomes",
      count: genomesAnalyzed,
      unitMax: nests(genomesAnalyzed, f.genomes_with_records),
      leftOut: capped
        ? `passing genomes with no unconflicted reported resistant or susceptible call on a drug of one of their species' datasets, then genomes past the cap of ${int(cap)} a species (${int(speciesAtCap)} of ${int(manifest.species.length)} species ${speciesAtCap === 1 ? "sits" : "sit"} at ${int(cap)}); the index counts neither`
        : `passing genomes with no unconflicted reported resistant or susceptible call on a drug of one of their species' datasets${waitingDatasets > 0 ? ", plus species not analyzed in this run" : ""}; the index does not count them`,
    },
    {
      key: "datasets_reported",
      step: "Datasets with a report",
      unit: "Datasets",
      count: datasetsReported,
      unitMax: nests(datasetsReported, datasetsCurated),
      leftOut:
        unreported === 0
          ? "none: every curated dataset has a report"
          : waitingDatasets === 0
            ? `${int(unreported)} of ${int(datasetsCurated)} datasets ${unreported === 1 ? "has" : "have"} one class under ${minWord ?? "the class minimum"} among the genomes kept`
            : `${int(unreported)} of ${int(datasetsCurated)} datasets have no report (${int(waitingDatasets)} in species not analyzed${
                unreported > waitingDatasets
                  ? `, ${int(unreported - waitingDatasets)} with one class under ${minWord ?? "the class minimum"} among the genomes kept`
                  : ""
              })`,
    },
  ];
};

const PathTable = ({ steps }: { steps: readonly PathStep[] }) => (
  <div id="m-path" data-curation-path="">
    <h3 tabIndex={-1} className="mb-2 text-control font-semibold text-ink">
      Curation path
    </h3>
    <ScrollFrame label="Curation path">
      <table className="w-full border-collapse text-control">
        <thead>
          <tr className="border-b border-hair">
            <th scope="col" className="py-2 pl-0 pr-2 text-left text-tiny font-medium text-blue">
              Step
            </th>
            <th scope="col" className="hidden px-2 py-2 text-left text-tiny font-medium text-blue sm:table-cell">
              Unit
            </th>
            <th scope="col" className="px-2 py-2 text-right text-tiny font-medium text-blue">
              Count
            </th>
            <th scope="col" className="hidden w-32 px-2 py-2 text-left text-tiny font-medium text-blue md:table-cell">
              Within unit
            </th>
            <th scope="col" className="py-2 pl-2 pr-0 text-left text-tiny font-medium text-blue">
              Left out
            </th>
          </tr>
        </thead>
        <tbody>
          {steps.map((s, i) => (
            <tr key={s.key} data-path-step={s.key} data-path-count={s.count} data-path-unit={s.unit} className="border-b border-hair">
              <th scope="row" className="py-2 pl-0 pr-2 text-left align-top font-medium text-ink">
                <span className="flex items-baseline gap-2">
                  <span aria-hidden="true" className="w-4 shrink-0 text-tiny tabular-nums text-muted-text">
                    {i + 1}
                  </span>
                  <span>
                    {s.step}
                    <span className="block text-tiny font-normal text-muted-text sm:hidden">{s.unit}</span>
                  </span>
                </span>
              </th>
              <td className="hidden whitespace-nowrap px-2 py-2 align-top text-tiny text-muted-text sm:table-cell">{s.unit}</td>
              <td data-path-value="" className="whitespace-nowrap px-2 py-2 text-right align-top font-semibold tabular-nums text-ink">
                {int(s.count)}
              </td>
              <td className="hidden w-32 px-2 py-2 align-top md:table-cell">
                {s.unitMax !== null ? (
                  <span className="flex h-[1.2rem] items-center">
                    <Bar value={s.count} max={s.unitMax} tone="kept" />
                  </span>
                ) : (
                  <span className="text-tiny text-muted-text" aria-hidden="true">
                    —
                  </span>
                )}
              </td>
              <td data-path-left-out="" className="py-2 pl-2 pr-0 align-top text-note text-muted-text">
                {s.leftOut === "" ? <span aria-label="None">—</span> : s.leftOut}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </ScrollFrame>
    <p className="type-note mt-2">
      Counts fall in four units (lab records, genomes, pairs, datasets); within one unit each later count nests inside the earlier one, and each bar compares a later step to
      the first count in its unit. Each count of genomes describes what laboratories deposited in BV-BRC; none says how common resistance is. Each species&apos; count of
      genomes with any label (in Species coverage below, intermediate calls and MIC-inferred labels included) is not a step of this path.
    </p>
  </div>
);

const FunnelSection = ({ manifest }: { manifest: Manifest }) => {
  const f = manifest.funnel;
  const failed = Object.values(f.genomes_failed_qc).reduce((a, b) => a + b, 0);
  const notPassed = failed + f.genomes_without_metadata;
  const nameDropped = f.records_class_label + f.records_not_antibiotic + f.records_unresolved_name;
  const genomesAdd = f.genomes_without_metadata + failed + f.genomes_passing_qc === f.genomes_with_records;
  const steps = curationPath(manifest);
  return (
    <DataSection id={DATA_ANCHORS.funnel} n={1} title="Curation funnel">
      <P>From BV-BRC lab records to datasets with a report in eight steps, in the order curation runs.</P>
      <PathTable steps={steps} />
      <Detail>
        Breakdowns of the left-out counts appear below in the same order: blue bars mark what goes on to the models, hatched blue marks MIC-inferred labels kept in records but
        excluded from training, and grey marks what is set aside.
      </Detail>
      <div className="grid grid-cols-1 gap-x-12 gap-y-6 xl:grid-cols-2">
        <div id="m-genomes">
          <h3 className="mb-2 text-control font-semibold text-ink">Genomes and QC outcomes</h3>
          <StageTable
            label="Genomes"
            whole={f.genomes_with_records}
            stages={[
              { label: "Passed QC", count: f.genomes_passing_qc, tone: "kept" },
              { label: "No genome metadata", count: f.genomes_without_metadata, tone: "dropped" },
              ...Object.entries(f.genomes_failed_qc)
                .sort((a, b) => b[1] - a[1])
                .map(([reason, count]): Stage => ({ label: `Failed QC: ${reason}`, count, tone: "dropped" })),
            ]}
          />
          <p className="type-note mt-2">
            {genomesAdd
              ? `${int(f.genomes_with_records)} genomes with records = ${int(f.genomes_passing_qc)} passing + ${int(f.genomes_without_metadata)} without metadata + ${int(failed)} failing QC.`
              : `${int(f.genomes_with_records)} genomes had records.`}
          </p>
        </div>
        <div id="m-records">
          <h3 className="mb-2 text-control font-semibold text-ink">Records set aside before pairing</h3>
          <StageTable
            label="Records set aside"
            whole={f.records}
            stages={[
              { label: "Name is a drug class, not a drug", count: f.records_class_label, tone: "dropped", note: "For example “beta-lactam”." },
              { label: "Name is not an antibiotic", count: f.records_not_antibiotic, tone: "dropped" },
              { label: "Name could not be resolved", count: f.records_unresolved_name, tone: "dropped" },
              { label: "No call and no MIC", count: f.records_no_measure, tone: "dropped", note: "Counted on passing genomes only." },
            ]}
          />
          <p className="type-note mt-2">
            Shares are of all {int(f.records)} BV-BRC records. The three name checks ({int(nameDropped)} records) read the records of every genome, passing or not; other
            records on the {int(notPassed)} {notPassed === 1 ? "genome" : "genomes"} that did not pass are skipped with no count; no-call/no-MIC records are counted on
            passing genomes only.
          </p>
        </div>
        <div id="m-pairs" className="xl:col-span-2">
          <h3 tabIndex={-1} className="mb-2 text-control font-semibold text-ink">Genome–antibiotic pairs</h3>
          <StageTable
            label="Pairs"
            whole={f.pairs}
            stages={[
              { label: "Resistant or susceptible as reported", count: f.pairs_reported, tone: "kept" },
              {
                label: "Labeled from MIC at the inferred breakpoint",
                count: f.pairs_labeled_from_mic,
                tone: "excluded",
                note: "An MIC but no reported call. Kept in the curated records; no model trains or is scored on these labels.",
              },
              { label: "Conflicting calls (both resistant and susceptible)", count: f.pairs_conflict, tone: "dropped" },
              { label: "Intermediate only", count: f.pairs_intermediate_only, tone: "dropped" },
              { label: "No label possible", count: f.pairs_unlabeled, tone: "dropped", note: "No call, and no MIC that falls clearly on one side of a breakpoint." },
            ]}
          />
          <p className="type-note mt-2">
            {int(f.pairs)} pairs = {int(f.pairs_reported)} reported + {int(f.pairs_labeled_from_mic)} from MIC + {int(f.pairs_conflict)} conflicting + {int(f.pairs_intermediate_only)} intermediate + {int(f.pairs_unlabeled)} unlabeled.
          </p>
        </div>
      </div>
    </DataSection>
  );
};

// ---- 2 isolate sources -----------------------------------------------------------------

const Sources = ({ corpus }: { corpus: Corpus }) => {
  const countries = corpus.countries.slice(0, 12);
  const totalCountry = corpus.countries.reduce((a, c) => a + c[1], 0);
  const years = [...corpus.years].sort((a, b) => a[0].localeCompare(b[0]));
  const top = corpus.countries[0];
  return (
    <>
      <p className="prose-measure text-note text-ink">
        BV-BRC isolates over-represent outbreaks, some countries, and some years.{" "}
        {top === undefined
          ? "No isolation country is recorded for the labeled genomes."
          : `${top[0]} supplies ${pct(top[1] / Math.max(totalCountry, 1))} of labeled genomes with a recorded country${years.length > 0 ? `; collection years run from ${years[0]![0]} to ${years[years.length - 1]![0]}` : ""}.`}
      </p>
      <div className="grid grid-cols-1 gap-x-12 gap-y-6 lg:grid-cols-2">
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Labeled genomes by isolation country</h3>
          <BarList rows={countries} label="Countries" />
          {corpus.countries.length > countries.length ? <p className="type-note mt-2">Top {countries.length} of {int(corpus.countries.length)} countries.</p> : null}
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Labeled genomes by collection year</h3>
          <YearColumns years={years} />
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Lab records by testing method</h3>
          <BarList rows={corpus.methods} label="Methods" />
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Lab records by testing standard</h3>
          <BarList rows={corpus.standards} label="Standards" />
        </div>
      </div>
    </>
  );
};

const SourcesSection = () => {
  const corpus = useCorpus();
  return (
    <DataSection id={DATA_ANCHORS.sources} n={2} title="Isolate sources">
      <P>
        An isolate is one sample of bacteria grown in a lab from one source, and each isolate here has one genome. Where and when the isolates were collected shapes
        every score.
      </P>
      {corpus.status === "loading" ? <PanelState kind="loading" title="Loading the corpus file" /> : null}
      {corpus.status === "error" ? <PanelState kind="error" title="The corpus file could not be read" message={loadErrorText(corpus.error)} /> : null}
      {corpus.status === "ready" ? <Sources corpus={corpus.data} /> : null}
    </DataSection>
  );
};

// ---- 3 species coverage ----------------------------------------------------------------

const SpeciesSection = ({ manifest, nav }: { manifest: Manifest; nav?: Navigator }) => (
  <DataSection id={DATA_ANCHORS.species} n={3} title="Species coverage">
    <P>
      Each row is one bacterial species in the curated corpus, with its reported datasets and selected genomes. Its count of genomes with any label (intermediate calls and
      MIC-inferred labels included) is shown for context and is not a step of the curation path.
    </P>
    <SpeciesCounts manifest={manifest} nav={nav} id="data-species-table" />
  </DataSection>
);

// ---- 4 labels and inclusion ------------------------------------------------------------

const LabelsSection = ({ manifest }: { manifest: Manifest }) => {
  const { cap, capped } = selectionCap(manifest);
  return (
    <DataSection id={DATA_ANCHORS.labels} n={4} title="Labels and inclusion">
      <P>
        The genomes and their lab results come from BV-BRC, a public database of bacterial genomes. A lab grows the bacterium next to the antibiotic. If the drug stops it
        growing, the lab reports it susceptible (S); if it keeps growing, resistant (R). That reported call is the label: what every model tries to predict and what every
        score is checked against.
      </P>
      <P>
        Some labs also report an MIC (minimum inhibitory concentration): the lowest amount of the drug, in mg/L, that stopped the bacterium growing. Labs test the drug at
        levels that double each time (0.5, 1, 2, 4 mg/L), and each level is one dilution step.
      </P>
      <Detail>
        Every learner trains and is scored on genomes with a reported resistant (R) or susceptible (S) lab call for the drug and no conflicting call, and the known-gene rules, panels, and MIC
        regression are scored on the same genomes. A genome and antibiotic pair with both an R and an S call is a conflict and gets no label; a pair with only intermediate calls is never
        an R or S label. A pair called susceptible and intermediate is labeled susceptible, and a pair called intermediate and resistant is labeled resistant.
      </Detail>
      <LabelRuleTable />
      <ClindamycinMethods manifest={manifest} />
      <Detail>
        A pair with an MIC but no reported call gets a label from an inferred breakpoint during curation, and no model trains or is scored on that label. A dataset report counts the labels
        it left out in its Method and coverage fold; a report written before that count was added shows it as not recorded. The species view&apos;s resistance profiles and cross-resistance
        use the consolidated labels, including those inferred from an MIC.
      </Detail>
      <Detail>
        Each inferred breakpoint is the dilution step whose calls agree best with the reported calls among isolates with both, fitted with at least 50 such isolates (10 of each call) and
        kept at 95% agreement or more; censored MICs that fall on neither side get no label. These are empirical thresholds from the data, not CLSI or EUCAST breakpoints.
      </Detail>
      <Detail>
        Antibiotic names are resolved to one drug vocabulary. Records naming a drug class such as “beta-lactam”, a substance that is not an antibiotic, or a name that cannot be resolved
        are set aside, as are records with neither a call nor an MIC.
      </Detail>
      <Detail>
        A dataset is one species and one antibiotic. {minimumSentence(classMinimum(manifest.datasets))}{" "}
        {capped
          ? `That minimum runs before and after the cap of ${int(cap)} genomes per species, and QC and FASTA availability reduce some cohorts further.`
          : "QC and FASTA availability decide which labeled genomes are selected."}
      </Detail>
    </DataSection>
  );
};

interface LabelRuleRow {
  key: string;
  calls: string;
  label: string;
  used: string;
  kept: boolean;
}

/** Consolidation of reported lab calls on one genome-and-drug pair under the
 *  release rule (`r-or-s-wins`, `src/curate.rs:77-89`, `src/curate.rs:753-793`). */
const LABEL_RULES: readonly LabelRuleRow[] = [
  { key: "s", calls: "S only", label: "Susceptible (S)", used: "Yes", kept: true },
  { key: "r", calls: "R only", label: "Resistant (R)", used: "Yes", kept: true },
  { key: "s-i", calls: "S and I", label: "Susceptible (S)", used: "Yes", kept: true },
  { key: "i-r", calls: "I and R", label: "Resistant (R)", used: "Yes", kept: true },
  { key: "i", calls: "I only", label: "Intermediate (I); no R or S label", used: "No", kept: false },
  { key: "conflict", calls: "R and S (with or without I)", label: "Conflict; no label", used: "No", kept: false },
  { key: "mic-only", calls: "No reported call, MIC only", label: "Inferred R or S at empirical breakpoint, else none", used: "No (curated records only)", kept: false },
];

const LabelRuleTable = () => (
  <div id="m-label-rule" className="max-w-[44rem]">
    <ScrollFrame label="Label rule">
      <table className="w-full border-collapse text-control">
        <thead>
          <tr className="border-b border-hair">
            <th scope="col" className="py-2 pl-0 pr-3 text-left text-tiny font-medium text-blue">
              Reported calls on pair
            </th>
            <th scope="col" className="px-3 py-2 text-left text-tiny font-medium text-blue">
              Label
            </th>
            <th scope="col" className="py-2 pl-3 pr-0 text-left text-tiny font-medium text-blue">
              Used by models
            </th>
          </tr>
        </thead>
        <tbody>
          {LABEL_RULES.map((row) => (
            <tr key={row.key} data-label-rule-row={row.key} className="border-b border-hair">
              <th scope="row" className={`py-1.5 pl-0 pr-3 text-left align-top font-medium ${row.kept ? "text-ink" : "text-muted-text"}`}>
                {row.calls}
              </th>
              <td className={`px-3 py-1.5 align-top ${row.kept ? "text-ink" : "text-muted-text"}`}>{row.label}</td>
              <td className={`py-1.5 pl-3 pr-0 align-top ${row.kept ? "font-medium text-ink" : "text-muted-text"}`}>{row.used}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </ScrollFrame>
  </div>
);

// ---- contents strip --------------------------------------------------------------------

const DATA_PARTS: readonly { title: string; text: string; to: string }[] = [
  { title: "Curation funnel", text: "The eight-step path from BV-BRC lab records to datasets with a report, plus QC, record, and pair breakdowns.", to: DATA_ANCHORS.funnel },
  { title: "Isolate sources", text: "Where and when isolates were collected, and which lab methods and standards tested them.", to: DATA_ANCHORS.sources },
  { title: "Species coverage", text: "Genomes analyzed, any-label pools, and reported datasets for each species.", to: DATA_ANCHORS.species },
  { title: "Labels and inclusion", text: "Reported R or S calls, MICs, conflicts, inferred breakpoints, and cohort rules.", to: DATA_ANCHORS.labels },
];

const DataMap = () => (
  <ol aria-label="Corpus sections" className="mt-8 grid grid-cols-1 gap-x-12 gap-y-4 lg:grid-cols-2">
    {DATA_PARTS.map((part, i) => (
      <li key={part.to} className="grid grid-cols-[1.5rem_minmax(0,1fr)] gap-x-2">
        <span aria-hidden="true" className="text-body font-semibold tabular-nums text-muted-text">
          {i + 1}
        </span>
        <div className="prose-measure min-w-0">
          <a
            href={`#${part.to}`}
            onClick={(e) => {
              e.preventDefault();
              showSection(part.to, true);
            }}
            className={`text-body font-semibold text-ink ${STATUS_LINK}`}
          >
            {part.title}
          </a>
          <p className="mt-0.5 text-note text-muted-text">{part.text}</p>
        </div>
      </li>
    ))}
  </ol>
);

export const DataView = ({ manifest, nav }: { manifest: Manifest; nav?: Navigator }) => {
  const [hash] = useState(() => (typeof window === "undefined" ? "" : decodeURIComponent(window.location.hash.slice(1))));
  useEffect(() => {
    if (hash !== "") showSection(hash);
  }, [hash]);
  const analyzed = manifest.species.reduce((a, s) => a + s.genomes_analyzed, 0);
  const datasetsCurated = manifest.species.reduce((a, s) => a + s.datasets_curated, 0);
  const activeSpecies = manifest.species.filter((s) => s.genomes_analyzed > 0).length;
  const f = manifest.funnel;
  const jump = (id: string) => (e: MouseEvent<HTMLAnchorElement>) => {
    e.preventDefault();
    showSection(id, true);
  };
  const summaryStats: BolsterStat[] = [
    { key: "records", label: "Public records", value: int(f.records), sub: "BV-BRC lab records", href: `#${DATA_ANCHORS.funnel}`, onJump: jump(DATA_ANCHORS.funnel) },
    { key: "pairs", label: "Reported pairs", value: int(f.pairs_reported), sub: `of ${int(f.pairs)} genome-drug pairs`, href: "#m-pairs", onJump: jump("m-pairs") },
    {
      key: "genomes",
      label: "Genomes analyzed",
      value: int(analyzed),
      sub: analyzed <= f.genomes_passing_qc ? `of ${int(f.genomes_passing_qc)} passing QC` : `across ${int(activeSpecies)} species`,
      href: `#${DATA_ANCHORS.species}`,
      onJump: jump(DATA_ANCHORS.species),
    },
    {
      key: "datasets",
      label: "Datasets",
      value: int(manifest.datasets.length),
      sub: datasetsCurated > manifest.datasets.length ? `of ${int(datasetsCurated)} curated, ${int(activeSpecies)} species` : `across ${int(activeSpecies)} species`,
      href: `#${DATA_ANCHORS.species}`,
      onJump: jump(DATA_ANCHORS.species),
    },
  ];
  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader
        title="Data"
        purpose="Where GRM's genomes and lab calls come from, what passes curation, and which records are set aside."
        status={<ReportStatus manifest={manifest} />}
      />
      <div className="mt-6">
        <BolsterStatStrip stats={summaryStats} ariaLabel="Corpus totals" />
      </div>
      <DataMap />
      <div className="mt-12 space-y-12">
        <FunnelSection manifest={manifest} />
        <SourcesSection />
        <SpeciesSection manifest={manifest} nav={nav} />
        <LabelsSection manifest={manifest} />
      </div>
    </div>
  );
};
