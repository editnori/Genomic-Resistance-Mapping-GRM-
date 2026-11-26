// Data: where GRM's genomes and lab calls come from, what passes curation,
// what is set aside, and which rules decide inclusion. Four linear sections
// put the numbers first (the curation funnel and isolate sources), then
// species coverage, then the label and inclusion rules. Every count comes
// from the manifest or the corpus file at render time.

import { useEffect, useState, type ReactNode } from "react";
import type { Corpus, Manifest } from "../data/schema";
import { int, pct } from "../data/format";
import { loadErrorText, useCorpus } from "../data/load";
import { classMinimum } from "../data/derive/overview";
import type { Navigator } from "../app/router";
import { PanelState, ScrollFrame, VIEW_SHELL } from "../ui/primitives";
import { ViewHeader } from "../ui/ViewHeader";
import { ReportStatus } from "../app/ReportStatus";
import { YearColumns } from "../charts/YearColumns";
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

/** A count's length on a track of the whole: blue for what was kept, grey
 *  for what was set aside. */
const Bar = ({ value, max, tone = "kept" }: { value: number; max: number; tone?: "kept" | "dropped" }) => (
  <span className="relative block h-2.5 w-full min-w-[6rem] overflow-hidden rounded-sm bg-surface-subtle" aria-hidden="true">
    <span
      className={`absolute inset-y-0 left-0 rounded-r-sm ${tone === "kept" ? "bg-[color:var(--series-s)]" : "bg-[color:var(--series-neutral)]"}`}
      style={{ width: `${max === 0 ? 0 : Math.max(0.5, (value / max) * 100)}%` }}
    />
  </span>
);

interface Stage {
  label: string;
  count: number;
  tone: "kept" | "dropped";
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

const FunnelSection = ({ manifest }: { manifest: Manifest }) => {
  const f = manifest.funnel;
  const failed = Object.values(f.genomes_failed_qc).reduce((a, b) => a + b, 0);
  const genomesAdd = f.genomes_without_metadata + failed + f.genomes_passing_qc === f.genomes_with_records;
  return (
    <DataSection id={DATA_ANCHORS.funnel} n={1} title="Curation funnel">
      <P>
        These tables follow the public records down to the genomes and antibiotics the models use. In the first two tables, the top row goes on to the models, with a blue
        bar on wider screens. Every other row is set aside, and its label says why. QC (quality control) checks each genome&apos;s assembled DNA, and each reason a genome failed has its own row.
      </P>
      <div className="grid grid-cols-1 gap-x-12 gap-y-6 xl:grid-cols-2">
        <div id="m-pairs" className="xl:col-span-2">
          <h3 className="mb-2 text-control font-semibold text-ink">Genome–antibiotic pairs</h3>
          <StageTable
            label="Pairs"
            whole={f.pairs}
            stages={[
              { label: "Resistant or susceptible as reported", count: f.pairs_reported, tone: "kept" },
              {
                label: "Labeled from MIC at the inferred breakpoint",
                count: f.pairs_labeled_from_mic,
                tone: "dropped",
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
              { label: "No call and no MIC", count: f.records_no_measure, tone: "dropped" },
            ]}
          />
          <p className="type-note mt-2">Shares are of all {int(f.records)} BV-BRC records. Records on genomes that failed QC are also set aside; they are counted by genome.</p>
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
    <P>Each row is one bacterial species in the curated corpus, with its reported datasets and selected genomes.</P>
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
          ? `Selection is capped at ${int(cap)} genomes per species when more are eligible, and QC and FASTA availability reduce some cohorts further.`
          : "QC and FASTA availability decide which labeled genomes are selected."}
      </Detail>
    </DataSection>
  );
};

// ---- contents strip --------------------------------------------------------------------

const DATA_PARTS: readonly { title: string; text: string; to: string }[] = [
  { title: "Curation funnel", text: "Pairs kept vs set aside, genomes passing QC, and records dropped before pairing.", to: DATA_ANCHORS.funnel },
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
  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader
        title="Data"
        purpose="Where GRM's genomes and lab calls come from, what passes curation, and which records are set aside."
        status={<ReportStatus manifest={manifest} />}
      />
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
