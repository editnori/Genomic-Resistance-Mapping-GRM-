// Methods: the evaluation audit, where every number comes from, how the
// learners are judged, and each limit beside the definition it qualifies.

import { useState, type ReactNode } from "react";
import type { Corpus, DrugNameRow, Manifest } from "../data/schema";
import { capitalize, int, median, mgL, pct, plural, speciesShort } from "../data/format";
import { loadErrorText, useCorpus } from "../data/load";
import { manifestAudit, nullFitWords } from "../data/derive/evaluation";
import { classMinimum, datasetRule } from "../data/derive/overview";
import { notScored } from "../data/derive/tests";
import { DataTable, type Column } from "../ui/DataTable";
import { FOCUS_RING, PanelState, ScrollFrame, Section, Segmented, WIDE_SHELL } from "../ui/primitives";
import { HoverTip } from "../ui/HoverTip";
import { YearColumns } from "../charts/YearColumns";
import { CorpusFigure } from "./findings";

/** The step that resolved one name component, in words. */
export const sourceWords = (source: string): string => {
  if (source === "sigdrug:vocabulary") return "SigDrug vocabulary";
  if (source.startsWith("sigdrug:")) return "SigDrug folds";
  if (source === "grm:vocabulary") return "GRM vocabulary";
  if (source === "grm:override") return "GRM override";
  const edit = /^grm:edit-distance-(\d+)$/.exec(source);
  if (edit !== null) return `Edit distance ${edit[1]}`;
  return source === "unresolved" ? "Unresolved" : source;
};

const KIND_WORD: Record<DrugNameRow["kind"], string> = { drug: "Drug", class: "Class label", excluded: "Not an antibiotic", unresolved: "Unresolved" };

const Bar = ({ value, max, tone = "kept" }: { value: number; max: number; tone?: "kept" | "dropped" }) => (
  <span className="relative block h-2.5 w-full min-w-[6rem] overflow-hidden rounded-sm bg-surface-subtle" aria-hidden="true">
    <span
      className={`absolute inset-y-0 left-0 rounded-r-sm ${tone === "kept" ? "bg-[color:var(--chart-mark)]" : "bg-[color:var(--series-neutral)]"}`}
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
            <th scope="row" className={`px-2 py-1.5 text-left font-medium ${s.tone === "kept" ? "text-ink" : "text-muted-text"}`}>
              {s.label}
              {s.note !== undefined ? <span className="block text-tiny font-normal text-muted-text">{s.note}</span> : null}
            </th>
            <td className="w-24 px-2 py-1.5 text-right tabular-nums">{int(s.count)}</td>
            <td className="w-16 px-2 py-1.5 text-right tabular-nums text-muted-text">{pct(whole === 0 ? 0 : s.count / whole)}</td>
            <td className="w-[40%] px-2 py-1.5">
              <Bar value={s.count} max={whole} tone={s.tone} />
            </td>
          </tr>
        ))}
      </tbody>
    </table>
  </ScrollFrame>
);

const Funnel = ({ manifest }: { manifest: Manifest }) => {
  const f = manifest.funnel;
  const failed = Object.values(f.genomes_failed_qc).reduce((a, b) => a + b, 0);
  const curated = manifest.species.reduce((a, s) => a + s.datasets_curated, 0);
  const genomesAdd = f.genomes_without_metadata + failed + f.genomes_passing_qc === f.genomes_with_records;
  return (
    <Section
      id="m-funnel"
      title="Curation funnel"
      answer={`${int(f.records)} BV-BRC lab records became ${int(f.pairs_reported)} genome and antibiotic pairs with a reported resistant or susceptible call on ${int(f.genomes_passing_qc)} genomes that passed QC, and ${plural(curated, "curated dataset")}, ${int(manifest.datasets.length)} analyzed so far.`}
    >
      <CorpusFigure manifest={manifest} />
      <p className="type-note mt-3">Each count is its own unit: records, genome and antibiotic pairs, genomes, and datasets. The tables below count each step in the unit it removes.</p>
      <div className="mt-8 grid gap-x-12 gap-y-8 xl:grid-cols-2">
        <div>
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
          <p className="type-note mt-2">Shares are of all {int(f.records)} records. Records on genomes that failed QC are also set aside; they are counted by genome below.</p>
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Genomes</h3>
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
        <div className="xl:col-span-2">
          <h3 className="mb-2 text-control font-semibold text-ink">Genome and antibiotic pairs</h3>
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
            {int(f.pairs)} pairs = {int(f.pairs_reported)} reported + {int(f.pairs_labeled_from_mic)} from MIC + {int(f.pairs_conflict)} conflicting + {int(f.pairs_intermediate_only)} intermediate + {int(f.pairs_unlabeled)} unlabeled. A dataset needs {datasetRule(classMinimum(manifest.datasets))}.
          </p>
        </div>
      </div>
    </Section>
  );
};

const nameColumns: Column<DrugNameRow>[] = [
  { id: "raw", header: "Name as submitted", sortValue: (r) => r.raw.toLowerCase(), firstDir: "asc", cell: (r) => <span className="whitespace-normal break-words">{r.raw}</span> },
  { id: "canonical", header: "Canonical name", sortValue: (r) => r.canonical, firstDir: "asc", cell: (r) => r.canonical },
  { id: "kind", header: "Kind", sortValue: (r) => r.kind, firstDir: "asc", cell: (r) => <span className={r.kind === "drug" ? "" : "text-muted-text"}>{KIND_WORD[r.kind]}</span> },
  { id: "records", header: "Records", align: "right", sortValue: (r) => r.records, cell: (r) => int(r.records) },
  {
    id: "steps",
    header: "Resolved by",
    tip: [{ term: "Resolved by", text: "The step that resolved each component: SigDrug vocabulary, SigDrug folds, GRM vocabulary, GRM override, edit distance, or unresolved." }],
    cell: (r) =>
      r.components.length === 0 ? (
        <span className="text-muted-text">Whole-name match</span>
      ) : (
        <span className="flex flex-col gap-0.5">
          {r.components.map((c, i) => (
            <span key={i} className="whitespace-normal">
              <span className="text-muted-text">{c.text}</span> → {c.canonical}{" "}
              <HoverTip focusable={false} lines={[{ term: sourceWords(c.source), source: c.source }]}>
                <span className={`text-tiny font-medium ${c.source === "unresolved" ? "text-coral-text" : "text-blue"}`}>{sourceWords(c.source)}</span>
              </HoverTip>
            </span>
          ))}
        </span>
      ),
  },
];

type KindFilter = "all" | DrugNameRow["kind"];

const Names = ({ corpus, manifest }: { corpus: Corpus; manifest: Manifest }) => {
  const [kind, setKind] = useState<KindFilter>("all");
  const [query, setQuery] = useState("");
  const names = corpus.drug_names;
  const q = query.trim().toLowerCase();
  const rows = names.filter((r) => (kind === "all" || r.kind === kind) && (q === "" || `${r.raw} ${r.canonical}`.toLowerCase().includes(q)));
  const total = names.reduce((a, r) => a + r.records, 0);
  const drugRecords = names.filter((r) => r.kind === "drug").reduce((a, r) => a + r.records, 0);
  const bySource = new Map<string, number>();
  names.forEach((r) => r.components.forEach((c) => bySource.set(sourceWords(c.source), (bySource.get(sourceWords(c.source)) ?? 0) + 1)));
  const count = (k: DrugNameRow["kind"]) => names.filter((r) => r.kind === k).length;
  return (
    <Section
      id="m-names"
      title="Antibiotic name resolution"
      answer={`${int(count("drug"))} of ${plural(names.length, "distinct name")} resolve to a drug, covering ${pct(total === 0 ? 0 : drugRecords / total)} of their records; ${int(count("unresolved"))} ${count("unresolved") === 1 ? "stays" : "stay"} unresolved.`}
    >
      <p className="type-note -mt-1 mb-3">
        Components resolved per step: {[...bySource.entries()].sort((a, b) => b[1] - a[1]).map(([s, n]) => `${s} ${int(n)}`).join(", ")}.
        {manifest.name_vocabulary !== null ? ` Reference vocabulary: ${manifest.name_vocabulary}.` : ""} Combination names split at “/”, and inhibitor partners sort last.
      </p>
      <div className="mb-3 flex flex-wrap items-center gap-3">
        <Segmented<KindFilter>
          ariaLabel="Kind"
          options={[
            { id: "all", label: "All" },
            { id: "drug", label: "Drugs" },
            { id: "class", label: "Class labels" },
            { id: "excluded", label: "Not antibiotics" },
            { id: "unresolved", label: "Unresolved" },
          ]}
          value={kind}
          onSelect={setKind}
        />
        <label className="min-w-[14rem] flex-1">
          <span className="sr-only">Search names</span>
          <input
            type="search"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            placeholder="Search names"
            className={`h-7 w-full max-w-sm rounded-full border border-hair bg-paper-pure px-3 text-control text-ink placeholder:text-muted-text hover:border-blue/40 ${FOCUS_RING}`}
          />
        </label>
        {q !== "" || kind !== "all" ? <span className="text-control text-muted-text">{int(rows.length)} matching</span> : null}
      </div>
      <DataTable label="Antibiotic names" rows={rows} columns={nameColumns} rowKey={(r) => r.raw} initialSort={{ id: "records", dir: "desc" }} showAllMatches={q !== ""} empty={<p className="type-note">No name matches.</p>} />
    </Section>
  );
};

const Breakpoints = ({ corpus }: { corpus: Corpus }) => {
  const rows = corpus.breakpoints;
  type B = (typeof rows)[number];
  const columns: Column<B>[] = [
    { id: "antibiotic", header: "Antibiotic", sortValue: (r) => r.antibiotic, firstDir: "asc", cell: (r) => capitalize(r.antibiotic) },
    { id: "species", header: "Species", sortValue: (r) => r.species, firstDir: "asc", cell: (r) => <i>{speciesShort(r.species)}</i> },
    { id: "bp", header: "Resistant above", align: "right", sortValue: (r) => r.r_above_log2, cell: (r) => `${mgL(r.r_above_log2)} mg/L` },
    { id: "agreement", header: "Agreement", align: "right", sortValue: (r) => r.agreement, cell: (r) => pct(r.agreement) },
    { id: "n", header: "Isolates compared", align: "right", sortValue: (r) => r.n, cell: (r) => int(r.n) },
    { id: "labels", header: "Labels from MIC", align: "right", sortValue: (r) => r.labels_from_mic, cell: (r) => int(r.labels_from_mic) },
  ];
  return (
    <Section
      id="m-breakpoints"
      title="Inferred MIC breakpoints"
      answer={
        rows.length === 0
          ? "No breakpoint could be inferred."
          : `${plural(rows.length, "breakpoint")} were inferred from isolates with both a reported call and an MIC; the median agrees with ${pct(median(rows.map((r) => r.agreement)))} of the calls it was fitted on.`
      }
    >
      <DataTable label="Breakpoints" rows={rows} columns={columns} rowKey={(r) => `${r.species}|${r.antibiotic}`} initialSort={{ id: "labels", dir: "desc" }} empty={null} />
      <p className="type-note mt-3">
        The breakpoint is the dilution step whose calls agree best with the reported resistant and susceptible calls among isolates with both. It is fitted only with at least 50 such isolates, 10 of each call, and kept only when it agrees with at least 95% of them. Censored MICs that do not fall clearly on one side count against agreement and receive no label. These are empirical thresholds from the data, not CLSI or EUCAST breakpoints. The labels they give stay in the curated records, counted in the Labels from MIC column, and no model trains or is scored on them.
      </p>
    </Section>
  );
};

const BarList = ({ rows, label }: { rows: [string, number][]; label: string }) => {
  const max = Math.max(0, ...rows.map((r) => r[1]));
  return (
    <ScrollFrame label={label}>
      <table className="w-full border-collapse text-control">
        <tbody>
          {rows.map(([k, v]) => (
            <tr key={k} className="border-b border-hair">
              <th scope="row" className="max-w-[12rem] px-2 py-1 text-left font-medium text-ink">
                {k}
              </th>
              <td className="w-20 px-2 py-1 text-right tabular-nums">{int(v)}</td>
              <td className="w-[45%] px-2 py-1">
                <Bar value={v} max={max} />
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </ScrollFrame>
  );
};

const Sources = ({ corpus }: { corpus: Corpus }) => {
  const countries = corpus.countries.slice(0, 12);
  const totalCountry = corpus.countries.reduce((a, c) => a + c[1], 0);
  const years = [...corpus.years].sort((a, b) => a[0].localeCompare(b[0]));
  const top = corpus.countries[0];
  return (
    <Section
      id="m-sources"
      title="Where the isolates come from"
      answer={
        top === undefined
          ? "No isolation country is recorded for the labeled genomes."
          : `${top[0]} supplies ${pct(top[1] / Math.max(totalCountry, 1))} of labeled genomes with a recorded country${years.length > 0 ? `; collection years run from ${years[0]![0]} to ${years[years.length - 1]![0]}` : ""}.`
      }
    >
      <p className="type-note mb-4">BV-BRC isolates over-represent outbreaks, some countries, and some years. Every score is agreement with lab calls for isolates like these and would change with a different population structure.</p>
      <div className="grid gap-x-12 gap-y-8 lg:grid-cols-2">
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Isolation country, labeled genomes</h3>
          <BarList rows={countries} label="Countries" />
          {corpus.countries.length > countries.length ? <p className="type-note mt-2">Top {countries.length} of {int(corpus.countries.length)} countries.</p> : null}
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Collection year, labeled genomes</h3>
          <YearColumns years={years} />
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Laboratory method, records</h3>
          <BarList rows={corpus.methods} label="Methods" />
        </div>
        <div>
          <h3 className="mb-2 text-control font-semibold text-ink">Testing standard, records</h3>
          <BarList rows={corpus.standards} label="Standards" />
        </div>
      </div>
    </Section>
  );
};

const Def = ({ term, children }: { term: string; children: ReactNode }) => (
  <div className="border-b border-hair py-3">
    <dt className="text-control font-semibold text-ink">{term}</dt>
    <dd className="prose-measure mt-1 text-note text-ink">{children}</dd>
  </div>
);

const Design = ({ manifest }: { manifest: Manifest }) => (
  <Section id="m-design" title="Learners and evaluation">
    <dl className="grid gap-x-12 lg:grid-cols-2">
      <Def term="Labels">
        Every learner trains and is scored on genomes with a reported resistant or susceptible lab call, and the catalog and panel comparisons are scored on the same genomes. A genome with an MIC but no reported call gets a label from the inferred breakpoint during curation, and no model trains or is scored on that label. A dataset report counts the labels it left out in its Method and coverage fold; a report written before that count was added shows it as not recorded. Calls come from different methods and standards across years, and breakpoints changed over time.
      </Def>
      <Def term="Features">
        {manifest.k === null
          ? "Each genome is counted into k-mers with GRM's own counter. The published reports do not all record the same k, so each dataset's Method and coverage fold states its own."
          : `Each genome is counted into its ${manifest.k}-mers with GRM's own counter.`}{" "}
        K-mers with the same presence pattern across a species&apos; genomes are one feature, kept when enough isolates carry it (each dataset&apos;s minimum is in its Method and coverage fold). Whether that filtering saw only training-fold isolates is stated by the evaluation audit above.
      </Def>
      <Def term="Set covering machine (SCM)">
        A short conjunction or disjunction of presence and absence rules, chosen greedily, the model Kover uses. It gives a yes or no call and is read as written.
      </Def>
      <Def term="Decision tree (CART)">A class-balanced Gini tree on presence and absence, its depth chosen by inner cross-validation.</Def>
      <Def term="Random forest and gradient boosting">200 randomized trees, and boosted trees with early stopping on a held-out slice; both score each isolate between 0 and 1.</Def>
      <Def term="L1 logistic regression (lasso)">A sparse linear model on screened features, its penalty chosen on a held-out slice.</Def>
      <Def term="Learner cross-check">
        On seven small test fixtures and on the 868-isolate <i>E. coli</i> ceftriaxone dataset from an earlier report, run with the same learner code as this report, each of the five learners matched a reference implementation of the same model (scikit-learn, XGBoost, or a reimplementation of its documented rules) exactly or within that reference&apos;s own run-to-run noise.
      </Def>
      <Def term="Known-gene catalog">
        The catalog-containment rule finds genes of the AMRFinderPlus catalog ({manifest.amrfinder_version === "" ? "version not recorded" : manifest.amrfinder_version}) by k-mer containment and applies a fixed rule per drug: resistant when the genome carries a gene of the classes or subclasses that act on that drug. It does not train. It does not test the point mutations AMRFinderPlus also describes, so mutational resistance counts against it.
      </Def>
      <Def term="Random split">Stratified 5-fold cross-validation, the scheme most published numbers use.</Def>
      <Def term="Lineage-held-out split">
        Whole lineages (MinHash clusters at a fixed Mash distance) go to one fold, so no test isolate has a near-identical relative in training. The difference from the random split estimates how much of the score depends on near-identical relatives in training, for this sampled population. It runs when a dataset has at least 10 lineages and the lineages can be split into folds that each hold both resistant and susceptible isolates; otherwise the report states why it is unavailable. A predictive k-mer can mark a lineage or a plasmid that travels with the cause; holding out lineages limits this association without removing it.
      </Def>
      <Def term="Balanced accuracy">The mean of sensitivity and specificity, over pooled out-of-fold predictions at a 0.5 threshold.</Def>
      <Def term="Null fits">{nullFitWords(manifest.evaluation_revision)}</Def>
      <Def term="VME and ME">
        Very major error is resistant called susceptible, over all resistant isolates; major error is susceptible called resistant, over all susceptible. Judged against FDA-style bars for susceptibility test devices, VME at most 1.5% and ME at most 3%. These are reference bars, not a regulatory claim.
      </Def>
      <Def term="Cannot-call rates">
        Not reported in this version. A conformal cannot-call rate needs a calibration set held apart from training, and the pipeline does not fit one yet; earlier reports computed them on pooled out-of-fold scores, so the dashboard does not show those.
      </Def>
      <Def term="MIC prediction">
        Censored gradient-boosted regression on log2 dilution steps, fitted on the MICs of genomes with a reported call; a genome whose label was inferred from its MIC is left out. Essential agreement is a prediction within one doubling dilution of an exact MIC. Categorical calls from the predicted MIC use each training fold&apos;s own empirical threshold and are scored against the reported resistant and susceptible labels.
      </Def>
      <Def term="Simulated PCR panels">
        A gene panel calls resistant when the assembled genome carries any of its targets, which follow the manufacturers&apos; published target lists; a probe assay calls resistant when a probed region differs from the wild-type sequence. Detection is in silico: target presence in genomes, not a PCR assay&apos;s analytical performance.
      </Def>
      <Def term="Feature annotation">
        A feature's k-mers are placed on a reference assembly that carries them and mapped to its gene annotation: the gene and codon range, or the nearest gene when the k-mers sit between genes, with the flanking sequence.
      </Def>
    </dl>
  </Section>
);

const EvaluationAudit = ({ manifest }: { manifest: Manifest }) => {
  const audit = manifestAudit(manifest);
  return (
  <Section id="m-audit" title="Evaluation audit">
      <p className="prose-measure text-note text-ink">
        {audit.total === 0
          ? "No dataset has a report yet."
          : audit.invalid > 0
            ? `${int(audit.invalid)} of ${plural(audit.total, "dataset evaluation")} do not pass the evaluation audit, so their metrics are provisional: ${audit.reasons.join(" ")}`
            : `Every dataset evaluation passes the internal evaluation audit. ${audit.scopes[0] ?? "That covers internal cross-validation against reported lab calls."}`}{" "}
        {audit.invalid > 0 ? "Passing audits cover internal cross-validation on these public genomes, not external or clinical validation." : ""}
      </p>
  </Section>
  );
};

export const MethodsView = ({ manifest }: { manifest: Manifest }) => {
  const corpus = useCorpus();
  const unsupportedRows = manifest.unsupported_panels ?? [];
  const unsupported = [...new Set(unsupportedRows.map((u) => u.panel))].map((panel) => {
    const datasets = unsupportedRows.filter((u) => u.panel === panel);
    return { ...datasets[0]!, datasets };
  });
  return (
    <div className={`${WIDE_SHELL} pb-16 pt-8`}>
      <h1 className="text-display font-semibold text-ink">Methods</h1>
      <div className="mt-8 space-y-10">
        <EvaluationAudit manifest={manifest} />
        <Funnel manifest={manifest} />
        {corpus.status === "loading" ? <PanelState kind="loading" title="Loading the corpus file" /> : null}
        {corpus.status === "error" ? <PanelState kind="error" title="The corpus file could not be read" message={loadErrorText(corpus.error)} /> : null}
        {corpus.status === "ready" ? (
          <>
            <Names corpus={corpus.data} manifest={manifest} />
            <Breakpoints corpus={corpus.data} />
            <Sources corpus={corpus.data} />
          </>
        ) : null}
        <Design manifest={manifest} />
        {unsupported.length > 0 ? (
          <Section id="m-unsupported" title="Assays not scored" answer={`${plural(unsupported.length, "assay")} ${unsupported.length === 1 ? "has" : "have"} no result on the datasets listed, so no sensitivity or specificity is reported for ${unsupported.length === 1 ? "it" : "them"}.`}>
            {unsupported.map((a) => (
              <article key={a.panel} className="border-t border-hair py-4">
                <h3 className="text-body font-semibold text-ink">{a.name}</h3>
                <p className="text-note text-muted-text">{a.datasets.map((d) => `${capitalize(d.antibiotic)} in ${speciesShort(d.species)}`).join(", ")}</p>
                <p className="mt-2 text-body text-ink">{notScored(a.reason)}</p>
              </article>
            ))}
          </Section>
        ) : null}
      </div>
    </div>
  );
};
