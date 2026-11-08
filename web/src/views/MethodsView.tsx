// Data & Methods: how the analysis works and where it can fail. Six plain
// steps open the view, in the order the work happens, each linking to its
// detail; then four linear sections. Which labels count and which genomes
// are in, with the curation funnel Overview's counts link to; the models and
// how they are evaluated; how to read each score; and the scope, limits, and
// audit of the evidence. Each section and definition opens with a plain line
// and keeps its technical detail after it. Every count comes from the
// manifest or the corpus file at render time; learner settings and fitted
// details stay in each dataset's report.

import { useEffect, useState, type ReactNode } from "react";
import type { Corpus, Manifest } from "../data/schema";
import { datasetTitle, int, pct, plural } from "../data/format";
import { BARS, barWords } from "../data/metrics";
import { loadErrorText, useCorpus } from "../data/load";
import { manifestAudit, nullFitWords } from "../data/derive/evaluation";
import { classMinimum } from "../data/derive/overview";
import { notScored } from "../data/derive/tests";
import { FOCUS_RING, PanelState, ScrollFrame, VIEW_SHELL } from "../ui/primitives";
import { ViewHeader } from "../ui/ViewHeader";
import { YearColumns } from "../charts/YearColumns";
import { minimumSentence, ReportLine, SectionHeading, selectionCap, showSection, STATUS_LINK } from "./findings";

/** The anchors other views link to. */
export const METHODS_ANCHORS = { funnel: "m-funnel", audit: "m-audit", limits: "m-limits" } as const;

const P = ({ children }: { children: ReactNode }) => <p className="prose-measure text-body text-ink">{children}</p>;

/** Technical detail after a section's plain lead, at the definitions' size. */
const Detail = ({ children }: { children: ReactNode }) => <p className="prose-measure text-note text-ink">{children}</p>;

const MethodsSection = ({ id, n, title, children }: { id: string; n: number; title: string; children: ReactNode }) => (
  <section id={id} aria-labelledby={`${id}-title`} className="scroll-mt-6 border-t border-hair pt-7">
    <SectionHeading id={id} n={n}>
      {title}
    </SectionHeading>
    <div className="mt-4 space-y-4">{children}</div>
  </section>
);

/** A part of a section, under its own heading; a link can land on it. */
const Part = ({ id, title, children }: { id?: string; title: string; children: ReactNode }) => (
  <div id={id} className="scroll-mt-6 space-y-3 pt-2">
    <h3 tabIndex={id === undefined ? undefined : -1} className={`rounded-sm text-body font-semibold text-ink ${FOCUS_RING}`}>
      {title}
    </h3>
    {children}
  </div>
);

/** A term: its plain meaning in one line, then its role in a sentence or
 *  two of technical detail. */
const Def = ({ term, plain, children }: { term: string; plain: string; children: ReactNode }) => (
  <div className="border-b border-hair py-3">
    <dt className="text-control font-semibold text-ink">{term}</dt>
    <dd className="prose-measure mt-1">
      <p className="text-body text-ink">{plain}</p>
      <p className="mt-1 text-note text-ink">{children}</p>
    </dd>
  </div>
);

const Defs = ({ children }: { children: ReactNode }) => <dl className="grid grid-cols-1 gap-x-12 border-t border-hair lg:grid-cols-2">{children}</dl>;

// ---- 1 labels and inclusion -----------------------------------------------------------

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

/** The curation funnel: pairs, genomes, and the records set aside before
 *  pairing, each table with its sum in words. */
const Funnel = ({ manifest }: { manifest: Manifest }) => {
  const f = manifest.funnel;
  const failed = Object.values(f.genomes_failed_qc).reduce((a, b) => a + b, 0);
  const genomesAdd = f.genomes_without_metadata + failed + f.genomes_passing_qc === f.genomes_with_records;
  return (
    <Part id={METHODS_ANCHORS.funnel} title="Curation funnel">
      <P>
        These tables follow the public records down to the genomes and antibiotics the models use. In the first two tables, the top row goes on to the models, with a blue
        bar on wider screens. Every other row is set aside, and its label says why. QC (quality control) checks each genome&apos;s assembled DNA, and each reason a genome failed has its own row.
      </P>
      <div className="grid grid-cols-1 gap-x-12 gap-y-6 xl:grid-cols-2">
        <div id="m-pairs" className="xl:col-span-2">
          <h4 className="mb-2 text-control font-semibold text-ink">Genome–antibiotic pairs</h4>
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
          <h4 className="mb-2 text-control font-semibold text-ink">Genomes and QC outcomes</h4>
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
          <h4 className="mb-2 text-control font-semibold text-ink">Records set aside before pairing</h4>
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
    </Part>
  );
};

const Labels = ({ manifest }: { manifest: Manifest }) => {
  const { cap, capped } = selectionCap(manifest);
  return (
    <MethodsSection id="m-labels" n={1} title="Labels and inclusion">
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
        Every learner trains and is scored on genomes with an unconflicted reported resistant (R) or susceptible (S) lab call for the drug, and the known-gene rules, panels, and MIC
        regression are scored on the same genomes. A genome and antibiotic pair with both an R and an S call is a conflict and gets no label; a pair with only intermediate calls is never
        an R or S label.
      </Detail>
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
      <Funnel manifest={manifest} />
    </MethodsSection>
  );
};

// ---- 2 models and evaluation -----------------------------------------------------------

const Models = ({ manifest }: { manifest: Manifest }) => (
  <MethodsSection id="m-models" n={2} title="Models and evaluation">
    <P>
      A learner is a method that studies genomes whose label is known and builds a model: a rule that predicts the label from the DNA alone. GRM reads each genome as
      yes-or-no answers (does it carry this piece of DNA?) and gives the same answers to every learner, so the learners can be compared. GRM compares five learners with a
      fixed rule built from a catalog of known resistance genes.
    </P>
    <P>
      Every score comes from genomes the model never trained on. A random split hides a random fifth of the genomes, so most hidden genomes still have close relatives
      (the same lineage) in training. The lineage split hides whole lineages, so no hidden genome has a close relative in training. The drop between the two shows how
      much a model relies on recognising relatives.
    </P>
    <Defs>
      <Def term="K-mer features" plain="Short overlapping pieces of DNA (k-mers), each one a yes-or-no question: does this genome carry it?">
        {manifest.k === null
          ? "Each genome is counted into k-mers with GRM's own counter. The published reports do not all record the same k, so each dataset's Method and coverage fold states its own."
          : `Each genome is counted into its ${manifest.k}-mers with GRM's own counter.`}{" "}
        K-mers with the same presence pattern across a species&apos; genomes are one feature, kept when enough isolates carry it.
      </Def>
      <Def term="Set covering machine (SCM)" plain="A learner that writes a few short yes-or-no rules a person can read.">A short set of presence and absence rules, chosen greedily as Kover does; it calls resistant when its rules hold and is read as written.</Def>
      <Def term="Decision tree (CART)" plain="A learner that asks one yes-or-no question after another, like a flowchart.">A class-balanced tree of presence and absence splits.</Def>
      <Def term="Random forest" plain="A learner that grows many different decision trees and lets them vote.">Many randomized trees whose votes score each isolate between 0 and 1.</Def>
      <Def term="Gradient-boosted trees" plain="A learner that adds small trees one at a time, each one fixing mistakes of the trees before it.">Trees added one at a time, each correcting the ones before; they score each isolate between 0 and 1.</Def>
      <Def term="L1 logistic regression (lasso)" plain="A learner that gives a weight to a small number of DNA pieces and ignores the rest.">A sparse linear model on screened features.</Def>
      <Def term="Known-gene rules" plain="A fixed rule with no training: call the genome resistant when it carries a known resistance gene that acts on the drug.">
        GRM finds acquired genes of the AMRFinderPlus catalog ({manifest.amrfinder_version === "" ? "version not recorded" : manifest.amrfinder_version}) in each genome with its own k-mer
        matching and applies its own fixed rule per drug: resistant when the genome carries a gene of a class or subclass that acts on the drug. The rule does not train. It does not run
        the AMRFinderPlus program or test the point mutations AMRFinderPlus also describes, so mutational resistance counts against it.
      </Def>
      <Def term="Simulated PCR panels" plain="Commercial lab gene tests, imitated on the genome sequence: a panel calls resistant when it finds any gene it looks for.">
        A gene panel calls resistant when the assembled genome carries any of its targets, which follow the manufacturers&apos; published target lists; a probe assay calls resistant when a
        probed region differs from the wild-type sequence. Detection is in silico, from assembled genomes.
      </Def>
      <Def term="MIC regression" plain="A model that predicts the MIC itself, the lowest amount of the drug that stops growth, instead of only R or S.">
        Censored gradient-boosted regression on log2 dilution steps, fitted on the MICs of genomes with a reported call; a genome whose label was inferred from its MIC is left out. It
        needs at least 100 isolates with an MIC, 50 of them exact.
      </Def>
      <Def term="Random split" plain="GRM cuts the genomes into five random parts (folds) and holds out each part once while the model trains on the other four.">Stratified 5-fold cross-validation, the scheme most published numbers use.</Def>
      <Def term="Lineage-held-out split" plain="The same test, but GRM always holds out close relatives (one lineage) together, so no held-out genome has a relative in training.">
        Whole lineages, clusters of genomes within a fixed Mash distance, go to one fold, so no test isolate has a near-identical relative in training. It runs when a dataset has at least
        10 lineages and the folds can each hold resistant and susceptible isolates; otherwise the report says why. The drop from the random split estimates how much a score depends on
        close relatives, for this sampled population. It does not remove every shortcut: a predictive k-mer can still mark a plasmid or a lineage that travels with the cause, and
        lineages absent from BV-BRC are not tested.
      </Def>
    </Defs>
    <p className="type-note">Each dataset&apos;s report gives its learner settings, the fitted rules and trees, and how its lineages were split.</p>
  </MethodsSection>
);

// ---- 3 how to read scores ----------------------------------------------------------------

const Scores = ({ manifest }: { manifest: Manifest }) => (
  <MethodsSection id="m-scores" n={3} title="How to read scores">
    <P>
      Every score compares a model&apos;s calls with the lab&apos;s calls on genomes the model never trained on. A call is right when it matches the lab. GRM counts the two
      kinds of wrong call separately, because calling a resistant genome susceptible is the more serious mistake.
    </P>
    <Defs>
      <Def term="Balanced accuracy" plain="The average of two shares: resistant genomes called resistant (sensitivity) and susceptible genomes called susceptible (specificity). 50% is no better than a coin toss.">The mean of sensitivity and specificity, over pooled out-of-fold predictions. A score of 0.5 or more calls resistant.</Def>
      <Def term="Best learner" plain="For each dataset, the learner that scored highest on the random split.">
        The learner with the highest random-split balanced accuracy among those that call every isolate on every split that ran. The known-gene rule is a baseline, never a candidate.
        Overview&apos;s figures use the best learner.
      </Def>
      <Def
        term="VME and ME"
        plain={`The two kinds of wrong call: a resistant genome called susceptible (very major error) and a susceptible genome called resistant (major error). Test devices are held to ${barWords(BARS.vme)} and ${barWords(BARS.me)}.`}
      >
        Very major error (VME) is resistant isolates called susceptible, over the resistant isolates called; major error (ME) is susceptible isolates called resistant, over the
        susceptible isolates called. They are read against FDA-style bars for susceptibility test devices, VME at most 1.5% and ME at most 3%, as a reference, not a regulatory claim.
      </Def>
      <Def term="Learner no-calls" plain="What happens when a fitted model learned nothing from the DNA, and how its genomes are counted.">{nullFitWords(manifest.evaluation_revision)}</Def>
      <Def term="Probe no-calls" plain="A probe test that finds none of its targets gives no answer, and that genome is left out of the test's own score.">
        A probe assay gives no call for an isolate whose genome carries none of its probe k-mers. That isolate is neither found nor missed: the panel&apos;s sensitivity and specificity
        cover the isolates with a probe call, and the learner&apos;s shares still count it.
      </Def>
      <Def
        term="MIC essential agreement"
        plain={`How often the predicted MIC is within one dilution step of the measured MIC. Test devices are held to ${barWords(BARS.mic)}.`}
      >
        A prediction within one doubling dilution of an exact MIC, a single measured value rather than a bound such as “≤0.25”; 90% is the bar used for susceptibility test devices.
        Calls from the predicted MIC use each training fold&apos;s own empirical threshold and are scored against the reported resistant and susceptible labels.
      </Def>
      <Def term="AUROC" plain="How well a model's scores rank resistant genomes above susceptible ones: 0.5 is chance, and 1 means every resistant genome ranks above every susceptible one.">
        The area under the receiver operating characteristic (ROC) curve, printed to three decimals in each dataset&apos;s learner table. A dash means the report has no
        AUROC for that learner.
      </Def>
      <Def term="What each Overview figure counts" plain="Which datasets each Overview figure includes. A median is the middle value when the datasets are lined up by score.">
        The lineage figure draws datasets with both splits; the known-gene figure, datasets with a rule and a lineage-held-out result; the PCR panel, datasets with that panel and a best
        learner; the MIC figure, datasets with an MIC model. Datasets with no best learner are left out of every median, and each figure states what it leaves out.
      </Def>
    </Defs>
  </MethodsSection>
);

// ---- 4 scope, limitations and audit -------------------------------------------------------

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

/** Where the labeled isolates come from: country, year, laboratory method,
 *  and testing standard, the sampling every score inherits. */
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
          <h4 className="mb-2 text-control font-semibold text-ink">Labeled genomes by isolation country</h4>
          <BarList rows={countries} label="Countries" />
          {corpus.countries.length > countries.length ? <p className="type-note mt-2">Top {countries.length} of {int(corpus.countries.length)} countries.</p> : null}
        </div>
        <div>
          <h4 className="mb-2 text-control font-semibold text-ink">Labeled genomes by collection year</h4>
          <YearColumns years={years} />
        </div>
        <div>
          <h4 className="mb-2 text-control font-semibold text-ink">Lab records by testing method</h4>
          <BarList rows={corpus.methods} label="Methods" />
        </div>
        <div>
          <h4 className="mb-2 text-control font-semibold text-ink">Lab records by testing standard</h4>
          <BarList rows={corpus.standards} label="Standards" />
        </div>
      </div>
    </>
  );
};

const NO_SCOPE = "Internal cross-validation against reported lab calls; not external or clinical validation.";

const Scope = ({ manifest }: { manifest: Manifest }) => {
  const corpus = useCorpus();
  const audit = manifestAudit(manifest);
  const { cap, capped } = selectionCap(manifest);
  const unsupportedRows = manifest.unsupported_panels ?? [];
  const unsupported = [...new Set(unsupportedRows.map((u) => u.panel))].map((panel) => {
    const datasets = unsupportedRows.filter((u) => u.panel === panel);
    return { ...datasets[0]!, datasets };
  });
  return (
    <MethodsSection id={METHODS_ANCHORS.limits} n={4} title="Scope, limitations and audit">
      <P>
        Every number in this report measures agreement with lab calls already stored in BV-BRC. The models were never tried on new samples, in another lab, or in a
        clinic.
      </P>
      <Part id={METHODS_ANCHORS.audit} title="Evaluation audit">
        <P>The evaluation audit is a check, run on each dataset&apos;s report, that its scores were computed under the current evaluation rules.</P>
        <P>
          {audit.total === 0
            ? "No dataset has a report yet."
            : audit.invalid > 0
              ? `${int(audit.invalid)} of ${plural(audit.total, "dataset evaluation")} do not pass the evaluation audit, so their metrics are provisional.`
              : `Every dataset evaluation passes the internal evaluation audit${manifest.evaluation_revision === undefined ? "" : ` (revision ${manifest.evaluation_revision})`}.`}
        </P>
        {audit.total === 0 ? null : audit.invalid > 0 ? (
          <ul className="prose-measure list-disc space-y-0.5 pl-5 text-note text-ink">
            {audit.reasons.map((r) => (
              <li key={r}>{r}</li>
            ))}
          </ul>
        ) : (
          <p className="prose-measure text-note text-ink">{audit.scopes[0] ?? NO_SCOPE}</p>
        )}
      </Part>
      <Part title="Why this is not clinical validation">
        <P>
          Every score is internal cross-validation: it measures agreement with the calls BV-BRC reports, on genomes from the same public collection. The genomes are a selection
          {capped ? `, at most ${int(cap)} per species in this run,` : ""} labeled by many laboratories with different methods and testing standards across years, while breakpoints
          changed over time. No external, prospective, or clinical evaluation was done, so no number here is the performance of a diagnostic test.
        </P>
      </Part>
      <Part id="m-sources" title="Where the isolates come from">
        <P>
          An isolate is one sample of bacteria grown in a lab from one source, and each isolate here has one genome. Where and when the isolates were collected shapes
          every score.
        </P>
        {corpus.status === "loading" ? <PanelState kind="loading" title="Loading the corpus file" /> : null}
        {corpus.status === "error" ? <PanelState kind="error" title="The corpus file could not be read" message={loadErrorText(corpus.error)} /> : null}
        {corpus.status === "ready" ? <Sources corpus={corpus.data} /> : null}
      </Part>
      <Part title="Learner cross-check">
        <P>GRM runs its own code for each learner, so GRM checked each one against a standard version of the same method.</P>
        <p className="prose-measure text-note text-ink">
          On seven small test fixtures and on the 868-isolate <i>E. coli</i> ceftriaxone dataset from an earlier report, run with the same learner code as this report, each of the five
          learners matched a reference implementation of the same model (scikit-learn, XGBoost, or a reimplementation of its documented rules) exactly or within that reference&apos;s own
          run-to-run noise.
        </p>
      </Part>
      {unsupported.length > 0 ? (
        <Part id="m-unsupported" title="Assays not scored">
          <p className="prose-measure text-note text-ink">
            {plural(unsupported.length, "assay")} {unsupported.length === 1 ? "has" : "have"} no result on the datasets listed, so no sensitivity or specificity is reported for{" "}
            {unsupported.length === 1 ? "it" : "them"}.
          </p>
          <ul className="space-y-3">
            {unsupported.map((a) => (
              <li key={a.panel} className="border-t border-hair pt-3">
                <p className="text-control font-semibold text-ink">{a.name}</p>
                <p className="text-note text-muted-text">{a.datasets.map(datasetTitle).join(", ")}</p>
                <p className="prose-measure mt-1 text-note text-ink">{notScored(a.reason)}</p>
              </li>
            ))}
          </ul>
        </Part>
      ) : null}
    </MethodsSection>
  );
};

// ---- the study in six steps ------------------------------------------------------------

/** How the study works, in the order the work happens: plain words only,
 *  each step linking to the section that holds its detail. */
const STUDY_STEPS: readonly { title: string; text: string; section: string; to: string }[] = [
  {
    title: "Collect genomes and lab results",
    text: "Genomes and their lab results come from BV-BRC, a public database. For each genome, a lab grew the bacterium next to one antibiotic and reported it resistant (R) if it kept growing or susceptible (S) if the drug stopped it.",
    section: "section 1, Labels and inclusion",
    to: "m-labels",
  },
  {
    title: "Train learners",
    text: "A learner studies genomes whose lab result is known and builds a model: a rule that predicts R or S from the DNA alone. GRM compares five learners with a fixed rule built from a catalog of known resistance genes.",
    section: "section 2, Models and evaluation",
    to: "m-models",
  },
  {
    title: "Test on genomes held out",
    text: "Every score comes from genomes the model never trained on. A random split holds out a random fifth of the genomes. A lineage split holds out whole families of close relatives, so a model cannot score well by recognising a relative.",
    section: "section 2, Models and evaluation",
    to: "m-models",
  },
  {
    title: "Score the calls",
    text: "Balanced accuracy is the average of the share of resistant genomes called resistant and the share of susceptible genomes called susceptible. 50% is no better than a coin toss.",
    section: "section 3, How to read scores",
    to: "m-scores",
  },
  {
    title: "Count the errors against limits",
    text: `A resistant genome called susceptible is a very major error (VME), and a susceptible genome called resistant is a major error (ME). Susceptibility test devices are held to a VME of ${barWords(BARS.vme)} and an ME of ${barWords(BARS.me)}.`,
    section: "section 3, How to read scores",
    to: "m-scores",
  },
  {
    title: "Read what the numbers support",
    text: "Every number is agreement with lab results in one public database. No external, prospective or clinical validation has been done.",
    section: "section 4, Scope, limitations and audit",
    to: METHODS_ANCHORS.limits,
  },
];

const StudySteps = () => (
  <ol aria-label="How the study works" className="mt-8 grid grid-cols-1 gap-x-12 gap-y-5 lg:grid-cols-2">
    {STUDY_STEPS.map((step, i) => (
      <li key={step.title} className="grid grid-cols-[1.5rem_minmax(0,1fr)] gap-x-2">
        <span aria-hidden="true" className="text-body font-semibold tabular-nums text-muted-text">
          {i + 1}
        </span>
        <div className="prose-measure min-w-0">
          <p className="text-body font-semibold text-ink">{step.title}</p>
          <p className="mt-1 text-body text-ink">{step.text}</p>
          <a
            href={`#${step.to}`}
            onClick={(e) => {
              e.preventDefault();
              showSection(step.to, true);
            }}
            className={`mt-1 inline-block text-note text-muted-text ${STATUS_LINK}`}
          >
            More in {step.section}
          </a>
        </div>
      </li>
    ))}
  </ol>
);

export const MethodsView = ({ manifest }: { manifest: Manifest }) => {
  // Read once per mount: a link from another view lands on one section.
  const [hash] = useState(() => (typeof window === "undefined" ? "" : decodeURIComponent(window.location.hash.slice(1))));
  useEffect(() => {
    if (hash !== "") showSection(hash);
  }, [hash]);
  return (
    <div className={`${VIEW_SHELL} pb-16 pt-8`}>
      <ViewHeader
        title="Data & Methods"
        purpose="How GRM turns public genomes and lab results into the scores in this report, and what those scores can and cannot show."
        status={
          <ReportLine
            manifest={manifest}
            scope={
              <a
                href={`#${METHODS_ANCHORS.limits}`}
                onClick={(e) => {
                  e.preventDefault();
                  showSection(METHODS_ANCHORS.limits, true);
                }}
                className={STATUS_LINK}
              >
                Internal cross-validation only
              </a>
            }
          />
        }
      />
      <StudySteps />
      <div className="mt-12 space-y-12">
        <Labels manifest={manifest} />
        <Models manifest={manifest} />
        <Scores manifest={manifest} />
        <Scope manifest={manifest} />
      </div>
    </div>
  );
};
