// The primer: five short steps a reader with no genomics or statistics
// background needs before the findings. What a dataset is, what a genome
// feature is, what the two test splits mean, how to read one result and its
// errors, and what the numbers can and cannot support. One real report
// runs through the steps; its numbers come from the manifest at render time,
// and the two drawings that are not data say so. Overview shows the steps
// as one section above the findings; Present shows one step per slide.

import type { ReactNode } from "react";
import { Check, X } from "lucide-react";
import type { CatalogRow, LearnerLite, Manifest } from "../data/schema";
import { capitalize, int, pct, speciesShort } from "../data/format";
import { BARS, barLabel } from "../data/metrics";
import { routeHref, type Navigator } from "../app/router";
import { FOCUS_RING } from "../ui/primitives";
import { inApp, LINK, showSection } from "./findings";

/** The report the primer reads from: the best random-split learner's scores
 *  on both splits, with every rate the steps print. */
export interface PrimerExample {
  row: CatalogRow;
  random: LearnerLite & { balanced_accuracy: number };
  lineage: LearnerLite & { balanced_accuracy: number; sensitivity: number; specificity: number };
}

/** The dataset the brief's walkthroughs used when this report has it. */
const PREFERRED = ["escherichia_coli__ceftriaxone", "staphylococcus_aureus__clindamycin"];

const usable = (row: CatalogRow): PrimerExample | null => {
  const best = row.learners.find((l) => l.learner === row.best_learner);
  const random = best?.random;
  const lineage = best?.lineage;
  if (random == null || lineage == null || random.balanced_accuracy === null) return null;
  if (lineage.balanced_accuracy === null || lineage.sensitivity === null || lineage.specificity === null) return null;
  return { row, random: { ...random, balanced_accuracy: random.balanced_accuracy }, lineage: { ...lineage, balanced_accuracy: lineage.balanced_accuracy, sensitivity: lineage.sensitivity, specificity: lineage.specificity } };
};

/** A preferred dataset when it has both splits, else the largest that does;
 *  null when no report has a lineage-held-out result. */
export const primerExample = (m: Manifest): PrimerExample | null => {
  for (const id of PREFERRED) {
    const row = m.datasets.find((d) => d.id === id);
    const ex = row === undefined ? null : usable(row);
    if (ex !== null) return ex;
  }
  return (
    [...m.datasets]
      .sort((a, b) => b.n - a.n || a.id.localeCompare(b.id))
      .map(usable)
      .find((e): e is PrimerExample => e !== null) ?? null
  );
};

/** "Ceftriaxone in E. coli". */
export const exampleName = (row: CatalogRow): string => `${capitalize(row.antibiotic)} in ${speciesShort(row.species)}`;

/** "ceftriaxone in E. coli", for mid-sentence. */
export const exampleNameInline = (row: CatalogRow): string => `${row.antibiotic} in ${speciesShort(row.species)}`;

/** The lineage split's very major errors as one sentence: which learner,
 *  which drug and species, and that resistant genomes were called
 *  susceptible. */
export const missedSentence = (ex: PrimerExample): string =>
  `With lineages held out, the best learner for ${exampleNameInline(ex.row)} called ${perHundred(ex).missed} of every 100 resistant genomes susceptible.`;

/** Per 100, rounded: how many resistant genomes were called susceptible and
 *  how many susceptible genomes were called resistant, on the lineage split. */
export const perHundred = (ex: PrimerExample): { missed: number; caught: number; flagged: number; cleared: number } => {
  const missed = Math.round((1 - ex.lineage.sensitivity) * 100);
  const flagged = Math.round((1 - ex.lineage.specificity) * 100);
  return { missed, caught: 100 - missed, flagged, cleared: 100 - flagged };
};

export interface PrimerStep {
  id: string;
  title: string;
  /** A slide's one-sentence headline. */
  lead: string;
  /** The page's fuller text: a few short sentences. */
  body: ReactNode;
  figure: ReactNode;
  /** What the figure shows, in a sentence or two. */
  caption: string;
}

// ---- figures ------------------------------------------------------------------

const FIG_LABEL = "text-control text-ink";
const FIG_NOTE = "text-control text-muted-text";

/** The example's genomes as one bar: resistant, then susceptible, each
 *  segment named with its count beneath, so colour is never the only cue. */
const CompositionFigure = ({ row }: { row: CatalogRow }) => {
  const s = row.n - row.n_resistant;
  const r = row.n === 0 ? 0 : row.n_resistant / row.n;
  return (
    <figure data-primer-figure="dataset" className="m-0 w-full max-w-[26rem]">
      <p className={FIG_NOTE}>{exampleName(row)}</p>
      <div className="mt-1.5 flex h-4 w-full gap-0.5" aria-hidden="true">
        <span className="h-full rounded-l-sm" style={{ width: `${r * 100}%`, background: "var(--series-r)" }} />
        <span className="h-full rounded-r-sm" style={{ width: `${(1 - r) * 100}%`, background: "var(--series-s)" }} />
      </div>
      <figcaption className="mt-1.5 flex justify-between gap-4">
        <span className={FIG_LABEL}>
          <span className="font-semibold tabular-nums">{int(row.n_resistant)}</span> resistant (R)
        </span>
        <span className={`${FIG_LABEL} text-right`}>
          <span className="font-semibold tabular-nums">{int(s)}</span> susceptible (S)
        </span>
      </figcaption>
    </figure>
  );
};

/** A drawing, not data: one 31-letter piece read from a stretch of DNA, and
 *  a yes-or-no table where each column is one feature. */
const SEQ = "TTGACCGATGCTGACGGCAGCGGTCAGCTGGATCGCTTAGCCAT";
const K_START = 6;
const K = 31;
const TABLE: { genome: string; call: "R" | "S"; has: [boolean, boolean, boolean] }[] = [
  { genome: "Genome 1", call: "R", has: [true, true, false] },
  { genome: "Genome 2", call: "R", has: [true, false, true] },
  { genome: "Genome 3", call: "S", has: [false, true, true] },
  { genome: "Genome 4", call: "S", has: [false, false, true] },
];

const FeatureFigure = () => (
  <figure data-primer-figure="feature" className="m-0 w-full max-w-[26rem]">
    <p className={FIG_NOTE}>Drawing, not data</p>
    <p className="mt-1.5 overflow-hidden whitespace-nowrap font-mono text-control tracking-normal text-muted-text" aria-hidden="true">
      {SEQ.slice(0, K_START)}
      <span className="rounded-sm bg-blue/[0.12] text-ink underline decoration-blue decoration-2 underline-offset-[3px]">{SEQ.slice(K_START, K_START + K)}</span>
      {SEQ.slice(K_START + K)}
    </p>
    <p className={`mt-1 ${FIG_NOTE}`}>One 31-letter piece (k-mer), underlined</p>
    <table className="mt-3 border-collapse text-control">
      <thead>
        <tr>
          <th scope="col" className="pr-4 text-left font-medium text-muted-text">
            <span className="sr-only">Genome</span>
          </th>
          {["Feature A", "B", "C"].map((h) => (
            <th key={h} scope="col" className="px-2 text-center font-medium text-muted-text">
              {h}
            </th>
          ))}
          <th scope="col" className="pl-4 text-left font-medium text-muted-text">
            Lab call
          </th>
        </tr>
      </thead>
      <tbody>
        {TABLE.map((t) => (
          <tr key={t.genome} className="border-t border-hair">
            <th scope="row" className="py-0.5 pr-4 text-left font-normal text-ink">
              {t.genome}
            </th>
            {t.has.map((h, i) => (
              <td key={i} className={`px-2 py-0.5 text-center ${h ? "font-semibold text-ink" : "text-muted-text"}`}>
                {h ? "yes" : "no"}
              </td>
            ))}
            <td className="py-0.5 pl-4 font-medium text-ink">{t.call}</td>
          </tr>
        ))}
      </tbody>
    </table>
  </figure>
);

/** A drawing of the two splits: five lineages of related genomes; filled
 *  squares are hidden from training. Random hides a scattered fifth, the
 *  lineage split one whole lineage. The example's scores sit beside each. */
const LINEAGES = [6, 4, 5, 3, 4];
const RANDOM_HIDDEN = new Set([2, 7, 12, 17, 20]);
const LINEAGE_HIDDEN = 2;

/** One split: its name and score on one line, then the lineages as one row
 *  of groups on fixed tracks, so both rows line up group for group at any
 *  width that holds them (about 290px); narrower, whole groups wrap. */
const SplitRow = ({ label, hidden, score }: { label: string; hidden: (lineage: number, index: number) => boolean; score?: string }) => {
  let index = 0;
  return (
    <div data-split-row="">
      <div className="flex items-baseline justify-between gap-x-4">
        <p className={FIG_LABEL}>{label}</p>
        {score === undefined ? null : <span className="text-headline font-semibold tabular-nums text-ink">{score}</span>}
      </div>
      <div className="mt-1 flex flex-wrap gap-x-2 gap-y-1" aria-hidden="true">
        {LINEAGES.map((size, l) => (
          <span key={l} data-lineage-group="" className="flex gap-0.5 border-b border-hair pb-1">
            {Array.from({ length: size }, () => {
              const i = index++;
              const h = hidden(l, i);
              return <span key={i} data-hidden={h ? "" : undefined} className={`block size-2.5 rounded-[2px] ${h ? "bg-ink" : "border border-ink/40"}`} />;
            })}
          </span>
        ))}
      </div>
    </div>
  );
};

const SplitFigure = ({ ex }: { ex: PrimerExample | null }) => (
  <figure data-primer-figure="splits" className="m-0 w-full max-w-[20rem] space-y-3">
    <SplitRow label="Random split" hidden={(_l, i) => RANDOM_HIDDEN.has(i)} score={ex === null ? undefined : pct(ex.random.balanced_accuracy)} />
    <SplitRow label="Lineages held out" hidden={(l) => l === LINEAGE_HIDDEN} score={ex === null ? undefined : pct(ex.lineage.balanced_accuracy)} />
    <figcaption className={`flex flex-wrap gap-x-4 gap-y-1 ${FIG_NOTE}`}>
      <span className="inline-flex items-center gap-1.5">
        <span aria-hidden="true" data-hidden-swatch="" className="block size-2.5 rounded-[2px] bg-ink" /> Hidden from training
      </span>
      <span className="inline-flex items-center gap-1.5">
        <span aria-hidden="true" className="block size-2.5 rounded-[2px] border border-ink/40" /> Trained on
      </span>
      <span>Groups are lineages{ex === null ? "" : `; scores are ${exampleNameInline(ex.row)}`}</span>
    </figcaption>
  </figure>
);

/** 100 genomes as a 10 by 10 grid; the errors are filled and come first, so
 *  their count reads by fill and position as well as colour. A filled square
 *  wears the colour of the call it got, as in the dataset bar. */
const Hundred = ({ errors, label, words, color }: { errors: number; label: string; words: string; color: string }) => (
  <div className="min-w-0">
    <p className={FIG_LABEL}>{label}</p>
    <div className="mt-1.5 grid w-max grid-cols-10 gap-[3px]" aria-hidden="true">
      {Array.from({ length: 100 }, (_, i) => (
        <span
          key={i}
          data-error={i < errors ? "" : undefined}
          className={`block size-[9px] rounded-[2px] ${i < errors ? "" : "border border-ink/25"}`}
          style={i < errors ? { background: color } : undefined}
        />
      ))}
    </div>
    <p className={`mt-1.5 ${FIG_LABEL}`}>
      <span className="font-semibold tabular-nums">{errors}</span> {words}
    </p>
  </div>
);

const ErrorFigure = ({ ex }: { ex: PrimerExample }) => {
  const h = perHundred(ex);
  return (
    <figure data-primer-figure="errors" className="m-0 flex flex-wrap gap-x-8 gap-y-4">
      <Hundred errors={h.missed} label="100 resistant genomes" words="called susceptible (VME)" color="var(--series-s)" />
      <Hundred errors={h.flagged} label="100 susceptible genomes" words="called resistant (ME)" color="var(--series-r)" />
    </figure>
  );
};

const SUPPORTS = [
  "Where genomes agree with reported lab calls in public BV-BRC data, and where they do not",
  "Which datasets lean on close relatives to score well",
  "Which markers and datasets are worth follow-up study",
];
const NOT_SUPPORTED = [
  "Choosing a treatment for a patient",
  "Performance in other labs, hospitals, test methods or years",
  "That a marked gene or k-mer causes resistance",
];

const ScopeFigure = () => (
  <figure data-primer-figure="scope" className="m-0 grid w-full max-w-[40rem] grid-cols-1 gap-x-8 gap-y-3 sm:grid-cols-2">
    <div>
      <p className="text-control font-semibold text-ink">Supports</p>
      <ul className="mt-1 space-y-1">
        {SUPPORTS.map((t) => (
          <li key={t} className="flex gap-2 text-control text-ink">
            <Check aria-hidden="true" size={14} strokeWidth={2.25} className="mt-[3px] shrink-0 text-live-text" />
            {t}
          </li>
        ))}
      </ul>
    </div>
    <div>
      <p className="text-control font-semibold text-ink">Does not support</p>
      <ul className="mt-1 space-y-1">
        {NOT_SUPPORTED.map((t) => (
          <li key={t} className="flex gap-2 text-control text-ink">
            <X aria-hidden="true" size={14} strokeWidth={2.25} className="mt-[3px] shrink-0 text-coral-text" />
            {t}
          </li>
        ))}
      </ul>
    </div>
  </figure>
);

// ---- steps --------------------------------------------------------------------

/** The five steps, with the example's numbers where the report has one. */
export const primerSteps = (m: Manifest, ex: PrimerExample | null = primerExample(m)): PrimerStep[] => {
  const name = ex === null ? null : exampleName(ex.row);
  const rule = ex?.row.scm_rules[0];
  const h = ex === null ? null : perHundred(ex);
  return [
    {
      id: "primer-dataset",
      title: "A dataset",
      lead:
        ex === null
          ? "A dataset is one species tested against one antibiotic."
          : `A dataset is one species tested against one antibiotic: ${name}, ${int(ex.row.n)} genomes, ${int(ex.row.n_resistant)} resistant.`,
      body: (
        <>
          <p>
            A lab grows the bacterium next to the antibiotic. If the drug stops it growing, the lab reports it susceptible (S); if it keeps growing, resistant (R). A dataset
            is every selected genome of one species with a reported R or S call for one antibiotic, and each dataset has one report.
          </p>
          {ex === null ? null : (
            <p className="mt-2">
              {name} has {int(ex.row.n)} genomes: {int(ex.row.n_resistant)} resistant and {int(ex.row.n - ex.row.n_resistant)} susceptible.
            </p>
          )}
        </>
      ),
      figure: ex === null ? null : <CompositionFigure row={ex.row} />,
      caption: "R and S are the calls the lab reported to BV-BRC. Every score in this report is agreement with those calls.",
    },
    {
      id: "primer-feature",
      title: "A genome feature",
      lead: "GRM reads each genome as yes-or-no answers: does it carry this piece of DNA?",
      body: (
        <>
          <p>
            GRM cuts every genome into overlapping 31-letter pieces of DNA (k-mers). Pieces found in exactly the same genomes merge into one feature, a yes-or-no column. A
            model learns which columns separate resistant from susceptible genomes.
          </p>
          {ex === null ? null : (
            <p className="mt-2">
              {name} has {int(ex.row.features)} features.
              {rule === undefined ? "" : ` The set covering machine's first rule is “${rule}”.`}
            </p>
          )}
        </>
      ),
      figure: <FeatureFigure />,
      caption: "A feature can travel with the cause of resistance, on the same plasmid or in the same lineage, without being the cause.",
    },
    {
      id: "primer-splits",
      title: "Two tests",
      lead:
        ex === null
          ? "Each model is scored twice: on a random split, and with whole lineages held out."
          : `${name} scores ${pct(ex.random.balanced_accuracy)} on random splits and ${pct(ex.lineage.balanced_accuracy)} with whole lineages held out.`,
      body: (
        <>
          <p>
            Every score comes from genomes the model never trained on. A random split hides a random fifth of the genomes, so most hidden genomes still have close relatives
            (the same lineage) in training. The lineage split hides whole lineages, so no hidden genome has a close relative in training. The drop between the two shows how
            much a model relies on recognising relatives.
          </p>
          <p className="mt-2">
            The score is balanced accuracy: the average of the share of resistant genomes called resistant (sensitivity) and the share of susceptible genomes called
            susceptible (specificity). 50% is no better than a coin toss.
          </p>
        </>
      ),
      figure: <SplitFigure ex={ex} />,
      caption: "Squares are genomes, grouped by lineage; the layout is a drawing. Both tests use the same public database, so neither is an external test.",
    },
    {
      id: "primer-errors",
      title: "Reading one result",
      lead:
        ex === null
          ? "Two kinds of error matter: a resistant genome called susceptible, and the reverse."
          : missedSentence(ex),
      body: (
        <>
          {h === null || ex === null ? null : (
            <p className="mb-2">
              With lineages held out, the best learner for {exampleNameInline(ex.row)} called {h.caught} of every 100 resistant genomes resistant and {h.cleared} of every 100
              susceptible genomes susceptible.
            </p>
          )}
          <p>
            A resistant genome called susceptible is a very major error (VME): in use, it would mean giving a drug that fails. A susceptible genome called resistant is a
            major error (ME): it would mean giving a broader drug than needed. Susceptibility test devices are held to a VME of {barLabel(BARS.vme)} and an ME of{" "}
            {barLabel(BARS.me)}; the Datasets table marks each value against those bars.
          </p>
        </>
      ),
      figure: ex === null ? null : <ErrorFigure ex={ex} />,
      caption: `Counts per 100, rounded, from the best random-split learner with lineages held out. Bars for test devices: VME ${barLabel(BARS.vme)}, ME ${barLabel(BARS.me)}.`,
    },
    {
      id: "primer-scope",
      title: "What the numbers support",
      lead: "These are research results on public data, not a validated test.",
      body: (
        <p>
          Every number is internal cross-validation on genomes and lab calls from one public database. No external, prospective or clinical validation has been done, and a
          marked feature shows association, not cause.
        </p>
      ),
      figure: <ScopeFigure />,
      caption: "Internal cross-validation only: both test splits come from the same BV-BRC records.",
    },
  ];
};

// ---- the Overview section ---------------------------------------------------------

/** The anchor of the primer on Overview. */
export const PRIMER_ANCHOR = "primer";

/** The primer as one section above the findings: a heading, the example
 *  report it reads from, then each step as text beside its figure. A link
 *  skips ahead to the findings for a reader who already knows the terms. */
export const PrimerSection = ({ manifest, nav, firstFinding }: { manifest: Manifest; nav: Navigator; firstFinding: string }) => {
  const ex = primerExample(manifest);
  const steps = primerSteps(manifest, ex);
  return (
    <section id={PRIMER_ANCHOR} aria-labelledby={`${PRIMER_ANCHOR}-title`} className="scroll-mt-6 border-t border-hair pt-7">
      <div className="flex flex-wrap items-baseline justify-between gap-x-6 gap-y-1">
        <h2 id={`${PRIMER_ANCHOR}-title`} tabIndex={-1} className={`rounded-sm text-headline font-semibold tracking-tight text-ink ${FOCUS_RING}`}>
          How to read this report
        </h2>
        <a
          href={`#${firstFinding}`}
          onClick={(e) => {
            e.preventDefault();
            showSection(firstFinding, true);
          }}
          className={`inline-flex min-h-11 items-center gap-1.5 text-control font-medium ${LINK}`}
        >
          Skip to the findings
          <span aria-hidden="true">↓</span>
        </a>
      </div>
      <p className="prose-measure mt-2 text-body text-ink">
        Five steps, no background needed.
        {ex === null ? null : (
          <>
            {" "}
            The example is one real report,{" "}
            <a href={routeHref({ ...nav.route, dataset: ex.row.id })} onClick={inApp(() => nav.openDataset(ex.row.id))} className={LINK}>
              {exampleName(ex.row)}
            </a>
            .
          </>
        )}
      </p>
      <ol className="mt-6">
        {steps.map((s) => (
          <li
            key={s.id}
            id={s.id}
            className="grid grid-cols-1 gap-x-10 gap-y-4 border-t border-hair py-6 first:border-t-0 first:pt-0 lg:grid-cols-[minmax(0,34rem)_minmax(0,1fr)]"
          >
            <div className="min-w-0">
              <h3 id={`${s.id}-title`} className="text-body font-semibold text-ink">
                {s.title}
              </h3>
              <div className="prose-measure mt-1.5 text-body text-ink">{s.body}</div>
            </div>
            {s.figure === null ? null : (
              <div className="min-w-0 self-center">
                {s.figure}
                <p className="prose-measure mt-3 text-note text-muted-text">{s.caption}</p>
              </div>
            )}
          </li>
        ))}
      </ol>
    </section>
  );
};
