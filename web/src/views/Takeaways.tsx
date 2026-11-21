// Executive summary at the top of Overview: one plain bottom-line lead,
// six scannable takeaway cards that link to their full sections below,
// the supports / does-not-support boundary, and a four-chapter jump bar.
// Every dynamic number is derived from the manifest, optional baselines
// summary, and bundled reference JSON files at render time.

import type { Manifest } from "../data/schema";
import { int, pct, plural, points } from "../data/format";
import { micBaselineSummary, type MicConstant } from "../data/baselines";
import { catalogSummary, micSummary, splitSummary } from "../data/derive/overview";
import { barSummary, learnerSummary, ruleBands, wordingSummary } from "../data/derive/evidence";
import { armById, inputSummary } from "../data/derive/inputs";
import { populationById, populationSummary } from "../data/derive/outside";
import { studySplitSummary } from "../data/derive/studySplit";
import { otherDatabaseSummary } from "../data/derive/otherDatabase";
import { FOCUS_RING } from "../ui/primitives";
import { LINK, showSection } from "./findings";
import { ScopeFigure } from "./primer";

export const TAKEAWAYS_ANCHOR = "takeaways";

export interface TakeawayCard {
  id: string;
  kicker: string;
  stat: string;
  body: string;
  target: string;
  linkLabel: string;
}

const CHAPTERS = [
  { id: "corpus", label: "Cohort & core scores" },
  { id: "catalog", label: "Rules, panels & MIC" },
  { id: "outside", label: "Studies & another database" },
  { id: "inputs", label: "Inputs & limits" },
] as const;

const median = (xs: readonly number[]): number | null => {
  if (xs.length === 0) return null;
  const s = [...xs].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 === 1 ? s[m]! : (s[m - 1]! + s[m]!) / 2;
};

/** The bottom-line summary lead above the six cards. */
export const bottomLineLead = (m: Manifest): string => {
  const s = splitSummary(m.datasets);
  const b = barSummary(m.datasets, "best", "lineage");
  const st = studySplitSummary();
  const randomAll = median(m.datasets.map((d) => d.best_random).filter((v): v is number => v !== null));
  if (m.datasets.length === 0) return "No dataset has a report yet.";
  const lineagePart =
    s.nLineage > 0 && s.medianLineage !== null
      ? ` and ${pct(s.medianLineage)} with whole lineages held out (${plural(s.nLineage, "dataset")})`
      : "";
  const screenPart =
    b.n > 0 ? ` With lineages held out, ${int(b.meets)} of ${plural(b.n, "dataset")} are shown to meet both limits of the adapted research screen.` : "";
  return `On public BV-BRC genomes, the best of five k-mer learners reaches a median balanced accuracy of ${pct(randomAll)} on random splits (${plural(m.datasets.length, "dataset")})${lineagePart}.${screenPart} Folds that hold whole studies out score worse than random folds for ${int(st.costShown)} of 5 learners. Every figure is retrospective. No external or clinical validation has been done.`;
};

/** Builds the six concise takeaway cards (1–2 sentences each), leaving full detail in each linked section. */
export const takeawayCards = (m: Manifest, constants?: ReadonlyMap<string, MicConstant>): TakeawayCard[] => {
  const rows = m.datasets;
  const s = splitSummary(rows);
  const b = barSummary(rows, "best", "lineage");
  const w = wordingSummary(rows, "best", "lineage");
  const ls = learnerSummary(rows);
  const scmRow = ls.rows.find((r) => r.learner === "scm");
  const cartRow = ls.rows.find((r) => r.learner === "cart");
  const rfRow = ls.rows.find((r) => r.learner === "rf");
  const cat = catalogSummary(rows);
  const bands = ruleBands(rows);
  const highBand = bands[bands.length - 1];
  const mic = micSummary(rows);
  const micBase = constants === undefined ? null : micBaselineSummary(rows, constants);
  const st = studySplitSummary();
  const outPrimary = populationSummary(populationById("primary"));
  const od = otherDatabaseSummary();
  const inp = inputSummary();
  const cap1000 = armById("cap-1000").all;
  const budget1000 = armById("budget-1000").all;

  const lineageVmes = rows
    .map((r) => r.learners.find((l) => l.learner === r.best_learner)?.lineage?.vme)
    .filter((v): v is number => v !== null && v !== undefined);
  const lineageMes = rows
    .map((r) => r.learners.find((l) => l.learner === r.best_learner)?.lineage?.me)
    .filter((v): v is number => v !== null && v !== undefined);
  const medVme = median(lineageVmes);
  const medMe = median(lineageMes);

  return [
    {
      id: "shift",
      kicker: "Unseen lineages and studies",
      stat: s.medianGap === null ? "–" : `${points(s.medianGap)} median loss`,
      body:
        s.nLineage === 0 || s.medianGap === null
          ? `In a retrospective split by study inside BV-BRC (${int(st.datasets)} datasets, ${int(st.species)} species), folds keeping BioProjects whole score worse than random folds for ${int(st.costShown)} of 5 learners.`
          : `Holding whole lineages out costs the best learner a median ${points(s.medianGap)} on ${int(s.nLineage)} datasets (${int(s.nLargeGap)} lose more than 5 points). Holding whole BioProjects out in BV-BRC (${int(st.datasets)} datasets, ${int(st.species)} species) worsens log loss and specificity for ${int(st.logLossWorse)} of 5 learners, and on ${int(outPrimary.frozen)} newly deposited unseen-study genomes SCM holds on ${int(outPrimary.holds.scm.holds)} of ${int(outPrimary.holds.scm.read)} datasets.`,
      target: "study-split",
      linkLabel: "See lineage and study splits",
    },
    {
      id: "screen",
      kicker: "Adapted research error screen",
      stat: b.n === 0 ? "–" : `${int(b.meets)} of ${int(b.n)} meet both`,
      body:
        b.n === 0
          ? "No dataset has both error rates with lineages held out yet."
          : `Against VME at most 1.5% and ME at most 3.0% by exact 95% intervals, the best learner with lineages held out is shown to meet both on ${int(b.meets)} of ${int(b.n)} datasets and shown to miss at least one on ${int(b.misses)} (${int(b.open)} not decided${medVme !== null && medMe !== null ? `; median VME ${pct(medVme)}, ME ${pct(medMe)}` : ""}). By the looser acceptance wording of FDA's 2009 guidance for growth-based susceptibility test systems (not an applicable criterion here), ${int(w.meets)} of ${int(w.n)} meet.`,
      target: "bars",
      linkLabel: "See error-rate screen",
    },
    {
      id: "learners",
      kicker: "Learners and inputs",
      stat: scmRow?.medianLineage == null ? "–" : `${pct(scmRow.medianLineage)} SCM alone`,
      body:
        scmRow?.medianLineage == null || cartRow?.medianLineage == null || rfRow?.medianLineage == null
          ? `In the 4-species input comparison (${int(inp.answered)} of ${int(inp.datasets)} datasets answered), no changed input was shown better than the default.`
          : `With lineages held out, median balanced accuracy is ${pct(scmRow.medianLineage)} for SCM and ${pct(cartRow.medianLineage)} for the decision tree against ${pct(rfRow.medianLineage)} for the random forest (${int(scmRow.n)} datasets each). Across ${int(inp.species)} species (${int(inp.answered)} of ${int(inp.datasets)} datasets answered), no input change was shown better than the default, while 1,000 genomes (${int(cap1000.worse)} of ${int(cap1000.datasets)}) and 1,000 features (${int(budget1000.worse)} of ${int(budget1000.datasets)}) read worse.`,
      target: "learners",
      linkLabel: "See learners and inputs",
    },
    {
      id: "known-gene",
      kicker: "Known-gene rules",
      stat: cat.nCompared === 0 ? "–" : `${int(cat.within.get(0.02) ?? 0)} of ${int(cat.nCompared)} within 2 pts`,
      body:
        cat.nCompared === 0
          ? "No dataset has both a known-gene rule and a lineage-held-out result yet."
          : `Fitted on nothing, the known-gene rule is ahead of the best learner or within 2 points with lineages held out on ${int(cat.within.get(0.02) ?? 0)} of ${int(cat.nCompared)} datasets${highBand && highBand.compared > 0 ? ` (${int(highBand.ruleHolds)} of ${int(highBand.compared)} where the rule itself scores 85% or more)` : ""}. The learner is ahead by more than 2 points on ${int(cat.regions.behind)} of ${int(cat.nCompared)}.`,
      target: "catalog",
      linkLabel: "See known-gene comparison",
    },
    {
      id: "other-db",
      kicker: "Another public database (NCBI)",
      stat: `${int(od.totals.scm.plan_meets)} of ${int(od.eightSpecies)} species (SCM)`,
      body: `Scored unchanged on ${int(od.totalGenomes)} NCBI Pathogen Detection genomes outside BV-BRC's labelled set (not external validation), SCM meets the plan's rule in ${int(od.totals.scm.plan_meets)} of ${int(od.eightSpecies)} species (${int(od.totals.scm.holds)} of ${int(od.totals.scm.estimable)} datasets; E. coli pilot: SCM ${int(od.pilot.scm.holds)} of ${int(od.pilot.scm.estimable)}, CART ${int(od.pilot.cart.holds)} of ${int(od.pilot.cart.estimable)}), CART in ${int(od.totals.cart.plan_meets)} of ${int(od.eightSpecies)}, and lasso in ${int(od.totals.lasso.plan_meets)} of ${int(od.totals.lasso.assessable)} species where it can be read (${int(od.totals.lasso.holds)} of its ${int(od.totals.lasso.estimable)} datasets across the eight). The rule reads recall and log loss only, and S. enterica and C. coli share BioProjects with training.`,
      target: "other-database",
      linkLabel: "See NCBI database test",
    },
    {
      id: "mic-data",
      kicker: "Quantitative MIC",
      stat: mic.n === 0 ? "–" : `${pct(mic.median)} within 1 dilution`,
      body:
        mic.n === 0
          ? "No dataset had enough MIC values for MIC regression."
          : micBase !== null && micBase.covered > 0
            ? `MIC models land within one doubling dilution on a median ${pct(micBase.medianModel)} of exact MICs (${int(micBase.covered)} datasets), against ${pct(micBase.medianConstant)} for a constant training-fold median MIC; the model is not above the constant on ${int(micBase.notAbove)} of ${int(micBase.covered)}. MIC folds are random, with no lineage split.`
            : `MIC models land within one doubling dilution on a median ${pct(mic.median)} of exact MICs across ${int(mic.n)} datasets (${int(mic.nAbove90)} reach 90% or more). MIC folds are random, with no lineage split.`,
      target: "mic",
      linkLabel: "See MIC results",
    },
  ];
};

export const Takeaways = ({ manifest, micConstants }: { manifest: Manifest; micConstants?: ReadonlyMap<string, MicConstant> }) => {
  const cards = takeawayCards(manifest, micConstants);
  return (
    <section id={TAKEAWAYS_ANCHOR} aria-labelledby={`${TAKEAWAYS_ANCHOR}-title`} className="scroll-mt-6 border-t border-hair pt-7">
      <div className="flex flex-wrap items-baseline justify-between gap-x-6 gap-y-2">
        <h2 id={`${TAKEAWAYS_ANCHOR}-title`} tabIndex={-1} className={`rounded-sm text-headline font-semibold tracking-tight text-ink ${FOCUS_RING}`}>
          At a glance
        </h2>
        <nav aria-label="Report chapters" className="flex flex-wrap items-center gap-x-4 gap-y-1 text-control">
          {CHAPTERS.map((ch) => (
            <a
              key={ch.id}
              href={`#${ch.id}`}
              onClick={(e) => {
                e.preventDefault();
                showSection(ch.id, true);
              }}
              className={`inline-flex min-h-11 items-center font-medium ${LINK}`}
            >
              {ch.label}
            </a>
          ))}
        </nav>
      </div>

      <p className="prose-measure mt-2 text-body text-ink">{bottomLineLead(manifest)}</p>

      <ol className="mt-6 grid grid-cols-1 gap-4 md:grid-cols-2 lg:grid-cols-3">
        {cards.map((c) => (
          <li key={c.id} data-takeaway-card={c.id} className="flex min-w-0 flex-col rounded-md border border-hair bg-paper-pure p-4">
            <span className="text-tiny font-medium uppercase tracking-wider text-muted-text">{c.kicker}</span>
            <span className="mt-1 text-headline font-semibold tabular-nums text-ink">{c.stat}</span>
            <p className="mt-2 text-note text-ink">{c.body}</p>
            <a
              href={`#${c.target}`}
              onClick={(e) => {
                e.preventDefault();
                showSection(c.target, true);
              }}
              className={`mt-auto inline-flex min-h-11 items-center gap-1 pt-2 text-control font-medium ${LINK}`}
            >
              {c.linkLabel}
              <span aria-hidden="true">↓</span>
            </a>
          </li>
        ))}
      </ol>

      <div className="mt-6 rounded-md border border-hair bg-paper-pure p-4">
        <ScopeFigure />
      </div>
    </section>
  );
};
