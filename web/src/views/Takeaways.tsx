// Above-the-fold forest-first summary at the top of Overview: one plain
// bottom-line sentence, a boxless bolster KPI strip (4 headline numbers with
// arrow deltas and denominators), six concise 1-sentence takeaway cards, a
// single supports / does-not-support boundary, and a 4-chapter jump bar.
// Every dynamic number is derived from the manifest, optional baselines
// summary, and bundled reference JSON files at render time.

import { useId, useState, type MouseEvent } from "react";
import { ChevronDown } from "lucide-react";
import type { Manifest } from "../data/schema";
import { int, pct, plural, points } from "../data/format";
import { micBaselineSummary, type MicConstant } from "../data/baselines";
import { catalogSummary, micSummary, speciesCount, splitSummary } from "../data/derive/overview";
import { barSummary, learnerSummary, ruleBands, wordingSummary } from "../data/derive/evidence";
import { armById, inputSummary } from "../data/derive/inputs";
import { populationById, populationSummary } from "../data/derive/outside";
import { studySplitSummary } from "../data/derive/studySplit";
import { otherDatabaseSummary } from "../data/derive/otherDatabase";
import { ADDITIONAL_STUDIES } from "../data/derive/additionalStudies";
import { BolsterStatStrip, type BolsterStat } from "../charts/amr";
import { FOCUS_RING } from "../ui/primitives";
import { LINK, selectionCap, showSection } from "./findings";
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
  { id: "corpus", label: "Cohort" },
  { id: "catalog", label: "Rules & MIC" },
  { id: "outside", label: "Studies & tools" },
  { id: "inputs", label: "Inputs & limits" },
] as const;

const median = (xs: readonly number[]): number | null => {
  if (xs.length === 0) return null;
  const s = [...xs].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 === 1 ? s[m]! : (s[m - 1]! + s[m]!) / 2;
};

/** The bottom-line summary lead above the bolster strip and cards. */
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
  return `On public BV-BRC genomes, the best of five k-mer learners reaches a median balanced accuracy of ${pct(randomAll)} on random splits (${plural(m.datasets.length, "dataset")})${lineagePart}.${screenPart} Folds that hold whole studies out score worse than random folds on log loss and specificity for ${int(st.logLossWorse)} of 5 learners. Every figure is retrospective. No external or clinical validation has been done.`;
};

/** Builds the four boxless headline stats for the bolster strip. */
const bolsterStats = (m: Manifest, constants?: ReadonlyMap<string, MicConstant>): BolsterStat[] => {
  const s = splitSummary(m.datasets);
  const b = barSummary(m.datasets, "best", "lineage");
  const st = studySplitSummary();
  const mic = micSummary(m.datasets);
  const micBase = constants === undefined ? null : micBaselineSummary(m.datasets, constants);
  const jump = (target: string) => (e: MouseEvent<HTMLAnchorElement>) => {
    e.preventDefault();
    showSection(target, true);
  };
  return [
    {
      key: "lineage",
      label: "Lineages held out",
      value: s.medianLineage === null ? "–" : pct(s.medianLineage),
      delta: s.medianGap === null ? undefined : `▼ ${(Math.abs(s.medianGap) * 100).toFixed(1)} pts`,
      tone: s.medianGap !== null && s.medianGap >= 0.05 ? "coral" : "ink",
      sub:
        s.nLineage === 0
          ? "No lineage split yet"
          : `from ${pct(s.medianRandom)} random · paired median drop · ${int(s.nLineage)} of ${int(m.datasets.length)} datasets`,
      href: "#lineage",
      onJump: jump("lineage"),
    },
    {
      key: "screen",
      label: "Adapted research error screen",
      value: b.n === 0 ? "–" : `${int(b.meets)} of ${int(b.n)}`,
      tone: b.n > 0 && b.meets === 0 ? "coral" : "ink",
      sub: b.n === 0 ? "No error rates yet" : `${int(b.misses)} miss ≥1 limit · ${int(b.open)} not decided`,
      href: "#bars",
      onJump: jump("bars"),
    },
    {
      key: "study",
      label: "Study splits",
      value: `${int(st.logLossWorse)} of 5 worse`,
      tone: "coral",
      sub: `on log loss & specificity · ${int(st.datasets)} datasets in ${int(st.species)} species`,
      href: "#study-split",
      onJump: jump("study-split"),
    },
    {
      key: "mic",
      label: "MIC (random folds)",
      value: mic.n === 0 ? "–" : pct(micBase !== null && micBase.covered > 0 ? micBase.medianModel : mic.median),
      sub:
        mic.n === 0
          ? "No MIC models yet"
          : micBase !== null && micBase.covered > 0
            ? `vs ${pct(micBase.medianConstant)} training-fold median · ${int(micBase.covered)} datasets`
            : `within 1 dilution · ${int(mic.n)} datasets`,
      href: "#mic",
      onJump: jump("mic"),
    },
  ];
};

/** Builds six concise takeaway cards, keeping every result caveat and 22 Nov wording intact. */
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
  const extra = ADDITIONAL_STUDIES;

  return [
    {
      id: "shift",
      kicker: "Unseen lineages, studies & dates",
      stat: s.medianGap === null ? "–" : `▼ ${(Math.abs(s.medianGap) * 100).toFixed(1)} pts median drop`,
      body:
        s.nLineage === 0 || s.medianGap === null
          ? `In a retrospective split by study inside BV-BRC (${int(st.datasets)} datasets, ${int(st.species)} species), folds keeping BioProjects whole score worse than random folds on log loss and specificity for ${int(st.logLossWorse)} of 5 learners.`
          : `Holding whole lineages out costs the best learner a median ${points(s.medianGap)} on ${int(s.nLineage)} datasets (${int(s.nLargeGap)} lose more than 5 points). Holding whole BioProjects out (${int(st.datasets)} datasets) or training only on earlier public dates (${int(extra.time_split.datasets)} datasets) worsens log loss for ${int(st.logLossWorse)} of 5 learners (recall worse for RF and GBT on study splits, and for ${int(extra.time_split.total_recall_worse)} of 5 on date splits); on ${int(outPrimary.frozen)} newly deposited unseen-study genomes SCM holds on ${int(outPrimary.holds.scm.holds)} of ${int(outPrimary.holds.scm.read)}.`,
      target: "study-split",
      linkLabel: "Study and date splits",
    },
    {
      id: "screen",
      kicker: "Adapted research error screen",
      stat: b.n === 0 ? "–" : `${int(b.meets)} of ${int(b.n)} meet both`,
      body:
        b.n === 0
          ? "No dataset has both error rates with lineages held out yet."
          : `Against VME ≤1.5% and ME ≤3.0% by exact 95% intervals, ${int(b.misses)} of ${int(b.n)} miss at least one limit and ${int(b.open)} are not decided; under the looser acceptance wording of FDA's 2009 guidance for growth-based susceptibility test systems (not an applicable criterion here), ${int(w.meets)} of ${int(w.n)} meet.`,
      target: "bars",
      linkLabel: "Error screen details",
    },
    {
      id: "learners",
      kicker: "Learners & inputs",
      stat: scmRow?.medianLineage == null ? "–" : `${pct(scmRow.medianLineage)} rule set`,
      body:
        scmRow?.medianLineage == null || cartRow?.medianLineage == null || rfRow?.medianLineage == null
          ? `Across ${int(inp.answered)} of ${int(inp.datasets)} datasets, no single input change beat the default.`
          : `With lineages held out, median balanced accuracy is ${pct(scmRow.medianLineage)} for the rule set (SCM) and ${pct(cartRow.medianLineage)} for the decision tree against ${pct(rfRow.medianLineage)} for the random forest (${int(scmRow.n)} datasets). Neither a dense neural network (${int(extra.dense_network.datasets)} datasets: ${extra.dense_network.outcome}) nor any single input change across ${int(inp.species)} species beat the default, while 1,000 genomes (${int(cap1000.worse)} of ${int(cap1000.datasets)}) and 1,000 features (${int(budget1000.worse)} of ${int(budget1000.datasets)}) read worse.`,
      target: "learners",
      linkLabel: "Learners and inputs",
    },
    {
      id: "known-gene",
      kicker: "Known genes & tools",
      stat: cat.nCompared === 0 ? "–" : `${int(cat.within.get(0.02) ?? 0)} of ${int(cat.nCompared)} within 2 pts`,
      body:
        cat.nCompared === 0
          ? "No dataset has both a known-gene rule and a lineage-held-out result yet."
          : `Fitted on nothing, GRM's known-gene rule is ahead or within 2 points of the best learner on ${int(cat.within.get(0.02) ?? 0)} of ${int(cat.nCompared)} datasets${highBand && highBand.compared > 0 ? ` (${int(highBand.ruleHolds)} of ${int(highBand.compared)} where the rule scores ≥85%)` : ""}. The full AMRFinderPlus 4.2.7 program leads GRM's rule on ${int(extra.amrfinderplus.program_vs_known_gene_rule.ahead)} of ${int(extra.amrfinderplus.datasets)} datasets and roughly ties the 5 learners on lineage splits, while TBProfiler 6.6.3 leads SCM on ${int(extra.tbprofiler.graded_vs_learners_lineage.scm.tool_higher)} of ${int(extra.tbprofiler.graded_endpoints)} graded M. tuberculosis drugs on assemblies.`,
      target: "comparisons",
      linkLabel: "Other comparisons",
    },
    {
      id: "other-db",
      kicker: "NCBI database",
      stat: `${int(od.totals.scm.plan_meets)} of ${int(od.eightSpecies)} species`,
      body: `Scored unchanged on ${int(od.totalGenomes)} NCBI Pathogen Detection genomes outside BV-BRC's labelled set (not external validation), SCM meets the plan's rule in ${int(od.totals.scm.plan_meets)} of ${int(od.eightSpecies)} species (${int(od.totals.scm.holds)} of ${int(od.totals.scm.estimable)} datasets; E. coli pilot: SCM ${int(od.pilot.scm.holds)} of ${int(od.pilot.scm.estimable)}, CART ${int(od.pilot.cart.holds)} of ${int(od.pilot.cart.estimable)}), CART in ${int(od.totals.cart.plan_meets)} of ${int(od.eightSpecies)}, and lasso in ${int(od.totals.lasso.plan_meets)} of ${int(od.totals.lasso.assessable)} species where it can be read (${int(od.totals.lasso.holds)} of its ${int(od.totals.lasso.estimable)} datasets across the eight). The rule reads recall and log loss only, and S. enterica and C. coli share BioProjects with training.`,
      target: "other-database",
      linkLabel: "NCBI database test",
    },
    {
      id: "mic-data",
      kicker: "MIC",
      stat: mic.n === 0 ? "–" : `${pct(mic.median)} within 1 dilution`,
      body:
        mic.n === 0
          ? "No dataset had enough MIC values for MIC regression."
          : micBase !== null && micBase.covered > 0
            ? `MIC models land within one doubling dilution on a median ${pct(micBase.medianModel)} of exact MICs (${int(micBase.covered)} datasets), against ${pct(micBase.medianConstant)} for a constant training-fold median MIC; the model is not above the constant on ${int(micBase.notAbove)} of ${int(micBase.covered)}. MIC folds are random, with no lineage split.`
            : `MIC models land within one doubling dilution on a median ${pct(mic.median)} of exact MICs across ${int(mic.n)} datasets (${int(mic.nAbove90)} reach 90% or more). MIC folds are random, with no lineage split.`,
      target: "mic",
      linkLabel: "MIC results",
    },
  ];
};

/** Expandable plain-English walkthrough of what GRM built, tested, and learned. */
const StartToFinishGuide = ({ manifest }: { manifest: Manifest }) => {
  const [open, setOpen] = useState(false);
  const panelId = useId();
  const analyzed = manifest.species.reduce((a, s) => a + s.genomes_analyzed, 0);
  const s = splitSummary(manifest.datasets);
  const b = barSummary(manifest.datasets, "best", "lineage");
  const st = studySplitSummary();
  const od = otherDatabaseSummary();
  const extra = ADDITIONAL_STUDIES;
  const { cap, capped } = selectionCap(manifest);

  return (
    <div className="mt-4">
      <button
        type="button"
        aria-expanded={open}
        aria-controls={open ? panelId : undefined}
        onClick={() => setOpen((v) => !v)}
        className={`inline-flex min-h-11 items-center gap-1.5 rounded-sm text-control font-medium text-blue transition-colors duration-fast ease-out hover:text-ink ${FOCUS_RING}`}
      >
        {open ? "Hide the guide" : "Read the guide"}
        <ChevronDown size={14} strokeWidth={1.75} aria-hidden="true" className={`transition-transform duration-fast ease-out ${open ? "rotate-180" : ""}`} />
      </button>
      {open ? (
        <div id={panelId} data-start-to-finish="" className="mt-3 grid grid-cols-1 gap-5 border-y border-hair py-5 md:grid-cols-2">
          <div>
            <h3 className="text-control font-semibold text-ink">1. Engine</h3>
            <p className="mt-1.5 text-note text-ink">
              GRM is a self-contained Rust engine that slices each bacterial genome into 31-letter DNA chunks (k-mers), merges identical presence patterns into features, and trains
              five models: a set covering machine (short yes-or-no rules), decision tree, random forest, boosted trees, and L1 logistic regression (lasso), beside a zero-training
              known-gene rule, simulated PCR panels, and quantitative MIC models.
            </p>
          </div>
          <div>
            <h3 className="text-control font-semibold text-ink">2. Data and splits</h3>
            <p className="mt-1.5 text-note text-ink">
              The main release covers {plural(manifest.datasets.length, "dataset")} across {plural(speciesCount(manifest), "species", "species")} using {int(analyzed)} genomes
              selected from {int(manifest.funnel.genomes_passing_qc)} QC-passing public BV-BRC genomes
              {capped ? ` (capped at ${int(cap)} per species in this run)` : ""}. Beyond random cross-validation, GRM refit the same five learners on unseen bacterial lineages (
              {int(s.nLineage)} datasets), held-out BioProject studies ({int(st.datasets)} datasets), and earlier-vs-later public dates ({int(extra.time_split.datasets)} datasets),
              and scored the released models unchanged on newly deposited BV-BRC genomes and {int(od.totalGenomes)} NCBI Pathogen Detection genomes outside BV-BRC&apos;s labelled
              set (not external validation).
            </p>
          </div>
          <div>
            <h3 className="text-control font-semibold text-ink">3. Findings</h3>
            <p className="mt-1.5 text-note text-ink">
              On lineage splits, SCM and CART had the highest medians (paired gaps from SCM to CART and to lasso do not clear zero). In separate comparisons, GRM&apos;s SCM led
              Kover 2.0.7 on {int(extra.kover.datasets)} small datasets, and a one-setting dense neural network established no gain over the five learners across{" "}
              {int(extra.dense_network.datasets)} datasets. Every harder split lowered scores: {int(b.meets)} of {int(b.n)} lineage-held-out datasets were shown to meet both
              limits of the adapted research error screen, holding out whole studies or earlier dates worsened log loss for all 5 learners (with recall worse for RF and GBT on
              study splits, and for {int(extra.time_split.total_recall_worse)} of 5 on date splits), AMRFinderPlus 4.2.7 led GRM&apos;s simplified known-gene rule on{" "}
              {int(extra.amrfinderplus.program_vs_known_gene_rule.ahead)} of {int(extra.amrfinderplus.datasets)} datasets while roughly tying the 5 lineage learners, and TBProfiler
              6.6.3 scored higher than SCM on {int(extra.tbprofiler.graded_vs_learners_lineage.scm.tool_higher)} of {int(extra.tbprofiler.graded_endpoints)} graded TB drugs on
              assemblies.
            </p>
          </div>
          <div>
            <h3 className="text-control font-semibold text-ink">4. Limits and navigation</h3>
            <p className="mt-1.5 text-note text-ink">
              Every figure is retrospective on public data, never clinical validation. Use the chapter links above to jump straight to any evidence section, open the 5-step visual
              guide right below for plate-level definitions, or switch to Datasets, Species, Data, or Methods in the top bar.
            </p>
          </div>
        </div>
      ) : null}
    </div>
  );
};

export const Takeaways = ({ manifest, micConstants }: { manifest: Manifest; micConstants?: ReadonlyMap<string, MicConstant> }) => {
  const cards = takeawayCards(manifest, micConstants);
  const stats = bolsterStats(manifest, micConstants);
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

      <div className="mt-5">
        <BolsterStatStrip stats={stats} ariaLabel="Headline results" />
      </div>

      <StartToFinishGuide manifest={manifest} />

      <ol className="mt-5 grid grid-cols-1 gap-x-8 gap-y-5 md:grid-cols-2 lg:grid-cols-3">
        {cards.map((c) => (
          <li key={c.id} data-takeaway-card={c.id} className="flex min-w-0 flex-col border-t border-hair pt-3">
            <span className="text-tiny font-medium text-muted-text">{c.kicker}</span>
            <span className="mt-0.5 text-headline font-semibold tabular-nums text-ink">{c.stat}</span>
            <p className="mt-1.5 text-note text-ink">{c.body}</p>
            <a
              href={`#${c.target}`}
              onClick={(e) => {
                e.preventDefault();
                showSection(c.target, true);
              }}
              className={`mt-auto inline-flex min-h-11 items-center gap-1 pt-1 text-control font-medium ${LINK}`}
            >
              {c.linkLabel}
              <span aria-hidden="true">↓</span>
            </a>
          </li>
        ))}
      </ol>

      <div className="mt-6 border-t border-hair pt-5">
        <ScopeFigure />
      </div>
    </section>
  );
};
