// Above-the-fold forest-first summary at the top of Overview: one plain
// bottom-line sentence, a boxless bolster KPI strip (4 headline numbers with
// arrow deltas and denominators), six concise 1-sentence takeaway cards, a
// single supports / does-not-support boundary, and a 4-chapter jump bar.
// Every dynamic number is derived from the manifest, optional baselines
// summary, and bundled reference JSON files at render time.

import { useId, useState, type MouseEvent } from "react";
import { ChevronDown } from "lucide-react";
import type { Manifest } from "../data/schema";
import { int, pct, plural } from "../data/format";
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
      ? ` to ${pct(s.medianLineage)} on unseen lineages (${plural(s.nLineage, "dataset")})`
      : "";
  const screenPart =
    b.n > 0 ? `, with ${int(b.meets)} of ${plural(b.n, "dataset")} shown to meet both error limits` : "";
  return `Median balanced accuracy drops from ${pct(randomAll)} on random splits (${plural(m.datasets.length, "dataset")})${lineagePart}${screenPart}, and holding whole studies out worsens ${int(st.logLossWorse)} of 5 learners. Every figure is retrospective. No external or clinical validation has been done.`;
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
      label: "Lineage drop",
      value: s.medianLineage === null ? "–" : pct(s.medianLineage),
      delta: s.medianGap === null ? undefined : `▼ ${(Math.abs(s.medianGap) * 100).toFixed(1)} pts`,
      tone: s.medianGap !== null && s.medianGap >= 0.03 ? "coral" : "ink",
      sub: s.nLineage === 0 ? "No lineage split yet" : `from ${pct(s.medianRandom)} random · ${int(s.nLineage)} of ${int(m.datasets.length)} datasets`,
      href: "#lineage",
      onJump: jump("lineage"),
    },
    {
      key: "screen",
      label: "Error screen",
      value: b.n === 0 ? "–" : `${int(b.meets)} of ${int(b.n)}`,
      tone: b.n > 0 && b.meets === 0 ? "coral" : "live",
      sub: b.n === 0 ? "No error rates yet" : `${int(b.misses)} miss ≥1 limit · ${int(b.open)} undecided`,
      href: "#bars",
      onJump: jump("bars"),
    },
    {
      key: "study",
      label: "Study splits",
      value: `${int(st.logLossWorse)} of 5 worse`,
      tone: "coral",
      sub: `${int(st.datasets)} datasets in ${int(st.species)} species`,
      href: "#study-split",
      onJump: jump("study-split"),
    },
    {
      key: "mic",
      label: "MIC agreement",
      value: mic.n === 0 ? "–" : pct(micBase !== null && micBase.covered > 0 ? micBase.medianModel : mic.median),
      sub:
        mic.n === 0
          ? "No MIC models yet"
          : micBase !== null && micBase.covered > 0
            ? `vs ${pct(micBase.medianConstant)} constant · ${int(micBase.covered)} datasets`
            : `within 1 dilution · ${int(mic.n)} datasets`,
      href: "#mic",
      onJump: jump("mic"),
    },
  ];
};

/** Builds six concise 1-sentence takeaway cards, leaving full detail in each linked section. */
export const takeawayCards = (m: Manifest, constants?: ReadonlyMap<string, MicConstant>): TakeawayCard[] => {
  const rows = m.datasets;
  const s = splitSummary(rows);
  const b = barSummary(rows, "best", "lineage");
  const w = wordingSummary(rows, "best", "lineage");
  const ls = learnerSummary(rows);
  const scmRow = ls.rows.find((r) => r.learner === "scm");
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
      kicker: "Unseen splits",
      stat: s.medianGap === null ? "–" : `▼ ${(Math.abs(s.medianGap) * 100).toFixed(1)} pts median`,
      body:
        s.nLineage === 0 || s.medianGap === null
          ? `Holding whole studies out lowers log loss and specificity for ${int(st.logLossWorse)} of 5 learners across ${int(st.datasets)} datasets.`
          : `${int(s.nLargeGap)} of ${int(s.nLineage)} datasets drop over 5 points on unseen lineages, and holding whole studies (${int(st.datasets)}) or earlier dates (${int(extra.time_split.datasets)}) out worsens all ${int(st.logLossWorse)} learners; on ${int(outPrimary.frozen)} unseen-study genomes the rule set holds on ${int(outPrimary.holds.scm.holds)} of ${int(outPrimary.holds.scm.read)}.`,
      target: "study-split",
      linkLabel: "Study and date splits",
    },
    {
      id: "screen",
      kicker: "Error screen",
      stat: b.n === 0 ? "–" : `${int(b.meets)} of ${int(b.n)} meet both`,
      body:
        b.n === 0
          ? "No dataset has both error rates with lineages held out yet."
          : `Against VME ≤1.5% and ME ≤3.0% by exact 95% intervals, ${int(b.misses)} of ${int(b.n)} miss at least one limit and ${int(b.open)} are undecided; under the looser wording of FDA's 2009 guidance for growth-based susceptibility test systems, ${int(w.meets)} of ${int(w.n)} meet.`,
      target: "bars",
      linkLabel: "Error screen details",
    },
    {
      id: "learners",
      kicker: "Learners & inputs",
      stat: scmRow?.medianLineage == null ? "–" : `${pct(scmRow.medianLineage)} rule set`,
      body:
        scmRow?.medianLineage == null
          ? `Across ${int(inp.answered)} of ${int(inp.datasets)} datasets, no single input change beat the default.`
          : `Neither a dense neural network (${int(extra.dense_network.datasets)} datasets) nor any single input tweak across ${int(inp.species)} species beat the default; cutting to 1,000 genomes (${int(cap1000.worse)} of ${int(cap1000.datasets)} worse) or 1,000 features (${int(budget1000.worse)} of ${int(budget1000.datasets)} worse) hurt.`,
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
          : `Zero-training known-gene rules match or beat the best learner on ${int(cat.within.get(0.02) ?? 0)} of ${int(cat.nCompared)} datasets${highBand && highBand.compared > 0 ? ` (${int(highBand.ruleHolds)} of ${int(highBand.compared)} when the rule scores ≥85%)` : ""}; full AMRFinderPlus leads GRM's rule on ${int(extra.amrfinderplus.program_vs_known_gene_rule.ahead)} of ${int(extra.amrfinderplus.datasets)}, and TBProfiler leads the rule set on ${int(extra.tbprofiler.graded_vs_learners_lineage.scm.tool_higher)} of ${int(extra.tbprofiler.graded_endpoints)} TB drugs.`,
      target: "comparisons",
      linkLabel: "Outside tool comparisons",
    },
    {
      id: "other-db",
      kicker: "NCBI database",
      stat: `${int(od.totals.scm.plan_meets)} of ${int(od.eightSpecies)} species`,
      body: `On ${int(od.totalGenomes)} NCBI Pathogen Detection genomes, the rule set meets the plan rule in ${int(od.totals.scm.plan_meets)} of ${int(od.eightSpecies)} species, the decision tree in ${int(od.totals.cart.plan_meets)} of ${int(od.eightSpecies)} (E. coli pilot ${int(od.pilot.cart.holds)} of ${int(od.pilot.cart.estimable)}), and lasso in ${int(od.totals.lasso.plan_meets)} of ${int(od.totals.lasso.assessable)} species where it can be read.`,
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
            ? `Models agree within one doubling dilution on a median ${pct(micBase.medianModel)} of exact MICs (${int(micBase.covered)} datasets) against ${pct(micBase.medianConstant)} for a constant training-fold median MIC, not above the constant on ${int(micBase.notAbove)} of ${int(micBase.covered)}.`
            : `Models agree within one doubling dilution on a median ${pct(mic.median)} of exact MICs across ${int(mic.n)} datasets (${int(mic.nAbove90)} at 90% or more) on random folds.`,
      target: "mic",
      linkLabel: "MIC results",
    },
  ];
};

/** Expandable plain-English start-to-finish walkthrough of what GRM built, tested, and learned. */
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
        {open ? "Hide start-to-finish guide" : "Read start-to-finish guide"}
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
              The release covers {plural(manifest.datasets.length, "dataset")} in {plural(speciesCount(manifest), "species", "species")} across {int(analyzed)} selected genomes
              from {int(manifest.funnel.genomes_passing_qc)} QC-passing public BV-BRC genomes
              {capped ? ` (capped at ${int(cap)} per species)` : ""}. Beyond random splits, GRM tests unseen bacterial lineages ({int(s.nLineage)} datasets), held-out studies (
              {int(st.datasets)}), earlier-vs-later dates ({int(extra.time_split.datasets)}), newly deposited BV-BRC genomes, and {int(od.totalGenomes)} NCBI Pathogen Detection genomes.
            </p>
          </div>
          <div>
            <h3 className="text-control font-semibold text-ink">3. Findings</h3>
            <p className="mt-1.5 text-note text-ink">
              Every stricter split lowers scores: {int(b.meets)} of {int(b.n)} unseen-lineage datasets are shown to meet both error limits, and holding out whole studies or earlier
              dates worsens all 5 models. The rule set leads Kover 2.0.7 on {int(extra.kover.datasets)} small datasets and a dense neural net adds no gain ({int(extra.dense_network.datasets)}{" "}
              datasets); full AMRFinderPlus 4.2.7 leads GRM&apos;s simplified gene rule on {int(extra.amrfinderplus.program_vs_known_gene_rule.ahead)} of{" "}
              {int(extra.amrfinderplus.datasets)} datasets while roughly tying the 5 learners, and TBProfiler 6.6.3 leads the rule set on{" "}
              {int(extra.tbprofiler.graded_vs_learners_lineage.scm.tool_higher)} of {int(extra.tbprofiler.graded_endpoints)} graded TB drugs.
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
