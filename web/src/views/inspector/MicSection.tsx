import type { DatasetReport, MicDataset, MicRun, Research } from "../../data/schema";
import type { Resource } from "../../data/load";
import { researchFor, type GuardPart, type MicResearchState } from "../../data/research";
import { int, mgL, pct, plural } from "../../data/format";
import { BARS } from "../../data/metrics";
import { MicHeatmap } from "../../charts/MicHeatmap";
import { verdictWords, PanelState, ScrollFrame, Section } from "../../ui/primitives";

export const MicSection = ({ r, research }: { r: DatasetReport; research: Resource<Research> }) => {
  const mic = r.mic;
  const withMic = r.samples.filter((s) => s.mic !== null);
  const exact = withMic.filter((s) => s.mic?.censor === "Exact").length;
  return (
    <Section
      id="ds-mic"
      title="MIC prediction"
      answer={
        mic === null
          ? undefined
          : `Predicted MIC is within one doubling dilution of the measured MIC for ${pct(mic.essential_agreement)} of the ${int(mic.n_exact)} isolates with an exact MIC.`
      }
    >
      {mic === null ? (
        <PanelState
          kind="empty"
          title="No MIC regression for this dataset"
          message={`MIC regression needs at least 100 isolates with an MIC, 50 of them exact. This dataset has ${int(withMic.length)} with an MIC, ${int(exact)} exact.`}
        />
      ) : (
        // grid-cols-1 keeps the narrow layout's one column at the dialog's width:
        // the heatmap sizes itself from its frame, and an auto column would grow
        // to the heatmap's own width and keep it there.
        <div className="grid grid-cols-1 gap-x-10 gap-y-6 lg:grid-cols-[minmax(0,auto)_minmax(0,1fr)]">
          <MicHeatmap confusion={mic.confusion} />
          <div className="space-y-3 text-body text-ink">
            <p>
              Counting censored values too, {pct(mic.essential_agreement_all)} of {int(mic.n)} isolates agree: a value reported as at most X agrees when the prediction is at most one step above X, and a value above X agrees when the prediction is at least X.
            </p>
            {mic.categorical !== null && mic.fold_breakpoints !== undefined ? (
              <p>
                Calling resistance from the predicted MIC at each training fold&apos;s own threshold (resistant above{" "}
                {[...new Set(mic.fold_breakpoints.flatMap((b) => (b === null ? [] : [b.r_above])))]
                  .sort((a, b) => a - b)
                  .map(mgL)
                  .join(" or ")}{" "}
                mg/L, fitted on that fold&apos;s training isolates) and scoring the calls against the reported resistant and susceptible labels gives{" "}
                {verdictWords(mic.categorical.vme, BARS.vme)} and {verdictWords(mic.categorical.me, BARS.me)} over {int(mic.categorical.n)} isolates with an MIC.
              </p>
            ) : mic.categorical !== null && mic.breakpoint !== null ? (
              <p>
                This report predates training-fold thresholds: its categorical errors{" "}
                {verdictWords(mic.categorical.vme, BARS.vme)} and {verdictWords(mic.categorical.me, BARS.me)} used the cohort breakpoint (resistant above{" "}
                {mgL(mic.breakpoint.r_above)} mg/L), which was fitted on all isolates including the ones scored.
              </p>
            ) : mic.fold_breakpoints !== undefined ? (
              <p>No categorical calls: at least one training fold had too few isolates to fit a threshold, and the pipeline does not score a subset.</p>
            ) : null}
            <p className="type-note">
              Censored gradient-boosted regression on log2 dilution steps, predicted out of fold on the random split.
            </p>
          </div>
        </div>
      )}
      <MicResearch r={r} research={research} />
    </Section>
  );
};

/** GRM fits an MIC model only with at least this many exact MICs; a rerun
 *  scored on fewer carries the exact-minimum limit. */
const EXACT_MINIMUM = 50;

/** What each guard part covers, as a mismatch names it. */
const PART_WORDS: Record<GuardPart, string> = {
  samples: "isolates",
  sample_ids: "isolates",
  labels: "resistant and susceptible labels",
  lineages: "lineage clusters",
  mic_values: "measured MICs",
  mic_predicted: "saved MIC predictions",
  scores: "held-out scores",
};

/** "a and b", "a, b, and c". */
const andList = (items: readonly string[]): string =>
  items.length < 3 ? items.join(" and ") : `${items.slice(0, -1).join(", ")}, and ${items[items.length - 1]}`;

/** Why the research for this dataset id gives this report none of its
 *  numbers. When only one side has an MIC model, the predictions part stands
 *  for that model, so the sentence names the model and no predictions. */
const unusedText = (state: Extract<MicResearchState, { status: "mismatch" | "no_model" }>): string => {
  const model =
    state.status === "no_model"
      ? "it scored GRM's MIC model and this report has none"
      : state.researchModel
        ? null
        : "this report has an MIC model and the revision-5 report it read had none";
  const parts = state.failed.filter((p) => model === null || p !== "mic_predicted");
  const reasons = [
    ...(model === null ? [] : [model]),
    ...(parts.length === 0 ? [] : [`this report's ${andList([...new Set(parts.map((p) => PART_WORDS[p]))])} differ from the ones it read`]),
  ];
  return `Revision-5 MIC research exists for this dataset, but ${reasons.join(", and ")}, so none of its numbers apply to this report.`;
};

/** GRM's MIC model beside the no-genome baseline from research.json. Numbers
 *  show only for an entry whose guard matches this report, which has an MIC
 *  model; otherwise the block says why none apply, and a withheld model shows
 *  only its reason. Nothing shows while research.json loads or when it is
 *  missing or invalid. */
const MicResearch = ({ r, research }: { r: DatasetReport; research: Resource<Research> }) => {
  if (research.status !== "ready") return null;
  const state = researchFor(research, r).mic;
  if (state.status === "loading" || state.status === "absent") return null;
  return (
    <div className="mt-8">
      <h3 className="text-control font-semibold text-ink">No-genome baseline</h3>
      {state.status === "current" ? (
        <BaselineComparison entry={state.entry} limits={research.data.mic.limits} />
      ) : (
        <p className="prose-measure mt-1 text-note text-ink">{state.status === "withheld" ? state.reason : unusedText(state)}</p>
      )}
    </div>
  );
};

const SPLITS = [
  { scheme: "random", label: "Random split" },
  { scheme: "lineage", label: "Lineages held out" },
] as const;

interface SplitRow {
  scheme: (typeof SPLITS)[number]["scheme"];
  label: string;
  model: MicRun;
  baseline: MicRun | null;
}

const agreement = (run: MicRun): string => (run.ea_exact === null ? "No exact MIC scored" : `${pct(run.ea_exact)} of ${int(run.n_exact)}`);

/** The splits where the no-genome baseline has higher essential agreement
 *  than GRM's model, as one sentence; null on a dataset where it never does.
 *  parseResearch makes the two score the same exact MICs on each split. */
const baselineWins = (rows: readonly SplitRow[]): string | null => {
  const won = rows.flatMap((x) =>
    x.baseline !== null && x.baseline.ea_exact !== null && x.model.ea_exact !== null && x.baseline.ea_exact > x.model.ea_exact
      ? [{ scheme: x.scheme, n: x.model.n_exact, b: x.baseline.ea_exact, m: x.model.ea_exact }]
      : [],
  );
  const first = won[0];
  if (first === undefined) return null;
  const where = (scheme: SplitRow["scheme"]) => (scheme === "random" ? "on the random split" : "with lineages held out");
  // Shares one exact MIC apart can round alike once a split has over 1,000.
  const pair = (b: number, m: number) => {
    const digits = pct(b) === pct(m) ? 2 : 1;
    return `${pct(b, digits)} against ${pct(m, digits)}`;
  };
  const lead = "The no-genome baseline has higher essential agreement than GRM's model";
  if (won.length === 1) return `${lead} ${where(first.scheme)}: ${pair(first.b, first.m)} of the same ${int(first.n)} exact MICs.`;
  const shared = won.every((x) => x.n === first.n);
  const clauses = won.map((x) => `${pair(x.b, x.m)}${shared ? "" : ` of ${int(x.n)} exact MICs`} ${where(x.scheme)}`);
  return `${lead} on both splits: ${clauses.join(" and ")}${shared ? `, of the same ${int(first.n)} exact MICs` : ""}.`;
};

const BaselineComparison = ({ entry, limits }: { entry: MicDataset; limits: Research["mic"]["limits"] }) => {
  // A saved report model has the random split only and no baseline; its row
  // stays, with the baseline marked as not run.
  const rows: SplitRow[] = SPLITS.flatMap(({ scheme, label }) => {
    const model = entry.model[scheme];
    return model === null ? [] : [{ scheme, label, model, baseline: entry.baseline?.[scheme] ?? null }];
  });
  const removed = entry.zone_diameters.removed;
  const holdout = entry.model.lineage?.holdout_random_folds ?? 0;
  const wins = baselineWins(rows);
  const notes = [
    limits.scope,
    ...(entry.baseline !== null ? [limits.baseline] : []),
    ...(entry.model.source === "saved_report" ? [limits.saved_only] : []),
    ...(rows.some((x) => x.model.n_exact < EXACT_MINIMUM) ? [limits.exact_minimum] : []),
    ...(removed > 0 ? [limits.zone_fix.replace("{removed}", plural(removed, "disk-diffusion zone diameter"))] : []),
    ...(holdout > 0
      ? [`In ${plural(holdout, "lineage-split fold")}, the early-stopping holdout (training isolates set aside to decide when boosting stops) fell back to random isolates instead of whole lineages.`]
      : []),
  ];
  return (
    <>
      <p className="prose-measure mt-1 text-note text-ink">
        {entry.model.source === "rerun"
          ? "The revision-5 MIC research refit GRM's MIC model on these isolates and scored the no-genome baseline on the same folds."
          : "The revision-5 MIC research scored this report's saved MIC predictions."}{" "}
        Essential agreement is the share of exact MICs predicted within one doubling dilution.
      </p>
      <div className="mt-3">
        <ScrollFrame label="Essential agreement against a no-genome baseline">
          <table className="w-full max-w-xl border-collapse text-control">
            <thead>
              <tr className="border-b border-hair">
                {["Split", "GRM model", "No-genome baseline"].map((h, i) => (
                  <th key={h} scope="col" className={`px-2 py-2 text-tiny font-medium text-blue ${i === 0 ? "text-left" : "text-right"}`}>
                    {h}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => (
                <tr key={row.scheme} className="border-b border-hair">
                  <th scope="row" className="px-2 py-1.5 text-left font-medium text-ink">
                    {row.label}
                  </th>
                  <td className="whitespace-nowrap px-2 py-1.5 text-right tabular-nums">{agreement(row.model)}</td>
                  <td className="whitespace-nowrap px-2 py-1.5 text-right tabular-nums">{row.baseline === null ? <span className="text-muted-text">Not run</span> : agreement(row.baseline)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </ScrollFrame>
      </div>
      {wins === null ? null : <p className="prose-measure mt-3 text-note text-ink">{wins}</p>}
      <ul className="type-note mt-3 list-disc space-y-0.5 pl-5">
        {notes.map((n) => (
          <li key={n}>{n}</li>
        ))}
      </ul>
    </>
  );
};
