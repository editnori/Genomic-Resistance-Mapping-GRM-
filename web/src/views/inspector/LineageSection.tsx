import { useState } from "react";
import type { DatasetReport } from "../../data/schema";
import { distance, int, LEARNER_NAME, pct, plural, stLabel } from "../../data/format";
import { loadErrorText, useSpecies } from "../../data/load";
import { bestLearner, oofErrors } from "../../data/derive/dataset";
import { ScatterMap, type MapPoint } from "../../charts/ScatterMap";
import { Dot, LegendItem, PanelState, Section } from "../../ui/primitives";

/** Share of resistant isolates held by the three lineages with the most. */
const topThreeShare = (r: DatasetReport): { lineages: number; share: number } => {
  const counts = new Map<number, number>();
  r.samples.forEach((s) => {
    if (s.resistant) counts.set(s.lineage, (counts.get(s.lineage) ?? 0) + 1);
  });
  const sorted = [...counts.values()].sort((a, b) => b - a);
  const top = sorted.slice(0, 3).reduce((a, b) => a + b, 0);
  return { lineages: counts.size, share: r.counts.resistant === 0 ? 0 : top / r.counts.resistant };
};

export const LineageSection = ({ r }: { r: DatasetReport }) => {
  const species = useSpecies(r.species_slug);
  const [showErrors, setShowErrors] = useState(false);
  const best = bestLearner(r);
  const errors = best === undefined ? [] : oofErrors(r, best);
  const nErrors = errors.filter(Boolean).length;
  const spread = topThreeShare(r);

  const sampleIndex = new Map(r.samples.map((s, i) => [s.id, i]));
  const points: MapPoint[] =
    species.status !== "ready"
      ? []
      : species.data.genomes.map((g) => {
          const i = sampleIndex.get(g.id);
          const s = i === undefined ? undefined : r.samples[i];
          return {
            x: g.x,
            y: g.y,
            layer: s === undefined ? "context" : "mark",
            tone: s === undefined ? "neutral" : s.resistant ? "r" : "s",
            flag: showErrors && i !== undefined && errors[i] === true,
            tip: () => [
              { term: g.id, text: s === undefined ? `Not tested for ${r.antibiotic}.` : `${s.resistant ? "Resistant" : "Susceptible"}${best !== undefined && i !== undefined ? `; ${LEARNER_NAME[best.learner]} scored ${pct(best.scores[i] ?? 0, 0)} resistant out of fold` : ""}.` },
              { source: `Lineage ${g.lineage}${g.mlst !== null ? ` · ${stLabel(g.mlst)}` : ""}${g.country !== null ? ` · ${g.country}` : ""}${g.year !== null ? ` ${g.year}` : ""}` },
            ],
          };
        });

  return (
    <Section
      id="ds-lineages"
      title="Lineages"
      answer={`Resistant isolates fall in ${plural(spread.lineages, "lineage")} of ${int(r.counts.lineages)}; the three with the most hold ${pct(spread.share)} of them.`}
      right={
        best !== undefined ? (
          <label className="flex min-h-11 cursor-pointer items-center gap-2 text-control text-ink">
            <input type="checkbox" checked={showErrors} onChange={(e) => setShowErrors(e.target.checked)} className="size-4 accent-[rgb(var(--blue))]" />
            {LEARNER_NAME[best.learner]} errors ({int(nErrors)})
          </label>
        ) : undefined
      }
    >
      <div className="mb-3 flex flex-wrap gap-x-5 gap-y-1">
        <LegendItem swatch={<Dot color="var(--series-r)" />} name="Resistant" value={int(r.counts.resistant)} />
        <LegendItem swatch={<Dot color="var(--series-s)" />} name="Susceptible" value={int(r.counts.susceptible)} />
        <LegendItem swatch={<Dot color="var(--series-neutral)" size={6} />} name="Not in this dataset" />
        {showErrors ? <LegendItem swatch={<Dot color="rgb(var(--ink))" ring size={12} />} name="Out-of-fold error, random split" value={int(nErrors)} /> : null}
      </div>
      {species.status === "loading" ? <PanelState kind="loading" title="Loading the species map" /> : null}
      {species.status === "error" ? <PanelState kind="error" title="The species map could not be read" message={loadErrorText(species.error)} /> : null}
      {species.status === "ready" ? (
        <>
          <ScatterMap points={points} label={`Map of ${int(species.data.genomes.length)} ${r.species} genomes by Mash distance, this dataset's isolates colored by phenotype.`} />
          <p className="type-note mt-2">
            Principal coordinates of Mash distances between genomes; two coordinates only approximate those distances. Lineages are clusters at Mash distance {distance(species.data.lineage_threshold)}, computed from the distances directly.
          </p>
        </>
      ) : null}
    </Section>
  );
};
