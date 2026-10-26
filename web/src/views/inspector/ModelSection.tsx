import { useState } from "react";
import type { CatalogRow, DatasetReport, FeatureCard, Rule } from "../../data/schema";
import { int, pct, plural } from "../../data/format";
import { featureIndex, featureLabel, featureWhere, ruleWords } from "../../data/derive/dataset";
import { TreeDiagram } from "../../charts/TreeDiagram";
import { Section, Segmented } from "../../ui/primitives";
import { Seq } from "../../ui/Seq";

const RuleBlock = ({ index, rule, card, r }: { index: number; rule: Rule; card: FeatureCard | undefined; r: DatasetReport }) => {
  const features = featureIndex(r);
  const step = r.models.scm.steps[index];
  const importance = r.models.scm.importance[index];
  const place = featureWhere(card, rule.feature);
  const locus = card?.locus ?? null;
  return (
    <article className="border-t border-hair py-4">
      <h4 className="flex flex-wrap items-baseline gap-x-3 gap-y-1 text-body font-semibold text-ink">
        <span className="text-control font-medium tabular-nums text-muted-text">Rule {index + 1}</span>
        {ruleWords(rule, features)}
        {card !== undefined && card.catalog_families.length > 0 ? (
          <span className="text-note font-normal text-muted-text">
            Catalog {card.catalog_families.length === 1 ? "family" : "families"} {card.catalog_families.join(", ")}
          </span>
        ) : null}
      </h4>
      <div className="mt-2 grid grid-cols-1 gap-x-8 gap-y-1 text-note text-ink md:grid-cols-2">
        <p>
          {place ?? "Not placed on a reference assembly."}
          {locus !== null && locus.gene !== null && locus.product !== null && locus.product !== locus.gene ? <span className="text-muted-text"> · {locus.product}</span> : null}
        </p>
        {card !== undefined ? (
          <p>
            Carried by {int(card.n_resistant)} of {int(r.counts.resistant)} resistant and {int(card.n_susceptible)} of {int(r.counts.susceptible)} susceptible isolates.
          </p>
        ) : null}
        {importance !== undefined ? <p>Dropping this rule raises training balanced error by {(importance * 100).toFixed(1)} points.</p> : null}
        {step !== undefined && step.n_equivalent > 0 ? (
          <p>
            {plural(step.n_equivalent, "other feature")} tie this rule exactly
            {step.equivalent.length > 0 ? `, including ${step.equivalent.slice(0, 3).map((e) => featureLabel(features.get(e.feature), e.feature)).join(", ")}` : ""}.
          </p>
        ) : null}
      </div>
      {locus !== null ? (
        <div className="mt-3">
          <p className="text-tiny text-muted-text">
            {plural(card?.n_kmers ?? 0, "k-mer")}; {int(locus.kmers_placed)} placed on genome {locus.genome_id}, contig {locus.contig}, {int(locus.start)}–{int(locus.end)}
            {locus.gene_strand !== null ? `, gene on the ${locus.gene_strand} strand` : ""}. Span in capitals.
          </p>
          <Seq context={locus.context} className="mt-1" />
        </div>
      ) : card !== undefined ? (
        <p className="mt-2 text-tiny text-muted-text">
          {int(card.n_kmers)} k-mers, for example {card.kmer}.
        </p>
      ) : null}
    </article>
  );
};

export const ModelSection = ({ r, row }: { r: DatasetReport; row?: CatalogRow }) => {
  const features = featureIndex(r);
  const scm = r.models.scm;
  const [deep, setDeep] = useState(false);
  const shallow = Math.min(3, r.models.cart_depth);
  const root = r.models.cart.nodes[0];
  const split = root?.split ?? null;
  const present = split === null ? undefined : r.models.cart.nodes[split.present];
  const absent = split === null ? undefined : r.models.cart.nodes[split.absent];

  return (
    <Section
      id="ds-model"
      title="The interpretable model"
      answer={
        scm.rules.length === 0
          ? "The set covering machine found no rule that beats calling every isolate one way."
          : `The set covering machine uses ${plural(scm.rules.length, "rule")}. It calls resistant when ${scm.kind === "disjunction" ? "any" : "every"} rule holds.`
      }
    >
      <div>
        {scm.rules.map((rule, i) => (
          <RuleBlock key={`${rule.feature}-${i}`} index={i} rule={rule} card={features.get(rule.feature)} r={r} />
        ))}
        <p className="type-note border-t border-hair pt-3">
          Rules were chosen greedily on all isolates for the final model, with the trade-off p = {scm.p} and the rule count picked by inner cross-validation (error {pct(scm.cv_error)}).
        </p>
      </div>

      <div className="mt-8">
        <div className="flex flex-wrap items-start justify-between gap-3">
          <div>
            <h3 className="text-control font-semibold text-ink">Decision tree (CART), depth {r.models.cart_depth}</h3>
            {split !== null && present !== undefined && absent !== undefined ? (
              <p className="mt-1 text-note text-ink">
                The first split is on {featureLabel(features.get(split.feature), split.feature)}: {pct(present.n === 0 ? 0 : present.n_resistant / present.n)} of the {int(present.n)} isolates carrying it are resistant, against {pct(absent.n === 0 ? 0 : absent.n_resistant / absent.n)} of the {int(absent.n)} without.
              </p>
            ) : (
              <p className="mt-1 text-note text-ink">The tree has no split: one leaf calls every isolate.</p>
            )}
          </div>
          {r.models.cart_depth > shallow ? (
            <Segmented
              ariaLabel="Tree depth shown"
              options={[
                { id: "shallow", label: `First ${shallow} levels` },
                { id: "all", label: `All ${r.models.cart_depth} levels` },
              ]}
              value={deep ? "all" : "shallow"}
              onSelect={(v) => setDeep(v === "all")}
            />
          ) : null}
        </div>
        <div className="mt-3">
          <TreeDiagram nodes={r.models.cart.nodes} depth={deep ? r.models.cart_depth : shallow} features={features} />
        </div>
      </div>
      {row?.kover == null ? null : (
        <p className="type-note mt-6">Legacy Kover: {pct(row.kover.balanced_accuracy)} balanced accuracy, {pct(row.kover.sensitivity)} sensitivity and {pct(row.kover.specificity)} specificity on {int(row.kover.n_examples)} examples, with {plural(row.kover.n_rules, "rule")}.</p>
      )}
    </Section>
  );
};
