"""Checks that README.md's Revision 14 lead paragraph and findings table match
committed reference files, and that docs/INDEX.md lists every page under docs/.
"""

import csv
import json
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
DOCS = ROOT / "docs"
INDEX = DOCS / "INDEX.md"
REV14 = ROOT / "reference" / "release-summary" / "rev14-274-20261112"


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_lead_paragraph():
    """Rebuild the README's Revision 14 lead paragraph from committed files."""
    unc = load_json(REV14 / "uncertainty.json")
    study_first = load_tsv(ROOT / "reference" / "study-split" / "first.tsv")
    study_ll_worse = sum(
        1 for r in study_first if r["contrast"] == "study" and r["measure"] == "log_loss" and r["reading"] == "worse"
    )
    m = unc["medians"]
    eb = unc["error_bars"]["rows"]["best:lineage"]
    return (
        f"On public BV-BRC genomes, the best of five k-mer learners reaches a median balanced accuracy of "
        f"{m['best:random']['value'] * 100:.1f}% on random splits ({m['best:random']['n']} datasets) and "
        f"{m['best:lineage']['value'] * 100:.1f}% with whole lineages held out ({m['best:lineage']['n']} datasets). "
        f"With lineages held out, {eb['shown_to_meet_both']} of {eb['n']} datasets are shown to meet both limits "
        f"of the adapted research screen. Folds that hold whole studies out score worse than random folds for "
        f"{study_ll_worse} of 5 learners. Every figure is retrospective. No external or clinical validation has been done."
    )


def build_findings_table():
    """Rebuild the README's 14-row Revision 14 summary table from committed files."""
    summary = load_json(REV14 / "summary.json")
    unc = load_json(REV14 / "uncertainty.json")
    cat = load_json(REV14 / "catalog.json")
    guid = load_json(REV14 / "guidance.json")

    study_sum = {r["key"]: r["value"] for r in load_tsv(ROOT / "reference" / "study-split" / "summary.tsv")}
    study_first = load_tsv(ROOT / "reference" / "study-split" / "first.tsv")
    study_ll_worse = sum(
        1 for r in study_first if r["contrast"] == "study" and r["measure"] == "log_loss" and r["reading"] == "worse"
    )
    study_sp_worse = sum(
        1 for r in study_first if r["contrast"] == "study" and r["measure"] == "specificity" and r["reading"] == "worse"
    )
    study_rec_worse = {
        r["learner"] for r in study_first if r["contrast"] == "study" and r["measure"] == "recall" and r["reading"] == "worse"
    }
    if not (study_ll_worse == 5 and study_sp_worse == 5 and study_rec_worse == {"rf", "gbt"}):
        raise AssertionError(
            f"unexpected study-split readings: ll={study_ll_worse}, sp={study_sp_worse}, rec={study_rec_worse}"
        )

    outside = load_json(ROOT / "web" / "src" / "data" / "outside-study.json")
    primary = next(p for p in outside["populations"] if p["id"] == "primary")
    primary_frozen = sum(s["frozen"] for s in primary["species"])
    read_states = {"below", "consistent", "above"}

    def outside_holds(learner):
        rows = [
            r for r in primary["rows"]
            if r["learner"] == learner and r["sensitivity"] in read_states and r["specificity"] in read_states
        ]
        holds = sum(1 for r in rows if r["sensitivity"] != "below" and r["specificity"] != "below")
        return holds, len(rows)

    scm_h, scm_n = outside_holds("scm")
    cart_h, cart_n = outside_holds("cart")
    lasso_h, lasso_n = outside_holds("lasso")

    od_species = load_tsv(ROOT / "reference" / "other-database" / "species.tsv")
    od_genomes = sum(int(r["scored"]) for r in od_species)
    od_totals = {r["learner"]: r for r in load_tsv(ROOT / "reference" / "other-database" / "totals.tsv")}
    od_learners = load_tsv(ROOT / "reference" / "other-database" / "learners.tsv")
    ecoli_scm = next(r for r in od_learners if r["species"] == "escherichia_coli" and r["learner"] == "scm")
    ecoli_cart = next(r for r in od_learners if r["species"] == "escherichia_coli" and r["learner"] == "cart")
    lasso_assessable = int(od_totals["lasso"]["species"]) - int(od_totals["lasso"]["plan_not_assessable"])

    inputs = load_json(ROOT / "web" / "src" / "data" / "input-comparison.json")
    arms = {a["id"]: a["all"] for a in inputs["arms"]}
    answered = sum(s["answered"] for s in inputs["species"])

    time_sum = {r["name"]: r["value"] for r in load_tsv(ROOT / "reference" / "time-split" / "summary.tsv")}

    kover_cohort = {r["measure"]: r["chosen"] for r in load_tsv(ROOT / "reference" / "kover-comparison" / "cohort.tsv")}
    kover_ba = next(
        r for r in load_tsv(ROOT / "reference" / "kover-comparison" / "summary.tsv") if r["measure"] == "balanced_accuracy"
    )

    dense_totals = {r["key"]: r["value"] for r in load_tsv(ROOT / "reference" / "dense-network" / "totals.tsv")}
    dense_outcome = next(
        r["value"] for r in load_tsv(ROOT / "reference" / "dense-network" / "rule.tsv") if r["statistic"] == "outcome"
    )
    dense_best5 = next(
        r
        for r in load_tsv(ROOT / "reference" / "dense-network" / "summary.tsv")
        if r["split"] == "lineage" and r["comparator"] == "best5" and r["measure"] == "balanced_accuracy"
    )

    amrf_sum = load_tsv(ROOT / "reference" / "amrfinderplus-comparison" / "summary.tsv")
    amrf_rule = next(
        r
        for r in amrf_sum
        if r["set"] == "primary-new-235" and r["measure"] == "ba" and r["contrast"] == "program_minus_known_gene_rule"
    )
    amrf_scm = next(
        r
        for r in amrf_sum
        if r["set"] == "primary-new-235" and r["measure"] == "ba" and r["contrast"] == "learner_minus_program:lineage:scm"
    )
    amrf_figs = {r["key"]: r["value"] for r in load_tsv(ROOT / "reference" / "amrfinderplus-comparison" / "figures.tsv")}

    tbp_sum = {r["name"]: r["value"] for r in load_tsv(ROOT / "reference" / "tbprofiler-comparison" / "summary.tsv")}
    tbp_scm = next(
        r
        for r in load_tsv(ROOT / "reference" / "tbprofiler-comparison" / "comparisons.tsv")
        if r["way"] == "who_assoc" and r["set"] == "graded" and r["caller"] == "scm:lineage"
    )

    m = unc["medians"]
    g = unc["gaps"]
    pl = unc["paired_lineage"]
    eb = unc["error_bars"]["rows"]["best:lineage"]
    gb = guid["rows"]["best:lineage"]
    c_all = cat["bands"]["all"]
    c_hi = cat["bands"]["0.85 and up"]
    learner_ahead_over_2 = c_all["with_a_lineage_score"] - c_all["rule_ahead_or_within_2pt_lineage"]

    return [
        "| Question | Revision 14 finding | Details |",
        "| --- | --- | --- |",
        f"| Coverage | {summary['datasets']} datasets across {summary['species']} species have a report "
        f"({summary['best']['n_both']} with a lineage split) | [release-summary.md](docs/release-summary.md) |",
        f"| Random vs. lineage split (best learner) | Median balanced accuracy {m['best:random']['value']:.4f} random "
        f"({m['best:random']['n']}) vs. {m['best:lineage']['value']:.4f} lineage ({m['best:lineage']['n']}); "
        f"paired median drop {g['best:gap']['value']:.4f}, {summary['best']['drop_over_5pt']} of "
        f"{summary['best']['n_both']} lose >5 points | [release-summary.md](docs/release-summary.md) |",
        f"| Five learners on unseen lineages | SCM {m['scm:lineage']['value']:.4f}, CART {m['cart:lineage']['value']:.4f}, "
        f"lasso {m['lasso:lineage']['value']:.4f}, GBT {m['gbt:lineage']['value']:.4f}, RF {m['rf:lineage']['value']:.4f} "
        f"({m['scm:lineage']['n']} datasets each); SCM's paired gap to CART ({pl['scm-cart:lineage']['value']:.4f}) "
        f"and to lasso ({pl['scm-lasso:lineage']['value']:.4f}) does not clear zero | [release-summary.md](docs/release-summary.md) |",
        f"| Adapted research error screen (lineage, best learner) | {eb['shown_to_meet_both']} of {eb['n']} shown to meet "
        f"VME ≤1.5% and ME ≤3% by 95% intervals ({eb['shown_to_miss_one']} shown to miss, {eb['not_decided']} not decided); "
        f"{gb['meets_wording']} of {gb['n']} meet FDA 2009 guidance wording | [release-summary.md](docs/release-summary.md) |",
        f"| Known-gene rule vs. best learner (lineage) | Rule ahead or within 2 points on "
        f"{c_all['rule_ahead_or_within_2pt_lineage']} of {c_all['with_a_lineage_score']} "
        f"({c_hi['rule_ahead_or_within_2pt_lineage']} of {c_hi['with_a_lineage_score']} where rule ≥0.85); "
        f"learner ahead by >2 points on {learner_ahead_over_2} of {c_all['with_a_lineage_score']} "
        f"| [release-summary.md](docs/release-summary.md) |",
        f"| Study split inside BV-BRC | Holding BioProjects whole scores worse than random folds for {study_ll_worse} of 5 "
        f"learners on log loss and specificity across {study_sum['datasets']} datasets in {study_sum['species']} species "
        f"(recall worse for RF and GBT) | [study-split.md](docs/study-split.md) |",
        f"| Newly deposited BV-BRC studies | On {primary_frozen:,} unseen-study genomes (rev {outside['evaluation_revision']} "
        f"reports), lineage figures hold on {scm_h} of {scm_n} datasets for SCM, {cart_h} of {cart_n} for CART, "
        f"and {lasso_h} of {lasso_n} for lasso | [outside-study-results.md](docs/outside-study-results.md) |",
        f"| Another public database (NCBI) | On {od_genomes:,} genomes ({len(od_species)} species), SCM meets the plan's "
        f"rule in {od_totals['scm']['plan_meets']} of {od_totals['scm']['species']} species "
        f"({od_totals['scm']['holds']} of {od_totals['scm']['estimable']} datasets; plus {ecoli_scm['holds']} of "
        f"{ecoli_scm['estimable']} in *E. coli* pilot), CART in {od_totals['cart']['plan_meets']} of "
        f"{od_totals['cart']['species']} ({od_totals['cart']['holds']} of {od_totals['cart']['estimable']}; plus "
        f"{ecoli_cart['holds']} of {ecoli_cart['estimable']} in *E. coli* pilot), "
        f"lasso in {od_totals['lasso']['plan_meets']} of the {lasso_assessable} species where it can be read "
        f"({od_totals['lasso']['holds']} of its {od_totals['lasso']['estimable']} datasets across the eight) "
        f"| [other-database-test.md](docs/other-database-test.md) |",
        f"| Input settings ({len(inputs['species'])} species, H-1) | No single change shown better than default across "
        f"{answered} answered datasets (cap-4000 better {arms['cap-4000']['better']} of {arms['cap-4000']['datasets']}; "
        f"cap-1000 worse {arms['cap-1000']['worse']} of {arms['cap-1000']['datasets']}; "
        f"budget-1000 worse {arms['budget-1000']['worse']} of {arms['budget-1000']['datasets']}) "
        f"| [input-comparison-runner.md](docs/input-comparison-runner.md) |",
        f"| First-public-day split | Training on earlier public dates scores worse than a size-matched random draw on "
        f"log loss and balanced accuracy for {time_sum['total_log_loss_worse']} of {time_sum['learners']} learners "
        f"across {time_sum['datasets']} datasets in {time_sum['species']} species "
        f"(recall worse for {time_sum['total_recall_worse']} of {time_sum['learners']}) "
        f"| [time-split.md](docs/time-split.md) |",
        f"| Kover 2.0.7 vs. GRM SCM | On the {kover_cohort['datasets']} smallest lineage-split datasets "
        f"({int(kover_cohort['genome_dataset_entries']):,} genome-dataset entries), Kover minus GRM SCM balanced "
        f"accuracy is {float(kover_ba['value']):.4f} (lower on {kover_ba['datasets_lower']} of {kover_cohort['datasets']}, "
        f"higher on {kover_ba['datasets_higher']}, equal on {kover_ba['datasets_equal']}) "
        f"| [kover-comparison.md](docs/kover-comparison.md) |",
        f"| Dense neural network | Across {dense_totals['datasets']} datasets in {dense_totals['species']} species "
        f"(one setting on 4,000 screened features), {dense_outcome}: median "
        f"{float(dense_best5['difference_median']):.4f} balanced accuracy vs. hindsight best of 5 learners "
        f"(higher on {dense_best5['positive']} of {dense_totals['datasets']}) | [dense-network.md](docs/dense-network.md) |",
        f"| AMRFinderPlus 4.2.7 program | Across {amrf_rule['datasets']} datasets in {amrf_rule['species']} species, "
        f"the full program leads GRM's known-gene rule on {amrf_rule['ahead']} of {amrf_rule['datasets']} "
        f"({amrf_figs['program_ahead_by_0.10_or_more']} by >=10 points, {amrf_figs['of_them_quinolone']} quinolones); "
        f"SCM minus program median balanced accuracy is {float(amrf_scm['median']):+.4f} "
        f"(95% species CI {float(amrf_scm['median_species_lo']):+.4f} to {float(amrf_scm['median_species_hi']):+.4f}) "
        f"| [amrfinderplus-comparison.md](docs/amrfinderplus-comparison.md) |",
        f"| TBProfiler 6.6.3 (*M. tuberculosis*) | On 2,000 assemblies across {tbp_sum['graded_endpoints']} graded TB drugs "
        f"({tbp_sum['by_name_endpoints']} named endpoints, {tbp_sum['by_name_rows_without_call']} no-calls), TBProfiler "
        f"scores higher than SCM on {tbp_scm['tool_higher']} of {tbp_scm['endpoints']} "
        f"(SCM minus TBProfiler median {float(tbp_scm['median_minus_tool']):.4f}) "
        f"| [tbprofiler-comparison.md](docs/tbprofiler-comparison.md) |",
    ]


class ReadmeAndDocsIndex(unittest.TestCase):
    def test_lead_paragraph_figures_match_committed_files(self):
        text = README.read_text(encoding="utf-8")
        self.assertIn(build_lead_paragraph() + "\n", text)

    def test_findings_table_matches_committed_files(self):
        text = README.read_text(encoding="utf-8")
        expected = "\n".join(build_findings_table()) + "\n"
        self.assertIn(expected, text)
        self.assertIn(
            "They are an adapted research screen. They are not an applicable FDA criterion, "
            "and meeting them does not make a learner fit for clinical use.",
            text,
        )

    def test_docs_index_covers_every_markdown_page_without_broken_links_or_em_dashes(self):
        index_text = INDEX.read_text(encoding="utf-8")
        self.assertNotIn("\u2014", index_text)
        actual = sorted(p.relative_to(DOCS).as_posix() for p in DOCS.rglob("*.md") if p != INDEX)
        linked = sorted(set(re.findall(r"\]\(((?!\.\./|https?://)[^)#]+\.md)\)", index_text)))
        self.assertEqual(linked, actual)
        for rel in linked:
            self.assertTrue((DOCS / rel).is_file(), f"missing target: docs/{rel}")

    def test_historical_snapshots_marked(self):
        index_text = INDEX.read_text(encoding="utf-8")
        self.assertIn("[research/genotype-to-phenotype.md](research/genotype-to-phenotype.md) | **Historical snapshot**", index_text)
        g2p = (DOCS / "research" / "genotype-to-phenotype.md").read_text(encoding="utf-8")
        self.assertIn("> **Historical snapshot (October 2026, evaluation revision 7).**", g2p)

    def test_ci_and_docs_wire_this_suite(self):
        ci_yml = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
        self.assertIn("          - suite: tests/readme_summary\n            top: .\n", ci_yml)
        n = len(unittest.defaultTestLoader.getTestCaseNames(ReadmeAndDocsIndex))
        ci_md = (DOCS / "ci.md").read_text(encoding="utf-8")
        self.assertIn(f"`tests/readme_summary` {n}", ci_md)


if __name__ == "__main__":
    unittest.main()
