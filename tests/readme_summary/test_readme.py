"""Checks that README.md's Revision 14 findings table matches the committed
reference files and that docs/INDEX.md lists every page under docs/.
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


def build_findings_table():
    """Rebuild the README's Revision 14 summary table from committed files."""
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
    assert study_ll_worse == 5 and study_sp_worse == 5 and study_rec_worse == {"rf", "gbt"}

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
    lasso_assessable = int(od_totals["lasso"]["species"]) - int(od_totals["lasso"]["plan_not_assessable"])

    inputs = load_json(ROOT / "web" / "src" / "data" / "input-comparison.json")
    arms = {a["id"]: a["all"] for a in inputs["arms"]}
    answered = sum(s["answered"] for s in inputs["species"])

    m = unc["medians"]
    g = unc["gaps"]
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
        f"({m['scm:lineage']['n']} datasets each) | [release-summary.md](docs/release-summary.md) |",
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
        f"{od_totals['cart']['species']} ({od_totals['cart']['holds']} of {od_totals['cart']['estimable']}), "
        f"lasso in {od_totals['lasso']['plan_meets']} of {lasso_assessable} assessable "
        f"({od_totals['lasso']['holds']} of {od_totals['lasso']['estimable']}) "
        f"| [other-database-test.md](docs/other-database-test.md) |",
        f"| Input settings ({len(inputs['species'])} species, H-1) | No single change shown better than default across "
        f"{answered} answered datasets (cap-4000 better {arms['cap-4000']['better']} of {arms['cap-4000']['datasets']}; "
        f"cap-1000 worse {arms['cap-1000']['worse']} of {arms['cap-1000']['datasets']}; "
        f"budget-1000 worse {arms['budget-1000']['worse']} of {arms['budget-1000']['datasets']}) "
        f"| [input-comparison-runner.md](docs/input-comparison-runner.md) |",
    ]


class ReadmeAndDocsIndex(unittest.TestCase):
    def test_findings_table_matches_committed_files(self):
        text = README.read_text(encoding="utf-8")
        expected = "\n".join(build_findings_table()) + "\n"
        self.assertIn(expected, text)

    def test_docs_index_covers_every_markdown_page_without_broken_links(self):
        index_text = INDEX.read_text(encoding="utf-8")
        actual = sorted(p.relative_to(DOCS).as_posix() for p in DOCS.rglob("*.md") if p != INDEX)
        linked = sorted(set(re.findall(r"\]\(((?!\.\./|https?://)[^)#]+\.md)\)", index_text)))
        self.assertEqual(linked, actual)
        for rel in linked:
            self.assertTrue((DOCS / rel).is_file(), f"missing target: docs/{rel}")

    def test_historical_snapshots_marked(self):
        index_text = INDEX.read_text(encoding="utf-8")
        self.assertIn("[work-program-20261027.md](work-program-20261027.md) — **Historical snapshot**", index_text)
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
