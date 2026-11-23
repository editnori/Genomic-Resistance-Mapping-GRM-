"""scripts/mechanism/page_tables.py: the gates and the tables of docs/mechanism-split-result.md.

Every species, dataset, antibiotic, class, family, count and figure here is
SYNTHETIC: three made-up species, six datasets of A, one of them in M, laid
out as run_species.sh's copies would hold them, and combine.py run over them.
"""

import contextlib
import hashlib
import io
import tempfile
import unittest
from pathlib import Path

from tests.mechanism import common as c

page = c.load("page_tables")
combine = c.load("combine")
strata = c.load("strata")
history = c.load("history")

COMMIT = "c" * 40
FIRST_RUN = "f" * 40
SPECIES = {"synthia_alba": "Synthia alba", "synthia_borea": "Synthia borea", "synthia_cana": "Synthia cana"}
DRUGS = ("synthomycin", "fictocillin")
IN_M = {"synthia_cana__fictocillin"}
NOT_IN_A = {"synthia_alba__nullomycin"}
LISTS_HEADER = ["dataset", "species_slug", "species", "antibiotic", "in_M", "in_A", "labels_depend_on"]
CONTRASTS = ["dataset", "species", "learner", "contrast", "d_recall", "d_log_loss", "d_specificity",
             "d_balanced_accuracy"]
# resistant and susceptible genomes per stratum of a dataset; all made up
RES = {"point_only": 11, "gene_only": 12, "both": 3, "neither": 10}
SUS = {"point_only": 2, "gene_only": 1, "both": 0, "neither": 40}
GUARD_P = {"bioproject": "0.500000000000", "pooled": "", "largest_out": "0.001000000000",
           "flagged_out": "0.300000000000", "lineage": "0.200000000000", "lineage_both": "0.300000000000",
           "countries_both": "0.300000000000", "intact": "0.300000000000"}


def lists_rows():
    out = []
    for slug, name in SPECIES.items():
        for drug in DRUGS:
            d = f"{slug}__{drug}"
            out.append({"dataset": d, "species_slug": slug, "species": name, "antibiotic": drug,
                        "in_M": int(d in IN_M), "in_A": 1, "labels_depend_on": "the maintainer's rule" if d in IN_M
                        else ""})
    out.append({"dataset": "synthia_alba__nullomycin", "species_slug": "synthia_alba", "species": "Synthia alba",
                "antibiotic": "nullomycin", "in_M": 0, "in_A": 0, "labels_depend_on": ""})
    return out


def recall(called, res):
    return f"{called / res:.12f}" if res >= 10 else ""


def species_tables(slug, k):
    """The ten result tables of one SYNTHETIC species run, as rows."""
    name = SPECIES[slug]
    datasets, strata_rows, pooled, nocall = [], [], [], []
    for drug in DRUGS:
        d = f"{slug}__{drug}"
        dset = "M" if d in IN_M else "primary"
        base = {"dataset": d, "species_slug": slug, "species": name, "antibiotic": drug, "drug_class": "synthoids",
                "set": dset, "labels_depend_on": "the maintainer's rule" if dset == "M" else ""}
        datasets.append(dict(base, genomes=100, with_table=99, rule_check=1, hard_check=1, learner_rows_equal=12,
                             learner_rows=12, program_check=1, record_check=1, entered=1, why=""))
        res = dict(RES)
        if d == "synthia_borea__fictocillin":
            res["point_only"] = 9
        for split in ("lineage", "random"):
            for i, row in enumerate(strata.ROWS):
                # gene only called a little more often than point only, more so for some rows
                called = {"gene_only": res["gene_only"] - (i % 3), "point_only": res["point_only"] - (i % 3) - k,
                          "both": res["both"], "neither": res["neither"] // 2}
                for s in strata.STRATA:
                    strata_rows.append(dict(base, split=split, row=row, stratum=s, kind="stratum",
                                            resistant=res[s], resistant_called=max(called[s], 0),
                                            susceptible=SUS[s], susceptible_called=0,
                                            recall=recall(max(called[s], 0), res[s])))
                for b in strata.BESIDE:
                    n = sum(res.values()) if b == "all" else 1
                    strata_rows.append(dict(base, split=split, row=row, stratum=b, kind="beside", resistant=n,
                                            resistant_called=n // 2, susceptible=0, susceptible_called=0,
                                            recall=recall(n // 2, n)))
                if row in strata.LEARNERS:
                    nocall.append({"dataset": d, "split": split, "learner": row, "no_call_resistant": 0,
                                   "no_call_susceptible": 0})
        for split in ("lineage", "random"):
            for row in strata.ROWS:
                pooled.append({"species_slug": slug, "species": name, "drug_class": "synthoids", "set": dset,
                               "split": split, "row": row, "stratum": "both", "datasets": 1, "resistant_pairs": 3,
                               "resistant_called": 3, "susceptible_pairs": 0, "susceptible_called": 0, "recall": ""})
    cells, families = [], []
    for cls, passes in (("QUINOLONE", True), ("TETRACYCLINE", False)):
        for periods in ("primary", "secondary"):
            base = {"species_slug": slug, "species": name, "cell": "species_class", "class": cls, "periods": periods,
                    "genomes": 99, "with_year_share": "0.900000000000", "earlier_genomes": 40, "later_genomes": 50,
                    "earlier_units": 4 if passes else 2, "later_units": 5 if passes else 1,
                    "earlier_point_disrupt": 0, "later_point_disrupt": 1}
            for guard in history.GUARDS:
                row = dict(base, guard=guard, passes="" if guard == "pooled" else int(passes))
                if passes:
                    row.update(earlier_point="0.200000000000", later_point=f"{0.3 + k / 100:.12f}",
                               d_point=f"{0.1 + k / 100:.12f}", earlier_gene="0.500000000000",
                               later_gene="0.600000000000", d_gene="0.100000000000", tvd="0.050000000000",
                               p_point=GUARD_P[guard], p_tvd=GUARD_P[guard])
                else:
                    row["note"] = "under 3 voting BioProjects in a period, no share written"
                cells.append(row)
            if passes:
                families.append({"species_slug": slug, "species": name, "class": cls, "periods": periods,
                                 "family": "qnrS", "earlier_genomes": 10, "later_genomes": 12,
                                 "earlier_mix": "0.800000000000", "later_mix": "0.700000000000"})
                families.append({"species_slug": slug, "species": name, "class": cls, "periods": periods,
                                 "family": "other families", "earlier_genomes": 3, "later_genomes": 2,
                                 "earlier_mix": "0.200000000000", "later_mix": "0.300000000000"})
    for drug in DRUGS:
        d = f"{slug}__{drug}"
        if d in IN_M:
            continue
        for periods in ("primary", "secondary"):
            for guard in ("bioproject", "pooled"):
                cells.append({"species_slug": slug, "species": name, "cell": "dataset_resistant", "class": "",
                              "dataset": d, "antibiotic": drug, "periods": periods, "guard": guard, "genomes": 36,
                              "earlier_genomes": 10, "later_genomes": 12, "earlier_units": 1, "later_units": 2,
                              "passes": "" if guard == "pooled" else 0, "earlier_point_disrupt": 0,
                              "later_point_disrupt": 0, "note": "under 3 voting BioProjects in a period"})
    years = [{"species_slug": slug, "species": name, "years": b, "genomes": n}
             for b, n in (("no year", 9), ("before 2000", 0), ("2000 to 2009", 10), ("2010 to 2014", 40),
                          ("2015 to 2019", 50), ("2020 on", 0), ("no table", 1))]
    lineage = [{"species_slug": slug, "species": name, "reports": 2, "report_pairs": 1, "differing_genome_pairs": 0,
                "given": 1}]
    d = f"{slug}__synthomycin"
    ts = [{"dataset": d, "species_slug": slug, "species": name, "antibiotic": "synthomycin", "reading": reading,
           "resistant_later_test": 8, "with_determinant": 6, "with_unseen_determinant": 2, "no_determinant": 2,
           "share_unseen": f"{2 / 6:.12f}"} for reading in ("point_as_one", "point_by_symbol")]
    return {
        "m1/datasets.tsv": (strata.DATASETS_HEADER, datasets),
        "m1/strata.tsv": (strata.STRATA_HEADER, strata_rows),
        "m1/pooled.tsv": (strata.POOLED_HEADER, pooled),
        "m1/nocall.tsv": (strata.NOCALL_HEADER, nocall),
        "m2/m2-cells.tsv": (history.CELLS_HEADER, cells),
        "m2/m2-families.tsv": (history.FAMILIES_HEADER, families),
        "m2/m2-years.tsv": (history.YEARS_HEADER, years),
        "m2/m2-lineage.tsv": (history.LINEAGE_HEADER, lineage),
        "m2/timesplit.tsv": (history.TIMESPLIT_HEADER, ts),
    }


def sums(results):
    names = sorted(p.relative_to(results).as_posix() for p in results.rglob("*")
                   if p.is_file() and p.name != "SHA256SUMS")
    text = "".join(f"{hashlib.sha256((results / n).read_bytes()).hexdigest()}  ./{n}\n" for n in names)
    (results / "SHA256SUMS").write_text(text)
    return hashlib.sha256(text.encode()).hexdigest(), len(names)


def lay_run(run, slug, tables, commit=COMMIT):
    res = run / "results"
    for name, (header, rows) in tables.items():
        c.write_tsv(res / name, header, rows)
    (res / "m1" / "summary.json").write_text("{}\n")
    digest, n = sums(res)
    (run / "DONE").write_text(f"DONE {slug} 2026-11-23T07:31:00Z commit {commit} results {n} SHA256SUMS {digest}\n")


def combine_into(reading, top):
    argv = ["--contrasts", str(top / "contrasts.tsv"), "--out", str(reading / "combined")]
    for slug in sorted(p.name for p in (reading / "runs").iterdir()):
        argv += ["--run", str(reading / "runs" / slug)]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        combine.main(argv)


def lay_reading(top, tamper=None, species=None, commit_of=None, after=None):
    """A SYNTHETIC reading directory under top/reading, with lists.tsv and contrasts.tsv beside it.
    tamper(slug, tables) changes a run's tables before they are summed; commit_of(slug) gives a run's commit;
    after(reading) changes files once the combined tables are written."""
    top = Path(top)
    c.write_tsv(top / "lists.tsv", LISTS_HEADER, lists_rows())
    con = []
    for slug in SPECIES:
        for i, learner in enumerate(strata.LEARNERS):
            for contrast in ("total", "study"):
                con.append({"dataset": f"{slug}__synthomycin", "species": slug, "learner": learner,
                            "contrast": contrast, "d_recall": str(-0.01 * (i + 1)), "d_log_loss": str(0.1 * i)})
    c.write_tsv(top / "contrasts.tsv", CONTRASTS, con)
    reading = top / "reading"
    for k, slug in enumerate(species or SPECIES):
        tables = species_tables(slug, k)
        if tamper:
            tamper(slug, tables)
        lay_run(reading / "runs" / slug, slug, tables, commit=(commit_of or (lambda s: COMMIT))(slug))
    combine_into(reading, top)
    (reading / "reading.txt").write_text("SYNTHETIC reading commands\n")
    (reading / "checks.txt").write_text("SYNTHETIC checks\n")
    if after:
        after(reading)
    return reading


def write(top, reading):
    top = Path(top)
    return page.write(reading, top / "lists.tsv", top / "contrasts.tsv", COMMIT, top / "out")


class Gates(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)

    def refused(self, reading, gate):
        with self.assertRaises(page.Stop) as cm:
            write(self.top, reading)
        self.assertIn(f"gate {gate}", str(cm.exception))
        return str(cm.exception)

    def test_a_clean_reading_passes_every_gate_and_writes_the_record(self):
        files = write(self.top, lay_reading(self.top))
        want = {"README.md", "gates.txt", "reading.txt", "checks.txt"} | {f"combined/{n}" for n in page.COMBINED} \
            | {"page-species.tsv", "page-datasets.tsv", "page-m1-datasets.tsv", "page-m1-strata.tsv",
               "page-m1-pooled.tsv", "page-m2-cells.tsv", "page-m2-families.tsv", "page-m2-years.tsv",
               "page-m2-lineage.tsv", "page-timesplit.tsv"}
        self.assertEqual(set(files), want)
        out = self.top / "out"
        for line in (out / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ")
            self.assertEqual(hashlib.sha256((out / name).read_bytes()).hexdigest(), digest)
        for name in page.COMBINED:
            self.assertEqual((out / "combined" / name).read_bytes(),
                             (self.top / "reading" / "combined" / name).read_bytes())
        self.assertEqual(len((out / "gates.txt").read_text().splitlines()), len(page.GATES))

    def test_a_species_missing_from_the_reading_is_refused(self):
        reading = lay_reading(self.top, species=["synthia_alba", "synthia_borea"])
        msg = self.refused(reading, 2)
        self.assertIn("synthia_cana", msg)
        # a species more than A is refused too
        sub = self.top / "more"
        sub.mkdir()
        reading = lay_reading(sub)
        (reading / "runs" / "synthia_dubia").mkdir()
        with self.assertRaises(page.Stop) as cm:
            write(sub, reading)
        self.assertIn("gate 2", str(cm.exception))
        # combine.py's tables missing, or a file it does not write, stop at gate 1
        for case in ("missing", "extra"):
            sub = self.top / case
            sub.mkdir()
            reading = lay_reading(sub)
            if case == "missing":
                (reading / "combined" / "m2-summary.tsv").unlink()
            else:
                (reading / "combined" / "extra.tsv").write_text("x\n")
            with self.assertRaises(page.Stop) as cm:
                write(sub, reading)
            self.assertIn("gate 1", str(cm.exception), case)

    def test_a_dataset_of_m_marked_primary_is_refused(self):
        def tamper(slug, tables):
            for name in ("m1/datasets.tsv", "m1/strata.tsv", "m1/pooled.tsv"):
                for r in tables[name][1]:
                    if r.get("dataset", "") in IN_M or (slug == "synthia_cana" and name == "m1/pooled.tsv"
                                                        and r["set"] == "M"):
                        r["set"], r["labels_depend_on"] = "primary", ""
        self.refused(lay_reading(self.top, tamper=tamper), 4)

        # strata.tsv alone puts the dataset of M in the primary set; datasets.tsv keeps it in M
        def strata_only(slug, tables):
            for r in tables["m1/strata.tsv"][1]:
                if r["dataset"] in IN_M:
                    r["set"], r["labels_depend_on"] = "primary", ""
        sub = self.top / "strata-only"
        sub.mkdir()
        with self.assertRaises(page.Stop) as cm:
            write(sub, lay_reading(sub, tamper=strata_only))
        self.assertIn("gate 4", str(cm.exception))

        # a dataset of M that waits for another label rule than lists.tsv gives
        def other_rule(slug, tables):
            for r in tables["m1/datasets.tsv"][1]:
                if r["dataset"] in IN_M:
                    r["labels_depend_on"] = "another rule"
        sub = self.top / "other-rule"
        sub.mkdir()
        with self.assertRaises(page.Stop) as cm:
            write(sub, lay_reading(sub, tamper=other_rule))
        self.assertIn("gate 4", str(cm.exception))

    def test_the_primary_reading_holds_no_dataset_of_m(self):
        write(self.top, lay_reading(self.top))
        m1 = {(r["split"], r["set"], r["row"]): r for r in c.read_tsv(self.top / "out" / "combined" / "m1-reading.tsv")}
        # primary: 5 datasets, of which synthia_borea__fictocillin has 9 point only resistant genomes
        self.assertEqual(m1[("lineage", "primary", "scm")]["datasets"], "4")
        self.assertEqual(m1[("lineage", "M", "scm")]["datasets"], "1")
        rows = c.read_tsv(self.top / "out" / "page-m1-datasets.tsv")
        self.assertEqual({r["dataset"] for r in rows if r["set"] == "M"}, IN_M)

    def test_a_first_run_copy_is_refused(self):
        reading = lay_reading(self.top, commit_of=lambda s: FIRST_RUN if s == "synthia_borea" else COMMIT)
        self.assertIn("synthia_borea", self.refused(reading, 3))
        # a DONE that names another species is refused
        sub = self.top / "other-done"
        sub.mkdir()
        reading = lay_reading(sub)
        done = reading / "runs" / "synthia_alba" / "DONE"
        done.write_text(done.read_text().replace("DONE synthia_alba ", "DONE synthia_borea "))
        with self.assertRaises(page.Stop) as cm:
            write(sub, reading)
        self.assertIn("gate 3", str(cm.exception))

    def test_a_table_changed_after_its_run_is_refused(self):
        def after(reading):
            p = reading / "runs" / "synthia_alba" / "results" / "m2" / "m2-years.tsv"
            p.write_text(p.read_text().replace("\t40\n", "\t41\n"))
        self.assertIn("m2/m2-years.tsv", self.refused(lay_reading(self.top, after=after), 3))

    def test_a_share_in_a_cell_under_its_floor_is_refused(self):
        def tamper(slug, tables):
            if slug == "synthia_alba":
                for r in tables["m2/m2-cells.tsv"][1]:
                    if r["class"] == "TETRACYCLINE" and r["guard"] == "bioproject":
                        r["earlier_point"] = "0.250000000000"
        self.assertIn("synthia_alba/m2/m2-cells.tsv", self.refused(lay_reading(self.top, tamper=tamper), 5))

    def test_a_recall_under_ten_resistant_is_refused(self):
        def tamper(slug, tables):
            if slug == "synthia_cana":
                for r in tables["m1/strata.tsv"][1]:
                    if r["stratum"] == "both":
                        r["recall"] = "1.000000000000"
        self.assertIn("synthia_cana/m1/strata.tsv", self.refused(lay_reading(self.top, tamper=tamper), 5))

    def test_combined_tables_not_from_the_runs_are_refused(self):
        def after(reading):
            p = reading / "combined" / "m1-reading.tsv"
            p.write_text(p.read_text().replace("gene only higher", "point only higher"))
        reading = lay_reading(self.top, after=after)
        if "point only higher" not in (reading / "combined" / "m1-reading.tsv").read_text():
            (reading / "combined" / "m1-reading.tsv").write_text(
                (reading / "combined" / "m1-reading.tsv").read_text() + "\n")
        self.refused(reading, 6)

    def test_a_full_output_folder_is_refused(self):
        (self.top / "out").mkdir()
        (self.top / "out" / "x").write_text("x")
        with self.assertRaises(page.Stop):
            write(self.top, lay_reading(self.top))


class Tables(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.top = Path(cls.tmp.name)
        write(cls.top, lay_reading(cls.top))
        cls.out = cls.top / "out"

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def rows(self, name):
        return c.read_tsv(self.out / name)

    def test_species_and_datasets_are_named_in_full(self):
        self.assertEqual([r["species"] for r in self.rows("page-species.tsv")], list(SPECIES.values()))
        for name in ("page-datasets.tsv", "page-m1-datasets.tsv", "page-m2-cells.tsv", "page-m2-families.tsv",
                     "page-m2-years.tsv", "page-timesplit.tsv", "page-m1-pooled.tsv"):
            got = {r["species"] for r in self.rows(name)}
            self.assertEqual(got, set(SPECIES.values()), name)
        ds = self.rows("page-datasets.tsv")
        self.assertEqual(len(ds), 6)
        self.assertNotIn("synthia_alba__nullomycin", {r["dataset"] for r in ds})
        self.assertEqual({r["antibiotic"] for r in ds}, set(DRUGS))
        printed = "\n".join(page.tables(self.out)["datasets"])
        self.assertIn("| *Synthia cana* | fictocillin | M | yes |", printed)

    def test_a_lineage_guard_row_keeps_its_own_figures(self):
        combined = {(r["species_slug"], r["cell"], r["class"], r["dataset"], r["periods"], r["guard"], r["statistic"]): r
                    for r in c.read_tsv(self.out / "combined" / "m2-reading.tsv")}
        cells = self.rows("page-m2-cells.tsv")
        lineage = [r for r in cells if r["guard"] == "lineage" and r["passes"] == "1"]
        self.assertEqual(len(lineage), 6)
        for r in cells:
            key = ("synthia_" + r["species"].split()[1], r["cell"], r["class"], r["dataset"], r["periods"], r["guard"])
            for stat in ("point", "tvd"):
                self.assertEqual(r[f"holm_p_{stat}"], combined[key + (stat,)]["holm_p"], (key, stat))
                self.assertEqual(r[f"moves_{stat}"], combined[key + (stat,)]["moves"], (key, stat))
        bp = {(r["species"], r["class"], r["periods"]): r["holm_p_point"] for r in cells
              if r["guard"] == "bioproject" and r["passes"] == "1"}
        for r in lineage:
            self.assertNotEqual(r["holm_p_point"], bp[(r["species"], r["class"], r["periods"])])
        # with the mix's p apart from the point share's, each statistic keeps its own Holm p and moves
        def tamper(slug, tables):
            for r in tables["m2/m2-cells.tsv"][1]:
                if r["guard"] == "bioproject" and r.get("p_tvd"):
                    r["p_tvd"] = "0.001000000000"
        with tempfile.TemporaryDirectory() as tmp:
            write(tmp, lay_reading(tmp, tamper=tamper))
            dec = [r for r in c.read_tsv(Path(tmp) / "out" / "page-m2-cells.tsv")
                   if (r["cell"], r["periods"], r["guard"]) == page.DECIDING and r["passes"] == "1"]
            self.assertEqual(len(dec), 3)
            self.assertEqual({(r["moves_point"], r["moves_tvd"]) for r in dec}, {("0", "1")})
            self.assertTrue(all(r["holm_p_tvd"] != r["holm_p_point"] for r in dec))

    def test_the_deciding_count_reads_the_bioproject_vote_alone(self):
        # the BioProject vote's p is 0.5 in every cell and the largest-out guard's 0.001: only a guard moves
        sp = self.rows("page-species.tsv")
        self.assertEqual([(r["m2_classes"], r["m2_pass_primary"], r["m2_moves_point"], r["m2_moves_tvd"]) for r in sp],
                         [("2", "1", "0", "0")] * 3)
        largest = [r for r in self.rows("page-m2-cells.tsv") if r["guard"] == "largest_out" and r["passes"] == "1"]
        self.assertTrue(largest and all(r["moves_point"] == "1" and r["decides"] == "0" for r in largest))
        printed = page.tables(self.out)["m2_cells"]
        self.assertEqual(len(printed), 2 + 3)
        self.assertTrue(all(line.endswith("| no |") for line in printed[2:]))
        summary = page.tables(self.out)["m2_summary"]
        self.assertIn("| 2010 to 2014 against 2015 to 2019 | one vote per BioProject (decides) | point | 6 | 3 | 0 | yes |",
                      summary)
        self.assertIn("| 2010 to 2014 against 2015 to 2019 | largest BioProject left out | point | 6 | 3 | 3 | no |",
                      summary)

    def test_a_cell_that_does_not_pass_prints_no_share(self):
        printed = "\n".join(page.tables(self.out)["m2_cells"])
        self.assertNotIn("TETRACYCLINE", printed)
        tet = [r for r in self.rows("page-m2-cells.tsv") if r["class"] == "TETRACYCLINE"]
        self.assertTrue(tet and all(r["earlier_point"] == r["later_point"] == r["tvd"] == "" for r in tet))

    def test_gene_less_point_needs_both_strata_at_ten(self):
        rows = {(r["dataset"], r["split"], r["row"]): r for r in self.rows("page-m1-datasets.tsv")}
        short = rows[("synthia_borea__fictocillin", "lineage", "scm")]
        self.assertEqual((short["point_only_resistant"], short["point_only_recall"], short["gene_less_point"]),
                         ("9", "", ""))
        full = rows[("synthia_alba__synthomycin", "lineage", "scm")]
        self.assertEqual(full["gene_less_point"], f"{round(12 / 12 - 11 / 11, 12):.12f}")
        cart = rows[("synthia_borea__synthomycin", "lineage", "cart")]
        self.assertEqual(cart["gene_less_point"], f"{round(11 / 12 - 9 / 11, 12):.12f}")

    def test_strata_counts_sum_over_entered_datasets_with_no_pooled_recall(self):
        rows = {(r["set"], r["split"], r["row"], r["stratum"]): r for r in self.rows("page-m1-strata.tsv")}
        gene = rows[("primary", "lineage", "program", "gene_only")]
        self.assertEqual((gene["datasets"], gene["with_resistant"], gene["ten_or_more"], gene["resistant"]),
                         ("5", "5", "5", "60"))
        point = rows[("primary", "lineage", "program", "point_only")]
        self.assertEqual((point["ten_or_more"], point["resistant"]), ("4", str(11 * 4 + 9)))
        self.assertNotIn("recall", page.M1_STRATA_COLUMNS)
        self.assertEqual(rows[("M", "lineage", "program", "gene_only")]["datasets"], "1")
        # neither holds exactly 10 resistant genomes in each dataset, which meets the floor of 10
        self.assertEqual(rows[("primary", "lineage", "program", "neither")]["ten_or_more"], "5")
        # a species' count needs both strata at 10: synthia_borea__fictocillin's point only stratum holds 9
        sp = {r["species"]: r for r in self.rows("page-species.tsv")}
        self.assertEqual(sp["Synthia borea"]["m1_primary_both_ten"], "1")

    def test_the_reading_table_gives_the_five_learners_then_the_two_callers_beside(self):
        printed = page.tables(self.out)["m1_reading"]
        self.assertEqual(len(printed), 2 + 7)
        self.assertTrue(printed[2].startswith("| set covering machine (`scm`) | 4 | 3 |"))
        self.assertTrue(printed[7].startswith("| GRM's known-gene rule |"))
        self.assertTrue(printed[8].endswith("| beside | decides nothing |"))

    def test_print_reads_the_record_alone(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(page.main(["print", str(self.out)]), 0)
        for name in ("m1_reading", "m1_beside", "m1_neither", "m1_strata", "m2_summary", "m2_secondary", "m2_cells",
                     "species", "datasets"):
            self.assertIn(f"[{name}]", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
