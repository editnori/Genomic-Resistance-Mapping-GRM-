"""The committed outside-study results (reference/outside-study) and the document that reads them.

The suite scores nothing. It checks that the document's tables are the ones the
output files give, that the files are the listed ones, and that each run's
record names the frozen kit, its own population and the frozen genomes.
"""

import csv
import gzip
import hashlib
import importlib.util
import json
import pathlib
import shutil
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reference" / "outside-study"
DOCUMENT = ROOT / "docs" / "outside-study-results.md"
FROZEN = ROOT / "scripts" / "outside-study" / "frozen"
KIT = "e3e1c345a55412ac8f7d07b691ebb56573d2d3d8"
NAME = "a retrospective test on newly deposited studies from the same database"
FILES = sorted(
    [
        "datasets.tsv",
        "duplicates.tsv",
        "features.tsv",
        "genomes.tsv",
        "groups.tsv",
        "leave_one_study_out.tsv",
        "metrics.tsv",
        "paired.tsv",
        "reading_counts.tsv",
        "reading.txt",
        "results.json.gz",
    ]
)
# docs/outside-study-test.md, item 1: the value of `population` in the frozen
# candidates and the flag in the frozen dataset table, per population.
POPULATIONS = {"primary": ("primary (other-study)", "eligible_primary"), "same-study": ("same-study", "eligible_same_study")}
STAGES = ["count", "fetch", "known-gene", "predict", "setup", "species-check"]


def load_tables():
    spec = importlib.util.spec_from_file_location("outside_study_results_tables", ROOT / "scripts" / "outside-study-results" / "tables.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_tables()


def tsv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def record(run):
    with gzip.open(run / "results.json.gz", "rt", encoding="utf-8") as f:
        return json.load(f)


def runs():
    return sorted(p for pop in POPULATIONS for p in (RESULTS / pop).iterdir() if p.is_dir())


class Tables(unittest.TestCase):
    def test_the_document_holds_the_tables_the_files_give(self):
        self.assertTrue(T.check(RESULTS, DOCUMENT))

    def test_a_changed_cell_in_the_document_fails_the_check(self):
        text = DOCUMENT.read_text(encoding="utf-8")
        begin, end = text.index(T.BEGIN), text.index(T.END)
        block = text[begin:end]
        self.assertIn(" below ", block)
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "changed.md"
            changed.write_text(text[:begin] + block.replace(" below ", " above ", 1) + text[end:], encoding="utf-8")
            self.assertFalse(T.check(RESULTS, changed))
            unmarked = pathlib.Path(tmp) / "unmarked.md"
            unmarked.write_text(text.replace(T.BEGIN, ""), encoding="utf-8")
            self.assertFalse(T.check(RESULTS, unmarked))

    def test_a_changed_reading_in_an_output_file_fails_the_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            copy = pathlib.Path(tmp) / "outside-study"
            shutil.copytree(RESULTS, copy)
            metrics = sorted(copy.glob("primary/*/metrics.tsv"))[0]
            text = metrics.read_text(encoding="utf-8")
            self.assertIn("\tbelow\t", text)
            metrics.write_text(text.replace("\tbelow\t", "\tabove\t", 1), encoding="utf-8")
            self.assertFalse(T.check(copy, DOCUMENT))

    def test_the_two_populations_are_printed_apart(self):
        text = T.build(RESULTS)
        self.assertEqual(text.count("**Primary population**"), 1)
        self.assertEqual(text.count("**Same-study population**"), 1)
        self.assertLess(text.index("**Primary population**"), text.index("**Same-study population**"))
        primary = text[: text.index("**Same-study population**")]
        frozen = sum(T.frozen_and_scored((run / "reading.txt").read_text(encoding="utf-8"))[0] for run in runs() if run.parent.name == "primary")
        self.assertIn(f"| All {len(list((RESULTS / 'primary').iterdir()))} | {frozen:,} |", primary)

    def test_a_dataset_in_both_populations_is_counted_once(self):
        text = T.build(RESULTS)
        block = text.split("**Datasets scored in both populations**", 1)[1].split("The same rows by the estimate alone", 1)
        self.assertEqual(len(block), 2)
        tables = []
        for part in block:
            table = part.split("**Rows cited in the text**", 1)[0]
            got = [[c.strip() for c in l.strip("|").split("|")] for l in table.split("\n") if l.startswith("| ") and not l.startswith("| Learner")]
            self.assertEqual([g[0] for g in got], [T.NAMES[l] for l in T.LEARNERS])
            tables.append(got)
        # The four ways a row can hold add up to the rows read in both.
        for g in tables[0]:
            self.assertEqual(int(g[1]), sum(int(x) for x in g[2:]), g)
        # Each measure's three counts add up to the same rows.
        for g, first in zip(tables[1], tables[0]):
            self.assertEqual(g[1], first[1])
            for cell in g[2:]:
                self.assertEqual(int(g[1]), sum(int(x) for x in cell.split(" / ")), g)

    def test_a_cited_row_that_is_missing_is_refused(self):
        data = T.load(RESULTS)
        with self.assertRaises(SystemExit):
            T.cited_row(data, "primary", "staphylococcus_aureus", "no_such_drug", "scm")

    def test_every_cited_figure_is_the_kit_s_own_text(self):
        # reading.txt gives each estimate in plain words. The percentages this
        # script formats from metrics.tsv must be the ones written there.
        data = T.load(RESULTS)
        for pop, species, drug, learner in T.CITED:
            row = T.cited_row(data, pop, species, drug, learner)
            lines = data[pop][species][0].split("\n")
            starts = [i for i, l in enumerate(lines) if l.startswith(f"{row['dataset']}: ")]
            self.assertEqual(len(starts), 1, (pop, species, drug))
            # A dataset's block is its header line and the indented lines under it.
            block = []
            for l in lines[starts[0] + 1 :]:
                if not l.startswith(" "):
                    break
                block.append(l)
            line = [l for l in block if l.startswith(f"  {learner}: ")]
            self.assertEqual(len(line), 1, (pop, species, drug, learner))
            for k in ("sensitivity", "specificity"):
                want = f"{k.capitalize()} {T.pct(row[k])} [{T.pct(row[k + '_lo'])}, {T.pct(row[k + '_hi'])}]"
                self.assertIn(want, line[0], (pop, species, drug, learner))
            self.assertIn(f"TP {row['tp']}, FN {row['fn']}, FP {row['fp']}, TN {row['tn']}.", line[0])


class Files(unittest.TestCase):
    def test_every_run_folder_holds_the_eleven_output_files(self):
        self.assertTrue(runs())
        for run in runs():
            self.assertEqual(sorted(p.name for p in run.iterdir()), FILES, run)

    def test_no_folder_but_the_two_populations(self):
        self.assertEqual(sorted(p.name for p in RESULTS.iterdir()), ["SHA256SUMS", "primary", "same-study"])

    def test_sha256sums_lists_every_file_with_its_hash(self):
        listed = {}
        for line in (RESULTS / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            self.assertNotIn(name, listed)
            listed[name] = digest
        on_disk = sorted(str(p.relative_to(RESULTS)) for p in RESULTS.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        self.assertEqual(sorted(listed), on_disk)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RESULTS / name).read_bytes()).hexdigest(), digest, name)


class Records(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.candidates = tsv(FROZEN / "candidates.tsv.gz")
        cls.datasets = tsv(FROZEN / "datasets.tsv")

    def test_each_record_names_the_kit_and_its_own_run(self):
        binaries = set()
        for run in runs():
            r = record(run)
            self.assertEqual(r["kit"]["commit"], KIT, run)
            self.assertEqual(r["kit"]["expect_commit"], KIT, run)
            self.assertEqual(r["name"], NAME, run)
            self.assertEqual((r["population"], r["species"]), (run.parent.name, run.name))
            self.assertIs(r["rehearsal"], False, run)
            self.assertIs(r["table_pins_checked"], True, run)
            self.assertEqual(r["evaluation_revision"], 12, run)
            self.assertEqual(sorted(r["attempts"]), STAGES, run)
            binaries.add(r["grm_sha256"])
        self.assertEqual(len(binaries), 1, "every run used one grm binary")

    def test_each_run_covers_the_frozen_genomes_and_datasets(self):
        for run in runs():
            pop, species = run.parent.name, run.name
            value, flag = POPULATIONS[pop]
            frozen = sorted(c["genome_id"] for c in self.candidates if c["species_slug"] == species and c["population"] == value and c["on_eligible_dataset"] == "yes")
            eligible = sorted(d["dataset_id"] for d in self.datasets if d["species_slug"] == species and d[flag] == "yes")
            self.assertTrue(frozen and eligible, run)
            self.assertEqual(sorted(g["genome_id"] for g in tsv(run / "genomes.tsv")), frozen, run)
            self.assertEqual(sorted(d["dataset"] for d in tsv(run / "datasets.tsv")), eligible, run)
            self.assertEqual(record(run)["genomes_frozen"], len(frozen), run)

    def test_no_species_of_a_frozen_population_is_missing(self):
        for pop, (value, flag) in POPULATIONS.items():
            want = sorted({d["species_slug"] for d in self.datasets if d[flag] == "yes"})
            self.assertEqual(sorted(p.name for p in (RESULTS / pop).iterdir()), want, pop)

    def test_the_counts_in_metrics_are_the_record_s(self):
        for run in runs():
            by_key = {(m["dataset"], m["learner"]): m for m in record(run)["metrics"]}
            rows = tsv(run / "metrics.tsv")
            self.assertEqual(len(rows), len(by_key), run)
            for row in rows:
                m = by_key[(row["dataset"], row["learner"])]
                for k in ("tp", "fn", "fp", "tn"):
                    self.assertEqual(int(row[k]), m[k], (run, row["dataset"], row["learner"], k))
                self.assertEqual(row["reading_sensitivity"], m["reading_sensitivity"])
                self.assertEqual(row["reading_specificity"], m["reading_specificity"])


if __name__ == "__main__":
    unittest.main()
