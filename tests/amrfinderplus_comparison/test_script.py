"""scripts/amrfinderplus-comparison/record.py on SYNTHETIC inputs: what it writes, what it refuses and where it
stops."""

import gzip
import hashlib
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction

from tests.amrfinderplus_comparison import synthetic
from tests.amrfinderplus_comparison.common import ROOT, SCRIPT, load_script, read_tsv

record = load_script()


def run(*args):
    done = subprocess.run([sys.executable, "-I", str(SCRIPT)] + [str(a) for a in args], cwd=ROOT,
                          capture_output=True, text=True)
    return done.returncode, done.stderr


def files(directory):
    d = pathlib.Path(directory)
    return {p.name: p.read_bytes() for p in sorted(d.iterdir())}


def edit_gz(path, old, new, count=1):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        text = fh.read()
    if old not in text:
        raise AssertionError("%r not in %s" % (old, path))
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write(text.replace(old, new, count))


def edit(path, old, new):
    text = pathlib.Path(path).read_text(encoding="utf-8")
    if old not in text:
        raise AssertionError("%r not in %s" % (old, path))
    pathlib.Path(path).write_text(text.replace(old, new, 1), encoding="utf-8")


class Script(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.top = pathlib.Path(cls.tmp) / "comparison"
        cls.final = synthetic.build(cls.top)
        cls.out = pathlib.Path(cls.tmp) / "record"
        code, err = run("--source", cls.final, "--out", cls.out)
        if code:
            raise AssertionError(err)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def copy_source(self, name):
        top = pathlib.Path(self.tmp) / name
        shutil.copytree(self.top, top)
        return top, top / "results" / "final"

    def test_writes_every_file(self):
        self.assertEqual(set(files(self.out)), set(record.ROW_FILES + record.DERIVED + ["README.md", "SHA256SUMS"]))
        for name in record.ROW_FILES + record.DERIVED:
            head = (self.out / name).read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
            self.assertEqual(head, record.COLUMNS[name], name)

    def test_same_bytes_twice_and_at_other_hash_seeds(self):
        again = pathlib.Path(self.tmp) / "again"
        self.assertEqual(run("--source", self.final, "--out", again)[0], 0)
        self.assertEqual(files(again), files(self.out))
        # Without -I, which would drop PYTHONHASHSEED; -s and -P keep site-packages and the cwd out.
        seen = set()
        for seed in ("0", "4242"):
            other = pathlib.Path(self.tmp) / ("seed" + seed)
            env = {"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"}
            subprocess.run([sys.executable, "-s", "-P", "-B", str(SCRIPT), "--source", str(self.final), "--out",
                            str(other)], check=True, env=env, cwd=ROOT, capture_output=True)
            self.assertEqual(files(other), files(self.out), seed)
            seen.add(subprocess.run([sys.executable, "-s", "-P", "-c", "print(hash('amrfinderplus comparison'))"],
                                    env=env, check=True, capture_output=True, text=True).stdout)
        self.assertEqual(len(seen), 2, "the two seeds gave one hash")

    def test_check_passes_on_its_own_output(self):
        self.assertEqual(run("--check", self.out), (0, ""))
        self.assertEqual(run("--source", self.final, "--check", self.out), (0, ""))

    def test_check_fails_on_a_changed_cell_a_dropped_file_and_an_extra_file(self):
        for name, change in (("callers.tsv", "cell"), ("summary.tsv", "cell"), ("intervals.tsv", "cell"),
                             ("figures.tsv", "drop"), ("README.md", "drop"), ("disagreements.tsv", "extra")):
            copy = pathlib.Path(self.tmp) / ("copy-%s-%s" % (name, change))
            shutil.copytree(self.out, copy)
            path = copy / name
            if change == "cell":
                lines = path.read_text(encoding="utf-8").split("\n")
                cells = lines[1].split("\t")
                cells[-1] = cells[-1] + "1"
                lines[1] = "\t".join(cells)
                path.write_text("\n".join(lines), encoding="utf-8")
            elif change == "drop":
                path.unlink()
            else:
                path.write_text("set\tcount\nx\t1\n", encoding="utf-8")
            code, err = run("--check", copy)
            self.assertEqual(code, 1, (name, change))
            self.assertTrue(err, (name, change))

    def test_out_refuses_a_directory_with_other_files(self):
        other = pathlib.Path(self.tmp) / "busy"
        other.mkdir()
        (other / "notes.txt").write_text("not a record\n", encoding="utf-8")
        code, err = run("--source", self.final, "--out", other)
        self.assertEqual(code, 1)
        self.assertIn("not an earlier record", err)
        self.assertEqual(sorted(p.name for p in other.iterdir()), ["notes.txt"])
        # An empty directory and an earlier record of the script are written over.
        empty = pathlib.Path(self.tmp) / "empty"
        empty.mkdir()
        self.assertEqual(run("--source", self.final, "--out", empty)[0], 0)
        self.assertEqual(run("--source", self.final, "--out", empty)[0], 0)
        self.assertEqual(files(empty), files(self.out))
        mixed = pathlib.Path(self.tmp) / "mixed"
        shutil.copytree(self.out, mixed)
        (mixed / "extra.tsv").write_text("a\n", encoding="utf-8")
        self.assertEqual(run("--source", self.final, "--out", mixed)[0], 1)

    def test_out_refuses_a_symlink_in_the_record(self):
        rec = pathlib.Path(self.tmp) / "with-link"
        shutil.copytree(self.out, rec)
        target = pathlib.Path(self.tmp) / "elsewhere.tsv"
        target.write_text("outside\n", encoding="utf-8")
        (rec / "callers.tsv").unlink()
        (rec / "callers.tsv").symlink_to(target)
        code, err = run("--source", self.final, "--out", rec)
        self.assertEqual(code, 1, err)
        self.assertIn("not an earlier record", err)
        self.assertTrue((rec / "callers.tsv").is_symlink())
        self.assertEqual(target.read_text(encoding="utf-8"), "outside\n")

    def test_out_refuses_a_symlink_as_the_record(self):
        rec = pathlib.Path(self.tmp) / "linked-record"
        shutil.copytree(self.out, rec)
        link = pathlib.Path(self.tmp) / "link"
        link.symlink_to(rec, target_is_directory=True)
        code, err = run("--source", self.final, "--out", link)
        self.assertEqual(code, 1, err)
        self.assertIn("a symlink", err)

    def test_one_datasets_measures_by_hand(self):
        r = record.rates({"tp": "3", "fn": "1", "tn": "1", "fp": "3"})
        self.assertEqual((r["sens"], r["spec"], r["ba"]), (Fraction(3, 4), Fraction(1, 4), Fraction(1, 2)))
        # An empty class gives 0.
        self.assertEqual(record.rates({"tp": "0", "fn": "0", "tn": "2", "fp": "0"})["ba"], Fraction(1, 2))

    def test_one_summary_cell_by_hand(self):
        d = record.describe([Fraction(1, 10), Fraction(1, 5), Fraction(0), Fraction(3, 10)])
        self.assertEqual((d["ahead"], d["level"], d["behind"]), (3, 1, 0))
        self.assertAlmostEqual(d["median"], 0.15, places=12)
        self.assertAlmostEqual(d["mean"], 0.15, places=12)
        self.assertAlmostEqual(d["q1"], 0.075, places=12)
        self.assertAlmostEqual(d["q3"], 0.225, places=12)
        self.assertEqual(d["sign_test_p"], "0.25")  # 3 ahead, 0 behind: 2 * (1/8)

    def test_out_needs_source(self):
        code, err = run("--out", pathlib.Path(self.tmp) / "nowhere")
        self.assertEqual(code, 2)
        self.assertIn("--out needs --source", err)

    def test_primary_reading_and_renamed_columns(self):
        diffs = read_tsv(self.out / "differences.tsv")
        self.assertEqual({r["reading"] for r in diffs}, {"tool_default"})
        with gzip.open(self.final / "differences.tsv.gz", "rt", encoding="utf-8") as fh:
            source = [line for line in fh.read().split("\n")[1:] if "\ttool_default\t" in line]
        self.assertEqual(len(diffs), len(source))
        for name in record.ROW_FILES + record.DERIVED:
            head = (self.out / name).read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
            self.assertNotIn("variant", head, name)
            self.assertNotIn("within_0.02", head, name)
            self.assertNotIn("share_within_0.02", head, name)
        for r in read_tsv(self.out / "summary.tsv"):
            if r["contrast"] == record.PR:
                self.assertEqual((r["learner_at_most_0.02_ahead"], r["within_0.02_either_way"]), ("", ""))
            else:
                self.assertNotEqual(r["learner_at_most_0.02_ahead"], "")

    def test_intervals_follow_the_species_count(self):
        sets = {r["set"]: r for r in read_tsv(self.out / "sets.tsv")}
        self.assertGreaterEqual(int(sets["primary-new-235"]["species"]), record.FEWEST_SPECIES)
        self.assertLess(int(sets["rule-M-quinolone-new-26"]["species"]), record.FEWEST_SPECIES)
        for r in read_tsv(self.out / "summary.tsv"):
            self.assertEqual(r["median_species_lo"] == "", int(r["species"]) < record.FEWEST_SPECIES, r["set"])

    def test_record_holds_no_seed_string_and_no_path(self):
        for name, data in files(self.out).items():
            text = data.decode("utf-8")
            self.assertNotIn(synthetic.SEED_STRING.split(":<")[0], text, name)
            self.assertNotIn(self.tmp, text, name)

    def test_source_stops(self):
        cases = [
            ("a summary cell", "results/final/cohort-summary.tsv.gz", None),
            ("an interval end", "results/final/cohort-summary.tsv.gz", "interval"),
            ("a plan list", "results/plan-lists/carbapenem-20.txt", None),
            ("the plan's hash", "PLAN.md", None),
            ("the seed string", "PLAN.md", "seed"),
            ("a genome id in a cell", "results/final/by-dataset-callers.tsv.gz", "genome"),
            ("a base run in a cell", "results/final/by-dataset-callers.tsv.gz", "bases"),
            ("a short genome id in a cell", "results/final/by-dataset-callers.tsv.gz", "short"),
            ("a rate above 1", "results/final/by-dataset-callers.tsv.gz", "rate"),
            ("the program version", "results/final/run-manifest.json", None),
            ("the final scoring's start", "results/logs/analyze-final.log", None),
        ]
        for label, rel, kind in cases:
            top, final = self.copy_source("stop-" + label.replace(" ", "-").replace("'", ""))
            path = top / rel
            if rel.endswith("cohort-summary.tsv.gz"):
                rows = read_tsv_gz(path)
                row = [r for r in rows if r["set"] == "primary-new-235" and r["contrast"] == record.PR
                       and r["measure"] == "ba"][0]
                col = "median_species_hi" if kind == "interval" else "first_ahead"
                old = "\t".join(row[c] for c in synthetic.COHORT_COLS)
                # An interval end moves by 0.01; a count gains a digit.
                row[col] = "%.6f" % (float(row[col]) + 0.01) if kind == "interval" else row[col] + "1"
                edit_gz(path, old, "\t".join(row[c] for c in synthetic.COHORT_COLS))
            elif rel.endswith(".txt"):
                text = path.read_text(encoding="utf-8")
                path.write_text("\n".join(text.split("\n")[1:]), encoding="utf-8")
            elif rel == "PLAN.md" and kind == "seed":
                text = path.read_text(encoding="utf-8").replace("seed string `" + synthetic.SEED_STRING, "seed `x")
                path.write_text(text, encoding="utf-8")
                sums = top / "PLAN.sha256"
                digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
                sums.write_text("%s  PLAN.md\n" % digest, encoding="utf-8")
            elif rel == "PLAN.md":
                edit(path, "A made-up plan", "A made-up plan, edited")
            elif kind == "genome":
                edit_gz(path, "\t0.100000\t", "\t1280.12345\t")
            elif kind == "bases":
                edit_gz(path, "\t0.100000\t", "\t" + "acgt" * 6 + "a\t")
            elif kind == "short":
                edit_gz(path, "\t0.100000\t", "\t1280.12\t")
            elif kind == "rate":
                edit_gz(path, "\t0.900000\t", "\t1.900000\t")
            elif rel.endswith(".json"):
                edit(path, "4.2.7", "4.2.8")
            else:
                edit(path, "scope final", "scope interim")
            out = pathlib.Path(self.tmp) / ("out-" + top.name)
            code, err = run("--source", final, "--out", out)
            self.assertEqual(code, 1, label)
            self.assertIn("stopped", err, label)
            self.assertFalse(out.exists(), label)
            if kind in ("genome", "bases", "short"):
                self.assertIn("base run", err, label)
            if kind == "rate":
                self.assertIn("outside 0 to 1", err, label)


def read_tsv_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        lines = fh.read().rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


if __name__ == "__main__":
    unittest.main()
