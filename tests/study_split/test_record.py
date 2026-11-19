"""reference/study-split: checksums whole, derived files equal derive(), each file's cells agree with the others, sizes and forms."""

import hashlib
import os
import pathlib
import re
import shutil
import tempfile
import unittest

from tests.study_split import agree
from tests.study_split.common import RECORD, load_script, read_tsv

record = load_script()
LIMIT = 4 * 1024 * 1024
# The files neither derive() nor agree.problems rebuilds whole, each pinned by its sha256. A sweep that changed one
# cell at three rows of every column of the row files found a change neither caught in nine of them; settings.tsv
# (its design keys), planfiles.tsv and inputs.sha256 are checked in part only. A change to one of these files is a
# change here too, and the README's Guards table names them from the script's PINNED.
PINNED = {
    "checks.tsv": "19b60b3e0d9c25cbd15a5a9c3733aa332eecdc7a4899606e90480cce05f6a58e",
    "distances.tsv": "eebeaf3cbd1e16849f4bef95d5d04850ba164c2772100ec3b2f027d5158c96f9",
    "eligibility.tsv": "b2119953453e75ce0590e0bd59fcad01971e653f4fbf3ef8ba3e6abcba95074f",
    "folds.tsv": "181ce61b5770226f84fb128fb00e856d0e02fdef16719845349011d6ba0887f5",
    "inputs.sha256": "1a008e7e1a51e4e7716f921e0b99f06690452cda215758a5f66fdde61dea1437",
    "lineage.tsv": "c3d5afd3ffd64e38d32ec3bdf6161f26c9ff05de8f3a0cef762508b19cf7d191",
    "order.tsv": "9f809da0fb7a2e4343f28dec0d8178652c2f359848d66c8791a77a616d44075e",
    "planfiles.tsv": "a22816705cdcaaba65f68b2ca223dac9d693d081fc71b4ac442efb70da3b8785",
    "plans.tsv": "2e200fef246ba0731674db16636926c34fe728bcc5ea89e7fe601f4a9df750cb",
    "settings.tsv": "1715e72c54e5774d39198c87b4b1504e62b651d20a3b05df6d9e87da063bd90b",
    "variants.tsv": "7feeedef407fb0dc6124ee6b5992f3cb7e014620cd34a125331cafd10f64ca9e",
    "variation.tsv": "3824c476e2ac346dab1a4a6f9f0e5be51bf82efbad09becd5efe2693284bdf3c",
}
# derive() computes about 46,000 exact intervals; worker processes share them out.
JOBS = max(1, min(4, os.cpu_count() or 1))


def rewrite_sum(directory, name):
    """Rewrite one file's line of SHA256SUMS, as a careful forger would."""
    d = pathlib.Path(directory)
    digest = hashlib.sha256((d / name).read_bytes()).hexdigest()
    lines = (d / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    lines = [f"{digest}  {name}" if x.split("  ", 1)[1] == name else x for x in lines]
    (d / "SHA256SUMS").write_text("".join(x + "\n" for x in lines), encoding="utf-8")


def change_cell(directory, name, key, column, new):
    """Set one cell of a TSV of the record, found by the first columns' values."""
    path = pathlib.Path(directory, name)
    lines = path.read_text(encoding="utf-8").split("\n")
    head = lines[0].split("\t")
    for i, line in enumerate(lines[1:-1], 1):
        cells = line.split("\t")
        if cells[:len(key)] == list(key):
            old = cells[head.index(column)]
            cells[head.index(column)] = new(old)
            assert cells[head.index(column)] != old
            lines[i] = "\t".join(cells)
            path.write_text("\n".join(lines), encoding="utf-8")
            return
    raise AssertionError(f"{name}: no row {key}")


class Record(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files, cls.tables = agree.load(RECORD)
        cls.text = {k: v.decode("utf-8") for k, v in cls.files.items()}
        rows = {name: cls.tables[name] for name in record.ROW_FILES}
        cls.derived = {name: record.tsv(name, r) for name, r in record.derive(rows, JOBS).items()}

    def failures(self, directory):
        """agree.problems, and each derived file that differs from derive() of the committed row files."""
        out = agree.problems(directory)
        for name in record.DERIVED:
            if pathlib.Path(directory, name).read_bytes() != self.derived[name]:
                out.append(f"{name}: differs from derive()")
        return out

    def test_checksums_whole(self):
        self.assertEqual(agree.checksum_problems(self.files), [])
        self.assertEqual(len(self.text["SHA256SUMS"].splitlines()), len(self.files) - 1)

    def test_under_four_megabytes(self):
        self.assertLess(sum(len(v) for v in self.files.values()), LIMIT)

    def test_files_are_the_scripts(self):
        names = set(record.COLUMNS) | {"README.md", "SHA256SUMS", "inputs.sha256"}
        self.assertEqual(set(self.files), names)
        self.assertEqual(set(record.ROW_FILES) | set(record.DERIVED), set(record.COLUMNS))
        self.assertEqual(set(record.DESCRIBE), set(record.COLUMNS))

    def test_columns_and_order(self):
        for name in record.COLUMNS:
            with self.subTest(name=name):
                self.assertEqual(read_tsv(RECORD / name, header_only=True), record.COLUMNS[name])
                rows = self.tables[name]
                self.assertTrue(rows)
                cols = record.COLUMNS[name]
                self.assertEqual(self.files[name], record.tsv(name, [[r[c] for c in cols] for r in rows]))
                if name in record.KEYS:
                    keys = [tuple(r[k] for k in record.KEYS[name]) for r in rows]
                    self.assertEqual(len(keys), len(set(keys)), "a key is written twice")

    def test_derived_files_equal_derive(self):
        for name in record.DERIVED:
            with self.subTest(name=name):
                self.assertEqual(self.files[name], self.derived[name])

    def test_cells_agree(self):
        self.assertEqual(agree.problems(RECORD), [])

    def test_pinned_files(self):
        self.assertEqual(sorted(PINNED), sorted(record.PINNED))
        for name, digest in PINNED.items():
            self.assertEqual(hashlib.sha256(self.files[name]).hexdigest(), digest, name)

    def test_guards_name_every_file(self):
        """The README's Guards table: every file under at least one guard, and Quoted the files the page reads."""
        guarded = set(record.DERIVED) | set(record.AGREED) | set(record.QUOTED) | set(record.PINNED)
        self.assertEqual(guarded | {"README.md", "SHA256SUMS"}, set(self.files))
        self.assertEqual(set(record.AGREED) & set(record.PINNED), set())
        self.assertLessEqual(set(record.ROW_FILES), set(record.AGREED) | set(record.PINNED))
        source = (pathlib.Path(__file__).parent / "lines.py").read_text(encoding="utf-8")
        self.assertEqual(set(record.QUOTED), set(re.findall(r'rd\(d / "([\w.]+\.tsv)"\)', source)))
        for name in record.PINNED:
            self.assertIn(f"`{name}`", self.text["README.md"].split("| Pinned |")[1].split("\n")[0])

    def test_readme_is_what_the_script_writes(self):
        settings = {r["key"]: r["value"] for r in self.tables["settings.tsv"]}
        self.assertEqual(self.files["README.md"], record.readme(self.files, settings))

    def test_plans_by_sha256(self):
        """The plans stay in the inputs: planfiles.tsv names each and gives its line of inputs.sha256."""
        inputs = dict(reversed(x.split("  ", 1)) for x in self.text["inputs.sha256"].splitlines())
        for exp in ("first", "second"):
            row = next(r for r in self.tables["planfiles.tsv"] if r["experiment"] == exp)
            self.assertEqual((row["file"], row["sha256"]), (f"{exp}/PLAN.md", inputs[f"{exp}/PLAN.md"]))
        self.assertFalse(any(name.startswith("plan/") for name in self.files))
        rule = next(r for r in self.tables["order.tsv"] if r["event"] == "stated_rule_choice")
        self.assertEqual(rule["source"], f"text of the second plan, sha256 {inputs['second/PLAN.md']}")

    def test_order_sources_inside_the_inputs(self):
        """Every source of order.tsv names a file of the input set and nothing outside it."""
        inputs = {x.split("  ", 1)[1] for x in self.text["inputs.sha256"].splitlines()}
        for r in self.tables["order.tsv"]:
            with self.subTest(event=r["event"]):
                named = re.findall(r"the input set's (\S+)", r["source"])
                if r["event"] == "stated_rule_choice":
                    self.assertEqual(named, [])
                else:
                    self.assertEqual(len(named), 1)
                    self.assertIn(named[0], inputs)
                    self.assertTrue(named[0].startswith(r["experiment"] + "/"))

    def test_holds_no_genome_label_or_score(self):
        barred = {"genome", "genomes_list", "label", "labels", "score", "scores", "prediction", "model", "feature",
                  "kmer", "importance", "weight", "study", "bioproject"}
        for name, cols in record.COLUMNS.items():
            self.assertFalse(barred & set(cols), name)
        # The largest row file has one row per dataset, plan and learner.
        self.assertEqual(len(self.tables["fits.tsv"]), len(self.tables["datasets.tsv"]) * 15 * 5)
        for name in record.COLUMNS:
            for r in self.tables[name]:
                for c in r.values():
                    self.assertNotRegex(c, r"^PRJ[NED][A-Z]\d+$", name)

    def test_no_genome_id_of_the_inputs(self):
        """The script counts the input tables' genome ids found as a token in the record, and writes the count."""
        ck = {r["check"]: r["count"] for r in self.tables["checks.tsv"]}
        self.assertEqual(ck["genome_ids_in_record"], "0")

    def test_a_changed_cell_fails_with_its_checksum_rewritten(self):
        d0 = self.tables["datasets.tsv"][0]["dataset"]
        bump = lambda x: str(int(x) + 1)
        nudge = lambda x: repr(float(x) + 1e-6)
        cases = [
            ("fits.tsv", (d0, "study-d0", "scm"), "tp", bump),
            ("fits.tsv", (d0, "random-d0", "rf"), "log_loss", nudge),
            ("lineage.tsv", (d0, "gbt"), "fp", bump),
            ("plans.tsv", (d0, "random-d2"), "lineage_mates", bump),
            ("plans.tsv", (d0, "matched-d1"), "resistant_0", bump),
            ("folds.tsv", (d0, "2", "lasso"), "random_log_loss", nudge),
            ("datasets.tsv", (d0,), "largest_study", bump),
            ("eligibility.tsv", (d0,), "with_study", lambda x: str(int(x) - 1)),
            ("checks.tsv", ("fold_files_checked",), "count", bump),
            ("checks.tsv", ("genome_ids_in_record",), "count", bump),
            ("order.tsv", ("second", "plan_builder"), "time", lambda x: x[:11] + "18:00:00Z"),
            ("summary.tsv", ("datasets",), "value", bump),
            ("contrasts.tsv", ("second", "gbt", "MR", "mean", "log_loss"), "reading", lambda x: "not shown"),
            ("draws.tsv", ("lasso", "MR", "matched", "cost"), "d4", bump),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for name, key, column, new in cases:
                with self.subTest(name=name, column=column):
                    copy = pathlib.Path(tmp, f"{name}-{column}")
                    shutil.copytree(RECORD, copy)
                    change_cell(copy, name, key, column, new)
                    rewrite_sum(copy, name)
                    files, _ = agree.load(copy)
                    self.assertEqual(agree.checksum_problems(files), [])
                    self.assertNotEqual(self.failures(copy), [])
                    shutil.rmtree(copy)

    def test_an_edited_plan_sha256_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            copy = pathlib.Path(tmp, "plan")
            shutil.copytree(RECORD, copy)
            change_cell(copy, "planfiles.tsv", ("second",), "sha256", lambda x: "0" + x[1:] if x[0] != "0" else "1" + x[1:])
            rewrite_sum(copy, "planfiles.tsv")
            self.assertIn("planfiles.tsv second: the plan's sha256 differs from inputs.sha256", agree.problems(copy))


if __name__ == "__main__":
    unittest.main()
