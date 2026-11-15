"""The committed record under reference/tie-rule: its files against each other and against SHA256SUMS.

The reports and traces it was written from are not in the repository, so a
row file cannot be checked against them here; `record.py --check` does that
where they are unpacked. What is checked here is that every derived file is
what the script derives from the row files, and that the counts agree.
"""

import hashlib
import importlib.util
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
DIR = ROOT / "reference" / "tie-rule"
SCRIPT = ROOT / "scripts" / "tie-rule" / "record.py"


def load():
    spec = importlib.util.spec_from_file_location("tie_rule_record_committed", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load()
HEX = re.compile(r"[0-9a-f]{64}")
SORT = {
    "datasets.tsv": lambda r: r["dataset"],
    "measures.tsv": lambda r: (r["dataset"], r["split"], R.LEARNERS.index(r["learner"])),
    "seconds.tsv": lambda r: (r["dataset"], R.LEARNERS.index(r["learner"])),
    "traces.tsv": lambda r: (r["dataset"], R.KINDS.index(r["kind"])),
    "reproduction.tsv": lambda r: r["dataset"],
    "splits.tsv": lambda r: (r["dataset"], r["split"]),
    "summary.tsv": lambda r: r["key"],
}


def text(name):
    return (DIR / name).read_text(encoding="utf-8")


class Record(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = R.read_rows(DIR)
        cls.derived = R.derive(R.read_rows(DIR))
        cls.summary = {r["key"]: r["value"] for r in R.read_table(DIR / "summary.tsv")}
        cls.run_ids = [r["dataset"] for r in cls.rows["datasets"] if r["run"]]

    def test_the_checksums_cover_every_file_whole(self):
        lines = text("SHA256SUMS").splitlines()
        names = [line.split("  ", 1)[1] for line in lines]
        self.assertEqual(names, sorted(p.name for p in DIR.iterdir() if p.name != "SHA256SUMS"))
        for line in lines:
            digest, name = line.split("  ", 1)
            self.assertEqual(hashlib.sha256((DIR / name).read_bytes()).hexdigest(), digest, name)

    def test_the_files_are_the_script_s_and_small(self):
        want = set(R.COLUMNS) | {"PLAN.md", "PLAN.sha256", "SHA256SUMS", "README.md"}
        self.assertEqual({p.name for p in DIR.iterdir()}, want)
        self.assertLess(sum(p.stat().st_size for p in DIR.iterdir()), 3 * 1024 * 1024)
        for p in DIR.iterdir():
            data = p.read_text(encoding="utf-8")
            self.assertTrue(data.endswith("\n"), p.name)
            self.assertNotIn("\r", data, p.name)

    def test_each_table_has_its_columns_and_order(self):
        for name, cols in R.COLUMNS.items():
            if cols is None:
                continue
            lines = text(name).splitlines()
            self.assertEqual(tuple(lines[0].split("\t")), cols, name)
            self.assertTrue(all(len(line.split("\t")) == len(cols) for line in lines[1:]), name)
            if name in SORT:
                rows = R.read_table(DIR / name)
                keys = [SORT[name](r) for r in rows]
                self.assertEqual(keys, sorted(keys), name)
                self.assertEqual(len(keys), len(set(keys)), name)

    def test_the_derived_files_are_derived_from_the_row_files(self):
        for name in R.DERIVED_FILES:
            self.assertEqual(R.table(name, self.derived[name[: -len(".tsv")]]), text(name), name)

    def test_the_settings_are_the_script_s(self):
        self.assertEqual(R.table("settings.tsv", [dict(key=k, value=v) for k, v in R.SETTINGS]), text("settings.tsv"))

    def test_the_plan_has_its_sha256(self):
        digest, name = text("PLAN.sha256").split()
        self.assertEqual((name, hashlib.sha256((DIR / "PLAN.md").read_bytes()).hexdigest()), ("PLAN.md", digest))

    def test_the_run_datasets_are_the_same_in_every_row_file(self):
        run = self.run_ids
        self.assertEqual(len(run), 169)
        self.assertEqual(len(self.rows["datasets"]), 274)
        self.assertEqual(sorted({r["dataset"] for r in self.rows["measures"]}), run)
        self.assertEqual([r["dataset"] for r in self.rows["reproduction"]], run)
        self.assertEqual([r["dataset"] for r in self.rows["seconds"]], [d for d in run for _ in R.LEARNERS])
        self.assertEqual([(r["dataset"], r["kind"]) for r in self.rows["traces"]], [(d, k) for d in run for k in R.KINDS])
        want = sorted([f"{a}/{d}.json" for d in run for a in "ABD"] + [f"{a}/{d}.trace.jsonl" for d in run for a in "CE"])
        self.assertEqual([p for p, _ in self.rows["inputs"]], want)
        self.assertTrue(all(HEX.fullmatch(h) for _, h in self.rows["inputs"]))
        self.assertEqual(self.summary["inputs.files"], 5 * len(run))

    def test_every_arm_reproduced(self):
        for r in self.rows["reproduction"]:
            self.assertEqual([r[c] for c in R.COLUMNS["reproduction.tsv"][1:]], [True] * 4, r["dataset"])

    def test_each_learner_has_every_split_a_dataset_has(self):
        have = {}
        for r in self.rows["measures"]:
            have.setdefault((r["dataset"], r["split"]), []).append(r["learner"])
        self.assertTrue(all(v == list(R.LEARNERS) for v in have.values()))
        self.assertEqual(sum(1 for d, s in have if s == "lineage"), 163)
        self.assertEqual(sum(1 for d, s in have if s == "random"), 169)

    def test_the_changes_are_the_differences(self):
        for r in self.rows["measures"]:
            for m in R.PAIRED_MEASURES:
                self.assertEqual(r[f"{m}_change"], r[f"{m}_exact"] - r[f"{m}_present"], (r["dataset"], r["learner"], m))
            if not r["scores_changed"]:
                self.assertEqual((r["calls_changed"], r["log_loss_change"], r["metrics_changed"]), (0, 0.0, False), r["dataset"])
            if r["calls_changed"] and not r["metrics_changed"]:
                # two calls of one class that change in opposite directions leave the metrics as they were
                self.assertEqual((r["calls_changed"] % 2, r["correct_exact"]), (0, r["correct_present"]), r["dataset"])
            self.assertLessEqual(abs(r["correct_exact"] - r["correct_present"]), r["calls_changed"])

    def test_the_counts_per_learner_agree_with_the_measures(self):
        either = {}
        for r in self.rows["measures"]:
            if r["metrics_changed"]:
                either.setdefault(r["learner"], set()).add(r["dataset"])
        counts = {k: len(either.get(k, ())) for k in R.LEARNERS}
        self.assertEqual(counts, {"scm": 0, "cart": 61, "rf": 166, "gbt": 149, "lasso": 0})
        for k in R.LEARNERS:
            self.assertEqual(self.summary[f"{k}.either.metrics_changed"], counts[k])
            for s in R.SPLITS:
                rs = [r for r in self.rows["measures"] if (r["learner"], r["split"]) == (k, s)]
                self.assertEqual(self.summary[f"{k}.{s}.calls_changed"], sum(r["calls_changed"] for r in rs))
                self.assertEqual(self.summary[f"{k}.{s}.genomes"], sum(r["genomes"] for r in rs))
        self.assertEqual(self.summary["trees.calls_changed"],
                         sum(r["calls_changed"] for r in self.rows["measures"] if r["learner"] in R.TREES))
        self.assertEqual(self.summary["trees.either.log_loss_changed"], 169)

    def test_the_trace_sums_agree_with_the_traces(self):
        for kind in R.KINDS:
            rs = [r for r in self.rows["traces"] if r["kind"] == kind]
            for c in R.COLUMNS["traces.tsv"][2:]:
                self.assertEqual(self.summary[f"trace.{kind}.{c}"], sum(r[c] for r in rs), (kind, c))
        for r in self.rows["traces"]:
            self.assertEqual(r["fit_mismatches"], 0, r["dataset"])
            differ = r["differ_tie_present"] + r["differ_order_present"] + r["differ_partition_present"]
            self.assertLessEqual(differ, r["searched_present"])
            self.assertLessEqual(r["exact_tied_present"], r["searched_present"])
            self.assertLessEqual(r["trees_differ"], r["tree_positions"])
            if r["kind"] in ("cart", "rf"):
                self.assertEqual((r["differ_order_present"], r["differ_partition_present"]), (0, 0), r["dataset"])

    def test_the_datasets_not_run_are_the_four_species(self):
        unrun = [r for r in self.rows["datasets"] if not r["run"]]
        self.assertEqual(len(unrun), 105)
        species = sorted({r["species"] for r in unrun})
        self.assertEqual(",".join(species), self.summary["species.not_run.names"])
        self.assertEqual(self.summary["species.both"], 0)
        self.assertEqual((len(species), self.summary["species.run"]), (4, 19))

    def test_the_strings_are_ids_species_learners_and_kinds(self):
        ids = {r["dataset"] for r in self.rows["datasets"]}
        species = {r["species"] for r in self.rows["datasets"]}
        self.assertTrue(all(re.fullmatch(r"[a-z]+(_[a-z]+)+", s) for s in species))
        allowed = {"dataset": ids, "species": species, "learner": set(R.LEARNERS), "kind": set(R.KINDS),
                   "split": set(R.SPLITS), "first": set(R.LEARNERS), "second": set(R.LEARNERS), "measure": set(R.PAIRED_MEASURES),
                   "best_present": set(R.LEARNERS), "best_exact": set(R.LEARNERS)}
        for name, cols in R.COLUMNS.items():
            if cols is None or name in ("summary.tsv", "settings.tsv"):
                continue
            for r in R.read_table(DIR / name):
                for c, v in r.items():
                    if isinstance(v, str):
                        self.assertIn(v, allowed.get(c, ()), (name, c, v))
        for k, v in self.summary.items():
            if isinstance(v, str):
                self.assertTrue(v in ids or all(s in species for s in v.split(",")), (k, v))


if __name__ == "__main__":
    unittest.main()
