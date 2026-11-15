"""The final models of the reports the outside-study runs read, in two revisions.

scripts/outside-study-results/models.py writes one digest of each final model
of a release's reports. reference/outside-study-models holds its output for
the revision the runs pinned and for the revision the site serves. The suite
fails if a model's digest differs between the two files. It ties each file to
its reports by their hashes, and it rebuilds from the committed runs the
features of each model a run scored. The reports are not in the repository,
so the digest itself is tested on a made-up report.
"""

import contextlib
import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import pathlib
import re
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reference" / "outside-study"
MODELS = ROOT / "reference" / "outside-study-models"
OLDER = MODELS / "rev12-274-20261105.tsv"
LATER = MODELS / "rev14-274-20261112.tsv"
README = MODELS / "README.md"
FIGURES = ROOT / "reference" / "outside-study-figures" / "rev14-274-20261112.tsv"
# The release's own checksum lines for the reports of the served revision.
LATER_LINES = FIGURES.parent / "rev14-274-20261112.reports.sha256"
# The kit's pins: the reports the runs of record were allowed to read.
PINS = ROOT / "scripts" / "outside-study" / "reports.sha256"
RELEASES = ROOT / "reference" / "release-summary" / "README.md"
HEX = r"^[0-9a-f]{64}$"


def load_models():
    spec = importlib.util.spec_from_file_location("outside_study_results_models", ROOT / "scripts" / "outside-study-results" / "models.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_models()
R, T = M.revision, M.tables


def records():
    """The record of each run of record."""
    out = []
    for path in sorted(RESULTS.glob("*/*/results.json.gz")):
        with gzip.open(path, "rt", encoding="utf-8") as f:
            out.append(json.load(f))
    return out


def text_sha(text):
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def changed_copy(path, tmp, key, column, value):
    """A copy of a digest file with one cell changed."""
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    (row,) = [r for r in rows if (r["dataset"], r["learner"]) == key]
    row[column] = value
    out = pathlib.Path(tmp) / f"changed-{path.name}"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, M.COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return out


def run_main(*args):
    """(exit status, what the command printed)."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        status = M.main(["models.py", *map(str, args)])
    return status, out.getvalue()


class DigestFiles(unittest.TestCase):
    """The two committed files, one for each revision."""

    @classmethod
    def setUpClass(cls):
        cls.datasets = M.datasets(RESULTS)
        (cls.old_rev, cls.older), (cls.new_rev, cls.later) = M.load(OLDER), M.load(LATER)

    def test_each_file_has_six_rows_for_each_dataset_of_the_runs(self):
        want = {(d, learner) for d in self.datasets for learner in M.ROWS}
        self.assertEqual(M.ROWS, ("scm", "cart", "rf", "gbt", "lasso", "known_gene"))
        self.assertEqual((set(self.older), set(self.later)), (want, want))
        self.assertTrue(self.datasets)

    def test_the_files_are_of_the_revision_the_runs_pinned_and_the_one_the_site_serves(self):
        self.assertEqual({r["evaluation_revision"] for r in records()}, {self.old_rev})
        self.assertEqual(R.load_figures(FIGURES)[0], self.new_rev)
        self.assertLess(self.old_rev, self.new_rev)
        for path, revision in ((OLDER, self.old_rev), (LATER, self.new_rev)):
            self.assertTrue(path.name.startswith(f"rev{revision}-"), path.name)

    def test_no_final_model_differs_between_the_two_revisions(self):
        same, differ = M.compare(self.older, self.later)
        self.assertEqual(differ, [])
        self.assertEqual(len(same), len(M.ROWS) * len(self.datasets))
        status, printed = run_main(OLDER, LATER)
        self.assertEqual(status, 0)
        n = len(self.datasets)
        self.assertEqual(printed.split("\n")[0], f"revision {self.old_rev} and revision {self.new_rev}, {n} datasets")
        for learner in M.LEARNERS:
            self.assertIn(f"{learner}: the same final model on {n} of {n} datasets\n", printed)
        ruled = self.with_rule()
        self.assertIn(f"known_gene: the same rule on {len(ruled)} of {len(ruled)} datasets with a rule, and no rule in either release on {n - len(ruled)}\n", printed)
        self.assertIn(f"the five learners: {5 * n} of {5 * n} dataset and learner pairs\n", printed)
        self.assertNotIn("differs", printed)

    def with_rule(self):
        """The datasets whose reports publish a known-gene rule, by the digest files."""
        ruled = [{d for d in self.datasets if rows[d, "known_gene"]["model_sha256"] != M.no_rule()} for rows in (self.older, self.later)]
        self.assertEqual(ruled[0], ruled[1])
        return ruled[0]

    def test_the_runs_scored_the_known_gene_rule_where_a_report_publishes_one(self):
        scored = {r["dataset"] for runs in T.load(RESULTS).values() for _, _, metrics, _ in runs.values() for r in metrics if r["learner"] == "known_gene"}
        self.assertEqual(self.with_rule(), scored)
        self.assertTrue(0 < len(scored) < len(self.datasets))

    def test_each_dataset_has_one_input_fingerprint_in_both_revisions(self):
        """An id of a feature a report does not list stands for one row of one matrix in both."""
        self.assertEqual(M.same_input(self.older, self.later), self.datasets)
        n = len(self.datasets)
        self.assertTrue(run_main(OLDER, LATER)[1].endswith(f"input fingerprint: the same on {n} of {n} datasets\n"))
        for rows in (self.older, self.later):
            for row in rows.values():
                self.assertRegex(row["input_fingerprint"], r"^[0-9a-f]{16}$")
        key = (self.datasets[-1], "rf")
        with tempfile.TemporaryDirectory() as tmp:
            # A fingerprint is one value on the six rows of a dataset.
            with self.assertRaises(SystemExit):
                M.load(changed_copy(LATER, tmp, key, "input_fingerprint", "0" * 16))
            with open(LATER, newline="", encoding="utf-8") as f:
                text = f.read()
            changed = pathlib.Path(tmp) / "later.tsv"
            fingerprint = self.later[key]["input_fingerprint"]
            self.assertEqual(text.count(f"\t{fingerprint}\t"), len(M.ROWS))
            changed.write_text(text.replace(f"\t{fingerprint}\t", "\t" + "0" * 16 + "\t"), encoding="utf-8")
            status, printed = run_main(OLDER, changed)
            self.assertEqual(status, 1)
            self.assertIn(f"input fingerprint: the same on {n - 1} of {n} datasets\n", printed)
            self.assertTrue(printed.endswith(f"input fingerprint differs: {key[0]}\n"))
            # No digest differs there, so the count of models stands and the command still fails.
            self.assertIn(f"the five learners: {5 * n} of {5 * n} dataset and learner pairs\n", printed)

    def test_a_changed_digest_in_either_file_fails_the_comparison(self):
        key = (self.datasets[0], "lasso")
        with tempfile.TemporaryDirectory() as tmp:
            for path in (OLDER, LATER):
                changed = changed_copy(path, tmp, key, "model_sha256", "0" * 64)
                pair = (changed, LATER) if path == OLDER else (OLDER, changed)
                self.assertEqual(M.compare(M.load(pair[0])[1], M.load(pair[1])[1])[1], [key])
                status, printed = run_main(*pair)
                self.assertEqual(status, 1)
                n = len(self.datasets)
                self.assertIn(f"lasso: the same final model on {n - 1} of {n} datasets\n", printed)
                self.assertIn(f"the five learners: {5 * n - 1} of {5 * n} dataset and learner pairs\n", printed)
                self.assertIn(f"\ndiffers: {key[0]} lasso\n", printed)
                self.assertNotIn("input fingerprint differs", printed)

    def test_each_older_report_is_the_one_its_run_read_and_the_kit_pinned(self):
        read = {}
        for record in records():
            for dataset, digest in record["reports"].items():
                self.assertEqual(read.setdefault(dataset, digest), digest, dataset)
        pins = R.score.read_pins(PINS)
        self.assertEqual(sorted(read), self.datasets)
        for (dataset, learner), row in self.older.items():
            self.assertEqual(row["report_sha256"], read[dataset], (dataset, learner))
            self.assertEqual(row["report_sha256"], pins[f"{dataset}.json"], (dataset, learner))

    def test_each_later_report_is_its_line_in_the_release_s_checksum_list(self):
        lines = {}
        for line in LATER_LINES.read_text(encoding="utf-8").split("\n")[:-1]:
            digest, dataset = re.fullmatch(r"([0-9a-f]{64})  \./datasets/([a-z0-9_]+)\.json", line).groups()
            lines[dataset] = digest
        figures = R.load_figures(FIGURES)[1]
        self.assertEqual(sorted(lines), self.datasets)
        for (dataset, learner), row in self.later.items():
            self.assertEqual(row["report_sha256"], lines[dataset], (dataset, learner))
        # The figure file of the section Revision 14 was read from the same reports.
        for (dataset, _), figure in figures.items():
            self.assertEqual(figure["report_sha256"], self.later[dataset, "scm"]["report_sha256"], dataset)

    def test_no_report_is_in_both_releases(self):
        for dataset in self.datasets:
            self.assertNotEqual(self.older[dataset, "scm"]["report_sha256"], self.later[dataset, "scm"]["report_sha256"], dataset)

    def test_each_model_a_run_scored_has_the_features_the_run_recorded(self):
        """`features.tsv` of the runs holds what a run read of each model it scored."""
        recorded = M.run_features(RESULTS)
        scored = {(r["dataset"], r["learner"]) for runs in T.load(RESULTS).values() for _, _, metrics, _ in runs.values() for r in metrics if r["learner"] in R.score.MODEL_LEARNERS}
        self.assertEqual(set(recorded), scored)
        self.assertGreater(len(scored), 2 * len(self.datasets))
        for key, features in recorded.items():
            for name, rows in (("older", self.older), ("later", self.later)):
                self.assertEqual(rows[key]["read_features_sha256"], M.sha(features), (name, key))
        # A scored model reads at least one feature, and no two scored models of a dataset read the same ones.
        self.assertTrue(all(recorded.values()))

    def test_each_cell_has_the_form_the_script_writes(self):
        for name, rows in (("older", self.older), ("later", self.later)):
            for (dataset, learner), row in rows.items():
                self.assertRegex(row["report_sha256"], HEX)
                self.assertRegex(row["model_sha256"], HEX)
                if learner in R.score.MODEL_LEARNERS:
                    self.assertRegex(row["read_features_sha256"], HEX, (name, dataset, learner))
                else:
                    self.assertEqual(row["read_features_sha256"], "", (name, dataset, learner))
            # Each digest names its learner, so no two rows of a dataset share one.
            for dataset in self.datasets:
                self.assertEqual(len({rows[dataset, learner]["model_sha256"] for learner in M.ROWS}), len(M.ROWS), dataset)

    def test_a_file_is_read_back_as_written(self):
        for path in (OLDER, LATER):
            out = io.StringIO()
            w = csv.DictWriter(out, M.COLUMNS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            rows = M.load(path)[1]
            w.writerows(rows[d, learner] for d in self.datasets for learner in M.ROWS)
            self.assertEqual(out.getvalue().encode("utf-8"), path.read_bytes())

    def test_the_readme_names_each_file_and_its_release(self):
        text = README.read_text(encoding="utf-8")
        releases = RELEASES.read_text(encoding="utf-8").split("\n")
        for path, revision in ((OLDER, self.old_rev), (LATER, self.new_rev)):
            self.assertIn(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n", text)
            # The release is the one reference/release-summary names for this revision.
            (line,) = [x for x in releases if x.startswith(f"| `{path.stem}/` | ") and ".tar." in x]
            _, _, asset, sha, rev, _ = (c.strip().strip("`") for c in line.split("|")[:6])
            self.assertTrue(asset.endswith((".tar.gz", ".tar.zst")), asset)
            self.assertEqual(int(rev), revision)
            self.assertRegex(sha, HEX)
            self.assertIn(f"| `{path.name}` | `{asset}` | `{sha}` | {revision} | {len(self.datasets)} reports |", text)
        self.assertIn(f"{len(self.older)} rows", text)
        same, differ = M.compare(self.older, self.later)
        pairs = len(M.LEARNERS) * len(self.datasets)
        self.assertIn(
            f"On {len(same) - len(self.datasets)} of the {pairs} pairs of a dataset and a learner the two files have the same `model_sha256`, "
            f"and on {len(self.with_rule())} of the {len(self.with_rule())} datasets whose reports publish a known-gene rule, the rule has. "
            f"The other {len(self.datasets) - len(self.with_rule())} publish no rule in either release.",
            text,
        )
        self.assertIn(f"The two reports of a dataset have one `input_fingerprint` on {len(M.same_input(self.older, self.later))} of the {len(self.datasets)} datasets.", text)
        self.assertEqual(differ, [])
        self.assertIn(f"The runs scored {len(M.run_features(RESULTS))} of these models", text)

    def test_a_file_with_other_columns_rows_or_revisions_is_refused(self):
        text = OLDER.read_text(encoding="utf-8")
        lines = text.split("\n")
        first = self.datasets[0]
        self.assertTrue(lines[1].startswith(f"{first}\tscm\t") and lines[7].startswith(self.datasets[1] + "\t"))
        other_report = lines[2].replace(lines[2].split("\t")[2], "f" * 64, 1)
        other_input = lines[2].replace(lines[2].split("\t")[4], "f" * 16, 1)
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "digests.tsv"
            for name, changed in (
                ("other columns", text.replace("model_sha256", "sha256", 1)),
                ("a row twice", "\n".join(lines[:2] + lines[1:])),
                ("two revisions", "\n".join(lines[:1] + [lines[1].replace(f"\t{self.old_rev}\t", f"\t{self.new_rev}\t", 1)] + lines[2:])),
                ("a dataset lacks a learner", "\n".join(lines[:1] + lines[2:])),
                ("a dataset names two reports", "\n".join(lines[:2] + [other_report] + lines[3:])),
                ("a dataset names two fingerprints", "\n".join(lines[:2] + [other_input] + lines[3:])),
                ("another learner", "\n".join(lines[:7] + [lines[1].replace("\tscm\t", "\tsvm\t", 1)] + lines[7:])),
            ):
                path.write_text(changed, encoding="utf-8")
                with self.assertRaises(SystemExit, msg=name):
                    M.load(path)
            path.write_text(text, encoding="utf-8")
            self.assertEqual(M.load(path)[0], self.old_rev)

    def test_two_files_of_one_revision_or_of_other_rows_are_refused(self):
        with self.assertRaises(SystemExit) as stop:
            M.summary(OLDER, OLDER)
        self.assertIn("both files", str(stop.exception))
        fewer = {k: v for k, v in self.later.items() if k[0] != self.datasets[0]}
        with self.assertRaises(SystemExit):
            M.compare(self.older, fewer)

    def test_a_wrong_command_line_prints_the_usage(self):
        for args in ([], [OLDER], [OLDER, "--digests"], [OLDER, LATER, LATER], [OLDER, "--report", LATER], ["--fields", OLDER], [OLDER, "--reports"]):
            with self.assertRaises(SystemExit) as stop:
                M.main(["models.py", *map(str, args)])
            self.assertEqual(stop.exception.code, M.__doc__)


def made_up_report(dataset):
    """A report with the fields the digests read. The numbers are no dataset's."""

    def feature(i, kmer, n_kmers, bits):
        return {
            "feature": i,
            "kmer": kmer,
            "n_kmers": n_kmers,
            "bits": bits,
            "n_resistant": 4,
            "n_susceptible": 1,
            "catalog_families": ["made-up family"],
            "catalog_products": [],
            "locus": {"genome_id": "1.1", "gene": None, "product": "made-up product"},
        }

    return {
        "schema": R.score.REPORT_SCHEMA,
        "evaluation": {"revision": 12, "valid": True},
        "id": dataset,
        "input_fingerprint": "made-up-input",
        "seconds": 1.5,
        "splits": [{"scheme": "lineage", "learners": [{"learner": "scm", "metrics": {"n": 12}}]}],
        "catalog": {"version": "made-up-1", "rule": [["C:MADE-UP"]], "metrics": {"n": 12}, "groups": []},
        "features": [
            feature(3, "ACGTA", 7, "AAE="),
            feature(5, "CCGTA", 1, "AAI="),
            feature(8, "GGGTA", 2, "AAM="),
            feature(13, "TTGTA", 4, "AAQ="),
            feature(21, "AAGTC", 9, "AAU="),
            # Listed, and named by no model.
            feature(144, "CAGTC", 1, "AAY="),
        ],
        "models": {
            "scm": {
                "kind": "disjunction",
                "p": 1.0,
                "cv_error": 0.25,
                "rules": [{"feature": 3, "presence": True}, {"feature": 8, "presence": False}],
                "steps": [
                    # Feature 34 is an equivalent rule the report does not list.
                    {"covered": 5, "lost": 0, "utility": 5.0, "equivalent": [{"feature": 21, "presence": True}, {"feature": 34, "presence": True}], "n_equivalent": 2},
                    {"covered": 2, "lost": 1, "utility": 1.0, "equivalent": [], "n_equivalent": 0},
                ],
                "importance": [0.5, 0.125],
            },
            "cart": {
                "nodes": [
                    {"value": 0.5, "n": 12, "n_resistant": 6, "depth": 0, "split": {"feature": 5, "gain": 0.25, "absent": 1, "present": 2}},
                    # A split at the published depth: in the report, and not read by `grm predict`.
                    {"value": 0.25, "n": 8, "n_resistant": 2, "depth": 1, "split": {"feature": 13, "gain": 0.125, "absent": 3, "present": 4}},
                    {"value": 1.0, "n": 4, "n_resistant": 4, "depth": 1, "split": None},
                    {"value": 0.0, "n": 6, "n_resistant": 0, "depth": 2, "split": None},
                    {"value": 1.0, "n": 2, "n_resistant": 2, "depth": 2, "split": None},
                ]
            },
            "cart_depth": 1,
            # Feature 55 has a weight and is not listed.
            "lasso": {"intercept": -1.5, "coef": [[3, 2.0], [55, -0.75]], "lambda": 0.0625},
            "lasso_fit": {"path_points": 30, "selected_point": 4, "refit_convergence": {"irls_steps": 3, "converged": True}},
            "rf_importance": [{"feature": 5, "weight": 0.5}, {"feature": 89, "weight": 0.25}],
            "gbt_importance": [{"feature": 8, "weight": 3.0}],
        },
    }


class DigestsFromReports(unittest.TestCase):
    """The digests of made-up reports. The numbers are not any dataset's."""

    DATASET = "shigella_flexneri__ciprofloxacin"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        base = pathlib.Path(self.tmp.name)
        # One committed run with one dataset names the dataset to read. Its report here is made up.
        self.root, self.reports = base / "outside-study", base / "datasets"
        (self.root / "primary").mkdir(parents=True)
        (self.root / "primary" / "shigella_flexneri").symlink_to(RESULTS / "primary" / "shigella_flexneri")
        self.reports.mkdir()
        self.base = made_up_report(self.DATASET)

    def write(self, report, reports=None):
        raw = json.dumps(report).encode()
        (pathlib.Path(reports or self.reports) / f"{self.DATASET}.json").write_bytes(raw)
        return hashlib.sha256(raw).hexdigest()

    def rows(self, report):
        self.write(report)
        return {r["learner"]: r for r in csv.DictReader(io.StringIO(M.digests(self.root, self.reports)), delimiter="\t")}

    def digests(self, report):
        return {learner: row["model_sha256"] for learner, row in self.rows(report).items()}

    def changed(self, report):
        """The learners whose digest differs from the unchanged made-up report's."""
        before, after = self.digests(self.base), self.digests(report)
        return {learner for learner in M.ROWS if before[learner] != after[learner]}

    def test_the_run_names_one_dataset(self):
        self.assertEqual(M.datasets(self.root), [self.DATASET])

    def test_a_row_is_the_report_s_model_with_its_hash(self):
        digest = self.write(self.base)
        text = M.digests(self.root, self.reports)
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        self.assertEqual([r["learner"] for r in rows], list(M.ROWS))
        for r in rows:
            self.assertEqual((r["dataset"], r["report_sha256"], r["evaluation_revision"], r["input_fingerprint"]), (self.DATASET, digest, "12", "made-up-input"))
        by = {r["learner"]: r for r in rows}
        # The text a digest is taken of, written out for three rows.
        self.assertEqual(by["known_gene"]["model_sha256"], text_sha('{"features":[],"learner":"known_gene","model":{"rule":[["C:MADE-UP"]],"version":"made-up-1"}}'))
        self.assertEqual(
            by["gbt"]["model_sha256"],
            text_sha('{"features":[[8,"GGGTA",2,"AAM="]],"learner":"gbt","model":{"gbt_importance":[{"feature":8,"weight":3.0}]}}'),
        )
        self.assertEqual(
            by["rf"]["model_sha256"],
            text_sha('{"features":[[5,"CCGTA",1,"AAI="],[89,null,null,null]],"learner":"rf","model":{"rf_importance":[{"feature":5,"weight":0.5},{"feature":89,"weight":0.25}]}}'),
        )
        # The features the kit reads: every rule, the splits above the published depth, every weight.
        self.assertEqual(by["scm"]["read_features_sha256"], text_sha('[[3,"ACGTA",7],[8,"GGGTA",2]]'))
        self.assertEqual(by["cart"]["read_features_sha256"], text_sha('[[5,"CCGTA",1]]'))
        self.assertEqual(by["lasso"]["read_features_sha256"], text_sha('[[3,"ACGTA",7],[55,null,null]]'))
        self.assertEqual([by[learner]["read_features_sha256"] for learner in ("rf", "gbt", "known_gene")], ["", "", ""])
        # What it writes is a file this script reads.
        path = pathlib.Path(self.tmp.name) / "digests.tsv"
        path.write_text(text, encoding="utf-8")
        revision, loaded = M.load(path)
        self.assertEqual((revision, sorted(loaded)), (12, sorted((self.DATASET, learner) for learner in M.ROWS)))
        self.assertEqual(run_main(self.root, "--digests", self.reports), (0, text))

    def test_a_record_names_every_feature_the_fields_name(self):
        record = M.record(self.base, "scm")
        self.assertEqual(sorted(record), ["features", "learner", "model"])
        self.assertEqual(record["model"], {"scm": self.base["models"]["scm"]})
        # Both rules, the listed equivalent rule, and the equivalent rule the report does not list.
        self.assertEqual(record["features"], [[3, "ACGTA", 7, "AAE="], [8, "GGGTA", 2, "AAM="], [21, "AAGTC", 9, "AAU="], [34, None, None, None]])
        # Both splits, the one at the published depth too.
        self.assertEqual([f[0] for f in M.record(self.base, "cart")["features"]], [5, 13])
        self.assertEqual(M.record(self.base, "cart")["model"], {"cart": self.base["models"]["cart"], "cart_depth": 1})
        self.assertEqual(M.record(self.base, "lasso")["features"], [[3, "ACGTA", 7, "AAE="], [55, None, None, None]])
        self.assertEqual(sorted(M.record(self.base, "lasso")["model"]), ["lasso", "lasso_fit"])
        self.assertEqual([f[0] for f in M.record(self.base, "rf")["features"]], [5, 89])
        self.assertEqual(M.record(self.base, "known_gene"), {"learner": "known_gene", "model": {"rule": [["C:MADE-UP"]], "version": "made-up-1"}, "features": []})

    def test_each_field_of_a_model_changes_its_digest_and_no_other(self):
        def at(*path):
            def change(report, value):
                node = report
                for key in path[:-1]:
                    node = node[key]
                node[path[-1]] = value

            return change

        m = ("models",)
        cases = [
            ("scm", at(*m, "scm", "kind"), "conjunction", "model.scm.kind"),
            ("scm", at(*m, "scm", "p"), 0.5, "model.scm.p"),
            ("scm", at(*m, "scm", "cv_error"), 0.125, "model.scm.cv_error"),
            ("scm", at(*m, "scm", "rules", 1, "presence"), True, "model.scm.rules[1].presence"),
            ("scm", at(*m, "scm", "steps", 0, "covered"), 6, "model.scm.steps[0].covered"),
            ("scm", at(*m, "scm", "steps", 1, "n_equivalent"), 1, "model.scm.steps[1].n_equivalent"),
            ("scm", at(*m, "scm", "importance", 1), 0.25, "model.scm.importance[1]"),
            ("cart", at(*m, "cart", "nodes", 1, "value"), 0.375, "model.cart.nodes[1].value"),
            ("cart", at(*m, "cart", "nodes", 0, "split", "present"), 1, "model.cart.nodes[0].split.present"),
            ("cart", at(*m, "cart", "nodes", 0, "split", "gain"), 0.5, "model.cart.nodes[0].split.gain"),
            ("cart", at(*m, "cart", "nodes", 3, "n_resistant"), 1, "model.cart.nodes[3].n_resistant"),
            ("cart", at(*m, "cart_depth"), 2, "model.cart_depth"),
            ("lasso", at(*m, "lasso", "intercept"), -1.25, "model.lasso.intercept"),
            ("lasso", at(*m, "lasso", "coef", 1, 1), -0.5, "model.lasso.coef[1][1]"),
            # The next number a report can hold after 2.0: a weight is compared to its last digit.
            ("lasso", at(*m, "lasso", "coef", 0, 1), math.nextafter(2.0, 3.0), "model.lasso.coef[0][1]"),
            ("lasso", at(*m, "lasso", "lambda"), 0.125, "model.lasso.lambda"),
            ("lasso", at(*m, "lasso_fit", "selected_point"), 5, "model.lasso_fit.selected_point"),
            ("lasso", at(*m, "lasso_fit", "refit_convergence", "converged"), False, "model.lasso_fit.refit_convergence.converged"),
            # A report with no lasso names no feature for it either.
            ("lasso", at(*m, "lasso"), None, "features model.lasso"),
            ("rf", at(*m, "rf_importance", 1, "weight"), 0.125, "model.rf_importance[1].weight"),
            ("gbt", at(*m, "gbt_importance", 0, "weight"), 2.0, "model.gbt_importance[0].weight"),
            ("known_gene", at("catalog", "rule"), [["C:MADE-UP"], ["C:OTHER"]], "model.rule"),
            ("known_gene", at("catalog", "version"), "made-up-2", "model.version"),
            ("known_gene", at("catalog"), None, "model"),
        ]
        for learner, change, value, path in cases:
            report = copy.deepcopy(self.base)
            change(report, value)
            self.assertEqual(self.changed(report), {learner}, path)
            self.assertEqual(M.fields(self.base, report), [(learner, p) for p in path.split()], path)
        self.assertEqual(self.changed(copy.deepcopy(self.base)), set())

    def test_a_feature_id_or_what_it_stands_for_changes_each_model_that_names_it(self):
        def feature(report, i):
            (entry,) = [f for f in report["features"] if f["feature"] == i]
            return entry

        cases = [
            # (what changes, the learners that name the feature)
            (lambda r: feature(r, 3).update(kmer="ACGTC"), {"scm", "lasso"}),
            (lambda r: feature(r, 5).update(bits="AAA="), {"cart", "rf"}),
            (lambda r: feature(r, 8).update(n_kmers=3), {"scm", "gbt"}),
            # A split at the published depth, which only the tree names.
            (lambda r: feature(r, 13).update(kmer="TTGTC"), {"cart"}),
            # An equivalent rule of a step.
            (lambda r: feature(r, 21).update(kmer="AAGTA"), {"scm"}),
            # A feature the report stops listing.
            (lambda r: r["features"].remove(feature(r, 13)), {"cart"}),
            # A feature the report starts listing.
            (lambda r: r["features"].append({**feature(r, 3), "feature": 55}), {"lasso"}),
            # Another id in a field.
            (lambda r: r["models"]["gbt_importance"][0].update(feature=5), {"gbt"}),
            (lambda r: r["models"]["scm"]["steps"][0]["equivalent"][1].update(feature=35), {"scm"}),
            (lambda r: r["models"]["cart"]["nodes"][1]["split"].update(feature=21), {"cart"}),
        ]
        for change, learners in cases:
            report = copy.deepcopy(self.base)
            change(report)
            self.assertEqual(self.changed(report), learners)
            self.assertEqual({learner for learner, _ in M.fields(self.base, report)}, learners)

    def test_what_the_digest_does_not_cover_leaves_it(self):
        def others(r):
            for f in r["features"]:
                f["locus"]["product"] = "another made-up product"
                f["catalog_families"] = []
                f["n_resistant"] += 1
            # A listed feature no model names.
            r["features"][-1]["kmer"] = "CAGTA"
            r["seconds"] = 2.5
            r["splits"][0]["learners"][0]["metrics"]["n"] = 13
            r["catalog"]["metrics"]["n"] = 13
            r["features"].reverse()

        report = copy.deepcopy(self.base)
        others(report)
        self.assertEqual(report["features"][0]["feature"], 144)
        self.assertEqual(self.changed(report), set())
        self.assertEqual(M.fields(self.base, report), [])
        # Another revision changes the row's revision and its report hash, and no digest.
        later = copy.deepcopy(self.base)
        later["evaluation"]["revision"] = 14
        before, after = self.rows(self.base), self.rows(later)
        for learner in M.ROWS:
            self.assertEqual((before[learner]["evaluation_revision"], after[learner]["evaluation_revision"]), ("12", "14"))
            self.assertNotEqual(before[learner]["report_sha256"], after[learner]["report_sha256"])
            for column in ("model_sha256", "read_features_sha256"):
                self.assertEqual(before[learner][column], after[learner][column])

    def test_the_digest_is_of_the_values_and_not_of_the_file_s_layout(self):
        before = self.rows(self.base)
        raw = json.dumps(self.base, indent=2, sort_keys=True).encode()
        (self.reports / f"{self.DATASET}.json").write_bytes(raw)
        after = {r["learner"]: r for r in csv.DictReader(io.StringIO(M.digests(self.root, self.reports)), delimiter="\t")}
        for learner in M.ROWS:
            self.assertNotEqual(before[learner]["report_sha256"], after[learner]["report_sha256"])
            self.assertEqual(before[learner]["model_sha256"], after[learner]["model_sha256"])
        # A whole number written with and without a decimal point is two texts, and the digest tells them apart.
        report = copy.deepcopy(self.base)
        report["models"]["scm"]["steps"][0]["utility"] = 5
        self.assertEqual(self.changed(report), {"scm"})
        self.assertEqual(M.fields(self.base, report), [("scm", "model.scm.steps[0].utility")])
        self.assertEqual(M.sha({"b": 1, "a": [1.5, None]}), text_sha('{"a":[1.5,null],"b":1}'))
        with self.assertRaises(ValueError):
            M.sha(float("nan"))

    def test_the_features_the_kit_reads(self):
        self.assertEqual(M.read_features(self.base, "scm"), [[3, "ACGTA", 7], [8, "GGGTA", 2]])
        self.assertEqual(M.read_features(self.base, "cart"), [[5, "CCGTA", 1]])
        self.assertEqual(M.read_features(self.base, "lasso"), [[3, "ACGTA", 7], [55, None, None]])
        for learner in ("rf", "gbt", "known_gene"):
            self.assertIsNone(M.read_features(self.base, learner))
        report = copy.deepcopy(self.base)
        report["models"]["lasso"] = None
        self.assertIsNone(M.read_features(report, "lasso"))
        self.assertEqual(self.rows(report)["lasso"]["read_features_sha256"], "")
        # A lasso that keeps no weight is a model with no feature, and not a missing model.
        report["models"]["lasso"] = {"intercept": 0.5, "coef": [], "lambda": 1.0}
        self.assertEqual(M.read_features(report, "lasso"), [])
        self.assertEqual(self.rows(report)["lasso"]["read_features_sha256"], text_sha("[]"))
        # A deeper published tree reads the second split too.
        report["models"]["cart_depth"] = 2
        self.assertEqual(M.read_features(report, "cart"), [[5, "CCGTA", 1], [13, "TTGTA", 4]])

    def test_the_features_of_the_committed_run_are_read_from_its_files(self):
        recorded = M.run_features(self.root)
        with open(self.root / "primary" / "shigella_flexneri" / "features.tsv", newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(set(recorded), {(r["dataset"], r["learner"]) for r in rows})
        for key, features in recorded.items():
            want = sorted([int(r["feature"]), r["kmer"], int(r["n_kmers"])] for r in rows if (r["dataset"], r["learner"]) == key)
            self.assertEqual(features, want)
            self.assertTrue(features)

    def digest_file(self, report, name):
        self.write(report)
        path = pathlib.Path(self.tmp.name) / name
        path.write_text(M.digests(self.root, self.reports), encoding="utf-8")
        return path

    def test_a_file_is_checked_against_its_reports(self):
        path = self.digest_file(self.base, "digests.tsv")
        self.assertEqual(M.stale(path, self.reports), [])
        self.assertEqual(run_main(path, "--reports", self.reports), (0, ""))
        # The report changes under the file.
        report = copy.deepcopy(self.base)
        report["models"]["lasso"]["intercept"] = -1.25
        self.write(report)
        want = [(self.DATASET, learner, "report_sha256") for learner in M.ROWS]
        want.insert(M.ROWS.index("lasso") + 1, (self.DATASET, "lasso", "model_sha256"))
        self.assertEqual(len(want), len(M.ROWS) + 1)
        self.assertEqual(M.stale(path, self.reports), want)
        status, printed = run_main(path, "--reports", self.reports)
        self.assertEqual(status, 1)
        self.assertIn(f"differs from its report: {self.DATASET} lasso model_sha256\n", printed)
        # The file changes over the report.
        self.write(self.base)
        for column in ("model_sha256", "read_features_sha256", "evaluation_revision"):
            changed = changed_copy(path, self.tmp.name, (self.DATASET, "cart"), column, "0" * 64 if column != "evaluation_revision" else "12")
            if column == "evaluation_revision":
                self.assertEqual(M.stale(changed, self.reports), [])
                continue
            self.assertEqual(M.stale(changed, self.reports), [(self.DATASET, "cart", column)])
            self.assertEqual(run_main(changed, "--reports", self.reports)[0], 1)

    def test_two_made_up_revisions_are_compared(self):
        older = self.digest_file(self.base, "rev12.tsv")
        later_report = copy.deepcopy(self.base)
        later_report["evaluation"]["revision"] = 14
        later = self.digest_file(later_report, "rev14.tsv")
        status, printed = run_main(older, later)
        self.assertEqual(status, 0)
        self.assertEqual(printed.split("\n")[0], "revision 12 and revision 14, 1 datasets")
        self.assertIn("the five learners: 5 of 5 dataset and learner pairs\n", printed)
        # A later release that changes one tree.
        later_report["models"]["cart"]["nodes"][2]["value"] = 0.75
        later = self.digest_file(later_report, "rev14.tsv")
        status, printed = run_main(older, later)
        self.assertEqual(status, 1)
        self.assertIn("cart: the same final model on 0 of 1 datasets\n", printed)
        self.assertIn("scm: the same final model on 1 of 1 datasets\n", printed)
        self.assertIn("known_gene: the same rule on 1 of 1 datasets with a rule, and no rule in either release on 0\n", printed)
        self.assertIn("the five learners: 4 of 5 dataset and learner pairs\n", printed)
        self.assertTrue(printed.endswith(f"input fingerprint: the same on 1 of 1 datasets\ndiffers: {self.DATASET} cart\n"))
        # A later release that changes the rule and no model.
        later_report = copy.deepcopy(self.base)
        later_report["evaluation"]["revision"] = 14
        later_report["catalog"]["version"] = "made-up-2"
        status, printed = run_main(older, self.digest_file(later_report, "rev14.tsv"))
        self.assertEqual(status, 1)
        self.assertIn("known_gene: the same rule on 0 of 1 datasets with a rule, and no rule in either release on 0\n", printed)
        self.assertIn("the five learners: 5 of 5 dataset and learner pairs\n", printed)
        # No rule in either release, and then a rule in the later one only.
        bare = copy.deepcopy(self.base)
        bare["catalog"] = None
        older_bare = self.digest_file(bare, "rev12-bare.tsv")
        self.assertEqual(M.load(older_bare)[1][self.DATASET, "known_gene"]["model_sha256"], M.no_rule())
        bare["evaluation"]["revision"] = 14
        status, printed = run_main(older_bare, self.digest_file(bare, "rev14.tsv"))
        self.assertEqual(status, 0)
        self.assertIn("known_gene: the same rule on 0 of 0 datasets with a rule, and no rule in either release on 1\n", printed)
        later_report = copy.deepcopy(self.base)
        later_report["evaluation"]["revision"] = 14
        status, printed = run_main(older_bare, self.digest_file(later_report, "rev14.tsv"))
        self.assertEqual(status, 1)
        self.assertIn("known_gene: the same rule on 0 of 1 datasets with a rule, and no rule in either release on 0\n", printed)
        self.assertIn(f"\ndiffers: {self.DATASET} known_gene\n", printed)
        # A later release with other inputs and the same published fields.
        later_report = copy.deepcopy(self.base)
        later_report["evaluation"]["revision"] = 14
        later_report["input_fingerprint"] = "another-made-up-input"
        status, printed = run_main(older, self.digest_file(later_report, "rev14.tsv"))
        self.assertEqual(status, 1)
        self.assertIn("the five learners: 5 of 5 dataset and learner pairs\n", printed)
        self.assertTrue(printed.endswith(f"input fingerprint: the same on 0 of 1 datasets\ninput fingerprint differs: {self.DATASET}\n"))

    def test_the_fields_that_differ_are_printed(self):
        base = pathlib.Path(self.tmp.name)
        (base / "later").mkdir()
        self.write(self.base)
        report = copy.deepcopy(self.base)
        report["models"]["lasso"]["coef"][1][1] = -0.5
        report["models"]["rf_importance"].append({"feature": 3, "weight": 0.0625})
        self.write(report, base / "later")
        older, later = self.reports / f"{self.DATASET}.json", base / "later" / f"{self.DATASET}.json"
        # A list of another length differs as a whole, and so does its list of features.
        self.assertEqual(run_main("--fields", older, later), (0, "rf\tfeatures\nrf\tmodel.rf_importance\nlasso\tmodel.lasso.coef[1][1]\n"))
        self.assertEqual(run_main("--fields", older, older), (0, ""))
        self.write({**self.base, "id": "shigella_flexneri__ampicillin"}, base / "later")
        with self.assertRaises(SystemExit):
            M.main(["models.py", "--fields", str(older), str(later)])

    def test_two_values_are_compared_field_by_field(self):
        self.assertEqual(M.walk({"a": 1, "b": [1, 2]}, {"a": 1, "b": [1, 2]}), [])
        self.assertEqual(M.walk({"a": 1, "b": [1, 2]}, {"a": 2, "b": [1, 3]}, "x"), ["x.a", "x.b[1]"])
        self.assertEqual(M.walk({"a": 1}, {"b": 1}), [".a", ".b"])
        self.assertEqual(M.walk([1, 2], [1, 2, 3], "x"), ["x"])
        self.assertEqual(M.walk(1, 1.0, "x"), ["x"])
        self.assertEqual(M.walk(True, 1, "x"), ["x"])
        self.assertEqual(M.walk(None, {"a": 1}, "x"), ["x"])

    def test_a_report_of_another_dataset_schema_or_an_invalid_evaluation_is_refused(self):
        self.write(self.base)
        self.assertTrue(M.digests(self.root, self.reports))
        for changed in (
            {**self.base, "id": "shigella_flexneri__ampicillin"},
            {**self.base, "schema": R.score.REPORT_SCHEMA + 1},
            {**self.base, "evaluation": {"revision": 12, "valid": False}},
            {**self.base, "evaluation": {"revision": 12}},
        ):
            self.write(changed)
            with self.assertRaises(SystemExit):
                M.digests(self.root, self.reports)
            path = self.digest_file(self.base, "digests.tsv")
            self.write(changed)
            with self.assertRaises(SystemExit):
                M.stale(path, self.reports)


if __name__ == "__main__":
    unittest.main()
