"""The committed record under reference/time-split: its files and sums, its columns and rows, its size, what it
must never hold, and the script's --check mode on it."""

import contextlib
import hashlib
import io
import re
import unittest

from tests.time_split.common import ARMS, LEARNERS, RECORD, load_script, pairs, read_tsv

record = load_script()
MAX_BYTES = 3 * 1024 * 1024
# The sha256 of SHA256SUMS as committed. Any change to a file of the record changes it.
PINNED = "093037620c3e4ce00c54f40a7cabb35d24f82fcdb1ec9b54cec197922a26613c"
DATASETS, SPECIES, CONTRASTS, MEASURES = 24, 4, 3, 4
ROWS = {"fits.tsv": DATASETS * len(ARMS) * len(LEARNERS), "datasets.tsv": DATASETS, "species.tsv": SPECIES,
        "review.tsv": len(LEARNERS), "order.tsv": 7, "contrasts.tsv": DATASETS * len(LEARNERS) * CONTRASTS,
        "readings.tsv": len(LEARNERS) * CONTRASTS * MEASURES,
        "species-medians.tsv": len(LEARNERS) * CONTRASTS * MEASURES * SPECIES, "rule.tsv": len(LEARNERS),
        "arms.tsv": len(LEARNERS) * len(ARMS), "mix.tsv": len(ARMS)}
# Column names that would mean a row per genome, a study named, or the content of a fit.
BANNED_COLUMN = re.compile(r"^(genome|sample|isolate|accession|assembly|bioproject|biosample|study|label|score|kmer|"
                           r"feature|marker|gene|mutation|variant|weight|importance|sequence)(_ids?|_names?|_accession)?$"
                           r"|(^|_)ids?$")
# Paths of inputs that may never feed the record.
BANNED_PATH = re.compile(r"(?i)seal[e]d|held-out|ks[0-9]{2}|vum[c]|clini[c]al|priv[a]te|ks[d]|sigdr[u]g")
INPUT_PATH = re.compile(r"(design\.json|studies\.json|(results|before)/[a-z0-9_]+__[a-z0-9_]+\.first-public-"
                        r"(past|random|study-out)\.json|plan/(PLAN\.md|PLAN\.sha256|tables\.sha256|scripts/(analyze|tables)\.py)|"
                        r"rerun/RERUN-PLAN\.(md|sha256)|logs/(first-run|rerun)\.log)")


def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = record.main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


class Files(unittest.TestCase):
    def test_file_list(self):
        self.assertEqual(sorted(p.name for p in RECORD.iterdir()), sorted(record.FILES))
        self.assertEqual(len(record.FILES), int(pairs(RECORD / "summary.tsv")["record_files"]))

    def test_sums(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text(encoding="utf-8").rstrip("\n").split("\n"):
            digest, name = line.split("  ")
            listed[name] = digest
        self.assertEqual(sorted(listed), sorted(n for n in record.FILES if n != "SHA256SUMS"))
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)

    def test_pinned(self):
        self.assertEqual(hashlib.sha256((RECORD / "SHA256SUMS").read_bytes()).hexdigest(), PINNED)

    def test_check_alone(self):
        code, out, err = run("--check", RECORD)
        self.assertEqual(code, 0, err)
        self.assertEqual(out, f"{len(record.FILES)} files agree with the row files and SHA256SUMS\n")

    def test_size(self):
        self.assertLess(sum(p.stat().st_size for p in RECORD.iterdir()), MAX_BYTES)

    def test_text_files_end_with_one_newline(self):
        for p in RECORD.iterdir():
            text = p.read_text(encoding="utf-8")
            self.assertTrue(text.endswith("\n") and not text.endswith("\n\n"), p.name)
            self.assertNotIn("\r", text, p.name)


class Tables(unittest.TestCase):
    def test_columns(self):
        for name, cols in record.COLUMNS.items():
            head = (RECORD / name).read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
            self.assertEqual(head, cols, name)

    def test_rows(self):
        for name, count in ROWS.items():
            self.assertEqual(len(read_tsv(RECORD / name)), count, name)
        self.assertEqual(len((RECORD / "inputs.sha256").read_text().rstrip("\n").split("\n")),
                         int(pairs(RECORD / "totals.tsv")["result_files"]) * 2 + 11)

    def test_one_row_per_dataset_arm_and_learner(self):
        datasets = [r["dataset"] for r in read_tsv(RECORD / "datasets.tsv")]
        fits = [(r["dataset"], r["arm"], r["learner"]) for r in read_tsv(RECORD / "fits.tsv")]
        self.assertEqual(fits, [(d, a, l) for d in datasets for a in ARMS for l in LEARNERS])
        for d in read_tsv(RECORD / "datasets.tsv"):
            self.assertRegex(d["dataset"], r"^[a-z]+_[a-z]+__[a-z0-9_]+$")
            self.assertEqual(d["dataset"], d["species"] + "__" + re.sub(r"[^a-z0-9]+", "_", d["antibiotic"].lower()))

    def test_checks(self):
        checks = {r["check"]: int(r["count"]) for r in read_tsv(RECORD / "checks.tsv")}
        files = DATASETS * len(ARMS)
        self.assertEqual(checks, {
            "result_files": files, "fold_file_hash_matches_design": files, "arms_test_the_same_genomes": files,
            "training_counts_match_design": files, "fits_match_evaluator_counts": files * len(LEARNERS),
            "before_result_files": files, "plans_match_their_hash": 2, "scripts_match_their_hash": 2,
            "rerun_runs_exit_0": files, "genome_ids_of_the_inputs_in_the_record": 0})

    def test_settings(self):
        s = pairs(RECORD / "settings.tsv")
        for k, v in record.FIXED_SETTINGS.items():
            self.assertEqual(s[k], v, k)
        self.assertEqual((s["seed"], s["screened_features"], s["inner_folds"]), ("42", "4000", "3"))
        self.assertEqual((s["min_test_class"], s["min_train_class"], s["min_reserve"], s["min_later_other_studies"]),
                         ("20", "25", "20", "20"))
        for k in ("plan_sha256", "analysis_sha256", "layout_sha256", "rerun_plan_sha256"):
            self.assertRegex(s[k], r"^[0-9a-f]{64}$")

    def test_cells_follow_from_each_other(self):
        fits = read_tsv(RECORD / "fits.tsv")
        for r in fits:
            c = {k: int(r[k]) for k in ("test", "test_resistant", "tp", "fn", "tn", "fp", "wrong_at_clip")}
            self.assertEqual(c["tp"] + c["fn"], c["test_resistant"])
            self.assertEqual(c["tn"] + c["fp"], c["test"] - c["test_resistant"])
            self.assertLessEqual(c["wrong_at_clip"], c["fn"] + c["fp"])
            if r["learner"] == "scm":
                self.assertEqual(int(r["scores_0_or_1"]), c["test"])
                self.assertEqual(c["wrong_at_clip"], c["fn"] + c["fp"])
        arms = {(r["learner"], r["arm"]): r for r in read_tsv(RECORD / "arms.tsv")}
        for (learner, arm), r in arms.items():
            mine = [f for f in fits if f["learner"] == learner and f["arm"] == arm]
            self.assertEqual(int(r["wrong_calls"]), sum(int(f["fn"]) + int(f["fp"]) for f in mine))
            self.assertEqual(int(r["tested_pairs"]), sum(int(f["test"]) for f in mine))


class NeverHeld(unittest.TestCase):
    def test_no_column_per_genome_or_study_or_feature(self):
        for name in record.COLUMNS:
            for c in (RECORD / name).read_text(encoding="utf-8").split("\n", 1)[0].split("\t"):
                self.assertIsNone(BANNED_COLUMN.search(c), (name, c))

    def test_no_text_of_a_banned_shape(self):
        for p in sorted(RECORD.iterdir()):
            text = p.read_text(encoding="utf-8")
            for shape in record.BANNED_SHAPES:
                if p.name == "inputs.sha256" and shape == r"[ACGT]{20,}":
                    continue
                self.assertIsNone(re.search(shape, text), (p.name, shape))

    def test_no_row_per_genome(self):
        pairs_ = int(pairs(RECORD / "summary.tsv")["tested_pairs"])
        genomes = int(pairs(RECORD / "totals.tsv")["distinct_test_genomes"])
        for p in sorted(RECORD.glob("*.tsv")):
            rows = len(read_tsv(p))
            self.assertLess(rows, genomes, p.name)
            self.assertLess(rows, pairs_, p.name)

    def test_inputs_are_the_named_files(self):
        lines = (RECORD / "inputs.sha256").read_text(encoding="utf-8").rstrip("\n").split("\n")
        for line in lines:
            digest, path = line.split("  ")
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
            self.assertRegex(path, INPUT_PATH)
            self.assertIsNone(BANNED_PATH.search(path), path)


if __name__ == "__main__":
    unittest.main()
