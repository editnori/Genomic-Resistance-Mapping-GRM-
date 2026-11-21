"""scripts/time-split/record.py on SYNTHETIC input sets: cells worked by hand, the plan's draw order, the same bytes
on every run, the faults that stop it, and its two --check modes."""

import contextlib
import io
import json
import math
import os
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.time_split import synthetic
from tests.time_split.common import LEARNERS, SCRIPT, load_script

record = load_script()


def run(*argv):
    """main() with its output caught: (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = record.main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def plan_interval(by_species, rng):
    """The interval as the plan's analysis script writes it: 10,000 draws, indices 250 and 9,749 of the sorted draws."""
    species = sorted(by_species)
    point = statistics.median([v for s in species for v in by_species[s]])
    draws = []
    for _ in range(10000):
        values = [v for s in rng.choices(species, k=len(species)) for v in by_species[s]]
        draws.append(statistics.median(values))
    draws.sort()
    return point, draws[int(0.025 * 10000)], draws[int(0.975 * 10000) - 1]


class Synthetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="time-split-"))
        cls.inputs = cls.tmp / "inputs"
        synthetic.write(cls.inputs)
        cls.out = cls.tmp / "record"
        code, _, err = run("--inputs", cls.inputs, "--out", cls.out)
        assert code == 0, err
        cls.files = {p.name: p.read_text(encoding="utf-8") for p in cls.out.iterdir()}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def copy_inputs(self):
        dest = Path(tempfile.mkdtemp(prefix="time-split-inputs-", dir=self.tmp))
        shutil.copytree(self.inputs, dest / "in")
        return dest / "in"

    def copy_record(self):
        dest = Path(tempfile.mkdtemp(prefix="time-split-record-", dir=self.tmp))
        shutil.copytree(self.out, dest / "rec")
        return dest / "rec"

    def stops(self, inputs, words):
        code, _, err = run("--inputs", inputs, "--out", inputs.parent / "out")
        self.assertEqual(code, 1, err)
        self.assertIn(words, err)
        self.assertFalse((inputs.parent / "out").exists(), "a stopped build wrote files")

    # ------------------------------------------------------------- cells by hand

    def test_a_fit_row_by_hand(self):
        truth = [True, True, False, False, True]
        scores = [1.0, 0.0, 0.2, 1.0, 0.5]
        r = record.fit_row("x__y", "x", "past", "cart", truth, scores, 10, 4, False)
        self.assertEqual((r["tp"], r["fn"], r["tn"], r["fp"]), (2, 1, 1, 1))
        self.assertEqual(r["wrong_at_clip"], 2)  # the resistant genome at 0.0 and the susceptible one at 1.0
        self.assertEqual(r["scores_0_or_1"], 3)
        wrong = -math.log(1e-4)
        loss = -math.log(1 - 1e-4) + wrong - math.log(0.8) - math.log(1 - (1 - 1e-4)) - math.log(0.5)
        self.assertAlmostEqual(r["loss_sum"], loss, places=12)
        self.assertAlmostEqual(r["clip_loss_sum"], wrong - math.log(1 - (1 - 1e-4)), places=12)
        m = record.measures(r)
        self.assertEqual(m["recall"], 2 / 3)
        self.assertEqual(m["specificity"], 1 / 2)
        self.assertEqual(m["balanced_accuracy"], (2 / 3 + 1 / 2) / 2)
        self.assertEqual(m["log_loss"], r["loss_sum"] / 5)

    def test_exact_interval_by_hand(self):
        # Two species: draws aa, ab, ba, bb give medians 1, 2, 2, 3.
        self.assertEqual(record.exact({"a": [1.0], "b": [3.0]}), (1.0, 3.0))
        # Three species of one dataset each: 27 ordered draws; the 1st and the 27th sorted medians.
        self.assertEqual(record.exact({"a": [0.0], "b": [1.0], "c": [5.0]}), (0.0, 5.0))

    def test_readings_by_hand(self):
        self.assertEqual(record.read(0.1, 0.3, "log_loss"), "worse")
        self.assertEqual(record.read(0.1, 0.3, "recall"), "better")
        self.assertEqual(record.read(-0.3, -0.1, "recall"), "worse")
        self.assertEqual(record.read(-0.3, 0.0, "recall"), "not shown")
        self.assertEqual(record.read(0.0, 0.2, "log_loss"), "not shown")

    def test_rule_by_hand(self):
        rows = []
        for l, rc, ll in (("scm", "worse", "not shown"), ("cart", "not shown", "worse"), ("rf", "worse", "worse"),
                          ("gbt", "not shown", "not shown"), ("lasso", "better", "worse")):
            rows += [{"learner": l, "contrast": "total", "measure": "recall", "reading": rc},
                     {"learner": l, "contrast": "total", "measure": "log_loss", "reading": ll}]
        per, answer, cost = record.rule(rows, 24, 4)
        self.assertEqual([r["rule"] for r in per], ["cost", "cost", "cost", "not shown", "mixed"])
        self.assertEqual((answer, cost), ("mixed", 3))
        calm = rows[:8] + [dict(rows[8], reading="not shown"), rows[9]]  # lasso: recall not shown, log loss worse
        self.assertEqual(record.rule(calm, 24, 4)[1:], ("cost shown", 4))
        self.assertEqual(record.rule(rows, 11, 4)[1], "descriptive only")
        self.assertEqual(record.rule(rows, 24, 2)[1], "descriptive only")

    def test_draws_in_the_plans_order(self):
        # The plan's script draws every interval from one generator, learner by learner, contrast by contrast,
        # measure by measure (recall, log loss, specificity, balanced accuracy).
        rows = [r for r in self.rows("fits.tsv")]
        got = {(r["learner"], r["contrast"], r["measure"]): r for r in self.rows("readings.tsv")}
        value = {(r["dataset"], r["learner"], r["arm"]): record.measures(r) for r in rows}
        datasets = list(dict.fromkeys(r["dataset"] for r in rows))
        species = {r["dataset"]: r["species"] for r in rows}
        rng = random.Random(20261114)
        n = 0
        for l in LEARNERS:
            for name, a, b in (("total", "past", "random"), ("study", "study-out", "random"), ("time", "past", "study-out")):
                for m in ("recall", "log_loss", "specificity", "balanced_accuracy"):
                    by = {}
                    for d in datasets:
                        by.setdefault(species[d], []).append(value[(d, l, a)][m] - value[(d, l, b)][m])
                    point, lo, hi = plan_interval(by, rng)
                    r = got[(l, name, m)]
                    self.assertEqual((r["median"], r["low"], r["high"]), (repr(point), repr(lo), repr(hi)))
                    n += 1
        self.assertEqual(n, 60)

    def rows(self, name):
        return record.parse_tsv(name, self.files[name]) if name in record.ROW_FILES else [
            dict(zip(self.files[name].split("\n")[0].split("\t"), line.split("\t")))
            for line in self.files[name].rstrip("\n").split("\n")[1:]]

    def test_synthetic_counts(self):
        fits = self.rows("fits.tsv")
        self.assertEqual(len(fits), 12 * 3 * 5)
        summary = dict((r["name"], r["value"]) for r in self.rows("summary.tsv"))
        self.assertEqual((summary["datasets"], summary["species"], summary["tested_pairs"]), ("12", "3", "96"))
        self.assertEqual((summary["ordered_species_draws"], summary["distinct_species_draws"]), ("27", "10"))
        checks = dict((r["check"], r["count"]) for r in self.rows("checks.tsv"))
        self.assertEqual(checks["genome_ids_of_the_inputs_in_the_record"], "0")
        self.assertEqual(checks["fits_match_evaluator_counts"], "180")
        review = self.rows("review.tsv")
        self.assertEqual([r["fits"] for r in review], [36] * 5)
        self.assertTrue(all(r["scores_changed"] > 0 for r in review), "the synthetic second run differs from the first")

    # ------------------------------------------------------------- the same bytes

    def test_same_bytes_on_a_second_run(self):
        code, out, err = run("--inputs", self.inputs, "--check", self.out)
        self.assertEqual(code, 0, err)
        self.assertIn("18 files rebuilt from the input set, byte for byte", out)

    def test_same_bytes_at_other_hash_seeds(self):
        # Without -I, which would drop PYTHONHASHSEED; -s and -P keep site-packages and the cwd out.
        seen = set()
        for seed in ("0", "4242"):
            out = self.tmp / f"seed-{seed}"
            env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "PYTHONHASHSEED": seed}
            subprocess.run([sys.executable, "-s", "-P", "-B", str(SCRIPT), "--inputs", str(self.inputs), "--out",
                            str(out)], check=True, env=env, capture_output=True)
            self.assertEqual({p.name: p.read_text(encoding="utf-8") for p in out.iterdir()}, self.files, seed)
            seen.add(subprocess.run([sys.executable, "-s", "-P", "-c", "import sys; print(sys.flags.ignore_environment, "
                                     "hash('time split'))"], env=env, check=True, capture_output=True, text=True).stdout)
        self.assertEqual(len(seen), 2, "the two seeds gave one hash")
        self.assertTrue(all(x.startswith("0 ") for x in seen), seen)

    def test_check_alone(self):
        code, out, err = run("--check", self.out)
        self.assertEqual(code, 0, err)
        self.assertIn("18 files agree with the row files and SHA256SUMS", out)

    def test_out_writes_only_into_its_own_record(self):
        other = self.tmp / "other"
        other.mkdir(exist_ok=True)
        (other / "notes.txt").write_text("not the record's\n")
        code, _, err = run("--inputs", self.inputs, "--out", other)
        self.assertEqual(code, 1)
        self.assertIn("holds files this script did not write", err)

    def test_out_rewrites_an_earlier_record(self):
        rec = self.copy_record()
        code, _, err = run("--inputs", self.inputs, "--out", rec)
        self.assertEqual(code, 0, err)
        self.assertEqual({p.name: p.read_text(encoding="utf-8") for p in rec.iterdir()}, self.files)

    def test_out_refuses_other_files_of_the_records_names(self):
        rec = self.copy_record()
        for p in rec.iterdir():
            p.write_text(f"not the record's {p.name}\n", encoding="utf-8")
        before = {p.name: p.read_bytes() for p in rec.iterdir()}
        code, _, err = run("--inputs", self.inputs, "--out", rec)
        self.assertEqual(code, 1, err)
        self.assertIn("SHA256SUMS", err)
        self.assertEqual({p.name: p.read_bytes() for p in rec.iterdir()}, before)

    def test_out_refuses_a_file_changed_since_its_sums(self):
        rec = self.copy_record()
        (rec / "mix.tsv").write_text("arm\nnot the record's\n", encoding="utf-8")
        code, _, err = run("--inputs", self.inputs, "--out", rec)
        self.assertEqual(code, 1, err)
        self.assertIn("changed since: mix.tsv", err)
        self.assertEqual((rec / "mix.tsv").read_text(encoding="utf-8"), "arm\nnot the record's\n")

    def test_out_refuses_a_symlink_in_the_record(self):
        rec = self.copy_record()
        target = rec.parent / "elsewhere.tsv"
        target.write_text(self.files["mix.tsv"], encoding="utf-8")
        (rec / "mix.tsv").unlink()
        (rec / "mix.tsv").symlink_to(target)
        code, _, err = run("--inputs", self.inputs, "--out", rec)
        self.assertEqual(code, 1, err)
        self.assertIn("mix.tsv: a symlink or not a plain file", err)
        self.assertTrue((rec / "mix.tsv").is_symlink())
        self.assertEqual(target.read_text(encoding="utf-8"), self.files["mix.tsv"])

    def test_out_refuses_a_symlink_as_the_record(self):
        rec = self.copy_record()
        link = rec.parent / "link"
        link.symlink_to(rec, target_is_directory=True)
        code, _, err = run("--inputs", self.inputs, "--out", link)
        self.assertEqual(code, 1, err)
        self.assertIn("a symlink", err)

    # ------------------------------------------------------------- faults that stop the build

    def test_fault_fold_file_hash(self):
        i = self.copy_inputs()
        p = i / "results" / "alpha_species__drug_one.first-public-random.json"
        body = json.loads(p.read_text())
        body["fold_file_sha256"] = "0" * 64
        p.write_text(json.dumps(body))
        self.stops(i, "not the design's fold file")

    def test_fault_arms_test_other_genomes(self):
        i = self.copy_inputs()
        p = i / "results" / "beta_species__drug_two.first-public-study-out.json"
        body = json.loads(p.read_text())
        body["samples"][0]["genome_id"] = "991111.0001"
        p.write_text(json.dumps(body))
        self.stops(i, "the arms test other genomes")

    def test_fault_evaluator_counts(self):
        i = self.copy_inputs()
        p = i / "results" / "gamma_species__drug_three.first-public-past.json"
        body = json.loads(p.read_text())
        body["split"]["learners"][2]["metrics"]["tp"] += 1
        p.write_text(json.dumps(body))
        self.stops(i, "confusion counts other than the evaluator's")

    def test_fault_a_missing_result_file(self):
        i = self.copy_inputs()
        (i / "before" / "alpha_species__drug_four.first-public-past.json").unlink()
        self.stops(i, "not in the input set")

    def test_fault_training_counts(self):
        i = self.copy_inputs()
        p = i / "results" / "alpha_species__drug_two.first-public-past.json"
        body = json.loads(p.read_text())
        body["folds"][0]["train"] += 1
        p.write_text(json.dumps(body))
        self.stops(i, "training counts")

    def test_fault_plan_hash(self):
        i = self.copy_inputs()
        (i / "plan" / "PLAN.md").write_text("# Synthetic plan, edited after its hash\n")
        self.stops(i, "not the hashed plan")

    def test_fault_analysis_script_hash(self):
        i = self.copy_inputs()
        (i / "plan" / "scripts" / "analyze.py").write_text("# edited\n")
        self.stops(i, "not the hashed analysis script")

    def test_fault_plan_after_first_result(self):
        i = self.copy_inputs()
        p = i / "rerun" / "RERUN-PLAN.sha256"
        p.write_text(p.read_text().replace("2026-01-02T09:00:00Z", "2026-01-02T09:30:00Z"))
        self.stops(i, "a plan was hashed after its first result")

    def test_fault_a_rerun_that_failed(self):
        i = self.copy_inputs()
        p = i / "logs" / "rerun.log"
        lines = p.read_text().split("\n")
        lines[3] = lines[3].replace(" rc=0 ", " rc=1 ")
        p.write_text("\n".join(lines))
        self.stops(i, "a rerun that did not exit 0")

    def test_fault_studies_disagree_with_the_design(self):
        i = self.copy_inputs()
        p = i / "studies.json"
        body = json.loads(p.read_text())
        body["per_dataset"]["beta_species__drug_one"]["tested_genomes"] += 1
        p.write_text(json.dumps(body))
        self.stops(i, "drawn studies or tested genomes other than the design's")

    def test_fault_studies_in_two_species_do_not_add_up(self):
        i = self.copy_inputs()
        p = i / "studies.json"
        body = json.loads(p.read_text())
        body["distinct_test_studies"] += 1
        p.write_text(json.dumps(body))
        self.stops(i, "studies in two species do not explain the sum")

    def test_fault_a_genome_id_in_a_cell(self):
        i = self.copy_inputs()
        p = i / "design.json"
        body = json.loads(p.read_text())
        body["datasets"][0]["antibiotic"] = "990000.1003"
        p.write_text(json.dumps(body))
        self.stops(i, "genome ids of the input set stand in the record")

    def test_fault_other_settings_in_the_two_runs(self):
        i = self.copy_inputs()
        p = i / "before" / "gamma_species__drug_one.first-public-random.json"
        body = json.loads(p.read_text())
        body["eval_config"]["screened_features"] = 2000
        p.write_text(json.dumps(body))
        self.stops(i, "the result files hold other settings")

    def test_fault_a_score_outside_0_to_1(self):
        i = self.copy_inputs()
        p = i / "results" / "beta_species__drug_three.first-public-random.json"
        body = json.loads(p.read_text())
        body["split"]["learners"][4]["scores"][0] = 1.5
        p.write_text(json.dumps(body))
        self.stops(i, "a score outside 0 to 1")

    # ------------------------------------------------------------- --check

    def test_check_with_inputs_finds_a_changed_score(self):
        i = self.copy_inputs()
        p = i / "results" / "alpha_species__drug_three.first-public-past.json"
        body = json.loads(p.read_text())
        learner = body["split"]["learners"][3]
        learner["scores"][0] = 0.25 if learner["scores"][0] >= 0.5 else 0.75
        t = [s["resistant"] for s in body["samples"]]
        sc = learner["scores"]
        learner["metrics"].update({
            "tp": sum(1 for a, s in zip(t, sc) if a and s >= 0.5), "fn": sum(1 for a, s in zip(t, sc) if a and s < 0.5),
            "tn": sum(1 for a, s in zip(t, sc) if not a and s < 0.5), "fp": sum(1 for a, s in zip(t, sc) if not a and s >= 0.5)})
        m = learner["metrics"]
        m["sensitivity"], m["specificity"] = m["tp"] / (m["tp"] + m["fn"]), m["tn"] / (m["tn"] + m["fp"])
        p.write_text(json.dumps(body))
        code, _, err = run("--inputs", i, "--check", self.out)
        self.assertEqual(code, 1)
        self.assertIn("differs from the input set's record", err)
        self.assertIn("fits.tsv", err)

    def test_check_alone_finds_a_changed_derived_cell(self):
        r = self.copy_record()
        p = r / "readings.tsv"
        p.write_text(p.read_text().replace("\tworse\t", "\tnot shown\t", 1) if "\tworse\t" in p.read_text()
                     else p.read_text().replace("\tnot shown\t", "\tworse\t", 1))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("readings.tsv", err)

    def test_check_alone_finds_a_changed_row_cell(self):
        r = self.copy_record()
        p = r / "fits.tsv"
        lines = p.read_text().split("\n")
        cells = lines[5].split("\t")
        cells[14] = repr(float(cells[14]) + 0.5)  # loss_sum
        lines[5] = "\t".join(cells)
        p.write_text("\n".join(lines))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("contrasts.tsv", err)

    def test_check_alone_finds_a_changed_count(self):
        r = self.copy_record()
        p = r / "fits.tsv"
        lines = p.read_text().split("\n")
        cells = lines[7].split("\t")
        cells[8] = str(int(cells[8]) + 1)  # tp
        lines[7] = "\t".join(cells)
        p.write_text("\n".join(lines))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("the confusion counts do not add up", err)

    def test_check_alone_finds_a_dropped_file(self):
        r = self.copy_record()
        (r / "mix.tsv").unlink()
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("mix.tsv", err)

    def test_check_alone_finds_an_added_file(self):
        r = self.copy_record()
        (r / "genomes.tsv").write_text("genome\n")
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("genomes.tsv", err)

    def test_check_alone_finds_a_changed_unbuilt_file(self):
        r = self.copy_record()
        p = r / "inputs.sha256"
        p.write_text(p.read_text().replace("0", "1", 1))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("SHA256SUMS", err)

    def test_check_alone_finds_an_added_row(self):
        r = self.copy_record()
        p = r / "fits.tsv"
        lines = p.read_text().split("\n")
        lines.insert(2, lines[1])
        p.write_text("\n".join(lines))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("not one row per dataset, arm and learner", err)

    def test_check_alone_finds_a_genome_column(self):
        r = self.copy_record()
        p = r / "species.tsv"
        lines = p.read_text().split("\n")
        lines[0] += "\tgenome"
        lines[1:-1] = [x + "\t990000.1000" for x in lines[1:-1]]
        p.write_text("\n".join(lines))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("species.tsv: the columns are not", err)

    def test_check_alone_finds_a_long_float(self):
        r = self.copy_record()
        p = r / "fits.tsv"
        lines = p.read_text().split("\n")
        cells = lines[3].split("\t")
        cells[15] = cells[15] + "0"  # clip_loss_sum, no longer in its shortest form
        lines[3] = "\t".join(cells)
        p.write_text("\n".join(lines))
        code, _, err = run("--check", r)
        self.assertEqual(code, 1)
        self.assertIn("not in its shortest form", err)


if __name__ == "__main__":
    unittest.main()
