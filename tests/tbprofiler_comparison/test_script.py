"""scripts/tbprofiler-comparison/record.py on a SYNTHETIC input set: cells worked by hand, the lineage interval
redrawn by a second, plain reading of the evaluator's choices, the same bytes on every run, the faults that stop
it, and its two --check modes."""

import contextlib
import gzip
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.tbprofiler_comparison import synthetic
from tests.tbprofiler_comparison.common import SCRIPT, load_script

record = load_script()
FILES = len(record.FILES)


def run(*argv):
    """main() with its output caught: (exit code, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = record.main([str(a) for a in argv])
    return code, out.getvalue(), err.getvalue()


def plain_groups(reports):
    """Genome to group number, as the evaluator numbers them, by a plain search over shared lineage ids."""
    links = {}
    for d in reports:
        by = {}
        for g, lin in zip(d["sid"], d["lineage"]):
            by.setdefault(lin, []).append(g)
        for gs in by.values():
            for g in gs:
                links.setdefault(g, set()).update(gs)
    seen, comps = set(), []
    for g in links:
        if g in seen:
            continue
        stack, comp = [g], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(links[x] - comp)
        seen |= comp
        comps.append(comp)
    comps.sort(key=min)
    return {g: i for i, c in enumerate(comps) for g in c}, len(comps)


class Synthetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="tbprofiler-comparison-"))
        cls.inputs = cls.tmp / "inputs"
        cls.counts = synthetic.write(cls.inputs)
        cls.out = cls.tmp / "record"
        code, _, err = run("--inputs", cls.inputs, "--out", cls.out)
        assert code == 0, err
        cls.files = {p.name: p.read_text(encoding="utf-8") for p in cls.out.iterdir()}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def copy_inputs(self):
        dest = Path(tempfile.mkdtemp(prefix="inputs-", dir=self.tmp))
        shutil.copytree(self.inputs, dest / "in")
        return dest / "in"

    def copy_record(self):
        dest = Path(tempfile.mkdtemp(prefix="record-", dir=self.tmp))
        shutil.copytree(self.out, dest / "rec")
        return dest / "rec"

    def rows(self, name):
        lines = self.files[name].rstrip("\n").split("\n")
        return [dict(zip(lines[0].split("\t"), line.split("\t"))) for line in lines[1:]]

    def stops(self, inputs, words):
        code, _, err = run("--inputs", inputs, "--out", inputs.parent / "out")
        self.assertEqual(code, 1, err)
        self.assertIn(words, err)
        self.assertFalse((inputs.parent / "out").exists(), "a stopped build wrote files")

    def edit_compact(self, inputs, change):
        p = inputs / "rev14_compact.jsonl.gz"
        with gzip.open(p, "rt", encoding="utf-8") as f:
            reports = [json.loads(line) for line in f]
        change(reports)
        with gzip.open(p, "wt", encoding="utf-8") as f:
            for d in reports:
                f.write(json.dumps(d) + "\n")

    def edit_calls(self, inputs, change):
        p = inputs / "calls.tsv"
        lines = p.read_text(encoding="utf-8").split("\n")
        head = lines[0].split("\t")
        rows = [dict(zip(head, x.split("\t"))) for x in lines[1:-1]]
        change(rows)
        p.write_text("\n".join([lines[0]] + ["\t".join(r[c] for c in head) for r in rows]) + "\n", encoding="utf-8")

    # ------------------------------------------------------------- cells by hand

    def test_rates_by_hand(self):
        self.assertEqual(record.rates(3, 1, 2, 1), (0.75, 2 / 3, (0.75 + 2 / 3) / 2))
        self.assertEqual(record.rates(0, 0, 5, 0), (0.0, 1.0, 0.5))
        self.assertEqual(record.rates(0, 0, 5, 0, strict=True), (None, None, None))
        self.assertEqual(record.confusion([True, True, False, False], [True, False, True, False]), (1, 1, 1, 1))

    def test_interval_and_median_by_hand(self):
        self.assertEqual(record.interval(list(range(100))), (2, 97, 100))
        self.assertEqual(record.interval([None] * 10 + list(range(39))), (None, None, 39))
        self.assertEqual(record.interval(list(range(2000))), (50, 1949, 2000))  # the 51st and the 1,950th
        self.assertEqual(record.median([3.0, 1.0, 2.0]), 2.0)
        self.assertEqual(record.median([4.0, 1.0, 2.0, 3.0]), 2.5)

    def test_groups_by_hand(self):
        reports = [{"sid": ["b", "a", "c", "d"], "lineage": [1, 1, 2, 3]}, {"sid": ["d", "c", "e"], "lineage": [5, 5, 6]}]
        group, n = record.lineage_groups(reports)
        self.assertEqual((group, n), ({"a": 0, "b": 0, "c": 1, "d": 1, "e": 2}, 3))

    def test_synthetic_counts(self):
        totals = dict((r["name"], r["value"]) for r in self.rows("totals.tsv"))
        self.assertEqual((totals["reports_of_the_species"], totals["primary_reports"], totals["proxy_reports"],
                          totals["not_named_reports"], totals["lineage_groups"], totals["distinct_genomes"]),
                         ("6", "4", "1", "1", "19", "160"))
        summary = dict((r["name"], r["value"]) for r in self.rows("summary.tsv"))
        self.assertEqual((summary["by_name_endpoints"], summary["graded_endpoints"], summary["ungraded_in_by_name"]),
                         ("4", "3", "1"))
        checks = dict((r["check"], r["count"]) for r in self.rows("checks.tsv"))
        self.assertEqual(checks["genome_ids_of_the_inputs_in_the_record"], "0")
        tool = [r for r in self.rows("callers.tsv")
                if (r["way"], r["dataset"], r["caller"]) == ("who_assoc", "mycobacterium_tuberculosis__betamycin", "tool")]
        self.assertEqual([int(tool[0][k]) for k in ("tp", "fp", "tn", "fn")], self.counts)
        ungraded = [r for r in self.rows("callers.tsv") if r["dataset"].endswith("__deltazid") and r["caller"] == "tool"
                    and r["way"] == "who_assoc"]
        self.assertEqual((ungraded[0]["tp"], ungraded[0]["fp"]), ("0", "0"))

    def test_proxies_enter_no_set(self):
        for r in self.rows("endpoints.tsv"):
            if r["kind"] == "proxy":
                self.assertEqual((r["by_name"], r["graded"]), ("no", "no"))
        med = [r for r in self.rows("medians.tsv") if r["set"] == "by_name"]
        self.assertEqual({r["endpoints"] for r in med}, {"4"})

    # ------------------------------------------------------------- the interval, read a second way

    def test_interval_redrawn_plainly(self):
        with gzip.open(self.inputs / "rev14_compact.jsonl.gz", "rt", encoding="utf-8") as f:
            reports = [d for d in map(json.loads, f) if d["species_slug"] == "mycobacterium_tuberculosis"]
        group, n = plain_groups(reports)
        calls = {}
        for line in (self.inputs / "calls.tsv").read_text(encoding="utf-8").split("\n")[1:-1]:
            cells = line.split("\t")
            calls[(cells[0], cells[1])] = cells[4]
        by_name = [r["dataset"] for r in self.rows("endpoints.tsv") if r["way"] == "who_assoc" and r["by_name"] == "yes"]
        per = {}
        for d in reports:
            if d["id"] not in by_name:
                continue
            counts = {}
            for g, y in zip(d["sid"], d["y"]):
                c = calls[(d["id"], g)]
                if c == "no_call":
                    continue
                tp, fn, s = counts.get(group[g], (0, 0, 0))
                counts[group[g]] = (tp + (y and c == "resistant"), fn + (y and c == "susceptible"), s + (not y))
            per[d["id"]] = counts
        rng = random.Random(20261114)
        meds, first = [], []
        for _ in range(2000):
            draws = [rng.randrange(n) for _ in range(n)]
            sens = []
            for d in by_name:
                tp, fn, s = (sum(per[d].get(c, (0, 0, 0))[j] for c in draws) for j in range(3))
                value = tp / (tp + fn) if tp + fn and s else None  # no value without both classes
                if value is not None:
                    sens.append(value)
                if d == by_name[0]:
                    first.append(value)
            sens.sort()
            k = len(sens) // 2
            meds.append(sens[k] if len(sens) % 2 else (sens[k - 1] + sens[k]) / 2)
        meds.sort()
        first = sorted(x for x in first if x is not None)
        iv = {(r["kind"], r["set"], r["dataset"], r["caller"], r["measure"]): r for r in self.rows("intervals.tsv")}
        got = iv[("median", "by_name", "-", "tool", "sensitivity")]
        self.assertEqual((got["low"], got["high"]), (repr(meds[50]), repr(meds[1949])))
        got = iv[("endpoint", "-", by_name[0], "tool", "sensitivity")]
        k = int(0.025 * len(first))
        self.assertEqual((got["low"], got["high"]), (repr(first[k]), repr(first[len(first) - 1 - k])))

    # ------------------------------------------------------------- the same bytes

    def test_same_bytes_on_a_second_run(self):
        code, out, err = run("--inputs", self.inputs, "--check", self.out)
        self.assertEqual(code, 0, err)
        self.assertIn(f"{FILES} files rebuilt from the input set, byte for byte", out)

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
                                     "hash('tbprofiler comparison'))"], env=env, check=True, capture_output=True,
                                    text=True).stdout)
        self.assertEqual(len(seen), 2, "the two seeds gave one hash")
        self.assertTrue(all(x.startswith("0 ") for x in seen), seen)

    def test_check_alone(self):
        code, out, err = run("--check", self.out)
        self.assertEqual(code, 0, err)
        self.assertIn(f"{FILES} files agree with the row files and SHA256SUMS", out)

    # ------------------------------------------------------------- --out

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
        (rec / "screen.tsv").write_text("way\nnot the record's\n", encoding="utf-8")
        code, _, err = run("--inputs", self.inputs, "--out", rec)
        self.assertEqual(code, 1, err)
        self.assertIn("changed since: screen.tsv", err)
        self.assertEqual((rec / "screen.tsv").read_text(encoding="utf-8"), "way\nnot the record's\n")

    def test_out_refuses_a_symlink_in_the_record(self):
        rec = self.copy_record()
        target = rec.parent / "elsewhere.tsv"
        target.write_text(self.files["screen.tsv"], encoding="utf-8")
        (rec / "screen.tsv").unlink()
        (rec / "screen.tsv").symlink_to(target)
        code, _, err = run("--inputs", self.inputs, "--out", rec)
        self.assertEqual(code, 1, err)
        self.assertIn("screen.tsv: a symlink or not a plain file", err)
        self.assertTrue((rec / "screen.tsv").is_symlink())
        self.assertEqual(target.read_text(encoding="utf-8"), self.files["screen.tsv"])

    def test_out_refuses_a_symlink_as_the_record(self):
        rec = self.copy_record()
        link = rec.parent / "link"
        link.symlink_to(rec, target_is_directory=True)
        code, _, err = run("--inputs", self.inputs, "--out", link)
        self.assertEqual(code, 1, err)
        self.assertIn("a symlink", err)

    def test_out_needs_inputs(self):
        code, _, err = run("--out", self.tmp / "nowhere")
        self.assertEqual(code, 1)
        self.assertIn("--out needs --inputs", err)

    # ------------------------------------------------------------- faults that stop the build

    def test_fault_plan_hash(self):
        i = self.copy_inputs()
        (i / "PLAN.md").write_text("# SYNTHETIC plan, edited after its hash\n")
        self.stops(i, "not the hashed plan")

    def test_fault_addendum_after_scoring(self):
        i = self.copy_inputs()
        p = i / "PLAN.sha256"
        p.write_text(p.read_text().replace("2026-01-01T10:05:00Z", "2026-01-01T11:30:00Z"))
        self.stops(i, "a plan was hashed after scoring began")

    def test_fault_a_missing_log_step(self):
        i = self.copy_inputs()
        p = i / "finish.log"
        p.write_text(p.read_text().replace("evaluator done", "evaluator stopped"))
        self.stops(i, "finish.log: no line 'evaluator done'")

    def test_fault_another_tool_version(self):
        i = self.copy_inputs()
        plan = b"# SYNTHETIC plan\n\nTool: TBProfiler 6.6.4 with its database, commit `3ed1cb99`.\n"
        (i / "PLAN.md").write_bytes(plan)
        p = i / "PLAN.sha256"
        p.write_text(p.read_text().replace(p.read_text()[:64], synthetic.sha(plan), 1))
        self.stops(i, "the plan names another tool version")

    def test_fault_a_missing_input(self):
        i = self.copy_inputs()
        (i / "exploratory.json").unlink()
        self.stops(i, "exploratory.json: not in the input set")

    def test_fault_call_file_genomes(self):
        i = self.copy_inputs()
        self.edit_calls(i, lambda rows: rows[3].update(genome_id="99006.9999"))
        self.stops(i, "the call file's genomes are not the report's")

    def test_fault_a_call_of_no_kind(self):
        i = self.copy_inputs()
        self.edit_calls(i, lambda rows: rows[5].update(tool_default="maybe"))
        self.stops(i, "a call that is not resistant, susceptible or no_call")

    def test_fault_a_resistant_call_without_a_counted_grade(self):
        i = self.copy_inputs()

        def change(rows):
            r = next(r for r in rows if r["who_assoc"] == "resistant")
            r["top_grade"] = "-"
        self.edit_calls(i, change)
        self.stops(i, "a resistant call with no counted grade")

    def test_fault_a_score_outside_0_to_1(self):
        i = self.copy_inputs()
        self.edit_compact(i, lambda reports: reports[1]["splits"]["lineage"]["scores"]["rf"].__setitem__(0, 1.5))
        self.stops(i, "the rf lineage scores")

    def test_fault_counted_entries(self):
        i = self.copy_inputs()
        p = i / "database-grades-by-drug.tsv"
        lines = p.read_text().split("\n")
        cells = lines[1].split("\t")
        cells[-1] = str(int(cells[-1]) + 1)
        lines[1] = "\t".join(cells)
        p.write_text("\n".join(lines))
        self.stops(i, "the counted entries are not the two counted grades")

    def test_fault_set_aside_counts(self):
        i = self.copy_inputs()
        p = i / "exploratory.json"
        body = json.loads(p.read_text())
        touched = body["rule_that_sets_an_entry_aside"]["primary_endpoints_touched"]
        next(iter(touched.values()))["counts_as_run"][0] += 1
        p.write_text(json.dumps(body))
        self.stops(i, "set-aside.tsv: the counts as run are not callers.tsv's")

    def test_fault_a_genome_id_in_a_cell(self):
        i = self.copy_inputs()
        self.edit_compact(i, lambda reports: reports[5].update(antibiotic=synthetic.genome_id(0)))
        self.stops(i, "genome ids of the input set stand in the record")

    def test_fault_a_base_run_in_a_cell(self):
        i = self.copy_inputs()
        self.edit_compact(i, lambda reports: reports[5].update(antibiotic="acgt" * 6 + "a"))
        self.stops(i, "text of the shape")

    # ------------------------------------------------------------- --check

    def test_check_with_inputs_finds_a_changed_call(self):
        i = self.copy_inputs()

        def change(rows):
            r = next(r for r in rows if r["kind"] == "primary" and r["who_assoc"] == "susceptible")
            r.update(who_assoc="no_call", who_assoc_only="no_call", tool_default="no_call", who_assoc_small="no_call")
        self.edit_calls(i, change)
        code, _, err = run("--inputs", i, "--check", self.out)
        self.assertEqual(code, 1)
        self.assertIn("differs from the input set's record", err)
        self.assertIn("callers.tsv", err)

    def check_alone_fails(self, rec, words):
        code, _, err = run("--check", rec)
        self.assertEqual(code, 1, err)
        self.assertIn(words, err)

    def edit_record(self, rec, name, change):
        p = rec / name
        lines = p.read_text().split("\n")
        change(lines)
        p.write_text("\n".join(lines))

    def test_check_alone_finds_a_changed_derived_cell(self):
        r = self.copy_record()
        self.edit_record(r, "screen.tsv", lambda lines: lines.__setitem__(1, lines[1][:-1] + "9"))
        self.check_alone_fails(r, "screen.tsv")

    def test_check_alone_finds_a_changed_interval_end(self):
        r = self.copy_record()

        def change(lines):
            at = next(k for k, x in enumerate(lines) if x.startswith("median\tby_name\t-\ttool\tsensitivity\t"))
            cells = lines[at].split("\t")
            cells[5] = repr(float(cells[5]) / 2)
            lines[at] = "\t".join(cells)
        self.edit_record(r, "intervals.tsv", change)
        self.check_alone_fails(r, "medians.tsv")

    def test_check_alone_finds_a_changed_count(self):
        r = self.copy_record()

        def change(lines):
            cells = lines[4].split("\t")
            cells[7] = str(int(cells[7]) + 1)  # tp
            lines[4] = "\t".join(cells)
        self.edit_record(r, "callers.tsv", change)
        self.check_alone_fails(r, "the confusion counts do not add up")

    def test_check_alone_finds_a_dropped_file(self):
        r = self.copy_record()
        (r / "screen.tsv").unlink()
        self.check_alone_fails(r, "screen.tsv")

    def test_check_alone_finds_an_added_file(self):
        r = self.copy_record()
        (r / "genomes.tsv").write_text("genome\n")
        self.check_alone_fails(r, "genomes.tsv")

    def test_check_alone_finds_a_changed_unbuilt_file(self):
        r = self.copy_record()
        p = r / "inputs.sha256"
        p.write_text(p.read_text().replace("0", "1", 1))
        self.check_alone_fails(r, "SHA256SUMS")

    def test_check_alone_finds_an_added_row(self):
        r = self.copy_record()
        self.edit_record(r, "callers.tsv", lambda lines: lines.insert(2, lines[1]))
        self.check_alone_fails(r, "not one row per way, report and caller")

    def test_check_alone_finds_a_row_per_genome(self):
        r = self.copy_record()

        def change(lines):
            cells = lines[1].split("\t")
            cells[2] = "genome"
            lines.insert(2, "\t".join(cells))
        self.edit_record(r, "callers.tsv", change)
        self.check_alone_fails(r, "not one row per way, report and caller")

    def test_check_alone_finds_a_proxy_moved_into_a_set(self):
        r = self.copy_record()
        self.edit_record(r, "reports.tsv", lambda lines: lines.__setitem__(
            next(k for k, x in enumerate(lines) if "\tproxy\t" in x),
            next(x for x in lines if "\tproxy\t" in x).replace("\tproxy\t", "\tprimary\t")))
        self.check_alone_fails(r, "reports.tsv")

    def test_check_alone_finds_a_genome_column(self):
        r = self.copy_record()

        def change(lines):
            lines[0] += "\tgenome"
            lines[1:-1] = [x + "\t" + synthetic.genome_id(0) for x in lines[1:-1]]
        self.edit_record(r, "reports.tsv", change)
        self.check_alone_fails(r, "reports.tsv: the columns are not")

    def test_check_alone_finds_a_long_float(self):
        r = self.copy_record()

        def change(lines):
            cells = lines[1].split("\t")
            cells[5] = cells[5] + "0"
            lines[1] = "\t".join(cells)
        self.edit_record(r, "intervals.tsv", change)
        self.check_alone_fails(r, "not in its shortest form")


if __name__ == "__main__":
    unittest.main()
