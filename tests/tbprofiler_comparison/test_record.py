"""The committed record under reference/tbprofiler-comparison: its files and sums, its columns and rows, its size,
the two sets and the proxies, what it must never hold, and the script's --check mode on it."""

import collections
import contextlib
import hashlib
import io
import pathlib
import re
import shutil
import statistics
import tempfile
import unittest

from tests.tbprofiler_comparison.common import LEARNERS, RECORD, SPLITS, WAY, load_script, pairs, read_tsv

record = load_script()
MAX_BYTES = 3 * 1024 * 1024
# The sha256 of SHA256SUMS as committed. Any change to a file of the record changes it.
PINNED = "a98a5a6eea8ec2e61141a008d60e7b24911109f00380fe3097bb4636afd1af2a"
REPORTS, NAMED, PROXIES, WAYS, TOOL_DRUGS = 22, 20, 4, 4, 18
CALLERS = 1 + len(LEARNERS) * len(SPLITS)
SETS, MEASURES = 2, 3
ROWS = {"reports.tsv": REPORTS, "callers.tsv": WAYS * NAMED * CALLERS,
        "intervals.tsv": NAMED * (CALLERS * MEASURES + CALLERS - 1) + SETS * (CALLERS * MEASURES + CALLERS - 1),
        "disagreements.tsv": NAMED, "grades.tsv": TOOL_DRUGS, "order.tsv": 5, "measures.tsv": WAYS * NAMED * CALLERS,
        "endpoints.tsv": WAYS * NAMED, "medians.tsv": WAYS * SETS * CALLERS * MEASURES,
        "comparisons.tsv": WAYS * SETS * (CALLERS - 1), "screen.tsv": WAYS * SETS, "bounds.tsv": WAYS * SETS * 4}
# Column names that would mean a row per genome, or the content of a call or a catalogue entry. A column named
# genomes is a count, and test_cells_follow_from_each_other holds it to the other counts of its row.
BANNED_COLUMN = re.compile(r"^(genome|sample|isolate|accession|assembly|bioproject|biosample|study|label|score|kmer|"
                           r"feature|marker|gene|mutation|variant|weight|importance|sequence|call)(_ids?|_names?)?$"
                           r"|^(samples|isolates|accessions|assemblies|labels|scores|kmers|features|markers|genes|"
                           r"mutations|variants|sequences)$|(^|_)ids?$")
BANNED_PATH = re.compile(r"(?i)seal[e]d|held-out|ks[0-9]{2}|vum[c]|clini[c]al|priv[a]te|ks[d]|sigdr[u]g")
INPUT_NAMES = ["calls.tsv", "rev14_compact.jsonl.gz", "database-grades-by-drug.tsv", "PLAN.md", "PLAN-ADDENDUM.md",
               "PLAN.sha256", "finish.log", "exploratory.json"]


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
    @classmethod
    def setUpClass(cls):
        cls.reports = read_tsv(RECORD / "reports.tsv")
        cls.callers = read_tsv(RECORD / "callers.tsv")
        cls.endpoints = {(r["way"], r["dataset"]): r for r in read_tsv(RECORD / "endpoints.tsv")}

    def test_columns(self):
        for name, cols in record.COLUMNS.items():
            head = (RECORD / name).read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
            self.assertEqual(head, cols, name)

    def test_rows(self):
        for name, count in ROWS.items():
            self.assertEqual(len(read_tsv(RECORD / name)), count, name)
        self.assertEqual(len((RECORD / "inputs.sha256").read_text().rstrip("\n").split("\n")), len(INPUT_NAMES))

    def test_one_row_per_way_report_and_caller(self):
        named = [r for r in self.reports if r["kind"] != "not named"]
        want = [(w, r["dataset"], c) for w in record.WAYS for r in named
                for c in ("tool",) + record.LEARNER_CALLERS]
        self.assertEqual([(r["way"], r["dataset"], r["caller"]) for r in self.callers], want)
        for r in self.reports:
            self.assertRegex(r["dataset"], r"^mycobacterium_tuberculosis__[a-z0-9_]+$")
            self.assertEqual(r["dataset"], "mycobacterium_tuberculosis__" + re.sub(r"[^a-z0-9]+", "_", r["antibiotic"]))

    def test_kinds(self):
        kinds = {}
        for r in self.reports:
            kinds.setdefault(r["kind"], []).append(r["antibiotic"])
        self.assertEqual({k: len(v) for k, v in kinds.items()}, {"primary": 16, "proxy": PROXIES, "not named": 2})
        self.assertEqual(kinds["not named"], ["amoxicillin", "nicotinamide"])
        self.assertEqual({r["antibiotic"]: r["genomes"] for r in self.reports if r["kind"] == "not named"},
                         {"amoxicillin": "57", "nicotinamide": "166"})

    def test_the_plans_16_and_the_addendums_14(self):
        by_name = [d for (w, d), e in self.endpoints.items() if w == WAY and e["by_name"] == "yes"]
        graded = [d for (w, d), e in self.endpoints.items() if w == WAY and e["graded"] == "yes"]
        kind = {r["dataset"]: r["kind"] for r in self.reports}
        self.assertEqual(len(by_name), 16)
        self.assertEqual(len(graded), 14)
        self.assertTrue(set(graded) <= set(by_name))
        self.assertEqual({kind[d] for d in by_name}, {"primary"})
        self.assertEqual(sorted(d.split("__")[1] for d in set(by_name) - set(graded)),
                         ["cycloserine", "para_aminosalicylic_acid"])
        readme = (RECORD / "README.md").read_text(encoding="utf-8")
        self.assertIn("`by_name` holds the scored primary reports, the plan's set; `graded` holds those whose tool drug "
                      "the database grades, the addendum's set. Proxy reports enter no set.", readme)

    def test_proxies_enter_no_median(self):
        kind = {r["dataset"]: r["kind"] for r in self.reports}
        for (w, d), e in self.endpoints.items():
            if kind[d] == "proxy":
                self.assertEqual((e["by_name"], e["graded"]), ("no", "no"), (w, d))
        measures = {(r["way"], r["dataset"], r["caller"]): r for r in read_tsv(RECORD / "measures.tsv")}
        for r in read_tsv(RECORD / "medians.tsv"):
            members = [d for (w, d), e in self.endpoints.items() if w == r["way"] and e[r["set"]] == "yes"]
            self.assertEqual(int(r["endpoints"]), len(members))
            values = [float(measures[(r["way"], d, r["caller"])][r["measure"]]) for d in members]
            self.assertEqual(float(r["median"]), statistics.median(values), (r["way"], r["set"], r["caller"]))

    def test_checks(self):
        checks = {r["check"]: int(r["count"]) for r in read_tsv(RECORD / "checks.tsv")}
        self.assertEqual(checks, {"reports_whose_genomes_equal_the_call_file": NAMED, "plans_match_their_hash": 2,
                                  "call_rows_read": 23968, "genome_ids_of_the_inputs_in_the_record": 0})

    def test_settings(self):
        s = pairs(RECORD / "settings.tsv")
        for k, v in record.FIXED_SETTINGS.items():
            self.assertEqual(s[k], v, k)
        self.assertEqual((s["resamples"], s["seed"], s["low_rank"], s["high_rank"]), ("2000", "20261114", "51", "1950"))
        for k in ("plan_sha256", "addendum_sha256", "compact_sha256"):
            self.assertRegex(s[k], r"^[0-9a-f]{64}$")

    def test_cells_follow_from_each_other(self):
        tool = {}
        for r in self.callers:
            c = {k: int(r[k]) for k in record.COUNT_COLS + ["caller_only_right", "tool_only_right"]}
            self.assertEqual(c["tp"] + c["fp"] + c["tn"] + c["fn"], c["with_call"])
            self.assertEqual(c["with_call"] + c["no_call_resistant"] + c["no_call_susceptible"], c["genomes"])
            if r["caller"] == "tool":
                tool[(r["way"], r["dataset"])] = c
            else:
                t = tool[(r["way"], r["dataset"])]
                self.assertEqual((c["tp"] + c["fn"], c["tn"] + c["fp"]), (t["tp"] + t["fn"], t["tn"] + t["fp"]))
                self.assertEqual(c["tp"] + c["tn"] - t["tp"] - t["tn"], c["caller_only_right"] - c["tool_only_right"])


class NeverHeld(unittest.TestCase):
    def test_no_column_per_genome(self):
        for name in record.COLUMNS:
            for c in (RECORD / name).read_text(encoding="utf-8").split("\n", 1)[0].split("\t"):
                self.assertIsNone(BANNED_COLUMN.search(c), (name, c))
                if c == "genomes":
                    self.assertTrue(all(re.fullmatch(r"0|[1-9]\d*", r[c]) for r in read_tsv(RECORD / name)), name)

    def test_no_text_of_a_banned_shape(self):
        for p in sorted(RECORD.iterdir()):
            text = p.read_text(encoding="utf-8")
            for shape in record.BANNED_SHAPES:
                self.assertIsNone(re.search(shape, text), (p.name, shape))

    def test_no_row_per_genome(self):
        # No file has more rows than the record's own shape allows: a row per way, report and caller, or per report
        # and set, caller and measure with the difference; far fewer than the genomes.
        genomes = int(pairs(RECORD / "totals.tsv")["distinct_genomes"])
        named = sum(1 for r in read_tsv(RECORD / "reports.tsv") if r["kind"] != "not named")
        callers = 1 + len(record.LEARNER_CALLERS)
        shape = max(len(record.WAYS) * named * callers,
                    (named + len(record.SETS)) * callers * (len(record.MEASURES) + 1))
        self.assertLess(shape, genomes)
        for p in sorted(RECORD.glob("*.tsv")):
            self.assertLessEqual(len(read_tsv(p)), shape, p.name)
        self.assertEqual(len(read_tsv(RECORD / "callers.tsv")), len(record.WAYS) * named * callers)

    def test_inputs_are_the_named_files(self):
        lines = (RECORD / "inputs.sha256").read_text(encoding="utf-8").rstrip("\n").split("\n")
        names = []
        for line in lines:
            digest, name = line.split("  ")
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
            self.assertIsNone(BANNED_PATH.search(name), name)
            names.append(name)
        self.assertEqual(names, INPUT_NAMES)


CATALOG_NAMES = {"catalog_endpoints", "catalog_rows", "catalog_tool_only_right", "catalog_only_right",
                 "tool_minus_catalog_median"}
# The groups of 1 the record gives outside the catalog row: the genomes without a call of 3 reports, by label, under
# each way in callers.tsv and under the primary way in disagreements.tsv. The public release names no such group by
# itself. They are counted here, so that no other group of 1 enters unseen.
GROUPS_OF_ONE = {"tool genomes without a call, by label": 30}


def small_cells(root=RECORD):
    """Every count of the record with the group the record gives it. A cell describes one genome when the group it
    counts has one member: a count of 1, or of 0, in a group of 1. Returns the count of such cells by kind."""
    found = collections.Counter()

    def check(kind, count, klass):
        if klass == 1:
            found[kind] += 1

    for r in read_tsv(root / "callers.tsv"):
        c = {k: int(r[k]) for k in record.COUNT_COLS + ["caller_only_right", "tool_only_right"]}
        res, sus = c["tp"] + c["fn"], c["tn"] + c["fp"]
        who = "tool" if r["caller"] == "tool" else "learner"
        for k, klass in (("tp", res), ("fn", res), ("tn", sus), ("fp", sus)):
            check(f"{who} {k}", c[k], klass)
        if who == "tool":
            check("tool genomes with a call", c["with_call"], c["genomes"])
            check("tool resistant genomes with a call", res, res + c["no_call_resistant"])
            check("tool susceptible genomes with a call", sus, sus + c["no_call_susceptible"])
            for k in ("no_call_resistant", "no_call_susceptible"):
                check("tool genomes without a call, by label", c[k], c["no_call_resistant"] + c["no_call_susceptible"])
        else:
            # A pair count is drawn from the genomes the tool called, which the learner's row is on.
            for k in ("caller_only_right", "tool_only_right"):
                check(f"learner {k}", c[k], c["with_call"])
    for r in read_tsv(root / "disagreements.tsv"):
        c = {k: int(v) for k, v in r.items() if k != "dataset"}
        for a, b in (("tool_r_label_s_assoc", "tool_r_label_s_interim"), ("tool_s_label_r_other", "tool_s_label_r_none"),
                     ("both_r_assoc", "both_r_interim"), ("both_s_other", "both_s_none")):
            for k in (a, b):
                check("grade split of a call", c[k], c[a] + c[b])
        for k in ("no_call_label_r", "no_call_label_s"):
            check("tool genomes without a call, by label", c[k], c["no_call_label_r"] + c["no_call_label_s"])
    for r in read_tsv(root / "comparisons.tsv"):
        for k in ("caller_only_right", "tool_only_right"):
            check(f"comparison {k}", int(r[k]), int(r["caller_only_right"]) + int(r["tool_only_right"]))
    s = pairs(root / "summary.tsv")
    for key in ("by_name", "graded"):
        check("rows with a call", int(s[f"{key}_rows_with_call"]), int(s[f"{key}_rows"]))
        check("rows without a call", int(s[f"{key}_rows_without_call"]), int(s[f"{key}_rows"]))
        check("rows without a call, resistant", int(s[f"{key}_rows_without_call_label_resistant"]),
              int(s[f"{key}_rows_without_call"]))
    for k in ("catalog_calls_agree", "catalog_tool_only_right", "catalog_only_right"):
        check(f"catalog {k}", int(s[k]), int(s["catalog_rows"]))
    a = pairs(root / "set-aside.tsv")
    for name in ("as_run", "counted"):
        tp, fp, tn, fn = (int(a[f"{k}_{name}"]) for k in ("tp", "fp", "tn", "fn"))
        for k, v, klass in (("tp", tp, tp + fn), ("fn", fn, tp + fn), ("tn", tn, tn + fp), ("fp", fp, tn + fp)):
            check(f"set-aside {k}", v, klass)
    check("set-aside entry", int(a["tn_counted"]), int(a["tn_as_run"]))
    return found


class SmallCells(unittest.TestCase):
    def test_no_catalog_row_per_report(self):
        # The catalog row stands only as figures over the plan's set: no caller row of it in any file.
        for p in sorted(RECORD.glob("*.tsv")):
            rows = read_tsv(p)
            if rows and "caller" in rows[0]:
                self.assertEqual({r["caller"] for r in rows} - {"tool", "-"}, set(record.LEARNER_CALLERS), p.name)
        for name in ("totals.tsv", "summary.tsv"):
            found = {k for k in pairs(RECORD / name) if "catalog" in k}
            self.assertEqual(found, CATALOG_NAMES | ({"catalog_calls_agree"} if name == "summary.tsv" else set()), name)
        for p in sorted(RECORD.iterdir()):
            if p.name in ("totals.tsv", "summary.tsv", "README.md", "reports.tsv"):
                continue
            self.assertNotIn("catalog", p.read_text(encoding="utf-8"), p.name)
        self.assertEqual({r["catalog_row"] for r in read_tsv(RECORD / "reports.tsv")}, {"yes", "no"})
        # The median is given to 3 decimals: in full it would show whether one genome had a call.
        self.assertRegex(pairs(RECORD / "totals.tsv")["tool_minus_catalog_median"], r"^-?\d\.\d{3}$")
        self.assertEqual(pairs(RECORD / "summary.tsv")["tool_minus_catalog_median"],
                         pairs(RECORD / "totals.tsv")["tool_minus_catalog_median"])
        readme = (RECORD / "README.md").read_text(encoding="utf-8")
        self.assertIn("The catalog row appears only as six-report totals, in `totals.tsv` and again in `summary.tsv`, "
                      "because on three reports a per-report count of it would describe one genome.", readme)

    def test_no_cell_describes_one_genome(self):
        self.assertEqual(dict(small_cells()), GROUPS_OF_ONE)

    def test_the_small_cell_check_finds_its_cells(self):
        # The check itself, on cells made up here: a count of 1 in a group of 1.
        tmp = tempfile.mkdtemp(prefix="small-cells-")
        try:
            dest = pathlib.Path(shutil.copytree(RECORD, f"{tmp}/rec"))
            p = dest / "summary.tsv"
            text = p.read_text(encoding="utf-8")
            s = pairs(p)
            text = text.replace(f"by_name_rows_without_call\t{s['by_name_rows_without_call']}\n",
                                "by_name_rows_without_call\t1\n")
            text = text.replace(f"by_name_rows_without_call_label_resistant\t"
                                f"{s['by_name_rows_without_call_label_resistant']}\n",
                                "by_name_rows_without_call_label_resistant\t1\n")
            p.write_text(text, encoding="utf-8")
            self.assertEqual(small_cells(dest)["rows without a call, resistant"], 1)
        finally:
            shutil.rmtree(tmp)


if __name__ == "__main__":
    unittest.main()
