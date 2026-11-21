"""The committed record under reference/tbprofiler-comparison: its files and sums, its columns and rows, its size,
the two sets and the proxies, what it must never hold, and the script's --check mode on it."""

import contextlib
import hashlib
import io
import re
import statistics
import unittest

from tests.tbprofiler_comparison.common import LEARNERS, RECORD, SPLITS, WAY, load_script, pairs, read_tsv

record = load_script()
MAX_BYTES = 3 * 1024 * 1024
# The sha256 of SHA256SUMS as committed. Any change to a file of the record changes it.
PINNED = "57f2bb8fbcc6cfcc1af955ca2006f683676d4e9444af7424bf8fac899ed03ee5"
REPORTS, NAMED, PROXIES, CATALOG, WAYS, TOOL_DRUGS = 22, 20, 4, 8, 4, 18
CALLERS = 1 + len(LEARNERS) * len(SPLITS)
SETS, MEASURES = 2, 3
ROWS = {"reports.tsv": REPORTS, "callers.tsv": WAYS * (NAMED * CALLERS + CATALOG),
        "intervals.tsv": NAMED * (CALLERS * MEASURES + CALLERS - 1) + SETS * (CALLERS * MEASURES + CALLERS - 1),
        "disagreements.tsv": NAMED, "grades.tsv": TOOL_DRUGS, "order.tsv": 5, "measures.tsv": WAYS * (NAMED * CALLERS + CATALOG),
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
                for c in ("tool",) + record.LEARNER_CALLERS + (("catalog",) if r["catalog_row"] == "yes" else ())]
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
        genomes = int(pairs(RECORD / "totals.tsv")["distinct_genomes"])
        for p in sorted(RECORD.glob("*.tsv")):
            self.assertLess(len(read_tsv(p)), genomes, p.name)

    def test_inputs_are_the_named_files(self):
        lines = (RECORD / "inputs.sha256").read_text(encoding="utf-8").rstrip("\n").split("\n")
        names = []
        for line in lines:
            digest, name = line.split("  ")
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
            self.assertIsNone(BANNED_PATH.search(name), name)
            names.append(name)
        self.assertEqual(names, INPUT_NAMES)


if __name__ == "__main__":
    unittest.main()
