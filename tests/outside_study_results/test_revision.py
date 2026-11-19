"""The outside-study estimates read against a later release's figures.

scripts/outside-study-results/revision.py reads reference/outside-study against
a figure file under reference/outside-study-figures. The suite checks the
committed figure file against the runs' own rows, the reading rule against the
kit's readings, and the tables of the page's section against the files. The
reports a figure file is read from are not in the repository, so `--figures`
is tested on made-up reports.
"""

import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import pathlib
import re
import shutil
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reference" / "outside-study"
FIGURES = ROOT / "reference" / "outside-study-figures" / "rev14-274-20261112.tsv"
# The release's own checksum lines for the reports the figure file was read from.
REPORT_LINES = FIGURES.parent / "rev14-274-20261112.reports.sha256"
README = FIGURES.parent / "README.md"
RELEASES = ROOT / "reference" / "release-summary" / "README.md"
DOCUMENT = ROOT / "docs" / "outside-study-results.md"


def load_revision():
    spec = importlib.util.spec_from_file_location("outside_study_results_revision", ROOT / "scripts" / "outside-study-results" / "revision.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_revision()
T = R.tables


def run_rows(data):
    """Every row of every run's metrics.tsv, read or not."""
    return [r for runs in data.values() for _, _, metrics, _ in runs.values() for r in metrics]


# The least and the greatest figure strictly between 0 and 1.
INSIDE_LOW, INSIDE_HIGH = R.INSIDE


def metrics(n, tp, fn, fp, tn):
    return {"n": n, "tp": tp, "fn": fn, "fp": fp, "tn": tn, "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp)}


class FigureFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = T.load(RESULTS)
        cls.revision, cls.figures = R.load_figures(FIGURES)

    def test_one_row_for_each_dataset_and_learner_of_the_runs(self):
        self.assertEqual(set(self.figures), {(r["dataset"], r["learner"]) for r in run_rows(self.data)})

    def test_each_dataset_has_one_report_of_one_revision(self):
        self.assertEqual(self.revision, 14)
        digests = {}
        for (dataset, _), f in self.figures.items():
            self.assertRegex(f["report_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(digests.setdefault(dataset, f["report_sha256"]), f["report_sha256"], dataset)
        # No two datasets share a report.
        self.assertEqual(len(set(digests.values())), len(digests))

    def test_each_report_hash_is_its_line_in_the_release_s_checksum_list(self):
        lines = REPORT_LINES.read_text(encoding="utf-8").split("\n")
        self.assertEqual(lines.pop(), "")
        listed = {}
        for line in lines:
            found = re.fullmatch(r"([0-9a-f]{64})  \./datasets/([a-z0-9_]+)\.json", line)
            self.assertIsNotNone(found, line)
            self.assertNotIn(found.group(2), listed)
            listed[found.group(2)] = found.group(1)
        self.assertEqual(listed, {dataset: f["report_sha256"] for (dataset, _), f in self.figures.items()})

    def test_each_figure_is_its_counts(self):
        with_figure = 0
        for key, f in self.figures.items():
            if f["lineage_note"]:
                self.assertEqual([f[c] for c in R.COLUMNS[4:11]], [""] * 7, key)
                continue
            with_figure += 1
            n, tp, fn, fp, tn = (int(f[f"lineage_{k}"]) for k in R.COUNTS)
            self.assertEqual(n, tp + fn + fp + tn, key)
            self.assertAlmostEqual(float(f["lineage_sensitivity"]), tp / (tp + fn), places=12, msg=key)
            self.assertAlmostEqual(float(f["lineage_specificity"]), tn / (tn + fp), places=12, msg=key)
        self.assertGreater(with_figure, 0)

    def test_a_dataset_with_no_lineage_split_has_the_kit_s_note_in_both(self):
        noted = {key: f["lineage_note"] for key, f in self.figures.items() if f["lineage_note"]}
        self.assertEqual(noted, {(r["dataset"], r["learner"]): r["lineage_note"] for r in run_rows(self.data) if r["lineage_note"]})
        self.assertTrue(noted)

    def test_a_dataset_in_both_populations_has_one_figure_in_the_runs(self):
        seen = {}
        for r in run_rows(self.data):
            self.assertEqual(seen.setdefault((r["dataset"], r["learner"]), R.counts(r)), R.counts(r))

    def test_the_readme_names_the_release_and_the_file(self):
        text = README.read_text(encoding="utf-8")
        for path in (FIGURES, REPORT_LINES):
            self.assertIn(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n", text)
        # The release is the one reference/release-summary names for this revision.
        (row,) = [line for line in RELEASES.read_text(encoding="utf-8").split("\n") if line.startswith("| `rev14-274-20261112/` | `rev14-274-20261112.tar.zst`")]
        asset = row.split("|")[3].strip().strip("`")
        self.assertRegex(asset, r"^[0-9a-f]{64}$")
        self.assertIn(asset, text)
        self.assertIn(f"{len({d for d, _ in self.figures})} reports", text)
        self.assertIn(f"{len(self.figures)} rows", text)


class Rereading(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = T.load(RESULTS)
        cls.revision, cls.figures = R.load_figures(FIGURES)
        cls.later = R.later(cls.data, cls.figures)

    def test_the_rule_on_the_printed_interval_gives_the_kit_s_reading_on_every_row_read(self):
        """The control of the method: the run's own figure, read from the printed interval."""
        checked = 0
        for rows in self.later.values():
            for r in rows:
                tp, fn, fp, tn = (int(r[f"lineage_{k}"]) for k in ("tp", "fn", "fp", "tn"))
                own = {"sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp)}
                for m in R.MEASURES:
                    got = R.score.reading(own[m], float(r[m + "_lo"]), float(r[m + "_hi"]))
                    self.assertEqual(got, r["reading_" + m], (r["dataset"], r["learner"], m))
                    # The script's own reading stops on none of them and gives the same.
                    self.assertEqual(R.reread(r, m, own[m]), r["reading_" + m], (r["dataset"], r["learner"], m))
                    checked += 1
        self.assertEqual(checked, 2 * sum(len(rows) for rows in self.later.values()))
        self.assertGreater(checked, 0)

    def test_a_run_s_figure_of_exactly_0_or_1_at_a_printed_end_is_read_as_the_kit_read_it(self):
        """The runs' own rows where a report figure of exactly 0 or 1 sits on an end printed as 0.000000 or 1.000000.

        In each the new estimate is at its limit too, 0 of n or n of n, so the
        kit set that end itself and the reading is certain."""
        met = 0
        for rows in self.later.values():
            for r in rows:
                for m, hit, miss in (("sensitivity", "tp", "fn"), ("specificity", "tn", "fp")):
                    figure = int(r["lineage_" + hit]) / (int(r["lineage_" + hit]) + int(r["lineage_" + miss]))
                    lo, hi = float(r[m + "_lo"]), float(r[m + "_hi"])
                    if figure not in (0.0, 1.0) or min(abs(figure - lo), abs(figure - hi)) >= R.MARGIN:
                        continue
                    met += 1
                    self.assertEqual(R.exact_ends(r, m), (figure == 0.0, figure == 1.0), (r["dataset"], r["learner"], m))
                    self.assertEqual((R.reread(r, m, figure), r["reading_" + m]), ("consistent", "consistent"), (r["dataset"], r["learner"], m))
        self.assertGreater(met, 0)

    def test_no_run_scored_a_class_large_enough_for_a_computed_end_to_print_at_the_limit(self):
        """From `least` genomes of a class on, the kit's computed end for 1 of n prints as 0.000000."""
        fmt, interval = R.score.fmt, R.score.clopper_pearson
        least = 50636
        self.assertEqual((fmt(interval(1, least)[0]), fmt(interval(1, least - 1)[0])), ("0.000000", "0.000001"))
        self.assertEqual((fmt(interval(least - 1, least)[1]), fmt(interval(least - 2, least - 1)[1])), ("1.000000", "0.999999"))
        largest = max(int(r[k]) for rows in self.later.values() for r in rows for k in ("scored_resistant", "scored_susceptible"))
        self.assertLess(largest, least)
        for rows in self.later.values():
            for r in rows:
                for m in R.MEASURES:
                    at_zero, at_one = R.exact_ends(r, m)
                    # An end printed at the limit is one the kit set, on every row read.
                    self.assertEqual((float(r[m + "_lo"]) == 0.0, float(r[m + "_hi"]) == 1.0), (at_zero, at_one), (r["dataset"], r["learner"], m))

    def test_the_rows_read_are_the_ones_the_pages_read(self):
        for pop, runs in self.data.items():
            want = [(r["dataset"], r["learner"]) for _, _, rows, _ in runs.values() for r in rows if r["learner"] in T.LEARNERS and T.holds(r) is not None]
            self.assertEqual([(r["dataset"], r["learner"]) for r in self.later[pop]], want)

    def test_a_row_with_the_run_s_counts_keeps_the_run_s_readings(self):
        same = [r for rows in self.later.values() for r in rows if not r["differs"]]
        self.assertTrue(same)
        for r in same:
            self.assertEqual(r["later"], {m: r["reading_" + m] for m in R.MEASURES})
            self.assertEqual(R.holds_later(r), T.holds(r))

    def test_a_row_whose_counts_differ_is_read_on_the_later_figure(self):
        differing = [r for rows in self.later.values() for r in rows if r["differs"]]
        self.assertTrue(differing)
        for r in differing:
            for m in R.MEASURES:
                figure = float(r["figure"]["lineage_" + m])
                lo, hi = float(r[m + "_lo"]), float(r[m + "_hi"])
                want = "below" if figure > hi else "above" if figure < lo else "consistent"
                self.assertEqual(r["later"][m], want, (r["dataset"], m))
                # Far from the six decimals an interval end is printed to.
                self.assertGreater(min(abs(figure - lo), abs(figure - hi)), 1000 * R.MARGIN)

    def test_a_figure_at_a_printed_interval_end_is_refused(self):
        r = next(r for rows in self.later.values() for r in rows if r["differs"])
        hi = float(r["sensitivity_hi"])
        # The kit computed this end: the estimate is not at its limit.
        self.assertEqual(R.exact_ends(r, "sensitivity"), (False, False))
        self.assertTrue(0.0 < hi - 2e-6 and hi + 2e-6 < 1.0)
        # Six printed decimals leave an end known to half a millionth.
        for figure in (hi, hi + 5e-7, hi - 5e-7, hi + 9e-7):
            with self.assertRaises(SystemExit):
                R.reread(r, "sensitivity", figure)
        self.assertEqual(R.reread(r, "sensitivity", hi + 2e-6), "below")
        self.assertEqual(R.reread(r, "sensitivity", hi - 2e-6), "consistent")

    def test_a_row_read_by_a_run_with_no_later_figure_is_refused(self):
        key = next((r["dataset"], r["learner"]) for rows in self.later.values() for r in rows)
        blank = {**self.figures[key], "lineage_sensitivity": ""}
        for figures in ({**self.figures, key: blank}, {k: f for k, f in self.figures.items() if k != key}):
            with self.assertRaises(SystemExit):
                R.later(self.data, figures)


class PrintedEnds(unittest.TestCase):
    """The stop at a printed interval end, on made-up rows. The counts are no dataset's.

    Each row is one measure of a made-up `metrics.tsv` row: some genomes of a
    class called right, and the kit's own interval for that count as the kit
    prints it. A class here has more genomes than any run scored, because only
    then does a computed end print as 0.000000 or 1.000000."""

    N = 60000
    MANY = 8000000

    def row(self, measure, right, n):
        """(the made-up row, the kit's interval before it is printed)."""
        lo, hi = R.score.clopper_pearson(right, n)
        hit, miss = ("tp", "fn") if measure == "sensitivity" else ("tn", "fp")
        row = {"dataset": "made_up__drug", "learner": "lasso", hit: str(right), miss: str(n - right)}
        row[measure + "_lo"], row[measure + "_hi"] = R.score.fmt(lo), R.score.fmt(hi)
        return row, (lo, hi)

    def kit(self, figure, interval):
        """The kit's reading on the interval it computed, before printing."""
        return R.score.reading(figure, *interval)

    def test_a_figure_of_exactly_0_is_read_from_the_counts_on_both_sides(self):
        for m in R.MEASURES:
            # The kit set the lower end to 0: the figure is not under it.
            exact, interval = self.row(m, 0, self.N)
            self.assertEqual((exact[m + "_lo"], interval[0]), ("0.000000", 0.0))
            self.assertEqual((R.reread(exact, m, 0.0), self.kit(0.0, interval)), ("consistent", "consistent"))
            # The kit computed the lower end, above 0, and it prints the same six decimals: the figure is under it.
            computed, interval = self.row(m, 1, self.N)
            self.assertEqual(computed[m + "_lo"], "0.000000")
            self.assertGreater(interval[0], 0.0)
            self.assertEqual((R.reread(computed, m, 0.0), self.kit(0.0, interval)), ("above", "above"))
            # The rule on the printed digits alone would read the second row as it reads the first.
            self.assertEqual(R.score.reading(0.0, float(computed[m + "_lo"]), float(computed[m + "_hi"])), "consistent")

    def test_a_figure_of_exactly_1_is_read_from_the_counts_on_both_sides(self):
        for m in R.MEASURES:
            exact, interval = self.row(m, self.N, self.N)
            self.assertEqual((exact[m + "_hi"], interval[1]), ("1.000000", 1.0))
            self.assertEqual((R.reread(exact, m, 1.0), self.kit(1.0, interval)), ("consistent", "consistent"))
            computed, interval = self.row(m, self.N - 1, self.N)
            self.assertEqual(computed[m + "_hi"], "1.000000")
            self.assertLess(interval[1], 1.0)
            self.assertEqual((R.reread(computed, m, 1.0), self.kit(1.0, interval)), ("below", "below"))
            self.assertEqual(R.score.reading(1.0, float(computed[m + "_lo"]), float(computed[m + "_hi"])), "consistent")

    def test_a_figure_at_its_limit_is_never_beyond_the_far_end(self):
        """A computed upper end that prints as 0.000000, and a computed lower end that prints as 1.000000."""
        for m in R.MEASURES:
            low, interval = self.row(m, 0, self.MANY)
            self.assertEqual((low[m + "_hi"], interval[1] > 0.0), ("0.000000", True))
            self.assertEqual((R.reread(low, m, 0.0), self.kit(0.0, interval)), ("consistent", "consistent"))
            self.assertEqual((R.reread(low, m, 1.0), self.kit(1.0, interval)), ("below", "below"))
            high, interval = self.row(m, self.MANY, self.MANY)
            self.assertEqual((high[m + "_lo"], interval[0] < 1.0), ("1.000000", True))
            self.assertEqual((R.reread(high, m, 1.0), self.kit(1.0, interval)), ("consistent", "consistent"))
            self.assertEqual((R.reread(high, m, 0.0), self.kit(0.0, interval)), ("above", "above"))

    def test_the_stop_stays_where_the_printed_digits_could_hide_the_reading(self):
        for m in R.MEASURES:
            # A computed lower end of about 4.2e-07, printed 0.000000.
            low, interval = self.row(m, 1, self.N)
            self.assertTrue(4e-7 < interval[0] < 5e-7)
            # The kit reads one of these two figures above and the other consistent, and the printed digits are the same.
            self.assertEqual((self.kit(3e-7, interval), self.kit(4.5e-7, interval)), ("above", "consistent"))
            for figure in (3e-7, 4.5e-7, 9e-7, INSIDE_LOW):
                with self.assertRaises(SystemExit) as stop:
                    R.reread(low, m, figure)
                self.assertIn("within one millionth", str(stop.exception))
            self.assertEqual((R.reread(low, m, 2e-6), self.kit(2e-6, interval)), ("consistent", "consistent"))
            # The mirror: a computed upper end printed 1.000000.
            high, interval = self.row(m, self.N - 1, self.N)
            self.assertEqual((self.kit(1 - 3e-7, interval), self.kit(1 - 4.5e-7, interval)), ("below", "consistent"))
            for figure in (1 - 3e-7, 1 - 4.5e-7, 1 - 9e-7, INSIDE_HIGH):
                with self.assertRaises(SystemExit):
                    R.reread(high, m, figure)
            self.assertEqual((R.reread(high, m, 1 - 2e-6), self.kit(1 - 2e-6, interval)), ("consistent", "consistent"))

    def test_no_stop_at_an_end_the_kit_set(self):
        """An end the kit set is printed exactly, so a figure beside it is read with no margin."""
        for m in R.MEASURES:
            low, interval = self.row(m, 0, self.N)
            for figure in (INSIDE_LOW, 3e-7, 9e-7):
                self.assertEqual((R.reread(low, m, figure), self.kit(figure, interval)), ("consistent", "consistent"))
            high, interval = self.row(m, self.N, self.N)
            for figure in (INSIDE_HIGH, 1 - 3e-7, 1 - 9e-7):
                self.assertEqual((R.reread(high, m, figure), self.kit(figure, interval)), ("consistent", "consistent"))
            # The other end of the same row is computed, and the stop holds there.
            with self.assertRaises(SystemExit):
                R.reread(low, m, float(low[m + "_hi"]) + 3e-7)
            with self.assertRaises(SystemExit):
                R.reread(high, m, float(high[m + "_lo"]) - 3e-7)

    def test_a_printed_end_that_is_not_the_kit_s_for_the_counts_is_refused(self):
        for m in R.MEASURES:
            low, _ = self.row(m, 0, self.N)
            high, _ = self.row(m, self.N, self.N)
            for row, end, wrong in ((low, "_lo", "0.000001"), (high, "_hi", "0.999999")):
                # A report figure of one half is over the first row's interval and under the second's.
                self.assertEqual(R.reread(row, m, 0.5), "below" if row is low else "above")
                with self.assertRaises(SystemExit) as stop:
                    R.reread({**row, m + end: wrong}, m, 0.5)
                self.assertIn("not the end the kit sets", str(stop.exception))

    def test_the_counts_of_a_measure_say_which_ends_the_kit_set(self):
        self.assertEqual(R.exact_ends({"tp": "0", "fn": "7", "fp": "3", "tn": "4"}, "sensitivity"), (True, False))
        self.assertEqual(R.exact_ends({"tp": "7", "fn": "0", "fp": "3", "tn": "4"}, "sensitivity"), (False, True))
        self.assertEqual(R.exact_ends({"tp": "0", "fn": "7", "fp": "3", "tn": "4"}, "specificity"), (False, False))
        self.assertEqual(R.exact_ends({"tp": "2", "fn": "7", "fp": "3", "tn": "0"}, "specificity"), (True, False))
        self.assertEqual(R.exact_ends({"tp": "2", "fn": "7", "fp": "0", "tn": "4"}, "specificity"), (False, True))
        for n in (1, 25, 1000, self.N):
            self.assertEqual((R.score.clopper_pearson(0, n)[0], R.score.clopper_pearson(n, n)[1]), (0.0, 1.0))
            lo, hi = R.score.clopper_pearson(1, n + 1)[0], R.score.clopper_pearson(n, n + 1)[1]
            self.assertTrue(0.0 < lo and hi < 1.0, n)


class Tables(unittest.TestCase):
    def changed_figures(self, tmp):
        """A figure file where the first row with `called_metrics` has those counts as its figure."""
        with open(FIGURES, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        row = next(r for r in rows if r["called_n"] != "")
        for k in R.COUNTS:
            row[f"lineage_{k}"] = row[f"called_{k}"]
        row["lineage_sensitivity"] = repr(int(row["called_tp"]) / (int(row["called_tp"]) + int(row["called_fn"])))
        row["lineage_specificity"] = repr(int(row["called_tn"]) / (int(row["called_tn"]) + int(row["called_fp"])))
        path = pathlib.Path(tmp) / "figures.tsv"
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, R.COLUMNS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        return path

    def test_the_document_holds_the_tables_the_files_give(self):
        self.assertTrue(R.check(RESULTS, FIGURES, DOCUMENT))
        self.assertEqual(R.main(["revision.py", str(RESULTS), str(FIGURES), "--check", str(DOCUMENT)]), 0)

    def test_the_figure_file_is_read_back_as_written(self):
        """The csv reader and writer used here leave the committed file byte for byte."""
        with tempfile.TemporaryDirectory() as tmp:
            with open(FIGURES, newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            copy = pathlib.Path(tmp) / "copy.tsv"
            with open(copy, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, R.COLUMNS, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            self.assertEqual(copy.read_bytes(), FIGURES.read_bytes())

    def test_a_changed_figure_fails_the_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            changed = self.changed_figures(tmp)
            self.assertFalse(R.check(RESULTS, changed, DOCUMENT))
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                self.assertEqual(R.main(["revision.py", str(RESULTS), str(changed), "--check", str(DOCUMENT)]), 1)
            self.assertIn("differ", err.getvalue())
            # The row is no longer one that differs, so the tables list two lines fewer.
            lines = lambda path: [x for x in R.build(RESULTS, path).split("\n") if x.startswith("| ")]
            self.assertLess(len(lines(changed)), len(lines(FIGURES)))

    def test_a_changed_cell_in_the_document_fails_the_check(self):
        text = DOCUMENT.read_text(encoding="utf-8")
        begin, end = text.index(R.BEGIN), text.index(R.END)
        block = text[begin:end]
        self.assertIn(" consistent ", block)
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "changed.md"
            changed.write_text(text[:begin] + block.replace(" consistent ", " above ", 1) + text[end:], encoding="utf-8")
            self.assertFalse(R.check(RESULTS, FIGURES, changed))
            unmarked = pathlib.Path(tmp) / "unmarked.md"
            unmarked.write_text(text.replace(R.BEGIN, ""), encoding="utf-8")
            self.assertFalse(R.check(RESULTS, FIGURES, unmarked))

    def test_a_changed_reading_in_a_listed_row_fails_the_check(self):
        """The tables print the run's reading of each row whose figure differs."""
        _, figures = R.load_figures(FIGURES)
        later = R.later(T.load(RESULTS), figures)
        pop, r = next((pop, r) for pop, rows in later.items() for r in rows if r["differs"] and "below" in (r["reading_sensitivity"], r["reading_specificity"]))
        with tempfile.TemporaryDirectory() as tmp:
            copy = pathlib.Path(tmp) / "outside-study"
            shutil.copytree(RESULTS, copy)
            path = copy / pop / r["species"] / "metrics.tsv"
            lines = path.read_text(encoding="utf-8").split("\n")
            (i,) = [i for i, line in enumerate(lines) if line.startswith(f"{r['dataset']}\t{r['learner']}\t")]
            self.assertIn("\tbelow\t", lines[i])
            lines[i] = lines[i].replace("\tbelow\t", "\tabove\t", 1)
            path.write_text("\n".join(lines), encoding="utf-8")
            self.assertFalse(R.check(copy, FIGURES, DOCUMENT))

    def test_each_differing_row_is_listed_measure_by_measure(self):
        _, figures = R.load_figures(FIGURES)
        differing = [r for rows in R.later(T.load(RESULTS), figures).values() for r in rows if r["differs"]]
        block = R.build(RESULTS, FIGURES)
        detail = block[block.index("**Rows read whose report figure differs**") : block.index("**Genomes of null-fit folds")]
        self.assertEqual(len([x for x in detail.split("\n") if x.startswith("| ") and not x.startswith("| Population")]), 2 * len(differing))
        for r in differing:
            for m in R.MEASURES:
                self.assertIn(f"| {m} | {T.estimate(r, m)} | {R.figure_cell(r, m)} | {r['reading_' + m]} | {R.figure_cell(r['figure'], m)} | {r['later'][m]} |", detail)

    def test_the_genomes_of_null_fit_folds_are_the_difference_of_the_two_figures(self):
        _, figures = R.load_figures(FIGURES)
        later = R.later(T.load(RESULTS), figures)
        rows = R.null_folds(later)
        kept = {(r["dataset"], r["learner"]) for read in later.values() for r in read if r["figure"]["called_n"] != ""}
        self.assertEqual(sorted((r["dataset"], r["learner"]) for r, *_ in rows), sorted(kept))
        self.assertTrue(rows)
        block = R.build(RESULTS, FIGURES)
        table = block[block.index("**Genomes of null-fit folds") :]
        self.assertEqual(len([x for x in table.split("\n") if x.startswith("| *")]), len(rows))
        for r, n, null, resistant, susceptible in rows:
            f = r["figure"]
            every, fitted = ({k: int(f[f"{part}_{k}"]) for k in R.COUNTS} for part in ("lineage", "called"))
            self.assertEqual((n, null), (every["n"], every["n"] - fitted["n"]))
            self.assertGreater(null, 0)
            self.assertEqual(resistant, every["tp"] + every["fp"] - fitted["tp"] - fitted["fp"])
            self.assertEqual(susceptible, every["fn"] + every["tn"] - fitted["fn"] - fitted["tn"])
            self.assertEqual(resistant + susceptible, null)
            self.assertIn(
                f"| {T.title(r['species'])} | {T.drug(r['dataset'])} | {T.NAMES[r['learner']]} | {T.num(n)} | {T.num(null)} "
                f"| {T.pct(null / n)} | {T.num(resistant)} | {T.num(susceptible)} |",
                table,
            )
        shares = [null / n for _, n, null, _, _ in rows]
        self.assertEqual(shares, sorted(shares))

    def test_a_second_figure_that_is_no_part_of_the_first_is_refused(self):
        _, figures = R.load_figures(FIGURES)
        later = R.later(T.load(RESULTS), figures)
        pop, r = next((pop, r) for pop, read in later.items() for r in read if r["figure"]["called_n"] != "")
        f = r["figure"]
        self.assertEqual(len(R.null_folds({pop: [r]})), 1)
        for changed in ({"called_n": str(int(f["lineage_n"]) + 1)}, {"called_tp": str(int(f["lineage_tp"]) + 1)}, {"called_n": str(int(f["called_n"]) - 1)}):
            with self.assertRaises(SystemExit):
                R.null_folds({pop: [{**r, "figure": {**f, **changed}}]})

    def test_a_wrong_command_line_prints_the_usage(self):
        for args in ([], [str(RESULTS)], [str(RESULTS), "--figures"], [str(RESULTS), str(FIGURES), "--check"], [str(RESULTS), str(FIGURES), "--chek", str(DOCUMENT)]):
            with self.assertRaises(SystemExit) as stop:
                R.main(["revision.py", *args])
            self.assertEqual(stop.exception.code, R.__doc__)

    def test_a_figure_file_with_other_columns_or_a_repeated_row_is_refused(self):
        text = FIGURES.read_text(encoding="utf-8")
        lines = text.split("\n")
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "figures.tsv"
            for changed in (text.replace("report_sha256", "sha256", 1), "\n".join(lines[:2] + lines[1:]), "\n".join(lines[:2] + [lines[2].replace("\t14\t", "\t15\t", 1)] + lines[3:])):
                path.write_text(changed, encoding="utf-8")
                with self.assertRaises(SystemExit):
                    R.load_figures(path)


class FiguresFromReports(unittest.TestCase):
    """`--figures` on made-up reports. The numbers are not any dataset's."""

    DATASET = "shigella_flexneri__ciprofloxacin"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        base = pathlib.Path(self.tmp.name)
        self.root, self.reports = base / "outside-study", base / "datasets"
        shutil.copytree(RESULTS / "primary" / "shigella_flexneri", self.root / "primary" / "shigella_flexneri")
        self.reports.mkdir()

    def write(self, report):
        raw = json.dumps(report).encode()
        (self.reports / f"{self.DATASET}.json").write_bytes(raw)
        return hashlib.sha256(raw).hexdigest()

    def rows(self):
        text = R.figures(self.root, self.reports)
        return {r["learner"]: r for r in csv.DictReader(io.StringIO(text), delimiter="\t")}

    def report(self, lineage):
        return {
            "id": self.DATASET,
            "schema": R.score.REPORT_SCHEMA,
            "evaluation": {"revision": 99, "valid": True},
            "splits": [{"scheme": "random", "learners": []}, lineage],
        }

    def test_a_row_is_the_report_s_lineage_row_with_its_hash(self):
        learners = [
            {"learner": "scm", "metrics": metrics(60, 20, 10, 5, 25)},
            {"learner": "lasso", "metrics": metrics(60, 28, 2, 20, 10), "called_metrics": metrics(22, 8, 2, 2, 10)},
            {"learner": "catalog", "metrics": metrics(60, 30, 0, 30, 0)},
        ]
        digest = self.write(self.report({"scheme": "lineage", "learners": learners}))
        rows = self.rows()
        self.assertEqual(sorted(rows), ["cart", "known_gene", "lasso", "scm"])
        for r in rows.values():
            self.assertEqual((r["dataset"], r["report_sha256"], r["evaluation_revision"]), (self.DATASET, digest, "99"))
        self.assertEqual(R.counts(rows["scm"]), ("60", "20", "10", "5", "25"))
        self.assertEqual((rows["scm"]["lineage_sensitivity"], rows["scm"]["lineage_specificity"]), (repr(20 / 30), repr(25 / 30)))
        self.assertEqual([rows["scm"][f"called_{k}"] for k in R.COUNTS], [""] * 5)
        # The known-gene rule is the report's `catalog` row.
        self.assertEqual(R.counts(rows["known_gene"]), ("60", "30", "0", "30", "0"))
        self.assertEqual([rows["lasso"][f"called_{k}"] for k in R.COUNTS], ["22", "8", "2", "2", "10"])
        self.assertEqual(R.counts(rows["lasso"]), ("60", "28", "2", "20", "10"))
        # A learner with no row has the kit's note and no figure.
        self.assertEqual(rows["cart"]["lineage_note"], "no cart row in the lineage split")
        self.assertEqual(rows["cart"]["lineage_n"], "")
        # What it writes is a figure file this script reads.
        path = pathlib.Path(self.tmp.name) / "figures.tsv"
        path.write_text(R.figures(self.root, self.reports), encoding="utf-8")
        revision, figures = R.load_figures(path)
        self.assertEqual((revision, sorted(l for _, l in figures)), (99, ["cart", "known_gene", "lasso", "scm"]))

    def test_an_unavailable_lineage_split_is_noted_for_every_learner(self):
        self.write(self.report({"scheme": "lineage", "unavailable_reason": "made up"}))
        for r in self.rows().values():
            self.assertEqual((r["lineage_note"], r["lineage_n"], r["called_n"]), ("lineage split unavailable: made up", "", ""))

    def test_a_report_of_another_dataset_is_refused(self):
        self.write({**self.report({"scheme": "lineage", "learners": []}), "id": "shigella_flexneri__ampicillin"})
        with self.assertRaises(SystemExit):
            R.figures(self.root, self.reports)

    def test_a_report_of_another_schema_or_an_invalid_evaluation_is_refused(self):
        good = self.report({"scheme": "lineage", "learners": []})
        self.write(good)
        self.assertTrue(R.figures(self.root, self.reports))
        for changed in (
            {**good, "schema": R.score.REPORT_SCHEMA + 1},
            {k: v for k, v in good.items() if k != "schema"},
            {**good, "evaluation": {"revision": 99, "valid": False}},
            {**good, "evaluation": {"revision": 99}},
        ):
            self.write(changed)
            with self.assertRaises(SystemExit):
                R.figures(self.root, self.reports)

    def test_the_command_prints_the_figure_file(self):
        self.write(self.report({"scheme": "lineage", "learners": [{"learner": "scm", "metrics": metrics(60, 20, 10, 5, 25)}]}))
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(R.main(["revision.py", str(self.root), "--figures", str(self.reports)]), 0)
        self.assertEqual(out.getvalue(), R.figures(self.root, self.reports))


if __name__ == "__main__":
    unittest.main()
