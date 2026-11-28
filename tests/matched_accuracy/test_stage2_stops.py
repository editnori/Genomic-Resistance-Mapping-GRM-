"""The stops of scripts/matched-accuracy/stage2.py: one made-up input per check, each one printed line and exit 1.

Before the line, a stop writes nothing; after it, the line and one stop line, and no export. Every input here is a
SYNTHETIC release of made_up2.py, changed in one place."""

import json
import unittest
from pathlib import Path
from unittest import mock

from tests.matched_accuracy import made_up2, stage2_common as c
from tests.matched_accuracy.stage2_common import stage2

OTHER = "0" * 64


def first_line_of(run, reading, dataset):
    return next(i for i, line in enumerate(run.row_list(reading).read_text().splitlines())
                if line.startswith(dataset + "\t"))


def set_field(dataset, field, value, reading="A"):
    """An `after` that sets one field of the first listed row of a dataset in one row list."""
    def after(run):
        i = first_line_of(run, reading, dataset)

        def change(lines):
            f = lines[i][:-1].split("\t")
            f[stage2.ROW_FIELDS.index(field)] = value
            lines[i] = "\t".join(f) + "\n"
            return lines
        run.rewrite_list(reading, change)
    return after


class Sums(c.Case, unittest.TestCase):
    def test_a_report_that_differs_from_its_line_of_the_sums_list_stops(self):
        run = self.run_stage2(after=lambda r: r.rewrite_report("synthia_kappa__drug_a",
                                                                lambda d: d.update(seconds=3.5), sums=False))
        self.assert_stops(run, "differs from its line of the release's sums list")

    def test_a_row_list_with_another_sha256_stops(self):
        for reading in ("A", "B"):
            with self.subTest(reading=reading):
                run = self.run_stage2(pins=lambda own, reading=reading: {
                    "ROW_LIST_SHA256": dict(own["ROW_LIST_SHA256"], **{reading: OTHER})})
                self.assert_stops(run, f"the row list of reading {reading} differs from the sha256 stage 1 printed")

    def test_a_count_table_with_another_sha256_stops(self):
        run = self.run_stage2(pins={"COUNTS_SHA256": OTHER})
        self.assert_stops(run, "stage 1's committed count table differs from its sha256")

    def test_a_genome_table_with_another_sha256_stops(self):
        run = self.run_stage2(pins={"HC20_SHA256": OTHER})
        self.assert_stops(run, "the genome table, 5 Nov, differs from its sha256")

    def test_a_report_whose_id_is_not_its_file_name_stops(self):
        # Stage 1 ran on kappa's report; its file then takes another name, with its line of the sums list.
        def after(run):
            old = run.report("synthia_kappa__drug_a")
            old.rename(old.with_name("synthia_kappa__drug_z.json"))
            run.rewrite_sums()
        self.assert_stops(self.run_stage2(after=after),
                          "synthia_kappa__drug_z.json gives the id synthia_kappa__drug_a, not its file name")

    def test_a_report_left_out_whose_id_is_not_its_file_name_stops(self):
        # A made-up report of a ninth species holds none of the five learners, so stage 1 leaves it out. After stage
        # 1 ran, its file takes another name, with its line of the sums list.
        def before(release):
            rows = [{"id": made_up2.made_up.gid("ch", n), "resistant": n % 2 == 0, "fold": n % 5, "lineage": 30 + n,
                     "mlst": "ST-1", "from_mic": False} for n in range(1, 8)]
            release["reports"].append(made_up2.made_up.report("synthia_chi__drug_a", "Synthia chi", "synthia_chi",
                                                              "drug_a", rows, None, learners=()))

        def after(run):
            old = run.report("synthia_chi__drug_a")
            old.rename(old.with_name("synthia_chi__drug_z.json"))
            run.rewrite_sums()
        self.assert_stops(self.run_stage2(before=before, after=after),
                          "synthia_chi__drug_z.json gives the id synthia_chi__drug_a, not its file name")

    def test_a_missing_row_list_stops(self):
        def argv(a):
            i = a.index("--rows-b")
            return a[:i + 1] + [a[i + 1] + "-gone"] + a[i + 2:]
        run = self.run_stage2(argv=argv)
        self.assert_stops(run, "the row list of reading B is missing")


class Lines(c.Case, unittest.TestCase):
    def test_a_repeated_line_stops(self):
        run = self.run_stage2(after=lambda r: r.rewrite_list("A", lambda lines: lines[:3] + lines[2:]))
        self.assert_stops(run, "a line of the row list of reading A is repeated")

    def test_lines_out_of_stage_1s_sort_stop(self):
        run = self.run_stage2(after=lambda r: r.rewrite_list("B", lambda lines: [lines[1], lines[0]] + lines[2:]))
        self.assert_stops(run, "the lines of the row list of reading B are out of stage 1's sort")

    def test_a_genome_listed_twice_for_one_dataset_stops(self):
        def change(lines):
            f = lines[0][:-1].split("\t")
            f[4] = f[4] + "9"   # the same genome again, with another BioProject text after it in the sort
            return [lines[0], "\t".join(f) + "\n"] + lines[1:]
        run = self.run_stage2(after=lambda r: r.rewrite_list("A", change))
        self.assert_stops(run, "gives one genome of one dataset twice")

    def test_a_line_without_eight_fields_stops(self):
        # The last line, so that a space in place of its first tab keeps the sort.
        run = self.run_stage2(after=lambda r: r.rewrite_list("A", lambda lines: lines[:-1]
                                                             + [lines[-1].replace("\t", " ", 1)]))
        self.assert_stops(run, "a line of the row list of reading A does not hold 8 fields")

    def test_a_list_without_a_last_line_break_stops(self):
        run = self.run_stage2(after=lambda r: r.row_list("A").write_bytes(r.row_list("A").read_bytes()[:-1]))
        self.assert_stops(run, "the row list of reading A does not end with a line break")

    def test_an_empty_list_stops(self):
        run = self.run_stage2(after=lambda r: r.row_list("A").write_bytes(b""))
        self.assert_stops(run, "the row list of reading A holds no line")

    def test_a_fold_that_is_not_an_integer_stops(self):
        run = self.run_stage2(after=set_field("synthia_kappa__drug_a", "fold", "01"))
        self.assert_stops(run, "a fold or lineage cluster that is not an integer")


class Rows(c.Case, unittest.TestCase):
    def test_a_listed_row_its_report_does_not_hold_stops(self):
        def change(lines):
            f = lines[0][:-1].split("\t")
            f[1] = f[1] + "x"   # a made-up genome id no report holds
            return sorted(lines + ["\t".join(f) + "\n"])
        run = self.run_stage2(after=lambda r: r.rewrite_list("A", change))
        self.assert_stops(run, "reading A lists a row its report does not hold")

    def test_a_listed_dataset_no_report_gives_stops(self):
        def change(lines):
            f = lines[0][:-1].split("\t")
            f[0] = "synthia_absent__drug_a"
            return sorted(lines + ["\t".join(f) + "\n"])
        run = self.run_stage2(after=lambda r: r.rewrite_list("B", change))
        self.assert_stops(run, "reading B lists a dataset no report gives")

    def test_a_listed_row_whose_fold_differs_from_its_reports_stops(self):
        run = self.run_stage2(after=set_field("synthia_kappa__drug_a", "fold", "3"))
        self.assert_stops(run, "reading A lists a row whose fold differs from its report's")

    def test_a_listed_row_whose_lineage_cluster_differs_from_its_reports_stops(self):
        run = self.run_stage2(after=set_field("synthia_kappa__drug_a", "lineage", "77", reading="B"))
        self.assert_stops(run, "reading B lists a row whose lineage cluster differs from its report's")


class Learners(c.Case, unittest.TestCase):
    MESSAGE = "lists a report without scm, cart, rf, gbt, lasso in that order under the lineage split"

    def test_a_listed_report_with_the_five_in_another_order_stops(self):
        def change(doc):
            split = next(s for s in doc["splits"] if s["scheme"] == "lineage")
            split["learners"][0], split["learners"][1] = split["learners"][1], split["learners"][0]
        run = self.run_stage2(after=lambda r: r.rewrite_report("synthia_rho__drug_b", change))
        self.assert_stops(run, "reading A " + self.MESSAGE)

    def test_a_listed_report_without_the_five_stops(self):
        def change(doc):
            split = next(s for s in doc["splits"] if s["scheme"] == "lineage")
            split["learners"] = []
            split["unavailable_reason"] = "made-up"
            del split["test_fold"]
        run = self.run_stage2(after=lambda r: r.rewrite_report("synthia_nu__drug_a", change))
        self.assert_stops(run, "reading A " + self.MESSAGE)


class Groups(c.Case, unittest.TestCase):
    def test_a_group_with_one_lab_outcome_stops(self):
        # mu's group 2 is one R row and two S rows; drop its R row from reading A.
        def after(run):
            doc = json.loads(run.report("synthia_mu__drug_a").read_bytes())
            r_row = next(s["id"] for s in doc["samples"] if s["resistant"] and s["lineage"] == 11)
            run.rewrite_list("A", lambda lines: [x for x in lines if x.split("\t")[1] != r_row])
        run = self.run_stage2(after=after)
        self.assert_stops(run, "reading A: a group rebuilt from the row list holds one lab outcome")

    def test_counts_that_differ_from_stage_1s_count_table_at_any_level_stop(self):
        cases = [("dataset", "synthia_rho__drug_b", "A", "pairs"), ("species", "synthia_xi", "B", "groups"),
                 ("overall", "all", "A", "species"), ("dataset", "synthia_kappa__drug_a", "B", "datasets"),
                 ("species", "synthia_tau", "A", "rows"), ("overall", "all", "B", "genomes")]
        for level, name, reading, column in cases:
            with self.subTest(level=level, column=column):
                def after(run, level=level, name=name, reading=reading, column=column):
                    lines = run.counts().read_text().split("\n")
                    head = lines[0].split("\t")
                    for i, line in enumerate(lines[1:-1], 1):
                        f = line.split("\t")
                        if (f[0], f[1], f[2]) == (level, name, reading):
                            f[head.index(column)] = str(int(f[head.index(column)]) + 1)
                            lines[i] = "\t".join(f)
                    run.counts().write_text("\n".join(lines))
                run = self.run_stage2(after=after)
                self.assert_stops(run, f"reading {reading}: the counts rebuilt from the row list differ from stage 1's "
                                       f"count table at the {level} level")

    def test_a_count_table_without_a_listed_dataset_stops(self):
        def after(run):
            lines = run.counts().read_text().split("\n")
            run.counts().write_text("\n".join(x for x in lines if not x.startswith("dataset\tsynthia_tau__drug_a\t")))
        run = self.run_stage2(after=after)
        self.assert_stops(run, "reading A: the counts rebuilt from the row list differ from stage 1's count table at "
                               "the dataset level")


class Typing(c.Case, unittest.TestCase):
    """The support of the MLST read and of the MLST and HC20 read, rebuilt from each list, against stage 1's
    committed typing table."""

    def test_a_typing_table_with_another_sha256_stops(self):
        run = self.run_stage2(pins={"TYPING_TABLE_SHA256": OTHER})
        self.assert_stops(run, "stage 1's committed typing table differs from its sha256")

    def test_a_committed_typing_table_without_stage_1s_columns_stops(self):
        def after(run):
            text = run.typing_table().read_text()
            run.typing_table().write_text(text.replace("primary_pairs", "primary_other", 1))
        self.assert_stops(self.run_stage2(after=after), "stage 1's committed typing table does not have stage 1's "
                                                        "columns")

    def test_typing_support_that_differs_from_stage_1s_typing_table_at_any_level_stops(self):
        cases = [("dataset", "synthia_kappa__drug_a", "A", "mlst", "pairs"),
                 ("species", "synthia_rho", "B", "mlst_hc20", "mlst_groups"),
                 ("overall", "all", "A", "mlst_hc20", "rows_lost_from_primary"),
                 ("dataset", "synthia_omicron__drug_a", "B", "mlst", "lacking_value_of_eligible_rows"),
                 ("species", "synthia_lambda", "A", "mlst", "genomes"),
                 ("overall", "all", "B", "mlst", "primary_pairs")]
        for level, name, reading, step, column in cases:
            with self.subTest(level=level, step=step, column=column):
                def change(head, f, level=level, name=name, reading=reading, step=step, column=column):
                    if f[:4] == [level, name, reading, step]:
                        f[head.index(column)] = str(int(f[head.index(column)]) + 1)
                run = self.run_stage2(after=lambda r: r.rewrite_table(r.typing_table(), change))
                self.assert_stops(run, f"reading {reading}: the {step} support rebuilt from the row list differs "
                                       f"from stage 1's typing table at the {level} level")

    def test_a_typing_table_without_a_listed_dataset_stops(self):
        def after(run):
            lines = run.typing_table().read_text().split("\n")
            run.typing_table().write_text("\n".join(
                x for x in lines if not x.startswith("dataset\tsynthia_tau__drug_a\tA\tmlst\t")))
        self.assert_stops(self.run_stage2(after=after), "reading A: the mlst support rebuilt from the row list "
                                                        "differs from stage 1's typing table at the dataset level")

    def test_an_hc20_group_other_than_the_one_stage_1_read_stops(self):
        # The genome table given to stage 2, with its own sha256, moves one resistant row of kappa's first group
        # to an HC20 group of its own. The count table still agrees; the MLST and HC20 read does not.
        def after(run):
            doc = json.loads(run.report("synthia_kappa__drug_a").read_bytes())
            genome = next(s["id"] for s in doc["samples"] if s["resistant"] and s["lineage"] == 10)
            path = run.top / "genomes.jsonl"
            lines = []
            for line in path.read_text().splitlines():
                g = json.loads(line)
                if g.get("genome_id") == genome:
                    g["cgmlst_hc20"] = "hc-9"
                lines.append(json.dumps(g) + "\n")
            path.write_text("".join(lines))
        # Levels are checked overall first, as the count table is.
        self.assert_stops(self.run_stage2(after=after), "reading A: the mlst_hc20 support rebuilt from the row list "
                                                        "differs from stage 1's typing table at the overall level")

    def test_the_rows_lacking_a_value_among_the_rows_entering_the_groups_are_not_rebuilt(self):
        # That column counts rows outside every eligible group, which no row list holds.
        def change(head, f):
            if f[:4] == ["overall", "all", "A", "mlst"]:
                i = head.index("lacking_value_of_rows_entering_groups")
                f[i] = str(int(f[i]) + 7)
        self.assert_passes(self.run_stage2(after=lambda r: r.rewrite_table(r.typing_table(), change)))


class Setup(c.Case, unittest.TestCase):
    def test_numpy_other_than_1_26_4_stops_before_the_line(self):
        with mock.patch.object(stage2.numpy, "__version__", "1" + ".26" + ".5"):
            run = self.run_stage2()
        self.assert_stops(run, "numpy is not the version the resampling is written for")

    def test_an_output_directory_that_holds_the_line_of_an_earlier_run_stops(self):
        def after(run):
            run.out.mkdir()
            (run.out / stage2.OPENED_FILE).write_text(stage2.OPENED + "\n")
        run = self.run_stage2(after=after)
        self.assertEqual(run.code, 1, run.lines)
        self.assertEqual(run.lines, [c.STOP + "the output directory already holds export, local or the line of an "
                                              "earlier run"])
        self.assertFalse(run.wrote_export())

    def test_an_output_directory_holding_an_export_or_stage_1s_local_directory_stops(self):
        for name in ("export", "local"):
            with self.subTest(name=name):
                def after(run, name=name):
                    (run.out / name).mkdir(parents=True)
                run = self.run_stage2(after=after)
                self.assertEqual(run.code, 1, run.lines)
                self.assertEqual(run.lines, [c.STOP + "the output directory already holds export, local or the line "
                                                      "of an earlier run"])
                self.assertEqual(sorted(p.name for p in run.out.iterdir()), [name])
                self.assertEqual(list((run.out / name).iterdir()), [])

    def test_a_missing_genome_table_stops_in_one_line_today(self):
        # Today the missing table stops in one line before the line. The Python error comes only with the stop cut.
        def argv(a):
            i = a.index("--hc20")
            return a[:i + 1] + [a[i + 1] + "-gone"] + a[i + 2:]
        self.assert_stops(self.run_stage2(argv=argv), "the genome table, 5 Nov, is missing")

    def test_a_committed_count_table_without_stage_1s_columns_stops(self):
        def after(run):
            text = run.counts().read_text()
            run.counts().write_text(text.replace("out_intermediate", "out_other", 1))
        run = self.run_stage2(after=after)
        self.assert_stops(run, "stage 1's committed count table does not have stage 1's columns")


def kappa(change):
    """An `after` that changes kappa's report, with its line of the sums list."""
    return lambda run: run.rewrite_report("synthia_kappa__drug_a", change)


def listed_row(doc):
    """The index of a row of kappa's first group, a listed row in both readings."""
    return next(i for i, s in enumerate(doc["samples"]) if s["lineage"] == 10)


def fit_sizes(learner, value):
    def change(doc):
        split = next(s for s in doc["splits"] if s["scheme"] == "lineage")
        next(e for e in split["learners"] if e["learner"] == learner)["fold_complexity"] = value
    return change


def listed_indices(run, dataset):
    """The report indices of a dataset's rows listed in either row list, in the order of the report."""
    genomes = {line.split("\t")[1] for reading in ("A", "B")
               for line in run.row_list(reading).read_text(encoding="utf-8").splitlines()
               if line.split("\t")[0] == dataset}
    doc = json.loads(run.report(dataset).read_bytes())
    return [i for i, s in enumerate(doc["samples"]) if s["id"] in genomes]


def calls_at(dataset, value, pick):
    """An `after` that sets the catalogue call of the listed rows pick(indices) gives, and keeps the indices on the
    run as run.listed."""
    def after(run):
        run.listed = listed_indices(run, dataset)

        def change(doc):
            for i in pick(run.listed):
                doc["samples"][i]["catalog_resistant"] = value
        run.rewrite_report(dataset, change)
    return after


# Values a fit size by fold may not take, each not a finite number: one for each learner and fold in turn.
NOT_FINITE = (None, float("nan"), "3", True, float("inf"), 10 ** 400, [3.0])


class Values(c.Case, unittest.TestCase):
    """The values of listed reports that the second loader reads, other than scores, checked before the line by
    count_bad: a catalogue call, the fit sizes and the names."""

    def test_a_catalogue_call_that_is_not_true_false_or_null_stops(self):
        for value in ("yes", 1, 0, 1.0, [True]):
            with self.subTest(value=value):
                def change(doc, value=value):
                    doc["samples"][listed_row(doc)]["catalog_resistant"] = value
                self.assert_stops(self.run_stage2(after=kappa(change)), "synthia_kappa__drug_a: listed rows whose "
                                  "catalog_resistant is not true, false or null: 1")

    def test_a_catalogue_call_on_an_unlisted_row_is_not_read(self):
        # The fill row of fold 2 leaves under the one-outcome rule, so no list holds it.
        def change(doc):
            next(s for s in doc["samples"] if s["lineage"] == 92)["catalog_resistant"] = "yes"
        self.assert_passes(self.run_stage2(after=kappa(change)))

    def test_a_fit_size_list_that_is_not_five_finite_numbers_stops_and_names_the_fit_size(self):
        for value in ([3.0] * 4, [3.0] * 6, None, [3.0] * 4 + [None], [3.0] * 4 + [float("nan")],
                      [3.0] * 4 + [float("inf")], [3.0] * 4 + [True], [3.0] * 4 + ["3"], [3.0] * 4 + [10 ** 400]):
            with self.subTest(value=value):
                self.assert_stops(self.run_stage2(after=kappa(fit_sizes("gbt", value))), "synthia_kappa__drug_a: the "
                                  "fit size by fold (fold_complexity) of gbt is not a list of 5 finite numbers")

    def test_a_bad_catalogue_call_on_any_listed_row_stops_and_each_is_counted(self):
        # Kappa's listed rows, by their place among the listed rows: the first, one in the middle and the last,
        # each alone, then every listed row at once.
        d = "synthia_kappa__drug_a"
        for where, pick in (("first", lambda ix: ix[:1]), ("middle", lambda ix: ix[len(ix) // 2:len(ix) // 2 + 1]),
                            ("last", lambda ix: ix[-1:]), ("every", lambda ix: ix)):
            with self.subTest(row=where):
                run = self.run_stage2(after=calls_at(d, "yes", pick))
                self.assertEqual(len(run.listed), 9)
                n = len(pick(run.listed))
                self.assert_stops(run, f"{d}: listed rows whose catalog_resistant is not true, false or null: {n}")

    def test_a_bad_fit_size_in_any_of_the_five_learners_at_any_fold_stops_and_names_the_learner(self):
        # In tau's report, each of the five learners at each of the 5 folds, the first and the last included.
        d = "synthia_tau__drug_a"
        for j, learner in enumerate(stage2.FIVE):
            for fold in range(5):
                value = NOT_FINITE[(j + fold) % len(NOT_FINITE)]
                with self.subTest(learner=learner, fold=fold, value=repr(value)[:12]):
                    def change(doc, learner=learner, fold=fold, value=value):
                        split = next(s for s in doc["splits"] if s["scheme"] == "lineage")
                        entry = next(e for e in split["learners"] if e["learner"] == learner)
                        entry["fold_complexity"][fold] = value
                    run = self.run_stage2(after=lambda r, change=change: r.rewrite_report(d, change))
                    self.assert_stops(run, f"{d}: the fit size by fold (fold_complexity) of {learner} is not a list "
                                           f"of 5 finite numbers")

    def test_fit_sizes_of_five_integers_and_a_null_fit_pass(self):
        self.assert_passes(self.run_stage2(after=kappa(fit_sizes("lasso", [3, 0, 3, 3, 3]))))

    def test_a_dataset_name_with_a_full_stop_stops(self):
        def before(release):
            release["reports"][0]["id"] = "synthia_kappa__drug" + ".q"
        self.assert_stops(self.run_stage2(before=before), "synthia_kappa__drug.q: its dataset or species name holds "
                                                          "a full stop, which no export file may hold")

    def test_a_species_name_with_a_full_stop_stops(self):
        def before(release):
            release["reports"][1]["species_slug"] = "synthia" + ".lambda"
        self.assert_stops(self.run_stage2(before=before), "synthia_lambda__drug_a: its dataset or species name holds "
                                                          "a full stop, which no export file may hold")

    def test_a_listed_row_that_moved_after_stage_1s_loader_read_it_stops(self):
        # Stage 1's loader is given kappa's report with its first two rows in each other's place, as if the file
        # changed between that read and the next. count_bad reads the checked bytes and finds each id elsewhere.
        real = stage2.stage1.load_reports

        def swapped(files):
            entering, left_out = real(files)
            for rep in entering:
                if rep["id"] == "synthia_kappa__drug_a":
                    s, f = rep["samples"], rep["lineage_split"]["test_fold"]
                    s[0], s[1], f[0], f[1] = s[1], s[0], f[1], f[0]
            return entering, left_out
        with mock.patch.object(stage2.stage1, "load_reports", side_effect=swapped):
            run = self.run_stage2()
        self.assert_stops(run, "synthia_kappa__drug_a: listed rows that moved after they were checked: 2")

    def test_a_report_changed_after_its_check_stops_before_the_line(self):
        real = stage2.load_reports

        def changed(files):
            out = real(files)
            path = dict(files)["synthia_rho__drug_b"]
            path.write_bytes(path.read_bytes() + b" ")
            return out
        with mock.patch.object(stage2, "load_reports", side_effect=changed):
            run = self.run_stage2()
        self.assert_stops(run, "synthia_rho__drug_b.json changed after it was checked")

    def test_what_count_bad_returns_is_integers_alone(self):
        run = self.run_stage2(after=kappa(fit_sizes("gbt", [3.0] * 4 + [float("nan")])))
        doc = json.loads(run.report("synthia_kappa__drug_a").read_bytes())
        listed = [(s["id"], i) for i, s in enumerate(doc["samples"])]
        got = stage2.count_bad(run.report("synthia_kappa__drug_a"), c.sha(run.report("synthia_kappa__drug_a")
                                                                           .read_bytes()), listed)
        self.assertEqual(list(got), [k for k, _ in stage2.BAD])
        self.assertEqual(set(type(v) for v in got.values()), {int})
        self.assertEqual(got["fold_complexity of gbt"], 1)
        self.assertEqual(sum(got.values()), 1)


class AfterTheLine(c.Case, unittest.TestCase):
    """A score the second loader cannot use stops the run after the line, with no export."""

    def score(self, value):
        def change(doc):
            split = next(s for s in doc["splits"] if s["scheme"] == "lineage")
            gbt = next(e for e in split["learners"] if e["learner"] == "gbt")
            i = next(i for i, s in enumerate(doc["samples"]) if s["lineage"] == 10)
            gbt["scores"][i] = value
        return lambda run: run.rewrite_report("synthia_kappa__drug_a", change)

    def test_a_null_score_stops(self):
        self.assert_stops_after_the_line(self.run_stage2(after=self.score(None)),
                                         "a learner's score is null, true or false, or not a number")

    def test_a_boolean_score_stops(self):
        self.assert_stops_after_the_line(self.run_stage2(after=self.score(True)),
                                         "a learner's score is null, true or false, or not a number")

    def test_a_text_score_stops(self):
        self.assert_stops_after_the_line(self.run_stage2(after=self.score("0" + ".5")),
                                         "a learner's score is null, true or false, or not a number")

    def test_a_score_not_finite_stops(self):
        for value in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(value=value):
                self.assert_stops_after_the_line(self.run_stage2(after=self.score(value)),
                                                 "a learner's score is not finite")

    def test_a_listed_report_that_cannot_be_read_after_the_line_stops_in_one_line_naming_it(self):
        # Finding 4 of stage 2's second review: a removed report, or a directory in its place, once gave a
        # Python error after the line.
        real = stage2.open_scores
        for how, error in (("removed", "FileNotFoundError"), ("a directory", "IsADirectoryError")):
            def opening(out, how=how):
                opened = real(out)
                path = Path(out).parent / "in" / "reports" / "synthia_rho__drug_b.json"
                path.unlink()
                if how == "a directory":
                    path.mkdir()
                return opened
            with self.subTest(how=how), mock.patch.object(stage2, "open_scores", side_effect=opening):
                run = self.run_stage2()
                self.assertEqual(run.report("synthia_rho__drug_b").parent, run.out.parent / "in" / "reports")
                self.assert_stops_after_the_line(run, f"synthia_rho__drug_b.json could not be read after the line "
                                                      f"({error})")
                self.assertEqual(run.lines[1], c.STOP + "synthia_rho__drug_b.json could not be read after the "
                                                        f"line ({error})")


if __name__ == "__main__":
    unittest.main()
