"""scripts/c3-inputs/page_tables.py: the gates and the tables of docs/input-arms-result.md.

Every reading, list row, species, antibiotic, dataset, binary and build here
is SYNTHETIC, made up in tests/c3_inputs/reading.py: three made-up species,
seven made-up datasets, one of them in M. One reading is built and written
once for the suite; each refusal test changes a copy of it and names the
gate that must stop it, by the words that gate prints.
"""

import contextlib
import hashlib
import io
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import reading as R

page = R.page
score = R.score


def rows_of(folder, name):
    return page.load(Path(folder) / name)


class Built(unittest.TestCase):
    """One SYNTHETIC reading, written once."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.top = Path(cls.tmp.name)
        try:
            cls.reading, cls.lists = R.build(cls.top / "clean")
            cls.out = cls.top / "out"
            cls.files = page.write(cls.reading, cls.lists, cls.out)
        except BaseException:
            cls.tmp.cleanup()
            raise

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def copy(self):
        """A copy of the clean reading to change, and its lists."""
        top = Path(tempfile.mkdtemp(dir=self.top))
        shutil.copytree(self.reading, top / "reading")
        return top / "reading"

    def lrows(self):
        return page.load(self.lists)

    def stops(self, reading, words, stage=None):
        """check_reading (or `stage`) on `reading` stops with words that match `words`."""
        stage = stage or (lambda: page.check_reading(reading, self.lrows()))
        with self.assertRaisesRegex(page.Stop, words):
            stage()


class Tables(Built):
    def test_a_clean_reading_passes_every_gate_and_writes_its_files(self):
        want = sorted(["README.md", "gates.txt", "reading.txt", "count/count.tsv", "count/count-no-constant.tsv",
                       "count/conditions.tsv"] + list(page.TABLE_COLUMNS))
        self.assertEqual(self.files, want)
        self.assertEqual(len((self.out / "gates.txt").read_text().splitlines()), len(page.GATES))

    def test_the_count_files_are_the_scorers_copied_and_summed(self):
        for name in page.COUNT_FILES:
            self.assertEqual((self.out / "count" / name).read_bytes(), (self.reading / "count" / name).read_bytes())
        sums = (self.out / "SHA256SUMS").read_text().splitlines()
        self.assertEqual([line[66:] for line in sums], self.files)
        for line in sums:
            self.assertEqual(hashlib.sha256((self.out / line[66:]).read_bytes()).hexdigest(), line[:64])
        readme = (self.out / "README.md").read_text()
        self.assertIn("They are score.py's own tables", readme)
        self.assertIn("copied them here as they are, and sums them in `SHA256SUMS`", readme)

    def test_every_species_and_antibiotic_is_named_in_full(self):
        species = rows_of(self.out, "page-species.tsv")
        self.assertEqual([r["species"] for r in species], ["Fakeria maxima", "Mockella minor", "Synthia coli"])
        for r in rows_of(self.out, "page-datasets.tsv"):
            self.assertEqual((r["species"], r["antibiotic"]), R.NAMES[r["dataset"]])
        self.assertEqual({r["species"] for r in rows_of(self.out, "page-species-count.tsv")},
                         {"Fakeria maxima", "Mockella minor", "Synthia coli"})

    def test_each_arm_stands_on_its_own_datasets(self):
        arms = {r["arm"]: r for r in rows_of(self.out, "page-arms.tsv")}
        for arm in page.PROGRAM_ARMS:
            self.assertEqual((arms[arm]["set"], arms[arm]["datasets"], arms[arm]["species"]), ("A", "3", "2"))
            self.assertEqual(arms[arm]["study_datasets"], "1")
        for arm in ("k31", "unitigs", "k31-counts", "k31-aa", "twin"):
            self.assertEqual((arms[arm]["set"], arms[arm]["datasets"], arms[arm]["species"]), ("L", "6", "3"))
        program = {r["dataset"] for r in rows_of(self.out, "page-datasets.tsv") if r["pair"] == "program-k31"}
        self.assertEqual(program, set(R.IN_A))
        builds = rows_of(self.out, "page-builds.tsv")
        self.assertEqual(sorted((r["species_slug"], r["kind"]) for r in builds),
                         sorted((s, k) for s in R.SPECIES for k in ("aa", "counts", "unitigs")))

    def test_the_twin_stands_apart_and_is_never_an_arm(self):
        conds = rows_of(self.out, "page-conditions.tsv")
        self.assertEqual(sorted({r["arm"] for r in conds}), sorted(page.ARMS))
        self.assertNotIn("twin", {r["arm"] for r in conds})
        twin = [r for r in rows_of(self.out, "page-datasets.tsv") if r["pair"] == "twin-k31"]
        self.assertEqual({r["split"] for r in twin}, {"lineage"})
        self.assertEqual(len(twin), len(R.FIRST) * len(R.LEARNERS))
        arms = {r["arm"]: r for r in rows_of(self.out, "page-arms.tsv")}
        self.assertEqual(arms["twin"]["study_datasets"], "0")

    def test_a_study_row_stays_on_the_study_split(self):
        # K31 plus counts with the forest: +0.05 and -0.05 on the lineage split, -0.05 on drug a's study run.
        rows = {(r["dataset"], r["split"]): r for r in rows_of(self.out, "page-datasets.tsv")
                if r["pair"] == "k31-counts-k31" and r["learner"] == "rf"}
        self.assertEqual(rows[("synthia_coli__drug_a", "lineage")]["difference"], "0.05")
        self.assertEqual(rows[("synthia_coli__drug_a", "study")]["difference"], "-0.05")
        counts = {(r["species_slug"], r["split"]): r for r in rows_of(self.out, "page-species-count.tsv")
                  if r["pair"] == "k31-counts-k31" and r["learner"] == "rf"}
        lineage, study = counts[("synthia_coli", "lineage")], counts[("synthia_coli", "study")]
        self.assertEqual((lineage["datasets"], lineage["better"], lineage["worse"], lineage["median"]), ("2", "1", "1", "0.0"))
        self.assertEqual((study["datasets"], study["better"], study["worse"]), ("1", "0", "1"))

    def test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero(self):
        rows = [r for r in rows_of(self.out, "page-datasets.tsv")
                if r["dataset"] == "mockella_minor__drug_e" and r["split"] == "lineage" and r["pair"] != "twin-k31"]
        self.assertTrue(rows)
        for r in rows:
            self.assertEqual((r["arm_value"], r["default_value"], r["difference"]), ("", "", ""))
            self.assertEqual((r["arm_draws_scored"], r["draws"]), ("2", "3"))
        counts = [r for r in rows_of(self.out, "page-species-count.tsv")
                  if r["species_slug"] == "mockella_minor" and r["split"] == "lineage" and r["pair"] == "unitigs-k31"]
        self.assertEqual({(r["datasets"], r["entered"], r["not_scored"]) for r in counts}, {("2", "1", "1")})

    def test_a_constant_model_row_is_kept_and_marked(self):
        rows = [r for r in rows_of(self.out, "page-datasets.tsv")
                if r["dataset"] == "synthia_coli__drug_a" and r["split"] == "study" and r["pair"] == "program-k31"]
        self.assertEqual(len(rows), len(R.LEARNERS))
        for r in rows:
            self.assertEqual((r["constant_draws"], r["study_refused"], r["arm_value"]), ("1", "1", "0.5"))
            self.assertNotEqual(r["difference"], "")
        arms = {r["arm"]: r for r in rows_of(self.out, "page-arms.tsv")}
        self.assertEqual(arms["program"]["study_refused"], "1")
        self.assertEqual(arms["k31-program"]["study_refused"], "0")
        study = [r for r in rows_of(self.out, "page-count.tsv") if r["split"] == "study" and r["pair"] == "program-k31"
                 and r["measure"] == page.PRIMARY]
        kept = {r["file"]: r for r in study if r["learner"] == "rf"}
        self.assertEqual((kept["count.tsv"]["constant"], kept["count.tsv"]["entered"]), ("1", "1"))
        self.assertEqual((kept["count-no-constant.tsv"]["constant"], kept["count-no-constant.tsv"]["entered"]), ("0", "0"))

    def test_the_stricter_reading_decides_nothing(self):
        conds = rows_of(self.out, "page-conditions.tsv")
        self.assertEqual({(r["reading"], r["decides"]) for r in conds}, {("ruling 15", "1"), ("stricter", "0")})
        # Per arm, ruling 15 reads conditions "1 and 2", 3, 4, 5 and 6, and the stricter reading 1 to 6.
        self.assertEqual(len(conds), 11 * len(page.ARMS))
        self.assertEqual(sum(r["reading"] == "ruling 15" for r in conds), 5 * len(page.ARMS))
        count = rows_of(self.out, "page-count.tsv")
        self.assertEqual({(r["file"], r["decides"]) for r in count}, {("count.tsv", "1"), ("count-no-constant.tsv", "0")})

    def test_one_binary_is_named_with_its_code(self):
        readme = (self.out / "README.md").read_text()
        self.assertIn(f"One binary: every species' `binary.txt` gives the code sha256 `{R.CODE}`, 3 of 3 species.", readme)
        self.assertEqual({r["code"] for r in rows_of(self.out, "page-binaries.tsv")}, {R.CODE})

    def test_two_binaries_are_each_named_with_their_species(self):
        reading = self.copy()
        other = "f" * 64
        run = reading / "run" / "mockella_minor"
        for name in ("binary.txt", "DONE"):
            (run / name).write_text((run / name).read_text().replace(R.CODE, other))
        data = page.check_reading(reading, self.lrows())
        tables = page.tables_of(data)
        page.check_tables(data, tables)
        self.assertEqual(page.binaries_sentence(data),
                         f"2 binaries by `binary.txt`'s code sha256, each named with its species in `page-binaries.tsv`: "
                         f"`{R.CODE}`: Fakeria maxima, Synthia coli; `{other}`: Mockella minor.")
        self.assertEqual({(r["species_slug"], r["code"]) for r in tables["page-binaries.tsv"]},
                         {("fakeria_maxima", R.CODE), ("synthia_coli", R.CODE), ("mockella_minor", other)})

    def test_the_count_rows_carry_conditions_one_and_two_on_the_lineage_primary_rows(self):
        count = [r for r in rows_of(self.out, "page-count.tsv") if r["file"] == "count.tsv"]
        lineage = [r for r in count if r["split"] == "lineage" and r["measure"] == page.PRIMARY]
        self.assertEqual(len(lineage), 30)
        self.assertTrue(all(r["meets_1"] in ("0", "1") and r["meets_2"] in ("0", "1") for r in lineage))
        self.assertTrue(all(r["meets_1"] == "" for r in count if r not in lineage))
        # Unitigs call 4 more resistant genomes on every dataset: better on each of the 5 lineage datasets scored.
        unitigs = [r for r in lineage if r["pair"] == "unitigs-k31"]
        self.assertEqual({(r["better"], r["entered"], r["median_difference"]) for r in unitigs}, {("5", "5", "0.05")})

    def test_print_prints_every_table_from_the_files_alone(self):
        printed = page.tables(self.out)
        self.assertEqual(list(printed), ["conditions", "stricter", "count", "no-constant", "beside", "arms", "species",
                                         "builds", "binaries"])
        self.assertEqual(len(printed["conditions"]), 2 + 5 * len(page.ARMS))
        self.assertEqual(len(printed["stricter"]), 2 + 6 * len(page.ARMS))
        self.assertEqual(len(printed["count"]), 2 + 30)
        self.assertEqual(len(printed["species"]), 2 + 3)
        self.assertIn("| *Synthia coli* | 2 | 2 | 1 | 1 |", "\n".join(printed["species"]))
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(page.main(["print", str(self.out)]), 0)
        self.assertIn("[binaries]", out.getvalue())

    def test_write_refuses_an_out_folder_that_holds_files(self):
        with self.assertRaisesRegex(page.Stop, "is not empty"):
            page.write(self.reading, self.lists, self.out)

    def test_no_file_holds_a_sequence(self):
        for name in self.files:
            self.assertNotRegex((self.out / name).read_text(), r"[ACGT]{31}")


class Gates(Built):
    def test_a_reading_without_its_commands_is_refused(self):
        reading = self.copy()
        (reading / "reading.txt").unlink()
        self.stops(reading, "gate 1 reading: no reading.txt")

    def test_a_folder_outside_the_reading_is_refused(self):
        reading = self.copy()
        (reading / "pilot").mkdir()
        self.stops(reading, "gate 1 reading: folders that are no part of the reading")

    def test_a_species_missing_from_the_reading_is_refused(self):
        reading = self.copy()
        shutil.rmtree(reading / "score" / "fakeria_maxima")
        self.stops(reading, r"gate 2 species: score/ has 2 species of 3; missing \['fakeria_maxima'\]")

    def test_a_species_twice_is_refused(self):
        reading = self.copy()
        shutil.copytree(reading / "score" / "fakeria_maxima", reading / "score" / "fakeria_maxima_again")
        self.stops(reading, "gate 2 species: score/ has 4 species of 3")

    def test_a_species_whose_tables_hold_another_species_is_refused(self):
        reading = self.copy()
        rows = rows_of(reading / "score" / "synthia_coli", "draws.tsv")
        more = rows_of(reading / "score" / "fakeria_maxima", "draws.tsv")
        score.write_tsv(reading / "score" / "synthia_coli" / "draws.tsv", rows + more, score.DRAW_COLUMNS)
        self.stops(reading, "gate 2 species: score/synthia_coli/draws.tsv holds rows of another species")

    def test_a_twin_draw_at_a_draw_seed_is_refused(self):
        reading = self.copy()

        def change(rows):
            for r in rows:
                if r["arm"] == "twin" and r["seed"] == "45":
                    r["seed"] = "42"
            return rows
        R.rewrite_species(reading, "synthia_coli", change)
        self.stops(reading, "gate 3 twin: score/synthia_coli/draws.tsv: a twin row of .* is no lineage draw at a twin seed")

    def test_conditions_that_read_the_twin_as_an_arm_are_refused(self):
        reading = self.copy()
        conds = rows_of(reading / "count", "conditions.tsv")
        conds += [dict(r, arm="twin") for r in conds if r["arm"] == "unitigs"]
        score.write_tsv(reading / "count" / "conditions.tsv", conds, score.CONDITION_COLUMNS)
        self.stops(reading, "gate 3 twin: conditions.tsv does not read the five arms alone")

    def test_a_study_run_on_the_lineage_split_is_refused(self):
        reading = self.copy()

        def change(rows):
            for r in rows:
                if r["role"] == "study" and r["dataset"] == "synthia_coli__drug_a" and r["arm"] == "unitigs":
                    r["split"] = "lineage"
            return rows
        R.rewrite_species(reading, "synthia_coli", change)
        self.stops(reading, "gate 4 study: score/synthia_coli/draws.tsv: a study row of synthia_coli__drug_a stands on "
                            "the lineage split")

    def test_a_dataset_of_m_is_refused(self):
        reading = self.copy()
        rows = rows_of(reading / "score" / "synthia_coli", "draws.tsv")
        more = [dict(r, dataset="synthia_coli__drug_m", antibiotic="drug m") for r in rows
                if r["dataset"] == "synthia_coli__drug_a"]
        R.rewrite_species(reading, "synthia_coli", lambda rs: rs + more)
        self.stops(reading, "gate 5 M: score/synthia_coli/draws.tsv holds 1 datasets of M")

    def test_a_dataset_named_otherwise_than_in_the_lists_is_refused(self):
        reading = self.copy()

        def change(rows):
            for r in rows:
                if r["dataset"] == "fakeria_maxima__drug_c":
                    r["antibiotic"] = "drug c sodium"
            return rows
        R.rewrite_species(reading, "fakeria_maxima", change)
        self.stops(reading, "gate 5 M: score/fakeria_maxima/draws.tsv: fakeria_maxima__drug_c is named otherwise")

    def test_a_program_arm_outside_a_is_refused(self):
        reading = self.copy()
        for slug in R.SPECIES:
            def change(rows, slug=slug):
                more = []
                for d in R.NOT_A:
                    if d.startswith(slug + "__"):
                        for arm in ("program", "k31-program"):
                            more += [dict(r, arm=arm) for r in rows if r["dataset"] == d and r["arm"] == "unitigs"]
                return rows + more
            R.rewrite_species(reading, slug, change)
        self.stops(reading, "gate 6 own datasets: score/fakeria_maxima/datasets.tsv holds other rows than its arms' "
                            "own datasets")

    def test_a_matrix_arm_missing_on_a_dataset_is_refused(self):
        reading = self.copy()
        R.rewrite_species(reading, "mockella_minor", lambda rows: [
            r for r in rows if not (r["dataset"] == "mockella_minor__drug_f" and r["arm"] == "k31-aa")])
        self.stops(reading, "gate 6 own datasets: score/mockella_minor/datasets.tsv holds other rows")

    def test_a_value_for_a_dataset_with_a_draw_not_scored_is_refused(self):
        reading = self.copy()
        path = reading / "score" / "mockella_minor" / "datasets.tsv"
        rows = rows_of(path.parent, path.name)
        for r in rows:
            if r["dataset"] == "mockella_minor__drug_e" and r["split"] == "lineage" and r["arm"] == "unitigs":
                r["scored"], r[page.PRIMARY] = "1", "0.0"
        score.write_tsv(path, rows, score.DATASET_COLUMNS)
        self.stops(reading, r"gate 7 draws: score/mockella_minor: mockella_minor__drug_e lineage unitigs \w+: a draw not "
                            "scored is read as a value")

    def test_a_dataset_figure_that_is_not_its_draws_mean_is_refused(self):
        reading = self.copy()
        path = reading / "score" / "fakeria_maxima" / "datasets.tsv"
        rows = rows_of(path.parent, path.name)
        rows[3]["recall"] = "0.123"
        score.write_tsv(path, rows, score.DATASET_COLUMNS)
        self.stops(reading, "gate 7 draws: score/fakeria_maxima: .* recall is not the mean of its draws")

    def test_a_constant_model_draw_off_an_arms_study_run_is_refused(self):
        reading = self.copy()

        def change(rows):
            for r in rows:
                if r["dataset"] == "synthia_coli__drug_b" and r["arm"] == "unitigs" and r["seed"] == "43":
                    r["constant"], r["why"] = "1", page.CONSTANT_WHY
            return rows
        R.rewrite_species(reading, "synthia_coli", change)
        self.stops(reading, "gate 8 constant: score/synthia_coli/draws.tsv: synthia_coli__drug_b unitigs lineage: a "
                            "constant-model draw that is no arm's scored study run")

    def test_a_constant_model_draw_without_its_reason_is_refused(self):
        reading = self.copy()

        def change(rows):
            for r in rows:
                if r["constant"] == "1":
                    r["why"] = ""
            return rows
        R.rewrite_species(reading, "synthia_coli", change)
        self.stops(reading, "gate 8 constant: .* a constant-model draw and its reason do not go together")

    def test_a_failed_check_is_refused(self):
        reading = self.copy()
        path = reading / "score" / "mockella_minor" / "checks.tsv"
        rows = rows_of(path.parent, path.name)
        rows[0]["passed"] = "0"
        score.write_tsv(path, rows, ["file", "check", "passed", "detail"])
        self.stops(reading, "gate 9 checks: score/mockella_minor/checks.tsv: 1 of 2 checks failed")

    def test_a_build_record_beyond_its_counts_is_refused(self):
        reading = self.copy()
        path = reading / "run" / "fakeria_maxima" / "builds" / "unitigs.build.json"
        rec = json.loads(path.read_text())
        rec["genomes_held"] = ["made-up genome"]
        path.write_text(json.dumps(rec))
        self.stops(reading, "gate 10 builds: run/fakeria_maxima/builds/unitigs.build.json holds more or less than its "
                            "counts")

    def test_a_missing_build_is_refused(self):
        reading = self.copy()
        (reading / "run" / "synthia_coli" / "builds" / "aa.build.json").unlink()
        self.stops(reading, "gate 10 builds: run/synthia_coli/builds holds")

    def test_a_done_line_naming_another_code_is_refused(self):
        reading = self.copy()
        path = reading / "run" / "synthia_coli" / "DONE"
        path.write_text(path.read_text().replace(R.CODE, "f" * 64))
        self.stops(reading, "gate 11 binary: run/synthia_coli/DONE does not name synthia_coli and binary.txt's binary "
                            "and code")

    def test_a_page_that_leaves_a_species_binary_out_is_refused(self):
        data = page.check_reading(self.reading, self.lrows())
        tables = page.tables_of(data)
        tables["page-binaries.tsv"] = [r for r in tables["page-binaries.tsv"] if r["species_slug"] != "mockella_minor"]
        self.stops(None, "gate 11 binary: page-binaries.tsv does not name each species with its binary",
                   lambda: page.check_tables(data, tables))

    def test_a_page_that_writes_the_stricter_reading_as_deciding_is_refused(self):
        data = page.check_reading(self.reading, self.lrows())
        tables = page.tables_of(data)
        for r in tables["page-conditions.tsv"]:
            r["decides"] = 1
        self.stops(None, "gate 12 conditions: page-conditions.tsv: .*the stricter reading's condition 1 is written as "
                         "deciding", lambda: page.check_tables(data, tables))

    def test_conditions_that_differ_from_the_count_are_refused(self):
        reading = self.copy()
        conds = rows_of(reading / "count", "conditions.tsv")
        for r in conds:
            if r["arm"] == "unitigs" and r["reading"] == "stricter" and r["condition"] == "1":
                r["learners"] = str(int(r["learners"]) + 1)
        score.write_tsv(reading / "count" / "conditions.tsv", conds, score.CONDITION_COLUMNS)
        data = page.check_reading(reading, self.lrows())
        self.stops(None, "gate 12 conditions: unitigs: the stricter reading's condition 1 read here off count.tsv",
                   lambda: page.check_tables(data, page.tables_of(data)))

    def test_a_count_that_is_not_the_scorers_is_refused(self):
        reading = self.copy()
        count = rows_of(reading / "count", "count.tsv")
        for r in count:
            if r["split"] == "lineage" and r["pair"] == "k31-aa-k31" and r["learner"] == "lasso":
                r["better"] = str(int(r["better"]) + 1)
        score.write_tsv(reading / "count" / "count.tsv", count, score.COUNT_COLUMNS)
        data = page.check_reading(reading, self.lrows())
        self.stops(None, "gate 13 count: count/count.tsv is not what score.py --count writes",
                   lambda: page.check_count(data))


class GatesHeld(Built):
    """Gates 4, 7, 10 and 11 on what each must refuse, where no test held them (SYNTHETIC)."""

    def test_a_study_run_at_a_draw_seed_is_refused(self):
        reading = self.copy()

        def change(rows):
            for r in rows:
                if r["role"] == "study" and r["dataset"] == "synthia_coli__drug_a" and r["arm"] == "unitigs":
                    r["seed"] = "43"
            return rows
        R.rewrite_species(reading, "synthia_coli", change)
        self.stops(reading, "gate 4 study: score/synthia_coli/draws.tsv: a study run of synthia_coli__drug_a at seed 43")

    def test_a_dataset_figure_a_thousandth_off_its_draws_mean_is_refused(self):
        reading = self.copy()
        path = reading / "score" / "fakeria_maxima" / "datasets.tsv"
        rows = rows_of(path.parent, path.name)
        row = next(r for r in rows if r["recall"] != "")
        row["recall"] = repr(float(row["recall"]) + 0.001)
        score.write_tsv(path, rows, score.DATASET_COLUMNS)
        self.stops(reading, "gate 7 draws: score/fakeria_maxima: .* recall is not the mean of its draws")

    def test_a_build_record_of_another_species_is_refused(self):
        reading = self.copy()
        path = reading / "run" / "synthia_coli" / "builds" / "unitigs.build.json"
        rec = json.loads(path.read_text())
        rec["species"] = "fakeria_maxima"
        path.write_text(json.dumps(rec))
        self.stops(reading, "gate 10 builds: run/synthia_coli/builds/unitigs.build.json is not the unitigs build of "
                            "synthia_coli")

    def test_a_done_line_naming_another_species_is_refused(self):
        reading = self.copy()
        path = reading / "run" / "synthia_coli" / "DONE"
        path.write_text(path.read_text().replace("DONE synthia_coli ", "DONE fakeria_maxima "))
        self.stops(reading, "gate 11 binary: run/synthia_coli/DONE does not name synthia_coli")


class Withdrawn(unittest.TestCase):
    def test_a_withdrawn_arm_has_no_rows_no_build_and_its_five_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            reading, lists = R.build(Path(tmp) / "w", withdrawn=("k31-aa",))
            out = Path(tmp) / "out"
            page.write(reading, lists, out)
            arms = {r["arm"]: r for r in rows_of(out, "page-arms.tsv")}
            self.assertEqual((arms["k31-aa"]["withdrawn"], arms["k31-aa"]["datasets"]), ("1", "0"))
            self.assertFalse([r for r in rows_of(out, "page-datasets.tsv") if r["pair"] == "k31-aa-k31"])
            self.assertFalse([r for r in rows_of(out, "page-builds.tsv") if r["kind"] == "aa"])
            mine = [r for r in rows_of(out, "page-count.tsv") if r["file"] == "count.tsv" and r["pair"] == "k31-aa-k31"
                    and r["split"] == "lineage" and r["measure"] == page.PRIMARY]
            self.assertEqual({(r["withdrawn"], r["sign_p"]) for r in mine}, {("1", "1.0")})


if __name__ == "__main__":
    unittest.main()
