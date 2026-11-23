"""scripts/c2-cap/check_folds.py and study.py: held-out genomes stay out of every wider arm's training (docs/cap-arms.md).

Every genome, lineage, BioProject, fold plan, cohort plan, training trace
and fold file here is SYNTHETIC, in the shapes grm arm and grm learn write.
The study split is proved on these SYNTHETIC fold plans only.
"""

import json
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

folds = common.load("check_folds")
study = common.load("study")

D = f"{common.SPECIES}__drug_a"
DRAWS, TWIN = 2, 2


class Lineage(unittest.TestCase):
    def run_check(self, breach=None, edit=None, none_added=0, keyed=None, refused=()):
        """(exit code, printed text, rows written) of check_folds.py lineage on a SYNTHETIC run."""
        sp = common.Species()
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            sp.structure(t / "structure.json")
            sp.cohorts(t / "cohorts")
            common.write_json(t / "roots.json", {"fits": {"cap-all": "cap-all"}, "same_as": {"cap-4000": "cap-all"},
                                                 "genomes": {"cap-2000": 60, "cap-4000": 90, "cap-all": 90}})
            common.lineage_fits(sp, t / "fits", D, DRAWS + TWIN, DRAWS, breach=breach, keyed=keyed,
                                refused=refused)
            if edit:
                edit(t / "fits")
            (t / "added-labelled.tsv").write_text(
                "dataset\tspecies\tantibiotic\tarm\tlabelled\tadded\tadded_labelled\tadded_resistant\t"
                f"added_susceptible\tnone_added\n{D}\tSynthia coli\tdrug a\tcap-all\t90\t30\t"
                f"{0 if none_added else 30}\t0\t{0 if none_added else 30}\t{int(bool(none_added))}\n", encoding="utf-8")
            rc, text = common.quiet(folds.main, ["check_folds.py", "lineage", str(t / "fits"), D, str(t / "cohorts"),
                                                 str(t / "structure.json"), str(t / "roots.json"), str(DRAWS),
                                                 str(TWIN), str(t / "out.tsv"), str(t / "added-labelled.tsv")])
            rows = (t / "out.tsv").read_text(encoding="utf-8").splitlines()
        return rc, text, rows

    def test_a_correct_run_passes_and_counts_every_fold(self):
        rc, text, rows = self.run_check()
        self.assertEqual(rc, 0, text)
        self.assertNotIn("breach", text)
        self.assertEqual(rows[0].split("\t"), folds.COLUMNS)
        body = [r.split("\t") for r in rows[1:]]
        self.assertEqual(len(body), common.K * (DRAWS + TWIN) + common.K * DRAWS)
        arm = [r for r in body if r[1] == "cap-all"]
        # Each fold keeps out the added genomes of the lineages it scores.
        self.assertTrue(all(int(r[8]) > 0 for r in arm))
        self.assertEqual(sum(int(r[5]) for r in arm), 60 * DRAWS)

    def test_a_replaced_draw_is_read_at_its_slot_not_its_keyed_repeat(self):
        # grm arm replaces a refused default lineage draw by the next keyed repeat past the slots: the plan and
        # the cohort plan hold that repeat, the trace line the slot (src/ablate/arm.rs, arm_plans).
        rc, text, rows = self.run_check(keyed={1: DRAWS + TWIN, 3: DRAWS + TWIN + 1})
        self.assertEqual(rc, 0, text)
        self.assertNotIn("breach", text)
        self.assertEqual(len(rows) - 1, common.K * (DRAWS + TWIN) + common.K * DRAWS)

    def test_a_cohort_plan_of_another_keyed_repeat_breaches(self):
        def edit(fits):
            p = fits / "cap-all" / D / "plans" / "lineage-1.cohort.json"
            doc = json.loads(p.read_text(encoding="utf-8"))
            doc["repeat"] = DRAWS + TWIN
            p.write_text(json.dumps(doc), encoding="utf-8")
        rc, text, _ = self.run_check(edit=edit)
        self.assertEqual(rc, 1)
        self.assertIn("breach: cap-all: a cohort plan made from another plan than the default's slot: 1", text)

    def test_an_unavailable_lineage_split_is_not_produced_and_does_not_stop_the_run(self):
        # grm arm refused the default's slot-1 draw past its two replacements: every lineage run of the dataset is
        # unavailable and not produced, and the run goes on (docs/input-comparison-runner.md, A refused repeat).
        rc, text, rows = self.run_check(refused={1})
        self.assertEqual(rc, 0, text)
        self.assertIn("lineage split unavailable: 1 default slot plans refused", text)
        self.assertNotIn("breach", text)
        body = [r.split("\t") for r in rows[1:]]
        # cap-all traced slot 0 before it reached the refused slot; those folds are still checked.
        self.assertEqual(sorted({(r[1], r[3]) for r in body}), [("cap-all", "0")])
        self.assertEqual(len(body), common.K)

    def test_a_cohort_plan_or_a_default_trace_line_of_an_unavailable_split_breaches(self):
        def edit(fits):
            common.write_json(fits / "cap-all" / D / "plans" / "lineage-1.cohort.json", {"repeat": 1})
            p = fits / "cap-2000" / D / "training.jsonl"
            line = json.dumps({"split": "lineage", "repeat": 0, "fold": 0, "train": []})
            p.write_text(p.read_text(encoding="utf-8") + line + "\n", encoding="utf-8")
        rc, text, _ = self.run_check(refused={1}, edit=edit)
        self.assertEqual(rc, 1)
        self.assertIn("breach: cap-all: a cohort plan made from a refused default plan: 1", text)
        self.assertIn("breach: a trace line of an unavailable lineage split: 1", text)

    def test_an_added_genome_of_a_held_out_lineage_in_training_breaches(self):
        rc, text, _ = self.run_check("lineage")
        self.assertEqual(rc, 1)
        self.assertIn("breach: an added genome of a held-out lineage in training: 1", text)

    def test_an_outer_test_fold_that_differs_between_arms_breaches(self):
        rc, text, _ = self.run_check("differs")
        self.assertEqual(rc, 1)
        self.assertIn("breach: an outer test fold that differs between arms: 1", text)

    def test_an_added_genome_in_an_outer_test_fold_breaches(self):
        rc, text, _ = self.run_check("added-scored")
        self.assertEqual(rc, 1)
        self.assertIn("breach: an added genome in an outer test fold: 1", text)

    def test_an_added_copy_of_a_scored_genome_in_training_breaches(self):
        rc, text, _ = self.run_check("copy")
        self.assertEqual(rc, 1)
        self.assertIn("breach: an added copy of a scored genome in training: 1", text)

    def test_a_cohort_plan_made_from_another_default_plan_breaches(self):
        def edit(fits):
            p = fits / "cap-all" / D / "plans" / "lineage-1.cohort.json"
            doc = json.loads(p.read_text(encoding="utf-8"))
            doc["reference_sha256"] = "plan0"
            p.write_text(json.dumps(doc), encoding="utf-8")
        rc, text, _ = self.run_check(edit=edit)
        self.assertEqual(rc, 1)
        self.assertIn("a cohort plan made from another plan than the default's slot", text)

    def test_a_missing_slot_breaches(self):
        rc, text, _ = self.run_check(edit=lambda fits: (fits / "cap-2000" / D / "plans" / "lineage-3.json").unlink())
        self.assertEqual(rc, 1)
        self.assertIn("breach: a default slot with no plan: 1", text)
        rc, text, _ = self.run_check(edit=lambda fits: (fits / "cap-all" / D / "plans" / "lineage-0.cohort.json").unlink())
        self.assertEqual(rc, 1)
        self.assertIn("breach: cap-all: a slot with no cohort plan: 1", text)

    def test_a_balance_drop_breaches(self):
        def edit(fits):
            p = fits / "cap-all" / D / "plans" / "lineage-0.cohort.json"
            doc = json.loads(p.read_text(encoding="utf-8"))
            doc["balance_dropped"][2] = 1
            p.write_text(json.dumps(doc), encoding="utf-8")
        rc, text, _ = self.run_check(edit=edit)
        self.assertEqual(rc, 1)
        self.assertIn("breach: a balance drop: 1", text)


    def test_a_none_added_arm_is_not_checked_and_a_fit_of_it_breaches(self):
        rc, text, _ = self.run_check(none_added=1)
        self.assertEqual(rc, 1)
        self.assertIn("breach: cap-all: a none-added arm fitted: 1", text)
        rc, text, rows = self.run_check(none_added=1, edit=lambda fits: __import__("shutil").rmtree(fits / "cap-all"))
        self.assertEqual(rc, 0, text)
        self.assertEqual({r.split("\t")[1] for r in rows[1:]}, {"cap-2000"})


class Study(unittest.TestCase):
    def setUp(self):
        self.sp = common.Species()
        self.f0 = common.fold_file(self.sp, D)
        self.bp = common.bioprojects(self.sp, self.f0)

    def test_place_trains_only_genomes_of_new_bioprojects_that_copy_no_tested_genome(self):
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        self.assertEqual(sorted(train), sorted(self.sp.added[8:10] + self.sp.added[12:20]))
        self.assertEqual((len(unused["bioproject"]), len(unused["unknown"]), len(unused["duplicate"])), (8, 10, 2))

    def test_an_added_genome_of_a_held_out_bioproject_stays_out_of_training(self):
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        self.assertFalse({g for g in train if self.bp[g] == "PRJ0"})
        self.assertEqual(sorted(unused["bioproject"]), sorted(self.sp.added[:8]))

    def test_bioprojects_join_several_files_and_a_genome_they_disagree_on_is_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            common.meta(t / "a.jsonl.gz", {"1.1": "PRJA", "2.1": "PRJB"})
            common.meta(t / "b.jsonl.gz", {"2.1": "PRJB", "3.1": "PRJC", "1.1": "PRJZ"})
            got = folds.bioprojects(f"{t / 'a.jsonl.gz'},{t / 'b.jsonl.gz'}")
        self.assertEqual(got, {"2.1": "PRJB", "3.1": "PRJC"})

    def test_an_added_genome_of_no_known_bioproject_stays_out_of_training(self):
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        self.assertFalse(set(train) & set(self.sp.added[20:]))
        self.assertEqual(sorted(unused["unknown"]), sorted(self.sp.added[20:]))

    def write_run(self, t, arm_file=None):
        """plans/, cohorts/, the meta rows and the grouping of a SYNTHETIC study run; the arm's fold file from
        study.py unless arm_file is given."""
        self.sp.structure(t / "structure.json")
        self.sp.cohorts(t / "cohorts")
        common.meta(t / "meta.jsonl.gz", self.bp)
        common.write_ids(t / "plans" / "study.txt", [D])
        common.write_json(t / "plans" / "reviewed" / "folds" / f"{D}.study-whole.json", self.f0)
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        common.write_json(t / "plans" / "cap-all" / f"{D}.study-whole.json",
                          arm_file if arm_file is not None else study.fold_file(self.f0, train, unused))

    def check(self, t):
        return common.quiet(folds.main, ["check_folds.py", "study", str(t / "plans"), str(t / "cohorts"),
                                         str(t / "meta.jsonl.gz"), str(t / "structure.json"), "cap-4000 cap-all",
                                         str(t / "out.tsv")])

    def test_study_files_keep_the_defaults_test_folds_and_pass_the_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            self.write_run(t)
            rc, text = self.check(t)
            rows = (t / "out.tsv").read_text(encoding="utf-8").splitlines()
        self.assertEqual(rc, 0, text)
        self.assertEqual(rows[1].split("\t"), [D, "cap-all", "5", str(sum(map(len, self.f0["test"]))), "0", "10", "20"])

    def test_the_study_check_breaches_on_a_held_out_bioproject_in_training(self):
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        bad = study.fold_file(self.f0, train + [self.sp.added[0]], unused)
        bad["unused"].remove(self.sp.added[0])
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            self.write_run(t, bad)
            rc, text = self.check(t)
        self.assertEqual(rc, 1)
        self.assertIn("breach: an added genome of a held-out BioProject in training: 1", text)

    def test_the_study_check_breaches_on_other_test_folds_and_an_added_tested_genome(self):
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        bad = study.fold_file(self.f0, train, unused)
        bad["test"] = [list(f) for f in bad["test"]]
        bad["test"][1].append(train[0])
        bad["train"].remove(train[0])
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            self.write_run(t, bad)
            rc, text = self.check(t)
        self.assertEqual(rc, 1)
        self.assertIn("breach: an outer test fold that differs between arms: 1", text)
        self.assertIn("breach: an added genome in an outer test fold: 1", text)

    def test_the_study_check_breaches_on_an_unknown_bioproject_a_copy_and_a_moved_default_genome(self):
        train, unused = study.place(self.f0, set(self.sp.added), self.bp, self.sp.digest)
        bad = study.fold_file(self.f0, train + [self.sp.added[25], self.sp.added[10]], unused)
        bad["unused"] = [g for g in bad["unused"] if g not in (self.sp.added[25], self.sp.added[10])]
        bad["unused"].remove(self.f0["unused"][0])
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            self.write_run(t, bad)
            rc, text = self.check(t)
        self.assertEqual(rc, 1)
        self.assertIn("breach: an added genome of no known BioProject in training: 1", text)
        self.assertIn("breach: an added copy of a tested genome in training: 1", text)
        self.assertIn("breach: a default genome moved: 1", text)

    def test_study_main_writes_a_file_only_where_a_genome_trains(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            self.write_run(t)
            other = f"{common.SPECIES}__drug_b"
            f1 = dict(self.f0, dataset=other)
            common.write_json(t / "plans" / "reviewed" / "folds" / f"{other}.study-whole.json", f1)
            common.write_ids(t / "plans" / "study.txt", [D, other])
            common.write_json(t / "curated" / "datasets.json",
                              [{"id": D, "species": "Synthia coli", "antibiotic": "drug a"},
                               {"id": other, "species": "Synthia coli", "antibiotic": "drug b"}])
            ph = [{"genome_id": g, "species": "Synthia coli", "antibiotic": "drug a", "reported": "R" if i % 3 else "S"}
                  for i, g in enumerate(self.sp.base + self.sp.added)]
            # drug b: only the added genomes of held-out BioProjects carry a call, so none trains.
            ph += [{"genome_id": g, "species": "Synthia coli", "antibiotic": "drug b", "reported": "S"}
                   for g in self.sp.base + self.sp.added[:8]]
            ph += [{"genome_id": self.sp.added[15], "species": "Synthia coli", "antibiotic": "drug b",
                    "reported": "R", "conflict": True}]
            (t / "curated" / "phenotypes.jsonl").write_text("".join(json.dumps(p) + "\n" for p in ph), encoding="utf-8")
            out = t / "arm"
            _, text = common.quiet(study.main, ["study.py", str(t / "plans" / "reviewed" / "folds"),
                                                str(t / "plans" / "study.txt"), str(t / "cohorts"), "cap-all",
                                                str(t / "curated"), str(t / "meta.jsonl.gz"), str(t / "structure.json"),
                                                str(out), str(t / "study.tsv")])
            written = sorted(p.name for p in out.iterdir())
            rows = [r.split("\t") for r in (t / "study.tsv").read_text(encoding="utf-8").splitlines()]
            doc = json.loads((out / f"{D}.study-whole.json").read_text(encoding="utf-8"))
        self.assertEqual(written, [f"{D}.study-whole.json"])
        self.assertEqual(rows[0], study.COLUMNS)
        self.assertEqual(rows[1:], [[D, "cap-all", "30", "10", "8", "10", "2", "1"], [other, "cap-all", "8", "0", "8", "0", "0", "0"]])
        self.assertEqual(doc["test"], self.f0["test"])
        self.assertIn("fold files written 1", text)


if __name__ == "__main__":
    unittest.main()
