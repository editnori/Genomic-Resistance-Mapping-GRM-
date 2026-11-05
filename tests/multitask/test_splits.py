"""splits.py on seeded synthetic genomes: one fold per genome for every task,
whole components, eligibility, rolling-origin folds and the lockbox guard.

usage: python3 -I -B -m unittest discover -s tests/multitask -t .
"""
import copy
import importlib.util
from pathlib import Path
import random
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]


def load(name):
    """Load scripts/multitask/<name>.py under the key splits.py uses for its sibling."""
    key = "grm_multitask_" + name
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, ROOT / "scripts" / "multitask" / (name + ".py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


lb = load("lockbox")
sp = load("splits")

CUTOFFS = ["2008-12-31", "2012-12-31", "2016-12-31", "2019-12-31"]


def synthetic(n=360, seed=11):
    """Uneven lineages, duplicate groups across lineages, partial and missing dates,
    two 0/1-labelled tasks, one names-only task and one rare task."""
    rng = random.Random(seed)
    lineages = ["L%02d" % i for i in range(18)]
    weights = [rng.randint(1, 8) for _ in lineages]
    genomes = []
    for i in range(n):
        r = rng.random()
        year = rng.randint(2003, 2023)
        date = None if r < 0.08 else str(year) if r < 0.15 else "%d-%02d-%02d" % (year, rng.randint(1, 12), rng.randint(1, 28))
        tasks = {}
        if rng.random() < 0.9:
            tasks["INH"] = int(rng.random() < 0.3)
        if rng.random() < 0.8:
            tasks["RIF"] = int(rng.random() < 0.25)
        if rng.random() < 0.5:
            tasks["EMB"] = None
        if i % 60 == 0:
            tasks["PZA"] = int(i % 120 == 0)
        genomes.append({"genome_id": "G%04d" % i, "lineage": rng.choices(lineages, weights)[0],
                        "duplicate_group": None, "collection_date": date, "tasks": tasks})
    free = list(range(n))
    rng.shuffle(free)
    for group in range(30):
        size = rng.choice([2, 3])
        if group < 3:
            members = free[:size]
        else:
            lineage = genomes[free[0]]["lineage"]
            members = [i for i in free if genomes[i]["lineage"] == lineage][:size]
        for i in members:
            genomes[i]["duplicate_group"] = "D%02d" % group
            free.remove(i)
    return genomes


def groups_by(genomes, field):
    groups = {}
    for g in genomes:
        if g[field] is not None:
            groups.setdefault(g[field], set()).add(g["genome_id"])
    return groups


def lockbox_and_development(genomes):
    box = lb.define_lockbox(genomes, recent_after="2020-12-31", access_log="/synthetic/never-written/access.jsonl",
                            holdout_lineages=["L03"], min_per_task=3, min_class=1)
    return box, lb.development_genomes(genomes, box)


class GlobalFolds(unittest.TestCase):
    def setUp(self):
        self.genomes = synthetic()

    def test_one_fold_per_genome_shared_by_every_task(self):
        folds = sp.global_folds(self.genomes, 5, 7, "random", lockbox=None)
        self.assertEqual(sorted(folds), sorted(g["genome_id"] for g in self.genomes))
        self.assertEqual(set(folds.values()), set(range(5)))
        by_task = {}
        for g in self.genomes:
            for task in g["tasks"]:
                by_task.setdefault(task, set()).add(g["genome_id"])
        sets = sp.fold_sets(folds)
        for fold in sets:
            self.assertFalse(set(fold["train"]) & set(fold["test"]))
            self.assertEqual(set(fold["train"]) | set(fold["test"]), set(folds))
        for task, labelled in by_task.items():
            tested = [gid for fold in sets for gid in fold["test"] if gid in labelled]
            self.assertEqual(sorted(tested), sorted(labelled), task)
            for fold in sets:
                self.assertEqual({g for g in labelled if folds[g] == fold["fold"]}, set(fold["test"]) & labelled)

    def test_duplicate_groups_and_lineages_stay_whole(self):
        duplicates, lineages = groups_by(self.genomes, "duplicate_group"), groups_by(self.genomes, "lineage")
        for scheme, whole in (("random", [duplicates]), ("lineage", [duplicates, lineages])):
            for seed in (1, 2, 3):
                folds = sp.global_folds(self.genomes, 5, seed, scheme, lockbox=None)
                for groups in whole:
                    for name, members in groups.items():
                        self.assertEqual(len({folds[g] for g in members}), 1, (scheme, seed, name))
        random_folds = sp.global_folds(self.genomes, 5, 1, "random", lockbox=None)
        self.assertTrue(any(len({random_folds[g] for g in members}) > 1 for members in lineages.values()))

    def test_duplicate_group_joins_lineages_into_one_unit(self):
        genomes = [{"genome_id": gid, "lineage": lineage, "duplicate_group": dup, "tasks": []}
                   for gid, lineage, dup in (("a", "L1", None), ("b", "L1", "D"), ("c", "L2", "D"), ("d", "L2", None),
                                             ("e", "L3", None), ("f", "L4", "E"), ("g", "L5", "E"))]
        self.assertEqual(sp.units(genomes, "lineage"), [["a", "b", "c", "d"], ["e"], ["f", "g"]])
        self.assertEqual(sp.units(genomes, "random"), [["a"], ["b", "c"], ["d"], ["e"], ["f", "g"]])
        with self.assertRaisesRegex(ValueError, "3 lineage units cannot fill 4 folds"):
            sp.global_folds(genomes, 4, 0, "lineage", lockbox=None)

    def test_eligibility_counts_match_a_recount(self):
        folds = sp.global_folds(self.genomes, 5, 7, "lineage", lockbox=None)
        report = sp.eligibility_counts(self.genomes, folds, tasks=["STR"], min_test=10, min_class=2)
        self.assertEqual(sorted(report), ["EMB", "INH", "PZA", "RIF", "STR"])
        for task, entry in report.items():
            self.assertEqual([f["fold"] for f in entry["folds"]], list(range(5)))
            for f in entry["folds"]:
                for side, members in (("test", lambda g: folds[g] == f["fold"]), ("train", lambda g: folds[g] != f["fold"])):
                    rows = [g for g in self.genomes if task in g["tasks"] and members(g["genome_id"])]
                    self.assertEqual(f[side]["n"], len(rows), (task, side))
                    if task in ("INH", "RIF", "PZA"):
                        self.assertEqual(f[side]["classes"], {c: sum(str(g["tasks"][task]) == c for g in rows) for c in "01"})
                    else:
                        self.assertIsNone(f[side]["classes"])
            bad = [f["fold"] for f in entry["folds"] if f["status"] == "insufficient_evidence"]
            self.assertEqual(entry["insufficient_folds"], bad)
            self.assertEqual(entry["status"], "insufficient_evidence" if bad else "ok")
        self.assertEqual(report["PZA"]["insufficient_folds"], list(range(5)))
        self.assertIn("test has 0 labelled genomes, below 10", report["STR"]["folds"][0]["reasons"])
        self.assertEqual(report["INH"]["status"], "ok")
        some = report["PZA"]["folds"][0]["reasons"]
        self.assertTrue(any(r.startswith("test has") for r in some) and any(r.startswith("train has") for r in some))

    def test_minimums_decide_status(self):
        folds = sp.global_folds(self.genomes, 5, 7, "random", lockbox=None)
        loose = sp.eligibility_counts(self.genomes, folds, min_test=1, min_class=1)
        self.assertEqual(loose["INH"]["status"], "ok")
        self.assertEqual(loose["EMB"]["status"], "ok")
        strict = sp.eligibility_counts(self.genomes, folds, min_test={"INH": 1000, "RIF": 1, "EMB": 1, "PZA": 1}, min_class=1)
        self.assertEqual(strict["INH"]["insufficient_folds"], list(range(5)))
        self.assertEqual(strict["RIF"]["status"], "ok")

    def test_minimums_below_one_and_missing_tasks_are_refused(self):
        # The review's G case: min_class=0 marked ok a fold whose training side lacked a class.
        folds = sp.global_folds(self.genomes, 5, 7, "random", lockbox=None)
        for min_test, min_class, message in ((0, 1, "min_test"), (10, 0, "min_class"), (-1, 1, "min_test"),
                                             (True, 1, "min_test"), (10, 1.5, "min_class"),
                                             ({"INH": 0, "RIF": 1, "EMB": 1, "PZA": 1, "STR": 1}, 1, "min_test\\['INH'\\]"),
                                             ({"INH": 5}, 1, "no minimum for task EMB, PZA, RIF"),
                                             ({"INH": 5, "RIF": 5, "EMB": 5, "PZA": 5}, 1, "no minimum for task STR")):
            with self.subTest(min_test=min_test, min_class=min_class), self.assertRaisesRegex(ValueError, message):
                sp.eligibility_counts(self.genomes, folds, tasks=["STR"], min_test=min_test, min_class=min_class)
        for build in (lambda: sp.global_manifest(self.genomes, 5, 7, "random", lockbox=None, min_class=0),
                      lambda: sp.rolling_origin_manifest(self.genomes, CUTOFFS, lockbox=None, min_test=0)):
            with self.assertRaises(ValueError):
                build()

    def test_fold_summary_gives_eligible_tasks_and_largest_unit(self):
        manifest = sp.global_manifest(self.genomes, 5, 7, "lineage", lockbox=None, tasks=["STR"])
        biggest = max(len(unit) for unit in sp.units(self.genomes, "lineage"))
        summary = manifest["fold_summary"]
        self.assertEqual([record["fold"] for record in summary], list(range(5)))
        self.assertEqual([record["genomes"] for record in summary], manifest["balance"]["genomes"]["counts"])
        self.assertEqual(max(record["largest_unit"] for record in summary), biggest)
        self.assertEqual(max(record["largest_unit_share"] for record in summary), round(biggest / len(self.genomes), 6))
        # Per fold, not only the maximum over folds (the re-check's I1).
        for record in summary:
            sizes = [len(unit) for unit in sp.units(self.genomes, "lineage") if manifest["folds"][unit[0]] == record["fold"]]
            self.assertEqual((record["genomes"], record["largest_unit"], record["largest_unit_share"]),
                             (sum(sizes), max(sizes), round(max(sizes) / len(self.genomes), 6)))
        self.assertTrue(any(record["largest_unit"] < biggest for record in summary))
        for record in summary:
            expected = sorted(task for task, entry in manifest["eligibility"].items()
                              if entry["folds"][record["fold"]]["status"] == "ok")
            self.assertEqual(sp.eligible_tasks(record), expected)
            self.assertNotIn("PZA", expected)
            self.assertNotIn("STR", expected)
        self.assertIn("INH", sp.eligible_tasks(summary[0]))
        rolling = sp.rolling_origin_manifest(self.genomes, CUTOFFS, lockbox=None, min_test=5, min_class=1)
        self.assertEqual(sp.eligible_tasks(rolling["fold_summary"][0]),
                         sorted(t for t, e in rolling["eligibility"].items() if e["folds"][0]["status"] == "ok"))
        with self.assertRaises(ValueError):
            sp.eligible_tasks(manifest["eligibility"])

    def test_stratification_balances_every_large_stratum(self):
        manifest = sp.global_manifest(self.genomes, 5, 7, "random", lockbox=None)
        for stratum in ("genomes", "task=INH", "task=RIF", "task=EMB", "task=INH,class=1", "task=RIF,class=1"):
            self.assertLess(manifest["balance"][stratum]["max_relative_deviation"], 0.1, stratum)
        # One lineage unit outweighs an even share: it should sit alone, with the rest spread evenly.
        lineage = sp.global_manifest(self.genomes, 5, 7, "lineage", lockbox=None)
        biggest = max(len(unit) for unit in sp.units(self.genomes, "lineage"))
        counts = sorted(lineage["balance"]["genomes"]["counts"])
        self.assertGreater(biggest, len(self.genomes) / 5)
        self.assertEqual(counts[-1], biggest)
        self.assertLessEqual(counts[-2] - counts[0], 2)

    def test_deterministic_and_order_free(self):
        manifest = sp.global_manifest(self.genomes, 5, 7, "lineage", lockbox=None)
        shuffled = copy.deepcopy(self.genomes)
        random.Random(5).shuffle(shuffled)
        for g in shuffled:
            g["tasks"] = dict(reversed(list(g["tasks"].items())))
        self.assertEqual(sp.global_manifest(shuffled, 5, 7, "lineage", lockbox=None), manifest)
        self.assertEqual(sp.verify_manifest(manifest), manifest["sha256"])
        self.assertNotEqual(sp.global_folds(self.genomes, 5, 8, "random", lockbox=None),
                            sp.global_folds(self.genomes, 5, 7, "random", lockbox=None))
        tampered = copy.deepcopy(manifest)
        tampered["folds"]["G0000"] = (tampered["folds"]["G0000"] + 1) % 5
        with self.assertRaisesRegex(ValueError, "sha256"):
            sp.verify_manifest(tampered)
        other_kind = {key: value for key, value in manifest.items() if key != "sha256"}
        other_kind["kind"] = "grm.multitask.other"
        other_kind["sha256"] = lb.sha256_json(other_kind)
        with self.assertRaisesRegex(ValueError, "not a GRM split manifest"):
            sp.verify_manifest(other_kind)

    def test_rejects_bad_parameters(self):
        for k, seed, scheme in ((1, 0, "random"), (True, 0, "random"), (5, "0", "random"), (5, 0, "species")):
            with self.subTest(k=k, seed=seed, scheme=scheme), self.assertRaises(ValueError):
                sp.global_folds(self.genomes, k, seed, scheme, lockbox=None)
        with self.assertRaises(TypeError):
            sp.global_folds(self.genomes, 5, 0, "random")


class RollingOrigin(unittest.TestCase):
    def setUp(self):
        self.genomes = synthetic()
        self.folds = sp.rolling_origin_dev_folds(self.genomes, CUTOFFS, lockbox=None)

    def test_monotone_windows(self):
        self.assertEqual([(f["fold"], f["cutoff"], f["test_until"]) for f in self.folds],
                         [(0, CUTOFFS[0], CUTOFFS[1]), (1, CUTOFFS[1], CUTOFFS[2]), (2, CUTOFFS[2], CUTOFFS[3]),
                          (3, CUTOFFS[3], None)])
        for before, after in zip(self.folds, self.folds[1:]):
            self.assertLess(set(before["train"]), set(after["train"]))
            self.assertEqual(set(after["train"]) - set(before["train"]), set(before["test"]))
        tested = [gid for f in self.folds for gid in f["test"]]
        self.assertEqual(len(tested), len(set(tested)))
        by_id = {g["genome_id"]: g for g in self.genomes}
        for f in self.folds:
            self.assertTrue(f["train"] and f["test"])
            self.assertFalse(set(f["train"]) & set(f["test"]))
            for gid in f["train"]:
                if by_id[gid]["collection_date"]:
                    self.assertLessEqual(lb.latest_day(by_id[gid]["collection_date"]).isoformat(), f["cutoff"])

    def test_duplicates_stay_whole_and_undated_are_listed(self):
        duplicates = groups_by(self.genomes, "duplicate_group")
        for f in self.folds:
            for members in duplicates.values():
                self.assertFalse(members & set(f["train"]) and members & set(f["test"]))
        manifest = sp.rolling_origin_manifest(self.genomes, CUTOFFS, lockbox=None)
        placed = set(self.folds[0]["train"]).union(*(f["test"] for f in self.folds))
        self.assertEqual(placed | set(manifest["undated"]), {g["genome_id"] for g in self.genomes})
        self.assertFalse(placed & set(manifest["undated"]))
        undated_dup = [members for members in duplicates.values()
                       if any(self.genomes[int(g[1:])]["collection_date"] is None for g in members)
                       and any(self.genomes[int(g[1:])]["collection_date"] for g in members)]
        self.assertTrue(undated_dup)
        self.assertTrue(all(members <= placed for members in undated_dup))

    def test_deterministic_and_validated(self):
        shuffled = copy.deepcopy(self.genomes)
        random.Random(9).shuffle(shuffled)
        manifest = sp.rolling_origin_manifest(self.genomes, CUTOFFS, lockbox=None, min_test=5, min_class=1)
        self.assertEqual(sp.rolling_origin_manifest(shuffled, CUTOFFS, lockbox=None, min_test=5, min_class=1), manifest)
        self.assertEqual(manifest["folds"], self.folds)
        self.assertEqual(sp.verify_manifest(manifest), manifest["sha256"])
        self.assertEqual(manifest["eligibility"]["PZA"]["status"], "insufficient_evidence")
        for bad in ([], ["2010-01-01", "2010-01-01"], ["2012-01-01", "2010-01-01"], ["2010"], "2010-01-01"):
            with self.subTest(cutoffs=bad), self.assertRaises(ValueError):
                sp.rolling_origin_dev_folds(self.genomes, bad, lockbox=None)


class LockboxExclusion(unittest.TestCase):
    def setUp(self):
        self.genomes = synthetic()
        self.box, self.development = lockbox_and_development(self.genomes)

    def test_split_functions_refuse_lockbox_genomes(self):
        self.assertGreater(self.box["counts"]["lockbox_genomes"], 20)
        calls = (lambda g, box: sp.global_folds(g, 5, 7, "lineage", lockbox=box),
                 lambda g, box: sp.global_manifest(g, 5, 7, "random", lockbox=box),
                 lambda g, box: sp.rolling_origin_dev_folds(g, CUTOFFS, lockbox=box),
                 lambda g, box: sp.rolling_origin_manifest(g, CUTOFFS, lockbox=box))
        lockbox_id = self.box["genome_ids"][0]
        planted = [copy.deepcopy(next(g for g in self.genomes if g["genome_id"] == lockbox_id))]
        sibling = dict(copy.deepcopy(self.development[0]), genome_id="NEW1", collection_date="2001-01-01",
                       duplicate_group=self.box["duplicate_groups"][0])
        tampered = dict(self.box, genome_ids=self.box["genome_ids"][1:])
        for call in calls:
            for genomes, box, message in ((self.genomes, self.box, "cannot enter development"),
                                          (self.development + planted, self.box, lockbox_id),
                                          (self.development + [sibling], self.box, "NEW1 \\(duplicate_group\\)"),
                                          (self.development, tampered, "tampered")):
                with self.assertRaisesRegex(lb.LockboxError, message):
                    call(genomes, box)

    def test_development_splits_never_hold_a_lockbox_genome(self):
        manifest = sp.global_manifest(self.development, 5, 7, "lineage", lockbox=self.box)
        self.assertEqual(manifest["lockbox_seal"], self.box["seal"])
        lb.assert_untouched(self.box, manifest["folds"])
        rolling = sp.rolling_origin_manifest(self.development, CUTOFFS, lockbox=self.box)
        self.assertEqual(rolling["lockbox_seal"], self.box["seal"])
        used = {gid for f in rolling["folds"] for side in ("train", "test") for gid in f[side]}
        lb.assert_untouched(self.box, used | set(rolling["undated"]))
        with self.assertRaises(lb.LockboxError):
            lb.assert_untouched(self.box, used | {self.box["genome_ids"][-1]})
        self.assertIsNone(sp.global_manifest(self.genomes, 5, 7, "lineage", lockbox=None)["lockbox_seal"])


if __name__ == "__main__":
    unittest.main()
