"""Synthetic-record tests for scripts/temporal/splits.py."""
import copy
import datetime as dt
import hashlib
import itertools
import random
import unittest

from scripts.temporal.splits import (
    EXCLUSION_REASONS,
    NOTES,
    AmbiguousLabel,
    GenomeRecord,
    SplitPolicy,
    canonical_json,
    input_sha256,
    label_as_of,
    manifest_sha256,
    rolling_origin,
    split_at,
)

DAY = dt.timedelta(days=1)
T = dt.date(2022, 1, 1)
EVAL = dt.date(2024, 12, 31)
HORIZON = 180
EMBARGO = 30
ANTIBIOTICS = ("INH", "RIF")


def event(recorded_at, value="R", antibiotic="INH", supersedes=None, source="r1"):
    return {"antibiotic": antibiotic, "value": value, "recorded_at": recorded_at,
            "source_release": source, "supersedes": supersedes}


def genome(genome_id, labels, collection=None, released=None, features=None,
           lineage="L1", group=None):
    return {"genome_id": genome_id, "collection_date": collection, "released_at": released,
            "features_available_at": features, "lineage": lineage,
            "duplicate_group": group or f"D-{genome_id}", "labels": labels}


def random_records(seed, count=300):
    """Synthetic genomes spread around T, with revisions and a few unknown dates."""
    rng = random.Random(seed)
    start = dt.date(2019, 1, 1)
    rows = []
    for index in range(count):
        collection = start + rng.randrange(0, 1800) * DAY
        released = collection + rng.randrange(0, 240) * DAY
        features = None
        if rng.random() < 0.3:
            features = released + rng.randrange(-40, 40) * DAY
        labels = []
        for antibiotic in ANTIBIOTICS:
            recorded = collection + rng.randrange(0, 500) * DAY
            previous = None
            for _ in range(rng.randrange(0, 4)):
                labels.append(event(recorded.isoformat(), rng.choice("RS"), antibiotic,
                                    previous and previous.isoformat(), f"r{rng.randrange(9)}"))
                previous, recorded = recorded, recorded + rng.randrange(1, 400) * DAY
        if rng.random() < 0.1:
            labels.append(event(None, rng.choice("RS"), rng.choice(ANTIBIOTICS)))
        rows.append(genome(
            f"G{index:04d}", labels,
            None if rng.random() < 0.04 else collection.isoformat(),
            None if rng.random() < 0.04 else released.isoformat(),
            features and features.isoformat(),
            rng.choice(["L1", "L2", f"L{collection.year}", f"L{collection.year}b", None]),
            None if rng.random() < 0.03 else
            rng.choice([f"D{rng.randrange(20)}", f"D{collection.year}-{rng.randrange(40)}"])))
    return rows


def split(rows, **overrides):
    args = {"T": T, "horizon": HORIZON, "embargo": EMBARGO, "eval_date": EVAL,
            "antibiotic": "INH", "policy": SplitPolicy()}
    args.update(overrides)
    return split_at(rows, **args)


def by_id(rows):
    return {row["genome_id"]: GenomeRecord.from_mapping(row) for row in rows}


def after_cutoff(row, days):
    """An ISO date `days` after both T and the row's collection date."""
    collection = row["collection_date"]
    floor = max(T, dt.date.fromisoformat(collection)) if collection else T
    return (floor + days * DAY).isoformat()


def place(manifest, genome_id):
    for arm in ("train", "test"):
        if genome_id in manifest[arm]:
            return f"{arm}={manifest[arm][genome_id]}"
    for reason, ids in manifest["excluded"].items():
        if genome_id in ids:
            return reason
    raise AssertionError(f"{genome_id} was not placed")


def leaf_paths(value, path=()):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from leaf_paths(child, path + (key,))
    elif isinstance(value, list) and value:
        for index, child in enumerate(value):
            yield from leaf_paths(child, path + (index,))
    else:
        yield path


class LabelAsOfTest(unittest.TestCase):
    record = GenomeRecord.from_mapping(genome("G1", [
        event("2021-03-01", "S"),
        event("2021-09-01", "R", supersedes="2021-03-01"),
        event("2022-06-01", "S", supersedes="2021-09-01"),
        event(None, "S"),
        event("2020-01-01", "R", antibiotic="RIF"),
    ]))

    def test_returns_latest_dated_event_recorded_on_or_before_cutoff(self):
        self.assertIsNone(label_as_of(self.record, "INH", "2021-02-28"))
        self.assertEqual(label_as_of(self.record, "INH", "2021-03-01").value, "S")
        self.assertEqual(label_as_of(self.record, "INH", "2021-08-31").value, "S")
        self.assertEqual(label_as_of(self.record, "INH", T).value, "R")
        self.assertEqual(label_as_of(self.record, "INH", EVAL).value, "S")
        self.assertIsNone(label_as_of(self.record, "EMB", EVAL))

    def test_same_day_events_collapse_or_are_ambiguous_only_where_read(self):
        agreeing = GenomeRecord.from_mapping(genome("G2", [
            event("2021-01-01", "R", source="b"), event("2021-01-01", "R", source="a")]))
        self.assertEqual(label_as_of(agreeing, "INH", T).source_release, "a")
        tied = GenomeRecord.from_mapping(genome("G3", [
            event("2021-01-01", "R", source="a"), event("2021-01-01", "S", source="b"),
            event("2021-06-01", "S", supersedes="2021-01-01")]))
        self.assertIsNone(label_as_of(tied, "INH", "2020-12-31"))
        with self.assertRaises(AmbiguousLabel):
            label_as_of(tied, "INH", "2021-01-01")
        with self.assertRaises(AmbiguousLabel):
            label_as_of(tied, "INH", "2021-05-31")
        self.assertEqual(label_as_of(tied, "INH", "2021-06-01").value, "S")
        self.assertTrue(issubclass(AmbiguousLabel, ValueError))

    def test_agreeing_tie_pick_is_the_first_in_canonical_order_for_any_input_order(self):
        earlier = [event("2021-01-01", source="a"), event("2021-01-15", source="a")]
        for tied, expected in (
                ([event("2021-02-01", source="b"),
                  event("2021-02-01", source="a", supersedes="2021-01-15"),
                  event("2021-02-01", source="a"),
                  event("2021-02-01", source="a", supersedes="2021-01-01")],
                 ("a", None)),
                ([event("2021-02-01", source="a", supersedes="2021-01-15"),
                  event("2021-02-01", source="b"),
                  event("2021-02-01", source="a", supersedes="2021-01-01")],
                 ("a", dt.date(2021, 1, 1)))):
            picks = set()
            for order in itertools.permutations(tied):
                record = GenomeRecord.from_mapping(genome("G", earlier + list(order)))
                pick = label_as_of(record, "INH", T)
                picks.add(pick)
                self.assertEqual((pick.source_release, pick.supersedes), expected)
                shuffled = list(reversed(earlier + list(order)))
                random.Random(len(picks)).shuffle(shuffled)
                record = GenomeRecord.from_mapping(genome("G", shuffled))
                self.assertEqual(label_as_of(record, "INH", T), pick)
            self.assertEqual(len(picks), 1)


class AvailabilityTest(unittest.TestCase):
    def test_no_train_genome_has_features_after_cutoff_minus_embargo(self):
        checked = 0
        for seed in range(5):
            rows = by_id(random_records(seed))
            for embargo in (0, 30, 365):
                manifest = split(list(rows.values()), embargo=embargo)
                for genome_id in manifest["train"]:
                    record = rows[genome_id]
                    for day in (record.released_at, record.features_available_at):
                        if day is not None:
                            self.assertLessEqual(day, T - embargo * DAY)
                    checked += 1
        self.assertGreater(checked, 500)

    def test_availability_is_the_later_of_release_and_features(self):
        rows = [genome("features_later", [event("2021-01-01")], released="2021-01-01",
                       features="2021-12-31"),
                genome("released_only", [event("2021-01-01")], released="2021-12-31"),
                genome("features_only", [event("2021-01-01")], features="2021-01-01"),
                genome("features_earlier", [event("2021-06-01")], collection="2020-06-01",
                       released="2023-01-01", features="2021-01-01")]
        manifest = split(rows)
        self.assertEqual(manifest["train"], {"features_only": "R"})
        self.assertEqual(manifest["excluded"]["embargo"], ["features_later", "released_only"])
        self.assertEqual(manifest["excluded"]["features_after_cutoff"], ["features_earlier"])
        self.assertEqual(manifest["flags"]["inconsistent_dates"], ["features_earlier"])
        self.assertEqual(manifest["counts"]["inconsistent_dates"], 1)

    def test_features_on_the_release_day_are_consistent(self):
        rows = [genome("same_day", [event("2021-06-01")], released="2021-01-01",
                       features="2021-01-01"),
                genome("day_before", [event("2021-06-01")], released="2021-01-02",
                       features="2021-01-01")]
        self.assertEqual(split(rows)["flags"]["inconsistent_dates"], ["day_before"])

    def test_extreme_horizon_embargo_or_cutoff_raise_value_error(self):
        for overrides in ({"horizon": 10**9, "eval_date": T},
                          {"T": "9999-12-31", "horizon": 1, "embargo": 0,
                           "eval_date": "9999-12-31"},
                          {"horizon": 1, "embargo": 10**6, "eval_date": T + DAY},
                          {"T": "0001-01-01", "embargo": 1, "eval_date": "0001-12-31"}):
            with self.assertRaises(ValueError, msg=overrides):
                split([], **overrides)

    def test_boundary_days(self):
        rows = [
            genome("feat_eq_T_minus_E", [event("2021-06-01")],
                   released=(T - EMBARGO * DAY).isoformat()),
            genome("feat_eq_T_minus_E_plus1", [event("2021-06-01")],
                   released=(T - (EMBARGO - 1) * DAY).isoformat()),
            genome("feat_eq_T", [event("2021-06-01")], released=T.isoformat()),
            genome("feat_eq_T_plus1", [event("2021-06-01")], released=(T + DAY).isoformat()),
            genome("label_eq_T", [event(T.isoformat())], released="2021-01-01"),
            genome("label_eq_T_plus1", [event((T + DAY).isoformat())], collection="2022-01-02",
                   released="2022-01-02"),
            genome("col_eq_T", [event("2022-02-01")], collection=T.isoformat(),
                   released="2021-12-01"),
            genome("col_eq_T_plus_H", [event("2022-08-01")],
                   collection=(T + HORIZON * DAY).isoformat(), released="2022-07-01"),
            genome("col_eq_T_plus_H_plus1", [event("2022-08-01")],
                   collection=(T + (HORIZON + 1) * DAY).isoformat(), released="2022-07-01"),
            genome("label_eq_eval", [event(EVAL.isoformat())], collection="2022-02-01",
                   released="2022-02-01"),
            genome("label_eq_eval_plus1", [event((EVAL + DAY).isoformat())],
                   collection="2022-02-01", released="2022-02-01"),
            genome("feat_eq_eval", [event("2022-03-01")], collection="2022-02-01",
                   released=EVAL.isoformat()),
            genome("feat_eq_eval_plus1", [event("2022-03-01")], collection="2022-02-01",
                   released=(EVAL + DAY).isoformat()),
            genome("bf_label_eq_T_plus_H", [event((T + HORIZON * DAY).isoformat())],
                   collection="2021-01-01", released="2021-01-01"),
            genome("bf_label_eq_T_plus_H_plus1", [event((T + (HORIZON + 1) * DAY).isoformat())],
                   collection="2021-01-01", released="2021-01-01"),
        ]
        expected = {
            "feat_eq_T_minus_E": "train=R", "feat_eq_T_minus_E_plus1": "embargo",
            "feat_eq_T": "embargo", "feat_eq_T_plus1": "features_after_cutoff",
            "label_eq_T": "train=R", "label_eq_T_plus1": "test=R", "col_eq_T": "test=R",
            "col_eq_T_plus_H": "test=R", "col_eq_T_plus_H_plus1": "after_horizon",
            "label_eq_eval": "test=R", "label_eq_eval_plus1": "no_label",
            "feat_eq_eval": "test=R", "feat_eq_eval_plus1": "features_after_eval",
            "bf_label_eq_T_plus_H": "test=R",
            "bf_label_eq_T_plus_H_plus1": "backfill_outside_horizon",
        }
        manifest = split(rows, policy=SplitPolicy(backfill="test"))
        self.assertEqual({genome_id: place(manifest, genome_id) for genome_id in expected},
                         expected)
        self.assertIn("col_eq_T", manifest["flags"]["backfill"])
        self.assertEqual(place(split(rows), "col_eq_T"), "backfill")
        self.assertEqual(split(rows, embargo=0)["train"]["feat_eq_T"], "R")

    def test_eval_date_may_equal_cutoff_plus_horizon(self):
        rows = [genome("G", [event("2022-03-01")], collection="2022-02-01",
                       released="2022-02-01")]
        manifest = split(rows, eval_date=(T + HORIZON * DAY).isoformat())
        self.assertEqual(manifest["test"], {"G": "R"})
        with self.assertRaises(ValueError):
            split(rows, eval_date=(T + (HORIZON - 1) * DAY).isoformat())

    def test_train_labels_equal_label_as_of_cutoff_and_test_labels_as_of_eval(self):
        for seed in range(5):
            rows = by_id(random_records(seed))
            for antibiotic in ANTIBIOTICS:
                manifest = split(list(rows.values()), antibiotic=antibiotic)
                self.assertTrue(manifest["train"] and manifest["test"])
                for genome_id, value in manifest["train"].items():
                    self.assertEqual(value, label_as_of(rows[genome_id], antibiotic, T).value)
                for genome_id, value in manifest["test"].items():
                    self.assertIsNone(label_as_of(rows[genome_id], antibiotic, T))
                    self.assertEqual(value, label_as_of(rows[genome_id], antibiotic, EVAL).value)

    def test_revision_after_cutoff_does_not_change_train(self):
        rows = [genome("G1", [event("2021-06-01", "S")], released="2021-06-01")]
        before = split(rows)
        rows[0]["labels"].append(event("2022-03-01", "R", supersedes="2021-06-01"))
        after = split(rows)
        self.assertEqual(before["train"], {"G1": "S"})
        self.assertEqual(after["train"], {"G1": "S"})
        self.assertEqual(after["test"], {})

    def test_every_record_lands_in_exactly_one_place(self):
        for policy in (SplitPolicy(), SplitPolicy("release", "test")):
            rows = random_records(3)
            manifest = split(rows, policy=policy)
            places = list(manifest["train"]) + list(manifest["test"])
            for ids in manifest["excluded"].values():
                places.extend(ids)
            self.assertEqual(sorted(places), sorted(row["genome_id"] for row in rows))
            self.assertEqual(sum(manifest["counts"]["excluded"].values())
                             + manifest["counts"]["train"] + manifest["counts"]["test"],
                             manifest["counts"]["records"])

    def test_window_policy_picks_collection_or_release_date(self):
        rows = [genome("G1", [event("2022-03-01")], collection="2021-12-01",
                       released="2022-02-01"),
                genome("G2", [event("2022-03-01")], collection="2022-02-01",
                       released="2022-09-01")]
        collection = split(rows)
        self.assertEqual(collection["test"], {"G2": "R"})
        self.assertEqual(collection["flags"]["backfill"], ["G1"])
        release = split(rows, policy=SplitPolicy(window="release"))
        self.assertEqual(release["test"], {"G1": "R"})
        self.assertEqual(release["excluded"]["after_horizon"], ["G2"])

    def test_label_known_at_cutoff_never_enters_test(self):
        rows = [genome(f"P{index}", [event("2021-10-01")], released="2022-03-01")
                for index in range(5)]
        rows.append(genome("revised", [event("2021-06-01", "S"),
                                       event("2022-03-01", "R", supersedes="2021-06-01")],
                           collection="2021-01-01", released="2021-01-01"))
        manifest = split(rows, policy=SplitPolicy("release", "test"))
        self.assertEqual(manifest["train"], {"revised": "S"})
        self.assertEqual(manifest["test"], {})
        self.assertEqual(manifest["counts"]["excluded"]["features_after_cutoff"], 5)
        self.assertEqual(manifest["notes"], list(NOTES))


class PostCutoffHistoryTest(unittest.TestCase):
    """Reviewer M-1: label history recorded after T must not change train."""

    def test_reviewer_repro_post_cutoff_conflict_and_eval_date(self):
        clean = [genome("G", [event("2021-06-01", "R")], released="2021-01-01")]
        leaky = [genome("G", [event("2021-06-01", "R"), event("2023-05-01", "R", source="a"),
                              event("2023-05-01", "S", source="b")], released="2021-01-01")]
        self.assertEqual(split(clean)["train"], {"G": "R"})
        self.assertEqual(split(leaky)["train"], {"G": "R"})
        self.assertEqual(split(leaky, eval_date="2023-04-30")["train"], {"G": "R"})

    def test_undated_event_does_not_change_train(self):
        rows = [genome("G", [event("2021-06-01", "R"), event(None, "S")],
                       released="2021-01-01")]
        self.assertEqual(split(rows)["train"], {"G": "R"})

    def test_conflict_read_at_cutoff_excludes_and_conflict_read_at_eval_blocks_test(self):
        rows = [
            genome("tie_at_T", [event("2021-06-01", "R", source="a"),
                                event("2021-06-01", "S", source="b")], released="2021-01-01"),
            genome("tie_resolved_before_T", [
                event("2021-03-01", "R", source="a"), event("2021-03-01", "S", source="b"),
                event("2021-06-01", "S", supersedes="2021-03-01")], released="2021-01-01"),
            genome("tie_at_eval", [event("2022-03-01", "R", source="a"),
                                   event("2022-03-01", "S", source="b")],
                   collection="2022-02-01", released="2022-02-01"),
            genome("tie_agreeing", [event("2021-06-01", "S", source="b"),
                                    event("2021-06-01", "S", source="a")],
                   released="2021-01-01"),
        ]
        manifest = split(rows)
        self.assertEqual(manifest["train"], {"tie_resolved_before_T": "S", "tie_agreeing": "S"})
        self.assertEqual(manifest["excluded"]["ambiguous_label"], ["tie_at_T", "tie_at_eval"])
        self.assertEqual(manifest["test"], {})

    def test_post_cutoff_junk_and_eval_date_never_change_train(self):
        for seed in range(4):
            rows = random_records(seed)
            baseline = split(rows)["train"]
            rng = random.Random(seed)
            for row in rows:
                day = after_cutoff(row, rng.randrange(1, 900))
                row["labels"] += [event(day, "R", source="x"), event(day, "S", source="y"),
                                  event(None, rng.choice("RS"), rng.choice(ANTIBIOTICS))]
            for eval_date in (EVAL, "2022-06-30", "2030-01-01"):
                self.assertEqual(split(rows, eval_date=eval_date)["train"], baseline)


class BackfillTest(unittest.TestCase):
    rows = [
        genome("old_late_label", [event("2022-03-01", "S")], collection="2021-01-01",
               released="2021-02-01"),
        genome("old_very_late_label", [event("2023-03-01", "S")], collection="2021-01-01",
               released="2021-02-01"),
        genome("new", [event("2022-03-01")], collection="2022-02-01", released="2022-02-15"),
    ]

    def test_backfill_is_excluded_from_train_and_flagged(self):
        manifest = split(self.rows)
        self.assertEqual(manifest["train"], {})
        self.assertEqual(manifest["test"], {"new": "R"})
        self.assertEqual(manifest["excluded"]["backfill"],
                         ["old_late_label", "old_very_late_label"])
        self.assertEqual(manifest["flags"]["backfill"],
                         ["old_late_label", "old_very_late_label"])
        self.assertEqual(manifest["counts"]["backfill"], 2)

    def test_test_policy_admits_backfill_labelled_within_horizon(self):
        manifest = split(self.rows, policy=SplitPolicy(backfill="test"))
        self.assertEqual(manifest["train"], {})
        self.assertEqual(manifest["test"], {"new": "R", "old_late_label": "S"})
        self.assertEqual(manifest["excluded"]["backfill_outside_horizon"],
                         ["old_very_late_label"])
        self.assertEqual(manifest["flags"]["backfill"],
                         ["old_late_label", "old_very_late_label"])


class UnknownDatesTest(unittest.TestCase):
    def test_unknown_dates_are_excluded_and_counted_by_reason(self):
        rows = [
            genome("only_undated", [event(None)], released="2021-01-01"),
            genome("undated_and_late", [event(None), event("2025-01-01")],
                   released="2021-01-01"),
            genome("no_features", [event("2021-01-01")]),
            genome("no_collection", [event("2022-03-01")], released="2022-02-01"),
            genome("no_collection_train", [event("2021-01-01")], released="2021-01-01"),
            genome("unlabelled", [event("2021-01-01", antibiotic="RIF")],
                   released="2021-01-01"),
            genome("no_events", [], released="2021-01-01"),
            genome("label_after_eval", [event("2025-01-01")], released="2021-01-01"),
        ]
        manifest = split(rows)
        self.assertEqual(manifest["excluded"]["unknown_label_recorded_at"],
                         ["only_undated", "undated_and_late"])
        self.assertEqual(manifest["excluded"]["no_label"],
                         ["label_after_eval", "no_events", "unlabelled"])
        counts = manifest["counts"]["excluded"]
        self.assertEqual(counts["unknown_features_available_at"], 1)
        self.assertEqual(counts["unknown_collection_date"], 1)
        self.assertEqual(manifest["train"], {"no_collection_train": "R"})
        self.assertEqual(set(counts), set(EXCLUSION_REASONS))

        release = split(rows, policy=SplitPolicy(window="release"))
        self.assertEqual(release["test"], {"no_collection": "R"})
        missing_release = [genome("G", [event("2022-03-01")], collection="2022-02-01",
                                  features="2022-02-01")]
        self.assertEqual(split(missing_release, policy=SplitPolicy(window="release"))
                         ["excluded"]["unknown_released_at"], ["G"])

    def test_random_records_report_every_unknown_reason(self):
        totals = dict.fromkeys(("unknown_label_recorded_at", "unknown_features_available_at",
                                "unknown_collection_date", "no_label"), 0)
        for seed in range(5):
            counts = split(random_records(seed))["counts"]["excluded"]
            for reason in totals:
                totals[reason] += counts[reason]
        for reason, total in totals.items():
            self.assertGreater(total, 0, reason)


class RollingOriginTest(unittest.TestCase):
    cutoffs = ["2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01", "2022-07-01"]

    def test_train_ids_grow_monotonically(self):
        rows = random_records(5)
        splits = rolling_origin(rows, self.cutoffs, HORIZON, EMBARGO, EVAL, "INH")
        self.assertEqual([item["parameters"]["T"] for item in splits], self.cutoffs)
        sizes = [len(item["train"]) for item in splits]
        self.assertEqual(sizes, sorted(sizes))
        self.assertLess(sizes[0], sizes[-1])
        for earlier, later in zip(splits, splits[1:]):
            self.assertLessEqual(set(earlier["train"]), set(later["train"]))

    def test_only_a_newly_ambiguous_label_leaves_train(self):
        rows = [genome("G", [event("2020-01-01", "R"), event("2021-03-01", "R", source="a"),
                             event("2021-03-01", "S", source="b")], released="2020-01-01"),
                genome("H", [event("2020-01-01", "S")], released="2020-01-01")]
        splits = rolling_origin(rows, self.cutoffs, HORIZON, EMBARGO, EVAL, "INH")
        self.assertEqual([place(item, "G") for item in splits],
                         ["train=R", "train=R", "ambiguous_label", "ambiguous_label",
                          "ambiguous_label"])
        for earlier, later in zip(splits, splits[1:]):
            self.assertLessEqual(set(earlier["train"]),
                                 set(later["train"]) | set(later["excluded"]["ambiguous_label"]))

    def test_rolling_origin_is_deterministic_and_order_free(self):
        rows = random_records(5)
        first = rolling_origin(rows, self.cutoffs, HORIZON, EMBARGO, EVAL, "INH")
        shuffled = list(reversed(rows))
        random.Random(1).shuffle(shuffled)
        second = rolling_origin(shuffled, self.cutoffs, HORIZON, EMBARGO, EVAL, "INH")
        self.assertEqual([item["sha256"] for item in first],
                         [item["sha256"] for item in second])
        self.assertEqual(first, second)
        self.assertEqual(first[3], split(rows))

    def test_cutoffs_must_strictly_increase(self):
        for cutoffs in ([], ["2021-01-01", "2020-01-01"], ["2021-01-01", "2021-01-01"]):
            with self.assertRaises(ValueError):
                rolling_origin([], cutoffs, HORIZON, EMBARGO, EVAL, "INH")


class NovelGroupTest(unittest.TestCase):
    train = dict(labels=[event("2021-01-01")], released="2021-01-01")
    new = dict(labels=[event("2022-03-01")], collection="2022-02-01", released="2022-02-01")

    def test_novel_arm_requires_unseen_lineage_and_duplicate_group(self):
        rows = [
            genome("train_a", lineage="L1", group="D1", **self.train),
            genome("train_b", lineage="L2", group="D2", **self.train),
            genome("seen_lineage", lineage="L1", group="D9", **self.new),
            genome("seen_group", lineage="L9", group="D2", **self.new),
            genome("both_seen", lineage="L2", group="D1", **self.new),
            genome("seen_lineage_no_group", lineage="L1", group=None, **self.new),
            genome("novel", lineage="L9", group="D9", **self.new),
            genome("unknown", lineage=None, group="D8", **self.new),
        ]
        rows[5]["duplicate_group"] = None
        manifest = split(rows)
        self.assertEqual(manifest["novel_group_test"], ["novel"])
        self.assertEqual(manifest["novel_group_excluded"], {
            "lineage_seen": ["both_seen", "seen_lineage", "seen_lineage_no_group"],
            "duplicate_group_seen": ["seen_group"],
            "group_unknown": ["unknown"],
            "novel_unverified": [],
        })
        self.assertEqual(manifest["counts"]["train_group_unknown"], 0)

    def test_null_train_group_makes_unseen_test_genomes_unverified(self):
        rows = [
            genome("train_blind", lineage=None, **self.train),
            genome("train_L1_nogroup", lineage="L1", **self.train),
            genome("test_L9_D9", lineage="L9", group="D9", **self.new),
            genome("test_L1_space", lineage=" L1 ", group="D7", **self.new),
            genome("test_nolineage_D9", lineage=None, group="D9", **self.new),
        ]
        rows[0]["duplicate_group"] = rows[1]["duplicate_group"] = None
        manifest = split(rows)
        self.assertEqual(manifest["novel_group_test"], [])
        self.assertEqual(manifest["novel_group_excluded"], {
            "lineage_seen": ["test_L1_space"],
            "duplicate_group_seen": [],
            "group_unknown": ["test_nolineage_D9"],
            "novel_unverified": ["test_L9_D9"],
        })
        self.assertEqual(manifest["counts"]["train_group_unknown"], 2)
        rows[0]["lineage"] = "L5"
        rows[0]["duplicate_group"] = "D5"
        rows[1]["duplicate_group"] = "D1"
        self.assertEqual(split(rows)["novel_group_test"], ["test_L9_D9"])

    def test_train_null_group_alone_blocks_novel_and_seen_group_beats_unknown(self):
        rows = [genome("train_L1_null_group", lineage="L1", **self.train),
                genome("train_x", lineage="L2", group="Dx", **self.train),
                genome("test_unseen", lineage="L9", group="D9", **self.new),
                genome("test_null_lineage_seen_group", lineage=None, group="Dx", **self.new)]
        rows[0]["duplicate_group"] = None
        manifest = split(rows)
        self.assertEqual(manifest["novel_group_test"], [])
        self.assertEqual(manifest["novel_group_excluded"]["novel_unverified"], ["test_unseen"])
        self.assertEqual(manifest["novel_group_excluded"]["duplicate_group_seen"],
                         ["test_null_lineage_seen_group"])
        self.assertEqual(manifest["counts"]["train_group_unknown"], 1)

    def test_rule_6_strips_then_requires_printable_ascii_and_refuses_aliases(self):
        train = genome("train", lineage="L1", group="D1", **self.train)
        stripped = [genome(f"test_group_{index}", lineage=f"L{index + 10}", group=value,
                           **self.new)
                    for index, value in enumerate(("D1 ", "\tD1", "\u00a0D1\n"))]
        manifest = split([train] + stripped)
        self.assertEqual(manifest["novel_group_excluded"]["duplicate_group_seen"],
                         ["test_group_0", "test_group_1", "test_group_2"])
        self.assertEqual({record.duplicate_group for record in by_id(stripped).values()}, {"D1"})
        self.assertEqual(by_id([genome("kept", [], lineage="Lin A")])["kept"].lineage, "Lin A")
        for name in ("duplicate_group", "lineage"):
            for value in ("D\t1", "D1\u200b", "\u200bD1", "D\u00ad1", "D1\x00", "L\u00e91",
                          "D\u00a01", "\uff24\uff11", "L\u2081", "D1\x7f", "   "):
                bad = dict(genome("test_bad", **self.new), **{name: value})
                with self.assertRaisesRegex(ValueError, f"^test_bad: {name} ", msg=repr(value)):
                    split([train, bad])
        for name, value, pair in (("duplicate_group", "d1", "'D1' and 'd1'"),
                                  ("lineage", "l1", "'L1' and 'l1'")):
            alias = dict(genome("test_alias", lineage="L9", group="D9", **self.new),
                         **{name: value})
            with self.assertRaisesRegex(ValueError, f"^{name} values {pair} collide", msg=value):
                split([train, alias])
        cross_field = genome("test_cross", lineage="d1", group="l1", **self.new)
        self.assertEqual(split([train, cross_field])["test"], {"test_cross": "R"})

    def test_rule_6_strips_every_python_whitespace_character_from_both_ends(self):
        stripped = (0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x1C, 0x1D, 0x1E, 0x1F, 0x20, 0x85, 0xA0,
                    0x1680, *range(0x2000, 0x200B), 0x2028, 0x2029, 0x202F, 0x205F, 0x3000)
        self.assertEqual(stripped, tuple(code for code in range(0x110000) if chr(code).isspace()))
        for code in stripped:
            space = chr(code)
            for value in (space + "D1", "D1" + space, space + "D1" + space):
                record = by_id([genome("G", [], lineage=value, group=value)])["G"]
                self.assertEqual((record.lineage, record.duplicate_group), ("D1", "D1"),
                                 msg=hex(code))

    def test_rule_6_accepts_tilde_and_refuses_del_inside_a_value(self):
        for name in ("lineage", "duplicate_group"):
            kept = by_id([dict(genome("G", []), **{name: "D" + chr(0x7E) + "1"})])["G"]
            self.assertEqual(getattr(kept, name), "D~1")
            with self.assertRaisesRegex(ValueError, f"^G: {name} "):
                by_id([dict(genome("G", []), **{name: "D" + chr(0x7F) + "1"})])

    def test_novel_arm_holds_on_random_records(self):
        for seed in range(5):
            rows = by_id(row for row in random_records(seed)
                         if row["lineage"] and row["duplicate_group"])
            manifest = split(list(rows.values()))
            self.assertEqual(split(list(reversed(rows.values()))), manifest)
            self.assertEqual(manifest["counts"]["train_group_unknown"], 0)
            lineages = {rows[genome_id].lineage for genome_id in manifest["train"]}
            groups = {rows[genome_id].duplicate_group for genome_id in manifest["train"]}
            self.assertTrue(manifest["novel_group_test"])
            for genome_id in manifest["novel_group_test"]:
                self.assertIn(genome_id, manifest["test"])
                self.assertNotIn(rows[genome_id].lineage, lineages)
                self.assertNotIn(rows[genome_id].duplicate_group, groups)
        with_nulls = split(random_records(0))
        self.assertGreater(with_nulls["counts"]["train_group_unknown"], 0)
        self.assertEqual(with_nulls["novel_group_test"], [])


class ManifestTest(unittest.TestCase):
    rows = [
        genome("G1", [event("2021-03-01", "S")], collection="2021-01-01",
               released="2021-02-01", lineage="L1", group="D1"),
        genome("G2", [event("2022-03-01")], collection="2022-02-01", released="2022-02-15",
               lineage="L2", group="D2"),
        genome("G3", [event("2022-03-01")], collection="2021-06-01", released="2021-07-01"),
    ]

    def test_canonical_json_sorts_keys_without_spaces(self):
        self.assertEqual(canonical_json({"b": [1, {"d": 2, "c": None}], "a": "x"}),
                         '{"a":"x","b":[1,{"c":null,"d":2}]}')

    def test_canonical_json_escapes_non_ascii_and_refuses_nan(self):
        text = canonical_json({"a": "\u00e9"})
        self.assertEqual(text, '{"a":"\\u00e9"}')
        self.assertEqual(hashlib.sha256(text.encode("utf-8")).hexdigest(),
                         "57e512a1ff02f75e939e93d2196b8d29ad3721b70ce25a601402d26410f453a2")
        with self.assertRaises(ValueError):
            canonical_json({"a": float("nan")})

    def test_hash_is_stable_across_runs_and_input_order(self):
        manifest = split(self.rows)
        self.assertEqual(manifest["sha256"], manifest_sha256(manifest))
        self.assertEqual(split(list(reversed(self.rows)))["sha256"], manifest["sha256"])
        self.assertEqual(
            manifest["sha256"],
            "92cb0612fa6267aebfcaa60198dce6f119ddb001534571e0449aa37195a6d031")
        self.assertEqual(
            manifest["input_sha256"],
            "0507e32a04dfaeb7d352168a0fe1da498aca6a64d1f9db7600c204b7eb28bfa2")

    def test_input_hash_sorts_labels_that_differ_only_in_supersedes(self):
        earlier = [event("2021-01-01"), event("2021-01-15")]
        pair = [event("2021-02-01", supersedes="2021-01-15"),
                event("2021-02-01", supersedes="2021-01-01")]
        forward = genome("G", earlier + pair, released="2021-01-01")
        backward = genome("G", earlier + pair[::-1], released="2021-01-01")
        self.assertEqual(input_sha256([forward]), input_sha256([backward]))
        self.assertEqual([row["supersedes"] for row in
                          by_id([backward])["G"].to_mapping()["labels"]],
                         [None, None, "2021-01-01", "2021-01-15"])

    def test_input_hash_covers_every_input_field(self):
        rows = random_records(2)
        base = input_sha256(rows)
        shuffled = copy.deepcopy(list(reversed(rows)))
        for row in shuffled:
            row["labels"].reverse()
        self.assertEqual(input_sha256(shuffled), base)
        train_id = next(iter(split(rows)["train"]))
        moved = copy.deepcopy(rows)
        row = next(row for row in moved if row["genome_id"] == train_id)
        row["released_at"] = "2018-01-01"
        self.assertEqual(split(moved)["train"], split(rows)["train"])
        self.assertNotEqual(split(moved)["input_sha256"], base)
        self.assertNotEqual(split(moved)["sha256"], split(rows)["sha256"])
        for key, value in (("lineage", "L-other"), ("duplicate_group", "D-other"),
                           ("collection_date", None), ("features_available_at", "2030-01-01")):
            edited = copy.deepcopy(rows)
            edited[0][key] = value
            self.assertNotEqual(input_sha256(edited), base, key)
        edited = copy.deepcopy(rows)
        edited[1]["labels"].append(event(None, antibiotic="EMB"))
        self.assertNotEqual(input_sha256(edited), base)

    def test_any_edit_changes_the_hash(self):
        manifest = split(random_records(2))
        paths = list(leaf_paths({k: v for k, v in manifest.items() if k != "sha256"}))
        self.assertGreater(len(paths), 100)
        for path in paths:
            edited = copy.deepcopy(manifest)
            parent = edited
            for key in path[:-1]:
                parent = parent[key]
            leaf = parent[path[-1]]
            if isinstance(leaf, list):
                parent[path[-1]] = ["added"]
            elif isinstance(leaf, int):
                parent[path[-1]] = leaf + 1
            elif isinstance(leaf, dict):
                parent[path[-1]] = {"added": 1}
            else:
                parent[path[-1]] = leaf + "x"
            self.assertNotEqual(manifest_sha256(edited), manifest["sha256"], path)
        renamed = copy.deepcopy(manifest)
        genome_id = next(iter(renamed["train"]))
        renamed["train"][genome_id + "x"] = renamed["train"].pop(genome_id)
        self.assertNotEqual(manifest_sha256(renamed), manifest["sha256"])

    def test_parameters_change_the_hash(self):
        base = split(self.rows)["sha256"]
        for overrides in ({"horizon": 181}, {"embargo": 31}, {"eval_date": "2025-01-01"},
                          {"T": "2022-01-02"}, {"antibiotic": "RIF"},
                          {"policy": SplitPolicy(window="release")},
                          {"policy": SplitPolicy(backfill="test")}):
            self.assertNotEqual(split(self.rows, **overrides)["sha256"], base, overrides)


class TestRecordIsolationTest(unittest.TestCase):
    def test_editing_test_records_cannot_change_train(self):
        rng = random.Random(99)
        for seed in range(3):
            rows = random_records(seed)
            baseline = split(rows)
            self.assertTrue(baseline["test"])
            test_ids = set(baseline["test"])
            for trial in range(40):
                edited = copy.deepcopy(rows)
                for row in edited:
                    if row["genome_id"] not in test_ids:
                        continue
                    row["collection_date"] = rng.choice([None, "2019-01-01", "2022-01-02"])
                    row["released_at"] = rng.choice([None, "2019-01-01", "2023-06-06"])
                    row["features_available_at"] = rng.choice([None, "2019-01-01"])
                    row["lineage"] = rng.choice([None, "L1", "L99"])
                    row["duplicate_group"] = rng.choice([None, "D1", "D999"])
                    for item in row["labels"]:
                        item["value"] = rng.choice("RS")
                        item["recorded_at"] = rng.choice(
                            [None, (T + rng.randrange(1, 900) * DAY).isoformat()])
                        item["supersedes"] = None
                    row["labels"].append(event((T + DAY).isoformat(), rng.choice("RS")))
                if trial % 4 == 0:
                    edited = [row for row in edited if row["genome_id"] not in test_ids]
                self.assertEqual(split(edited)["train"], baseline["train"])


class ValidationTest(unittest.TestCase):
    def test_rejects_malformed_inputs(self):
        good = genome("G1", [event("2021-01-01")], collection="2020-12-01",
                      released="2021-01-01")
        no_labels = {key: value for key, value in good.items() if key != "labels"}
        no_source = {key: value for key, value in event("2021-01-01").items()
                     if key != "source_release"}
        cases = {
            "slash date": [dict(good, released_at="2021/01/01")],
            "compact date": [dict(good, released_at="20210101")],
            "impossible date": [dict(good, released_at="2021-02-30")],
            "unknown key": [dict(good, features_available="2021-01-01")],
            "missing labels": [no_labels],
            "missing genome_id": [{"labels": []}],
            "labels not a list": [dict(good, labels="R")],
            "record not a mapping": [["G1"]],
            "missing source_release": [dict(good, labels=[no_source])],
            "empty source_release": [dict(good, labels=[event("2021-01-01", source="")])],
            "unknown event key": [dict(good, labels=[dict(event("2021-01-01"), extra=1)])],
            "bad value": [dict(good, labels=[event("2021-01-01", "I")])],
            "empty lineage": [dict(good, lineage="  ")],
            "empty group": [dict(good, duplicate_group=" \t")],
            "undated event supersedes missing day": [dict(good, labels=[
                event("2021-01-01"), event(None, "S", supersedes="2020-06-01")])],
            "supersedes itself": [dict(good, labels=[
                event("2021-01-01", supersedes="2021-01-01")])],
            "supersedes missing day": [dict(good, labels=[
                event("2021-01-01"), event("2021-06-01", "S", supersedes="2020-12-25")])],
            "supersedes other antibiotic": [dict(good, labels=[
                event("2021-01-01", antibiotic="RIF"),
                event("2021-06-01", "S", supersedes="2021-01-01")])],
            "label before collection": [dict(good, labels=[event("2020-11-30")])],
            "duplicate id": [good, good],
        }
        for name, rows in cases.items():
            with self.assertRaises(ValueError, msg=name):
                split(rows)
        for name in ("lineage", "duplicate_group"):
            with self.assertRaisesRegex(ValueError, f"^G1: {name} "):
                split([dict(good, **{name: "   "})])
        self.assertEqual(split([good])["train"], {"G1": "R"})
        for overrides in ({"horizon": 0}, {"horizon": True}, {"embargo": -1},
                          {"embargo": 1.5}, {"eval_date": "2022-06-29"}, {"T": None},
                          {"antibiotic": ""}, {"policy": "collection"},
                          {"T": dt.datetime(2022, 1, 1)}):
            with self.assertRaises(ValueError, msg=overrides):
                split([good], **overrides)
        with self.assertRaises(ValueError):
            SplitPolicy(window="deposit")
        with self.assertRaises(ValueError):
            SplitPolicy(backfill="train")

    def test_accepts_contract_edge_cases(self):
        rows = [genome("empty", []),
                genome("spaced", [event("2021-01-01")], released="2021-01-01",
                       lineage="\tL4.2 "),
                genome("undated_supersedes", [event("2021-01-01"),
                                              event(None, "S", supersedes="2021-01-01")],
                       released="2021-01-01")]
        records = by_id(rows)
        self.assertEqual(records["empty"].labels, ())
        self.assertEqual(records["spaced"].lineage, "L4.2")
        self.assertEqual(split(rows)["train"], {"spaced": "R", "undated_supersedes": "R"})


if __name__ == "__main__":
    unittest.main()
