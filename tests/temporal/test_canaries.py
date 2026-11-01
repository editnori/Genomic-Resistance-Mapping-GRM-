"""canaries.py on synthetic records: every canary passes its clean control
and catches its planted leak, with the offending ids, features or artifacts.

The cutoff is 2026-03-01. Genomes:
    g01 INH S, revised to R after the cutoff
    g02 INH R, never revised
    g03 collected before the cutoff, INH recorded only after it (backfill)
    g04 released and labelled on the cutoff day itself (in training)
    g05 collected before, released after the cutoff (test period)
    g06 released early, features available only after the cutoff
    g07 INH S revised to R before the cutoff; RIF R

usage: python3 -I -B -m unittest discover -s tests/temporal -t .
"""
import collections
import dataclasses
import datetime as dt
import functools
import hashlib
import importlib.util
import itertools
import json
import os
import random
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def load(name):
    """Load scripts/temporal/<name>.py by path under the key canaries.py uses."""
    key = f"grm_temporal_{name}"
    module = sys.modules.get(key)
    if module is None:
        spec = importlib.util.spec_from_file_location(key, ROOT / "scripts" / "temporal" / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return module


A = load("artifacts")
C = load("canaries")
T = "2026-03-01"
CUTOFF = dt.date(2026, 3, 1)


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


RELEASES = [
    {"kind": "catalog", "release_id": "cat-v1", "released_at": "2026-01-01", "sha256": sha("cat-v1"),
     "vocabulary": ["katG", "rpoB", "fabG1"]},
    {"kind": "catalog", "release_id": "cat-v2", "released_at": "2026-03-01", "sha256": sha("cat-v2"),
     "vocabulary": ["katG", "rpoB", "inhA"]},
    {"kind": "catalog", "release_id": "cat-v3", "released_at": "2026-06-01", "sha256": sha("cat-v3"),
     "vocabulary": ["katG", "rpoB", "inhA", "embB"]},
    {"kind": "annotation", "release_id": "ann-1", "released_at": "2025-12-01", "sha256": sha("ann-1"),
     "vocabulary": ["geneA", "geneB"]},
    {"kind": "annotation", "release_id": "ann-2", "released_at": "2026-03-02", "sha256": sha("ann-2"),
     "vocabulary": ["geneA", "geneB", "geneC"]},
]


def event(antibiotic, value, recorded_at, supersedes=None, source="src-1"):
    return {"antibiotic": antibiotic, "value": value, "recorded_at": recorded_at,
            "source_release": source, "supersedes": supersedes}


def genome(gid, collected, released, labels, available=None, lineage="L4", group=None):
    row = {"genome_id": gid, "collection_date": collected, "released_at": released,
           "lineage": lineage, "duplicate_group": group, "labels": labels}
    if available is not None:
        row["features_available_at"] = available
    return row


RECORDS = [
    genome("g01", "2025-06-01", "2025-07-01",
           [event("INH", "S", "2025-08-01"), event("INH", "R", "2026-04-01", "2025-08-01", "src-2")]),
    genome("g02", "2025-05-01", "2025-06-01", [event("INH", "R", "2025-09-01")], lineage="L2"),
    genome("g03", "2025-10-01", "2025-11-01", [event("INH", "R", "2026-05-01", source="src-2")]),
    genome("g04", "2025-11-01", "2026-03-01", [event("INH", "S", "2026-03-01")], group="d1"),
    genome("g05", "2026-02-01", "2026-04-01", [event("INH", "R", "2026-04-02", source="src-2")], group="d1"),
    genome("g06", "2024-12-01", "2025-01-01", [event("INH", "S", "2025-02-01")], available="2026-05-01"),
    genome("g07", "2025-02-01", "2025-03-01",
           [event("INH", "S", "2025-04-01"), event("INH", "R", "2026-01-15", "2025-04-01"),
            event("RIF", "R", "2025-04-01")], lineage=None),
]

FEATURES = {
    "g01": ["katG", "geneA"], "g02": ["katG", "rpoB", "geneB"], "g03": ["inhA"],
    "g04": ["rpoB", "geneA"], "g05": ["embB", "geneC", "katG"], "g06": ["inhA", "geneB"],
    "g07": ["katG", "inhA"],
}


def registry():
    return A.ArtifactRegistry.from_records(RELEASES)


def records():
    return C.load_records(RECORDS)


def clean_training():
    return C.training_labels_as_of(RECORDS, "INH", T)


def clean_history(record):
    """M-1's leak: no undated INH event and no same-day INH conflict anywhere,
    including after the cutoff."""
    events = [e for e in record.labels if e.antibiotic == "INH"]
    days = {}
    for e in events:
        days.setdefault(e.recorded_at, set()).add(e.value)
    return all(e.recorded_at is not None for e in events) and all(len(v) == 1 for v in days.values())


def make_fit(label="as_of", vocabulary="registry", gate=True, noisy=False, history="as_of",
             availability="max"):
    """A small numpy fit; each keyword plants one leak when changed."""
    reg = registry()
    noise = itertools.count()

    def fit(recs, features, cutoff):
        if vocabulary == "registry":
            vocab = sorted(reg.vocabulary(cutoff))
        else:  # planted: vocabulary from every genome's features
            vocab = sorted(frozenset().union(*features.values()))
        index = {gene: i for i, gene in enumerate(vocab)}
        rows = []
        for record in recs:
            if availability == "max":
                available = record.available_at
            else:  # planted: the pre-rule-5 reading, features_available_at overriding released_at
                available = record.features_available_at or record.released_at
            if gate and (available is None or available > cutoff):
                continue
            if history == "all" and not clean_history(record):  # planted: M-1
                continue
            if label == "as_of":
                found = C.label_as_of(record, "INH", cutoff)
            else:  # planted: the latest label ever recorded
                found = ([e for e in record.labels if e.antibiotic == "INH"] or [None])[-1]
            if found is not None:
                rows.append((record.genome_id, found.value))
        x = np.zeros((len(rows), len(vocab)))
        for i, (gid, _) in enumerate(rows):
            for gene in features.get(gid, ()):
                if gene in index:
                    x[i, index[gene]] = 1.0
        y = np.array([value == "R" for _, value in rows], dtype=float)
        bias = float(y.mean()) if rows else 0.0
        if noisy:
            bias += next(noise) * 1e-12
        weights = x.T @ (y - bias) / max(len(rows), 1)
        selected = [vocab[j] for j in range(len(vocab)) if x[:, j].any()]
        return {"vocabulary": vocab, "encoder": index, "selected_features": selected,
                "model": {"weights": weights, "bias": bias}, "training_genomes": [gid for gid, _ in rows]}
    # The registry vocabulary at T does not depend on the training rows (NN-1).
    fit.constant_keys = ("vocabulary", "encoder") if vocabulary == "registry" else ()
    return fit


class Records(unittest.TestCase):
    def test_load_and_defaults(self):
        recs = records()
        self.assertEqual([r.genome_id for r in recs], sorted(FEATURES))
        g01 = recs[0]
        self.assertIsNone(g01.features_available_at)
        self.assertEqual(g01.available_at, dt.date(2025, 7, 1))
        self.assertEqual(recs[5].available_at, dt.date(2026, 5, 1))
        self.assertEqual(g01.labels[1].supersedes, dt.date(2025, 8, 1))
        # Missing nullable fields are null.
        bare = C.load_record({"genome_id": "x", "labels": []})
        self.assertIsNone(bare.released_at)
        self.assertFalse(C.knowable(bare, T))

    def test_label_as_of_boundaries(self):
        recs = {r.genome_id: r for r in records()}
        self.assertEqual(C.label_as_of(recs["g04"], "INH", "2026-03-01").value, "S")
        self.assertIsNone(C.label_as_of(recs["g04"], "INH", "2026-02-28"))
        self.assertEqual(C.label_as_of(recs["g01"], "INH", "2026-03-31").value, "S")
        self.assertEqual(C.label_as_of(recs["g01"], "INH", "2026-04-01").value, "R")
        self.assertEqual(C.label_as_of(recs["g07"], "INH", T).value, "R")
        self.assertIsNone(C.label_as_of(recs["g07"], "EMB", T))

    def test_label_order_is_canonical(self):
        row = dict(RECORDS[6], labels=list(reversed(RECORDS[6]["labels"])))
        self.assertEqual(C.load_record(row), C.load_record(RECORDS[6]))

    def test_contract_violations_are_refused(self):
        no_source = {k: v for k, v in event("INH", "R", "2025-09-01").items() if k != "source_release"}
        no_recorded_key = {k: v for k, v in event("INH", "R", "2025-09-01").items() if k != "recorded_at"}
        bad = [
            dict(RECORDS[1], extra=1),
            {k: v for k, v in RECORDS[1].items() if k != "labels"},  # rule 2
            dict(RECORDS[1], labels="INH R"),
            dict(RECORDS[1], genome_id=""),
            dict(RECORDS[1], released_at="2025/06/01"),
            dict(RECORDS[1], lineage=4),
            dict(RECORDS[1], lineage="   "),
            dict(RECORDS[1], labels=[event("INH", "I", "2025-09-01")]),
            dict(RECORDS[1], labels=[no_recorded_key]),
            dict(RECORDS[1], labels=[event("INH", "R", "2025-09-01", source="")]),  # rule 2
            dict(RECORDS[1], labels=[no_source]),  # rule 2
            dict(RECORDS[1], labels=[dict(event("INH", "R", "2025-09-01"), note="x")]),
            dict(RECORDS[1], labels=[event("INH", "R", "2025-09-01", "2025-08-01")]),  # rule 4
            dict(RECORDS[1], labels=[event("INH", "R", "2025-09-01", "2025-09-01")]),
            dict(RECORDS[1], labels=[event("INH", "R", "2025-09-01"),
                                     event("RIF", "S", "2025-10-01", "2025-09-01")]),  # rule 4
            dict(RECORDS[1], labels=[event("INH", "R", None), event("INH", "S", "2025-10-01", "2025-09-01")]),
            dict(RECORDS[1], labels=[event("INH", "R", "2025-04-30")]),  # rule 5: before collection
        ]
        for row in bad:
            with self.subTest(row=row), self.assertRaises(C.RecordError):
                C.load_record(row)
        with self.assertRaisesRegex(C.RecordError, "duplicate genome_id"):
            C.load_records(RECORDS + [RECORDS[0]])


class ContractRules(unittest.TestCase):
    """The five binding rulings, one test per rule."""

    def test_rule1_undated_events_load_and_are_never_knowable(self):
        rec = C.load_record(dict(RECORDS[1], labels=[event("INH", "S", None), event("INH", "R", "2025-09-01"),
                                                     event("RIF", "R", None)]))
        self.assertEqual(C.label_as_of(rec, "INH", "2099-01-01").value, "R")
        self.assertIsNone(C.label_as_of(rec, "RIF", "2099-01-01"))
        self.assertEqual(C.training_labels_as_of([rec], "RIF", T), [])
        # An undated event may revise a dated one; it still never counts.
        rec = C.load_record(dict(RECORDS[1], labels=[event("INH", "R", "2025-09-01"),
                                                     event("INH", "S", None, "2025-09-01")]))
        self.assertEqual(C.label_as_of(rec, "INH", "2099-01-01").value, "R")

    def test_rule2_empty_labels_load(self):
        rec = C.load_record(dict(RECORDS[1], labels=[]))
        self.assertEqual(rec.labels, ())
        self.assertIsNone(C.label_as_of(rec, "INH", T))

    def test_rule3_same_day_events(self):
        agree = [event("INH", "R", "2025-09-01", source="src-b"), event("INH", "R", "2025-09-01", source="src-a")]
        rec = C.load_record(dict(RECORDS[1], labels=agree))
        self.assertEqual(len(rec.labels), 2)
        self.assertEqual(C.label_as_of(rec, "INH", T).source_release, "src-a")
        conflict = [event("INH", "R", "2025-09-01", source="src-a"), event("INH", "S", "2025-09-01", source="src-b")]
        rec = C.load_record(dict(RECORDS[1], labels=conflict))
        self.assertIsNone(C.label_as_of(rec, "INH", "2025-08-31"))
        with self.assertRaises(C.AmbiguousLabel):
            C.label_as_of(rec, "INH", "2025-09-01")
        self.assertEqual(C.training_labels_as_of([rec], "INH", T), [])
        # A later dated event in the window supersedes the conflict.
        rec = C.load_record(dict(RECORDS[1], labels=conflict + [event("INH", "S", "2025-12-01", "2025-09-01")]))
        with self.assertRaises(C.AmbiguousLabel):
            C.label_as_of(rec, "INH", "2025-11-30")
        self.assertEqual(C.label_as_of(rec, "INH", "2025-12-01").value, "S")
        # A conflict after the cutoff never reaches the train-side read.
        rec = C.load_record(dict(RECORDS[1], labels=[event("INH", "R", "2025-09-01")] + [
            event("INH", v, "2026-05-01", source=f"src-{v}") for v in "RS"]))
        self.assertEqual(C.training_labels_as_of([rec], "INH", T)[0]["value"], "R")
        with self.assertRaises(C.AmbiguousLabel):
            C.label_as_of(rec, "INH", "2026-05-01")

    def test_rule5_availability_is_the_later_date(self):
        late_release = C.load_record(dict(RECORDS[1], released_at="2026-04-01", features_available_at="2026-01-01"))
        self.assertEqual(late_release.available_at, dt.date(2026, 4, 1))
        self.assertTrue(late_release.inconsistent_dates)
        self.assertFalse(C.knowable(late_release, T))
        late_features = C.load_record(dict(RECORDS[1], features_available_at="2026-04-01"))
        self.assertEqual(late_features.available_at, dt.date(2026, 4, 1))
        self.assertFalse(late_features.inconsistent_dates)
        only_features = C.load_record(dict(RECORDS[1], released_at=None, features_available_at="2025-07-01"))
        self.assertEqual(only_features.available_at, dt.date(2025, 7, 1))
        self.assertEqual(C.inconsistent_dates([late_release, dict(RECORDS[0])]), ["g02"])
        # Labels on the collection day load; undated labels skip the check.
        C.load_record(dict(RECORDS[1], labels=[event("INH", "R", "2025-05-01"), event("RIF", "R", None)]))

    def test_rule6_lineage_and_duplicate_group_are_stripped(self):
        rec = C.load_record(dict(RECORDS[1], lineage=" 4.2 ", duplicate_group="D1 "))
        self.assertEqual((rec.lineage, rec.duplicate_group), ("4.2", "D1"))
        self.assertEqual(C.load_record(dict(RECORDS[1], duplicate_group="\td1\n")).duplicate_group, "d1")  # case kept
        # Rule 6 v2: str.strip() first, so characters it counts as whitespace go
        # (MN-1); interior spaces and every printable ASCII character stay.
        for raw, kept in (("\x1fL2", "L2"), ("L2\x85", "L2"), ("\x0bL2\x0c", "L2"), ("L 4.2", "L 4.2"),
                          ("!~", "!~"), ("a-b_c/d", "a-b_c/d")):
            with self.subTest(raw=raw):
                self.assertEqual(C.load_record(dict(RECORDS[1], lineage=raw)).lineage, kept)
        for field in ("lineage", "duplicate_group"):
            for bad in ("", "   ", "\t\n", "D\t1", "D1​", "​D1", "D­1", "D1\x00", "L\x1f2",
                        "D1\x7f", "Lé1", "Ｄ１", "D 1", "D 1", "Ｄ1"):
                with self.subTest(field=field, value=bad), self.assertRaisesRegex(
                        ValueError, f"g02 {field}: (expected a non-empty string|.* outside printable ASCII)"):
                    C.load_record(dict(RECORDS[1], **{field: bad}))

    def test_rule6_alias_collisions_in_one_input(self):
        def rows(*values, field="duplicate_group"):
            key = "lineage" if field == "lineage" else "group"
            return [genome(f"a{i}", "2025-01-01", "2025-02-01", [], **{key: v}) for i, v in enumerate(values)]
        # 'D1 ' strips to 'D1' itself, so it is the same value, not an alias.
        self.assertEqual({r.duplicate_group for r in C.load_records(rows("D1", "D1 ", " D1"))}, {"D1"})
        for field in ("duplicate_group", "lineage"):
            # Case is the alias rule 6 v2 leaves to the collision check; full-width,
            # a zero-width space and a tab are refused per genome before it.
            for alias, reason in (("d1", f"{field} alias collision: 'D1' \\(a0\\) and 'd1' \\(a1\\)"),
                                  ("Ｄ１", f"a1 {field}: .* outside printable ASCII"),
                                  ("D1​", f"a1 {field}: .* outside printable ASCII"),
                                  ("D\t1", f"a1 {field}: .* outside printable ASCII")):
                with self.subTest(field=field, alias=alias), self.assertRaisesRegex(ValueError, reason):
                    C.load_records(rows("D1", alias, field=field))
            with self.assertRaisesRegex(ValueError, f"{field} alias collision: 'L 4' \\(a0\\) and 'l 4' \\(a1\\)"):
                C.load_records(rows("L 4", " l 4 ", field=field))
        # Groups and lineages are separate namespaces; distinct values pass.
        C.load_records([genome("a0", "2025-01-01", "2025-02-01", [], lineage="D1", group="d1"),
                        genome("a1", "2025-01-01", "2025-02-01", [], lineage="L2", group="D2")])
        # Every canary reads its records through load_records.
        with self.assertRaisesRegex(ValueError, "alias collision"):
            C.backfill(rows("D1", "d1"), T, [])

    def test_ruling_b_tie_break_ignores_input_order(self):
        # Same day, same value, same source_release; only supersedes differs.
        history = [event("INH", "S", "2025-06-01"), event("INH", "S", "2025-07-01"),
                   event("INH", "R", "2025-09-01", "2025-07-01", "src-a"),
                   event("INH", "R", "2025-09-01", "2025-06-01", "src-a"),
                   event("INH", "R", "2025-09-01", None, "src-a"), event("INH", "R", "2025-09-01", None, "src-b")]
        rng = random.Random(7)
        for orders in range(12):
            labels = list(reversed(history)) if orders == 0 else rng.sample(history, len(history))
            got = C.label_as_of(C.load_record(dict(RECORDS[1], labels=labels)), "INH", T)
            self.assertEqual((got.value, got.source_release, got.supersedes), ("R", "src-a", None))
        without_null = [e for e in history if e["supersedes"] is not None or e["value"] == "S"]
        for labels in (without_null, list(reversed(without_null))):
            got = C.label_as_of(C.load_record(dict(RECORDS[1], labels=labels)), "INH", T)
            self.assertEqual((got.source_release, got.supersedes), ("src-a", dt.date(2025, 6, 1)))

    def test_reference_training_labels(self):
        self.assertEqual([(r["genome_id"], r["value"]) for r in clean_training()],
                         [("g01", "S"), ("g02", "R"), ("g04", "S"), ("g07", "R")])


class FutureMarker(unittest.TestCase):
    def test_clean(self):
        used = sorted(registry().vocabulary(T)) + [A.UNKNOWN]
        result = C.future_marker(registry(), T, used)
        self.assertTrue(result.passed, result)
        self.assertEqual(result.offenders, ())

    def test_leak(self):
        result = C.future_marker(registry(), T, ["katG", "embB", "geneC", "fabG1", "madeUp"])
        self.assertFalse(result.passed)
        self.assertEqual(result.offenders, (
            ("embB", "catalog:cat-v3@2026-06-01", "only in a release after the cutoff"),
            ("fabG1", "catalog:cat-v1@2026-01-01", "retired before the cutoff"),
            ("geneC", "annotation:ann-2@2026-03-02", "only in a release after the cutoff"),
            ("madeUp", "-", "in no registered release"),
        ))

    def test_feature_is_clean_from_its_release_day(self):
        self.assertFalse(C.future_marker(registry(), "2026-03-01", ["geneC"]).passed)
        self.assertTrue(C.future_marker(registry(), "2026-03-02", ["geneC"]).passed)

    def test_n6_generator_input_is_counted(self):
        result = C.future_marker(registry(), T, (f for f in ["katG", "rpoB", "katG"]))
        self.assertEqual(result.detail, "all 2 feature(s) are in the vocabulary as of 2026-03-01")
        leak = C.future_marker(registry(), T, (f for f in ["katG", "embB"]))
        self.assertEqual([o[0] for o in leak.offenders], ["embB"])

    def test_unresolvable_artifacts_fail(self):
        reg = registry()
        reg.add({"kind": "catalog", "release_id": "cat-v4", "released_at": "2026-02-01", "sha256": None,
                 "vocabulary": ["katG"]})
        result = C.future_marker(reg, "2026-02-15", ["katG"])
        self.assertFalse(result.passed)
        self.assertIn("catalog:cat-v4@2026-02-01: sha256 is missing", result.offenders[0][2])


class ArtifactManifest(unittest.TestCase):
    def test_clean(self):
        result = C.artifact_manifest(registry(), T, registry().manifest(T))
        self.assertTrue(result.passed, result)

    def test_leak(self):
        declared = registry().manifest(T)
        declared["artifacts"]["catalog"] = registry().manifest("2026-06-01")["artifacts"]["catalog"]
        result = C.artifact_manifest(registry(), T, declared)
        self.assertFalse(result.passed)
        self.assertEqual({where for where, _ in result.offenders}, {"catalog"})
        self.assertIn("release_id declared 'cat-v3', resolved 'cat-v2'", [p for _, p in result.offenders])


class Backfill(unittest.TestCase):
    def test_clean(self):
        result = C.backfill(RECORDS, T, clean_training())
        self.assertTrue(result.passed, result)

    def test_leak(self):
        leaked = clean_training() + [{"genome_id": "g03", "antibiotic": "INH", "value": "R"}]
        result = C.backfill(RECORDS, T, leaked)
        self.assertFalse(result.passed)
        self.assertEqual(result.offenders, (
            ("g03", "INH", "first recorded 2026-05-01, after the cutoff; collected 2025-10-01"),))

    def test_other_unknown_labels(self):
        rows = [("g07", "EMB", "R"), ("nope", "INH", "S")]
        self.assertEqual(C.backfill(RECORDS, T, rows).offenders, (
            ("g07", "EMB", "no label event at all"), ("nope", "INH", "genome not in the records")))

    def test_undated_label_leak(self):
        recs = RECORDS + [genome("g08", "2025-01-01", "2025-02-01", [event("INH", "R", None)])]
        self.assertEqual(C.backfill(recs, T, [("g08", "INH", "R")]).offenders, (
            ("g08", "INH", "only undated label events, never knowable"),))

    def test_n1_unknowable_genome_leak(self):
        # The reviewer's B3a: x2 released 2026-04-01, INH recorded 2026-02-20, T 2026-03-01.
        clean = C.training_labels_as_of(REVIEW_RECORDS, "INH", T)
        self.assertTrue(C.backfill(REVIEW_RECORDS, T, clean).passed)
        leak = clean + [{"genome_id": "x2", "antibiotic": "INH", "value": "R", "recorded_at": "2026-02-20"}]
        self.assertEqual(C.backfill(REVIEW_RECORDS, T, leak).offenders, (
            ("x2", "INH", "features available 2026-04-01, after the cutoff"),))
        never = REVIEW_RECORDS + [genome("x3", "2025-01-01", None, [event("INH", "R", "2025-02-01")])]
        self.assertEqual(C.backfill(never, T, [("x3", "INH", "R")]).offenders, (
            ("x3", "INH", "features available never, after the cutoff"),))

    def test_ambiguous_label_is_not_backfill(self):
        recs = RECORDS + [genome("g08", "2025-01-01", "2025-02-01",
                                 [event("INH", v, "2025-03-01", source=f"src-{v}") for v in "RS"])]
        self.assertTrue(C.backfill(recs, T, [("g08", "INH", "R")]).passed)

    def test_bad_rows_are_refused(self):
        for row in ({"genome_id": "g01", "antibiotic": "INH"}, ("g01", "INH", "X"),
                    {"genome_id": "g01", "antibiotic": "INH", "value": "S", "weight": 1}):
            with self.subTest(row=row), self.assertRaises(C.RecordError):
                C.backfill(RECORDS, T, [row])


class LabelRevision(unittest.TestCase):
    def test_clean(self):
        result = C.label_revision(RECORDS, T, clean_training())
        self.assertTrue(result.passed, result)
        # A revision recorded before the cutoff is the label as of the cutoff.
        self.assertTrue(C.label_revision(RECORDS, T, [("g07", "INH", "R", "2026-01-15")]).passed)

    def test_leak(self):
        leaked = [dict(row, value="R") if row["genome_id"] == "g01" else row for row in clean_training()]
        result = C.label_revision(RECORDS, T, leaked)
        self.assertFalse(result.passed)
        self.assertEqual(result.offenders, (
            ("g01", "INH", "R from the revision recorded 2026-04-01; as of the cutoff S, recorded 2025-08-01"),))

    def test_declared_later_event_and_corrupt_value(self):
        rows = [("g02", "INH", "R", "2026-04-01"), ("g02", "INH", "S")]
        self.assertEqual(C.label_revision(RECORDS, T, rows).offenders, (
            ("g02", "INH", "S is in no later revision; as of the cutoff R, recorded 2025-09-01"),
            ("g02", "INH", "declares the event recorded 2026-04-01; as of the cutoff R, recorded 2025-09-01"),
        ))

    def test_backfilled_rows_are_left_to_backfill(self):
        self.assertTrue(C.label_revision(RECORDS, T, [("g03", "INH", "R")]).passed)

    def test_n1_unknowable_genome_leak(self):
        # The reviewer's B3b: the same x2 row, which matches its label as of T.
        clean = C.training_labels_as_of(REVIEW_RECORDS, "INH", T)
        self.assertTrue(C.label_revision(REVIEW_RECORDS, T, clean).passed)
        leak = clean + [("x2", "INH", "R", "2026-02-20")]
        self.assertEqual(C.label_revision(REVIEW_RECORDS, T, leak).offenders, (
            ("x2", "INH", "features available 2026-04-01, after the cutoff"),))
        # g06: label recorded 2025-02-01, features only from 2026-05-01.
        self.assertEqual(C.label_revision(RECORDS, T, [("g06", "INH", "S")]).offenders, (
            ("g06", "INH", "features available 2026-05-01, after the cutoff"),))

    def test_ambiguous_and_undated_leaks(self):
        recs = RECORDS + [
            genome("g08", "2025-01-01", "2025-02-01", [event("INH", v, "2025-03-01", source=f"src-{v}") for v in "RS"]),
            genome("g09", "2025-01-01", "2025-02-01", [event("INH", "S", "2025-03-01"), event("INH", "R", None)]),
        ]
        self.assertEqual(C.label_revision(recs, T, [("g08", "INH", "R"), ("g09", "INH", "R")]).offenders, (
            ("g08", "INH", "R and S both recorded 2025-03-01, the latest day by the cutoff"),
            ("g09", "INH", "R from an undated event, never knowable; as of the cutoff S, recorded 2025-03-01"),
        ))

    def test_conflict_after_the_cutoff_leaves_training_alone(self):
        # M-1: same-day conflicts and undated events after T change neither rows nor verdicts.
        recs = [dict(row, labels=row["labels"] + [event("INH", v, "2026-05-01", source=f"tie-{v}") for v in "RS"]
                     + [event("INH", "R", None, source="undated")]) for row in RECORDS]
        self.assertEqual(C.training_labels_as_of(recs, "INH", T), clean_training())
        self.assertTrue(C.label_revision(recs, T, clean_training()).passed)
        self.assertTrue(C.backfill(recs, T, clean_training()).passed)


# The reviewer's fixture (grm-review-task-n-20261030/probes.py, B section):
# t1-t4 train at T; x1 public before T but labelled after it (ruling A: not a
# training row); x2 released after T with a label recorded before T, sharing
# t2's lineage and duplicate group.
REVIEW_RECORDS = [
    genome("t1", "2025-10-01", "2025-11-01", [event("INH", "S", "2025-12-01")], lineage="L2", group="dgA"),
    genome("t2", "2025-10-01", "2025-11-01", [event("INH", "R", "2025-12-15")], lineage="L4", group="dgB"),
    genome("t3", "2025-10-01", "2025-11-01",
           [event("INH", "S", "2026-01-10"), event("INH", "R", "2026-04-01", "2026-01-10", "s2")],
           lineage="L1", group="dgC"),
    genome("t4", "2025-10-01", "2025-11-01", [event("INH", "R", "2026-01-20")], lineage="L3", group="dgD"),
    genome("x1", "2026-01-01", "2026-02-01", [event("INH", "R", "2026-05-01")], lineage="L5", group="dgE"),
    genome("x2", "2026-01-15", "2026-04-01", [event("INH", "R", "2026-02-20")], lineage="L4", group="dgB"),
]
REVIEW_FEATURES = {"t1": ["katG", "rpoB"], "t2": ["katG"], "t3": ["inhA"], "t4": ["katG", "inhA"],
                   "x1": ["katG", "embB"], "x2": ["katG", "pncA"]}


def honest(recs, feats, cutoff):
    return {"train": sorted((r["genome_id"], r["value"]) for r in C.training_labels_as_of(recs, "INH", cutoff))}


def vocab_all(recs, feats, cutoff):
    return {"vocab": sorted(set().union(*feats.values()))}


def vocab_knowable(recs, feats, cutoff):
    """Ruling A leak: vocabulary over every knowable genome, labelled or not."""
    return {"vocab": sorted(set().union(*(feats.get(r.genome_id, ()) for r in recs if C.knowable(r, cutoff))))}


def vocab_training(recs, feats, cutoff):
    """Clean under ruling A: vocabulary over the training rows only."""
    train = {r["genome_id"] for r in C.training_labels_as_of(recs, "INH", cutoff)}
    return {"vocab": sorted(set().union(*(feats.get(gid, ()) for gid in train)))}


def make_cached(inner):
    """The reviewer's B1b: memoise on the cutoff."""
    cache = {}

    def fit(recs, feats, cutoff):
        if cutoff not in cache:
            cache[cutoff] = inner(recs, feats, cutoff)
        return cache[cutoff]
    return fit


def make_lru_leak():
    """The reviewer's B1c: functools.lru_cache keyed on the cutoff."""
    last = {}

    @functools.lru_cache(maxsize=None)
    def cached_vocab(cutoff):
        return tuple(sorted(set().union(*last["feats"].values())))

    def fit(recs, feats, cutoff):
        last["feats"] = feats
        return {"vocab": list(cached_vocab(cutoff))}
    return fit


def make_dedup(field, with_values=True):
    """The reviewer's B2a/B2b: drop training genomes that share a lineage or
    duplicate group with a genome outside the training rows. with_values=False
    is the reviewer's exact fit (ids only), which the positive control rejects."""
    def fit(recs, feats, cutoff):
        rows = C.training_labels_as_of(recs, "INH", cutoff)
        train = {r["genome_id"] for r in rows}
        others = {getattr(r, field) for r in recs if r.genome_id not in train} - {None}
        by = {r.genome_id: r for r in recs}
        kept = [r for r in rows if getattr(by[r["genome_id"]], field) not in others]
        return {"train": sorted((r["genome_id"], r["value"]) if with_values else r["genome_id"] for r in kept)}
    return fit


class TestDataInvariance(unittest.TestCase):
    def run_canary(self, fit, **kwargs):
        kwargs.setdefault("constant_keys", getattr(fit, "constant_keys", ()))
        return C.invariance_to_test_data(RECORDS, FEATURES, T, fit, antibiotics=("INH",), **kwargs)

    def run_review(self, fit, **kwargs):
        return C.invariance_to_test_data(REVIEW_RECORDS, REVIEW_FEATURES, T, fit, antibiotics=("INH",), **kwargs)

    def test_clean(self):
        result = self.run_canary(make_fit())
        self.assertTrue(result.passed, result)
        self.assertEqual(make_fit()(records(), C._features(FEATURES), CUTOFF)["training_genomes"],
                         ["g01", "g02", "g04", "g07"])
        self.assertEqual(C.training_ids(RECORDS, T, ("INH",)), {"g01", "g02", "g04", "g07"})
        for fit in (honest, vocab_training):
            with self.subTest(fit.__name__):
                self.assertTrue(self.run_review(fit).passed, self.run_review(fit))

    def test_latest_label_leak(self):
        result = self.run_canary(make_fit(label="latest"))
        self.assertFalse(result.passed)
        # Its g03 row (public at T, labelled after it) also breaks ruling A, so
        # features and new_genomes fire too.
        self.assertTrue({"labels", "label_ties", "features", "new_genomes"} <= {o[0] for o in result.offenders})
        self.assertIn(("labels", "model", "changed by test-period records"), result.offenders)

    def test_m1_history_filter_leak(self):
        # Dropping a training genome over a conflict or undated event recorded
        # after T (M-1) is invisible to the uniform label flips.
        result = self.run_canary(make_fit(history="all"))
        self.assertFalse(result.passed)
        self.assertEqual({o[0] for o in result.offenders}, {"label_ties", "undated_labels"})
        self.assertIn(("label_ties", "training_genomes", "changed by test-period records"), result.offenders)
        self.assertIn(("undated_labels", "training_genomes", "changed by test-period records"), result.offenders)

    def test_pre_rule5_availability_leak(self):
        # Only the added genome with features dated before its later release
        # (inconsistent_dates) separates the two readings in this fixture.
        result = self.run_canary(make_fit(availability="features_first"))
        self.assertFalse(result.passed)
        self.assertEqual({o[0] for o in result.offenders}, {"new_genomes"})
        self.assertIn(("new_genomes", "training_genomes", "changed by test-period records"), result.offenders)

    def test_vocabulary_from_test_genomes_leak(self):
        result = self.run_canary(make_fit(vocabulary="all"))
        self.assertFalse(result.passed)
        changed = {(o[0], o[1]) for o in result.offenders}
        self.assertTrue({("features", "vocabulary"), ("new_genomes", "vocabulary"),
                         ("features", "encoder"), ("new_genomes", "encoder")} <= changed, changed)
        self.assertNotIn("labels", {o[0] for o in result.offenders})

    def test_ruling_a_knowable_but_unlabelled_genomes_leak(self):
        # x1 was public at T but labelled after it: its genes must not reach the vocabulary.
        result = self.run_review(vocab_knowable)
        self.assertFalse(result.passed)
        self.assertIn(("features", "vocab", "changed by test-period records"), result.offenders)
        self.assertIn(("new_genomes", "vocab", "changed by test-period records"), result.offenders)

    def test_ungated_genomes_leak(self):
        result = self.run_canary(make_fit(gate=False))
        self.assertFalse(result.passed)
        changed = {(o[0], o[1]) for o in result.offenders}
        # g06 (features after the cutoff) is trained on, so its perturbed genes
        # move the model; the added genomes carry labels recorded by the cutoff.
        self.assertTrue({("features", "model"), ("new_genomes", "training_genomes")} <= changed, changed)

    def test_nondeterministic_fit_is_not_blamed_on_test_data(self):
        result = self.run_canary(make_fit(noisy=True))
        self.assertFalse(result.passed)
        self.assertEqual(result.offenders, (("refit", "model", "differs on unchanged inputs"),))

    def test_n2_cached_leaks_fail_the_positive_control(self):
        self.assertFalse(self.run_review(vocab_all).passed)

        def insensitive(*keys):
            return tuple(("control", key, "insensitive to its own training data; not certified") for key in keys)
        for name, fit, keys in (("dict cache", make_cached(vocab_all), ["vocab"]),
                                ("lru_cache", make_lru_leak(), ["vocab"]),
                                ("cached clean fit", make_cached(make_fit()),
                                 ["encoder", "model", "selected_features", "training_genomes", "vocabulary"]),
                                ("constant", lambda *_: {"k": 1}, ["k"]), ("constant non-mapping", lambda *_: 1,
                                                                           ["<fit>"])):
            with self.subTest(name):
                self.assertEqual(self.run_review(fit).offenders, insensitive(*keys))

    def test_nn1_control_applies_per_key(self):
        # The re-checker's G1: honest training-row keys plus one key computed
        # from every genome and cached on the cutoff.
        def make_one_key_cached():
            cache = {}

            def fit(recs, feats, cutoff):
                out = vocab_training(recs, feats, cutoff)
                out.update(honest(recs, feats, cutoff))
                if cutoff not in cache:
                    cache[cutoff] = sorted(set().union(*feats.values()))
                out["screened_genes"] = cache[cutoff]
                return out
            return fit
        result = self.run_review(make_one_key_cached())
        self.assertEqual(result.offenders,
                         (("control", "screened_genes", "insensitive to its own training data; not certified"),))
        # Uncached, the same key is caught as a leak by the perturbations instead.
        def uncached(recs, feats, cutoff):
            return dict(make_one_key_cached()(recs, feats, cutoff))
        offenders = self.run_review(uncached).offenders
        self.assertIn(("features", "screened_genes", "changed by test-period records"), offenders)
        self.assertNotIn("control", {o[0] for o in offenders})
        # Each key needs its own control: n_train moves only when a row is dropped.
        def counted(recs, feats, cutoff):
            return dict(vocab_training(recs, feats, cutoff), n_train=len(honest(recs, feats, cutoff)["train"]))
        self.assertTrue(self.run_review(counted).passed, self.run_review(counted))
        with mock.patch.dict(C.TRAINING_CONTROLS):
            del C.TRAINING_CONTROLS["training_rows"]
            offenders = self.run_review(counted).offenders
        self.assertEqual(offenders, (("control", "n_train", "insensitive to its own training data; not certified"),))

    def test_nn1_constant_keys(self):
        def with_constants(recs, feats, cutoff):
            return dict(vocab_training(recs, feats, cutoff), schema=2, registry=sorted(registry().vocabulary(cutoff)))
        self.assertEqual({o[1] for o in self.run_review(with_constants).offenders}, {"schema", "registry"})
        self.assertTrue(self.run_review(with_constants, constant_keys=("schema", "registry")).passed)
        result = self.run_review(with_constants, constant_keys=["schema", "registry", "vocab", "absent"])
        self.assertEqual(result.offenders,
                         (("control", "absent", "named in constant_keys but not fitted"),
                          ("control", "vocab", "named in constant_keys but changed by training data")))
        # A declared constant is still compared under every perturbation.
        def leaky_constant(recs, feats, cutoff):
            return dict(vocab_training(recs, feats, cutoff), n_genomes=len(recs))
        self.assertIn(("new_genomes", "n_genomes", "changed by test-period records"),
                      self.run_review(leaky_constant, constant_keys=("n_genomes",)).offenders)
        with self.assertRaisesRegex(ValueError, "constant_keys"):
            self.run_review(with_constants, constant_keys="schema")

    def test_nn1_registry_vocabulary_is_a_declared_constant(self):
        # The re-checker's C2: training-row keys plus the registry vocabulary at T
        # and each training genome encoded against it.
        def registry_fit(recs, feats, cutoff):
            vocab = registry().vocabulary(cutoff)
            train = {r["genome_id"] for r in C.training_labels_as_of(recs, "INH", cutoff)}
            return dict(honest(recs, feats, cutoff), registry_vocabulary=sorted(vocab),
                        encoded={gid: A.encode_genes(feats.get(gid, ()), vocab) for gid in sorted(train)})
        self.assertEqual(self.run_review(registry_fit).offenders,
                         (("control", "registry_vocabulary", "insensitive to its own training data; not certified"),))
        result = self.run_review(registry_fit, constant_keys=["registry_vocabulary"])
        self.assertTrue(result.passed, result)

    def test_nn2_dates_of_test_genomes(self):
        # The re-checker's G2 and G3: dates of genomes outside the training rows.
        def age_vs_newest(recs, feats, cutoff):
            train = {r["genome_id"] for r in C.training_labels_as_of(recs, "INH", cutoff)}
            newest = max(r.released_at for r in recs if r.released_at)
            return dict(honest(recs, feats, cutoff),
                        age={r.genome_id: (newest - r.released_at).days for r in recs if r.genome_id in train})

        def scaler(recs, feats, cutoff):
            days = [r.available_at.toordinal() for r in recs if r.available_at]
            return dict(honest(recs, feats, cutoff), scaler=[min(days), max(days)])

        def latest_collection(recs, feats, cutoff):
            return dict(honest(recs, feats, cutoff),
                        last=max((r.collection_date for r in recs if r.collection_date), default=None))
        def own_dates_vs_cutoff(recs, feats, cutoff):
            # Clean: each training genome's own dates, measured against T.
            train = {r["genome_id"] for r in C.training_labels_as_of(recs, "INH", cutoff)}
            return dict(honest(recs, feats, cutoff), own={
                r.genome_id: [None if d is None else (cutoff - d).days
                              for d in (r.released_at, r.features_available_at, r.available_at, r.collection_date)]
                for r in recs if r.genome_id in train})
        for fixture in ((REVIEW_RECORDS, REVIEW_FEATURES), (RECORDS, FEATURES)):
            result = C.invariance_to_test_data(*fixture, T, own_dates_vs_cutoff, antibiotics=("INH",))
            self.assertTrue(result.passed, result)
        for fit, key in ((age_vs_newest, "age"), (scaler, "scaler")):
            with self.subTest(fit.__name__):
                offenders = self.run_review(fit).offenders
                self.assertIn(("dates", key, "changed by test-period records"), offenders)
                self.assertIn(("new_genomes", key, "changed by test-period records"), offenders)
        recs, feats = C.load_records(REVIEW_RECORDS), C._features(REVIEW_FEATURES)
        training = C.training_ids(recs, CUTOFF, ("INH",))
        moved = {r.genome_id: r for r in C.perturb_dates(recs, feats, CUTOFF, training)[0]}
        before = {r.genome_id: r for r in recs}
        for gid in training:  # training rows keep every date, collection_date included
            self.assertEqual(moved[gid], before[gid])
            self.assertIsNotNone(moved[gid].collection_date)
        self.assertEqual({gid for gid in moved if moved[gid] != before[gid]}, {"x1", "x2"})
        # x1 (public at T, labelled after) moves earlier; x2 (released after T) moves later.
        self.assertEqual((moved["x1"].released_at, moved["x1"].collection_date), (dt.date(2024, 12, 28), None))
        self.assertEqual(moved["x2"].released_at, dt.date(2027, 5, 6))
        self.assertTrue(C.knowable(moved["x1"], CUTOFF))
        self.assertFalse(C.knowable(moved["x2"], CUTOFF))
        self.assertIn(("dates", "last", "changed by test-period records"), self.run_review(latest_collection).offenders)
        # g06: released before T, features only after it; each date moves away from T on its own side.
        g06 = {r.genome_id: r for r in C.perturb_dates(records(), {}, CUTOFF, frozenset())[0]}["g06"]
        self.assertEqual((g06.released_at, g06.features_available_at), (dt.date(2023, 11, 28), dt.date(2027, 6, 5)))
        self.assertFalse(C.knowable(g06, CUTOFF))
        self.assertEqual((C._shift(dt.date.max, CUTOFF), C._shift(dt.date(1, 1, 2), CUTOFF), C._shift(None, CUTOFF)),
                         (dt.date.max, dt.date.min, None))

    def test_nn5_ambiguous_label_is_not_a_training_row(self):
        # The re-checker's t4: its only INH day at T holds R and S.
        rows = REVIEW_RECORDS + [genome("t9", "2025-10-01", "2025-11-01",
                                        [event("INH", "R", "2026-01-05"), event("INH", "S", "2026-01-05", source="s0")],
                                        lineage="L9", group="dgZ")]
        recs = C.load_records(rows)
        t9 = {r.genome_id: r for r in recs}["t9"]
        self.assertTrue(C.knowable(t9, CUTOFF))
        self.assertFalse(C.is_training_row(t9, CUTOFF, ["INH"]))
        self.assertEqual(C.training_ids(rows, T, ("INH",)), {"t1", "t2", "t3", "t4"})
        # A fit that reads t9 through labels recorded by T only (TL-2: the
        # re-check's NL6 shape, the conflict counted as R) is caught by
        # drop_non_training and training_view_only, and by nothing else.
        def prior_ties_as_r(recs, feats, cutoff):
            values = []
            for r in recs:
                known = C.known_events(r, "INH", cutoff) if C.knowable(r, cutoff) else []
                if known:
                    day = known[-1].recorded_at
                    values.append("R" if any(e.value == "R" for e in known if e.recorded_at == day) else "S")
            return dict(honest(recs, feats, cutoff), prior_r=values.count("R") / len(values))
        result = C.invariance_to_test_data(rows, REVIEW_FEATURES, T, prior_ties_as_r, antibiotics=("INH",))
        self.assertEqual(result.offenders, (("drop_non_training", "prior_r", "changed by test-period records"),
                                            ("training_view_only", "prior_r", "changed by test-period records")))

    def test_n2_controls_apply_one_at_a_time(self):
        # Flipping labels and complementing genes together cancels in make_fit's
        # centred linear score; each control alone moves it.
        recs, feats = records(), C._features(FEATURES)
        training = C.training_ids(recs, CUTOFF, ("INH",))
        fit = make_fit()
        base = C._fingerprints(fit(recs, feats, CUTOFF))
        for name, control in C.TRAINING_CONTROLS.items():
            with self.subTest(name):
                self.assertIn("model", C._changed(base, C._fingerprints(fit(*control(recs, feats, CUTOFF, training),
                                                                            CUTOFF))))

    def test_n2_no_training_rows_is_not_certified(self):
        result = C.invariance_to_test_data(RECORDS, FEATURES, "2025-01-01", honest, antibiotics=("INH",))
        self.assertEqual(result.offenders,
                         (("control", "*", "no training rows at the cutoff; invariance not certified"),))

    def test_n3_dedup_against_test_genomes_leaks(self):
        self.assertEqual(make_dedup("duplicate_group")(C.load_records(REVIEW_RECORDS), {}, CUTOFF)["train"],
                         [("t1", "S"), ("t3", "S"), ("t4", "R")])
        for field in ("duplicate_group", "lineage"):
            with self.subTest(field):
                result = self.run_review(make_dedup(field))
                self.assertFalse(result.passed)
                # new_genomes fires too: the added genomes copy training metadata.
                self.assertEqual({o[0] for o in result.offenders},
                                 {"metadata_shared", "metadata_fresh", "new_genomes", "drop_non_training",
                                  "training_view_only"})
                self.assertIn(("metadata_fresh", "train", "changed by test-period records"), result.offenders)
                self.assertIn(("metadata_shared", "train", "changed by test-period records"), result.offenders)
                # The reviewer's exact fit returns ids only. The drop-row control
                # moves it now (NN-1), so the metadata perturbations catch it.
                ids_only = self.run_review(make_dedup(field, False)).offenders
                self.assertIn(("metadata_shared", "train", "changed by test-period records"), ids_only)
                self.assertNotIn("control", {o[0] for o in ids_only})

    def test_n3_metadata_perturbations(self):
        recs = C.load_records(REVIEW_RECORDS)
        training = C.training_ids(recs, CUTOFF, ("INH",))
        self.assertEqual(training, {"t1", "t2", "t3", "t4"})
        shared = {r.genome_id: (r.lineage, r.duplicate_group, r.collection_date)
                  for r in C.perturb_metadata_shared(recs, {}, CUTOFF, training)[0]}
        self.assertEqual(shared["x1"], ("L2", "dgA", dt.date(2026, 1, 1)))
        self.assertEqual(shared["x2"], ("L4", "dgB", dt.date(2026, 1, 15)))
        fresh = {r.genome_id: (r.lineage, r.duplicate_group, r.collection_date)
                 for r in C.perturb_metadata_fresh(recs, {}, CUTOFF, training)[0]}
        self.assertEqual(fresh["x1"], ("canary-lineage-ec31682fde561917", "canary-group-ec31682fde561917", None))
        for gid in training:
            self.assertEqual(shared[gid], fresh[gid])
        added = [r for r in C.add_genomes(recs, {}, CUTOFF, training)[0] if r.genome_id.startswith("canary-new")]
        self.assertEqual([(r.lineage, r.duplicate_group) for r in added],
                         [("L2", "dgA"), ("L4", "dgB"), ("L1", "dgC"), ("L3", "dgD"), ("L2", "dgA"), ("L4", "dgB")])
        # The fifth is public at T with its only label after T (ruling A).
        self.assertTrue(C.knowable(added[4], CUTOFF))
        self.assertFalse(C.is_training_row(added[4], CUTOFF, ["INH"]))
        # The sixth is dated a day after every date in the records (NN-2): x1's INH label on 2026-05-01.
        self.assertEqual((added[5].collection_date, added[5].released_at, added[5].labels[0].recorded_at),
                         (dt.date(2026, 5, 2),) * 3)

    def test_perturbation_that_touches_training_data_is_a_canary_bug(self):
        def bad_features(recs, features, cutoff, training):
            return recs, dict(features, g02=frozenset({"inhA"}))

        def bad_membership(recs, features, cutoff, training):  # g03 gains a label at T
            return tuple(r if r.genome_id != "g03" else dataclasses.replace(
                r, labels=r.labels + (C.LabelEvent("INH", "S", cutoff, "x"),)) for r in recs), features

        def bad_cutoff_day(recs, features, cutoff, training):  # T-4: g04's label recorded on T flipped
            return tuple(r if r.genome_id != "g04" else dataclasses.replace(r, labels=tuple(
                dataclasses.replace(e, value=C._flip(e.value)) if e.recorded_at == cutoff else e for e in r.labels))
                for r in recs), features
        self.assertIn("g04", C.training_ids(RECORDS, T, ("INH",)))
        for bad in (bad_features, bad_membership, bad_cutoff_day):
            with self.subTest(bad.__name__), self.assertRaisesRegex(
                    ValueError, f"perturbation '{bad.__name__}' changed the training rows"):
                self.run_canary(make_fit(), perturbations={bad.__name__: bad})

    def test_builtin_perturbations_change_only_test_period_data(self):
        for recs, feats in ((records(), C._features(FEATURES)),
                            (C.load_records(REVIEW_RECORDS), C._features(REVIEW_FEATURES))):
            training = C.training_ids(recs, CUTOFF)
            view = C._training_view(recs, feats, CUTOFF, training)
            for name, perturb in C.PERTURBATIONS.items():
                new_recs, new_feats = perturb(recs, feats, CUTOFF, training)
                new_recs, new_feats = C.load_records(new_recs), C._features(new_feats)
                with self.subTest(name):
                    self.assertNotEqual((new_recs, new_feats), (recs, feats))
                    self.assertEqual(C.training_ids(new_recs, CUTOFF), training)
                    self.assertEqual(C._training_view(new_recs, new_feats, CUTOFF, training), view)

    def test_features_perturbation_covers_knowable_unlabelled_genomes(self):
        recs, feats = C.load_records(REVIEW_RECORDS), C._features(REVIEW_FEATURES)
        new = C.perturb_features(recs, feats, CUTOFF, C.training_ids(recs, CUTOFF))[1]
        self.assertEqual({gid for gid in feats if new[gid] != feats[gid]}, {"x1", "x2"})

    def test_label_perturbation_flips_later_events_and_adds_revisions(self):
        new = {r.genome_id: r for r in C.perturb_labels(records(), C._features(FEATURES), CUTOFF, frozenset())[0]}
        # NLK-1: g01's 2026-04-01 R from src-2, superseding 2025-08-01, becomes S
        # 400 days later from canary-resourced with no supersedes.
        self.assertEqual([(e.antibiotic, e.value, e.recorded_at.isoformat(), e.source_release,
                           e.supersedes and e.supersedes.isoformat()) for e in new["g01"].labels],
                         [("INH", "S", "2025-08-01", "src-1", None),
                          ("INH", "S", "2027-05-06", "canary-resourced", None),
                          ("INH", "R", "2027-05-07", "canary-later-release", "2027-05-06"),
                          ("RIF", "R", "2026-03-02", "canary-later-release", None)])
        # An undated event stays undated; with no supersedes it gains the latest dated event by T.
        undated = C.load_records([dict(RECORDS[1], labels=RECORDS[1]["labels"] + [event("INH", "R", None)])])
        rewritten = [e for e in C.perturb_labels(undated, {}, CUTOFF, frozenset())[0][0].labels if e.recorded_at is None]
        self.assertEqual([(e.value, e.source_release, e.supersedes) for e in rewritten],
                         [("S", "canary-resourced", dt.date(2025, 9, 1))])
        # The date cap keeps a slot after it: 9999-12-29 moves to 9999-12-30, which stays.
        # 9998-11-26 is the one date whose shift lands on 9999-12-31 without overflowing.
        # 9998-11-27 and 9999-01-01 overflow, and the fallback is the cap, not a day later (LEAK-2).
        self.assertEqual([C._later_date(dt.date(9998, 11, 26)), C._later_date(dt.date(9998, 11, 27)),
                          C._later_date(dt.date(9999, 1, 1)), C._later_date(dt.date(9999, 12, 29)),
                          C._later_date(dt.date(9999, 12, 30)), C._later_date(dt.date(9999, 12, 31))],
                         [dt.date(9999, 12, 30)] * 5 + [dt.date(9999, 12, 31)])
        self.assertEqual([(e.value, e.recorded_at.isoformat()) for e in new["g07"].labels if e.antibiotic == "RIF"],
                         [("R", "2025-04-01"), ("S", "2026-03-02")])

    def test_null_supersedes_points_at_latest_same_antibiotic_event_by_t(self):
        # T5-1: two dated INH events by T, the later on T itself, and a RIF event
        # by T that is earlier than INH's latest. max, <= cutoff and the
        # antibiotic each matter (NR5-13, NR5-15, NR5-14).
        rec = C.load_records([genome("s1", "2025-01-01", "2025-02-01", [
            event("INH", "S", "2025-06-01"), event("INH", "R", "2026-03-01", "2025-06-01"),
            event("INH", "S", "2026-05-01"), event("RIF", "S", "2025-07-01"), event("RIF", "R", "2026-04-01")])])
        labels = C.perturb_labels(rec, {}, CUTOFF, frozenset())[0][0].labels
        self.assertEqual(sorted((e.antibiotic, e.recorded_at.isoformat(), e.supersedes.isoformat())
                                for e in labels if e.source_release == C.LATER_SOURCE),
                         [("INH", "2027-06-05", "2026-03-01"), ("RIF", "2027-05-06", "2025-07-01")])

    def test_tie_and_undated_perturbations(self):
        def rows(perturb, gid):
            new = {r.genome_id: r for r in perturb(records(), C._features(FEATURES), CUTOFF, frozenset())[0]}
            return [(e.antibiotic, e.value, e.recorded_at and e.recorded_at.isoformat(), e.source_release)
                    for e in new[gid].labels]
        self.assertEqual(rows(C.perturb_label_ties, "g01"), [
            ("INH", "S", "2025-08-01", "src-1"),
            ("INH", "R", "2026-04-01", "src-2"), ("INH", "S", "2026-04-01", "canary-tie"),
            ("INH", "R", "2026-04-02", "canary-tie-r"), ("INH", "S", "2026-04-02", "canary-tie-s"),
            ("RIF", "R", "2026-03-02", "canary-tie-r"), ("RIF", "S", "2026-03-02", "canary-tie-s"),
        ])
        # A genome collected after the cutoff gets its new events no earlier than collection.
        late = C.load_records([genome("g10", "2026-06-01", "2026-06-02", [])])
        new = C.perturb_label_ties(late, {}, CUTOFF, frozenset())[0][0]
        self.assertEqual({e.recorded_at for e in new.labels}, {dt.date(2026, 6, 1)})
        new = C.perturb_labels(late, {}, CUTOFF, frozenset())[0][0]
        self.assertEqual({e.recorded_at for e in new.labels}, {dt.date(2026, 6, 1)})
        undated = rows(C.perturb_undated_labels, "g02")
        self.assertEqual([r for r in undated if r[2] is None], [
            ("INH", "R", None, "canary-undated"), ("INH", "S", None, "canary-undated"),
            ("RIF", "R", None, "canary-undated"), ("RIF", "S", None, "canary-undated")])

    def test_fingerprint(self):
        self.assertEqual(C.fingerprint({"a": {1, 2}, "b": np.arange(3.0)}),
                         C.fingerprint({"b": np.arange(3.0), "a": {2, 1}}))
        self.assertNotEqual(C.fingerprint(np.arange(3.0)), C.fingerprint(np.arange(3)))
        self.assertNotEqual(C.fingerprint(0.1 + 0.2), C.fingerprint(0.3))
        self.assertNotEqual(C.fingerprint([1]), C.fingerprint([True]))
        with self.assertRaises(TypeError):
            C.fingerprint(object())


# The joint re-check's fixture (grm-recheck-mnp-20261030/temporal/temporal_n_new_fits.py).
# Training rows at T: t1, t2, t3, t5, t6, t7. a1 is ambiguous at T; u1 and u2
# have no availability date (never knowable) but INH labels before T; f1's
# features came after T; p1 is public at T, labelled after it.
def mnp(gid, collected, released, labels, lineage, group, available=None):
    return {"genome_id": gid, "collection_date": collected, "released_at": released,
            "features_available_at": available, "lineage": lineage, "duplicate_group": group, "labels": labels}


MNP_RECORDS = [
    mnp("t1", "2024-12-01", "2025-01-10", [event("INH", "S", "2025-02-01", source="r1"),
                                           event("INH", "R", "2026-05-01", "2025-02-01", "r2")], "LA", "DA"),
    mnp("t2", "2025-01-01", "2025-02-10", [event("INH", "R", "2025-03-01", source="r1")], "LB", "DB"),
    mnp("t3", "2025-02-01", "2025-03-10", [event("INH", "S", "2025-04-01", source="r1")], "LC", "DC"),
    mnp("t5", "2025-03-01", "2025-04-10", [event("INH", "S", "2025-05-01", source="r1"),
                                           event("INH", "R", "2026-02-01", "2025-05-01", "r1")], "LD", "DD"),
    mnp("t6", "2025-04-01", "2025-05-10", [event("INH", "R", "2025-06-01", source="r1"),
                                           event("INH", "R", "2025-06-01", source="r0")], "LE", "DE"),
    mnp("t7", None, "2025-06-10", [event("INH", "S", "2025-07-01", source="r1")], None, "DF", "2025-06-20"),
    mnp("a1", "2025-05-01", "2025-06-01", [event("INH", "R", "2025-08-01", source="r1"),
                                           event("INH", "S", "2025-08-01", source="r0")], "LG", "DG"),
    mnp("u1", "2025-06-01", None, [event("INH", "R", "2025-09-01", source="r1")], "LH", "DH"),
    mnp("u2", None, None, [event("INH", "R", "2025-10-01", source="r1")], "LI", "DI"),
    mnp("f1", "2025-07-01", "2025-08-01", [event("INH", "R", "2025-09-01", source="r1")], "LJ", "DJ", "2026-04-01"),
    mnp("p1", "2025-08-01", "2025-09-01", [event("INH", "S", "2026-04-01", source="r1")], "LK", "DK"),
    mnp("x1", "2026-02-01", "2026-06-01", [event("INH", "R", "2026-06-15", source="r1")], "LB", "DX"),
    mnp("x2", "2026-04-01", "2026-08-01", [event("INH", "S", "2026-08-05", source="r1")], "LX", "DB"),
    mnp("x3", "2026-05-01", "2026-12-01", [], "LY", "DY"),
]
MNP_FEATURES = {"t1": {"katG", "ahpC"}, "t2": {"katG", "rpoB"}, "t3": {"inhA"}, "t5": {"katG", "inhA"},
                "t6": {"rpoB", "ahpC"}, "t7": {"inhA", "ahpC"}, "a1": {"katG", "embB"}, "u1": {"rpoB"},
                "u2": {"katG"}, "f1": {"embB"}, "p1": {"geneQ"}, "x1": {"embB", "katG"}, "x2": {"rpoB_S450L"},
                "x3": {"geneZ"}}
# u1, u2 and f1 relabelled S: the R share over every genome labelled at T,
# availability ignored, is exactly 1/3, and stays 1/3 after new_genomes (TL-3).
MNP_THIRD = [dict(r, labels=[dict(e, value="S") for e in r["labels"]]) if r["genome_id"] in ("u1", "u2", "f1")
             else r for r in MNP_RECORDS]
MNP_REGISTRY = [
    {"kind": "catalog", "release_id": "c1", "released_at": "2025-11-01", "sha256": sha("c1"),
     "vocabulary": ["katG", "rpoB", "inhA", "ahpC"]},
    {"kind": "catalog", "release_id": "c2", "released_at": "2026-05-01", "sha256": sha("c2"),
     "vocabulary": ["katG", "rpoB", "inhA", "ahpC", "embB"]},
    {"kind": "annotation", "release_id": "a1", "released_at": "2025-10-01", "sha256": sha("a1"),
     "vocabulary": ["geneQ"]},
]


def mnp_train(recs, cutoff):
    return [r for r in recs if C.is_training_row(r, cutoff, ["INH"])]


def mnp_honest(recs, feats, cutoff):
    rows = mnp_train(recs, cutoff)
    labels = {r.genome_id: C.label_as_of(r, "INH", cutoff).value for r in rows}
    vocab = sorted(set().union(*(feats.get(r.genome_id, ()) for r in rows)))
    return {"vocabulary": vocab, "weights": {gene: sum(1 if labels[r.genome_id] == "R" else -1 for r in rows
                                                       if gene in feats.get(r.genome_id, ())) for gene in vocab}}


def share_r(values):
    values = list(values)
    return sum(v == "R" for v in values) / len(values) if values else None


def nl5_prior_none_is_available(recs, feats, cutoff):
    """TL-1: a falsy availability test, so a genome with no availability date counts."""
    values = []
    for r in recs:
        if r.available_at and r.available_at > cutoff:
            continue
        try:
            found = C.label_as_of(r, "INH", cutoff)
        except C.AmbiguousLabel:
            continue
        if found is not None:
            values.append(found.value)
    return dict(mnp_honest(recs, feats, cutoff), prior_r=share_r(values),
                class_weight_r=values.count("S") / max(1, values.count("R")))


def nl6_prior_ties_resolve_to_r(recs, feats, cutoff):
    """TL-2: a prior over knowable genomes; a same-day R/S conflict counts as R."""
    values = []
    for r in recs:
        known = C.known_events(r, "INH", cutoff) if C.knowable(r, cutoff) else []
        if known:
            day = known[-1].recorded_at
            values.append("R" if any(e.value == "R" for e in known if e.recorded_at == day) else "S")
    return dict(mnp_honest(recs, feats, cutoff), prior_r=share_r(values))


def nl7_prior_ignores_availability(recs, feats, cutoff):
    """TL-3: a prior over every genome labelled at T, availability ignored."""
    values = []
    for r in recs:
        try:
            found = C.label_as_of(r, "INH", cutoff)
        except C.AmbiguousLabel:
            continue
        if found is not None:
            values.append(found.value)
    return dict(mnp_honest(recs, feats, cutoff), prior_r=share_r(values))


def nc1_age_vs_cutoff(recs, feats, cutoff):
    return dict(mnp_honest(recs, feats, cutoff), age={
        r.genome_id: [None if d is None else (cutoff - d).days
                      for d in (r.available_at, r.collection_date, C.label_as_of(r, "INH", cutoff).recorded_at)]
        for r in mnp_train(recs, cutoff)})


def nc2_registry_vocab_declared(recs, feats, cutoff):
    vocab = A.ArtifactRegistry.from_records(MNP_REGISTRY).vocabulary(cutoff)
    return dict(mnp_honest(recs, feats, cutoff), registry_vocabulary=sorted(vocab),
                encoded={r.genome_id: A.encode_genes(feats.get(r.genome_id, ()), vocab)
                         for r in mnp_train(recs, cutoff)})


def nc3_n_train(recs, feats, cutoff):
    rows = mnp_train(recs, cutoff)
    return {"n_train": len(rows), "ids": sorted(r.genome_id for r in rows)}


def nc4_prior_training_only(recs, feats, cutoff):
    return dict(mnp_honest(recs, feats, cutoff),
                prior_r=share_r(C.label_as_of(r, "INH", cutoff).value for r in mnp_train(recs, cutoff)))


def nc5_m_train_embargo_30(recs, feats, cutoff):
    manifest = M.split_at([r.to_mapping() for r in recs], cutoff, 90, 30, cutoff + dt.timedelta(days=400), "INH")
    return {"train": manifest["train"],
            "vocab": sorted(set().union(*(feats.get(gid, ()) for gid in manifest["train"])))}


def nc6_hyperparameters_declared(recs, feats, cutoff):
    return dict(mnp_honest(recs, feats, cutoff), schema=3, l2=0.5)


class DropNonTraining(unittest.TestCase):
    """Round 4: TL-1 to TL-5, AGR-2 and T-3 on the joint re-check's fixture."""

    def run_mnp(self, fit, rows=MNP_RECORDS, **kwargs):
        return C.invariance_to_test_data(rows, MNP_FEATURES, T, fit, antibiotics=["INH"], **kwargs)

    def test_fixture(self):
        recs = C.load_records(MNP_RECORDS)
        self.assertEqual(C.training_ids(recs, T, ["INH"]), {"t1", "t2", "t3", "t5", "t6", "t7"})
        by = {r.genome_id: r for r in recs}
        self.assertEqual([C.knowable(by[g], CUTOFF) for g in ("a1", "u1", "u2", "f1", "p1")],
                         [True, False, False, False, True])

    def test_label_only_leaks_fail(self):
        # TL-1 (NL5), TL-2 (NL6) and TL-3 (NL7b): drop_non_training and training_view_only
        # alone catch each.
        for name, fit, rows, keys in (
                ("NL5", nl5_prior_none_is_available, MNP_RECORDS, {"prior_r", "class_weight_r"}),
                ("NL6", nl6_prior_ties_resolve_to_r, MNP_RECORDS, {"prior_r"}),
                ("NL7b", nl7_prior_ignores_availability, MNP_THIRD, {"prior_r"})):
            with self.subTest(name):
                result = self.run_mnp(fit, rows)
                self.assertFalse(result.passed)
                self.assertEqual({(o[0], o[1]) for o in result.offenders},
                                 {(p, k) for p in ("drop_non_training", "training_view_only") for k in keys})
                without = dict(C.PERTURBATIONS)
                del without["drop_non_training"], without["training_view_only"]
                self.assertTrue(self.run_mnp(fit, rows, perturbations=without).passed)  # the round-3 gap

    def test_tl3_not_yet_public_labels(self):
        # NL7 on the generic fixture is also caught by new_genomes; on MNP_THIRD
        # the added labels leave the 1/3 share unchanged, so only the drop catches it.
        offenders = self.run_mnp(nl7_prior_ignores_availability).offenders
        self.assertEqual({o[0] for o in offenders}, {"new_genomes", "drop_non_training", "training_view_only"})
        recs = C.load_records(MNP_THIRD)
        training = C.training_ids(recs, CUTOFF, ["INH"])
        feats = C._features(MNP_FEATURES)
        full = nl7_prior_ignores_availability(recs, feats, CUTOFF)["prior_r"]
        added = nl7_prior_ignores_availability(*C.add_genomes(recs, feats, CUTOFF, training), CUTOFF)["prior_r"]
        alone = nl7_prior_ignores_availability(*C.drop_non_training(recs, feats, CUTOFF, training), CUTOFF)["prior_r"]
        self.assertEqual((full, added, alone), (1 / 3, 1 / 3, 0.5))

    def test_clean_fits_pass(self):
        for name, fit, kwargs in (("NC1", nc1_age_vs_cutoff, {}),
                                  ("NC2", nc2_registry_vocab_declared, {"constant_keys": ["registry_vocabulary"]}),
                                  ("NC3", nc3_n_train, {}), ("NC4", nc4_prior_training_only, {}),
                                  ("NC5", nc5_m_train_embargo_30, {}),
                                  ("NC6", nc6_hyperparameters_declared, {"constant_keys": ["schema", "l2"]})):
            with self.subTest(name):
                if fit is nc5_m_train_embargo_30 and M is None:
                    self.skipTest("Task M's splits.py is neither merged nor in its worktree")
                result = self.run_mnp(fit, **kwargs)
                self.assertTrue(result.passed, result)

    def test_drop_non_training_keeps_exactly_the_training_rows(self):
        recs, feats = C.load_records(MNP_RECORDS), C._features(MNP_FEATURES)
        training = C.training_ids(recs, CUTOFF, ["INH"])
        kept, kept_features = C.drop_non_training(recs, dict(feats, orphan=frozenset({"katG"})), CUTOFF, training)
        self.assertEqual(tuple(kept), tuple(r for r in recs if r.genome_id in training))
        self.assertEqual(kept_features, {gid: feats[gid] for gid in training})

    def test_tl4_declared_constant_keys_are_reported(self):
        passed = self.run_mnp(nc6_hyperparameters_declared, constant_keys=["schema", "l2"])
        self.assertEqual(passed.constant_keys, ("l2", "schema"))
        self.assertEqual(passed.to_dict()["constant_keys"], ["l2", "schema"])
        self.assertIn("constant_keys declared by the caller, each a claim to check: ['l2', 'schema']", passed.detail)
        failed = self.run_mnp(nc6_hyperparameters_declared, constant_keys=["schema"])
        self.assertFalse(failed.passed)
        self.assertEqual(failed.to_dict()["constant_keys"], ["schema"])
        self.assertIn("['schema']", failed.detail)
        plain = self.run_mnp(nc3_n_train)
        self.assertEqual((plain.constant_keys, plain.to_dict()["constant_keys"]), ((), []))
        self.assertNotIn("constant_keys", plain.detail)

    def test_agr2_fresh_metadata_always_passes_rule6(self):
        # The re-check's ids: case-only twins, non-ASCII and control characters,
        # and an id whose old fresh value aliased a training lineage.
        def rows(*ids, training_lineage="LT"):
            out = [genome("t1", "2025-01-01", "2025-02-01", [event("INH", "S", "2025-03-01")], lineage=training_lineage),
                   genome("t2", "2025-01-01", "2025-02-01", [event("INH", "R", "2025-03-02")], lineage="LU")]
            return out + [genome(gid, "2026-04-01", "2026-05-01", [], lineage=None) for gid in ids]
        for name, recs in (("case", rows("x1", "X1")), ("e-acute", rows("gé1")), ("tab", rows("g\t1")),
                           ("zero-width", rows("g​1")), ("training alias", rows("y1", training_lineage="canary-lineage-Y1"))):
            with self.subTest(name):
                loaded = C.load_records(recs)
                fresh = C.load_records(C.perturb_metadata_fresh(loaded, {}, CUTOFF, C.training_ids(loaded, T))[0])
                for r in fresh:
                    if r.genome_id not in ("t1", "t2"):
                        self.assertRegex(r.lineage, r"^canary-lineage-[0-9a-f]{16}$")
                        self.assertRegex(r.duplicate_group, r"^canary-group-[0-9a-f]{16}$")
                self.assertTrue(C.invariance_to_test_data(recs, {}, T, honest, antibiotics=("INH",)).passed)
        self.assertEqual(C._id_digest("x1"), hashlib.sha256(b"x1").hexdigest()[:16])
        self.assertNotEqual(C._id_digest("x1"), C._id_digest("X1"))

    def test_t3_dates_on_the_cutoff_keep_their_side(self):
        recs = C.load_records([genome("t1", "2025-01-01", "2025-02-01", [event("INH", "S", "2025-03-01")]),
                               genome("t2", "2025-01-01", "2025-02-01", [event("INH", "R", "2025-03-02")]),
                               genome("u0", "2025-06-01", T, [], available=T),
                               genome("u9", "2025-06-01", "2026-03-02", [])])
        moved = {r.genome_id: r for r in C.perturb_dates(recs, {}, CUTOFF, C.training_ids(recs, T))[0]}
        self.assertEqual((moved["u0"].released_at, moved["u0"].features_available_at), (CUTOFF - C.DATE_SHIFT,) * 2)
        self.assertTrue(C.knowable(moved["u0"], CUTOFF))
        self.assertEqual(moved["u9"].released_at, dt.date(2026, 3, 2) + C.DATE_SHIFT)
        self.assertEqual(C._shift(CUTOFF, CUTOFF), CUTOFF - C.DATE_SHIFT)

    def test_nl9_fit_that_raises_under_a_perturbation_fails(self):
        def nl9_days_to_next_release(recs, feats, cutoff):  # the re-check's NL9 leak
            nxt = min(r.available_at for r in recs if r.available_at and r.available_at > cutoff)
            return dict(mnp_honest(recs, feats, cutoff),
                        to_next={r.genome_id: (nxt - r.available_at).days for r in mnp_train(recs, cutoff)})
        result = self.run_mnp(nl9_days_to_next_release)
        self.assertFalse(result.passed)
        self.assertIn(("drop_non_training", "*", "fit raised ValueError: min() iterable argument is empty"),
                      result.offenders)
        self.assertIn(("dates", "to_next", "changed by test-period records"), result.offenders)

        def long_message(recs, feats, cutoff):
            if "x3" not in {r.genome_id for r in recs}:  # only drop_non_training removes x3
                raise KeyError("k" * 500)
            return mnp_honest(recs, feats, cutoff)
        reason = [o[2] for o in self.run_mnp(long_message).offenders if o[0] == "drop_non_training"]
        self.assertEqual(reason, ["fit raised KeyError: '" + "k" * 199])

        def interrupted(recs, feats, cutoff):
            if "x3" not in {r.genome_id for r in recs}:  # only drop_non_training removes x3
                raise KeyboardInterrupt
            return mnp_honest(recs, feats, cutoff)
        with self.assertRaises(KeyboardInterrupt):
            self.run_mnp(interrupted, perturbations={"drop_non_training": C.drop_non_training})

        def broken(recs, feats, cutoff):  # the unchanged inputs and the controls stay loud
            raise RuntimeError("setup")
        with self.assertRaisesRegex(RuntimeError, "setup"):
            self.run_mnp(broken)
        calls = []

        def breaks_in_control(recs, feats, cutoff):
            calls.append(1)
            if len(calls) > 2:
                raise RuntimeError("control")
            return mnp_honest(recs, feats, cutoff)
        with self.assertRaisesRegex(RuntimeError, "control"):
            self.run_mnp(breaks_in_control)

    def test_tl5_label_on_the_last_date_is_refused(self):
        rows = REVIEW_RECORDS + [genome("q9", "2025-01-01", "2025-02-01", [event("INH", "S", "9999-12-31")])]
        for name in ("labels", "label_ties"):
            with self.subTest(name), self.assertRaisesRegex(ValueError, "q9: an event or cutoff on 9999-12-31"):
                C.invariance_to_test_data(rows, REVIEW_FEATURES, T, honest, antibiotics=("INH",),
                                          perturbations={name: C.PERTURBATIONS[name]})


# The round-5 re-check's fixture (grm-recheck-np-r2-20261030/nleak/nleak_r5_fits.py).
# Every training row in NLK_FULL has a dated INH event after T; in NLK_PARTIAL t6 has none.
class FitError(Exception):
    """A fit's own exception type, outside every built-in family (T5-2)."""


def nlk(gid, collected, released, labels, lineage, group, available=None):
    return {"genome_id": gid, "collection_date": collected, "released_at": released,
            "features_available_at": available, "lineage": lineage, "duplicate_group": group, "labels": labels}


def nev(value, recorded, source, supersedes=None, antibiotic="INH"):
    return event(antibiotic, value, recorded, supersedes, source)


NLK_TRAIN = [
    nlk("t1", "2025-05-01", "2025-06-01", [nev("R", "2025-08-01", "labA"), nev("S", "2026-04-15", "labA", "2025-08-01")],
        "L1", "D1"),
    nlk("t2", "2025-06-01", "2025-07-01", [nev("S", "2025-09-01", "labB"), nev("S", "2026-06-01", "labC")], "L2", "D2"),
    nlk("t3", "2025-04-01", "2025-05-01", [nev("S", "2025-06-01", "labB"), nev("R", "2025-10-01", "labA", "2025-06-01"),
                                           nev("R", "2026-03-20", "labA")], "L1", "D3"),
    nlk("t4", "2025-08-01", "2025-09-01", [nev("S", "2025-11-01", "labB"), nev("R", "2026-09-01", "labD", "2025-11-01")],
        "L3", "D4"),
    nlk("t5", "2025-03-01", "2025-04-01", [nev("R", "2025-12-01", "labC"), nev("S", "2026-05-05", "labC")], "L2", "D5"),
    nlk("t6", "2025-09-01", "2025-10-01", [nev("S", "2026-01-10", "labD"), nev("S", "2026-03-02", "labB")], "L3", "D6"),
]
NLK_NON_TRAIN = [
    nlk("n1", "2025-07-01", "2025-08-01", [nev("S", "2026-04-01", "labA")], "L1", "D7"),
    nlk("n2", "2025-12-01", "2026-04-01", [nev("R", "2026-02-01", "labB")], "L4", "D8"),
    nlk("n3", "2025-09-01", None, [nev("S", "2025-10-01", "labC")], "L2", "D9"),
    nlk("n4", "2025-10-01", "2025-11-01", [nev("R", "2026-02-15", "labA"), nev("S", "2026-02-15", "labB")], "L3", "D10"),
    nlk("n5", "2025-10-01", "2025-11-01", [nev("R", None, "labD")], "L4", "D11"),
    nlk("n6", "2025-10-01", "2025-12-01", [nev("R", "2026-01-01", "labA", antibiotic="RIF")], "L1", "D12"),
]
NLK_FULL = NLK_TRAIN + NLK_NON_TRAIN
NLK_PARTIAL = NLK_TRAIN[:5] + [nlk("t6", "2025-09-01", "2025-10-01", [nev("S", "2026-01-10", "labD")], "L3", "D6")] \
    + NLK_NON_TRAIN
NLK_FEATURES = {"t1": ["katG", "inhA"], "t2": ["rpoB"], "t3": ["katG", "ahpC", "fabG1"], "t4": ["rpoB", "embB"],
                "t5": ["katG", "inhA", "gyrA"], "t6": ["pncA"], "n1": ["katG", "rrs"], "n2": ["eis", "inhA"],
                "n3": ["katG"], "n4": ["rpoB", "rrs"], "n5": ["gyrA", "eis"], "n6": ["embB"], "o1": ["katG", "zz1"]}
NLK_REGISTRY = [
    {"kind": "catalog", "release_id": "c1", "released_at": "2025-11-01", "sha256": sha("c1"),
     "vocabulary": ["katG", "rpoB", "inhA", "ahpC", "fabG1"]},
    {"kind": "catalog", "release_id": "c2", "released_at": "2026-05-01", "sha256": sha("c2"),
     "vocabulary": ["katG", "rpoB", "inhA", "ahpC", "fabG1", "embB"]},
    {"kind": "annotation", "release_id": "a1", "released_at": "2025-10-01", "sha256": sha("a1"), "vocabulary": ["gyrA"]},
]


def nlk_rows(recs, cutoff):
    return [r for r in recs if C.is_training_row(r, cutoff, ["INH"])]


def nlk_later(r, cutoff):
    return sorted((e for e in r.labels if e.antibiotic == "INH" and e.recorded_at is not None and e.recorded_at > cutoff),
                  key=lambda e: e.recorded_at)


def nlk_model(rows, feats, cutoff, weight=None):
    """The re-check's honest weighted model: responds to labels, genes and rows."""
    weight = weight or {}
    w = {r.genome_id: float(weight.get(r.genome_id, 1.0)) for r in rows}
    y = {r.genome_id: 1.0 if C.label_as_of(r, "INH", cutoff).value == "R" else 0.0 for r in rows}
    vocab = sorted(set().union(*(feats.get(r.genome_id, ()) for r in rows))) if rows else []
    return {"prior_r": sum(w[k] * y[k] for k in w) / sum(w.values()), "vocabulary": vocab, "n_train": len(rows),
            "gene_weight": {g: sum(w[r.genome_id] * y[r.genome_id] for r in rows if g in feats.get(r.genome_id, ()))
                            for g in vocab}}


def make_lk7(weight_of):
    def fit(recs, feats, cutoff):
        rows = nlk_rows(recs, cutoff)
        return nlk_model(rows, feats, cutoff, {r.genome_id: weight_of(r, cutoff) for r in rows})
    return fit


def w_first_later(r, cutoff):  # LK7a: the date of the first report after T
    events = nlk_later(r, cutoff)
    return 1.0 / (1 + (events[0].recorded_at - cutoff).days) if events else 1.0


def w_same_source(r, cutoff):  # LK7b: the as-of-T label's own source reported again after T
    known = C.label_as_of(r, "INH", cutoff)
    return 2.0 if any(e.source_release == known.source_release for e in nlk_later(r, cutoff)) else 1.0


def w_superseded(r, cutoff):  # LK7d: the as-of-T label superseded after T
    known = C.label_as_of(r, "INH", cutoff)
    return 3.0 if any(e.supersedes == known.recorded_at for e in nlk_later(r, cutoff)) else 1.0


def w_undated_same_source(r, cutoff):  # the verifier's undated variant
    known = C.label_as_of(r, "INH", cutoff)
    return 2.0 if any(e.recorded_at is None and e.source_release == known.source_release for e in r.labels) else 1.0


def cf1_training_scaler(recs, feats, cutoff):
    rows = nlk_rows(recs, cutoff)
    counts = np.array([len(feats.get(r.genome_id, ())) for r in rows], dtype=float)
    mu, sd = counts.mean(), counts.std()
    return dict(nlk_model(rows, feats, cutoff), z=(counts - mu) / (sd if sd > 0 else 1.0), mu_sd=[float(mu), float(sd)])


def cf2_registry_declared(recs, feats, cutoff):
    rows = nlk_rows(recs, cutoff)
    vocab = sorted(A.ArtifactRegistry.from_records(NLK_REGISTRY).vocabulary(cutoff))
    return {"registry_vocabulary": vocab,
            "encoded": {r.genome_id: A.encode_genes(feats.get(r.genome_id, ()), set(vocab)) for r in rows},
            "labels": {r.genome_id: C.label_as_of(r, "INH", cutoff).value for r in rows}}


def cf3_n_train(recs, feats, cutoff):
    rows = nlk_rows(recs, cutoff)
    return {"n_train": len(rows), "ids": [r.genome_id for r in rows]}


def cf4_revisions_by_t(recs, feats, cutoff):
    rows = nlk_rows(recs, cutoff)
    return dict(nlk_model(rows, feats, cutoff), revisions_by_t={
        r.genome_id: sum(1 for e in r.labels if e.antibiotic == "INH" and e.recorded_at and e.recorded_at <= cutoff
                         and e.supersedes) for r in rows})


def cf5_bootstrap_constant_seed(recs, feats, cutoff):
    rows = nlk_rows(recs, cutoff)
    y = np.array([1.0 if C.label_as_of(r, "INH", cutoff).value == "R" else 0.0 for r in rows])
    draws = np.random.default_rng(0).integers(0, len(rows), size=(16, len(rows)))
    return {"boot_prior": y[draws].mean(axis=1), "n_train": len(rows)}


def cf6_seed_n_train(recs, feats, cutoff):
    out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
    rng = np.random.default_rng(out["n_train"])
    out["init"] = rng.normal(size=len(out["vocabulary"])) + np.array([out["gene_weight"][g] for g in out["vocabulary"]])
    return out


def flip_only_labels(recs, features, cutoff, training):
    """Round 4's perturb_labels: a later event's value flipped, its date, source and supersedes kept."""
    with mock.patch.object(C, "_rewrite_later", lambda e, record, c: dataclasses.replace(e, value=C._flip(e.value))):
        return C.perturb_labels(recs, features, cutoff, training)


class PostCutoffEvents(unittest.TestCase):
    """Round 5: NLK-1, NLK-2 and TM-N1 to TM-N5."""

    def run_nlk(self, fit, rows=NLK_FULL, **kwargs):
        return C.invariance_to_test_data(rows, NLK_FEATURES, T, fit, antibiotics=["INH"], **kwargs)

    def test_fixture(self):
        self.assertEqual(C.training_ids(NLK_FULL, T, ["INH"]), {"t1", "t2", "t3", "t4", "t5", "t6"})
        self.assertEqual(C.training_ids(NLK_PARTIAL, T, ["INH"]), {"t1", "t2", "t3", "t4", "t5", "t6"})

    def test_nlk1_later_report_fits_fail(self):
        for name, weight, rows in (("LK7a first later date", w_first_later, NLK_FULL),
                                   ("LK7b same source again", w_same_source, NLK_FULL),
                                   ("LK7b' same source on PARTIAL", w_same_source, NLK_PARTIAL),
                                   ("LK7d superseded after T", w_superseded, NLK_FULL)):
            with self.subTest(name):
                result = self.run_nlk(make_lk7(weight), rows)
                self.assertFalse(result.passed)
                kinds = {o[0] for o in result.offenders}
                self.assertTrue({"drop_post_cutoff_events", "labels"} <= kinds, kinds)
                # Round 4's defaults (labels flipped values only, no drop, no training view) passed each one.
                old = dict(C.PERTURBATIONS, labels=flip_only_labels)
                del old["drop_post_cutoff_events"], old["training_view_only"]
                self.assertTrue(self.run_nlk(make_lk7(weight), rows, perturbations=old).passed)
        undated = NLK_FULL + [nlk("t7", "2025-02-01", "2025-03-01", [nev("R", "2025-04-01", "labZ"),
                                                                       nev("R", None, "labZ")], "L1", "D13")]
        offenders = self.run_nlk(make_lk7(w_undated_same_source), undated).offenders
        self.assertIn(("drop_post_cutoff_events", "prior_r", "changed by test-period records"), offenders)

    def test_nlk1_clean_fits_pass(self):
        for name, fit, rows, kwargs in (("CF1", cf1_training_scaler, NLK_FULL, {}),
                                        ("CF1'", cf1_training_scaler, NLK_PARTIAL, {}),
                                        ("CF2", cf2_registry_declared, NLK_FULL,
                                         {"constant_keys": ["registry_vocabulary"]}),
                                        ("CF3", cf3_n_train, NLK_FULL, {}), ("CF4", cf4_revisions_by_t, NLK_FULL, {}),
                                        ("CF4'", cf4_revisions_by_t, NLK_PARTIAL, {}),
                                        ("CF5", cf5_bootstrap_constant_seed, NLK_FULL, {}),
                                        ("CF6", cf6_seed_n_train, NLK_FULL, {})):
            with self.subTest(name):
                result = self.run_nlk(fit, rows, **kwargs)
                self.assertTrue(result.passed, result)

    def test_drop_post_cutoff_events_keeps_only_dated_events_by_t(self):
        recs = C.load_records(NLK_FULL)
        dropped = {r.genome_id: r for r in C.drop_post_cutoff_events(recs, {}, CUTOFF, frozenset())[0]}
        for r in recs:
            self.assertEqual(dropped[r.genome_id].labels,
                             tuple(e for e in r.labels if e.recorded_at is not None and e.recorded_at <= CUTOFF))
        self.assertEqual(dropped["t1"].labels, (C.LabelEvent("INH", "R", dt.date(2025, 8, 1), "labA"),))
        self.assertEqual(dropped["n5"].labels, ())

    def test_nlk2_keys_that_print_alike_are_refused(self):
        def collide(first, second):
            def fit(recs, feats, cutoff):
                out = {first: sorted(set().union(*feats.values()))}  # vocabulary over every genome
                out.update(nlk_model(nlk_rows(recs, cutoff), feats, cutoff))
                out[second] = out["vocabulary"]
                return out
            return fit
        # Since MC3-1 every fitted key must be an exact str, so the non-str key
        # of each pair is refused before any two keys can collide.
        for first, second in ((0, "0"), ("0", 0), (None, "None"), (dt.date(2026, 3, 1), "2026-03-01"), (True, "True")):
            odd = first if type(first) is not str else second
            with self.subTest(first=first, second=second), self.assertRaises(TypeError) as caught:
                self.run_nlk(collide(first, second))
            self.assertIn(f"fitted key {odd!r} is a", str(caught.exception))
        with self.assertRaisesRegex(TypeError, "must each be an exact str"):
            C._fingerprints({0: 1, "a": 2})
        self.assertEqual(C._fingerprints({"0": 1, "a": 2}).keys(), {"0", "a"})

    def test_tm_n1_any_exception_from_the_perturbed_fit_is_an_offender(self):
        def raising(error):
            def fit(recs, feats, cutoff):
                if "n1" not in {r.genome_id for r in recs}:  # only drop_non_training removes n1
                    raise error
                return nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
            return fit
        for error, text in ((ZeroDivisionError("division by zero"), "ZeroDivisionError: division by zero"),
                            (TypeError("bad operand"), "TypeError: bad operand"),
                            (AttributeError("no lower"), "AttributeError: no lower"),
                            (RuntimeError("boom"), "RuntimeError: boom"),
                            (AssertionError("assert"), "AssertionError: assert"),
                            (OSError("disk"), "OSError: disk"), (StopIteration("done"), "StopIteration: done"),
                            (FitError("custom"), "FitError: custom")):
            with self.subTest(type(error).__name__):
                self.assertEqual(self.run_nlk(raising(error)).offenders,
                                 (("drop_non_training", "*", f"fit raised {text}"),
                                  ("training_view_only", "*", f"fit raised {text}")))

    def test_tm_n2_the_loop_continues_after_an_offender(self):
        def raises_under_dates(recs, feats, cutoff):  # dates is first in sorted order
            if any(r.genome_id == "n1" and r.released_at != dt.date(2025, 8, 1) for r in recs):
                raise TypeError("n1 moved")
            return nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
        self.assertEqual(sorted(C.PERTURBATIONS)[0], "dates")
        self.assertEqual(self.run_nlk(raises_under_dates).offenders, (("dates", "*", "fit raised TypeError: n1 moved"),))

        def raises_under_dates_reads_later(recs, feats, cutoff):  # T5-3: later perturbations still run
            later = tuple(e.value for r in recs for e in nlk_later(r, cutoff))
            return dict(raises_under_dates(recs, feats, cutoff), later_values=later)
        self.assertEqual(self.run_nlk(raises_under_dates_reads_later).offenders,
                         (("dates", "*", "fit raised TypeError: n1 moved"),
                          ("drop_non_training", "later_values", "changed by test-period records"),
                          ("drop_post_cutoff_events", "later_values", "changed by test-period records"),
                          ("label_ties", "later_values", "changed by test-period records"),
                          ("labels", "later_values", "changed by test-period records"),
                          ("new_genomes", "later_values", "changed by test-period records"),
                          ("training_view_only", "later_values", "changed by test-period records")))

        def counts_n1_raises_under_features(recs, feats, cutoff):
            if feats.get("n1") != frozenset({"katG", "rrs"}) and "n1" in feats:
                raise TypeError("n1 genes")
            rows = nlk_rows(recs, cutoff)
            return dict(nlk_model(rows, feats, cutoff), rows_and_n1=len(rows) + sum(r.genome_id == "n1" for r in recs))
        self.assertEqual(self.run_nlk(counts_n1_raises_under_features).offenders,
                         (("drop_non_training", "rows_and_n1", "changed by test-period records"),
                          ("features", "*", "fit raised TypeError: n1 genes"),
                          ("training_view_only", "rows_and_n1", "changed by test-period records")))

    def test_tm_n3_invalid_records_from_a_perturbation_propagate(self):
        def bad_rule6(recs, features, cutoff, training):  # non-training n1 gets an alias of t1's lineage
            return tuple(dataclasses.replace(r, lineage="l1") if r.genome_id == "n1" else r for r in recs), features
        with self.assertRaisesRegex(ValueError, "lineage alias collision: 'l1' \\(n1\\) and 'L1'"):
            self.run_nlk(make_lk7(lambda r, c: 1.0), perturbations={"bad_rule6": bad_rule6})

        def bad_features(recs, features, cutoff, training):
            return recs, dict(features, n1="katG")
        with self.assertRaisesRegex(ValueError, "n1: features must be a collection of strings"):
            self.run_nlk(make_lk7(lambda r, c: 1.0), perturbations={"bad_features": bad_features})

    def test_tm_n4_fresh_digest_encoding(self):
        self.assertEqual((C._id_digest("gé1"), C._id_digest("gè1")), ("c386226d0a03c5fa", "f4855fd1e9fc7225"))
        self.assertEqual(C._id_digest("g\ud800"), hashlib.sha256(b"g\xed\xa0\x80").hexdigest()[:16])
        rows = [genome("t1", "2025-01-01", "2025-02-01", [event("INH", "S", "2025-03-01")]),
                genome("t2", "2025-01-01", "2025-02-01", [event("INH", "R", "2025-03-02")]),
                genome("gé1", "2026-04-01", "2026-05-01", []), genome("gè1", "2026-04-01", "2026-05-01", []),
                genome("g\ud800", "2026-04-01", "2026-05-01", [])]
        self.assertTrue(C.invariance_to_test_data(rows, {}, T, honest, antibiotics=("INH",)).passed)

    def test_tm_n5_last_date_boundary(self):
        def with_label(day):
            return REVIEW_RECORDS + [genome("q9", "2025-01-01", "2025-02-01", [event("INH", "S", day)])]
        self.assertTrue(C.invariance_to_test_data(with_label("9999-12-30"), REVIEW_FEATURES, T, honest,
                                                  antibiotics=("INH",)).passed)
        with self.assertRaisesRegex(ValueError, "q9: an event or cutoff on 9999-12-31"):
            C.invariance_to_test_data(with_label("9999-12-31"), REVIEW_FEATURES, T, honest, antibiotics=("INH",))


def latest_date(recs):
    dates = [d for r in recs for d in (r.collection_date, r.released_at, r.features_available_at) if d is not None]
    return max(dates + [e.recorded_at for r in recs for e in r.labels if e.recorded_at is not None])


def vr1_horizon(recs, feats, cutoff):
    # LEAK-1 (the verifier's VR1): label ages from min(T, latest date anywhere in the input).
    rows = nlk_rows(recs, cutoff)
    ref = min(cutoff, latest_date(recs))
    return dict(nlk_model(rows, feats, cutoff),
                age={r.genome_id: (ref - C.label_as_of(r, "INH", cutoff).recorded_at).days for r in rows})


def lr2_embargo(recs, feats, cutoff):
    # LEAK-1 (the leak lens's LR2): a 60-day embargo only when the input has a date after T.
    embargo = dt.timedelta(days=60) if latest_date(recs) > cutoff else dt.timedelta(0)
    rows = [r for r in nlk_rows(recs, cutoff) if C.label_as_of(r, "INH", cutoff).recorded_at <= cutoff - embargo]
    out = nlk_model(rows, feats, cutoff) if rows else {"prior_r": None}
    out.pop("n_train", None)  # the training_rows control drops t6, which the embargo already left out
    return out


def lr_min_prior(recs, feats, cutoff):
    # LEAK-1 (fin_lr_min.py): the prior is shrunk when any date in the input is after T.
    out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
    if latest_date(recs) > cutoff:
        out["prior_r"] = out["prior_r"] * 0.9 + 0.01  # no fixed point at the fixture's prior of 0.5
    return out


def clean_horizon(recs, feats, cutoff):
    # Clean: ages from the latest date in the fit's own training view.
    view = C._training_view(recs, feats, cutoff, {r.genome_id for r in nlk_rows(recs, cutoff)})
    ref = latest_date([record for record, _ in view.values()])
    rows = nlk_rows(recs, cutoff)
    return dict(nlk_model(rows, feats, cutoff),
                age={r.genome_id: (ref - C.label_as_of(r, "INH", cutoff).recorded_at).days for r in rows})


class TrainingViewOnly(unittest.TestCase):
    """Round 6: LEAK-1, the joint presence of both kinds of test-period data."""

    def run_nlk(self, fit, rows=NLK_FULL, **kwargs):
        return C.invariance_to_test_data(rows, NLK_FEATURES, T, fit, antibiotics=["INH"], **kwargs)

    def test_is_exactly_the_training_view(self):
        recs, feats = C.load_records(NLK_FULL), C._features(dict(NLK_FEATURES))
        del feats["t2"]  # a training row without features stays without them
        training = C.training_ids(recs, CUTOFF, ["INH"])
        view_records, view_features = C.training_view_only(recs, feats, CUTOFF, training)
        view = C._training_view(recs, feats, CUTOFF, training)
        self.assertEqual(view_records, tuple(record for record, _ in view.values()))
        self.assertEqual([r.genome_id for r in view_records], ["t1", "t2", "t3", "t4", "t5", "t6"])
        self.assertEqual(sorted(view_features), ["t1", "t3", "t4", "t5", "t6"])
        self.assertTrue(all(e.recorded_at is not None and e.recorded_at <= CUTOFF for r in view_records for e in r.labels))
        self.assertEqual(C._training_view(C.load_records(view_records), view_features, CUTOFF, training), view)
        self.assertIs(C.PERTURBATIONS["training_view_only"], C.training_view_only)

    def test_fixture_has_both_kinds(self):
        recs = C.load_records(NLK_FULL)
        training = C.training_ids(recs, CUTOFF, ["INH"])
        for name in ("drop_non_training", "drop_post_cutoff_events"):
            self.assertGreater(latest_date(C.PERTURBATIONS[name](recs, {}, CUTOFF, training)[0]), CUTOFF, name)
        self.assertLess(latest_date(C.training_view_only(recs, {}, CUTOFF, training)[0]), CUTOFF)

    def test_leak1_presence_fits_fail(self):
        for name, fit, key in (("VR1 horizon", vr1_horizon, "age"), ("LR2 embargo", lr2_embargo, "prior_r"),
                               ("fin_lr_min prior", lr_min_prior, "prior_r")):
            with self.subTest(name):
                result = self.run_nlk(fit)
                self.assertFalse(result.passed)
                self.assertIn(("training_view_only", key, "changed by test-period records"), result.offenders)
                self.assertEqual({o[0] for o in result.offenders}, {"training_view_only"})
                # Round 5's defaults passed each one: every other perturbation leaves a date after T.
                old = dict(C.PERTURBATIONS)
                del old["training_view_only"]
                self.assertTrue(self.run_nlk(fit, perturbations=old).passed)

    def test_r6_l1_undated_events_on_training_rows_are_dropped(self):
        # The re-check's RK-01 fixture: t1's only test-period data is an undated
        # event, and no event anywhere is dated after T.
        rows = [genome("t1", None, "2025-01-01", [event("INH", "S", "2025-02-01", source="labA"),
                                                  event("INH", "R", None, source="labU")], lineage="L1", group="D1"),
                genome("t2", None, "2025-03-01", [event("INH", "R", "2025-04-01", source="labB")], lineage="L2",
                       group="D2"),
                genome("t3", None, "2025-05-01", [event("INH", "R", "2025-06-01", source="labA")], lineage="L1",
                       group="D3"),
                genome("x1", None, "2026-04-01", [event("INH", "S", "2025-07-01", source="labC")], lineage="L3",
                       group="D4")]
        feats = {"t1": ["katG"], "t2": ["katG", "rpoB"], "t3": ["inhA"], "x1": ["eis"]}
        recs = C.load_records(rows)
        training = C.training_ids(recs, CUTOFF, ["INH"])
        self.assertEqual(training, {"t1", "t2", "t3"})
        view_records = {r.genome_id: r for r in C.training_view_only(recs, C._features(feats), CUTOFF, training)[0]}
        self.assertEqual(view_records["t1"].labels, (C.LabelEvent("INH", "S", dt.date(2025, 2, 1), "labA"),))

        def any_test_period(recs, feats, cutoff):
            rows_by_t = nlk_rows(recs, cutoff)
            ids = {r.genome_id for r in rows_by_t}
            present = (any(r.genome_id not in ids for r in recs)
                       or any(C._later(e, cutoff) for r in recs for e in r.labels))
            out = nlk_model(rows_by_t, feats, cutoff)
            out["prior_r"] = (out["prior_r"], present)  # LEAK-1's shape: folded into a fitted key
            return out

        result = C.invariance_to_test_data(rows, feats, T, any_test_period, antibiotics=["INH"])
        self.assertEqual(result.offenders, (("training_view_only", "prior_r", "changed by test-period records"),))
        honest_fit = C.invariance_to_test_data(
            rows, feats, T, lambda r, f, c: nlk_model(nlk_rows(r, c), f, c), antibiotics=["INH"])
        self.assertTrue(honest_fit.passed, honest_fit)

        # MC-2 (MN-2, MN-3): t1's only test-period data is one RIF event, an
        # antibiotic t1 has no dated label for by T, undated or after T.
        def one_col_per_antibiotic(recs, feats, cutoff):
            return dict(nlk_model(nlk_rows(recs, cutoff), feats, cutoff),
                        cols=sorted({e.antibiotic for r in recs for e in r.labels}))

        for name, recorded in (("RIF undated", None), ("RIF after T", "2026-04-01")):
            with self.subTest(name):
                rif = [dict(rows[0], labels=[event("INH", "S", "2025-02-01", source="labA"),
                                             event("RIF", "R", recorded, source="labU")])] + rows[1:]
                recs = C.load_records(rif)
                view_t1 = C.training_view_only(recs, C._features(feats), CUTOFF, training)[0][0]
                self.assertEqual(view_t1.labels, (C.LabelEvent("INH", "S", dt.date(2025, 2, 1), "labA"),))
                self.assertEqual(C.invariance_to_test_data(rif, feats, T, any_test_period, antibiotics=["INH"]).offenders,
                                 (("training_view_only", "prior_r", "changed by test-period records"),))
                cols = C.invariance_to_test_data(rif, feats, T, one_col_per_antibiotic, antibiotics=["INH"])
                self.assertIn(("training_view_only", "cols", "changed by test-period records"), cols.offenders)
                self.assertTrue(C.invariance_to_test_data(
                    rif, feats, T, lambda r, f, c: nlk_model(nlk_rows(r, c), f, c), antibiotics=["INH"]).passed)

    def test_r6_l2_a_key_only_the_perturbed_fit_returns_is_an_offender(self):
        self.assertEqual(C._changed({"a": "1"}, {"a": "1", "b": "2"}), ["b"])
        self.assertEqual(C._changed({"a": "1", "b": "2"}, {"a": "1"}), ["b"])
        # The re-check's RK-06 fixture and VR2's threshold with one key per gene:
        # on the training view the threshold drops to 1 and only adds keys.
        rows = [genome("t1", None, "2025-01-01", [event("INH", "R", "2025-02-01", source="labA"),
                                                  event("INH", "S", "2026-05-01", "2025-02-01", "labA")],
                       lineage="L1", group="D1"),
                genome("t2", None, "2025-03-01", [event("INH", "R", "2025-04-01", source="labB")], lineage="L2",
                       group="D2"),
                genome("t3", None, "2025-05-01", [event("INH", "S", "2025-06-01", source="labA")], lineage="L1",
                       group="D3"),
                genome("x1", None, "2026-04-01", [event("INH", "S", "2025-07-01", source="labC")], lineage="L3",
                       group="D4")]
        feats = {"t1": ["katG", "rpoB"], "t2": ["katG", "inhA"], "t3": ["katG", "gyrA"], "x1": ["eis"]}

        def per_gene_mincount(recs, feats, cutoff):
            rows_by_t = nlk_rows(recs, cutoff)
            threshold = 2 if latest_date(recs) > cutoff else 1
            seen = collections.Counter(g for r in rows_by_t for g in feats.get(r.genome_id, ()))
            y = {r.genome_id: C.label_as_of(r, "INH", cutoff).value for r in rows_by_t}
            out = {"prior_r": sum(v == "R" for v in y.values()) / len(rows_by_t)}
            for gene in sorted(g for g, n in seen.items() if n >= threshold):
                carriers = [y[r.genome_id] for r in rows_by_t if gene in feats.get(r.genome_id, ())]
                out[f"w:{gene}"] = (carriers.count("R"), carriers.count("S"))
            return out

        result = C.invariance_to_test_data(rows, feats, T, per_gene_mincount, antibiotics=["INH"])
        self.assertEqual(result.offenders, tuple(("training_view_only", f"w:{g}", "changed by test-period records")
                                                 for g in ("gyrA", "inhA", "rpoB")))

    def test_documented_limits_i1_i2(self):
        # I-1: a bare value is the key "<fit>". Since MC2-1 a tuple and a list differ.
        self.assertNotEqual(C.fingerprint((1, 2)), C.fingerprint([1, 2]))
        # Since MC3-1 the bare value's key is reserved, so the two forms never meet.
        self.assertEqual(C._fingerprints(5).keys(), {"<fit>"})
        with self.assertRaisesRegex(ValueError, "reserved"):
            C._fingerprints({"<fit>": 5})
        # I-2: a cutoff of 0001-01-01 raises OverflowError, never a pass.
        rows = [genome("t1", None, "0001-01-01", [event("INH", "S", "0001-01-01")])]
        honest_fit = lambda r, f, c: nlk_model(nlk_rows(r, c), f, c)  # noqa: E731
        with self.assertRaises(OverflowError):
            C.invariance_to_test_data(rows, {"t1": ["katG"]}, "0001-01-01", honest_fit, antibiotics=["INH"])
        # MC-3: with no training row it fails without raising, and without
        # new_genomes it runs, so an honest fit passes.
        unlabelled = [genome("t1", None, "0001-01-01", [])]
        self.assertEqual(C.invariance_to_test_data(unlabelled, {"t1": ["katG"]}, "0001-01-01",
                                                   lambda r, f, c: {"n": len(nlk_rows(r, c))},
                                                   antibiotics=["INH"]).offenders,
                         (("control", "*", "no training rows at the cutoff; invariance not certified"),))
        two = rows + [genome("t2", None, "0001-01-01", [event("INH", "R", "0001-01-01")])]
        without = {k: v for k, v in C.PERTURBATIONS.items() if k != "new_genomes"}
        self.assertTrue(C.invariance_to_test_data(two, {"t1": ["katG"], "t2": ["rpoB"]}, "0001-01-01", honest_fit,
                                                  perturbations=without, antibiotics=["INH"]).passed)

    def test_mc1_masked_and_structured_arrays(self):
        # Pairs the round-6 fingerprint could not tell apart.
        data = np.array([1.0, 2.0, 3.0])
        self.assertNotEqual(C.fingerprint(np.zeros(2, dtype=[("katG", "f8")])),
                            C.fingerprint(np.zeros(2, dtype=[("rpoB", "f8")])))
        self.assertNotEqual(C.fingerprint(np.array(5.0)), C.fingerprint(np.array([5.0])))
        self.assertNotEqual(C.fingerprint(np.zeros(1, dtype=[("katG", "f8")])[0]),
                            C.fingerprint(np.zeros(1, dtype=[("rpoB", "f8")])[0]))
        # Any ndarray subclass is refused, numpy.ma.MaskedArray included since
        # MC4, here and through the canary.
        class Tagged(np.ndarray):
            pass
        matrix = data.reshape(1, 3).view(np.matrix)  # a view: np.asmatrix warns under -W error
        for value in (matrix, data.view(Tagged), np.ma.MaskedArray(data, mask=[False, True, False])):
            with self.subTest(type(value).__name__), self.assertRaisesRegex(TypeError, "an ndarray subclass"):
                C.fingerprint(value)
        with self.assertRaisesRegex(TypeError, "an ndarray subclass"):
            C.invariance_to_test_data(NLK_FULL, NLK_FEATURES, T, lambda r, f, c: {"m": matrix},
                                      antibiotics=["INH"])

        def training_ids_of(recs, cutoff):
            return {r.genome_id for r in nlk_rows(recs, cutoff)}

        def masked(leak):
            def fit(recs, feats, cutoff):
                rows = nlk_rows(recs, cutoff)
                out = nlk_model(rows, feats, cutoff)
                genes = out["vocabulary"]
                counts = np.array([out["gene_weight"][g] for g in genes], dtype=float)
                if leak:  # screening on genomes outside the training rows
                    ids = training_ids_of(recs, cutoff)
                    other = set().union(*(g for gid, g in feats.items() if gid not in ids))
                    mask = [g in other for g in genes]
                else:  # honest: masked where no training row is R
                    mask = [c == 0 for c in counts]
                return dict(out, w=np.ma.MaskedArray(counts, mask=mask))
            return fit

        def structured(leak):
            def fit(recs, feats, cutoff):
                rows = nlk_rows(recs, cutoff)
                out = nlk_model(rows, feats, cutoff)
                pool = feats if leak else {r.genome_id: feats.get(r.genome_id, ()) for r in rows}
                last = max(set().union(*pool.values()))  # the field is named after one gene
                return dict(out, top=np.array([(out["prior_r"],)], dtype=[(last, "f8")]))
            return fit

        # The masked fit is refused, leaking and honest (MC4); the structured
        # leak fails, and its honest twin passes.
        for leak in (True, False):
            with self.subTest("masked", leak=leak), self.assertRaisesRegex(TypeError, "an ndarray subclass"):
                self.run_nlk(masked(leak))
        result = self.run_nlk(structured(True))
        self.assertFalse(result.passed)
        self.assertIn(("training_view_only", "top", "changed by test-period records"), result.offenders)
        self.assertEqual({o[1] for o in result.offenders}, {"top"})
        self.assertTrue(self.run_nlk(structured(False)).passed)

    def test_mc2_2_masked_and_structured_details(self):
        # MX-1 and MX-2: since MC4 a masked array is refused, whatever it hides.
        for value in (np.ma.MaskedArray([1.0, 2.0], mask=[False, True]), np.ma.MaskedArray([1.0, 1.0], mask=[True, False]),
                      np.ma.masked_array([1.0], fill_value=-1.0), np.ma.MaskedArray([1.0], hard_mask=True)):
            with self.subTest(repr(value)), self.assertRaisesRegex(TypeError, "an ndarray subclass"):
                C.fingerprint(value)
        # MX-3: each field's values, not just the names and dtypes.
        self.assertNotEqual(C.fingerprint(np.zeros(2, dtype=[("katG", "f8")])),
                            C.fingerprint(np.ones(2, dtype=[("katG", "f8")])))
        # Titles and offsets are part of a structured dtype.
        self.assertNotEqual(C.fingerprint(np.zeros(1, dtype=[(("zz1", "w"), "f8")])),
                            C.fingerprint(np.zeros(1, dtype=[("w", "f8")])))
        self.assertNotEqual(C.fingerprint(np.zeros(1, dtype={"names": ["w"], "formats": ["f8"], "offsets": [8],
                                                             "itemsize": 16})),
                            C.fingerprint(np.zeros(1, dtype={"names": ["w"], "formats": ["f8"], "offsets": [0],
                                                             "itemsize": 16})))

        # MX-4 and MC2-1: a MaskedArray subclass, a masked array over an
        # ndarray subclass, and one with an added attribute are refused too.
        class SubMasked(np.ma.MaskedArray):
            pass

        class Tagged(np.ndarray):
            pass

        data = np.array([1.0, 2.0, 3.0])
        tagged_attr = np.ma.MaskedArray(data)
        tagged_attr.screened = ["zz1"]
        for name, value in (("MaskedArray subclass", SubMasked(data)),
                            ("masked over Tagged", np.ma.MaskedArray(data.view(Tagged))),
                            ("masked over matrix", np.ma.MaskedArray(data.reshape(1, 3).view(np.matrix))),
                            ("masked with an attribute", tagged_attr)):
            with self.subTest(name), self.assertRaisesRegex(TypeError, "cannot fingerprint"):
                C.fingerprint(value)

    def test_mc2_1_fingerprint_fails_closed(self):
        # The exact type is part of the form.
        forms = [[1], (1,), {0: 1}, {1}, frozenset([1]), 1, True, 1.0, "1", b"1", np.int64(1), np.float64(1.0)]
        self.assertEqual(len({C.fingerprint(v) for v in forms}), len(forms))
        # Floats keep every bit, NaN sign and payload included.
        nan = float("nan")
        neg_nan = struct.unpack(">d", bytes.fromhex("fff8000000000000"))[0]
        payload = struct.unpack(">d", bytes.fromhex("7ff8000000000005"))[0]
        self.assertEqual(len({C.fingerprint(v) for v in (nan, neg_nan, payload)}), 3)
        self.assertEqual(len({C.fingerprint(np.float64(v)) for v in (nan, neg_nan, payload)}), 3)
        self.assertNotEqual(C.fingerprint(0.0), C.fingerprint(-0.0))
        # numpy scalars are read by their bytes, not item(): datetime64 is
        # accepted, and a float32 NaN keeps its own payload.
        self.assertNotEqual(C.fingerprint(np.datetime64("2026-01-01")), C.fingerprint(np.datetime64("2026-01-02")))
        # Scalars whose item() is a refused type still fingerprint, by their bytes.
        for value in (np.complex64(1 + 2j), np.datetime64("2026-01-01T00:00:00"), np.timedelta64(5, "s")):
            with self.subTest(type(value).__name__):
                self.assertEqual(len(C.fingerprint(value)), 64)
        f32 = np.frombuffer(bytes.fromhex("7fc00001"), dtype=">f4").astype("<f4")[0]
        f32b = np.frombuffer(bytes.fromhex("7fc00002"), dtype=">f4").astype("<f4")[0]
        self.assertNotEqual(C.fingerprint(f32), C.fingerprint(f32b))
        # Huge ints fingerprint; longdouble is refused with TypeError.
        self.assertNotEqual(C.fingerprint(10 ** 5000), C.fingerprint(10 ** 5000 + 1))
        with self.assertRaisesRegex(TypeError, "longdouble"):
            C.fingerprint(np.longdouble(1))

        @dataclasses.dataclass
        class Point:
            x: float

        Top = collections.namedtuple("Top", ["zz1"])

        class Attr(dict):
            pass

        refused = {"namedtuple": Top(0.5), "defaultdict": collections.defaultdict(float, {"a": 1.0}),
                   "OrderedDict": collections.OrderedDict(a=1.0), "dict subclass": Attr(a=1.0),
                   "dataclass": Point(1.0), "matrix": np.zeros((1, 2)).view(np.matrix), "datetime": dt.datetime(2026, 1, 1),
                   "complex": 1j, "object": object(), "nested namedtuple": [{"k": Top(0.5)}]}
        for name, value in refused.items():
            with self.subTest(name), self.assertRaisesRegex(TypeError, "cannot fingerprint"):
                C.fingerprint(value)

    def test_mc3_1_keys_masked_flags_and_dtype_metadata(self):
        def outside(recs, feats, cutoff):
            ids = {r.genome_id for r in nlk_rows(recs, cutoff)}
            return sorted(set().union(*(set(g) for gid, g in feats.items() if gid not in ids)))

        class Quiet:  # a hashable key whose str() hides what it carries
            def __init__(self, genes):
                self.genes = tuple(genes)

            def __str__(self):
                return "screen"

            def __hash__(self):
                return 0

            def __eq__(self, other):
                return isinstance(other, Quiet)

        def lying_key(leak):
            def fit(recs, feats, cutoff):
                out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
                out[Quiet(outside(recs, feats, cutoff) if leak else [])] = 1
                return out
            return fit

        def typed_key(leak):  # the key's type alone carries the bit
            def fit(recs, feats, cutoff):
                out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
                out[7 if leak and outside(recs, feats, cutoff) else "7"] = 1
                return out
            return fit

        def shared_mask(leak, attr):
            def fit(recs, feats, cutoff):
                out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
                w = np.ma.MaskedArray(np.array([out["prior_r"]]), mask=[False])
                genes = outside(recs, feats, cutoff) if leak else []
                setattr(w, attr, genes if attr != "_sharedmask_bool" else bool(genes))
                return dict(out, w=w)
            return fit

        def metadata(leak, structured):
            def fit(recs, feats, cutoff):
                out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
                dtype = np.dtype("f8", metadata={"screened": outside(recs, feats, cutoff) if leak else []})
                if structured:
                    dtype = np.dtype([("p", dtype)])
                return dict(out, w=np.zeros(2, dtype))
            return fit

        for name, make in (("lying __str__ key", lying_key), ("_sharedmask from the genes",
                                                              lambda leak: shared_mask(leak, "_sharedmask")),
                           ("_isfield from the genes", lambda leak: shared_mask(leak, "_isfield")),
                           ("dtype metadata", lambda leak: metadata(leak, False)),
                           ("field dtype metadata", lambda leak: metadata(leak, True))):
            for leak in (True, False):
                with self.subTest(name, leak=leak), self.assertRaises(TypeError):
                    self.run_nlk(make(leak))
        # A key of the wrong type is refused even when the honest twin's key is a str.
        with self.assertRaisesRegex(TypeError, "exact str"):
            self.run_nlk(typed_key(True))
        self.assertTrue(self.run_nlk(typed_key(False), constant_keys=["7"]).passed)
        # MC4-1: a masked array's state (flags, fill_value, _basedict, _optinfo)
        # is no longer read at all: the array is refused.
        b = np.ma.MaskedArray(np.array([1.0]), mask=[False])
        b._optinfo["screened"] = ["zz1"]
        c = np.ma.MaskedArray(np.array([1.0]), mask=[False], fill_value=np.ma.MaskedArray([7.0]))
        for value in (b, c):
            with self.subTest(repr(value)), self.assertRaisesRegex(TypeError, "an ndarray subclass"):
                C.fingerprint(value)
        # Sub-array and nested fields with metadata are refused too.
        meta = np.dtype("f8", metadata={"k": 1})
        for dtype in (np.dtype([("s", meta, (2,))]), np.dtype([("o", [("i", meta)])])):
            with self.subTest(str(dtype)), self.assertRaisesRegex(TypeError, "dtype metadata"):
                C.fingerprint(np.zeros(1, dtype))
        # MC3-4: a bad _baseclass raises TypeError, not AttributeError.
        e = np.ma.MaskedArray(np.array([1.0]))
        e._baseclass = 5
        with self.assertRaisesRegex(TypeError, "cannot fingerprint"):
            C.fingerprint(e)

    def test_mc3_2_sets_and_scalars(self):
        # MY-3: set elements by their canonical form, refusals included.
        nan = float("nan")
        neg_nan = struct.unpack(">d", bytes.fromhex("fff8000000000000"))[0]
        self.assertNotEqual(C.fingerprint({nan}), C.fingerprint({neg_nan}))
        self.assertNotEqual(C.fingerprint({1.0}), C.fingerprint({np.float64(1.0)}))
        with self.assertRaisesRegex(TypeError, "cannot fingerprint"):
            C.fingerprint({collections.namedtuple("Top", ["a"])(1)})
        # MY-5: a numpy scalar's dtype is part of its form.
        self.assertNotEqual(C.fingerprint(np.int64(0)), C.fingerprint(np.float64(0.0)))
        # Same type and bytes, different unit: only the dtype tells them apart.
        self.assertNotEqual(C.fingerprint(np.datetime64(1, "D")), C.fingerprint(np.datetime64(1, "s")))
        self.assertNotEqual(C.fingerprint(np.int32(1)), C.fingerprint(np.frombuffer(np.int32(1).tobytes(), "f4")[0]))
        # numpy scalar aliases are different types, so they differ too.
        self.assertNotEqual(C.fingerprint(np.longlong(5)), C.fingerprint(np.int64(5)))
        # I-3, documented limits: memory layout is not part of an array's values.
        a = np.arange(6.0).reshape(2, 3)
        self.assertEqual(C.fingerprint(a), C.fingerprint(np.asfortranarray(a)))
        self.assertEqual(C.fingerprint(a), C.fingerprint(np.arange(12.0)[:6].reshape(2, 3)))
        x6 = np.arange(12.0)
        self.assertEqual(C.fingerprint(x6[::2]), C.fingerprint(x6[::2].copy()))  # strides
        frozen = a.copy()
        frozen.flags.writeable = False
        self.assertEqual(C.fingerprint(a), C.fingerprint(frozen))
        packed = np.dtype({"names": ["a", "b"], "formats": ["u1", "f8"], "offsets": [0, 8], "itemsize": 16})
        aligned = np.dtype([("a", "u1"), ("b", "f8")], align=True)
        self.assertNotEqual(packed.isalignedstruct, aligned.isalignedstruct)
        self.assertEqual(C.fingerprint(np.zeros(1, packed)), C.fingerprint(np.zeros(1, aligned)))
        padded = np.dtype({"names": ["w"], "formats": ["f8"], "offsets": [0], "itemsize": 16})
        x, y = np.zeros(1, padded), np.zeros(1, padded)
        y.view(np.uint8)[8] = 7  # a padding byte
        self.assertEqual(C.fingerprint(x), C.fingerprint(y))
        # MY-6: sets sorted, whatever order they iterate in.
        words = [f"g{i:03d}" for i in range(200)]
        forward, backward = frozenset(words), frozenset(reversed(words))
        if list(forward) == list(backward):  # make the iteration orders differ
            backward = frozenset(list(reversed(list(forward))))
        self.assertEqual(C.fingerprint(forward), C.fingerprint(backward))
        self.assertEqual(C.fingerprint(set(words)), C.fingerprint(set(reversed(words))))

    def test_mc4_dtype_aliases_void_scalars_and_keys(self):
        # MC4-2b: dtype aliases differ, on arrays, fields, sub-array fields and scalars.
        for a, b in (("q", "l"), ("Q", "L")):
            with self.subTest(a=a, b=b):
                self.assertEqual(np.dtype(a).str, np.dtype(b).str)
                self.assertNotEqual(C.fingerprint(np.zeros(2, a)), C.fingerprint(np.zeros(2, b)))
                self.assertNotEqual(C.fingerprint(np.zeros(1, [("w", a)])), C.fingerprint(np.zeros(1, [("w", b)])))
                self.assertNotEqual(C.fingerprint(np.zeros(1, [("w", a, (2,))])),
                                    C.fingerprint(np.zeros(1, [("w", b, (2,))])))
                self.assertNotEqual(C.fingerprint(np.zeros(1, a)[0]), C.fingerprint(np.zeros(1, b)[0]))
        # MC4-2a: an unstructured void scalar is refused, and metadata on a
        # structured void scalar is refused.
        with self.assertRaisesRegex(TypeError, "unstructured void scalar"):
            C.fingerprint(np.zeros(1, "V8")[0])
        with self.assertRaisesRegex(TypeError, "dtype metadata"):
            C.fingerprint(np.zeros(1, np.dtype([("p", "f8")], metadata={"screened": ["zz1"]}))[0])
        self.assertEqual(len(C.fingerprint(np.zeros(1, [("p", "f8")])[0])), 64)
        # MZ-7: metadata on a sub-array field's own dtype (its base and values carry none).
        sub = np.dtype(("f8", (2,)), metadata={"screened": ["zz1"]})
        with self.assertRaisesRegex(TypeError, "dtype metadata"):
            C.fingerprint(np.zeros(1, [("s", sub)]))
        # MZ-1: a str subclass is not an exact str key.
        class Key(str):
            pass
        with self.assertRaisesRegex(TypeError, "exact str"):
            C._fingerprints({Key("w"): 1})
        with self.assertRaisesRegex(TypeError, "exact str"):
            self.run_nlk(lambda r, f, c: dict(nlk_model(nlk_rows(r, c), f, c), **{Key("w"): 1}))

    def test_mc2_1_leaking_fits_of_refused_types(self):
        try:
            from scipy.sparse import dok_matrix
        except ImportError:  # the suite's environment has no scipy: a dict subclass with a shape stands in
            class dok_matrix(dict):  # noqa: N801
                def __init__(self, shape):
                    super().__init__()
                    self.shape = shape

        class Tagged(np.ndarray):
            pass

        def outside(recs, feats, cutoff):  # genes that genomes outside the training rows carry
            ids = {r.genome_id for r in nlk_rows(recs, cutoff)}
            return sorted(set().union(*(set(g) for gid, g in feats.items() if gid not in ids)))

        def gene_rows(recs, feats, cutoff, leak):
            out = nlk_model(nlk_rows(recs, cutoff), feats, cutoff)
            return out, (outside(recs, feats, cutoff) if leak else [])

        def dok(leak):
            def fit(recs, feats, cutoff):
                out, screened = gene_rows(recs, feats, cutoff, leak)
                m = dok_matrix((1, len(out["vocabulary"]) + len(screened)))
                m[0, 0] = out["prior_r"]
                return dict(out, w=m)
            return fit

        def named(leak):
            def fit(recs, feats, cutoff):
                out, screened = gene_rows(recs, feats, cutoff, leak)
                return dict(out, top=collections.namedtuple("Top", [(screened or out["vocabulary"])[-1]])(out["prior_r"]))
            return fit

        def default(leak):
            def fit(recs, feats, cutoff):
                out, screened = gene_rows(recs, feats, cutoff, leak)
                return dict(out, gene_weight=collections.defaultdict(functools.partial(float, len(screened)),
                                                                     out["gene_weight"]))
            return fit

        def masked_over_subclass(leak):
            def fit(recs, feats, cutoff):
                out, screened = gene_rows(recs, feats, cutoff, leak)
                base = np.array([out["prior_r"]]).view(Tagged)
                base.tag = screened
                return dict(out, w=np.ma.MaskedArray(base))
            return fit

        def top_level(kind):
            def make(leak):
                def fit(recs, feats, cutoff):
                    out, screened = gene_rows(recs, feats, cutoff, leak)
                    if kind is collections.defaultdict:
                        return collections.defaultdict(functools.partial(float, len(screened)), out)
                    return kind(out)
                return fit
            return make

        for name, make in (("dok_matrix", dok), ("namedtuple", named), ("defaultdict", default),
                           ("masked over a subclass", masked_over_subclass),
                           ("top-level defaultdict", top_level(collections.defaultdict)),
                           ("top-level OrderedDict", top_level(collections.OrderedDict))):
            for leak in (True, False):
                with self.subTest(name, leak=leak), self.assertRaisesRegex(TypeError, "cannot fingerprint"):
                    self.run_nlk(make(leak))

    def test_leak1_clean_fits_pass(self):
        for name, fit, rows in (("clean horizon", clean_horizon, NLK_FULL), ("clean horizon PARTIAL", clean_horizon,
                                                                            NLK_PARTIAL),
                                ("honest", lambda recs, feats, cutoff: nlk_model(nlk_rows(recs, cutoff), feats, cutoff),
                                 NLK_FULL)):
            with self.subTest(name):
                self.assertTrue(self.run_nlk(fit, rows).passed)


class UnseenVocabulary(unittest.TestCase):
    TEST = {gid: FEATURES[gid] for gid in ("g05", "g06")}

    def vocab(self):
        return registry().vocabulary(T)

    def test_clean_under_both_policies(self):
        for policy in ("unknown", "drop"):
            with self.subTest(policy):
                vocab = self.vocab()
                result = C.unseen_vocabulary(vocab, self.TEST, policy=policy)
                self.assertTrue(result.passed, result)
                self.assertEqual(vocab, registry().vocabulary(T))
        self.assertEqual(A.encode_genes(FEATURES["g05"], self.vocab()), ("katG", A.UNKNOWN))
        self.assertEqual(A.encode_genes(FEATURES["g05"], self.vocab(), "drop"), ("katG",))

    def test_vocabulary_growth_leak(self):
        def grows(genes, vocab):
            vocab.update(genes)
            return tuple(sorted(set(genes)))
        result = C.unseen_vocabulary(self.vocab(), self.TEST, grows)
        self.assertFalse(result.passed)
        self.assertIn(("<vocabulary>", "embB", "added while encoding test genomes"), result.offenders)
        self.assertIn(("<vocabulary>", "geneC", "added while encoding test genomes"), result.offenders)
        self.assertIn(("g05", "embB", "token outside the vocabulary"), result.offenders)

    def test_passthrough_and_order_leaks(self):
        def passthrough(genes, vocab):
            return tuple(genes)
        result = C.unseen_vocabulary(self.vocab(), self.TEST, passthrough)
        self.assertEqual(sorted(set(result.offenders)), sorted(result.offenders))
        self.assertIn(("g05", "*", "encoding depends on gene order or repetition"), result.offenders)
        self.assertIn(("g05", "geneC", "token outside the vocabulary"), result.offenders)
        self.assertIn(("g05", "*", f"2 unseen gene(s) not mapped to {A.UNKNOWN}"), result.offenders)
        # g06 has no unseen genes, but passing input order through is still order-dependent.
        self.assertEqual([o for o in result.offenders if o[0] == "g06"],
                         [("g06", "*", "encoding depends on gene order or repetition")])

    def test_policy_mismatch_leaks(self):
        drops = lambda genes, vocab: A.encode_genes(genes, vocab, "drop")  # noqa: E731
        result = C.unseen_vocabulary(self.vocab(), self.TEST, drops, policy="unknown")
        self.assertEqual(result.offenders, (("g05", "*", f"2 unseen gene(s) not mapped to {A.UNKNOWN}"),))
        maps = lambda genes, vocab: A.encode_genes(genes, vocab, "unknown")  # noqa: E731
        result = C.unseen_vocabulary(self.vocab(), self.TEST, maps, policy="drop")
        self.assertEqual(result.offenders, (("g05", "*", f"{A.UNKNOWN} emitted under policy 'drop' "
                                                         "with 2 unseen gene(s)"),))
        forgets = lambda genes, vocab: ()  # noqa: E731
        result = C.unseen_vocabulary(self.vocab(), self.TEST, forgets, policy="drop")
        self.assertIn(("g05", "katG", "known gene not encoded"), result.offenders)

    def test_n4_known_tokens_the_genome_does_not_carry(self):
        # The reviewer's B5a/B5b, under policy "drop".
        vocab = frozenset({"katG", "rpoB", "inhA"})

        def nearest(genes, voc):  # unseen allele remapped onto its known gene
            return tuple(sorted({g if g in voc else g.split("_")[0] for g in genes
                                 if g in voc or g.split("_")[0] in voc}))

        def phantom(genes, voc):  # a known gene added to every genome
            return tuple(sorted({g for g in genes if g in voc} | {"inhA"}))
        self.assertTrue(C.unseen_vocabulary(vocab, {"x9": ["rpoB_S450L"], "x2": ["rpoB"]}, policy="drop").passed)
        self.assertEqual(C.unseen_vocabulary(vocab, {"x9": ["rpoB_S450L"]}, nearest, policy="drop").offenders,
                         (("x9", "rpoB", "known token the genome does not carry"),))
        self.assertEqual(C.unseen_vocabulary(vocab, {"x2": ["rpoB"]}, phantom, policy="drop").offenders,
                         (("x2", "inhA", "known token the genome does not carry"),))
        # Under "unknown" the remap is caught twice: a phantom known token and no UNKNOWN.
        self.assertEqual({o[2] for o in C.unseen_vocabulary(vocab, {"x9": ["rpoB_S450L"]}, nearest).offenders},
                         {"known token the genome does not carry", f"1 unseen gene(s) not mapped to {A.UNKNOWN}"})

    def test_nn3_tokens_are_a_multiset(self):
        # The re-checker's G6: the unseen allele becomes a second copy of a gene the genome carries.
        vocab = frozenset({"katG", "rpoB", "inhA"})

        def duplicate(genes, voc):
            out = [g for g in sorted(set(genes)) if g in voc]
            return tuple(sorted(out + ["rpoB" for g in genes if g == "rpoB_S450L" and "rpoB" in genes]))

        def double_unknown(genes, voc):
            return A.encode_genes(genes, voc) + (A.UNKNOWN,)
        x4 = {"x4": ["rpoB", "rpoB_S450L"]}
        self.assertTrue(C.unseen_vocabulary(vocab, x4, policy="drop").passed)
        self.assertEqual(C.unseen_vocabulary(vocab, x4, duplicate, "drop").offenders, (("x4", "rpoB", "token repeated"),))
        self.assertEqual(C.unseen_vocabulary(vocab, x4, double_unknown).offenders,
                         (("x4", A.UNKNOWN, "token repeated"),))


M_SPLITS = (ROOT / "scripts" / "temporal" / "splits.py",
            Path("/home/layth/review/grm-task-m-temporal-splits-20261030/scripts/temporal/splits.py"))


def load_m():
    """Task M's splits.py, read-only and without writing bytecode beside it:
    the merged copy when present, else M's worktree."""
    path = next((p for p in M_SPLITS if p.is_file()), None)
    if path is None:
        return None, None
    key = "grm_temporal_splits_crosscheck"
    module = sys.modules.get(key)
    if module is None:
        spec = importlib.util.spec_from_file_location(key, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        saved, sys.dont_write_bytecode = sys.dont_write_bytecode, True
        try:
            spec.loader.exec_module(module)
        finally:
            sys.dont_write_bytecode = saved
    return module, path


M, M_PATH = load_m()
BASE = dt.date(2025, 1, 1)
ONE_DAY = dt.timedelta(days=1)
CROSS_CHECK = {}


def day(n):
    return (BASE + dt.timedelta(days=n)).isoformat()


def synthetic_rows(seed, n=60):
    """Seeded records that hit all five rulings: undated events, same-day
    agreeing and conflicting events, supersedes chains, both availability
    orders and omitted optional keys. Days come from a small pool per
    antibiotic so ties are common."""
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        collected = rng.choice([None, rng.randrange(0, 200)])
        labels = []
        for antibiotic in ("INH", "RIF"):
            pool = [rng.randrange(collected or 0, 600) for _ in range(3)]
            used, dated = set(), []
            for _ in range(rng.randrange(0, 6)):
                when = None if rng.random() < 0.15 else rng.choice(pool)
                source = f"src-{rng.choice('abcd')}"
                if (when, source) in used:
                    continue
                used.add((when, source))
                earlier = [d for d in dated if when is None or d < when]
                supersedes = rng.choice(earlier) if earlier and rng.random() < 0.4 else None
                labels.append({"antibiotic": antibiotic, "value": rng.choice("RS"),
                               "recorded_at": None if when is None else day(when), "source_release": source,
                               "supersedes": None if supersedes is None else day(supersedes)})
                if when is not None:
                    dated.append(when)
        rng.shuffle(labels)
        row = {"genome_id": f"s{seed}-{i:03d}", "released_at": rng.choice([None, day(rng.randrange(0, 500))]),
               "lineage": rng.choice([None, "L1", " L2 "]), "labels": labels}
        if collected is not None or rng.random() < 0.5:
            row["collection_date"] = None if collected is None else day(collected)
        if rng.random() < 0.6:
            row["features_available_at"] = rng.choice([None, day(rng.randrange(0, 500))])
        rows.append(row)
    return rows


def outcome(module, record, antibiotic, cutoff):
    try:
        event = module.label_as_of(record, antibiotic, cutoff)
    except module.AmbiguousLabel:
        return "ambiguous"
    if event is None:
        return None
    return (event.value, event.recorded_at.isoformat(), event.source_release,
            event.supersedes and event.supersedes.isoformat())


def m_train(rows, cutoff, antibiotic="INH"):
    cutoff = C._date(cutoff, "cutoff", False)
    return M.split_at(rows, cutoff, 90, 0, cutoff + dt.timedelta(days=365), antibiotic)["train"]


@unittest.skipIf(M is None, "Task M's splits.py is neither merged nor in its worktree")
class CrossCheckWithM(unittest.TestCase):
    """N's loader and label_as_of against M's splits.py on the same rows."""

    def test_label_as_of_matches_on_seeded_records(self):
        compared, kinds = 0, collections.Counter()
        for seed in range(10):
            rows = synthetic_rows(seed)
            mine = {r.genome_id: r for r in C.load_records(rows)}
            theirs = {r.genome_id: r for r in (M.GenomeRecord.from_mapping(row) for row in rows)}
            self.assertEqual(sorted(mine), sorted(theirs))
            for gid, record in mine.items():
                other = theirs[gid]
                self.assertEqual(record.available_at, other.available_at, gid)
                self.assertEqual(record.inconsistent_dates, other.inconsistent_dates, gid)
                self.assertEqual(record.lineage, other.lineage, gid)
                kinds["inconsistent_dates"] += record.inconsistent_dates
                kinds["undated_events"] += sum(e.recorded_at is None for e in record.labels)
                days = sorted({e.recorded_at for e in record.labels if e.recorded_at is not None})
                cutoffs = sorted({BASE + dt.timedelta(days=n) for n in range(-1, 620, 9)}
                                 | set(days) | {d - ONE_DAY for d in days})
                for antibiotic in ("INH", "RIF", "EMB"):
                    for cutoff in cutoffs:
                        got, want = outcome(C, record, antibiotic, cutoff), outcome(M, other, antibiotic, cutoff)
                        self.assertEqual(got, want, (gid, antibiotic, cutoff))
                        compared += 1
                        known = C.known_events(record, antibiotic, cutoff)
                        tied = [e for e in known if e.recorded_at == known[-1].recorded_at] if known else []
                        if got == "ambiguous":
                            kinds["ambiguous"] += 1
                        elif got is None:
                            kinds["none"] += 1
                        else:
                            kinds["label"] += 1
                            kinds["agreeing_tie_read"] += len(tied) > 1
                            earlier = {d for e in known for d in [e.recorded_at] if d < known[-1].recorded_at}
                            kinds["superseded_conflict_read"] += any(
                                len({e.value for e in known if e.recorded_at == d}) > 1 for d in earlier)
        for kind in ("ambiguous", "none", "label", "agreeing_tie_read", "superseded_conflict_read",
                     "undated_events", "inconsistent_dates"):
            self.assertGreater(kinds[kind], 0, (kind, kinds))
        self.assertGreater(compared, 10000)
        CROSS_CHECK["label_as_of"] = {"compared": compared, **kinds}

    def test_load_verdicts_match(self):
        base = genome("x", "2025-05-01", "2025-06-01", [event("INH", "R", "2025-09-01")])
        no_source = {k: v for k, v in event("INH", "R", "2025-09-01").items() if k != "source_release"}
        cases = {
            "undated event": dict(base, labels=[event("INH", "R", None)]),
            "empty labels": dict(base, labels=[]),
            "no labels key": {k: v for k, v in base.items() if k != "labels"},
            "empty source_release": dict(base, labels=[event("INH", "R", "2025-09-01", source="")]),
            "no source_release key": dict(base, labels=[no_source]),
            "same day agreeing": dict(base, labels=[event("INH", "R", "2025-09-01", source=s) for s in "ab"]),
            "same day conflict": dict(base, labels=[event("INH", v, "2025-09-01", source=v) for v in "RS"]),
            "supersedes missing day": dict(base, labels=[event("INH", "R", "2025-09-01", "2025-08-01")]),
            "supersedes other antibiotic": dict(base, labels=[event("INH", "R", "2025-09-01"),
                                                              event("RIF", "S", "2025-10-01", "2025-09-01")]),
            "supersedes undated only": dict(base, labels=[event("INH", "R", None),
                                                          event("INH", "S", "2025-10-01", "2025-09-01")]),
            "undated supersedes dated": dict(base, labels=[event("INH", "R", "2025-09-01"),
                                                           event("INH", "S", None, "2025-09-01")]),
            "supersedes not earlier": dict(base, labels=[event("INH", "R", "2025-09-01", "2025-09-01")]),
            "label before collection": dict(base, labels=[event("INH", "R", "2025-04-30")]),
            "label on collection day": dict(base, labels=[event("INH", "R", "2025-05-01")]),
            "features before release": dict(base, features_available_at="2025-01-01"),
            "whitespace lineage": dict(base, lineage="  "),
            "unknown key": dict(base, note="x"),
            "datetime date": dict(base, released_at="2025-06-01T00:00:00"),
        }
        verdicts = {}
        for name, row in cases.items():
            with self.subTest(name):
                results = []
                for load in (C.load_record, M.GenomeRecord.from_mapping):
                    try:
                        load(row)
                        results.append("accept")
                    except ValueError:
                        results.append("refuse")
                self.assertEqual(results[0], results[1], name)
                verdicts[name] = results[0]
        self.assertEqual(sum(v == "refuse" for v in verdicts.values()), 11, verdicts)
        CROSS_CHECK["load_verdicts"] = verdicts

    def test_train_matches_m_split(self):
        compared = 0
        for seed in range(10):
            rows = synthetic_rows(seed)
            for n in range(0, 600, 30):
                cutoff = day(n)
                mine = {r["genome_id"]: r["value"] for r in C.training_labels_as_of(rows, "INH", cutoff)}
                self.assertEqual(mine, m_train(rows, cutoff), (seed, cutoff))
                compared += 1
        self.assertEqual(m_train(RECORDS, T), {r["genome_id"]: r["value"] for r in clean_training()})
        CROSS_CHECK["train_splits"] = compared

    def test_m_split_passes_invariance(self):
        def fit(recs, features, cutoff):
            return {"train": m_train([r.to_mapping() for r in recs], cutoff)}
        results = [C.invariance_to_test_data(RECORDS, FEATURES, T, fit, antibiotics=("INH",)),
                   C.invariance_to_test_data(REVIEW_RECORDS, REVIEW_FEATURES, T, fit, antibiotics=("INH",))]
        results += [C.invariance_to_test_data(synthetic_rows(seed), {}, day(300), fit, antibiotics=("INH",))
                    for seed in range(3)]
        for result in results:
            self.assertTrue(result.passed, result)
        CROSS_CHECK["m_split_invariance"] = [r.passed for r in results]


def all_results(seed=None):
    """Every canary's clean control and planted leak, as JSON-ready dicts."""
    rows, feats, train = list(RECORDS), dict(FEATURES), clean_training()
    if seed is not None:
        rng = random.Random(seed)
        rng.shuffle(rows)
        rng.shuffle(train)
        feats = dict(rng.sample(list(feats.items()), len(feats)))
        feats = {gid: list(reversed(genes)) for gid, genes in feats.items()}
    reg = registry()
    leaked_rev = [dict(r, value="R") if r["genome_id"] == "g01" else r for r in train]
    test = {gid: feats[gid] for gid in ("g05", "g06")}
    grows = lambda genes, vocab: (vocab.update(genes), tuple(sorted(set(genes))))[1]  # noqa: E731
    review = list(REVIEW_RECORDS)
    review_train = C.training_labels_as_of(review, "INH", T)
    if seed is not None:
        rng.shuffle(review)
        rng.shuffle(review_train)

    def invariance(recs, fs, fit):
        return C.invariance_to_test_data(recs, fs, T, fit, antibiotics=("INH",),
                                         constant_keys=getattr(fit, "constant_keys", ()))
    results = [
        C.future_marker(reg, T, sorted(reg.vocabulary(T))),
        C.future_marker(reg, T, ["katG", "embB", "geneC", "fabG1", "madeUp"]),
        C.artifact_manifest(reg, T, reg.manifest(T)),
        C.artifact_manifest(reg, T, reg.manifest("2026-06-01")),
        C.backfill(rows, T, train),
        C.backfill(rows, T, train + [("g03", "INH", "R")]),
        C.backfill(review, T, review_train + [("x2", "INH", "R")]),  # N-1
        C.label_revision(rows, T, train),
        C.label_revision(rows, T, leaked_rev),
        C.label_revision(review, T, review_train + [("x2", "INH", "R")]),  # N-1
        invariance(rows, feats, make_fit()),
        invariance(rows, feats, make_fit(label="latest")),
        invariance(rows, feats, make_fit(vocabulary="all")),
        invariance(rows, feats, make_fit(gate=False)),
        invariance(rows, feats, make_fit(history="all")),
        invariance(rows, feats, make_fit(availability="features_first")),
        invariance(review, REVIEW_FEATURES, make_cached(vocab_all)),  # N-2
        invariance(review, REVIEW_FEATURES, make_dedup("duplicate_group")),  # N-3
        invariance(review, REVIEW_FEATURES, vocab_knowable),  # ruling A
        C.unseen_vocabulary(reg.vocabulary(T), test),
        C.unseen_vocabulary(reg.vocabulary(T), test, grows),
        C.unseen_vocabulary(reg.vocabulary(T), {"x9": ["rpoB", "katG_S315T"]},
                            lambda genes, voc: tuple(sorted({g.split("_")[0] for g in genes} & voc)), "drop"),  # N-4
    ]
    return [result.to_dict() for result in results]


class Determinism(unittest.TestCase):
    def test_every_canary_passes_clean_and_fails_leak(self):
        results = all_results()
        self.assertEqual([r["passed"] for r in results],
                         [True, False] * 2 + [True, False, False] * 2 + [True] + [False] * 8 + [True, False, False])
        for result in results:
            self.assertEqual(bool(result["offenders"]), not result["passed"], result)

    def test_input_order_does_not_matter(self):
        expected = all_results()
        for seed in range(3):
            self.assertEqual(all_results(seed), expected)

    def test_hash_seed_does_not_matter(self):
        # -I would imply -E and drop PYTHONHASHSEED (review N-5), so the child
        # runs with -s -P -B, an explicit environment and a cwd under .tmp/.
        code = ("import importlib.util, json, os, sys; "
                "print(json.dumps([os.environ['PYTHONHASHSEED'], hash('grm')])); "
                "spec = importlib.util.spec_from_file_location('canary_tests', sys.argv[1]); "
                "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); "
                "print(json.dumps(m.all_results(1), sort_keys=True))")
        scratch = ROOT / ".tmp"
        scratch.mkdir(exist_ok=True)
        seen, outputs = {}, set()
        with tempfile.TemporaryDirectory(dir=scratch) as cwd:
            for hash_seed in ("0", "1", "12345"):
                env = {"PYTHONHASHSEED": hash_seed, "PATH": os.environ.get("PATH", "/usr/bin:/bin"), "TMPDIR": cwd}
                done = subprocess.run([sys.executable, "-s", "-P", "-B", "-c", code, __file__], env=env, cwd=cwd,
                                      capture_output=True, text=True, check=True)
                header, body = done.stdout.split("\n", 1)
                reported_seed, grm_hash = json.loads(header)
                self.assertEqual(reported_seed, hash_seed)
                seen[hash_seed] = grm_hash
                outputs.add(body)
        self.assertEqual(len(set(seen.values())), 3, seen)  # every seed took effect
        self.assertEqual(len(outputs), 1)
        self.assertEqual(json.loads(outputs.pop()), json.loads(json.dumps(all_results(1), sort_keys=True)))


if __name__ == "__main__":
    unittest.main()
