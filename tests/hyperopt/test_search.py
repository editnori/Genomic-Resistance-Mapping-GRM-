"""Deterministic search, failure recording, resume and history round-trips."""
import collections
import copy
import hashlib
import json
import os
import tempfile
import unittest
from unittest import mock

from scripts.hyperopt import protocol as hp
from scripts.hyperopt import search

FOLD_PLAN = {"inner_folds": 3}
FOLD_PLAN_HASH = hp.fold_plan_hash(FOLD_PLAN)
ALL_FAMILIES = hp.LEARNERS + ("generic",)
TOP = 1 - 2**-53


def build(**kwargs):
    kwargs.setdefault("families", ALL_FAMILIES)
    kwargs.setdefault("trials_per_family", 6)
    return hp.build_protocol(FOLD_PLAN_HASH, **kwargs)


def rehash(protocol):
    protocol["protocol_hash"] = hp.compute_hash(protocol)
    hp.verify_protocol(protocol)
    return protocol


def objective(family, config, fold_plan, seed):
    """Scores derived only from the inputs, so reruns are exact."""
    key = hp.canonical_json([family, config, seed])
    scores = []
    for fold in reversed(range(fold_plan["inner_folds"])):
        byte = hashlib.sha256(f"{key}|{fold}".encode()).digest()[0]
        scores.append({"fold": fold, "metrics": {"balanced_accuracy": byte / 255}})
    return {"status": "ok", "inner_scores": scores, "fit_count": fold_plan["inner_folds"],
            "seconds": 0.25, "complexity": len(key)}


def is_default(family, config):
    return hp.canonical_json(config) == hp.canonical_json(hp.FAMILIES[family]["default"])


class SearchCase(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name

    def path(self, name="history.jsonl"):
        return os.path.join(self.dir, name)

    def run_search(self, protocol, fn=objective, path=None, fold_plan=FOLD_PLAN, **kwargs):
        return search.run_search(protocol, fn, fold_plan=fold_plan,
                                 history_path=path or self.path(), **kwargs)

    def read_bytes(self, path=None):
        with open(path or self.path(), "rb") as stream:
            return stream.read()

    def rewrite_line(self, index, edit, path=None):
        with open(path or self.path()) as stream:
            lines = stream.readlines()
        record = json.loads(lines[index])
        edit(record)
        lines[index] = hp.history_line(record)
        with open(path or self.path(), "w") as stream:
            stream.writelines(lines)


class DeterminismTests(SearchCase):
    def test_same_protocol_and_seed_write_identical_bytes(self):
        self.run_search(build(), path=self.path("a.jsonl"))
        self.run_search(build(), path=self.path("b.jsonl"))
        first = self.read_bytes(self.path("a.jsonl"))
        self.assertEqual(first, self.read_bytes(self.path("b.jsonl")))
        self.run_search(build(sampler_seed=7), path=self.path("c.jsonl"))
        self.assertNotEqual(first, self.read_bytes(self.path("c.jsonl")))

    def test_trial_zero_is_the_default_for_every_family(self):
        records = self.run_search(build())
        firsts = {r["family"]: r for r in records if r["trial_id"] == 0}
        self.assertEqual(set(firsts), set(ALL_FAMILIES))
        for family, record in firsts.items():
            self.assertTrue(record["is_default"])
            self.assertTrue(is_default(family, record["config"]), family)
        self.assertFalse(any(r["is_default"] for r in records if r["trial_id"]))

    def test_every_family_gets_the_same_budget(self):
        records = self.run_search(build(trials_per_family=5))
        counts = collections.Counter(r["family"] for r in records)
        self.assertEqual(counts, {family: 5 for family in ALL_FAMILIES})

    def test_run_order_is_trial_id_then_family(self):
        plan = search.plan_trials(build(trials_per_family=3))
        self.assertEqual([(t["trial_id"], t["family"]) for t in plan],
                         [(i, f) for i in range(3) for f in sorted(ALL_FAMILIES)])

    def test_sampled_configs_lie_inside_their_spaces(self):
        protocol = build(trials_per_family=300)
        for trial in search.plan_trials(protocol):
            if trial["trial_id"]:
                space = protocol["families"][trial["family"]]["space"]
                self.assertIsNone(hp.config_problem(trial["config"], space), trial)
                self.assertNotIn(hp.INTERNAL_SEARCH, trial["config"])

    def test_sample_value_reaches_both_bounds(self):
        for family in hp.FAMILIES.values():
            for spec in family["space"].values():
                with self.subTest(spec=spec):
                    if spec["type"] == "choice":
                        ends = (spec["values"][0], spec["values"][-1])
                    else:
                        ends = (spec["low"], spec["high"])
                    self.assertEqual(search.sample_value(spec, 0.0), ends[0])
                    self.assertEqual(search.sample_value(spec, TOP), ends[1])

    def test_linear_int_draws_every_value_equally(self):
        spec = {"type": "int", "low": 1, "high": 4, "log": False}
        counts = collections.Counter(search.sample_value(spec, i / 400) for i in range(400))
        self.assertEqual(counts, {1: 100, 2: 100, 3: 100, 4: 100})

    def test_adding_a_family_keeps_other_draws(self):
        small = search.plan_trials(build(families=("cart",)))
        large = [t for t in search.plan_trials(build()) if t["family"] == "cart"]
        self.assertEqual(small, large)

    def test_families_with_one_space_draw_independently(self):
        protocol = build(families=("generic",), trials_per_family=8)
        protocol["families"]["twin"] = copy.deepcopy(protocol["families"]["generic"])
        rehash(protocol)
        configs = collections.defaultdict(list)
        for trial in search.plan_trials(protocol):
            if trial["trial_id"]:
                configs[trial["family"]].append(trial["config"])
        self.assertNotEqual(configs["generic"], configs["twin"])


class ObjectiveCallTests(SearchCase):
    def test_objective_receives_the_fit_seed_the_record_states(self):
        seeds = []

        def spy(family, config, fold_plan, seed):
            seeds.append(seed)
            return objective(family, config, fold_plan, seed)

        records = self.run_search(build(sampler_seed=7, fit_seed=42, trials_per_family=2), spy)
        self.assertEqual(set(seeds), {42})
        self.assertEqual({r["seed"] for r in records}, {42})

    def test_objective_receives_private_copies_of_config_and_plan(self):
        def mutating(family, config, fold_plan, seed):
            result = objective(family, config, fold_plan, seed)
            config.clear()
            fold_plan["inner_folds"] = 99
            return result

        records = self.run_search(build(trials_per_family=2), mutating)
        self.assertTrue(all(r["status"] == "ok" and r["config"] for r in records))
        self.assertEqual(FOLD_PLAN, {"inner_folds": 3})

    def test_every_line_is_fsynced(self):
        with mock.patch("os.fsync", wraps=os.fsync) as fsync:
            records = self.run_search(build(trials_per_family=2))
        self.assertEqual(fsync.call_count, len(records))


class HistoryTests(SearchCase):
    def test_lines_have_exactly_the_documented_fields(self):
        protocol = build()
        records = self.run_search(protocol)
        with open(self.path()) as stream:
            lines = [json.loads(line) for line in stream]
        self.assertEqual(lines, records)
        for record in lines:
            self.assertEqual(set(record), hp.HISTORY_FIELDS)
            self.assertEqual(record["protocol_hash"], protocol["protocol_hash"])
            self.assertEqual(record["fold_plan_hash"], FOLD_PLAN_HASH)
            self.assertEqual([s["fold"] for s in record["inner_scores"]], [0, 1, 2])

    def test_jsonl_round_trips(self):
        protocol = build()
        records = self.run_search(protocol)
        read = hp.read_history(self.path(), protocol, [0, 1, 2])
        self.assertEqual(read, records)
        self.assertEqual("".join(map(hp.history_line, read)).encode(), self.read_bytes())

    def test_torn_last_line_is_refused_until_removed(self):
        self.run_search(build(), path=self.path("full.jsonl"))
        self.run_search(build())
        data = self.read_bytes()
        with open(self.path(), "wb") as stream:
            stream.write(data[:-40])
        with self.assertRaisesRegex(hp.ProtocolError, "delete everything after the last newline"):
            self.run_search(build())
        kept = data[:-40]
        with open(self.path(), "wb") as stream:
            stream.write(kept[:kept.rindex(b"\n") + 1])
        self.run_search(build())
        self.assertEqual(self.read_bytes(), self.read_bytes(self.path("full.jsonl")))

    def test_repeated_keys_and_non_canonical_lines_are_refused(self):
        protocol = build(families=("generic",), trials_per_family=2)
        self.run_search(protocol)
        original = self.read_bytes().decode()
        lines = original.splitlines(keepends=True)
        cases = {
            "repeated JSON keys": lines[1].replace("{", '{"status":"failed",', 1),
            "not canonical": lines[1].replace(",", ", ", 1),
        }
        for message, line in cases.items():
            with self.subTest(message):
                with open(self.path(), "w") as stream:
                    stream.write(lines[0] + line)
                with self.assertRaisesRegex(hp.ProtocolError, message):
                    hp.read_history(self.path(), protocol)

    def test_crlf_line_endings_are_refused(self):
        protocol = build(families=("generic",), trials_per_family=4)
        self.run_search(protocol)
        data = self.read_bytes()
        for name, bad in {"all CRLF": data.replace(b"\n", b"\r\n"),
                          "CRLF prefix": data.replace(b"\n", b"\r\n", 2).rsplit(b"\n", 3)[0]
                          + b"\n",
                          "bare CR": data.replace(b"\n", b"\r", 1)}.items():
            with self.subTest(name):
                with open(self.path(), "wb") as stream:
                    stream.write(bad)
                with self.assertRaisesRegex(hp.ProtocolError, ":1: carriage return"):
                    hp.read_history(self.path(), protocol)
                with self.assertRaisesRegex(hp.ProtocolError, "carriage return"):
                    self.run_search(protocol)
                self.assertEqual(self.read_bytes(), bad)

    def test_huge_integer_in_history_is_a_protocol_error(self):
        protocol = build(families=("generic",), trials_per_family=2)
        self.run_search(protocol)
        self.rewrite_line(1, lambda r: r.update(seconds=10**400))
        with self.assertRaisesRegex(hp.ProtocolError, ":2: seconds must be"):
            hp.read_history(self.path(), protocol)

    def test_non_ascii_history_is_a_protocol_error(self):
        protocol = build(families=("generic",), trials_per_family=2)
        self.run_search(protocol)
        with open(self.path(), "ab") as stream:
            stream.write(b"\xc3\xa9\n")
        with self.assertRaisesRegex(hp.ProtocolError, "not ASCII"):
            hp.read_history(self.path(), protocol)

    def test_repeated_trial_is_refused(self):
        protocol = build(families=("generic",), trials_per_family=2)
        self.run_search(protocol)
        with open(self.path()) as stream:
            first = stream.readline()
        with open(self.path(), "a") as stream:
            stream.write(first)
        with self.assertRaisesRegex(hp.ProtocolError, "repeats"):
            hp.read_history(self.path(), protocol)


class RecordCheckTests(SearchCase):
    def setUp(self):
        super().setUp()
        self.protocol = build(families=("generic",), trials_per_family=3)
        self.records = self.run_search(self.protocol)

    def assertRefused(self, record, message, **kwargs):
        with self.assertRaisesRegex(hp.ProtocolError, message):
            hp.validate_record(record, **kwargs)

    def test_written_records_pass(self):
        for record in self.records:
            hp.validate_record(record, self.protocol, [0, 1, 2])

    def test_is_default_must_match_trial_zero(self):
        for index in (0, 1):
            record = copy.deepcopy(self.records[index])
            record["is_default"] = not record["is_default"]
            self.assertRefused(record, "is_default")

    def test_seed_must_be_the_fit_seed(self):
        record = dict(self.records[1], seed=self.protocol["seeds"]["fit"] + 1)
        self.assertRefused(record, "seeds.fit", protocol=self.protocol)

    def test_trial_id_must_be_within_the_budget(self):
        record = dict(self.records[1], trial_id=self.protocol["trials_per_family"])
        self.assertRefused(record, "trials_per_family", protocol=self.protocol)

    def test_trial_zero_must_carry_the_default_config(self):
        record = dict(self.records[0], config=self.records[1]["config"])
        self.assertRefused(record, "default config", protocol=self.protocol)

    def test_ok_record_must_match_the_plan_folds(self):
        record = self.records[1]
        self.assertRefused(record, "exactly the plan's inner folds", folds=[0, 1, 2, 3])
        one = dict(record, inner_scores=record["inner_scores"][:1])
        self.assertRefused(one, "at least 2", folds=[0])

    def test_internal_search_config_is_only_for_trial_zero(self):
        protocol = build(families=("cart",), trials_per_family=2)
        records = self.run_search(protocol, path=self.path("cart.jsonl"))
        record = dict(records[1], config=hp.FAMILIES["cart"]["default"])
        self.assertRefused(record, "config", protocol=protocol)


class FoldPlanTests(SearchCase):
    def test_hash_matches_task_d_bytes(self):
        # Task D's canaries.fold_plan_hash({"inner_folds": 3}), recorded by the review.
        self.assertEqual(FOLD_PLAN_HASH,
                         "cff93e3b784c576c9416dce970165c96b2b0dd2010dd49835aa5fa32642d4183")
        self.assertEqual(hp.fold_plan_hash({"b": [1, 2.5, None, True], "a": "x"}),
                         hashlib.sha256(b'{"a":"x","b":[1,2.5,null,true]}').hexdigest())

    def test_plan_must_be_plain_json(self):
        for plan in ({1: "x"}, {"inner_folds": {0, 1}}, {"x": float("nan")}, {"x": object()}):
            with self.subTest(plan=plan):
                with self.assertRaisesRegex(hp.ProtocolError, "plain JSON"):
                    hp.fold_plan_hash(plan)

    def test_inner_fold_ids(self):
        self.assertEqual(hp.inner_fold_ids({"inner_folds": 3}), [0, 1, 2])
        self.assertEqual(hp.inner_fold_ids({"inner_folds": [2, 0, None, 1, 0]}), [0, 1, 2])
        for bad in (1, True, [0, 0, None], [0, -1], [0, "1"], "3", None, [True, False, None],
                    [0, 1, True]):
            with self.subTest(bad=bad):
                with self.assertRaises(hp.ProtocolError):
                    hp.inner_fold_ids({"inner_folds": bad})
        with self.assertRaises(hp.ProtocolError):
            hp.inner_fold_ids([0, 1])

    def test_plan_other_than_the_hashed_one_is_refused(self):
        other = {"inner_folds": 3, "outer_fold": 99}
        with self.assertRaisesRegex(hp.ProtocolError, "fold plan hashes to"):
            self.run_search(build(), fold_plan=other)
        self.assertFalse(os.path.exists(self.path()))
        trial = search.plan_trials(build())[0]
        with self.assertRaisesRegex(hp.ProtocolError, "fold plan hashes to"):
            search.run_trial(build(), objective, trial, other)

    def test_supplied_hash_must_also_match(self):
        with self.assertRaisesRegex(hp.ProtocolError, "fold plan hash"):
            self.run_search(build(), fold_plan_hash="e" * 64)
        self.assertFalse(os.path.exists(self.path()))
        self.run_search(build(trials_per_family=1), fold_plan_hash=FOLD_PLAN_HASH)

    def test_plan_with_one_inner_fold_is_refused_before_running(self):
        plan = {"inner_folds": 1}
        protocol = hp.build_protocol(hp.fold_plan_hash(plan), families=("generic",))
        with self.assertRaisesRegex(hp.ProtocolError, "at least 2"):
            self.run_search(protocol, fold_plan=plan)
        self.assertFalse(os.path.exists(self.path()))

    def test_per_sample_label_plan_runs(self):
        plan = {"inner_folds": [0, 1, None, 1, 0, None]}
        protocol = hp.build_protocol(hp.fold_plan_hash(plan), families=("generic",),
                                     trials_per_family=2)

        def labelled(family, config, fold_plan, seed):
            folds = sorted({f for f in fold_plan["inner_folds"] if f is not None})
            return objective(family, config, {"inner_folds": len(folds)}, seed)

        records = self.run_search(protocol, labelled, fold_plan=plan)
        self.assertEqual([r["status"] for r in records], ["ok", "ok"])


class FailureTests(SearchCase):
    def test_failed_and_unconverged_trials_are_recorded_with_reasons(self):
        def flaky(family, config, fold_plan, seed):
            result = objective(family, config, fold_plan, seed)
            trial = (family, is_default(family, config))
            if trial == ("scm", True):
                raise RuntimeError("solver blew up")
            if trial == ("cart", True):
                return dict(result, status="unconverged", reason="path hit iteration cap",
                            inner_scores=result["inner_scores"][:1])
            if trial == ("rf", True):
                return dict(result, status="failed", reason="no resistant samples",
                            complexity=None, inner_scores=[])
            if trial == ("gbt", True):
                result["inner_scores"][0]["metrics"]["balanced_accuracy"] = float("nan")
            if trial == ("lasso", True):
                result["inner_scores"][1]["metrics"] = {"auroc": 0.5}
            if trial == ("generic", True):
                return dict(result, reason="ok trials carry no reason")
            return result

        protocol = build(trials_per_family=2)
        records = self.run_search(protocol, flaky)
        self.assertEqual(len(records), 12)
        defaults = {r["family"]: r for r in records if r["is_default"]}
        expected = {
            "scm": ("failed", "objective raised RuntimeError: solver blew up"),
            "cart": ("unconverged", "path hit iteration cap"),
            "rf": ("failed", "no resistant samples"),
            "gbt": ("failed", "invalid objective result: inner_scores: fold 2 metrics"),
            "lasso": ("failed", "invalid objective result: fold 1 is missing metrics"),
            "generic": ("failed", "invalid objective result: reason must be present"),
        }
        for family, (status, reason) in expected.items():
            with self.subTest(family):
                record = defaults[family]
                self.assertEqual(record["status"], status)
                self.assertTrue(record["reason"].startswith(reason), record["reason"])
                self.assertEqual(set(record), hp.HISTORY_FIELDS | {"reason"})
        self.assertEqual(len(defaults["cart"]["inner_scores"]), 1)
        self.assertIsNone(defaults["scm"]["complexity"])
        self.assertTrue(all(r["status"] == "ok" for r in records if not r["is_default"]))
        self.assertEqual(hp.read_history(self.path(), protocol, [0, 1, 2]), records)

    def test_results_off_the_plan_folds_become_failed_trials(self):
        def ragged(family, config, fold_plan, seed):
            result = objective(family, config, fold_plan, seed)
            if family == "rf" and is_default(family, config):
                return dict(result, inner_scores=result["inner_scores"][1:])
            if family == "gbt" and is_default(family, config):
                return dict(result, inner_scores=result["inner_scores"][:1])
            if family == "lasso" and is_default(family, config):
                extra = {"fold": 7, "metrics": {"balanced_accuracy": 0.5}}
                return dict(result, status="unconverged", reason="cap",
                            inner_scores=result["inner_scores"][:1] + [extra])
            return result

        protocol = build(trials_per_family=2)
        records = self.run_search(protocol, ragged)
        defaults = {r["family"]: r for r in records if r["is_default"]}
        expected = {
            "rf": "exactly the plan's inner folds [0, 1, 2], got [0, 1]",
            "gbt": "at least 2 inner folds",
            "lasso": "folds [7] are not in the plan's inner folds",
        }
        for family, reason in expected.items():
            with self.subTest(family):
                record = defaults[family]
                self.assertEqual(record["status"], "failed")
                self.assertTrue(record["reason"].startswith("invalid objective result: "))
                self.assertIn(reason, record["reason"])
                self.assertEqual(record["inner_scores"], [])
        ok_folds = {tuple(s["fold"] for s in r["inner_scores"])
                    for r in records if r["status"] == "ok"}
        self.assertEqual(ok_folds, {(0, 1, 2)})
        self.assertEqual(hp.read_history(self.path(), protocol, [0, 1, 2]), records)

    def test_malformed_results_become_failed_trials(self):
        cases = {
            "not an object": lambda r: None,
            "missing fields": lambda r: {"status": "ok"},
            "unknown fields": lambda r: dict(r, extra=1),
            "unknown status": lambda r: dict(r, status="done", reason="x"),
            "empty ok scores": lambda r: dict(r, inner_scores=[]),
            "repeated fold": lambda r: dict(r, inner_scores=r["inner_scores"][:1] * 2),
            "negative seconds": lambda r: dict(r, seconds=-1.0),
            "boolean fit count": lambda r: dict(r, fit_count=True),
            "null complexity on ok": lambda r: dict(r, complexity=None),
            "integer metric key": lambda r: dict(r, inner_scores=[
                {"fold": e["fold"], "metrics": {**e["metrics"], 1: 0.5}}
                for e in r["inner_scores"]]),
        }
        cases.update(self.HUGE)
        protocol = build(families=("generic",))
        trial = search.plan_trials(protocol)[0]
        for name, mutate in cases.items():
            with self.subTest(name):
                record = search.run_trial(
                    protocol, lambda *args: mutate(objective(*args)), trial, FOLD_PLAN)
                self.assertEqual(record["status"], "failed")
                self.assertTrue(record["reason"].startswith("invalid objective result: "))
                hp.validate_record(record, protocol, [0, 1, 2])

    # Integers too large for a float: refused as invalid, not raised as OverflowError.
    HUGE = {
        "huge metric": lambda r: dict(r, inner_scores=[
            {"fold": e["fold"], "metrics": {"balanced_accuracy": 10**400}}
            for e in r["inner_scores"]]),
        "huge seconds": lambda r: dict(r, seconds=10**400),
        "huge complexity": lambda r: dict(r, complexity=10**400),
    }

    def test_huge_integers_become_failed_trials_with_reasons(self):
        expected = {"huge metric": "are not finite numbers",
                    "huge seconds": "seconds must be a non-negative number",
                    "huge complexity": "complexity must be a non-negative number"}
        protocol = build(families=("generic",), trials_per_family=2)

        def huge(family, config, fold_plan, seed):
            result = objective(family, config, fold_plan, seed)
            return self.HUGE[name](result) if is_default(family, config) else result

        for name, reason in expected.items():
            with self.subTest(name):
                path = self.path(f"{name}.jsonl")
                records = self.run_search(protocol, huge, path=path)
                self.assertEqual([r["status"] for r in records], ["failed", "ok"])
                self.assertIn(reason, records[0]["reason"])
                self.assertEqual(hp.read_history(path, protocol, [0, 1, 2]), records)


class RefusalTests(SearchCase):
    def test_mutated_protocol_is_refused_before_running(self):
        protocol = build()
        protocol["families"]["cart"]["space"]["max_depth"]["high"] = 12
        with self.assertRaisesRegex(hp.ProtocolError, "hash mismatch"):
            self.run_search(protocol)

    def test_resume_appends_the_remaining_trials_identically(self):
        self.run_search(build(), path=self.path("full.jsonl"))
        calls = []

        def interrupted(*args):
            if len(calls) == 10:
                raise KeyboardInterrupt
            calls.append(args)
            return objective(*args)

        with self.assertRaises(KeyboardInterrupt):
            self.run_search(build(), interrupted)
        self.assertEqual(len(hp.read_history(self.path())), 10)
        self.run_search(build())
        self.assertEqual(self.read_bytes(), self.read_bytes(self.path("full.jsonl")))

    def test_resumed_history_with_mixed_protocol_hashes_is_refused(self):
        first, second = build(), build(fit_seed=7)
        self.run_search(first)
        self.run_search(second, path=self.path("other.jsonl"))
        with open(self.path("other.jsonl")) as src, open(self.path(), "a") as dst:
            dst.write(src.readline())
        with self.assertRaisesRegex(hp.ProtocolError, "mixes protocol hashes"):
            self.run_search(first)

    def test_history_from_another_protocol_is_refused(self):
        self.run_search(build(fit_seed=7))
        with self.assertRaisesRegex(hp.ProtocolError, "protocol_hash"):
            self.run_search(build())

    def test_history_out_of_plan_order_is_refused(self):
        protocol = build()
        self.run_search(protocol)
        with open(self.path()) as stream:
            lines = stream.readlines()
        with open(self.path(), "w") as stream:
            stream.writelines(lines[1:])
        with self.assertRaisesRegex(hp.ProtocolError, "not planned trial"):
            self.run_search(protocol)

    def test_resume_rechecks_folds_against_the_plan(self):
        protocol = build(families=("generic",), trials_per_family=3)
        self.run_search(protocol)
        self.rewrite_line(1, lambda r: r.update(inner_scores=r["inner_scores"][1:]))
        hp.read_history(self.path(), protocol)
        with self.assertRaisesRegex(hp.ProtocolError,
                                    ":2: inner_scores: an ok trial must report exactly"):
            self.run_search(protocol)

    def test_resume_refuses_an_in_space_but_unplanned_config(self):
        protocol = build(families=("generic",), trials_per_family=3)
        records = self.run_search(protocol)
        self.rewrite_line(1, lambda r: r.update(config=records[2]["config"]))
        hp.read_history(self.path(), protocol)
        with self.assertRaisesRegex(hp.ProtocolError, "not planned trial generic:1"):
            self.run_search(protocol)

    def test_resume_compares_configs_by_canonical_json(self):
        protocol = build(families=("generic",), trials_per_family=8)
        protocol["families"]["generic"] = {
            "default": {"flag": True},
            "space": {"flag": {"type": "choice", "values": [True, 1]}}}
        rehash(protocol)
        plan = search.plan_trials(protocol)
        index = next(i for i, t in enumerate(plan)
                     if t["trial_id"] and t["config"]["flag"] is not True)
        self.run_search(protocol)
        self.rewrite_line(index, lambda r: r.update(config={"flag": True}))
        with self.assertRaisesRegex(hp.ProtocolError, "not planned trial"):
            self.run_search(protocol)
        self.rewrite_line(index, lambda r: r.update(config={"flag": 1}))
        self.rewrite_line(0, lambda r: r.update(config={"flag": 1}))
        with self.assertRaisesRegex(hp.ProtocolError, "default config"):
            self.run_search(protocol)


if __name__ == "__main__":
    unittest.main()
