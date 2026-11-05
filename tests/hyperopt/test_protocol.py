"""Protocol hashing, freezing and validation."""
import copy
import json
import os
import re
import tempfile
import unittest
from pathlib import Path

from scripts.hyperopt import protocol as hp

FOLD_PLAN = "f" * 64
RUST = Path(__file__).resolve().parents[2] / "src" / "learn"


def build(**kwargs):
    return hp.build_protocol(FOLD_PLAN, **kwargs)


class ProtocolHashTests(unittest.TestCase):
    def test_hash_is_stable_and_matches_content(self):
        first, second = build(), build()
        self.assertEqual(first["protocol_hash"], second["protocol_hash"])
        self.assertEqual(first["protocol_hash"], hp.compute_hash(first))
        self.assertRegex(first["protocol_hash"], "^[0-9a-f]{64}$")

    def test_default_protocol_hash_is_pinned(self):
        # Changes when any default, space, budget, seed or objective does. The
        # review's F1 dropped lasso lambda_min_ratio, which moved it from
        # 47edf88bfb507184812b28fb9247f7d17663576a94e2664b670a3c7edf526a39.
        self.assertEqual(build()["protocol_hash"],
                         "e5c93376219cbe5feea30deb282e4104623b753395a6d9f25a21cf4df9ca8a4b")

    def test_hash_ignores_key_order(self):
        protocol = build()
        reordered = {key: protocol[key] for key in reversed(list(protocol))}
        self.assertEqual(hp.compute_hash(reordered), protocol["protocol_hash"])

    def test_every_field_change_changes_the_hash(self):
        base = build(constraints=[{"metric": "complexity", "op": "<=", "bound": 50}])
        edits = {
            "schema_version": lambda p: p.update(schema_version=2),
            "fold_plan_hash": lambda p: p.update(fold_plan_hash="e" * 64),
            "objective.name": lambda p: p["objective"].update(name="auroc"),
            "objective.direction": lambda p: p["objective"].update(direction="minimize"),
            "constraint bound": lambda p: p["objective"]["constraints"][0].update(bound=51),
            "sampler seed": lambda p: p["seeds"].update(sampler=1),
            "fit seed": lambda p: p["seeds"].update(fit=1),
            "budget": lambda p: p.update(trials_per_family=33),
            "default": lambda p: p["families"]["rf"]["default"].update(max_depth=29),
            "internal grid": lambda p: p["families"]["scm"]["default"]["p"].pop(),
            "space bound": lambda p: p["families"]["cart"]["space"]["max_depth"].update(high=31),
            "space scale": lambda p: p["families"]["gbt"]["space"]["rounds"].update(log=False),
            "choice values": lambda p: p["families"]["scm"]["space"]["kind"]["values"].pop(),
            "float bound": lambda p: p["families"]["rf"]["space"]["subsample"].update(
                low=0.30000000000000004),
            "family removed": lambda p: p["families"].pop("lasso"),
        }
        seen = {base["protocol_hash"]}
        for name, edit in edits.items():
            with self.subTest(name):
                changed = copy.deepcopy(base)
                edit(changed)
                digest = hp.compute_hash(changed)
                self.assertNotIn(digest, seen)
                seen.add(digest)
                with self.assertRaises(hp.ProtocolError):
                    hp.verify_protocol(changed)

    def test_canonical_json_is_compact_sorted_and_finite(self):
        self.assertEqual(hp.canonical_json({"b": 0.1, "a": [1e-4, 1, 2.5e20]}),
                         '{"a":[0.0001,1,2.5e+20],"b":0.1}')
        for value in (float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                hp.canonical_json({"x": value})


class ProtocolFileTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = os.path.join(tmp.name, "protocol.json")

    def test_save_load_round_trip(self):
        protocol = build(families=hp.LEARNERS + ("generic",))
        hp.save_protocol(protocol, self.path)
        self.assertEqual(hp.load_protocol(self.path), protocol)

    def test_tampered_space_is_refused_on_load(self):
        hp.save_protocol(build(), self.path)
        with open(self.path) as stream:
            data = json.load(stream)
        data["families"]["cart"]["space"]["max_depth"]["high"] = 40
        with open(self.path, "w") as stream:
            json.dump(data, stream)
        with self.assertRaisesRegex(hp.ProtocolError, "hash mismatch"):
            hp.load_protocol(self.path)

    def test_rehashed_invalid_protocol_is_refused(self):
        protocol = build()
        protocol["families"]["cart"]["default"]["max_depth"] = [99]
        protocol["protocol_hash"] = hp.compute_hash(protocol)
        with open(self.path, "w") as stream:
            json.dump(protocol, stream)
        with self.assertRaisesRegex(hp.ProtocolError, "cart default"):
            hp.load_protocol(self.path)

    def test_repeated_json_keys_are_refused_on_load(self):
        hp.save_protocol(build(), self.path)
        with open(self.path) as stream:
            text = stream.read()
        with open(self.path, "w") as stream:
            stream.write(text.replace("{", '{"schema_version": 1,', 1))
        with self.assertRaisesRegex(hp.ProtocolError, "repeated JSON keys"):
            hp.load_protocol(self.path)

    def test_non_finite_json_is_refused(self):
        with open(self.path, "w") as stream:
            stream.write('{"schema_version": NaN}')
        with self.assertRaises(hp.ProtocolError):
            hp.load_protocol(self.path)


class ProtocolValidationTests(unittest.TestCase):
    def test_every_default_is_valid_for_its_space(self):
        for name, family in hp.FAMILIES.items():
            with self.subTest(name):
                self.assertIsNone(hp.default_problem(family["default"], family["space"]))

    def test_defaults_are_the_production_defaults(self):
        self.assertEqual(hp.FAMILIES["scm"]["default"], {
            "internal_search": True, "kind": ["conjunction", "disjunction"],
            "max_rules": [1, 2, 3, 4, 5, 6],
            "p": [0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0]})
        self.assertEqual(hp.FAMILIES["cart"]["default"], {
            "internal_search": True, "max_depth": list(range(1, 9)), "min_leaf": [1]})
        for name in ("rf", "gbt", "lasso", "generic"):
            with self.subTest(name):
                family = hp.FAMILIES[name]
                self.assertFalse(hp.is_internal_search(family["default"]))
                self.assertIsNone(hp.config_problem(family["default"], family["space"]))

    def test_lasso_searches_only_what_production_exposes(self):
        # fit_path_with_status takes n_lambdas; the smallest-penalty ratio is
        # fixed at 1e-3 inside lasso.rs (review F1), so it is not searched.
        self.assertEqual(hp.FAMILIES["lasso"], {
            "default": {"n_lambdas": 30},
            "space": {"n_lambdas": {"type": "int", "low": 10, "high": 100, "log": False}}})

    def test_bad_internal_search_defaults_are_refused(self):
        space = hp.FAMILIES["cart"]["space"]
        good = hp.FAMILIES["cart"]["default"]
        cases = {
            "differ": {"internal_search": True, "max_depth": [1, 2]},
            "must list": dict(good, min_leaf=1),
            "must list ": dict(good, min_leaf=[]),
            "repeats": dict(good, max_depth=[1, 1]),
            "outside": dict(good, max_depth=[0, 1]),
            "parameters": dict(good, internal_search=False),
            "parameters ": dict(good, internal_search=1),
        }
        for expected, default in cases.items():
            with self.subTest(expected):
                self.assertIn(expected.strip(), hp.default_problem(default, space))
        self.assertFalse(hp.is_internal_search(dict(good, internal_search=1)))

    def test_config_problem_names_the_offending_parameter(self):
        space = hp.FAMILIES["generic"]["space"]
        cases = {
            "n=9": {"mode": "a", "n": 9, "x": 0.0},
            "x=0": {"mode": "a", "n": 4, "x": 0},
            "mode='d'": {"mode": "d", "n": 4, "x": 0.0},
            "n=True": {"mode": "a", "n": True, "x": 0.0},
            "differ": {"mode": "a", "n": 4},
        }
        for expected, config in cases.items():
            with self.subTest(expected):
                self.assertIn(expected, hp.config_problem(config, space))

    def test_integers_too_large_for_a_float_are_refused(self):
        with self.assertRaisesRegex(hp.ProtocolError, "constraint"):
            build(constraints=[{"metric": "complexity", "op": "<=", "bound": 10**400}])
        protocol = build()
        protocol["objective"]["constraints"] = [
            {"metric": "complexity", "op": "<=", "bound": 10**400}]
        protocol["protocol_hash"] = hp.compute_hash(protocol)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "protocol.json")
            with open(path, "w") as stream:
                json.dump(protocol, stream)
            with self.assertRaisesRegex(hp.ProtocolError, "constraint"):
                hp.load_protocol(path)

    def test_parameter_names_must_be_non_empty_and_not_reserved(self):
        cases = {"": {"type": "int", "low": 1, "high": 2, "log": False},
                 "internal_search": {"type": "choice", "values": [False, True]}}
        for param, spec in cases.items():
            with self.subTest(param=param):
                protocol = build()
                protocol["families"]["rf"]["space"][param] = spec
                protocol["families"]["rf"]["default"][param] = False if param else 1
                protocol["protocol_hash"] = hp.compute_hash(protocol)
                with self.assertRaisesRegex(hp.ProtocolError, f"parameter '{param}' must be"):
                    hp.verify_protocol(protocol)

    def test_malformed_protocols_are_refused(self):
        cases = {
            "direction": lambda p: p["objective"].update(direction="up"),
            "constraint op": lambda p: p["objective"].update(
                constraints=[{"metric": "complexity", "op": "<", "bound": 1}]),
            "empty fold plan": lambda p: p.update(fold_plan_hash=""),
            "zero budget": lambda p: p.update(trials_per_family=0),
            "int bound on float": lambda p: p["families"]["rf"]["space"]["subsample"].update(
                low=0),
            "log from zero": lambda p: p["families"]["cart"]["space"]["max_depth"].update(
                low=0, log=True),
            "extra field": lambda p: p.update(note="x"),
            "schema_version true": lambda p: p.update(schema_version=True),
            "schema_version 1.0": lambda p: p.update(schema_version=1.0),
            "empty family name": lambda p: p["families"].update({"": p["families"].pop("rf")}),
            "fixed default outside space": lambda p: p["families"]["rf"]["default"].update(
                max_depth=31),
            "huge constraint bound": lambda p: p["objective"].update(
                constraints=[{"metric": "complexity", "op": "<=", "bound": 10**400}]),
            "huge choice value": lambda p: p["families"]["scm"]["space"]["kind"][
                "values"].append(10**400),
            "huge int default": lambda p: p["families"]["rf"]["default"].update(
                max_depth=10**400),
        }
        for name, edit in cases.items():
            with self.subTest(name):
                protocol = build()
                edit(protocol)
                protocol["protocol_hash"] = hp.compute_hash(protocol)
                with self.assertRaises(hp.ProtocolError):
                    hp.verify_protocol(protocol)


def rust_literal(text):
    """A Rust numeric literal as the JSON number the protocol writes."""
    return float(text) if re.fullmatch(r"\d+\.\d+", text) else int(text)


def rust_default(file, name):
    """The `field: value` lines of `impl Default for <name>` in src/learn."""
    text = (RUST / file).read_text()
    block = re.search(rf"impl Default for {name} \{{\s*fn default\(\) -> Self \{{"
                      rf"\s*{name} \{{\n(.*?)\n\s*\}}\n", text, re.S)
    if block is None:
        raise AssertionError(f"no `impl Default for {name}` in src/learn/{file}")
    return dict(re.findall(r"^\s*(\w+): (.+),$", block.group(1), re.M))


def rust_item(text, start):
    """The Rust item that begins at `start`: through the `}` that closes its
    body, or through its `;` if it has no body (a `const`). Brackets are
    counted, so a `;` inside `[f64; 5]` doesn't end it."""
    i = text.index(start)
    depth, body = 0, False
    for j in range(i, len(text)):
        c = text[j]
        if c in "([{":
            depth += 1
            body = body or (c == "{" and depth == 1)
        elif c in ")]}":
            depth -= 1
            if depth == 0 and c == "}" and body:
                return text[i:j + 1]
        elif c == ";" and depth == 0 and not body:
            return text[i:j + 1]
    raise AssertionError(f"the item at {start!r} doesn't end")


def rust_after(file, start, pattern):
    """The first match of `pattern` inside the item that begins at `start` in
    src/learn/<file>. The search stops at the item's end, so a default that
    stops being a literal fails here rather than matching a later item."""
    text = (RUST / file).read_text()
    match = re.search(pattern, rust_item(text, start))
    if match is None:
        raise AssertionError(f"{pattern!r} not found in the item at {start!r} in src/learn/{file}")
    return match.group(1)


class RustDefaultTests(unittest.TestCase):
    """Review F7: trial 0 is production's default, so each default is read
    back from the Rust source it copies, and a changed Rust default fails
    here instead of drifting from the protocol."""

    def assertSame(self, got, want):
        # Canonical JSON, so an int and a float of equal value differ.
        self.assertEqual(hp.canonical_json(got), hp.canonical_json(want))

    def test_fixed_defaults_are_the_rust_default_impls(self):
        for family, file, name in (("rf", "tree.rs", "ForestParams"),
                                   ("gbt", "gbt.rs", "GbtParams")):
            with self.subTest(family):
                fields = rust_default(file, name)
                default = hp.FAMILIES[family]["default"]
                self.assertSame({k: rust_literal(fields[k]) for k in default}, default)
        config = rust_default("eval.rs", "EvalConfig")
        self.assertEqual((config["forest"], config["gbt"]),
                         ("ForestParams::default()", "GbtParams::default()"))
        self.assertSame(hp.FAMILIES["lasso"]["default"],
                        {"n_lambdas": rust_literal(config["lasso_lambdas"])})
        # The path's smallest-penalty ratio is fixed, so it is not searched (F1).
        self.assertEqual(rust_after("lasso.rs", "fn path_problem(",
                                    r"let ratio: f64 = ([^;]+);"), "1e-3")

    def test_internal_search_defaults_are_the_rust_grids(self):
        config = rust_default("eval.rs", "EvalConfig")
        scm, cart = hp.FAMILIES["scm"]["default"], hp.FAMILIES["cart"]["default"]
        grid = rust_after("scm.rs", "pub const P_GRID", r"= \[([^\]]*)\];")
        self.assertSame(scm["p"], [rust_literal(v.strip()) for v in grid.split(",")])
        kinds = rust_after("scm.rs", "pub fn select(", r"for kind in \[([^\]]*)\]")
        self.assertEqual(scm["kind"], [k.strip().removeprefix("Kind::").lower()
                                       for k in kinds.split(",")])
        self.assertEqual(scm["max_rules"],
                         list(range(1, rust_literal(config["scm_max_rules"]) + 1)))
        self.assertEqual(cart["max_depth"],
                         list(range(1, rust_literal(config["cart_max_depth"]) + 1)))
        for file, start in (("eval.rs", "fn fit_cart("), ("tree.rs", "pub fn select_depth(")):
            with self.subTest(start):
                self.assertSame(cart["min_leaf"],
                                [rust_literal(rust_after(file, start, r"min_leaf: (\d+),"))])


if __name__ == "__main__":
    unittest.main()
