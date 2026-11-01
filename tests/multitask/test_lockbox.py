"""lockbox.py on hand-built synthetic fixtures: definition, identifiers, seal, guards and access log.

usage: python3 -I -B -m unittest discover -s tests/multitask -t .
"""
import copy
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import stat
import sys
import tempfile
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

RECENT_AFTER = "2020-06-30"
NOW = dt.datetime(2026, 10, 30, 12, 0, tzinfo=dt.timezone.utc)
UNWRITTEN_LOG = "/synthetic/never-written/access.jsonl"
# Rule 6 v2 against lineage L9: spellings stripped to L9, and spellings refused with the refusal's text.
CLAIMED = ("L9 ", " L9", "L9\t", "\tL9")
REFUSED = (("l9", "alias collision: 'L9' and 'l9'"), ("\U0000ff2c9", "U\\+FF2C"), ("L9\U0000200b", "U\\+200B"),
           ("L\t9", "U\\+0009"), ("L9\U0000034f", "U\\+034F"), ("L9\U0000fe0f", "U\\+FE0F"),
           ("L9\U000e0100", "U\\+E0100"), ("\U00003164L9", "U\\+3164"), ("L9\U0000115f", "U\\+115F"),
           ("\U000004c09", "U\\+04C0"), ("L9\x7f", "U\\+007F"))
# Every character Python's str.strip() removes, pinned literally (the re-check's AGR-3).
STRIP_SET = ("\t\n\x0b\x0c\r\x1c\x1d\x1e\x1f \x85\xa0          "
             "      　")


def genome(gid, lineage, dup, date, tasks):
    return {"genome_id": gid, "lineage": lineage, "duplicate_group": dup, "collection_date": date, "tasks": tasks}


def fixture():
    return [
        genome("a1", "L1", None, "2010-05-01", {"INH": 1, "RIF": 0}),
        genome("a2", "L1", "D1", "2012", {"INH": 0, "RIF": 0}),
        genome("b1", "L2", "D2", "2021-03-04", {"INH": 1, "RIF": 1}),   # recent
        genome("b2", "L2", None, "2014-06", {"INH": 0}),
        genome("c1", "L3", "D2", "2011-01-01", {"INH": 0, "RIF": 1}),   # duplicate of b1
        genome("c2", "L3", None, "2020", {"INH": 1}),                   # may be after 2020-06-30
        genome("c3", "L3", None, "2020-06", {"RIF": 0}),                # ends on the cutoff: not recent
        genome("d1", "L4", None, "2016-02-02", {"INH": 1}),             # held-out lineage
        genome("d2", "L4", "D1", None, {"RIF": 1}),                     # held-out lineage, pulls in a2
        genome("e1", "L5", None, None, {"INH": 0}),
    ]


def define(genomes=None, **overrides):
    params = dict(recent_after=RECENT_AFTER, access_log=UNWRITTEN_LOG, holdout_lineages=["L4", "L9"], tasks=["EMB"],
                  min_per_task=5, min_class=2)
    params.update(overrides)
    return lb.define_lockbox(fixture() if genomes is None else genomes, **params)


def reviewer_base():
    """The independent review's probe_leak.py fixture: 40 old genomes, recent r1, held-out h1."""
    def g(gid, lineage, dup, date, tasks=None):
        return genome(gid, lineage, dup, date, tasks or {"INH": 0})
    base = [g("a%02d" % i, "L%d" % (i % 4), None, "2010-01-01", {"INH": i % 2}) for i in range(40)]
    return base + [g("r1", "L1", None, "2022-01-01"), g("h1", "L9", None, "2011-01-01")], g


def release0():
    """The re-check's rc_p_probe_deep.py fixture: recent r1 and its old duplicate a in D1, held-out h,
    and a development duplicate group G = {d1, d2}."""
    filler = [genome("f%02d" % i, "L%d" % (i % 4), None, "2010-01-01", {"INH": i % 2}) for i in range(30)]
    return filler + [genome("r1", "L1", "D1", "2022-01-01", {"INH": 0}), genome("a", "L2", "D1", "2010-01-01", {"INH": 0}),
                     genome("h", "L9", None, "2011-01-01", {"INH": 0}), genome("d1", "L3", "G", "2010-01-01", {"INH": 0}),
                     genome("d2", "L3", "G", "2010-01-01", {"INH": 0})]


def release(base, *extra, drop=(), regroup=None):
    out = [copy.deepcopy(x) for x in base if x["genome_id"] not in drop]
    for x in out:
        x["duplicate_group"] = (regroup or {}).get(x["genome_id"], x["duplicate_group"])
    return out + [genome(gid, lineage, dup, date, {"INH": 0}) for gid, lineage, dup, date in extra]


def resealed(box, **changes):
    forged = copy.deepcopy(box)
    forged.pop("seal")
    forged.update(changes)
    return dict(forged, seal=lb.sha256_json(forged))


def near_miss(seal):
    """The seal with only its last hex digit changed."""
    return seal[:-1] + ("0" if seal[-1] != "0" else "1")


class Definition(unittest.TestCase):
    def test_union_of_recent_and_holdout_lineages_closed_over_duplicates(self):
        box = define()
        self.assertEqual(box["genome_ids"], ["a2", "b1", "c1", "c2", "d1", "d2"])
        self.assertEqual(box["duplicate_groups"], ["D1", "D2"])
        self.assertEqual(box["definition_groups"], {"D1": ["a2", "d2"], "D2": ["b1", "c1"]})
        self.assertEqual(box["definition"]["recent_after"], RECENT_AFTER)
        self.assertEqual(box["definition"]["holdout_lineages"], ["L4", "L9"])
        self.assertEqual(box["access_log"], UNWRITTEN_LOG)
        self.assertEqual(box["counts"], {
            "input_genomes": 10, "lockbox_genomes": 6, "development_genomes": 4, "recent": 2, "holdout_lineage": 2,
            "duplicate_closure": 2, "undated_input": 2, "holdout_lineages_absent": ["L9"]})

    def test_insufficient_evidence_per_task(self):
        tasks = define()["tasks"]
        self.assertEqual(tasks["INH"], {"n": 5, "classes": {"0": 2, "1": 3}, "reasons": [], "status": "ok"})
        self.assertEqual(tasks["RIF"]["status"], "insufficient_evidence")
        self.assertEqual(tasks["RIF"]["reasons"], ["lockbox has 4 labelled genomes, below 5",
                                                   "lockbox has 1 of class 0, below 2"])
        self.assertEqual(tasks["EMB"], {"n": 0, "classes": None, "status": "insufficient_evidence",
                                        "reasons": ["lockbox has 0 labelled genomes, below 5"]})

    def test_rejects_bad_input(self):
        rows = fixture()
        cases = [
            (rows + [genome("a1", "L1", None, None, {})], "duplicate genome_id"),
            (rows + [genome("z1", "L1", None, None, ["INH"])], "mixes"),
            ([genome("z1", "L1", None, "2020-13", {})], "month"),
            ([genome("z1", "L1", None, "20/01/2020", {})], "collection_date"),
            ([genome("z1", "L1", None, None, {"INH": 2})], "0 or 1"),
            ([genome("z1", None, None, None, {})], "lineage"),
            ([genome("z1", "L1", "", None, {})], "duplicate_group"),
            ([genome("z1", "L1", None, None, "INH")], "string"),
        ]
        for genomes, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                define(genomes)
        with self.assertRaisesRegex(ValueError, "YYYY-MM-DD"):
            define(recent_after="2020")
        with self.assertRaisesRegex(ValueError, "list of lineage"):
            define(holdout_lineages="L4")

    def test_access_log_must_be_an_absolute_path(self):
        # The re-check's L4: a relative log made the bound chain depend on the working directory.
        for access_log in ("", "access.jsonl", "./access.jsonl", "logs/access.jsonl", " /abs/access.jsonl",
                           "/abs/access.jsonl\n", None):
            with self.subTest(access_log=access_log), self.assertRaisesRegex(ValueError, "access_log must be an absolute path"):
                define(access_log=access_log)

    def test_minimums_below_one_are_refused(self):
        # The re-check's L5: min_per_task=0, min_class=0 marked a lockbox task with no positives ok.
        for name in ("min_per_task", "min_class"):
            for value in (0, -1, True, 1.5, "2", None):
                with self.subTest(name=name, value=value), \
                        self.assertRaisesRegex(ValueError, name + " must be an integer of at least 1"):
                    define(**{name: value})
        self.assertEqual(define(min_per_task=1, min_class=1)["definition"]["min_class"], 1)


class Identifiers(unittest.TestCase):
    """Rule 6 v2: strip, refuse empty or anything outside printable ASCII but an interior space,
    refuse NFKC-casefold aliases; otherwise keep as given."""

    def test_lineage_and_duplicate_ids_are_stripped_and_otherwise_kept(self):
        rows = lb.normalize_genomes([genome("x1", " L9 ", " D 1 ", None, []), genome("x2", "4.2.1", "Dx\n", None, []),
                                     genome("x3", "\tL8", None, None, []), genome("x4", "L9 x", "\tD2\r\n", None, [])])
        self.assertEqual([(r["lineage"], r["duplicate_group"]) for r in rows],
                         [("L9", "D 1"), ("4.2.1", "Dx"), ("L8", None), ("L9 x", "D2")])

    def test_refusals(self):
        cases = [(genome("x1", lineage, None, None, []), message) for lineage, message in REFUSED if "alias" not in message]
        cases += [(genome("x1", "L9", "D1" + lineage[2:], None, []), "duplicate_group .*" + message)
                  for lineage, message in REFUSED if lineage.startswith("L9")]
        cases += [
            (genome("x1", "\U0000200bL9", None, None, []), "lineage .* in genome 'x1' holds U\\+200B"),
            (genome("x1", "L\U000000ad9", None, None, []), "U\\+00AD"),
            (genome("x1", "L\U000000a09", None, None, []), "U\\+00A0"),
            (genome("x1", "L9", "D\n1", None, []), "U\\+000A"),
            (genome("x1", "   ", None, None, []), "lineage must be non-empty in genome 'x1'"),
            (genome("x1", "L9", " \t", None, []), "duplicate_group must be non-empty"),
            (genome("x1 ", "L9", None, None, []), "genome_id 'x1 ' has surrounding whitespace"),
            (genome("\tx1", "L9", None, None, []), "genome_id '\\\\tx1' has surrounding whitespace"),
            (genome("x\U0000200b1", "L9", None, None, []), "genome_id .* U\\+200B"),
            (genome("x1\U0000034f", "L9", None, None, []), "genome_id .* U\\+034F"),
            (genome("x\t1", "L9", None, None, []), "genome_id .* U\\+0009"),
        ]
        self.assertEqual(len(cases), 10 + 6 + 11)
        for record, message in cases:
            with self.subTest(record=record), self.assertRaisesRegex(ValueError, message):
                lb.normalize_genomes([record])

    def test_the_full_strip_set_is_stripped_from_lineages_and_groups_and_refused_around_genome_ids(self):
        self.assertEqual(len(STRIP_SET), 29)
        self.assertEqual(set(STRIP_SET), {chr(c) for c in range(0x110000) if chr(c).strip() == ""})
        for c in STRIP_SET:
            with self.subTest(char="U+%04X" % ord(c)):
                rows = lb.normalize_genomes([genome("x1", c + "L9", "D1" + c, None, []), genome("x2", "L8" + c, c + "D2", None, [])])
                self.assertEqual([(r["lineage"], r["duplicate_group"]) for r in rows], [("L9", "D1"), ("L8", "D2")])
                for gid in (c + "x1", "x1" + c):
                    with self.assertRaisesRegex(ValueError, "has surrounding whitespace"):
                        lb.normalize_genomes([genome(gid, "L9", None, None, [])])
                if c != " ":
                    with self.assertRaisesRegex(ValueError, "U\\+%04X" % ord(c)):
                        lb.normalize_genomes([genome("x1", "L" + c + "9", None, None, [])])

    def test_printable_ascii_bounds(self):
        # The re-check's T-2 and AGR-3: "!" (U+0021) and "~" (U+007E) are kept, DEL (U+007F) is refused,
        # each inside a value and at its edges.
        rows = lb.normalize_genomes([genome("x~1", "L~9", "D~1", None, []), genome("!x2~", "~L8!", "!D2~", None, [])])
        self.assertEqual([(r["genome_id"], r["lineage"], r["duplicate_group"]) for r in rows],
                         [("!x2~", "~L8!", "!D2~"), ("x~1", "L~9", "D~1")])
        for record in (genome("x1", "L\x7f9", None, None, []), genome("x1", "\x7fL9", None, None, []),
                       genome("x1", "L9", "D\x7f1", None, []), genome("x\x7f1", "L9", None, None, []),
                       genome("x1\x7f", "L9", None, None, [])):
            with self.subTest(record=record), self.assertRaisesRegex(ValueError, "U\\+007F"):
                lb.normalize_genomes([record])
        for c in map(chr, range(0x21, 0x7f)):
            with self.subTest(char=c):
                self.assertEqual(lb.identifier("a" + c + "b", "lineage"), "a" + c + "b")
                self.assertEqual(lb.identifier(c, "genome_id", strip=False), c)
        for c in map(chr, list(range(0x20)) + [0x7f, 0x80, 0xa0]):
            with self.subTest(char="U+%04X" % ord(c)), self.assertRaisesRegex(ValueError, "U\\+%04X" % ord(c)):
                lb.identifier("a" + c + "b", "lineage")

    def test_alias_collisions_name_both_values(self):
        cases = [
            ([genome("x1", "L9", None, None, []), genome("x2", "l9", None, None, [])], "lineage alias collision: 'L9' and 'l9'"),
            ([genome("x1", "L1", "D1", None, []), genome("x2", "L1", "d1", None, [])], "duplicate_group alias collision: 'D1' and 'd1'"),
            ([genome("X1", "L1", None, None, []), genome("x1", "L1", None, None, [])], "genome_id alias collision: 'X1' and 'x1'"),
            ([genome("x1", "L9", None, None, []), genome("x2", " l9\t", None, None, [])], "lineage alias collision: 'L9' and 'l9'"),
        ]
        for genomes, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                lb.normalize_genomes(genomes)
        with self.assertRaisesRegex(ValueError, "lineage alias collision: 'L4' and 'l4'"):
            define(holdout_lineages=["l4"])
        self.assertEqual(define(holdout_lineages=[" L4 ", "\tL9"])["definition"]["holdout_lineages"], ["L4", "L9"])

    def test_aliases_of_a_held_out_lineage_never_enter_development(self):
        base, g = reviewer_base()
        box = lb.define_lockbox(base, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        dev = lb.development_genomes(base, box)
        for alias in CLAIMED:
            with self.subTest(alias=alias):
                spelled = base + [g("alias", alias, None, "2009-01-01")]
                self.assertIn({"genome_id": "alias", "reasons": ["holdout_lineage"]}, lb.lockbox_conflicts(spelled, box))
                self.assertNotIn("alias", [r["genome_id"] for r in lb.development_genomes(spelled, box)])
                with self.assertRaisesRegex(lb.LockboxError, "alias \\(holdout_lineage\\)"):
                    sp.global_folds(dev + [g("alias", alias, None, "2009-01-01")], 5, 1, "random", lockbox=box)
        for alias, message in REFUSED:
            with self.subTest(alias=alias):
                with self.assertRaisesRegex(ValueError, message):
                    lb.development_genomes(base + [g("alias", alias, None, "2009-01-01")], box)
                with self.assertRaisesRegex(ValueError, message):
                    sp.global_folds(dev + [g("alias", alias, None, "2009-01-01")], 5, 1, "random", lockbox=box)
        kept = lb.development_genomes(base + [g("other", "L9 x", None, "2009-01-01")], box)
        self.assertEqual([r["lineage"] for r in kept if r["genome_id"] == "other"], ["L9 x"])

    def test_aliases_never_split_a_lineage_unit(self):
        base = [genome("y%02d" % i, "L%d" % (i % 6), None, None, {"INH": i % 2}) for i in range(30)]
        for alias in CLAIMED:
            with self.subTest(alias=alias):
                genomes = base + [genome("alias", alias, None, None, {"INH": 1}), genome("true9", "L9", None, None, {"INH": 0})]
                self.assertIn(["alias", "true9"], sp.units(genomes, "lineage"))
                folds = sp.global_folds(genomes, 5, 1, "lineage", lockbox=None)
                self.assertEqual(folds["alias"], folds["true9"])
        for alias, message in REFUSED:
            with self.subTest(alias=alias), self.assertRaisesRegex(ValueError, message):
                sp.global_folds(base + [genome("alias", alias, None, None, {"INH": 1}),
                                        genome("true9", "L9", None, None, {"INH": 0})], 5, 1, "lineage", lockbox=None)

    def test_stored_ids_collide_with_new_input(self):
        box = define()
        for record, message in ((genome("A2", "L6", None, None, {}), "genome_id alias collision: 'A2' and 'a2'"),
                                (genome("x1", "L6", "d2", None, {}), "duplicate_group alias collision: 'D2' and 'd2'"),
                                (genome("x1", "l4", None, None, {}), "lineage alias collision: 'L4' and 'l4'")):
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                lb.lockbox_conflicts([record], box)

    def test_stored_aliases_beyond_the_grouped_lockbox(self):
        # The re-check's LB-L1: a genome-id alias of h, a lockbox id in no group; and LB-L2: a group
        # alias of G, a definition group outside the lockbox. Each must be refused, not let through.
        r0 = release0()
        box = lb.define_lockbox(r0, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        self.assertEqual((box["genome_ids"], box["duplicate_groups"]), (["a", "h", "r1"], ["D1"]))
        self.assertNotIn("h", [gid for members in box["definition_groups"].values() for gid in members])
        cases = ((release(r0, ("H", "L0", None, "2010-01-01"), drop=("h",)), (), "genome_id alias collision: 'H' and 'h'"),
                 (r0[:30], [release(r0, ("H", "L0", None, "2010-01-01"), drop=("h",))], "genome_id alias collision: 'H' and 'h'"),
                 (release(r0, ("n", "L0", "g", "2010-01-01"), drop=("d1", "d2")), (), "duplicate_group alias collision: 'G' and 'g'"),
                 (release(r0, ("n", "L0", "g", "2010-01-01"), regroup={"d1": "D1", "d2": None}), (),
                  "duplicate_group alias collision: 'G' and 'g'"))
        for genomes, releases, message in cases:
            with self.subTest(message=message, releases=bool(releases)):
                for check in (lambda: lb.lockbox_conflicts(genomes, box, releases=releases),
                              lambda: lb.development_genomes(genomes, box, releases=releases),
                              lambda: sp.global_folds(genomes, 5, 1, "random", lockbox=box, releases=releases)):
                    with self.assertRaisesRegex(ValueError, message):
                        check()


class Seal(unittest.TestCase):
    def test_seal_is_sha256_of_canonical_json(self):
        box = define()
        body = {key: value for key, value in box.items() if key != "seal"}
        text = json.dumps(body, sort_keys=True, separators=(",", ":"))
        self.assertEqual(box["seal"], hashlib.sha256(text.encode("utf-8")).hexdigest())
        self.assertEqual(lb.verify_seal(box), box["seal"])
        self.assertEqual(lb.verify_seal(box, expected_seal=box["seal"]), box["seal"])

    def test_seal_is_stable_and_tracks_the_definition(self):
        rows = fixture()
        random.Random(3).shuffle(rows)
        self.assertEqual(define(rows, holdout_lineages=["L9", "L4", "L4"])["seal"], define()["seal"])
        self.assertNotEqual(define(recent_after="2020-07-01")["seal"], define()["seal"])
        self.assertNotEqual(define(min_class=1)["seal"], define()["seal"])
        self.assertNotEqual(define(access_log=UNWRITTEN_LOG + ".2")["seal"], define()["seal"])
        regrouped = fixture()
        regrouped[0]["duplicate_group"] = "D3"
        self.assertEqual(define(regrouped)["genome_ids"], define()["genome_ids"])
        self.assertNotEqual(define(regrouped)["seal"], define()["seal"])

    def test_tampered_manifest_is_refused_everywhere(self):
        box = define()
        tampered = copy.deepcopy(box)
        tampered["genome_ids"].remove("c1")
        for check in (lambda m: lb.verify_seal(m), lambda m: lb.assert_untouched(m, []),
                      lambda m: lb.lockbox_conflicts(fixture(), m), lambda m: lb.development_genomes(fixture(), m)):
            with self.assertRaisesRegex(lb.LockboxError, "tampered"):
                check(tampered)
        with self.assertRaisesRegex(lb.LockboxError, "not a GRM lockbox"):
            lb.verify_seal(dict(box, kind="other"))

    def test_expected_seal_is_compared_in_full(self):
        # The re-check's I2: a seal that differs only in its last digit must be refused.
        box = define()
        dev = lb.development_genomes(fixture(), box)
        near = near_miss(box["seal"])
        self.assertNotEqual(near, box["seal"])
        manifest = sp.global_manifest(dev, 2, 1, "random", lockbox=box, min_test=1, min_class=1)
        for check in (lambda: lb.verify_seal(box, expected_seal=near),
                      lambda: lb.require_development(dev, box, expected_seal=near),
                      lambda: lb.assert_untouched(box, [], expected_seal=near),
                      lambda: sp.global_folds(dev, 2, 1, "random", lockbox=box, expected_seal=near),
                      lambda: sp.verify_manifest(manifest, expected_seal=near)):
            with self.assertRaisesRegex(lb.LockboxError, "expected seal"):
                check()

    def test_resealed_forgery_is_refused_against_the_recorded_seal(self):
        # The review's case 5: h1 dropped and its lineage un-held, then resealed.
        base, g = reviewer_base()
        box = lb.define_lockbox(base, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        forged = resealed(box, genome_ids=["r1"], definition=dict(box["definition"], holdout_lineages=[]))
        self.assertEqual(lb.verify_seal(forged), forged["seal"])
        dev = base[:40] + [g("h1", "L9", None, "2011-01-01")]
        manifest = sp.global_manifest(dev, 5, 1, "random", lockbox=forged)
        self.assertIn("h1", manifest["folds"])
        for check in (lambda: lb.verify_seal(forged, expected_seal=box["seal"]),
                      lambda: lb.require_development(dev, forged, expected_seal=box["seal"]),
                      lambda: lb.development_genomes(dev, forged, expected_seal=box["seal"]),
                      lambda: lb.assert_untouched(forged, [], expected_seal=box["seal"]),
                      lambda: sp.verify_manifest(manifest, expected_seal=box["seal"]),
                      lambda: sp.global_folds(dev, 5, 1, "random", lockbox=forged, expected_seal=box["seal"]),
                      lambda: sp.global_manifest(dev, 5, 1, "random", lockbox=forged, expected_seal=box["seal"]),
                      lambda: sp.rolling_origin_dev_folds(dev, ["2009-06-30"], lockbox=forged, expected_seal=box["seal"]),
                      lambda: sp.rolling_origin_manifest(dev, ["2009-06-30"], lockbox=forged, expected_seal=box["seal"])):
            with self.assertRaisesRegex(lb.LockboxError, "expected seal"):
                check()
        self.assertEqual(sp.verify_manifest(manifest, expected_seal=forged["seal"]), manifest["sha256"])
        dev_true = lb.development_genomes(dev, box, expected_seal=box["seal"])
        self.assertNotIn("h1", sp.global_folds(dev_true, 5, 1, "random", lockbox=box, expected_seal=box["seal"]))
        with self.assertRaisesRegex(lb.LockboxError, "without a lockbox"):
            sp.global_folds(dev_true, 5, 1, "random", lockbox=None, expected_seal=box["seal"])

    def test_written_once_read_only_with_a_genesis_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, log = os.path.join(tmp, "lockbox.json"), os.path.join(tmp, "access.jsonl")
            box = define(access_log=log)
            written = lb.write_lockbox(box, path, actor="test", now=NOW)
            self.assertEqual(stat.S_IMODE(os.stat(path).st_mode), 0o444)
            self.assertEqual(lb.load_lockbox(path), box)
            self.assertEqual(lb.load_lockbox(path, expected_seal=box["seal"]), box)
            genesis = lb.read_access_log(log)
            self.assertEqual([(r["event"], r["seal"], r["prev_sha256"]) for r in genesis], [("seal", box["seal"], lb.LOG_GENESIS)])
            self.assertEqual(written, {"seal": box["seal"], "log_head": {"records": 1, "record_sha256": genesis[0]["record_sha256"]}})
            self.assertEqual(lb.verify_log(log, written["log_head"], seal=box["seal"]), written["log_head"])
            with self.assertRaises(FileExistsError):
                lb.write_lockbox(box, os.path.join(tmp, "again.json"), actor="test")
            self.assertFalse(os.path.exists(os.path.join(tmp, "again.json")))
            os.remove(log)
            with self.assertRaises(FileExistsError):
                lb.write_lockbox(box, path, actor="test")
            self.assertFalse(os.path.exists(log))
            os.chmod(path, 0o644)
            Path(path).write_text(Path(path).read_text().replace('"c1",', ""))
            with self.assertRaisesRegex(lb.LockboxError, "tampered"):
                lb.load_lockbox(path)
            forged = resealed(box, genome_ids=box["genome_ids"][1:], access_log=os.path.join(tmp, "forged.jsonl"))
            with self.assertRaisesRegex(lb.LockboxError, "expected seal"):
                lb.write_lockbox(forged, os.path.join(tmp, "other.json"), actor="test", expected_seal=box["seal"])
            self.assertFalse(os.path.exists(os.path.join(tmp, "other.json")))


class Guards(unittest.TestCase):
    def test_conflicts_name_each_rule(self):
        box = define()
        self.assertEqual([g["genome_id"] for g in lb.development_genomes(fixture(), box)], ["a1", "b2", "c3", "e1"])
        newcomers = [genome("x1", "L6", "D1", "2000-01-01", {}), genome("x2", "L4", None, "2001-01-01", {}),
                     genome("x3", "L6", None, "2024-01-01", {}), genome("x4", "L6", None, "2020-07", {}),
                     genome("x5", "L6", None, None, {})]
        self.assertEqual(lb.lockbox_conflicts(newcomers, box), [
            {"genome_id": "x1", "reasons": ["duplicate_group"]}, {"genome_id": "x2", "reasons": ["holdout_lineage"]},
            {"genome_id": "x3", "reasons": ["recent"]}, {"genome_id": "x4", "reasons": ["recent"]}])
        conflicts = {c["genome_id"]: c["reasons"] for c in lb.lockbox_conflicts(fixture(), box)}
        self.assertEqual(conflicts["d2"], ["lockbox_genome", "duplicate_group", "holdout_lineage"])
        self.assertEqual(conflicts["c2"], ["lockbox_genome", "recent"])
        with self.assertRaisesRegex(lb.LockboxError, "6 lockbox genomes cannot enter development"):
            lb.require_development(fixture(), box)
        self.assertEqual(lb.require_development(lb.development_genomes(fixture(), box), box), box["seal"])

    def test_duplicate_closure_over_a_new_release(self):
        # The review's H1 cases 1, 1c and 1d: each partner must be claimed, and every split refuses it.
        base, g = reviewer_base()
        box = lb.define_lockbox(base, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        self.assertEqual((box["genome_ids"], box["duplicate_groups"]), (["h1", "r1"], []))
        regrouped = copy.deepcopy(base)
        next(x for x in regrouped if x["genome_id"] == "r1")["duplicate_group"] = "D9"
        cases = {
            "r1copy": (regrouped + [g("r1copy", "L1", "D9", "2010-01-01")],
                       {"r1": ["lockbox_genome", "recent"], "r1copy": ["duplicate_closure"]}),
            "n_old": (base + [g("n_recent", "L2", "D7", "2023-05-05"), g("n_old", "L2", "D7", "2009-01-01")],
                      {"n_recent": ["recent"], "n_old": ["duplicate_closure"]}),
            "n_h_dup": (base + [g("n_h", "L9", "D8", "2009-01-01"), g("n_h_dup", "L2", "D8", "2009-01-01")],
                        {"n_h": ["holdout_lineage"], "n_h_dup": ["duplicate_closure"]}),
        }
        for leaked, (genomes, expected) in cases.items():
            with self.subTest(leaked=leaked):
                conflicts = {c["genome_id"]: c["reasons"] for c in lb.lockbox_conflicts(genomes, box)}
                for gid, reasons in expected.items():
                    self.assertEqual(conflicts[gid], reasons)
                dev = lb.development_genomes(genomes, box)
                self.assertNotIn(leaked, {x["genome_id"] for x in dev})
                self.assertEqual(len(dev), 40)
                for call in (lambda r: sp.global_folds(r, 5, 1, "random", lockbox=box),
                             lambda r: sp.global_manifest(r, 5, 1, "lineage", lockbox=box),
                             lambda r: sp.rolling_origin_dev_folds(r, ["2009-06-30"], lockbox=box)):
                    with self.assertRaisesRegex(lb.LockboxError, "%s \\(duplicate_closure\\)" % leaked):
                        call(genomes)
                sp.global_folds(dev, 5, 1, "random", lockbox=box)

    def assert_claimed(self, genomes, box, expected, releases=()):
        """expected {genome_id: reasons} are claimed and left out of development, and every split refuses the release."""
        conflicts = {c["genome_id"]: c["reasons"] for c in lb.lockbox_conflicts(genomes, box, releases=releases)}
        self.assertEqual({gid: conflicts.get(gid) for gid in expected}, expected)
        dev = lb.development_genomes(genomes, box, releases=releases)
        dev_ids = {x["genome_id"] for x in dev}
        self.assertFalse(dev_ids & set(expected))
        self.assertLessEqual({"f%02d" % i for i in range(30)}, dev_ids)
        shown = min(gid for gid, reasons in expected.items() if reasons[-1] in ("duplicate_closure", "earlier_release"))
        for call in (lambda r: sp.global_folds(r, 5, 1, "random", lockbox=box, releases=releases),
                     lambda r: sp.global_manifest(r, 3, 1, "lineage", lockbox=box, releases=releases),
                     lambda r: sp.rolling_origin_dev_folds(r, ["2009-06-30"], lockbox=box, releases=releases),
                     lambda r: sp.rolling_origin_manifest(r, ["2009-06-30"], lockbox=box, releases=releases)):
            with self.assertRaisesRegex(lb.LockboxError, "%s \\(%s\\)" % (shown, ", ".join(expected[shown]))):
                call(genomes)
            call(dev)

    def test_cross_release_chains_through_definition_groups(self):
        # The re-check's L1 cases e and f: an old development group G = {d1, d2} re-split by a later
        # release, with d1 moved in beside a lockbox genome; only the stored definition groups link d2.
        r0 = release0()
        box = lb.define_lockbox(r0, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        self.assertEqual(box["genome_ids"], ["a", "h", "r1"])
        self.assertEqual(box["definition_groups"], {"D1": ["a", "r1"], "G": ["d1", "d2"]})
        case_e = release(r0, ("e", "L0", "Y", "2010-01-01"), regroup={"d1": "D1", "d2": "Y"})
        self.assert_claimed(case_e, box, {"d1": ["duplicate_group"], "d2": ["duplicate_closure"], "e": ["duplicate_closure"]})
        case_f = release(r0, ("n", "L0", "N", "2024-01-01"), regroup={"d1": "N"})
        self.assert_claimed(case_f, box, {"n": ["recent"], "d1": ["duplicate_closure"], "d2": ["duplicate_closure"]})
        unstored = resealed(box, definition_groups={"D1": ["a", "r1"]})
        self.assertIn("d2", [x["genome_id"] for x in lb.development_genomes(case_f, unstored)])

    def test_cross_release_chains_through_earlier_releases(self):
        # The re-check's L1 case d: lockbox genome a moves in with b in R1, then leaves; R2 keeps b with c.
        r0 = release0()
        box = lb.define_lockbox(r0, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        r1 = release(r0, ("b", "L2", "E", "2010-01-01"), regroup={"a": "E"})
        self.assert_claimed(r1, box, {"a": ["lockbox_genome"], "b": ["duplicate_closure"]})
        r2 = release(r0, ("b", "L2", "E", "2010-01-01"), ("c", "L2", "E", "2010-01-01"), drop=("a",))
        self.assertTrue({"b", "c"} <= {x["genome_id"] for x in lb.development_genomes(r2, box)})  # R2 alone links nothing
        self.assert_claimed(r2, box, {"b": ["duplicate_closure"], "c": ["duplicate_closure"]}, releases=[r1])
        self.assertEqual(lb.lockbox_conflicts(r2, box, releases=iter([r1])), lb.lockbox_conflicts(r2, box, releases=[r1]))
        # A genome R1 claims as recent, added after definition, then dropped: its partners stay claimed.
        r1b = release(r0, ("n", "L0", "M", "2024-01-01"), ("x", "L0", "M", "2010-01-01"))
        r2b = release(r0, ("x", "L0", "M", "2010-01-01"), ("y", "L0", "M", "2010-01-01"))
        self.assert_claimed(r2b, box, {"x": ["duplicate_closure"], "y": ["duplicate_closure"]}, releases=[r1b])
        # The same genome back with a corrected, old date and no group: still claimed.
        r2c = release(r0, ("n", "L0", None, "2010-01-01"))
        self.assert_claimed(r2c, box, {"n": ["earlier_release"]}, releases=[r1, r1b])

    def test_releases_are_validated(self):
        r0 = release0()
        box = lb.define_lockbox(r0, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        for releases in (r0, r0[0], "r0", [r0, r0[0]]):
            with self.subTest(releases=type(releases)), self.assertRaisesRegex(ValueError, "list of releases"):
                lb.lockbox_conflicts(r0, box, releases=releases)
        for genomes, releases, message in (
                (release(r0, ("D1", "L0", None, None), drop=("d1",)), (), "genome_id alias collision: 'D1' and 'd1'"),
                (r0, [release(r0, regroup={"d1": "g", "d2": "g"})], "duplicate_group alias collision: 'G' and 'g'"),
                (r0, [release(r0, ("R1", "L0", None, None), drop=("r1",))], "genome_id alias collision: 'R1' and 'r1'"),
                (r0, [release(r0, ("z", "l9", None, None), drop=("h",))], "lineage alias collision: 'L9' and 'l9'"),
                (r0, [release(r0, ("z", "L\U0000034f9", None, None))], "U\\+034F")):
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                lb.lockbox_conflicts(genomes, box, releases=releases)
        with self.assertRaisesRegex(lb.LockboxError, "without a lockbox"):
            sp.global_folds(r0, 5, 1, "random", lockbox=None, releases=[r0])

    def test_assert_untouched_catches_a_planted_leak(self):
        box = define()
        clean = ["a1", "b2", "c3", "e1"]
        lb.assert_untouched(box, clean)
        lb.assert_untouched(box, iter(clean))
        with self.assertRaisesRegex(lb.LockboxError, r"1 lockbox genomes were used: c1$"):
            lb.assert_untouched(box, clean + ["c1"])
        with self.assertRaises(TypeError):
            lb.assert_untouched(box, "c1")

    def test_assert_untouched_applies_the_genome_id_rule(self):
        # The re-check's LB-M1: assert_untouched matched exactly, so variants of a lockbox id that every
        # other guard refuses passed. Lockbox ids a, h (ungrouped) and r1; d1 is a stored development member.
        r0 = release0()
        box = lb.define_lockbox(r0, recent_after="2020-12-31", access_log=UNWRITTEN_LOG, holdout_lineages=["L9"],
                                min_per_task=1, min_class=1)
        dev = [x["genome_id"] for x in lb.development_genomes(r0, box)]
        self.assertIn("d1", dev)
        lb.assert_untouched(box, dev)
        lb.assert_untouched(box, iter(dev), expected_seal=box["seal"])
        lb.assert_untouched(box, [])
        for leak in ("r1", "h", "a"):
            with self.subTest(leak=leak), self.assertRaisesRegex(lb.LockboxError, r"1 lockbox genomes were used: %s$" % leak):
                lb.assert_untouched(box, dev + [leak])
        # The re-check's PLB-L1: a lockbox id inside a one-shot iterator or generator, read once.
        for leak, used in (("r1", lambda: iter(dev + ["r1"])), ("h", lambda: (gid for gid in ["h"])),
                           ("a", lambda: map(str, dev + ["a"]))):
            with self.subTest(one_shot=leak), self.assertRaisesRegex(lb.LockboxError, r"1 lockbox genomes were used: %s$" % leak):
                lb.assert_untouched(box, used())
        with self.assertRaisesRegex(ValueError, "genome_id alias collision: 'R1' and 'r1'"):
            lb.assert_untouched(box, (gid for gid in ["d1", "R1"]))
        # PLB-L2: two used ids that alias each other, neither stored, as rule 6 v2 refuses in one input.
        for used in (["f01", "F01"], dev + ["dev07", "DEV07"], iter(["dev07", "DEV07"])):
            with self.subTest(used=used), self.assertRaisesRegex(ValueError, "genome_id alias collision: '(F01' and 'f01|DEV07' and 'dev07)'"):
                lb.assert_untouched(box, used)
        with self.assertRaisesRegex(ValueError, "genome_id alias collision"):
            lb.lockbox_conflicts([genome("dev07", "L0", None, None, {}), genome("DEV07", "L0", None, None, {})], box)
        for used in (b"r1", b""):
            with self.subTest(used=used), self.assertRaisesRegex(TypeError, "not one string"):
                lb.assert_untouched(box, used)
        cases = (("R1", "genome_id alias collision: 'R1' and 'r1'"), ("r1 ", "'r1 ' has surrounding whitespace"),
                 (" r1", "' r1' has surrounding whitespace"), ("r1​", "U\\+200B"), ("ｒ１", "U\\+FF52"),
                 (1, "genome_id must be a string, got 1"), (None, "genome_id must be a string"),
                 ("H", "genome_id alias collision: 'H' and 'h'"), ("h\t", "has surrounding whitespace"),
                 ("A", "genome_id alias collision: 'A' and 'a'"), ("D1", "genome_id alias collision: 'D1' and 'd1'"),
                 ("", "genome_id must be non-empty"), ("r\x7f1", "U\\+007F"))
        for used, message in cases:
            with self.subTest(used=used):
                with self.assertRaisesRegex(ValueError, message):
                    lb.assert_untouched(box, [used])
                with self.assertRaisesRegex(ValueError, message):
                    lb.assert_untouched(box, dev + [used])
                if isinstance(used, str):
                    with self.assertRaises(ValueError):
                        lb.development_genomes(release(r0, (used, "L0", None, "2010-01-01"), drop=("r1", "h", "a", "d1")), box)


class AccessLog(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.log = os.path.join(self.tmp.name, "access.jsonl")
        self.box = define(access_log=self.log)
        self.seal = self.box["seal"]
        self.genesis = lb.write_lockbox(self.box, os.path.join(self.tmp.name, "lockbox.json"), actor="test", now=NOW)
        self.head = self.genesis["log_head"]

    def unseal(self, manifest=None, approval="SYNTH-APPROVAL-1", actor="test-runner", expected=None, head=None):
        """unseal with the head this test recorded, which it updates after every recorded attempt."""
        try:
            granted = lb.unseal(self.box if manifest is None else manifest, "synthetic confirmation", approval,
                                actor=actor, expected_seal=self.seal if expected is None else expected,
                                expected_head=self.head if head is None else head, now=NOW)
        except lb.LockboxError as error:
            self.head = error.log_head or self.head
            raise
        self.head = granted["log_head"]
        return granted

    def assert_refused_unrecorded(self, message, **kwargs):
        before = Path(self.log).read_bytes() if os.path.exists(self.log) else None
        with self.assertRaisesRegex(lb.LockboxError, "not recorded: .*" + message) as caught:
            self.unseal(**kwargs)
        self.assertIsNone(caught.exception.log_head)
        self.assertEqual(Path(self.log).read_bytes() if os.path.exists(self.log) else None, before)

    def test_unseal_needs_approval_and_logs_every_attempt(self):
        for approval in (None, "", "  "):
            with self.assertRaisesRegex(lb.LockboxError, "missing approval_ref"):
                self.unseal(approval=approval)
        with self.assertRaisesRegex(lb.LockboxError, "missing actor"):
            self.unseal(actor="")
        with self.assertRaisesRegex(lb.LockboxError, "missing expected_seal"):
            self.unseal(expected="")
        with self.assertRaisesRegex(lb.LockboxError, "tampered"):
            self.unseal(dict(self.box, genome_ids=["a1"]))
        with self.assertRaisesRegex(lb.LockboxError, "expected seal") as caught:
            self.unseal(expected="0" * 64)
        self.assertEqual(caught.exception.log_head, {"records": 8, "record_sha256": lb.read_access_log(self.log)[-1]["record_sha256"]})
        granted = self.unseal()
        self.assertEqual(granted["genome_ids"], self.box["genome_ids"])
        self.assertEqual(granted["log_head"], lb.verify_log(self.log, seal=self.seal))
        records = lb.read_access_log(self.log)
        self.assertEqual(granted["log_head"], {"records": 9, "record_sha256": records[-1]["record_sha256"]})
        self.assertEqual([r["outcome"] for r in records], ["written"] + ["refused"] * 7 + ["granted"])
        last = records[-1]
        self.assertEqual({k: last[k] for k in ("event", "time", "actor", "purpose", "approval_ref", "seal", "reason")},
                         {"event": "unseal", "time": "2026-10-30T12:00:00+00:00", "actor": "test-runner",
                          "purpose": "synthetic confirmation", "approval_ref": "SYNTH-APPROVAL-1",
                          "seal": self.seal, "reason": None})
        self.assertEqual(records[1]["approval_ref"], None)
        self.assertIn("tampered", records[6]["reason"])
        self.assertIn("expected seal", records[7]["reason"])
        for before, after in zip(records, records[1:]):
            self.assertEqual(after["prev_sha256"], before["record_sha256"])

    def test_no_grant_without_the_lockbox_log(self):
        # The review's case 6: a deleted log, or a log that was never started, must not start a fresh chain.
        unwritten = define(access_log=os.path.join(self.tmp.name, "never.jsonl"))
        with self.assertRaisesRegex(lb.LockboxError, "not recorded: .*is missing"):
            lb.unseal(unwritten, "p", "SYNTH-1", actor="a", expected_seal=unwritten["seal"], expected_head=self.head)
        self.assertFalse(os.path.exists(unwritten["access_log"]))
        os.remove(self.log)
        self.assert_refused_unrecorded("is missing")
        self.assertFalse(os.path.exists(self.log))
        with self.assertRaises(TypeError):
            lb.unseal(self.box, "p", "SYNTH-1", actor="a", expected_seal=self.seal, expected_head=self.head,
                      log_path=self.log)

    def test_log_started_for_another_seal_is_refused(self):
        other = define(access_log=self.log, min_class=1)
        with self.assertRaisesRegex(lb.LockboxError, "not recorded: .*another seal"):
            lb.unseal(other, "p", "SYNTH-1", actor="a", expected_seal=other["seal"], expected_head=self.head)
        self.assertEqual(len(lb.read_access_log(self.log)), 1)

    def test_unseal_requires_the_recorded_head(self):
        # The re-check's L3: a re-created log, a log cut back to its first record and a hand-written
        # first record were all granted without a recorded head.
        with self.assertRaises(TypeError):
            lb.unseal(self.box, "p", "SYNTH-1", actor="a", expected_seal=self.seal)
        for bad in (None, {}, {"records": 0, "record_sha256": "0" * 64}, {"records": 1}, {"records": True,
                    "record_sha256": self.head["record_sha256"]}, dict(self.head, extra=1), [1, self.head["record_sha256"]]):
            with self.subTest(head=bad):
                self.assert_refused_unrecorded("expected_head", head=bad if bad is not None else {})
        first = self.unseal()["log_head"]
        second = self.unseal()["log_head"]
        self.assertEqual(second["records"], 3)
        self.assert_refused_unrecorded("1 records past the expected head", head=first)
        self.assert_refused_unrecorded("no longer holds the recorded head", head=dict(second, records=4))
        lines = Path(self.log).read_text().splitlines(keepends=True)

        Path(self.log).write_text(lines[0])                                         # cut back to the first record
        self.assertEqual(lb.verify_log(self.log, seal=self.seal), self.genesis["log_head"])
        self.assert_refused_unrecorded("no longer holds the recorded head")

        os.remove(self.log)                                                         # re-created by write_lockbox
        recreated = lb.write_lockbox(self.box, os.path.join(self.tmp.name, "again.json"), actor="test",
                                     now=NOW + dt.timedelta(hours=1))
        self.assertNotEqual(recreated["log_head"], self.genesis["log_head"])
        self.assert_refused_unrecorded("no longer holds the recorded head")
        self.assert_refused_unrecorded("no longer holds the recorded head", head=self.genesis["log_head"])

        os.remove(self.log)                                                         # a hand-written first record
        forged = {"event": "seal", "time": "2026-10-30T13:00:00+00:00", "actor": "forger", "purpose": "define lockbox",
                  "approval_ref": None, "seal": self.seal, "outcome": "written", "reason": None, "prev_sha256": lb.LOG_GENESIS}
        Path(self.log).write_text(lb.canonical_json(dict(forged, record_sha256=lb.sha256_json(forged))) + "\n")
        self.assertEqual(lb.verify_log(self.log, seal=self.seal)["records"], 1)
        self.assert_refused_unrecorded("no longer holds the recorded head")
        self.assert_refused_unrecorded("no longer holds the recorded head", head=self.genesis["log_head"])

    def test_recorded_head_detects_truncation_and_rewrite(self):
        first = self.unseal()["log_head"]
        second = self.unseal()["log_head"]
        self.assertEqual(lb.verify_log(self.log, first, seal=self.seal), second)
        lines = Path(self.log).read_text().splitlines(keepends=True)
        Path(self.log).write_text("".join(lines[:-1]))
        self.assertEqual(lb.verify_log(self.log, seal=self.seal), first)            # the chain alone reads clean
        with self.assertRaisesRegex(lb.LockboxError, "no longer holds the recorded head"):
            lb.verify_log(self.log, second, seal=self.seal)
        # The review's whole-chain rewrite with forged approvals: clean as a chain, caught by a recorded head.
        records, previous = lb.read_access_log(self.log), lb.LOG_GENESIS
        with open(self.log, "w") as stream:
            for record in records:
                body = {k: v for k, v in record.items() if k not in ("prev_sha256", "record_sha256")}
                if body["event"] == "unseal":
                    body["approval_ref"] = "FORGED"
                body["prev_sha256"] = previous
                previous = lb.sha256_json(body)
                stream.write(lb.canonical_json(dict(body, record_sha256=previous)) + "\n")
        self.assertEqual([r["approval_ref"] for r in lb.read_access_log(self.log)], [None, "FORGED"])
        self.assertEqual(lb.verify_log(self.log, self.genesis["log_head"], seal=self.seal), lb.verify_log(self.log, seal=self.seal))
        with self.assertRaisesRegex(lb.LockboxError, "no longer holds the recorded head"):
            lb.verify_log(self.log, first, seal=self.seal)
        for bad in ({"records": 0, "record_sha256": "x"}, {"records": 1}, {}):
            with self.assertRaises(ValueError):
                lb.verify_log(self.log, bad, seal=self.seal)

    def test_verify_log_requires_the_seal(self):
        # The re-check's I6: seal=None skipped the genesis binding.
        with self.assertRaises(TypeError):
            lb.verify_log(self.log)
        for bad in (None, "", self.seal[:12], self.seal.upper()):
            with self.subTest(seal=bad), self.assertRaisesRegex(ValueError, "seal must be"):
                lb.verify_log(self.log, seal=bad)
        with self.assertRaisesRegex(lb.LockboxError, "another seal"):
            lb.verify_log(self.log, seal=near_miss(self.seal))
        self.assertEqual(lb.verify_log(self.log, seal=self.seal), self.genesis["log_head"])

    def test_edits_inside_the_chain_block_unseal(self):
        self.unseal()
        first = Path(self.log).read_bytes()
        with self.assertRaises(lb.LockboxError):
            self.unseal(approval=None)
        self.unseal()
        grown = Path(self.log).read_bytes()
        self.assertTrue(grown.startswith(first))
        self.assertEqual(grown.count(b"\n"), 4)

        lines = grown.decode().splitlines(keepends=True)
        edited = lines[1].replace("SYNTH-APPROVAL-1", "SYNTH-APPROVAL-2")
        for broken, message in (([lines[0], edited] + lines[2:], "breaks the hash chain"),
                                (lines[2:], "breaks the hash chain"),
                                (lines[:2] + lines[3:], "breaks the hash chain"),
                                ([grown.decode()[:-1]], "breaks the hash chain"),
                                (lines[1:2], "breaks the hash chain"),
                                ([], "does not start with a seal record")):
            with self.subTest(lines=len(broken)):
                Path(self.log).write_text("".join(broken))
                with self.assertRaisesRegex(lb.LockboxError, message):
                    lb.read_access_log(self.log)
                self.assert_refused_unrecorded(message)


if __name__ == "__main__":
    unittest.main()
