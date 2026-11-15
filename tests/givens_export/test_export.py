"""GRM's export for Givens (reference/givens-export/grm.json) and the script that writes it.

The file must be exactly what scripts/givens-export/export.py writes from the committed sources. The suite then holds
the file to the rules Givens' own check applies to a module, to the documents' sentences, and to the export's limits:
a number is read from a file or quoted from a document, a missing count says so, and no dataset, species, antibiotic
or sequence is named. The planted faults run the script on a changed copy of its spec or of its source files. The
tests name a few species and antibiotics, because each of those names must stop the script.
"""

import contextlib
import copy
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
KIT = ROOT / "scripts" / "givens-export"
FILE = ROOT / "reference" / "givens-export" / "grm.json"
SUMS = ROOT / "reference" / "givens-export" / "SHA256SUMS"


def load_export():
    spec = importlib.util.spec_from_file_location("givens_export", KIT / "export.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load_export()
SPEC = json.loads((KIT / "spec.json").read_text(encoding="utf-8"))
MODULE = json.loads(FILE.read_text(encoding="utf-8"))
ENTRIES = {x["id"]: x for x in MODULE["nodes"] + MODULE["assumptions"]}
LINES = [(x["id"], line) for x in ENTRIES.values() for line in x["measured"]]


def text(path):
    return (ROOT / path).read_text(encoding="utf-8")


def line(entry, what):
    """The one measured line of an entry whose text holds `what`."""
    found = [x for x in ENTRIES[entry]["measured"] if what in x["what"]]
    assert len(found) == 1, (entry, what, len(found))
    return found[0]


def of(entry, what):
    x = line(entry, what)
    return x["count"], x["of"]


def changed(edit):
    """The spec after one edit."""
    spec = copy.deepcopy(SPEC)
    edit(spec)
    return spec


def node(spec, node_id):
    return next(n for n in spec["nodes"] if n["id"] == node_id)


def copied(tmp):
    """A directory that holds a copy of every committed file the script reads, and nothing else of the repository."""
    root = pathlib.Path(tmp)
    for p in list(SPEC["files"].values()) + MODULE["sources"]["documents"] + [E.VOCABULARY]:
        (root / p).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / p, root / p)
    (root / "tests" / "outside_study_results").mkdir(parents=True)
    (root / "reference" / "givens-export").mkdir(parents=True)
    return root


def data_names():
    """The names the screen takes from the files the export reads."""
    sources = E.Sources(ROOT, SPEC["files"])
    for alias in SPEC["files"]:
        sources.data(alias)
    return E._barred(ROOT, sources)


def given(spec, tail):
    return next(g for g in spec["assumptions"] if g["id"] == f"grm:given/{tail}")


class Committed(unittest.TestCase):
    def test_the_file_is_what_the_script_writes(self):
        self.assertEqual(FILE.read_text(encoding="utf-8"), E.build())
        self.assertEqual(E.main(["export.py", "--check", str(FILE)]), 0)
        self.assertEqual(E.main(["export.py", "--check"]), 0)
        self.assertEqual(MODULE, E.build_module(SPEC))

    def test_the_hash_list_is_the_file_s(self):
        # One line, as sha256sum prints it from the repository root, so Givens can store the list as it is.
        digest = hashlib.sha256(FILE.read_bytes()).hexdigest()
        self.assertEqual(SUMS.read_text(encoding="utf-8"), f"{digest}  reference/givens-export/grm.json\n")
        self.assertEqual(SUMS.read_text(encoding="utf-8"), E.sums(FILE.read_text(encoding="utf-8")))
        readme = text("reference/givens-export/README.md")
        self.assertIn("`SHA256SUMS`", readme)
        self.assertIsNone(re.search(r"[0-9a-f]{64}", readme), "the hash stands in SHA256SUMS alone")

    def test_every_source_file_is_listed_with_its_sha256(self):
        files = MODULE["sources"]["files"]
        self.assertEqual([f["path"] for f in files], sorted(set(SPEC["files"].values())))
        for f in files:
            with self.subTest(file=f["path"]):
                self.assertEqual(hashlib.sha256((ROOT / f["path"]).read_bytes()).hexdigest(), f["sha256"])
                self.assertTrue(f["path"].startswith("reference/release-summary/") or f["path"] == E.OUTSIDE)
        self.assertEqual(MODULE["sources"]["documents"], sorted({q["path"] for x in ENTRIES.values() for q in x["quotes"]}
                                                                | {x["source"] for _, x in LINES if x["basis"] == "document"}))

    def test_the_file_names_no_commit_and_no_time(self):
        self.assertEqual(list(MODULE), ["format", "id", "namespace", "name", "job", "repo", "release", "binding",
                                        "binding_note", "sources", "nodes", "edges", "assumptions", "import_rule",
                                        "left_out"])
        self.assertEqual((MODULE["format"], MODULE["binding"]), ("grm-givens-export/1", "typed_export"))
        self.assertIsNone(re.search(r"\b20\d\d-\d\d-\d\dT|\b[0-9a-f]{40}\b", FILE.read_text(encoding="utf-8")))


class Shape(unittest.TestCase):
    """The rules Givens' tools/validate.mjs applies to a module, so the file needs no repair there."""

    def test_nodes(self):
        self.assertEqual(len(MODULE["nodes"]), len({n["id"] for n in MODULE["nodes"]}))
        for n in MODULE["nodes"]:
            with self.subTest(node=n["id"]):
                self.assertTrue(n["id"].startswith("grm:") and not n["id"].startswith("grm:given/"))
                self.assertIn(n["kind"], E.KINDS)
                self.assertTrue(n["label"] and n["definition"] and n["source"] and n["quotes"])
                self.assertIn(n.get("truth_state", "asserted"), E.TRUTH_STATES)
                self.assertLessEqual(set(n), {"id", "kind", "label", "definition", "truth_state", "source", "quotes",
                                              "measured", "unmeasured"})

    def test_edges(self):
        ids = {n["id"] for n in MODULE["nodes"]}
        for e in MODULE["edges"]:
            with self.subTest(edge=e):
                self.assertEqual(set(e), {"from", "type", "to"})
                self.assertIn(e["type"], E.EDGE_TYPES)
                self.assertTrue(e["from"] in ids and e["to"] in ids)
        self.assertEqual(len(MODULE["edges"]), len({json.dumps(e, sort_keys=True) for e in MODULE["edges"]}))

    def test_assumptions(self):
        ids = {n["id"] for n in MODULE["nodes"]}
        for g in MODULE["assumptions"]:
            with self.subTest(given=g["id"]):
                self.assertTrue(g["id"].startswith("grm:given/"))
                self.assertTrue(g["statement"] and g["hidden_in"] and g["when_false"] and g["under"])
                self.assertLessEqual(set(g["under"]), ids)
                self.assertTrue(g["measured"] or g.get("unmeasured"))

    def test_measured_lines(self):
        for entry, x in LINES:
            with self.subTest(entry=entry, what=x["what"]):
                self.assertTrue(x["what"] and x["source"])
                self.assertLessEqual(set(x), {"count", "of", "what", "source", "basis", "evaluation_revision", "read",
                                              "quote", "rebuilt_by", "values", "note"})
                for n in (x["count"], x["of"]):
                    self.assertTrue(n is None or (type(n) is int and n >= 0))
                if x["count"] is not None and x["of"] is not None:
                    self.assertLessEqual(x["count"], x["of"])
                if x["of"] is None:
                    self.assertTrue(x.get("note"))
                self.assertTrue(x["count"] is not None or x.get("values"))


class Numbers(unittest.TestCase):
    def test_every_line_is_read_from_a_file_or_quoted_from_a_document(self):
        listed = {f["path"] for f in MODULE["sources"]["files"]}
        for entry, x in LINES:
            with self.subTest(entry=entry, what=x["what"]):
                if x["basis"] == "file":
                    self.assertLessEqual(set(x["source"].split("; ")), listed)
                    self.assertTrue(x["read"] and "quote" not in x)
                    for where in list(x["read"].values()) + [v["read"] for v in x.get("values", [])]:
                        self.assertTrue(any(p in where for p in listed), where)
                else:
                    self.assertEqual(x["basis"], "document")
                    self.assertIn(x["quote"], text(x["source"]))
                    stated = E._numbers(x["quote"])
                    self.assertIn(str(x["of"]), stated)
                    self.assertTrue(str(x["count"]) in stated or x["count"] == 48 + 25, x["count"])
                    self.assertNotIn("read", x)

    def test_a_line_is_of_the_revision_its_files_are_of(self):
        revision = {f["path"]: (12 if f["path"] == E.OUTSIDE else int(E.SUMMARY_DIR.match(f["path"]).group(1)))
                    for f in MODULE["sources"]["files"]}
        for path, n in revision.items():
            doc = json.loads(text(path))
            self.assertEqual(doc.get("evaluation_revision", doc.get("revision")), n, path)
        for entry, x in LINES:
            if x["basis"] != "file":
                continue
            with self.subTest(entry=entry, what=x["what"]):
                if "values" in x:
                    self.assertEqual([v["evaluation_revision"] for v in x["values"]], [12, 14])
                    for v in x["values"]:
                        self.assertEqual({revision[p] for p in revision if p in v["read"]}, {v["evaluation_revision"]})
                else:
                    self.assertEqual({revision[p] for p in x["source"].split("; ")}, {x["evaluation_revision"]})

    def test_the_counts_add_up(self):
        same = "grm:given/learners-scored-on-same-genomes"
        left_out, rows = of(same, "left genomes without a call")
        self.assertEqual((left_out, rows), (122, 2710))
        self.assertEqual(of(same, "boosted trees' rows")[0] + of(same, "lasso's rows")[0], left_out)
        self.assertEqual(of(same, "scored on every genome under the flag"), (left_out, rows))
        self.assertEqual(of(same, "still count fewer genomes"), (0, rows))
        self.assertEqual(of("grm:state/scored_at_prevalence", "carry the flag"), (left_out, rows))
        self.assertEqual(of("grm:learner_row", "learner rows in the release"), (rows, rows))
        random, datasets = of("grm:split/random", "random split")
        lineage = of("grm:split/lineage", "lineage split")[0]
        self.assertEqual(len(E.LEARNERS) * (random + lineage), rows)
        self.assertEqual(lineage + of("grm:state/split_unavailable", "no lineage split")[0], datasets)
        reported, curated = of("grm:dataset", "have a report")
        self.assertEqual((reported, reported + of("grm:state/no_report", "no report")[0]), (datasets, curated))
        read, held = of("grm:outside_population/primary", "can be read")
        self.assertEqual(read + of("grm:reading/not_readable", "cannot be read")[0], held)

    def test_the_counts_are_the_documents(self):
        new = "grm:given/figure-holds-on-new-studies"
        primary, same = "grm:outside_population/primary", "grm:outside_population/same_study"
        results = text("docs/outside-study-results.md")
        read, datasets = of(primary, "can be read")
        scm, cart, lasso = (of(new, f"{name} holds, {where}") for name, where in (
            ("set covering machine", "on genomes from studies the cohort does not hold"),
            ("decision tree", "same population"), ("lasso", "of the datasets it can be scored and read on, same")))
        above, below = of("grm:reading/below", "reads above"), of("grm:reading/below", "reads below")
        scored, frozen = of(same, "genomes scored")
        scm2, cart2, lasso2 = (of(new, f"{name} holds, {where}on new genomes from studies the cohort already holds")
                               for name, where in (("set covering machine", ""), ("decision tree", ""),
                                                   ("lasso", "of the datasets it can be scored and read on, ")))
        sentences = [
            f"The primary population is {of(primary, 'genomes scored')[1]:,} genomes on {datasets} datasets in "
            f"{of(primary, 'species')[0]} species.",
            f"so {read} datasets can be read.",
            f"The set covering machine holds on {scm[0]} of the {scm[1]} datasets and reads below on the other "
            f"{scm[1] - scm[0]}.",
            f"The decision tree holds on {cart[0]} of the {cart[1]} datasets.",
            f"The lasso can be scored and read on {lasso[1]} of the {read} datasets, and holds on {lasso[0]} of them.",
            f"sensitivity reads above the report's figure on {above[0]} of the {above[1]} datasets and below it on "
            f"{below[0]}.",
            f"The same-study population is {frozen:,} genomes on {of(same, 'can be read')[1]} datasets in "
            f"{of(same, 'species')[0]} species.",
            f"{scored:,} of the {frozen:,} genomes were scored.",
            f"The set covering machine holds on {scm2[0]} of the {scm2[1]} datasets.",
            f"The decision tree holds on {cart2[0]} of the {cart2[1]} datasets.",
            f"The lasso holds on {lasso2[0]} of the {lasso2[1]} same-study datasets read for it against the figure "
            f"over fitted folds",
        ]
        for sentence in sentences:
            with self.subTest(sentence=sentence):
                self.assertIn(sentence, results)
        same_genomes = "grm:given/learners-scored-on-same-genomes"
        flagged = of(same_genomes, "scored on every genome under the flag")[0]
        self.assertIn(f"In revision 14, the {flagged} learner rows that carry the flag are scored on every genome",
                      text("docs/release-summary.md"))
        for name in ("boosted trees", "lasso"):
            values = line(same_genomes, f"{name}' median" if name.endswith("s") else f"{name}'s median")["values"]
            self.assertNotEqual(values[0]["value"], values[1]["value"], name)
            self.assertTrue(all(type(v["value"]) is float and 0.5 < v["value"] < 1 for v in values))

    def test_a_count_follows_its_source_file(self):
        # A changed source file changes the line that reads it and its sha256, and nothing else.
        def flip(doc):
            row = next(r for r in doc["populations"][0]["rows"] if r["learner"] == "scm"
                       and r["sensitivity"] == r["specificity"] == "consistent")
            row["specificity"] = "below"

        def one_less(doc):
            doc["flagged_rows"] -= 1

        for path, edit, entry, what in ((E.OUTSIDE, flip, "grm:given/figure-holds-on-new-studies",
                                         "set covering machine holds, on genomes"),
                                        (SPEC["files"]["index14"], one_less, "grm:state/scored_at_prevalence",
                                         "carry the flag")):
            with self.subTest(file=path), tempfile.TemporaryDirectory() as tmp:
                root = copied(tmp)
                doc = json.loads((root / path).read_text(encoding="utf-8"))
                edit(doc)
                (root / path).write_text(json.dumps(doc), encoding="utf-8")
                built = E.build_module(SPEC, root)
                got = next(x for x in next(e for e in built["nodes"] + built["assumptions"] if e["id"] == entry)[
                    "measured"] if what in x["what"])
                self.assertEqual(got["count"], line(entry, what)["count"] - 1)
                before = {f["path"]: f["sha256"] for f in MODULE["sources"]["files"]}
                after = {f["path"]: f["sha256"] for f in built["sources"]["files"]}
                self.assertEqual([p for p in before if before[p] != after[p]], [path])


class Missing(unittest.TestCase):
    def test_what_grm_has_not_counted_says_so(self):
        mic = ENTRIES["grm:state/label_from_mic_excluded"]
        self.assertEqual(mic["measured"], [])
        self.assertIn("No committed file counts", mic["unmeasured"])
        folds = ENTRIES["grm:given/folds-hold-different-genomes"]
        self.assertIn("No committed file counts duplicate groups on the cohorts of revision 10 or later",
                      folds["unmeasured"])
        self.assertEqual({x["evaluation_revision"] for x in folds["measured"]}, {9})
        self.assertEqual(sorted(i for i, x in ENTRIES.items() if "unmeasured" in x),
                         ["grm:given/folds-hold-different-genomes", "grm:state/label_from_mic_excluded"])

    def test_a_number_stands_in_a_measured_line_only(self):
        # 5-fold, a 95% interval, and the revisions 9, 10, 12, 13 and 14 are the only numbers the words may hold.
        allowed = {"5", "95", "9", "10", "12", "13", "14"}
        for entry in ENTRIES.values():
            for key in ("label", "definition", "statement", "hidden_in", "when_false", "unmeasured"):
                with self.subTest(entry=entry["id"], field=key):
                    self.assertLessEqual(E._numbers(entry.get(key, "")), allowed)
        for sentence in MODULE["import_rule"] + MODULE["left_out"] + [MODULE["job"], MODULE["binding_note"]]:
            self.assertLessEqual(E._numbers(sentence), allowed | {"256"}, sentence)

    def test_a_document_line_names_the_suite_that_rebuilds_it_or_none(self):
        rebuilt = {(entry, x["count"], x["of"]) for entry, x in LINES if "rebuilt_by" in x}
        self.assertEqual(rebuilt, {("grm:given/figure-holds-on-new-studies", 8, 25),
                                   ("grm:given/figure-holds-on-new-studies", 7, 15)})
        for _, x in LINES:
            if "rebuilt_by" in x:
                self.assertEqual((x["basis"], x["rebuilt_by"]), ("document", "tests/outside_study_results"))
                self.assertTrue((ROOT / x["rebuilt_by"]).is_dir())

    def test_an_assumption_with_no_count_must_say_unmeasured(self):
        def strip(spec):
            g = given(spec, "folds-hold-different-genomes")
            g["measured"] = []
            del g["unmeasured"]

        with self.assertRaisesRegex(E.Fault, "says so in unmeasured"):
            E.build_module(changed(strip))
        with self.assertRaisesRegex(E.Fault, "unmeasured says what GRM has not counted"):
            E.build_module(changed(lambda s: given(s, "folds-hold-different-genomes").update(unmeasured=" ")))


class Limits(unittest.TestCase):
    """Each planted fault must stop the script."""

    def refused(self, message, edit):
        with self.assertRaisesRegex(E.Fault, message):
            E.build_module(changed(edit))

    def test_the_file_as_committed_passes_the_screen(self):
        barred = data_names()
        self.assertGreater(len(barred), 250)
        for name in ("escherichia coli", "e. coli", "klebsiella_pneumoniae", "meropenem", "carbapenem",
                     "escherichia_coli__amikacin"):
            self.assertIn(name, barred)
        vocabulary = E._vocabulary(ROOT)
        self.assertGreater(len(vocabulary - barred), 250)
        E.screen(MODULE, barred | vocabulary)

    def test_a_name_of_the_vocabulary_is_refused(self):
        # Names the files read do not hold: the screen has them from src/drugs.rs alone.
        only = E._vocabulary(ROOT) - data_names()
        cases = {
            "an antibiotic of no outside-study dataset": "tylosin",
            "a spelling the vocabulary folds": "rifampicin",
            "a brand": "synercid",
            "one component of a combination": "quinupristin",
            "a drug class": "glycopeptide",
            "a drug class in the plural": "glycopeptides",
            "a class label in the singular": "fluoroquinolone",
            "a class of several words": "siderophore cephalosporin",
        }
        for case, name in cases.items():
            with self.subTest(case=case):
                self.assertIn(name, only)
                self.refused("may not name", lambda s, name=name: s["left_out"].append(f"Rows for {name.title()}."))
                self.refused("may not name", lambda s, name=name: node(s, "grm:antibiotic_vocabulary").update(
                    definition=f"The names GRM resolves antibiotics through, such as {name}."))
        for name in sorted(only):
            with self.subTest(name=name):
                with self.assertRaisesRegex(E.Fault, "may not name"):
                    E.screen({"left_out": [f"Rows for {name}."]}, only)

    def test_a_name_has_its_spellings(self):
        self.assertEqual(E._forms("Trimethoprim/Sulfamethoxazole"),
                         {"trimethoprim/sulfamethoxazole", "trimethoprim", "trimethoprims", "sulfamethoxazole",
                          "sulfamethoxazoles"})
        self.assertEqual(E._forms("fluoroquinolones"), {"fluoroquinolones", "fluoroquinolone"})
        self.assertEqual(E._forms("escherichia_coli__amikacin"), {"escherichia_coli__amikacin",
                                                                  "escherichia coli  amikacin"})
        # A name is refused as a word of its own, and passes inside a longer word.
        E.screen({"left_out": ["An antipenicillinase row."]}, {"penicillin"})
        with self.assertRaisesRegex(E.Fault, "may not name 'penicillin'"):
            E.screen({"left_out": ["A penicillin-binding row."]}, {"penicillin"})
        # A name of three letters or fewer is too short to bar: it would stop the script on ordinary words.
        with tempfile.TemporaryDirectory() as tmp:
            root, source, head = copied(tmp), text(E.VOCABULARY), "const AMBIGUOUS: &[&str] = &["
            self.assertEqual(source.count(head), 1)
            (root / E.VOCABULARY).write_text(source.replace(head, head + '"abc", "abcd", '), encoding="utf-8")
            self.assertEqual(E._vocabulary(root) - E._vocabulary(ROOT), {"abcd", "abcds"})

    def test_the_vocabulary_must_be_readable(self):
        # A table the screen can no longer read stops the script. It never passes with fewer names.
        source = text(E.VOCABULARY)
        edits = [(f"const {table}: ", f"const {table}_OLD: ") for table in E.VOCABULARY_TABLES]
        edits.append((f"fn {E.VOCABULARY_CLASSES}(", f"fn {E.VOCABULARY_CLASSES}_old("))
        for old, new in edits:
            with self.subTest(renamed=old), tempfile.TemporaryDirectory() as tmp:
                root = copied(tmp)
                self.assertEqual(source.count(old), 1)
                (root / E.VOCABULARY).write_text(source.replace(old, new), encoding="utf-8")
                with self.assertRaisesRegex(E.Fault, "holds no name"):
                    E.build_module(SPEC, root)
        with tempfile.TemporaryDirectory() as tmp:
            root = copied(tmp)
            (root / E.VOCABULARY).unlink()
            with self.assertRaisesRegex(E.Fault, "src/drugs.rs: no such committed file"):
                E.build_module(SPEC, root)
            self.assertEqual(E.build_module(SPEC, copied(tempfile.mkdtemp(dir=tmp))), MODULE)

    def test_a_name_is_refused(self):
        cases = [
            ("a species", lambda s: node(s, "grm:dataset").update(definition="One species, such as Escherichia coli.")),
            ("a species, short", lambda s: node(s, "grm:split").update(label="Split of E. coli")),
            ("an antibiotic", lambda s: given(s, "figure-holds-on-new-studies")["measured"][0].update(
                what="datasets on which the set covering machine holds for meropenem")),
            ("a dataset", lambda s: given(s, "figure-holds-on-new-studies")["measured"][4].update(
                note="As on escherichia_coli__amikacin.")),
            ("a drug class", lambda s: s["left_out"].append("Rows for the carbapenem class.")),
            ("a species in a rule", lambda s: s["import_rule"].append("Staphylococcus aureus has its own rule.")),
        ]
        for name, edit in cases:
            with self.subTest(case=name):
                self.refused("may not name", edit)

    def test_a_quote_that_names_a_species_is_refused(self):
        quote = "The other 11 are *K. pneumoniae* genomes that failed the species check"
        self.assertIn(quote, text("docs/outside-study-results.md"))
        self.refused("may not name", lambda s: node(s, "grm:outside_population/same_study")["quotes"].append(
            {"path": "docs/outside-study-results.md", "quote": quote}))

    def test_a_sequence_is_refused(self):
        self.refused("holds no sequence", lambda s: node(s, "grm:duplicate_group").update(
            definition="Genomes that share ACGTACGTACGTACGTACGTACGTACGTACG."))

    def test_a_quote_must_be_in_its_document(self):
        def reword(spec):
            q = node(spec, "grm:dataset")["quotes"][0]
            q["quote"] = q["quote"].replace("one antibiotic", "one drug")

        self.refused("does not hold the sentence", reword)
        self.refused("a quote comes from README.md or a page under docs/",
                     lambda s: node(s, "grm:dataset")["quotes"][0].update(path="src/report.rs"))
        self.refused("at least four words", lambda s: node(s, "grm:dataset")["quotes"][0].update(quote="A dataset"))
        self.refused("quotes at least one sentence", lambda s: node(s, "grm:split").update(quotes=[]))

    def test_a_document_number_must_be_in_the_quote(self):
        self.refused("not a number of the quoted sentence",
                     lambda s: node(s, "grm:state/no_report")["measured"][0].update(count=4))
        self.refused("not a number of the quoted sentence",
                     lambda s: node(s, "grm:state/no_report")["measured"][0].update(of=276))
        self.refused("no suite", lambda s: node(s, "grm:state/no_report")["measured"][0].update(rebuilt_by="tests/none"))

    def test_a_number_is_never_written_in_the_spec(self):
        self.refused("never written in the spec",
                     lambda s: node(s, "grm:state/scored_at_prevalence")["measured"][0].update(count=122))

    def test_a_field_outside_the_list_is_refused(self):
        for at in (["per_species"], ["mic", "with_model"], ["catalog_vs_best", "random", "n"],
                   ["best_learner_counts", "scm"], ["per_learner", "scm", "median_random"],
                   ["error_bars_best_learner", "random", "both"]):
            with self.subTest(field=at):
                self.refused("not a field the export may read", lambda s: node(s, "grm:split/lineage")["measured"][
                    0].update(count={"file": "summary14", "at": at}))
        self.refused("not a field the export may read", lambda s: node(s, "grm:split/random")["measured"][0].update(
            count={"file": "guidance14", "at": ["rows", "catalog:random", "n"]}))
        # the list of datasets with no lineage split may be counted and never read
        self.refused("is not a number", lambda s: node(s, "grm:state/split_unavailable")["measured"][0].update(
            count={"file": "summary14", "at": ["no_lineage"]}))

    def test_a_file_outside_the_list_is_refused(self):
        for path in ("reference/release-summary/rev14-274-20261112/catalog.json",
                     "reference/release-summary/rev14-274-20261112/uncertainty.json",
                     "web/src/data/input-comparison.json", "reference/outside-study/primary/escherichia_coli/metrics.tsv",
                     "reference/release-summary/../../web/src/data/schema.ts"):
            with self.subTest(file=path):
                self.refused("may not read this file", lambda s: s["files"].update(index14=path))
        self.refused("is not named under files", lambda s: node(s, "grm:split/lineage")["measured"][0].update(
            count={"file": "release", "at": ["datasets"]}))
        self.refused("no line reads it", lambda s: s["files"].update(
            more="reference/release-summary/rev14-274-20261112/guidance.json"))

    def test_an_unknown_field_is_refused(self):
        cases = [
            lambda s: s.update(samples=[]),
            lambda s: node(s, "grm:dataset").update(features=["k-mer"]),
            lambda s: node(s, "grm:dataset")["measured"][0].update(dataset="one"),
            lambda s: given(s, "figure-holds-on-new-studies").update(models=[]),
            lambda s: s["edges"][0].update(relation="causes"),
        ]
        for i, edit in enumerate(cases):
            with self.subTest(case=i):
                self.refused("the export has no field", edit)

    def test_a_count_must_fit_its_denominator_and_its_revision(self):
        def swap(spec):
            x = node(spec, "grm:state/scored_at_prevalence")["measured"][0]
            x["count"], x["of"] = x["of"], x["count"]

        self.refused("not a count within its denominator", swap)
        self.refused("says why in note", lambda s: node(s, "grm:outside_population/primary")["measured"][2].pop("note"))
        self.refused("and the line says 14", lambda s: node(s, "grm:state/scored_at_prevalence")["measured"][0].update(
            count={"file": "index12", "at": ["unflagged_with_no_call"]}))
        self.refused("and the line says 14", lambda s: given(s, "figure-holds-on-new-studies")["measured"][0].update(
            revision=14))
        self.refused("and the line says 12", lambda s: given(s, "learners-scored-on-same-genomes")["measured"][7][
            "values"][0].update(value={"file": "summary14", "at": ["per_learner", "gbt", "median_lineage"]}))

    def test_the_graph_must_close(self):
        self.refused("both ends of an edge", lambda s: s["edges"][0].update(to="grm:model"))
        self.refused("unknown edge type", lambda s: s["edges"][0].update(type="causes"))
        self.refused("written twice", lambda s: s["edges"].append(dict(s["edges"][0])))
        self.refused("sits under at least one node", lambda s: given(s, "figure-holds-on-new-studies").update(
            under=["grm:model"]))
        self.refused("defined once", lambda s: s["nodes"].append(dict(s["nodes"][0])))
        self.refused("unknown kind", lambda s: node(s, "grm:dataset").update(kind="gene"))
        self.refused("unknown truth state", lambda s: node(s, "grm:dataset").update(truth_state="proven"))
        self.refused("typed_export", lambda s: s.update(binding="provisional_snapshot"))

    def test_the_command_reports_a_changed_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            other = pathlib.Path(tmp) / "grm.json"
            other.write_text(FILE.read_text(encoding="utf-8").replace('"count": 32,', '"count": 33,'), encoding="utf-8")
            self.assertNotEqual(other.read_text(encoding="utf-8"), FILE.read_text(encoding="utf-8"))
            with contextlib.redirect_stderr(io.StringIO()) as said:
                self.assertEqual(E.main(["export.py", "--check", str(other)]), 1)
            self.assertIn("the file differs", said.getvalue())

    def test_the_command_writes_and_checks_both_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = copied(tmp)
            export, sums = root / E.EXPORT, root / E.SUMS
            with contextlib.redirect_stderr(io.StringIO()) as said:
                self.assertEqual(E.main(["export.py", "--check"], root), 1)
            self.assertEqual(said.getvalue().count("the file differs"), 2)
            self.assertEqual(E.main(["export.py", "--write"], root), 0)
            self.assertEqual((export.read_bytes(), sums.read_bytes()), (FILE.read_bytes(), SUMS.read_bytes()))
            self.assertEqual(E.main(["export.py", "--check"], root), 0)
            for path, other in ((sums, sums.read_text(encoding="utf-8").replace("  ", " ")),
                                (export, export.read_text(encoding="utf-8").replace('"count": 32,', '"count": 33,'))):
                with self.subTest(changed=path.name):
                    kept = path.read_bytes()
                    path.write_text(other, encoding="utf-8")
                    self.assertNotEqual(path.read_bytes(), kept)
                    with contextlib.redirect_stderr(io.StringIO()) as said:
                        self.assertEqual(E.main(["export.py", "--check"], root), 1)
                    self.assertEqual(said.getvalue(), f"{path}: the file differs from what the committed sources give\n")
                    path.write_bytes(kept)
            with self.assertRaises(SystemExit):
                E.main(["export.py", "--write", str(export)], root)


class Wiring(unittest.TestCase):
    def test_the_suite_is_in_ci_and_its_document(self):
        self.assertIn("          - suite: tests/givens_export\n            top: .\n", text(".github/workflows/ci.yml"))
        tests = sum(len(unittest.defaultTestLoader.getTestCaseNames(case))
                    for case in (Committed, Shape, Numbers, Missing, Limits, Wiring))
        self.assertIn(f"`tests/givens_export` {tests}", text("docs/ci.md"))

    def test_the_readme_and_the_page_name_the_file(self):
        readme = text("README.md")
        part = readme[readme.index("\n## Givens\n"):]
        part = part[: part.index("\n## ", 4)]
        for name in ("reference/givens-export/grm.json", "scripts/givens-export/export.py", "docs/givens-export.md"):
            self.assertIn(name, part)
        self.assertIn("GRM reads nothing from Givens", part)
        page = text("docs/givens-export.md")
        documents = sum(x["basis"] == "document" for _, x in LINES)
        unmeasured = sum("unmeasured" in x for x in ENTRIES.values())
        for sentence in (
            f"It holds {len(MODULE['nodes'])} definitions, {len(MODULE['edges'])} links between them and "
            f"{len(MODULE['assumptions'])} assumptions.",
            f"Its {len(LINES)} measured lines split into {len(LINES) - documents} read from a file and {documents} "
            f"quoted from a document.",
            f"GRM has no count for {unmeasured} entries, and they say so.",
        ):
            with self.subTest(sentence=sentence):
                self.assertIn(sentence, page)
        self.assertEqual(sorted(re.findall(r"`(reference/release-summary/[^`]+|web/src/data/[^`]+)`", page)),
                         [f["path"] for f in MODULE["sources"]["files"]])

    def test_every_count_of_the_page_is_a_line_of_the_file(self):
        page = text("docs/givens-export.md")
        stated = {(int(a), int(b)) for a, b in re.findall(r"(\d+) of (\d+)", page)}
        documents = {(x["count"], x["of"]) for _, x in LINES if x["basis"] == "document"}
        self.assertEqual(stated, documents)
        self.assertLessEqual(E._numbers(re.sub(r"`[^`]*`|\d+ of \d+|\(https[^)]*\)", "", page)),
                             {"22", "24", "3", "36", "28", "8", "2", "9", "12", "14", "1", "3.12", "2026", "256"})


if __name__ == "__main__":
    unittest.main()
