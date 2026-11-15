"""collect.py list on a synthetic project: the rows, the classes, and what stops it."""

import os
import subprocess
import sys
import unittest

from tests.collect import kit


class TheList(unittest.TestCase):
    def test_every_file_has_its_size_sha256_and_class(self):
        d, proj, rules = kit.project(self)
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual((rc, err), (0, kit.MARKED))
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text())
        # 21 + 20 public; 21 + 0 + 16 + 10 internal; 9 input; 3 + 194 + 11 + 3 regenerable
        self.assertEqual(out, "class\tentries\tbytes\npublic-result\t2\t41\ninternal-record\t4\t47\n"
                              "public-input\t1\t9\nregenerable\t4\t211\n"
                              f"list: wrote 11 rows to {d / 'list.tsv'}\n")

    def test_the_list_does_not_depend_on_the_threads_or_on_where_the_project_lies(self):
        d, proj, rules = kit.project(self)
        other = kit.tmp(self)
        kit.write(other / "elsewhere", kit.FILES)
        texts = []
        for threads, top in ((1, proj), (8, proj), (3, other / "elsewhere")):
            rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={top}", "--out", d / "l.tsv",
                                 "--threads", threads)
            self.assertEqual((rc, err), (0, kit.MARKED))
            texts.append((d / "l.tsv").read_bytes())
        self.assertEqual(texts[0], texts[1])
        self.assertEqual(texts[0], texts[2])

    def test_a_file_over_one_block_is_hashed_to_its_end(self):
        d = kit.tmp(self)
        data = kit.noise(1, (1 << 20) + 4097)  # one block of the reader and a little more
        kit.write(d / "proj", {"big.bin": data})
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual((d / "l.tsv").read_text(encoding="utf-8"),
                         kit.list_text([("proj/big.bin", 1052673, data, "internal-record")]))

    def test_nothing_under_a_regenerable_directory_is_opened(self):
        d, proj, rules = kit.project(self)
        for rel in ("build/target/debug/app", "build/__pycache__/m.pyc", "env/lib/site.py"):
            os.chmod(proj / rel, 0)
        os.symlink("/nonexistent", proj / "build" / "target" / "dangling")
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual((rc, err), (0, kit.MARKED))
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text())

    def test_a_directory_is_regenerable_by_a_marker_file_and_not_by_a_name_like_target(self):
        """Each directory a marker made regenerable is named on stderr with its marker and its count of files: a
        marker can sit in a folder of records too, and a reader must be able to notice it."""
        d = kit.tmp(self)
        kit.write(d / "proj", {
            "a/target/plan.md": b"what we aim at",  # 14 bytes: a folder that is only named target
            "a/.cargo/config.toml": b"[build]\n",  # 8: a tracked setting, not a cache
            "b/target/CACHEDIR.TAG": b"Signature: 8a477f597d28d172789f06886806bc55\n", "b/target/x.o": b"o" * 6,
            "c/np/pyvenv.cfg": b"home = /usr/bin\n", "c/np/lib/site.py": b"s" * 5,  # 16 + 5
            "d/.pytest_cache/CACHEDIR.TAG": b"Signature: 8a477f597d28d172789f06886806bc55\n",  # 44
        })
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, "list: 3 directories are regenerable for a marker file they hold, so their "
                                        "files go into no archive:\n"
                                        "  proj/b/target/  CACHEDIR.TAG  2 files\n"
                                        "  proj/c/np/  pyvenv.cfg  2 files\n"
                                        "  proj/d/.pytest_cache/  CACHEDIR.TAG  1 files\n"))
        self.assertEqual((d / "l.tsv").read_text(encoding="utf-8"), kit.list_text([
            ("proj/a/.cargo/config.toml", 8, b"[build]\n", "internal-record"),
            ("proj/a/target/plan.md", 14, b"what we aim at", "internal-record"),
            ("proj/b/target/", 50, None, "regenerable"),
            ("proj/c/np/", 21, None, "regenerable"),
            ("proj/d/.pytest_cache/", 44, None, "regenerable"),
        ]))

    def test_records_beside_a_marker_are_named_with_the_marker_and_their_count(self):
        """A marker of any kind (an empty file, a directory, a dangling link) makes its directory regenerable, so
        the records beside it go into no archive; list says so on stderr."""
        d = kit.tmp(self)
        kit.write(d / "proj", {"analysis/results.tsv": b"r\t1\n", "analysis/plan.md": b"plan\n",
                               "analysis/pyvenv.cfg": b"", "runs/r1/out.json": b"{}\n", "notes/n.md": b"n\n"})
        (d / "proj" / "runs" / "r1" / "CACHEDIR.TAG").mkdir()
        os.symlink(d / "nowhere", d / "proj" / "notes" / "pyvenv.cfg")
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, out, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, "list: 3 directories are regenerable for a marker file they hold, so their "
                                        "files go into no archive:\n"
                                        "  proj/analysis/  pyvenv.cfg  3 files\n"
                                        "  proj/notes/  pyvenv.cfg  1 files\n"
                                        "  proj/runs/r1/  CACHEDIR.TAG  1 files\n"), out)

    def test_the_first_matching_rule_gives_the_class(self):
        d, proj, rules = kit.project(self)
        rules.write_text("internal-record\tproj/reports/beta.json\n" + kit.RULES, encoding="utf-8")
        rc, _, _ = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 0)
        rows = [(p, n, data, "internal-record" if p == "proj/reports/beta.json" else cls)
                for p, n, data, cls in kit.ROWS]
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text(rows))

    def test_a_file_no_rule_matches_stops_the_list_and_is_named(self):
        d, proj, rules = kit.project(self)
        rules.write_text(kit.RULES.replace("internal-record\t*\n", "internal-record\tproj/notes/*.md\n"),
                         encoding="utf-8")
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("gives a class to 3 files", err)
        for name in ("proj/notes/alpha-copy.json", "proj/notes/empty.txt", "proj/notes/messages/m1.txt"):
            self.assertIn(f"  {name}\n", err)
        self.assertFalse((d / "list.tsv").exists())

    def test_a_root_may_be_one_file_under_a_nested_label(self):
        d, proj, rules = kit.project(self)
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj/notes/ledger.md={proj / 'notes' / 'ledger.md'}",
                             "--root", f"proj/reports={proj / 'reports'}", "--out", d / "list.tsv")
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"),
                         kit.list_text([r for r in kit.ROWS if r[0] == "proj/notes/ledger.md" or "/reports/" in r[0]]))

    def test_two_roots_that_give_one_path_stop_the_list(self):
        d, proj, rules = kit.project(self)
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}",
                             "--root", f"proj/notes={proj / 'notes'}", "--out", d / "list.tsv")
        self.assertEqual(rc, 2)
        self.assertIn("4 paths are listed by two sources", err)
        self.assertFalse((d / "list.tsv").exists())

    def test_a_link_stops_the_list_unless_it_is_left_out(self):
        d, proj, rules = kit.project(self)
        os.symlink("ledger.md", proj / "notes" / "link.md")
        os.mkfifo(proj / "notes" / "pipe")
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 2)
        self.assertIn("2 links or special files", err)
        self.assertIn("  proj/notes/link.md\n  proj/notes/pipe\n", err)
        self.assertFalse((d / "list.tsv").exists())
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv",
                             "--skip-special")
        self.assertEqual(rc, 0)
        self.assertIn("left out 2 links or special files", err)
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text())

    def test_a_name_with_a_control_character_stops_the_list(self):
        d, proj, rules = kit.project(self)
        (proj / "notes" / "two\tcells.md").write_bytes(b"x")
        (proj / "notes" / "two\nlines.md").write_bytes(b"x")
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 2)
        self.assertIn("2 names a list cannot hold", err)
        self.assertIn(r"'proj/notes/two\tcells.md'", err)
        self.assertFalse((d / "list.tsv").exists())

    def test_rules_a_list_cannot_use(self):
        d, proj, rules = kit.project(self)
        for text, said in (("internal\t*\n", "'internal' is not a class a rule can give"),
                           ("refused\t*\n", "'refused' is not a class a rule can give"),
                           ("internal-record *\n", "want CLASS<TAB>PATTERN"),
                           ("internal-record\t\n", "want CLASS<TAB>PATTERN"),
                           ("internal-record\tproj/notes/\n", "only regenerable takes one")):
            with self.subTest(text=text):
                rules.write_text(text, encoding="utf-8")
                rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
                self.assertEqual((rc, out), (2, ""))
                self.assertIn(f"{rules}:1: ", err)
                self.assertIn(said, err)
                self.assertFalse((d / "list.tsv").exists())

    def test_a_root_that_is_not_there_and_no_root_at_all(self):
        d, proj, rules = kit.project(self)
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={d / 'absent'}", "--out", d / "list.tsv")
        self.assertEqual(rc, 2)
        self.assertIn("cannot read the root", err)
        rc, _, err = kit.run("list", "--rules", rules, "--out", d / "list.tsv")
        self.assertEqual(rc, 2)
        self.assertIn("give at least one --root", err)
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"../up={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 2)
        self.assertIn("want NAME=VALUE with a plain name", err)


class BarredWords(unittest.TestCase):
    """A path holding a barred word is named, never opened and never entered, and nothing is written."""

    def test_each_word_in_any_case_stops_the_list(self):
        for word in kit.BARRED:
            for spelt in (word, word.upper(), word.capitalize()):
                with self.subTest(word=spelt):
                    d, proj, rules = kit.project(self)
                    kit.write(proj, {f"notes/{spelt}-set/table.tsv": b"x", f"notes/list-{spelt}.txt": b"y"})
                    os.chmod(proj / "notes" / f"{spelt}-set", 0)  # entering it would fail with status 2
                    os.chmod(proj / "notes" / f"list-{spelt}.txt", 0)
                    rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "l.tsv")
                    self.assertEqual((rc, out), (3, ""))
                    self.assertIn("refused: 2 paths hold a barred word; none was opened:", err)
                    self.assertIn(f"  proj/notes/{spelt}-set  ({word})\n", err)
                    self.assertIn(f"  proj/notes/list-{spelt}.txt  ({word})\n", err)
                    self.assertIn("nothing was written", err)
                    self.assertFalse((d / "l.tsv").exists())

    def test_a_barred_word_in_another_form_or_split_by_a_separator_stops_the_list(self):
        refused = {"\uff2b\uff33\uff24-notes.txt": "ksd",  # full width
                   "\u212asd.txt": "ksd",  # the Kelvin sign
                   "\u017fealed.txt": "sealed",  # the long s
                   "Held_Out.tsv": "held-out", "held out.tsv": "held-out", "held.out.tsv": "held-out",
                   "heldout.tsv": "held-out", "KS-30.tsv": "ks30", "ks_50.tsv": "ks50", "ks - 30.tsv": "ks30",
                   "Sig-Drug.tsv": "sigdrug", "sig__drug.tsv": "sigdrug", "\uff33\uff29\uff27drug.tsv": "sigdrug",
                   "ks/50.tsv": "ks50", "held/out/x.tsv": "held-out"}
        kept = ("tasks/data.tsv", "ks-3-0.tsv", "sealing.txt", "sig-dru-g.tsv", "k-s-d.txt")
        d, proj, rules = kit.project(self)
        kit.write(proj / "notes", {name: b"x" for name in list(refused) + list(kept)})
        for name in refused:
            os.chmod(proj / "notes" / name, 0)  # opening it would fail with status 2
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "l.tsv")
        self.assertEqual((rc, out), (3, ""), err)
        self.assertIn(f"refused: {len(refused)} paths hold a barred word; none was opened:", err)
        for name, word in refused.items():
            shown = "proj/notes/held/out" if name.startswith("held/") else f"proj/notes/{name}"
            self.assertIn(f"  {shown}  ({word})\n", err)
        for name in kept:
            self.assertNotIn(name, err)
        for name in refused:
            os.chmod(proj / "notes" / name, 0o600)
            os.unlink(proj / "notes" / name)
        os.removedirs(proj / "notes" / "held" / "out")
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, kit.MARKED), out)
        listed = (d / "l.tsv").read_text(encoding="utf-8")
        for name in kept:
            self.assertIn(f"proj/notes/{name}\t1\t{kit.sha(b'x')}\tinternal-record\n", listed)

    def test_a_barred_name_under_a_regenerable_directory_stops_the_list_too(self):
        d, proj, rules = kit.project(self)
        kit.write(proj, {"build/target/debug/private/key.bin": b"k"})
        os.chmod(proj / "build" / "target" / "debug" / "private", 0)
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 3)
        self.assertIn("  proj/build/target/debug/private  (private)\n", err)
        self.assertFalse((d / "list.tsv").exists())

    def test_record_refused_writes_the_name_alone_and_goes_on(self):
        d, proj, rules = kit.project(self)
        kit.write(proj, {"notes/ks30/table.tsv": b"x", "build/target/vumc.o": b"12345"})
        os.chmod(proj / "notes" / "ks30", 0)
        os.chmod(proj / "build" / "target" / "vumc.o", 0)
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv",
                               "--record-refused")
        self.assertEqual(rc, 0)
        self.assertIn("2 paths hold a barred word; none was opened", err)
        rows = kit.ROWS + [("proj/build/target/vumc.o", "-", None, "refused"),
                           ("proj/notes/ks30", "-", None, "refused")]
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text(sorted(rows, key=lambda r: r[0])))
        self.assertIn("refused\t2\t0\n", out)  # the 5 bytes of vumc.o are not counted in target/ either

    def test_refused_paths_are_counted_by_top_folder_with_and_without_record_refused(self):
        d, proj, rules = kit.project(self)
        kit.write(d / "other", {"ksd-a.txt": b"a", "b/private.txt": b"b", "c.txt": b"c"})
        kit.write(proj, {"notes/ks30.txt": b"x"})
        for more in ((), ("--record-refused",)):
            with self.subTest(more=more):
                rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}",
                                       "--root", f"other={d / 'other'}", "--out", d / "l.tsv", *more)
                self.assertEqual(rc, 0 if more else 3, err)
                self.assertIn("3 paths hold a barred word; none was opened:\n", err)
                self.assertIn("list: refused paths by top folder; a refused path is in no archive, so look at it "
                              "by hand:\n  other  2\n  proj  1\n", err)

    def test_a_list_written_with_record_refused_is_taken_by_its_own_pack_and_summary(self):
        d, proj, rules = kit.project(self)
        kit.write(proj, {"notes/ks30/table.tsv": b"x", "notes/Sealed.md": b"y"})
        os.chmod(proj / "notes" / "ks30", 0)
        os.chmod(proj / "notes" / "Sealed.md", 0)
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv",
                               "--record-refused")
        self.assertEqual(rc, 0, err)
        rc, out, err = kit.run("pack", "--list", d / "list.tsv", "--class", "internal-record", "--name", "demo",
                               "--out", d / "out", "--root", f"proj={proj}")
        self.assertEqual((rc, err), (0, ""), out)
        rc, out, err = kit.run("check", "--parts", d / "out" / "demo.internal-record.parts.sha256")
        self.assertEqual((rc, err), (0, ""), out)
        rc, out, err = kit.run("summary", "--place", f"here={d / 'list.tsv'}")
        self.assertEqual((rc, err), (0, ""), out)
        self.assertIn("refused", out)

    def test_a_refused_path_whose_list_path_holds_no_word_stops_record_refused_too(self):
        d, proj, rules = kit.project(self)
        kit.write(d / "sigdrug-data", {"table.tsv": b"x"})
        os.symlink(d / "sigdrug-data", proj / "notes" / "tables")
        table = d / "table.tsv"
        table.write_text(f"path\tbytes\tsha256\nx/held-out/proj/a.tsv\t1\t{kit.sha(b'a')}\n", encoding="utf-8")
        cases = (("a link with a clean name", ("--root", f"proj={proj}"), "proj/notes/tables"),
                 ("a root with a clean label", ("--root", f"clean={d / 'sigdrug-data'}"), "clean"),
                 ("a word --strip-prefix cuts off", ("--listing", table, "--strip-prefix", "x/held-out/"),
                  "proj/a.tsv"))
        for name, args, path in cases:
            with self.subTest(name):
                rc, out, err = kit.run("list", "--rules", rules, *args, "--out", d / "l.tsv", "--record-refused")
                self.assertEqual((rc, out), (3, ""), err)
                self.assertIn("list: 1 of them hold no barred word in their list path (their place on disk does, "
                              "through a link or the root, or --strip-prefix cut it off), and a refused row must "
                              f"hold its word; nothing was written:\n  {path}\n", err)
                self.assertFalse((d / "l.tsv").exists())

    def test_a_root_or_a_link_that_leads_to_a_barred_place_is_refused(self):
        d, proj, rules = kit.project(self)
        kit.write(d / "sigdrug-data", {"table.tsv": b"x"})
        with self.subTest("the root's own place"):
            rc, _, err = kit.run("list", "--rules", rules, "--root", f"clean={d / 'sigdrug-data'}",
                                 "--out", d / "l.tsv")
            self.assertEqual(rc, 3)
            self.assertIn("  clean  (sigdrug)\n", err)
        with self.subTest("a link with a clean name"):
            os.symlink(d / "sigdrug-data", proj / "notes" / "tables")
            rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "l.tsv",
                                 "--skip-special")
            self.assertEqual(rc, 3)
            self.assertIn("  proj/notes/tables  (sigdrug)\n", err)
        self.assertFalse((d / "l.tsv").exists())

    def test_a_barred_word_in_an_argument_stops_every_command(self):
        d, proj, rules = kit.project(self)
        for args in (("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "ksd" / "list.tsv"),
                     ("list", "--rules", d / "Private" / "rules.tsv", "--root", f"proj={proj}", "--out", d / "l.tsv"),
                     ("list", "--rules", rules, "--listing", d / "held-out.tsv", "--out", d / "l.tsv"),
                     ("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "x", "--out", d / "vumc"),
                     ("pack", "--list", d / "ks50.tsv", "--class", "internal-record", "--name", "x", "--out", d / "o"),
                     ("check", d / "clinical" / "x.internal-record"),
                     ("check", "--parts", d / "sealed" / "x.internal-record.parts.sha256"),
                     ("summary", "--place", f"here={d / 'sigdrug.tsv'}"),
                     ("summary", "--place", f"here={d / 'l.tsv'}", "--hashes", f"rel={d / 'ks30.sha256'}"),
                     ("pack", "--list", d / "l.tsv", "--class", "public-result", "--name", "x", "--out", d / "o",
                      "--published", d / "ks50" / "release.sha256"),
                     ("summary", "--place", f"here={d / 'l.tsv'}", "--out", d / "vumc" / "all.tsv")):
            with self.subTest(args=args[:2]):
                rc, out, err = kit.run(*args)
                self.assertEqual((rc, out), (3, ""))
                self.assertIn("refused: the argument", err)
        rc, _, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "KSD-set",
                             "--out", d / "o")
        self.assertEqual(rc, 3)
        self.assertIn("the name KSD-set holds the barred word 'ksd'", err)
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["proj", "rules.tsv"])

    def test_a_barred_label_on_a_clean_place_stops_the_list(self):
        d = kit.tmp(self)
        kit.write(d / "proj", {"notes/a.md": b"synthetic\n"})
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, out, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"KS50={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, out), (3, ""), err)
        self.assertIn("refused: 1 paths hold a barred word; none was opened:\n  KS50  (ks50)\n", err)
        self.assertFalse((d / "l.tsv").exists())

    def test_a_barred_word_before_dot_dot_stops_every_command(self):
        """The path is judged as given: d/sealed/../l.tsv names d/l.tsv, but it holds a barred word. None of the
        barred folders below exists, and none is made."""
        d, proj, lst = kit.listed(self)
        published = kit.published(d, proj)
        rules = d / "rules.tsv"
        cut = d / "out" / "demo.internal-record"
        sums = d / "out" / "demo.internal-record.parts.sha256"

        def via(word, path):
            return f"{d}/{word}/../{path.relative_to(d)}"
        pack = ("pack", "--list", lst, "--class", "public-result", "--name", "demo", "--root", f"proj={proj}",
                "--published", published, "--out", d / "out")
        arguments = (
            (("list", "--rules", via("ks50", rules), "--root", f"proj={proj}", "--out", d / "l.tsv"), "ks50"),
            (("list", "--rules", rules, "--root", f"proj={proj}", "--out", via("Sealed", d / "l.tsv")), "sealed"),
            (("list", "--rules", rules, "--listing", via("vumc", lst), "--out", d / "l.tsv"), "vumc"),
            (("pack", "--list", via("KS30", lst)) + pack[3:], "ks30"),
            (pack[:-1] + (via("clinical", d / "out"),), "clinical"),
            (pack[:9] + ("--published", via("sigdrug", published), "--out", d / "out"), "sigdrug"),
            (("check", via("private", cut)), "private"),
            (("check", "--parts", via("held-out", sums)), "held-out"),
            (("summary", "--place", f"here={via('ksd', lst)}"), "ksd"),
            (("summary", "--place", f"here={lst}", "--hashes", f"rel={via('ks50', published)}"), "ks50"),
            (("summary", "--place", f"here={lst}", "--out", via("sealed", d / "s.tsv")), "sealed"),
        )
        for args, word in arguments:
            with self.subTest(args=args[:3], word=word):
                rc, out, err = kit.run(*args)
                self.assertEqual((rc, out), (3, ""), err)
                self.assertIn(f"refused: the argument {d}/", err)
                self.assertIn(f"holds the barred word {word!r}", err)
        with self.subTest("a root of list"):
            rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={via('held-out', proj)}",
                                   "--out", d / "l.tsv")
            self.assertEqual((rc, out), (3, ""), err)
            self.assertIn("refused: 1 paths hold a barred word; none was opened:\n  proj  (held-out)\n", err)
        for root, word in ((f"proj={via('ksd', proj)}", "ksd"), (f"spare={via('Private', proj)}", "private")):
            with self.subTest("a root of pack", root=root):
                rc, out, err = kit.run(*pack[:7], "--root", root, *pack[9:])
                self.assertEqual((rc, out), (3, ""), err)
                self.assertIn(f"pack: refused: the root {root} holds the barred word {word!r}", err)
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["list.tsv", "proj", "release.sha256", "rules.tsv"])


class TheFileSystemEncoding(unittest.TestCase):
    """Names are read through the file system encoding, so the tool runs only where it is UTF-8."""

    # Runs the tool with sys.getfilesystemencoding() giving another encoding; the locale of this machine is left alone.
    LATIN1 = ("import runpy, sys; sys.getfilesystemencoding = lambda: 'iso8859-1'; sys.argv = sys.argv[1:]; "
              "runpy.run_path(sys.argv[0], run_name='__main__')")

    def test_another_encoding_stops_every_command_before_it_reads_anything(self):
        d, proj, rules = kit.project(self)
        for args in (("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "l.tsv"),
                     ("check", d / "x.internal-record"), ("summary", "--place", f"here={d / 'l.tsv'}")):
            with self.subTest(command=args[0]):
                p = subprocess.run([sys.executable, "-I", "-B", "-W", "error", "-c", self.LATIN1, str(kit.TOOL),
                                    *map(str, args)], capture_output=True, text=True, stdin=subprocess.DEVNULL)
                self.assertEqual((p.returncode, p.stdout), (2, ""), p.stderr)
                self.assertEqual(p.stderr, "collect: the file system encoding is iso8859-1, not UTF-8, so names would "
                                           "not be read the same on every machine; nothing was done. Run it under a "
                                           "UTF-8 locale (LC_ALL=C.UTF-8) or as python3 -I -X utf8 "
                                           "scripts/collect/collect.py\n")
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["proj", "rules.tsv"])

    def test_utf_8_under_another_name_runs(self):
        d, proj, rules = kit.project(self)
        code = self.LATIN1.replace("'iso8859-1'", "'UTF8'")
        p = subprocess.run([sys.executable, "-I", "-B", "-W", "error", "-c", code, str(kit.TOOL), "list", "--rules",
                            str(rules), "--root", f"proj={proj}", "--out", str(d / "l.tsv")],
                           capture_output=True, text=True, stdin=subprocess.DEVNULL)
        self.assertEqual((p.returncode, p.stderr), (0, kit.MARKED), p.stdout)


class ATableInPlaceOfAWalk(unittest.TestCase):
    """--listing: the rows of a table made elsewhere, classed by the same rules. Only the table is opened."""

    HEADER = "sha256\tbytes\tmtime\tpath\n"

    def table(self, d, rows):
        text = self.HEADER + "".join(f"{kit.sha(data)}\t{len(data)}\t0\t/far/home/{p}\n" for p, data in rows)
        (d / "far.tsv").write_text(text, encoding="utf-8")
        return d / "far.tsv"

    def test_rows_take_their_class_from_the_rules(self):
        d, proj, rules = kit.project(self)
        rows = [("proj/reports/alpha.json", b"A"), ("proj/notes/n.md", b"note"),
                ("proj/web/node_modules/pkg/index.js", b"js"), ("proj/env/lib/x.py", b"py"), ("proj/a.tmp", b"t"),
                ("proj/b/target/debug/x.o", b"obj"), ("proj/b/target/CACHEDIR.TAG", b"Signature"),
                ("proj/c/target/plan.md", b"plan")]
        rc, out, err = kit.run("list", "--rules", rules, "--listing", self.table(d, rows), "--strip-prefix",
                               "/far/home/", "--out", d / "list.tsv")
        self.assertEqual((rc, err), (0, "list: 1 directories are regenerable for a marker file they hold, so their "
                                        "files go into no archive:\n  proj/b/target/  CACHEDIR.TAG  2 files\n"))
        want = [("proj/a.tmp", 1, b"t", "regenerable"),
                ("proj/b/target/CACHEDIR.TAG", 9, b"Signature", "regenerable"),  # the marker's own directory
                ("proj/b/target/debug/x.o", 3, b"obj", "regenerable"),
                ("proj/c/target/plan.md", 4, b"plan", "internal-record"),  # no marker: only named target
                ("proj/env/lib/x.py", 2, b"py", "regenerable"),
                ("proj/notes/n.md", 4, b"note", "internal-record"),
                ("proj/reports/alpha.json", 1, b"A", "public-result"),
                ("proj/web/node_modules/pkg/index.js", 2, b"js", "regenerable")]
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text(want))

    def test_a_barred_row_stops_the_list_or_is_recorded_by_name(self):
        d, proj, rules = kit.project(self)
        table = self.table(d, [("proj/notes/n.md", b"note"), ("proj/Held-Out/x.tsv", b"x")])
        rc, _, err = kit.run("list", "--rules", rules, "--listing", table, "--strip-prefix", "/far/home/",
                             "--out", d / "list.tsv")
        self.assertEqual(rc, 3)
        self.assertIn("  proj/Held-Out/x.tsv  (held-out)\n", err)
        self.assertFalse((d / "list.tsv").exists())
        rc, _, _ = kit.run("list", "--rules", rules, "--listing", table, "--strip-prefix", "/far/home/",
                           "--out", d / "list.tsv", "--record-refused")
        self.assertEqual(rc, 0)
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"),
                         "path\tbytes\tsha256\tclass\nproj/Held-Out/x.tsv\t-\t-\trefused\n"
                         f"proj/notes/n.md\t4\t{kit.sha(b'note')}\tinternal-record\n")

    def test_a_word_in_the_part_that_is_stripped_refuses_the_row(self):
        d, proj, rules = kit.project(self)
        (d / "far.tsv").write_text(self.HEADER + f"{kit.sha(b'x')}\t1\t0\t/far/ksd/proj/notes/n.md\n", encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", rules, "--listing", d / "far.tsv", "--strip-prefix", "/far/ksd/",
                             "--out", d / "list.tsv")
        self.assertEqual(rc, 3)
        self.assertIn("  proj/notes/n.md  (ksd)\n", err)

    def test_a_table_that_is_not_one_stops_the_list(self):
        d, proj, rules = kit.project(self)
        good = f"{kit.sha(b'x')}\t1\t0\t/far/home/proj/notes/n.md\n"
        for text, said in (("sha256\tbytes\tmtime\n" + good, "the first line must name path, bytes and sha256"),
                           (self.HEADER + good.rstrip("\n"), "the last line must end"),
                           (self.HEADER + good.replace("\t1\t", "\t1\t1\t"), "5 cells, want 4"),
                           (self.HEADER + good.replace("\t1\t", "\tone\t"), "a size and a sha256"),
                           (self.HEADER + good.replace(kit.sha(b"x"), "0" * 63), "a size and a sha256"),
                           (self.HEADER + good.replace("/far/home/", "/near/home/"), "want a path under"),
                           (self.HEADER + good.replace("proj/notes", "proj/../notes"), "names a list cannot hold")):
            with self.subTest(said=said):
                (d / "far.tsv").write_text(text, encoding="utf-8")
                rc, out, err = kit.run("list", "--rules", rules, "--listing", d / "far.tsv", "--strip-prefix",
                                       "/far/home/", "--out", d / "list.tsv")
                self.assertEqual((rc, out), (2, ""))
                self.assertIn(said, err)
                self.assertFalse((d / "list.tsv").exists())

    def test_a_table_and_a_walk_make_one_list(self):
        d, proj, rules = kit.project(self)
        table = self.table(d, [("far/notes/n.md", b"note")])
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--listing", table,
                             "--strip-prefix", "/far/home/", "--out", d / "list.tsv")
        self.assertEqual((rc, err), (0, kit.MARKED))
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"),
                         kit.list_text([("far/notes/n.md", 4, b"note", "internal-record")] + kit.ROWS))


if __name__ == "__main__":
    unittest.main()
