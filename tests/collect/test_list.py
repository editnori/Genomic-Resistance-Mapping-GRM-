"""collect.py list on a synthetic project: the rows, the classes, and what stops it."""

import os
import unittest

from tests.collect import kit


class TheList(unittest.TestCase):
    def test_every_file_has_its_size_sha256_and_class(self):
        d, proj, rules = kit.project(self)
        rc, out, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual((rc, err), (0, ""))
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
            self.assertEqual((rc, err), (0, ""))
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
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"), kit.list_text())

    def test_a_directory_is_regenerable_by_a_marker_file_and_not_by_a_name_like_target(self):
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
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual((d / "l.tsv").read_text(encoding="utf-8"), kit.list_text([
            ("proj/a/.cargo/config.toml", 8, b"[build]\n", "internal-record"),
            ("proj/a/target/plan.md", 14, b"what we aim at", "internal-record"),
            ("proj/b/target/", 50, None, "regenerable"),
            ("proj/c/np/", 21, None, "regenerable"),
            ("proj/d/.pytest_cache/", 44, None, "regenerable"),
        ]))

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
                     ("summary", "--place", f"here={d / 'l.tsv'}", "--hashes", f"rel={d / 'ks30.sha256'}")):
            with self.subTest(args=args[:2]):
                rc, out, err = kit.run(*args)
                self.assertEqual((rc, out), (3, ""))
                self.assertIn("refused: the argument", err)
        rc, _, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "KSD-set",
                             "--out", d / "o")
        self.assertEqual(rc, 3)
        self.assertIn("the name KSD-set holds the barred word 'ksd'", err)
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["proj", "rules.tsv"])


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
        self.assertEqual((rc, err), (0, ""))
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
        self.assertEqual((rc, err), (0, ""))
        self.assertEqual((d / "list.tsv").read_text(encoding="utf-8"),
                         kit.list_text([("far/notes/n.md", 4, b"note", "internal-record")] + kit.ROWS))


if __name__ == "__main__":
    unittest.main()
