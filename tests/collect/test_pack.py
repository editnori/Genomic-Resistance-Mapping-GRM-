"""collect.py pack on a synthetic project: one class per archive, parts under the size limit, and what stops it."""

import os
import unittest

from tests.collect import kit

INNER = "path\tbytes\tsha256\tclass\tpart\n"


class OneClassPerArchive(unittest.TestCase):
    def test_an_archive_holds_the_files_of_its_class_and_their_list(self):
        for compress, ext in (("gz", ".tar.gz"), ("xz", ".tar.xz"), ("none", ".tar")):
            with self.subTest(compress=compress):
                d, proj, sums = kit.packed(self, "internal-record", compress)
                self.assertEqual(sorted(p.name for p in (d / "out").iterdir()),
                                 [f"demo.internal-record.part001{ext}", "demo.internal-record.parts.sha256"])
                part = d / "out" / f"demo.internal-record.part001{ext}"
                self.assertEqual(sums.read_text(encoding="utf-8"), f"{kit.sha(part.read_bytes())}  {part.name}\n")
                inner = INNER + "".join(f"{p}\t{n}\t{kit.sha(data)}\t{cls}\t001\n"
                                        for p, n, data, cls in kit.ROWS if cls == "internal-record")
                want = {"demo.internal-record/LIST.tsv": inner.encode()}
                want.update({f"demo.internal-record/files/{p}": data
                             for p, n, data, cls in kit.ROWS if cls == "internal-record"})
                self.assertEqual(kit.members(part), want)

    def test_the_public_archive_names_no_internal_file_and_the_internal_one_no_public_file(self):
        d, proj, lst = kit.listed(self)
        for cls, key in (("public-result", ("--published", kit.published(d, proj))), ("internal-record", ()),
                         ("public-input", ())):
            rc, _, err = kit.run("pack", "--list", lst, "--class", cls, "--name", "demo", "--out", d / "out",
                                 "--root", f"proj={proj}", *key)
            self.assertEqual((rc, err), (0, ""))
        held = {cls: kit.members(d / "out" / f"demo.{cls}.part001.tar.gz")
                for cls in ("public-result", "internal-record", "public-input")}
        self.assertEqual(sorted(held["public-result"]), ["demo.public-result/LIST.tsv"] +
                         [f"demo.public-result/files/{p}" for p in kit.PUBLIC])
        self.assertEqual(sorted(held["public-input"]),
                         ["demo.public-input/LIST.tsv", "demo.public-input/files/proj/inputs/genome1.fna"])
        for cls, other in (("public-result", kit.INTERNAL), ("internal-record", kit.PUBLIC),
                           ("public-input", kit.INTERNAL + kit.PUBLIC)):
            everything = b"\n".join(name.encode() + b"\n" + data for name, data in held[cls].items())
            for path in other:
                self.assertNotIn(path.encode(), everything, f"{cls} names {path}")
        self.assertNotIn(b"a message", b"".join(held["public-result"].values()))

    def test_a_public_archive_needs_every_file_in_a_published_checksum_file(self):
        d, proj, lst = kit.listed(self)
        base = ("pack", "--list", lst, "--class", "public-result", "--name", "demo", "--out", d / "out",
                "--root", f"proj={proj}")
        rc, out, err = kit.run(*base)
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("public-result needs --published", err)
        (d / "partial.sha256").write_text(f"{kit.sha(kit.FILES['reports/alpha.json'])}  alpha.json\n", encoding="utf-8")
        rc, out, err = kit.run(*base, "--published", d / "partial.sha256")
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("1 files listed as public-result stand in no --published file with their sha256 and their "
                      "name", err)
        self.assertIn("  proj/reports/beta.json\n", err)
        self.assertFalse((d / "out").exists())
        (d / "rest.sha256").write_text(f"{kit.sha(kit.FILES['reports/beta.json'])} *beta.json\n", encoding="utf-8")
        rc, _, err = kit.run(*base, "--published", d / "partial.sha256", "--published", d / "rest.sha256")
        self.assertEqual((rc, err), (0, ""))

    def test_a_rule_that_calls_a_note_public_does_not_get_it_into_the_public_archive(self):
        d, proj, rules = kit.project(self)
        rules.write_text("public-result\tproj/notes/ledger.md\n" + kit.RULES, encoding="utf-8")
        rc, _, _ = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 0)
        rc, _, err = kit.run("pack", "--list", d / "list.tsv", "--class", "public-result", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={proj}", "--published", kit.published(d, proj))
        self.assertEqual(rc, 2)
        self.assertIn("  proj/notes/ledger.md\n", err)
        self.assertFalse((d / "out").exists())

    def test_regenerable_and_refused_are_never_packed(self):
        d, proj, lst = kit.listed(self)
        for cls, said in (("regenerable", "regenerable is never packed"), ("refused", "refused is never packed"),
                          ("notes", "unknown class 'notes'")):
            rc, out, err = kit.run("pack", "--list", lst, "--class", cls, "--name", "demo", "--out", d / "out",
                                   "--root", f"proj={proj}")
            self.assertEqual((rc, out), (2, ""))
            self.assertIn(said, err)
        self.assertFalse((d / "out").exists())

    def test_a_list_that_gives_a_barred_path_a_packed_class_stops_the_pack(self):
        d, proj, lst = kit.listed(self)
        kit.write(proj, {"notes/ksd.md": b"x"})
        lst.write_text(lst.read_text(encoding="utf-8") + f"proj/notes/ksd.md\t1\t{kit.sha(b'x')}\tinternal-record\n",
                       encoding="utf-8")
        rc, out, err = kit.run("pack", "--list", lst, "--class", "internal-record", "--name", "demo",
                               "--out", d / "out", "--root", f"proj={proj}")
        self.assertEqual((rc, out), (3, ""))
        self.assertIn("refused: proj/notes/ksd.md holds the barred word 'ksd'", err)
        self.assertFalse((d / "out").exists())

    def test_a_root_at_a_barred_place_stops_the_pack(self):
        d, proj, lst = kit.listed(self)
        os.rename(proj, d / "vumc-copy")
        rc, _, err = kit.run("pack", "--list", lst, "--class", "internal-record", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={d / 'vumc-copy'}")
        self.assertEqual(rc, 3)
        self.assertIn(f"pack: refused: the root proj={d / 'vumc-copy'} holds the barred word 'vumc'", err)
        self.assertFalse((d / "out").exists())

    def test_a_folder_of_the_root_that_leads_to_a_barred_place_stops_the_pack(self):
        d, proj, lst = kit.listed(self)
        os.rename(proj / "notes", d / "sigdrug-notes")
        os.symlink(d / "sigdrug-notes", proj / "notes")
        rc, _, err = kit.run("pack", "--list", lst, "--class", "internal-record", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={proj}")
        self.assertEqual(rc, 3)
        self.assertIn(f"pack: refused: proj/notes/alpha-copy.json is at {proj / 'notes' / 'alpha-copy.json'}, which "
                      f"holds the barred word 'sigdrug'", err)
        self.assertFalse((d / "out").exists())


class PublishedByNameAndContent(unittest.TestCase):
    """A public file needs its sha256 and its base name together in a --published file. notes/alpha-copy.json holds
    the bytes of reports/alpha.json under another name, and notes/empty.txt is empty."""

    def pack_public(self, rules_first, published_text):
        d, proj, rules = kit.project(self)
        rules.write_text(rules_first + kit.RULES, encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
        self.assertEqual(rc, 0, err)
        (d / "release.sha256").write_text(published_text, encoding="utf-8")
        rc, out, err = kit.run("pack", "--list", d / "list.tsv", "--class", "public-result", "--name", "demo",
                               "--out", d / "out", "--root", f"proj={proj}", "--published", d / "release.sha256")
        return d, rc, out, err

    RELEASE = "".join(f"{kit.sha(kit.FILES[f'reports/{n}'])}  {n}\n" for n in ("alpha.json", "beta.json"))

    def test_an_internal_copy_of_a_published_file_is_not_packed_as_public(self):
        for rule in ("public-result\tproj/notes/alpha-copy.json\n", "public-result\tproj/*.json\n"):
            with self.subTest(rule=rule):
                d, rc, out, err = self.pack_public(rule, self.RELEASE)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn("1 files listed as public-result stand in no --published file with their sha256 and "
                              "their name, so they are not packed as public:\n  proj/notes/alpha-copy.json\n", err)
                self.assertFalse((d / "out").exists())

    def test_an_empty_file_needs_its_name_too(self):
        d, rc, out, err = self.pack_public("public-result\tproj/notes/empty.txt\n",
                                           self.RELEASE + f"{kit.EMPTY_SHA256}  .keep\n")
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("  proj/notes/empty.txt\n", err)
        self.assertFalse((d / "out").exists())
        d, rc, out, err = self.pack_public("public-result\tproj/notes/empty.txt\n",
                                           self.RELEASE + f"{kit.EMPTY_SHA256}  data/empty.txt\n")
        self.assertEqual((rc, err), (0, ""))
        self.assertIn("demo.public-result/files/proj/notes/empty.txt",
                      kit.members(d / "out" / "demo.public-result.part001.tar.gz"))

    def test_every_name_a_published_file_gives_for_one_sha256_counts(self):
        alpha = kit.sha(kit.FILES["reports/alpha.json"])
        for first, second in (("alpha.json", "alpha-copy.json"), ("alpha-copy.json", "alpha.json")):
            with self.subTest(first=first):
                d, rc, out, err = self.pack_public(
                    "public-result\tproj/notes/alpha-copy.json\n",
                    f"{alpha}  {first}\n{alpha}  ./sub/{second}\n{kit.sha(kit.FILES['reports/beta.json'])}  beta.json\n")
                self.assertEqual((rc, err), (0, ""))
                self.assertEqual(sorted(kit.members(d / "out" / "demo.public-result.part001.tar.gz")),
                                 ["demo.public-result/LIST.tsv", "demo.public-result/files/proj/notes/alpha-copy.json",
                                  "demo.public-result/files/proj/reports/alpha.json",
                                  "demo.public-result/files/proj/reports/beta.json"])


class Parts(unittest.TestCase):
    """12 synthetic files of 40,000 bytes that do not compress. In a tar each takes 40,960 bytes: a 512-byte header
    and 79 blocks. With --max-bytes 150000 a part has room for three (122,880) and not for four (163,840), whatever
    the list and the tar end take of the rest, so the 12 files make 4 parts."""

    def big(self):
        d = kit.tmp(self)
        files = {f"data/f{i:02d}.bin": kit.noise(i, 40000) for i in range(12)}
        kit.write(d / "proj", files)
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, ""))
        return d, files

    def test_no_part_is_over_the_limit_and_every_file_is_in_the_part_the_list_names(self):
        for compress in ("gz", "xz", "none"):
            with self.subTest(compress=compress):
                d, files = self.big()
                rc, out, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                                       "--out", d / "out", "--root", f"proj={d / 'proj'}", "--max-bytes", 150000,
                                       "--compress", compress)
                self.assertEqual((rc, err), (0, ""))
                sums = d / "out" / "demo.internal-record.parts.sha256"
                parts = kit.parts_of(sums)
                self.assertEqual([p.name[:25] for p in parts], [f"demo.internal-record.part{i:03d}"[:25]
                                                                for i in (1, 2, 3, 4)])
                self.assertIn("12 files of class internal-record, 480000 bytes, in 4 parts", out)
                inner = INNER + "".join(f"proj/data/f{i:02d}.bin\t40000\t{kit.sha(files[f'data/f{i:02d}.bin'])}\t"
                                        f"internal-record\t{i // 3 + 1:03d}\n" for i in range(12))
                for n, part in enumerate(parts):
                    self.assertLessEqual(part.stat().st_size, 150000)
                    self.assertGreater(part.stat().st_size, 120000)
                    want = {"demo.internal-record/LIST.tsv": inner.encode()}
                    want.update({f"demo.internal-record/files/proj/data/f{i:02d}.bin": files[f"data/f{i:02d}.bin"]
                                 for i in range(3 * n, 3 * n + 3)})
                    self.assertEqual(kit.members(part), want)
                rc, out, _ = kit.run("check", "--parts", sums)
                self.assertEqual(rc, 0, out)
                self.assertIn("4 parts with the listed sha256, 12 files, 480000 bytes", out)
                rc, out, _ = kit.run("check", kit.unpack(sums, d / "un"))
                self.assertEqual(rc, 0, out)

    def test_a_part_leaves_room_for_data_that_grows_when_compressed(self):
        """Data that does not compress comes out of gzip a little larger: about 5 bytes more for every 16,383.

        The limit here is 64 MiB, 67,108,864 bytes. Two synthetic files of noise: big.bin with 67,075,072 bytes
        (67,075,584 in the tar with its header) and small.bin with 20,480 (20,992). The list takes 1,024 and the tar
        end is counted at 11,264. As plain tar bytes all of it comes to 67,108,864, exactly the limit, and gzip would
        add about 20,000 bytes to it. So the room of a part is first cut by what compression can add: 20,546 bytes
        at this limit. big.bin then still fits (446 bytes to spare) and small.bin goes to a second part."""
        d = kit.tmp(self)
        kit.write(d / "proj", {"big.bin": kit.noise("big", 67_075_072), "small.bin": kit.noise("small", 20_480)})
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, ""))
        rc, out, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                               "--out", d / "out", "--root", f"proj={d / 'proj'}", "--max-bytes", 67_108_864)
        self.assertEqual((rc, err), (0, ""), out)
        sums = d / "out" / "demo.internal-record.parts.sha256"
        one, two = kit.parts_of(sums)
        self.assertLessEqual(one.stat().st_size, 67_108_864)
        self.assertGreater(one.stat().st_size, 67_075_072)  # the noise did not shrink
        self.assertEqual(sorted(kit.members(two)), ["demo.internal-record/LIST.tsv",
                                                    "demo.internal-record/files/proj/small.bin"])
        rc, out, _ = kit.run("check", "--parts", sums)
        self.assertEqual(rc, 0, out)
        self.assertIn("2 parts with the listed sha256, 2 files, 67095552 bytes", out)  # 67,075,072 + 20,480

    def test_a_file_too_large_for_a_part_stops_the_pack_and_nothing_is_written(self):
        d, files = self.big()
        rc, out, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                               "--out", d / "out", "--root", f"proj={d / 'proj'}", "--max-bytes", 40000)
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("proj/data/f00.bin needs 40960 bytes with its header", err)
        self.assertIn("nothing was written", err)
        self.assertFalse((d / "out").exists())

    def test_the_same_files_pack_to_the_same_bytes(self):
        d, files = self.big()
        for out in ("one", "two"):
            rc, _, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                                 "--out", d / out, "--root", f"proj={d / 'proj'}", "--max-bytes", 150000)
            self.assertEqual((rc, err), (0, ""))
        self.assertEqual((d / "one" / "demo.internal-record.parts.sha256").read_bytes(),
                         (d / "two" / "demo.internal-record.parts.sha256").read_bytes())

    def test_existing_parts_are_not_written_over(self):
        d, proj, sums = kit.packed(self)
        before = sums.read_bytes()
        rc, out, err = kit.run("pack", "--list", d / "list.tsv", "--class", "internal-record", "--name", "demo",
                               "--out", d / "out", "--root", f"proj={proj}")
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("exists; nothing was written", err)
        self.assertEqual(sums.read_bytes(), before)


class TheSourceMustBeAsListed(unittest.TestCase):
    def pack(self, d, proj):
        return kit.run("pack", "--list", d / "list.tsv", "--class", "internal-record", "--name", "demo",
                       "--out", d / "out", "--root", f"proj={proj}")

    def test_a_file_changed_since_the_list_stops_the_pack_and_leaves_no_part(self):
        for name, change in (("same size", b"# Ledger\n\nrow 2\n"), ("longer", b"# Ledger\n\nrow 1\nrow 2\n"),
                             ("shorter", b"# Ledger\n")):
            with self.subTest(name):
                d, proj, lst = kit.listed(self)
                (proj / "notes" / "ledger.md").write_bytes(change)
                rc, out, err = self.pack(d, proj)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn("proj/notes/ledger.md changed since it was listed", err)
                self.assertEqual(list((d / "out").iterdir()), [])

    def test_a_file_that_is_gone_or_became_a_link_stops_the_pack(self):
        d, proj, lst = kit.listed(self)
        (proj / "notes" / "ledger.md").unlink()
        rc, _, err = self.pack(d, proj)
        self.assertEqual(rc, 2)
        self.assertIn("cannot open proj/notes/ledger.md", err)
        (d / "elsewhere.md").write_bytes(kit.FILES["notes/ledger.md"])  # the listed content, somewhere else
        os.symlink(d / "elsewhere.md", proj / "notes" / "ledger.md")
        rc, _, err = self.pack(d, proj)
        self.assertEqual(rc, 2)
        self.assertIn("cannot open proj/notes/ledger.md", err)
        self.assertEqual(list((d / "out").iterdir()), [])

    def test_a_path_with_no_root_stops_the_pack(self):
        d, proj, lst = kit.listed(self)
        rc, _, err = kit.run("pack", "--list", lst, "--class", "internal-record", "--name", "demo",
                             "--out", d / "out", "--root", f"other={proj}")
        self.assertEqual(rc, 2)
        self.assertIn("no --root LABEL=PATH holds proj/notes/alpha-copy.json", err)

    def test_long_names_and_names_outside_ascii_come_back_as_they_were(self):
        d = kit.tmp(self)
        long = "deep/" + "/".join(["a-folder-with-a-long-name"] * 6) + "/" + "n" * 120 + ".md"  # 281 characters
        files = {long: b"long", "café/résumé.md": b"accents", "with space/a b.md": b"space"}
        kit.write(d / "proj", files)
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        rc, _, err = kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        self.assertEqual((rc, err), (0, ""))
        rc, _, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={d / 'proj'}", "--max-bytes", 30000)
        self.assertEqual((rc, err), (0, ""))
        sums = d / "out" / "demo.internal-record.parts.sha256"
        got = kit.members(kit.parts_of(sums)[0])
        self.assertEqual({k: v for k, v in got.items() if "/files/" in k},
                         {f"demo.internal-record/files/proj/{p}": data for p, data in files.items()})
        self.assertLessEqual(kit.parts_of(sums)[0].stat().st_size, 30000)
        for args in (("--parts", sums), (kit.unpack(sums, d / "un"),)):
            rc, out, _ = kit.run("check", *args)
            self.assertEqual(rc, 0, out)

    def test_a_list_that_is_not_one_stops_the_pack(self):
        d, proj, lst = kit.listed(self)
        good = lst.read_text(encoding="utf-8")
        row = f"proj/notes/ledger.md\t16\t{kit.sha(kit.FILES['notes/ledger.md'])}\tinternal-record\n"
        self.assertIn(row, good)
        for text, said in ((good.replace("path\tbytes", "bytes\tpath"), "not a list"),
                           (good + row, "is listed twice"),
                           (good.replace(row, row.replace("\t16\t", "\t-\t")), "bytes is '-'"),
                           (good.replace(row, row.replace("internal-record", "notes")), "unknown class 'notes'"),
                           (good.replace(row, row.replace("proj/notes/ledger.md", "/etc/passwd")),
                            "is not a plain relative path"),
                           (good.replace(row, row.replace("proj/notes", "proj/../..")), "is not a plain relative path"),
                           (good.replace(row, row.replace("ledger.md", "ledger.md/")), "is a directory row"),
                           (good.replace(row, row.replace("internal-record", "refused")), "holds no barred word")):
            with self.subTest(said=said):
                lst.write_text(text, encoding="utf-8")
                rc, out, err = self.pack(d, proj)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn(said, err)
                self.assertFalse((d / "out").exists())


if __name__ == "__main__":
    unittest.main()
