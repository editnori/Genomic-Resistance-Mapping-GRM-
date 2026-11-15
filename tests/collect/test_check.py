"""collect.py check: an archive that is as its list says passes, and each planted difference is named.

Every archive here is packed from the synthetic project of kit.py. The internal-record archive holds 4 files and
47 bytes (21 + 0 + 16 + 10)."""

import gzip
import io
import lzma
import os
import shutil
import tarfile
import unittest

from tests.collect import kit

OK_FOLDER = ("4 files, 47 bytes, class internal-record, every size and sha256 as listed, "
             "nothing else under files/\n")


class AnUnpackedArchive(unittest.TestCase):
    def unpacked(self, cls="internal-record"):
        d, proj, sums = kit.packed(self, cls)
        return d, kit.unpack(sums, d / "un")

    def fails(self, top, *said, good):
        rc, out, err = kit.run("check", top)
        self.assertEqual((rc, err), (1, ""), out)
        for s in said:
            self.assertIn(s, out)
        self.assertTrue(out.endswith(f"check: FAILED {top}: {good} of 4 listed files are as listed\n"), out)

    def test_an_archive_as_listed_passes(self):
        d, top = self.unpacked()
        rc, out, err = kit.run("check", top)
        self.assertEqual((rc, err, out), (0, "", f"check: OK {top}: " + OK_FOLDER))
        d, top = self.unpacked("public-result")
        rc, out, err = kit.run("check", top, "--threads", 1)
        self.assertEqual((rc, out), (0, f"check: OK {top}: 2 files, 41 bytes, class public-result, every size and "
                                        "sha256 as listed, nothing else under files/\n"))

    def test_a_folder_named_for_another_class_is_found(self):
        d, top = self.unpacked()
        rc, out, err = kit.run("check", f"{top}/")
        self.assertEqual((rc, err, out), (0, "", f"check: OK {top}/: " + OK_FOLDER))
        renamed = top.parent / "demo.public-result"
        os.rename(top, renamed)
        self.fails(renamed, "the list is of class internal-record and the folder is named demo.public-result",
                   good=4)

    def test_one_changed_byte_of_the_same_size_is_found(self):
        d, top = self.unpacked()
        (top / "files/proj/notes/ledger.md").write_bytes(b"# Ledger\n\nrow 2\n")
        self.fails(top, "1 files have the listed size and another sha256:\n  proj/notes/ledger.md\n", good=3)

    def test_a_change_past_the_first_block_of_a_large_file_is_found(self):
        d = kit.tmp(self)
        data = kit.noise(2, (1 << 20) + 4097)
        kit.write(d / "proj", {"big.bin": data})
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        rc, _, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={d / 'proj'}")
        self.assertEqual((rc, err), (0, ""))
        top = kit.unpack(d / "out" / "demo.internal-record.parts.sha256", d / "un")
        rc, out, _ = kit.run("check", top)
        self.assertEqual(rc, 0, out)
        (top / "files/proj/big.bin").write_bytes(data[:-1] + bytes([data[-1] ^ 1]))  # the last byte, one bit
        rc, out, _ = kit.run("check", top)
        self.assertEqual(rc, 1)
        self.assertIn("1 files have the listed size and another sha256:\n  proj/big.bin\n", out)

    def test_a_shorter_or_longer_file_is_found(self):
        d, top = self.unpacked()
        (top / "files/proj/notes/ledger.md").write_bytes(b"# Ledger\n")
        (top / "files/proj/notes/empty.txt").write_bytes(b"x")
        self.fails(top, "2 files have another size:\n"
                        "  proj/notes/empty.txt (the list has 0 bytes, the file has 1)\n"
                        "  proj/notes/ledger.md (the list has 16 bytes, the file has 9)\n", good=2)

    def test_a_missing_file_is_found_with_its_part(self):
        d, top = self.unpacked()
        (top / "files/proj/notes/messages/m1.txt").unlink()
        self.fails(top, "1 listed files are missing (1 of part 001):\n  proj/notes/messages/m1.txt\n", good=3)

    def test_a_file_the_list_does_not_name_is_found(self):
        d, top = self.unpacked()
        (top / "files/proj/notes/added.md").write_bytes(b"new")
        (top / "files/proj/reports").mkdir()
        (top / "files/proj/reports/alpha.json").write_bytes(kit.FILES["reports/alpha.json"])
        self.fails(top, "2 files are not in the list:\n  proj/notes/added.md\n  proj/reports/alpha.json\n", good=4)

    def test_a_moved_file_is_missing_and_extra(self):
        d, top = self.unpacked()
        os.rename(top / "files/proj/notes/ledger.md", top / "files/proj/notes/Ledger.md")
        self.fails(top, "1 listed files are missing (1 of part 001):\n  proj/notes/ledger.md\n",
                   "1 files are not in the list:\n  proj/notes/Ledger.md\n", good=3)

    def test_two_files_with_their_contents_exchanged_are_found(self):
        d, top = self.unpacked()
        (top / "files/proj/notes/ledger.md").write_bytes(b"a message\n" + b"123456")  # 16 bytes
        (top / "files/proj/notes/messages/m1.txt").write_bytes(b"# Ledger\n\n")  # 10 bytes
        self.fails(top, "2 files have the listed size and another sha256:\n  proj/notes/ledger.md\n"
                        "  proj/notes/messages/m1.txt\n", good=2)

    def test_a_link_in_place_of_a_file_is_not_followed(self):
        d, top = self.unpacked()
        (d / "same.md").write_bytes(kit.FILES["notes/ledger.md"])
        (top / "files/proj/notes/ledger.md").unlink()
        os.symlink(d / "same.md", top / "files/proj/notes/ledger.md")
        self.fails(top, "1 listed files are missing", "1 links or special files, which an archive does not hold:\n"
                                                       "  proj/notes/ledger.md\n", good=3)

    def test_anything_beside_the_list_and_the_files_is_found(self):
        d, top = self.unpacked()
        (top / "NOTES.md").write_bytes(b"x")
        self.fails(top, "1 entries beside LIST.tsv and files/:\n  NOTES.md\n", good=4)

    def test_a_barred_name_is_found_and_not_opened(self):
        d, top = self.unpacked()
        kit.write(top / "files/proj/notes", {"Clinical/t.tsv": b"x", "ks50.txt": b"y"})
        os.chmod(top / "files/proj/notes/Clinical", 0)
        os.chmod(top / "files/proj/notes/ks50.txt", 0)
        self.fails(top, "2 names hold a barred word; none was opened:\n  proj/notes/Clinical\n"
                        "  proj/notes/ks50.txt\n", good=4)

    def test_a_list_edited_after_packing_is_found(self):
        d, top = self.unpacked()
        text = (top / "LIST.tsv").read_text(encoding="utf-8")
        row = f"proj/notes/ledger.md\t16\t{kit.sha(kit.FILES['notes/ledger.md'])}\tinternal-record\t001\n"
        self.assertIn(row, text)
        cases = (
            (text.replace(row, ""), 1, "1 files are not in the list:\n  proj/notes/ledger.md\n"),
            (text.replace(row, row.replace(kit.sha(kit.FILES["notes/ledger.md"]), "0" * 64)), 1,
             "1 files have the listed size and another sha256"),
            (text.replace(row, row.replace("\t16\t", "\t17\t")), 1, "the list has 17 bytes, the file has 16"),
            (text.replace(row, row.replace("internal-record", "public-result")), 1,
             "the list names the classes internal-record, public-result; an archive holds one packed class"),
            (text.replace(row, row.replace("\t001\n", "\t1\n")), 2, "part is '1'"),
            (text.replace(row, row.replace("\t001\n", "\n")), 2, "4 cells, want 5"),
            (text.replace(row, row.replace("ledger.md", "sealed.md")), 3, "holds the barred word 'sealed'"),
        )
        for new, status, said in cases:
            with self.subTest(said=said):
                (top / "LIST.tsv").write_text(new, encoding="utf-8")
                rc, out, err = kit.run("check", top)
                self.assertEqual(rc, status, out + err)
                self.assertIn(said, out + err)

    def test_a_list_of_a_class_that_is_never_packed_fails(self):
        d, top = self.unpacked()
        (top / "LIST.tsv").write_text("path\tbytes\tsha256\tclass\tpart\n"
                                      "proj/build/target/\t194\t-\tregenerable\t001\n", encoding="utf-8")
        rc, out, _ = kit.run("check", top)
        self.assertEqual(rc, 1)
        self.assertIn("the list names the classes regenerable; an archive holds one packed class", out)

    def test_a_folder_with_no_list_is_not_checked(self):
        d, top = self.unpacked()
        (top / "LIST.tsv").unlink()
        rc, out, err = kit.run("check", top)
        self.assertEqual((rc, out), (2, ""))
        self.assertIn("LIST.tsv is not there", err)
        rc, _, err = kit.run("check")
        self.assertEqual(rc, 2)
        rc, _, err = kit.run("check", top, "--parts", d / "out" / "demo.internal-record.parts.sha256")
        self.assertEqual(rc, 2)
        self.assertIn("give one unpacked folder, or --parts", err)


class ThePartsWhereTheyLie(unittest.TestCase):
    """check --parts reads the packed parts as streams. Nothing is unpacked."""

    OK = "check: OK demo.internal-record: 1 parts with the listed sha256, " + OK_FOLDER

    def fails(self, sums, *said):
        rc, out, err = kit.run("check", "--parts", sums)
        self.assertEqual((rc, err), (1, ""), out)
        for s in said:
            self.assertIn(s, out)
        self.assertIn("check: FAILED demo.internal-record", out)
        return out

    def test_parts_as_packed_pass_and_nothing_is_unpacked(self):
        for compress in ("gz", "xz", "none"):
            with self.subTest(compress=compress):
                d, proj, sums = kit.packed(self, "internal-record", compress)
                before = sorted(p.name for p in sums.parent.iterdir())
                rc, out, err = kit.run("check", "--parts", sums, cwd=d)
                self.assertEqual((rc, err, out), (0, "", self.OK))
                self.assertEqual(sorted(p.name for p in sums.parent.iterdir()), before)
                self.assertEqual(sorted(p.name for p in d.iterdir()), ["list.tsv", "out", "proj", "rules.tsv"])

    def test_a_part_with_one_changed_byte_is_found(self):
        d, proj, sums = kit.packed(self, "internal-record", "none")
        part = kit.parts_of(sums)[0]
        data = bytearray(part.read_bytes())
        at = data.index(b"row 1")  # inside notes/ledger.md
        data[at + 4] = ord("2")
        part.write_bytes(bytes(data))
        self.fails(sums, f"the part {part.name} has the sha256 {kit.sha(bytes(data))}, and {sums.name} has ",
                   "1 files have the listed size and another sha256:\n  proj/notes/ledger.md\n")

    def test_a_cut_part_and_a_missing_part_are_found(self):
        d, proj, sums = kit.packed(self)
        part = kit.parts_of(sums)[0]
        whole = part.read_bytes()
        part.write_bytes(whole[:len(whole) // 2])
        out = self.fails(sums, f"the part {part.name} has the sha256 ")
        self.assertIn(f"the part {part.name} fails when read to the end of its gzip stream: the stream stops before "
                      f"its end", out)
        part.unlink()
        self.fails(sums, f"the part {part.name} is not beside {sums.name}", "no list could be read")

    def test_a_change_hidden_by_a_new_parts_file_is_found_by_the_list_inside(self):
        d, proj, sums = kit.packed(self)
        part = kit.parts_of(sums)[0]
        kit.rewrite(part, lambda name, data: [(name, b"# Ledger\n\nrow 9\n" if name.endswith("ledger.md") else data)])
        kit.resum(sums)
        out = self.fails(sums, "1 files have the listed size and another sha256:\n  proj/notes/ledger.md\n")
        self.assertNotIn("has the sha256", out)

    def test_members_the_list_does_not_name_are_found(self):
        cases = (
            ("an added file", lambda n, data: [(n, data)] + ([("demo.internal-record/files/proj/notes/new.md", b"n")]
                                                             if n.endswith("LIST.tsv") else []),
             "1 files are not in the list:\n  proj/notes/new.md\n"),
            ("a file twice", lambda n, data: [(n, data)] * (2 if n.endswith("m1.txt") else 1),
             "proj/notes/messages/m1.txt is in the archive twice"),
            ("a file left out", lambda n, data: [] if n.endswith("m1.txt") else [(n, data)],
             "1 listed files are missing (1 of part 001):\n  proj/notes/messages/m1.txt\n"),
            ("a member outside the folder", lambda n, data: [(n, data)] + ([("elsewhere/x", b"x")]
                                                                           if n.endswith("LIST.tsv") else []),
             "1 entries beside LIST.tsv and files/:\n  elsewhere/x\n"),
            ("a barred name", lambda n, data: [(n, data)] + ([("demo.internal-record/files/proj/SigDrug/v.tsv", b"x")]
                                                             if n.endswith("LIST.tsv") else []),
             "1 names hold a barred word; none was opened:\n  proj/SigDrug/v.tsv\n"),
            ("no list", lambda n, data: [] if n.endswith("LIST.tsv") else [(n, data)],
             "holds no demo.internal-record/LIST.tsv"),
        )
        for name, change, said in cases:
            with self.subTest(name):
                d, proj, sums = kit.packed(self)
                kit.rewrite(kit.parts_of(sums)[0], change)
                kit.resum(sums)
                self.fails(sums, said)

    def test_a_link_inside_a_part_is_found(self):
        d, proj, sums = kit.packed(self, "internal-record", "none")
        part = kit.parts_of(sums)[0]
        old = kit.members(part)
        with tarfile.open(part, "w", format=tarfile.PAX_FORMAT) as tar:
            for name, data in old.items():
                info = tarfile.TarInfo(name)
                if name.endswith("ledger.md"):
                    info.type, info.linkname = tarfile.SYMTYPE, "/etc/hostname"
                    tar.addfile(info)
                else:
                    info.size = len(data)
                    tar.addfile(info, io.BytesIO(data))
        kit.resum(sums)
        self.fails(sums, "1 listed files are missing", "1 links or special files, which an archive does not hold:\n"
                                                       "  proj/notes/ledger.md\n")

    def test_the_part_of_another_class_under_this_name_is_found(self):
        d, proj, sums = kit.packed(self)
        rc, _, err = kit.run("pack", "--list", d / "list.tsv", "--class", "public-input", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={proj}")
        self.assertEqual((rc, err), (0, ""))
        part = kit.parts_of(sums)[0]
        shutil.copyfile(d / "out" / "demo.public-input.part001.tar.gz", part)
        kit.resum(sums)
        self.fails(sums, "holds no demo.internal-record/LIST.tsv", "no list could be read")


    def test_a_list_of_another_class_inside_this_archive_is_found(self):
        d, proj, sums = kit.packed(self)
        kit.rewrite(kit.parts_of(sums)[0], lambda n, data: [(n, data.replace(b"internal-record", b"public-input")
                                                             if n.endswith("LIST.tsv") else data)])
        kit.resum(sums)
        self.fails(sums, "the list is of class public-input and the archive is named demo.internal-record")

    def test_files_in_another_part_than_the_list_says_and_parts_that_disagree(self):
        d = kit.tmp(self)
        kit.write(d / "proj", {f"f{i}.bin": kit.noise(i, 40000) for i in range(6)})
        (d / "rules.tsv").write_text("internal-record\t*\n", encoding="utf-8")
        kit.run("list", "--rules", d / "rules.tsv", "--root", f"proj={d / 'proj'}", "--out", d / "l.tsv")
        rc, _, err = kit.run("pack", "--list", d / "l.tsv", "--class", "internal-record", "--name", "demo",
                             "--out", d / "out", "--root", f"proj={d / 'proj'}", "--max-bytes", 150000)
        self.assertEqual((rc, err), (0, ""))
        sums = d / "out" / "demo.internal-record.parts.sha256"
        one, two = kit.parts_of(sums)
        with self.subTest("the two parts under each other's name"):
            one.rename(d / "held")
            two.rename(one)
            (d / "held").rename(two)
            out = self.fails(sums, "6 files are in another part than listed:\n"
                                   "  proj/f0.bin (the list has part 001, the file is in part 002)\n")
            self.assertIn("has the sha256", out)
            kit.resum(sums)
            self.fails(sums, "6 files are in another part than listed")
            two.rename(d / "held")
            one.rename(two)
            (d / "held").rename(one)
            kit.resum(sums)
            self.assertEqual(kit.run("check", "--parts", sums)[0], 0)
        with self.subTest("a part whose list differs"):
            kit.rewrite(two, lambda n, data: [(n, data.replace(b"\t002\n", b"\t001\n") if n.endswith("LIST.tsv")
                                               else data)])
            kit.resum(sums)
            self.fails(sums, "the parts do not hold the same LIST.tsv")
        with self.subTest("a part the parts file no longer names"):
            sums.write_text(sums.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")
            self.fails(sums, "the list has files in the parts 002, which demo.internal-record.parts.sha256 does not "
                             "name", "3 listed files are missing (3 of part 002)")

    def test_a_parts_file_that_is_not_one_is_not_checked(self):
        d, proj, sums = kit.packed(self)
        line = sums.read_text(encoding="utf-8")
        for name, text, said in (("demo.internal-record.parts.sha256", line.replace("part001", "part1"),
                                  "is not one part of demo.internal-record"),
                                 ("demo.internal-record.parts.sha256", line.replace("demo.", "../demo."),
                                  "is not one part of demo.internal-record"),
                                 ("demo.internal-record.parts.sha256", "", "names no part"),
                                 ("demo.internal-record.parts.sha256", "not a line\n", "not a sha256sum line"),
                                 ("parts.txt", line, "--parts takes a NAME.CLASS.parts.sha256 file")):
            with self.subTest(said=said):
                (sums.parent / name).write_text(text, encoding="utf-8")
                rc, out, err = kit.run("check", "--parts", sums.parent / name)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn(said, err)


if __name__ == "__main__":
    unittest.main()


class EachPartIsReadToItsEnd(unittest.TestCase):
    """check --parts reads each part to the end of its gzip or xz stream, so the stream's own check is read, and
    after the tar end mark nothing but zero bytes may follow. The parts file is rewritten to match each change, as
    someone hiding it would, so only the reading of the stream can find it."""

    STREAM = {"gz": "gzip", "xz": "xz", "none": "tar"}

    def changed(self, compress, change):
        """Packs the project, changes its one part's bytes through change(bytes) -> bytes, and checks it."""
        d, proj, sums = kit.packed(self, "internal-record", compress)
        part = kit.parts_of(sums)[0]
        part.write_bytes(change(part.read_bytes()))
        kit.resum(sums)
        rc, out, err = kit.run("check", "--parts", sums)
        self.assertEqual((rc, err), (1, ""), out)
        self.assertIn(f"check: the part {part.name} fails when read to the end of its {self.STREAM[compress]} "
                      f"stream: ", out)
        self.assertIn("check: FAILED demo.internal-record", out)
        self.assertNotIn("has the sha256", out)
        return out

    @staticmethod
    def tar_inside(compress, data):
        return gzip.decompress(data) if compress == "gz" else lzma.decompress(data) if compress == "xz" else data

    @staticmethod
    def packed_again(compress, data):
        return (gzip.compress(data, mtime=0) if compress == "gz" else
                lzma.compress(data, format=lzma.FORMAT_XZ) if compress == "xz" else data)

    @staticmethod
    def end_mark_at(tar_bytes):
        with tarfile.open(fileobj=io.BytesIO(tar_bytes), mode="r|") as tar:
            for _ in tar:
                pass
            return tar.offset

    def test_a_damaged_check_at_the_end_of_the_stream_is_found(self):
        def flip(at):
            return lambda data: data[:at] + bytes([data[at] ^ 1]) + data[at + 1:]
        cases = (("gz", -8, "the gzip CRC"), ("gz", -4, "the gzip length"), ("xz", -12, "the xz stream footer"))
        for compress, at, name in cases:
            with self.subTest(name):
                self.changed(compress, flip(at))

    def test_bytes_after_the_compressed_stream_are_found(self):
        cases = (("gz", b"x", "a byte"), ("gz", gzip.compress(b"more", mtime=0), "a second gzip member"),
                 ("xz", b"x", "a byte"), ("xz", lzma.compress(b"more"), "a second xz stream"),
                 ("xz", b"\0" * 4, "xz stream padding"))
        for compress, more, name in cases:
            with self.subTest(compress=compress, case=name):
                out = self.changed(compress, lambda data: data + more)
                self.assertIn("bytes follow the end of the stream", out)

    def test_bytes_after_the_tar_end_mark_are_found(self):
        second = io.BytesIO()
        with tarfile.open(fileobj=second, mode="w", format=tarfile.PAX_FORMAT) as tar:
            info = tarfile.TarInfo("demo.internal-record/files/proj/notes/hidden.md")
            info.size = 6
            tar.addfile(info, io.BytesIO(b"hidden"))
        for compress in ("gz", "xz", "none"):
            for more, name in ((b"x", "a byte"), (second.getvalue(), "a second tar")):
                with self.subTest(compress=compress, case=name):
                    out = self.changed(compress, lambda data: self.packed_again(
                        compress, self.tar_inside(compress, data) + more))
                    self.assertIn("bytes that are not zero follow its last member", out)
                    self.assertNotIn("hidden.md", out)

    def test_a_missing_or_cut_tar_end_mark_is_found(self):
        for compress in ("gz", "xz", "none"):
            for keep, name in ((0, "no end mark"), (tarfile.BLOCKSIZE, "one zero block")):
                with self.subTest(compress=compress, case=name):
                    def cut(data):
                        tar = self.tar_inside(compress, data)
                        return self.packed_again(compress, tar[:self.end_mark_at(tar) + keep])
                    out = self.changed(compress, cut)
                    self.assertIn("the tar end mark is missing or cut short", out)
