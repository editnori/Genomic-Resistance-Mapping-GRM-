"""SYNTHETIC tests of the step that takes score.py's three tables out of a species' score.tar (scripts/c2-cap/take_tables.py).

Every archive here is made in the test, with made-up tables and a made-up score.log. Each read the step makes is
recorded, and no read may touch score.log's bytes; no file named score.log may be written.
"""

import contextlib
import gzip
import io
import tarfile
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

take = common.load("take_tables")
LOG = b"score.py's output, made up: 0 of 0\n" * 32  # 1,152 bytes
TABLE = {"checks.tsv": b"row\tcheck\tpassed\tdetail\n", "draws.tsv": b"dataset\tarm\n" * 50,
         "datasets.tsv": b"dataset\tsplit\n"}


class Audit(io.RawIOBase):
    """A file that records the byte ranges read from it."""

    def __init__(self, data):
        self.data, self.pos, self.reads = data, 0, []

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, pos, whence=0):
        self.pos = pos if whence == 0 else self.pos + pos if whence == 1 else len(self.data) + pos
        return self.pos

    def tell(self):
        return self.pos

    def read(self, n=-1):
        n = len(self.data) - self.pos if n is None or n < 0 else n
        out = self.data[self.pos:self.pos + n]
        self.reads.append((self.pos, self.pos + len(out)))
        self.pos += len(out)
        return out


def archive(members, fmt=tarfile.GNU_FORMAT):
    """Bytes of a tar holding members [(name, bytes or None for a folder)], and score.log's data range."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=fmt) as t:
        for name, data in members:
            info = tarfile.TarInfo(name)
            if data is None:
                info.type = tarfile.DIRTYPE
                t.addfile(info)
            else:
                info.size = len(data)
                t.addfile(info, io.BytesIO(data))
    raw = buf.getvalue()
    with tarfile.open(fileobj=io.BytesIO(raw)) as t:
        logs = [(m.offset_data, m.offset_data + m.size) for m in t.getmembers() if m.name.endswith("score.log")]
    return raw, logs


def driver(prefix="score/", log=LOG):
    """As `tar -cf score.tar score` makes it: the folder, then its files, score.log among them."""
    return [(prefix.rstrip("/"), None), (prefix + "checks.tsv", TABLE["checks.tsv"]), (prefix + "score.log", log),
            (prefix + "draws.tsv", TABLE["draws.tsv"]), (prefix + "datasets.tsv", TABLE["datasets.tsv"])]


class Take(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def taken(self, members):
        raw, logs = archive(members)
        f = Audit(raw)
        found = take.take(f)
        for lo, hi in logs:
            for a, b in f.reads:
                self.assertFalse(a < hi and lo < b, f"read {a} to {b} touches score.log's bytes {lo} to {hi}")
        return found

    def test_the_three_tables_are_taken_and_score_log_is_never_read_or_written(self):
        for log in (LOG, LOG[:1024], b""):
            with self.subTest(log_bytes=len(log)):
                found = self.taken(driver(log=log))
                self.assertEqual(found, TABLE)
        lines = take.write(found, self.tmp / "species")
        self.assertEqual(sorted(p.name for p in (self.tmp / "species" / "score").iterdir()), sorted(TABLE))
        self.assertEqual(len(lines), 3)
        self.assertEqual([p for p in self.tmp.rglob("*") if p.name == "score.log"], [])

    def test_members_named_with_a_dot_folder_are_taken(self):
        self.assertEqual(self.taken(driver(prefix="./score/")), TABLE)

    def test_no_other_member_is_taken(self):
        found = self.taken(driver() + [("score/other.tsv", b"x\n"), ("results/checks.tsv", b"y\n")])
        self.assertEqual(found, TABLE)

    def test_a_compressed_archive_is_refused(self):
        raw, _ = archive(driver())
        with self.assertRaises(take.Refuse):
            take.take(io.BytesIO(gzip.compress(raw)))

    def test_a_missing_table_is_refused(self):
        with self.assertRaises(take.Refuse):
            self.taken([m for m in driver() if not m[0].endswith("draws.tsv")])

    def test_a_table_given_twice_is_refused(self):
        with self.assertRaises(take.Refuse):
            self.taken(driver() + [("./score/checks.tsv", b"other\n")])

    def test_a_link_member_is_refused(self):
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w", format=tarfile.GNU_FORMAT) as t:
            for name, data in driver():
                if name.endswith("datasets.tsv"):
                    link = tarfile.TarInfo(name)
                    link.type, link.linkname = tarfile.SYMTYPE, "score.log"
                    t.addfile(link)
                elif data is not None:
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    t.addfile(info, io.BytesIO(data))
        with self.assertRaises(take.Refuse):
            take.take(io.BytesIO(buf.getvalue()))

    def test_an_archive_cut_short_is_refused(self):
        raw, _ = archive(driver())
        with self.assertRaises(take.Refuse):
            take.take(io.BytesIO(raw[:1536]))

    def test_a_score_folder_that_already_holds_a_table_is_refused(self):
        found = self.taken(driver())
        (self.tmp / "species" / "score").mkdir(parents=True)
        (self.tmp / "species" / "score" / "draws.tsv").write_bytes(b"old\n")
        with self.assertRaises(take.Refuse):
            take.write(found, self.tmp / "species")
        self.assertEqual((self.tmp / "species" / "score" / "draws.tsv").read_bytes(), b"old\n")

    def test_the_command_writes_the_three_tables_from_a_file(self):
        raw, _ = archive(driver())
        (self.tmp / "score.tar").write_bytes(raw)
        with contextlib.redirect_stdout(io.StringIO()) as printed:
            self.assertEqual(take.main([str(self.tmp / "score.tar"), str(self.tmp / "species")]), 0)
        self.assertEqual(len(printed.getvalue().splitlines()), 3)
        for name, data in TABLE.items():
            self.assertEqual((self.tmp / "species" / "score" / name).read_bytes(), data)
        self.assertFalse((self.tmp / "species" / "score" / "score.log").exists())


if __name__ == "__main__":
    unittest.main()
