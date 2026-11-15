"""docs/data-collection.md: its length, its classes and barred words, and its commands, run as written.

The commands run in a temporary folder on a small synthetic record tree. `/path/to/records` is that tree, and
`/path/to/release` a synthetic release folder beside it."""

import hashlib
import importlib.util
import os
import re
import shlex
import subprocess
import sys
import unittest

from tests.collect import kit

DOCUMENT = kit.ROOT / "docs" / "data-collection.md"


def tool():
    spec = importlib.util.spec_from_file_location("collect_tool", kit.TOOL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def blocks(text, kind):
    return re.findall(r"^```" + kind + r"\n(.*?)^```$", text, flags=re.S | re.M)


class TheDocument(unittest.TestCase):
    def setUp(self):
        self.text = DOCUMENT.read_text(encoding="utf-8")

    def test_it_is_under_400_words(self):
        self.assertLess(len(self.text.split()), 400)

    def test_it_names_the_tool_s_classes_barred_words_and_part_size(self):
        t = tool()
        rows = re.findall(r"^\| `([a-z-]+)` \|.*$", self.text, flags=re.M)
        self.assertEqual(tuple(rows), t.CLASSES)
        refused = next(line for line in self.text.split("\n") if line.startswith("| `refused` |"))
        self.assertEqual(tuple(re.findall(r"`([^`]+)`", refused)[1:]), t.BARRED)
        self.assertEqual(t.BARRED, kit.BARRED)
        self.assertEqual(t.PART_BYTES, 1_900_000_000)
        self.assertIn("each at most 1.9 GB (`--max-bytes`)", self.text)
        for cls in t.PACKABLE:
            self.assertRegex(self.text, rf"\| `{cls}` \| yes \|")
        named = self.text.split("names a `regenerable` directory. ")[1].split(" need no rule.")[0]
        self.assertEqual(tuple(re.findall(r"`([^`]+)`", named)), t.REGENERABLE_DIRS + t.REGENERABLE_MARKERS)

    def test_it_says_what_each_check_proves_and_what_the_guard_does_not_see(self):
        text = " ".join(self.text.split())
        for said in ("becomes a row in no archive: look at it by hand",
                     "The guard reads names only: a hard link with a clean name is read.",
                     "`--published`, the checksum file of a release already published, lists it by name and content.",
                     "Each class has its own folder", "Send `archive/public` alone.",
                     "Parts hold file times, so a second pack gives other bytes; the list repeats.",
                     "Check the parts first:", "It trusts the `parts.sha256` beside them, so keep a copy apart.",
                     "Then unpack every part:", "`check DIR` proves only that the folder matches the list inside it."):
            self.assertIn(said, text)
        self.assertLess(text.index("$C check --parts"), text.index("tar -xzf"))

    def test_its_commands_run_as_written(self):
        """The page's commands on a small synthetic tree, with one change: the internal archive is cut into parts of
        at most 17,000 bytes, so that its restore must unpack each of three parts."""
        d = kit.tmp(self)
        alpha, beta = b'{"dataset": "alpha"}\n', b'{"dataset": "beta"}\n'
        records = {"release/alpha.json": alpha, "release/beta.json": beta, "notes/ledger.md": b"# Ledger\n",
                   "notes/a.bin": kit.noise(1, 3000), "notes/b.bin": kit.noise(2, 3000),
                   "notes/c.bin": kit.noise(3, 3000), "web/node_modules/pkg/index.js": b"js"}
        kit.write(d / "records", records)
        sums = f"{kit.sha(alpha)}  alpha.json\n{kit.sha(beta)}  beta.json\n"
        kit.write(d / "release", {"SHA256SUMS": sums.encode()})
        internal = [v for k, v in records.items() if k.startswith("notes/")]
        public = [alpha, beta]
        rules, = blocks(self.text, "tsv")
        (d / "rules.tsv").write_text(rules, encoding="utf-8")
        os.symlink(kit.ROOT / "scripts", d / "scripts")
        script = "".join(blocks(self.text, "sh"))
        self.assertEqual(script.count("/path/to/records"), 3)
        self.assertEqual(script.count("/path/to/release"), 1)
        self.assertEqual(script.count("python3 -I scripts/collect/collect.py"), 1)
        self.assertEqual(script.count("--out archive/internal\n"), 1)
        script = script.replace("--out archive/internal\n", "--out archive/internal --max-bytes 17000\n")
        script = script.replace("/path/to/records", shlex.quote(str(d / "records")))
        script = script.replace("/path/to/release", shlex.quote(str(d / "release")))
        script = script.replace("python3 -I", f"{shlex.quote(sys.executable)} -I -B -W error")
        # No standard input: with a socket there and an empty environment, bash reads its startup files.
        p = subprocess.run(["bash", "-euc", script], cwd=d, capture_output=True, text=True,
                           stdin=subprocess.DEVNULL)
        self.assertEqual((p.returncode, p.stderr), (0, ""), p.stdout)
        self.assertEqual(sorted(f.name for f in (d / "archive").iterdir()), ["internal", "public"])
        self.assertEqual(sorted(f.name for f in (d / "archive" / "internal").iterdir()),
                         [f"grm-20261114.internal-record.part00{n}.tar.gz" for n in (1, 2, 3)] +
                         ["grm-20261114.internal-record.parts.sha256"])
        self.assertEqual(sorted(f.name for f in (d / "archive" / "public").iterdir()), [
            "grm-20261114.public-result.part001.tar.gz", "grm-20261114.public-result.parts.sha256"])
        files = f"{len(internal)} files, {sum(map(len, internal))} bytes, class internal-record"
        self.assertIn(f"check: OK grm-20261114.internal-record: 3 parts with the listed sha256, {files}", p.stdout)
        self.assertIn(f"check: OK restore/grm-20261114.internal-record: {files}", p.stdout)
        n = sum(map(len, public))
        self.assertIn(f"\npublic-result\t2\t{n}\t2\t{n}\t2\t{n}\t0\t0\n", p.stdout)
        n = len(records["web/node_modules/pkg/index.js"])
        self.assertIn(f"\nregenerable\t1\t{n}\t0\t0\t0\t0\t1\t{n}\n", p.stdout)
        members = kit.members(d / "archive" / "public" / "grm-20261114.public-result.part001.tar.gz")
        self.assertEqual(sorted(members), ["grm-20261114.public-result/LIST.tsv",
                                           "grm-20261114.public-result/files/records/release/alpha.json",
                                           "grm-20261114.public-result/files/records/release/beta.json"])

if __name__ == "__main__":
    unittest.main()
