"""docs/data-collection.md: its length, its classes and barred words, and its commands, run as written.

The commands run in a temporary folder on a small synthetic record tree. `/path/to/records` is that tree."""

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

    def test_its_commands_run_as_written(self):
        d = kit.tmp(self)
        alpha, beta = b'{"dataset": "alpha"}\n', b'{"dataset": "beta"}\n'
        kit.write(d / "records", {
            "release/alpha.json": alpha, "release/beta.json": beta,
            "release/SHA256SUMS": (f"{hashlib.sha256(alpha).hexdigest()}  alpha.json\n"
                                   f"{hashlib.sha256(beta).hexdigest()}  beta.json\n").encode(),
            "notes/ledger.md": b"# Ledger\n", "web/node_modules/pkg/index.js": b"js"})
        rules, = blocks(self.text, "tsv")
        (d / "rules.tsv").write_text(rules, encoding="utf-8")
        os.symlink(kit.ROOT / "scripts", d / "scripts")
        script = "".join(blocks(self.text, "sh"))
        self.assertEqual(script.count("/path/to/records"), 4)
        self.assertEqual(script.count("python3 -I scripts/collect/collect.py"), 1)
        script = script.replace("/path/to/records", shlex.quote(str(d / "records")))
        script = script.replace("python3 -I", f"{shlex.quote(sys.executable)} -I -B -W error")
        # No standard input: with a socket there and an empty environment, bash reads its startup files.
        p = subprocess.run(["bash", "-euc", script], cwd=d, capture_output=True, text=True,
                           stdin=subprocess.DEVNULL)
        self.assertEqual((p.returncode, p.stderr), (0, ""), p.stdout)
        self.assertEqual(sorted(f.name for f in (d / "archive").iterdir()), [
            "grm-20261114.internal-record.part001.tar.gz", "grm-20261114.internal-record.parts.sha256",
            "grm-20261114.public-result.part001.tar.gz", "grm-20261114.public-result.parts.sha256"])
        # internal-record: SHA256SUMS (153 bytes: lines of 64 + 2 + 10 + 1 and 64 + 2 + 9 + 1) and ledger.md (9)
        self.assertIn("check: OK grm-20261114.internal-record: 1 parts with the listed sha256, 2 files, 162 bytes, "
                      "class internal-record", p.stdout)
        self.assertIn("check: OK restore/grm-20261114.internal-record: 2 files, 162 bytes, class internal-record",
                      p.stdout)
        self.assertIn("\npublic-result\t2\t41\t2\t41\t2\t41\t0\t0\n", p.stdout)
        self.assertIn("\nregenerable\t1\t2\t0\t0\t0\t0\t1\t2\n", p.stdout)
        public = kit.members(d / "archive" / "grm-20261114.public-result.part001.tar.gz")
        self.assertEqual(sorted(public), ["grm-20261114.public-result/LIST.tsv",
                                          "grm-20261114.public-result/files/records/release/alpha.json",
                                          "grm-20261114.public-result/files/records/release/beta.json"])


if __name__ == "__main__":
    unittest.main()
