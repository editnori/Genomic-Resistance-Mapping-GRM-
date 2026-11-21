"""The tree itself: no tracked file holds a private detail, and the page that held most of them stays out."""

import unittest

from tests.public_words import rules


class Tree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.names = rules.tracked()
        cls.hits, cls.unread = rules.findings()

    def test_no_tracked_file_holds_a_private_detail(self):
        shown = "\n".join(f"{kind} {path}:{number}" for kind, path, number in self.hits[:60])
        self.assertEqual(len(self.hits), 0, f"{rules.counts(self.hits)}\n{shown}")

    def test_the_scan_reads_the_tree(self):
        for name in ("README.md", "Cargo.toml", "src/main.rs", "docs/ci.md", ".github/workflows/ci.yml",
                     rules.SUITE + "/rules.py"):
            self.assertIn(name, self.names)
        self.assertGreater(len(self.names), 500)

    def test_only_compressed_files_are_left_unread(self):
        self.assertEqual([name for name in self.unread if not name.endswith(".gz")], [])

    def test_the_removed_page_stays_out(self):
        self.assertFalse(rules.ROOT.joinpath(*rules.REMOVED_PAGE).exists())
        self.assertNotIn("/".join(rules.REMOVED_PAGE), self.names)


if __name__ == "__main__":
    unittest.main()
