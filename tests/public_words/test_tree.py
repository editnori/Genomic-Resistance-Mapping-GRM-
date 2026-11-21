"""The tree itself: no tracked file holds a private detail, and the page that held most of them stays out."""

import sys
import unittest

from tests.public_words import rules


class Tree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.names = rules.tracked()
        cls.checked = rules.checked_names()
        cls.hits, cls.unread = rules.findings(names=cls.checked)

    def test_no_tracked_file_holds_a_private_detail(self):
        shown = "\n".join(rules.listing(self.hits[:60], self.checked))
        self.assertEqual(len(self.hits), 0, f"{rules.counts(self.hits)}\n{shown}")

    def test_the_extra_words(self):
        """Says whether the extra words were checked. Without them the suite still passes."""
        print(rules.extra_state(), file=sys.stderr, flush=True)
        unusable = rules.extra_names()[1]
        self.assertEqual(unusable, 0, f"lines of {rules.EXTRA} that are not one run of letters and digits: {unusable}")

    def test_the_scan_reads_the_tree(self):
        for name in ("README.md", "Cargo.toml", "src/main.rs", "docs/ci.md", ".github/workflows/ci.yml",
                     rules.SUITE + "/rules.py"):
            self.assertIn(name, self.names)
        self.assertGreater(len(self.names), 500)

    def test_only_compressed_files_are_left_unread(self):
        left = [rules.shown(name, self.checked) for name in self.unread if not name.endswith(".gz")]
        self.assertEqual(left, [])

    def test_the_removed_page_stays_out(self):
        self.assertFalse(rules.ROOT.joinpath(*rules.REMOVED_PAGE).exists())
        self.assertNotIn("/".join(rules.REMOVED_PAGE), self.names)


if __name__ == "__main__":
    unittest.main()
