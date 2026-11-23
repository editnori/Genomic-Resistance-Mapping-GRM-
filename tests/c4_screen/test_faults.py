"""Each planted fault of the runners whose files C4 changes still applies to this tree: C3's
(scripts/c3-inputs/planted_faults.py), C9's (scripts/c9-pick/planted_faults.py) and C4's own. A runner stops at
the first fault whose text to change does not occur exactly once in its file, so a line C4 adds that repeats
another block's planted text stops that block's runner. Reads the committed runners and the files they name; reads
no data and plants nothing.
"""

import importlib.util
import unittest

from tests.c4_screen import common

RUNNERS = ("c3-inputs", "c9-pick", "c4-screen")


def runner(kit):
    spec = importlib.util.spec_from_file_location(f"c4_screen_faults_{kit.replace('-', '_')}",
                                                  common.ROOT / "scripts" / kit / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Faults(unittest.TestCase):
    def test_each_runners_fault_text_occurs_once_in_its_file(self):
        for kit in RUNNERS:
            faults = runner(kit).FAULTS
            self.assertTrue(faults, kit)
            for fault in faults:
                n, path, old = fault[0], fault[2], fault[3]
                with self.subTest(runner=kit, fault=n):
                    self.assertEqual((common.ROOT / path).read_text(encoding="utf-8").count(old), 1, path)


if __name__ == "__main__":
    unittest.main()
