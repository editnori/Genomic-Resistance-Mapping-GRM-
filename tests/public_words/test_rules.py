"""Each rule on made-up lines: what it must find and what it must leave alone. Every example is put together
from pieces here, so that this file holds no hit itself."""

import unittest

from tests.public_words import rules

S = "/"
B = chr(92)
MADE_UP = frozenset({rules.sha("xavier")})
PUBLIC = ".".join(["8"] * 4)   # an address outside the private ranges, put together here


def kinds(text, path="docs/page.md", names=MADE_UP):
    return [kind for kind, _ in rules.hits_in(path, text, names)]


class Paths(unittest.TestCase):
    def test_a_home_path_under_a_name(self):
        self.assertEqual(kinds("ran in " + S + "home" + S + "zed" + S + "review"), ["K1"])
        self.assertEqual(kinds("ran in " + S + "Users" + S + "zed" + S + "review"), ["K1"])

    def test_neutral_home_paths_pass(self):
        for name in sorted(rules.NEUTRAL_HOMES) + ["...", "$USER", "<name>"]:
            self.assertEqual(kinds("ran in " + S + "home" + S + name + S + "review"), [], name)

    def test_a_windows_user_path(self):
        self.assertEqual(kinds("C" + ":" + B + "Users" + B + "zed" + B + "x"), ["K2"])
        self.assertEqual(kinds("C" + ":" + B + B + "Users" + B + B + "zed"), ["K2"])
        self.assertEqual(kinds(S + "mnt" + S + "c" + S + "Users" + S + "zed" + S + "x"), ["K2"])

    def test_a_windows_path_outside_users_passes(self):
        self.assertEqual(kinds("C" + ":" + B + "grm-stage" + B + "x"), [])

    def test_the_line_number(self):
        text = "one\ntwo " + S + "home" + S + "zed\nthree"
        self.assertEqual(rules.hits_in("docs/page.md", text, MADE_UP), [("K1", 2)])


class LinksAndMachines(unittest.TestCase):
    def test_a_thread_link(self):
        self.assertEqual(kinds("see t3" + "-thread" + ":" + S + S + "v2/abc"), ["K3"])
        self.assertEqual(kinds("see t3" + "-context" + ":" + S + S + "v1/abc"), ["K3"])

    def test_a_box_id(self):
        self.assertEqual(kinds("box bx" + "_" + "a1b2c3d4 ran it"), ["K4"])

    def test_a_stub_id_passes(self):
        self.assertEqual(kinds("ID=bx" + "_" + "stub1"), [])
        self.assertEqual(kinds("ID=bx" + "_" + "a1b2c3d4e5"), [])

    def test_a_host_name(self):
        self.assertEqual(kinds("ssh3." + "vast" + ".ai port 2"), ["K5"])
        self.assertEqual(kinds("on box-a." + "ts" + ".net"), ["K5"])
        self.assertEqual(kinds("on DESKTOP" + "-AB12CDE"), ["K5"])

    def test_a_public_address(self):
        self.assertEqual(kinds("at " + PUBLIC + " now"), ["K5"])

    def test_other_addresses_pass(self):
        for text in ("127.0.0.1", "0.0.0.0", "10.1.2.3", "192.168.1.9", "172.20.0.1", "169.254.1.1", "192.0.2.7",
                     "198.51.100.7", "203.0.113.7", "999.1.1.1", "(EC 2.3.1.82)", "1.2.3", "1.2.3.4.5"):
            self.assertEqual(kinds("x " + text + " y"), [], text)

    def test_a_lock_file_may_hold_a_version_of_four_parts(self):
        self.assertEqual(kinds("version " + PUBLIC, path="Cargo.lock"), [])
        self.assertEqual(kinds("version " + PUBLIC, path="web/package-lock.json"), [])


class Names(unittest.TestCase):
    def test_a_barred_name_in_any_case(self):
        for text in ("Xavier ruled", "XAVIER ruled", "xavier ruled", "Xavier's ruling", "xavier_review", "xavier888",
                     "XavierSmith", "by-xavier"):
            self.assertEqual(kinds(text), ["K6"], text)

    def test_each_time_it_stands(self):
        self.assertEqual(kinds("Xavier and xavier\nand Xavier"), ["K6", "K6", "K6"])

    def test_other_words_pass(self):
        self.assertEqual(kinds("xaviera, saxavier and xavi"), [])

    def test_a_licence_and_an_author_line_may_hold_it(self):
        self.assertEqual(kinds("Copyright Xavier", path="LICENSE"), [])
        self.assertEqual(kinds("Copyright Xavier", path="web/LICENSE.md"), [])
        self.assertEqual(kinds('authors = ["Xavier"]', path="Cargo.toml"), [])
        self.assertEqual(kinds('  "author": "Xavier",', path="web/package.json"), [])
        self.assertEqual(kinds("written by Xavier", path="Cargo.toml"), ["K6"])

    def test_a_path_can_be_a_hit(self):
        self.assertEqual(kinds("", path="docs/xavier-notes.md"), [])
        self.assertEqual([k for k, _ in rules.hits_in("", "docs/xavier-notes.md", MADE_UP)], ["K6"])

    def test_the_names_are_held_as_sha256(self):
        self.assertEqual(len(rules.NAMES), 2)
        for value in rules.NAMES:
            self.assertRegex(value, r"^[0-9a-f]{64}$")
        self.assertEqual(rules.sha("XAVIER"), rules.sha("xavier"))


class Counts(unittest.TestCase):
    def test_hits_and_files_by_kind(self):
        got = rules.counts([("K1", "a", 1), ("K1", "a", 2), ("K1", "b", 1), ("K6", "a", 0)])
        self.assertEqual(got, {"K1": (3, 2), "K2": (0, 0), "K3": (0, 0), "K4": (0, 0), "K5": (0, 0), "K6": (1, 1)})


if __name__ == "__main__":
    unittest.main()
