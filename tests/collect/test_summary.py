"""collect.py summary: how many places hold each content, counted by sha256.

Three synthetic places. The contents are single letters repeated, so every size below can be counted by hand:

  content   bytes   machine-a                      machine-b            release.sha256
  A         100     notes/a.md, notes/a-copy.md    notes/a.md           -
  B         200     notes/b.md                     -                    -
  C         300     reports/c.json                 -                    named
  D         400     -                              reports/d.json       named
  E         500     notes/e.md (internal-record)   pub/e.json (public)  -
  F         600     -                              -                    named
  -         7000    web/node_modules/ (one row)    -                    -
  -         -       -                              notes/ks30 (refused) -
"""

import unittest

from tests.collect import kit

A, B, C, D, E, F = (letter.encode() * n for letter, n in zip("ABCDEF", (100, 200, 300, 400, 500, 600)))
HEAD = "path\tbytes\tsha256\tclass\n"


def row(path, data, cls):
    return f"{path}\t{len(data)}\t{kit.sha(data)}\t{cls}\n"


class Copies(unittest.TestCase):
    def places(self):
        d = kit.tmp(self)
        (d / "a.tsv").write_text(HEAD + row("notes/a-copy.md", A, "internal-record")
                                 + row("notes/a.md", A, "internal-record")
                                 + row("notes/b.md", B, "internal-record") + row("notes/e.md", E, "internal-record")
                                 + row("reports/c.json", C, "public-result")
                                 + "web/node_modules/\t7000\t-\tregenerable\n", encoding="utf-8")
        (d / "b.tsv").write_text(HEAD + row("notes/a.md", A, "internal-record") + "notes/ks30\t-\t-\trefused\n"
                                 + row("pub/e.json", E, "public-result") + row("reports/d.json", D, "public-result"),
                                 encoding="utf-8")
        (d / "release.sha256").write_text("".join(f"{kit.sha(x)}  {n}\n" for x, n in ((C, "c.json"), (D, "d.json"),
                                                                                     (F, "f.json"))), encoding="utf-8")
        return d

    def test_the_tables(self):
        d = self.places()
        rc, out, err = kit.run("summary", "--place", f"machine-a={d / 'a.tsv'}", "--place", f"machine-b={d / 'b.tsv'}",
                               "--hashes", f"release={d / 'release.sha256'}", "--out", d / "all.tsv")
        self.assertEqual((rc, err), (0, ""))
        tables = [t.split("\n") for t in out.split("\n\n")]
        self.assertEqual(tables[0][1:], [
            "class\tentries\tbytes\tcontents\tbytes_without_duplicates\tcontents_in_one_place\tbytes_in_one_place"
            "\tentries_not_hashed\tbytes_not_hashed",
            # public-result: C 300, E 500, D 400. All three are in two places (C and D in the release, E in a and b).
            "public-result\t3\t1200\t3\t1200\t0\t0\t0\t0",
            # internal-record: A three times (300), B 200, E 500 = 1000 in 5 entries; without duplicates A B E = 800.
            # Only B is in one place.
            "internal-record\t5\t1000\t3\t800\t1\t200\t0\t0",
            "regenerable\t1\t7000\t0\t0\t0\t0\t1\t7000",
            "refused\t1\t0\t0\t0\t0\t0\t1\t0",
            # all: 10 entries, 1200 + 1000 + 7000; contents A B C D E = 1500; in one place: B.
            "all\t10\t9200\t5\t1500\t1\t200\t2\t7000",
        ])
        self.assertEqual(tables[1][1:], [
            "class\tplaces\theld_by\tcontents\tbytes",
            "public-result\t2\tmachine-a+machine-b\t1\t500",
            "public-result\t2\tmachine-a+release\t1\t300",
            "public-result\t2\tmachine-b+release\t1\t400",
            "internal-record\t1\tmachine-a\t1\t200",
            "internal-record\t2\tmachine-a+machine-b\t2\t600",
            "all\t1\tmachine-a\t1\t200",
            "all\t2\tmachine-a+machine-b\t2\t600",
            "all\t2\tmachine-a+release\t1\t300",
            "all\t2\tmachine-b+release\t1\t400",
        ])
        self.assertEqual(tables[2][1:], [
            "place\tclass\tentries\tbytes\tbytes_without_duplicates",
            "machine-a\tpublic-result\t1\t300\t300",
            "machine-a\tinternal-record\t4\t900\t800",  # A twice: 100 + 100 + 200 + 500, and 100 less
            "machine-a\tregenerable\t1\t7000\t7000",
            "machine-b\tpublic-result\t2\t900\t900",
            "machine-b\tinternal-record\t1\t100\t100",
            "machine-b\trefused\t1\t0\t0",
        ])
        # The release names C, D and F; F is in no list.
        self.assertEqual(tables[3][1:], ["place\tsha256_named\tin_no_list", "release\t3\t1",
                                         f"summary: wrote 10 rows to {d / 'all.tsv'}", ""])
        self.assertEqual((d / "all.tsv").read_text(encoding="utf-8"), "".join(line + "\n" for line in (
            "place\tpath\tbytes\tsha256\tclass\tplaces\theld_by",
            f"machine-a\tnotes/a-copy.md\t100\t{kit.sha(A)}\tinternal-record\t2\tmachine-a+machine-b",
            f"machine-a\tnotes/a.md\t100\t{kit.sha(A)}\tinternal-record\t2\tmachine-a+machine-b",
            f"machine-a\tnotes/b.md\t200\t{kit.sha(B)}\tinternal-record\t1\tmachine-a",
            f"machine-a\tnotes/e.md\t500\t{kit.sha(E)}\tinternal-record\t2\tmachine-a+machine-b",
            f"machine-a\treports/c.json\t300\t{kit.sha(C)}\tpublic-result\t2\tmachine-a+release",
            "machine-a\tweb/node_modules/\t7000\t-\tregenerable\t-\t-",
            f"machine-b\tnotes/a.md\t100\t{kit.sha(A)}\tinternal-record\t2\tmachine-a+machine-b",
            "machine-b\tnotes/ks30\t-\t-\trefused\t-\t-",
            f"machine-b\tpub/e.json\t500\t{kit.sha(E)}\tpublic-result\t2\tmachine-a+machine-b",
            f"machine-b\treports/d.json\t400\t{kit.sha(D)}\tpublic-result\t2\tmachine-b+release",
        )))

    def test_one_place_alone_has_every_content_once(self):
        d = self.places()
        rc, out, _ = kit.run("summary", "--place", f"machine-a={d / 'a.tsv'}")
        self.assertEqual(rc, 0)
        self.assertIn("\nall\t6\t8200\t4\t1100\t4\t1100\t1\t7000\n", out)  # A B C E: 100 + 200 + 300 + 500
        self.assertNotIn("summary: wrote", out)

    def test_inputs_a_summary_cannot_use(self):
        d = self.places()
        for args, said in ((("--place", f"a={d / 'a.tsv'}", "--place", f"a={d / 'b.tsv'}"), "a is given twice"),
                           (("--place", f"a={d / 'a.tsv'}", "--hashes", f"a={d / 'release.sha256'}"),
                            "is both a --place and a --hashes name"),
                           (("--place", f"a/b={d / 'a.tsv'}",), "want NAME=VALUE with a plain name"),
                           (("--place", f"a={d / 'release.sha256'}",), "not a list"),
                           (("--place", f"a={d / 'a.tsv'}", "--hashes", f"r={d / 'b.tsv'}"), "not a sha256sum line"),
                           (("--place", f"a={d / 'absent.tsv'}",), "cannot read the list")):
            with self.subTest(said=said):
                rc, out, err = kit.run("summary", *args)
                self.assertEqual((rc, out), (2, ""))
                self.assertIn(said, err)

    def test_a_place_name_with_a_barred_word_is_refused(self):
        d = self.places()
        for args, said in ((("--place", f"KS50={d / 'a.tsv'}"), "the place name KS50 holds the barred word 'ks50'"),
                           (("--place", f"a={d / 'a.tsv'}", "--hashes", f"sealed-release={d / 'release.sha256'}"),
                            "the place name sealed-release holds the barred word 'sealed'")):
            with self.subTest(said=said):
                rc, out, err = kit.run("summary", *args)
                self.assertEqual((rc, out), (3, ""))
                self.assertEqual(err, f"summary: refused: {said}\n")


if __name__ == "__main__":
    unittest.main()
