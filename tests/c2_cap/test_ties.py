"""SYNTHETIC tests of the count of genomes tied at the 4,096th place of the selection rule (scripts/c2-cap/ties.py).

Every genome, call and id here is made up in the test.
"""

import json
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

ties = common.load("ties")
SPECIES = [line.split("\t")[:2] for line in (common.KIT / "species.tsv").read_text(encoding="utf-8").splitlines()[1:]
           if line.strip()]


def calls(species, resistant, susceptible, first=1):
    """Made-up phenotypes: `resistant` genomes resistant and `susceptible` susceptible to drug a."""
    out = []
    for n in range(resistant + susceptible):
        out.append({"genome_id": f"9000.{first + n}", "species": species, "antibiotic": "drug a",
                    "reported": "R" if n < resistant else "S", "conflict": False, "label_from_mic": False})
    return out


class Ties(unittest.TestCase):
    def test_the_genomes_sharing_the_score_at_the_cut_are_counted_kept_and_left(self):
        # 2,000 susceptible genomes score 1/2,000 and rank first; the 3,000 resistant ones tie at 1/3,000.
        phen = calls("Synthetica borealis", 3000, 2000)
        row = ties.tie_row(phen, [{"species": "Synthetica borealis", "antibiotic": "drug a"}], "Synthetica borealis",
                           "synthetica_borealis")
        self.assertEqual((row["ranked"], row["cut"], row["tied"], row["kept"], row["left"]), (5000, "4096", 3000, 2096, 904))

    def test_a_cut_on_a_score_held_by_one_genome_counts_one(self):
        phen = calls("Synthetica borealis", 1, 4095) + calls("Synthetica borealis", 0, 1, first=9001)
        phen += [dict(phen[-1], antibiotic="drug b", reported="R")]
        data = [{"species": "Synthetica borealis", "antibiotic": a} for a in ("drug a", "drug b")]
        row = ties.tie_row(phen, data, "Synthetica borealis", "synthetica_borealis")
        self.assertEqual(row["ranked"], 4097)
        self.assertEqual((row["tied"], row["kept"], row["left"]), (4095, 4094, 1))

    def test_a_species_ranking_4096_or_fewer_has_no_cut(self):
        phen = calls("Synthetica borealis", 10, 20)
        row = ties.tie_row(phen, [{"species": "Synthetica borealis", "antibiotic": "drug a"}], "Synthetica borealis",
                           "synthetica_borealis")
        self.assertEqual((row["ranked"], row["cut"], row["tied"]), (30, "no cut", ""))

    def test_the_table_holds_one_row_per_species_and_no_genome_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            cur = Path(tmp) / "curated"
            cur.mkdir()
            (cur / "species.json").write_text(json.dumps([{"slug": s, "species": n} for s, n in SPECIES]))
            (cur / "datasets.json").write_text(json.dumps([{"species": n, "antibiotic": "drug a"} for _, n in SPECIES]))
            phen = calls(SPECIES[0][1], 3000, 2000) + calls(SPECIES[1][1], 5, 5, first=7001)
            (cur / "phenotypes.jsonl").write_text("".join(json.dumps(p) + "\n" for p in phen))
            out = Path(tmp) / "ties.tsv"
            import contextlib
            import io
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(ties.main([str(cur), str(out)]), 0)
            text = out.read_text()
            rows = [dict(zip(ties.COLUMNS, line.split("\t"))) for line in text.splitlines()[1:]]
            self.assertEqual([r["species_slug"] for r in rows], [s for s, _ in SPECIES])
            self.assertEqual((rows[0]["tied"], rows[0]["kept"], rows[0]["left"]), ("3000", "2096", "904"))
            self.assertEqual((rows[1]["ranked"], rows[1]["cut"]), ("10", "no cut"))
            self.assertNotIn("9000.", text)
            self.assertNotIn("7001.", text)


if __name__ == "__main__":
    unittest.main()
