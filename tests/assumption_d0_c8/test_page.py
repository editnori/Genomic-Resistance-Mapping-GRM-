"""docs/refit-result.md against the committed tables of reference/refit-result/.

This is the one test file of the suite that reads committed result tables:
the aggregates the reading program wrote on its one run, by species,
antibiotic, learner and split. It fits nothing and reads no genome.
"""

import collections
import csv
import hashlib
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts", "assumption-d0-c8"))

import page_tables  # noqa: E402

RECORD = os.path.join(ROOT, "reference", "refit-result")
PAGE = os.path.join(ROOT, "docs", "refit-result.md")
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
FOUR = ["absent", "removed", "lineage", "unused"]


def rows(name):
    with open(os.path.join(RECORD, name), newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def text(name):
    with open(os.path.join(RECORD, name)) as handle:
        return handle.read()


def causes_of(row):
    return [] if row["cause"] == "no missed fold" else row["cause"].split("+")


class Page(unittest.TestCase):
    def setUp(self):
        with open(PAGE) as handle:
            self.text = handle.read()
        self.lines = self.text.split("\n")

    def test_every_printed_line_stands_in_the_page(self):
        tables = page_tables.tables(RECORD)
        self.assertEqual(
            sorted(tables),
            ["by_species", "causes", "checks", "family", "names", "order", "overlap", "shares", "tied"],
        )
        for name, printed in tables.items():
            if name == "names":
                self.assertEqual(len(printed), 23)
                for line in printed:
                    self.assertEqual(self.lines.count(line), 1, line)
                continue
            self.assertGreater(len(printed), 2, name)
            self.assertIn(printed[1], self.lines, name)
            for line in printed[:1] + printed[2:]:
                self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_order_is_the_one_the_tables_give(self):
        order = {r["cause"]: r for r in rows("page-d0-order.tsv")}
        self.assertEqual(sorted(order), ["absent", "chosen", "lineage", "removed", "unused"])
        self.assertEqual({r["miss_pairs"] for r in order.values()}, {"177"})
        pairs = {c: int(r["miss_pairs_in_cause"]) for c, r in order.items()}
        self.assertEqual(pairs, {"absent": 107, "removed": 15, "lineage": 1, "chosen": 28, "unused": 33})
        self.assertEqual(sum(pairs.values()), 184)
        self.assertEqual([c for c, r in order.items() if r["sets_the_order"] == "yes"], ["absent"])
        self.assertEqual(max(FOUR, key=pairs.get), "absent")
        self.assertGreaterEqual(pairs["absent"] / 177, 0.40)
        # Ruling 6 on the plan's 204 miss pairs, with the 27 unread pairs given to the next cause.
        nxt = max(pairs[c] for c in FOUR if c != "absent")
        self.assertEqual(nxt, 33)
        self.assertGreaterEqual(pairs["absent"] / 204, 0.40)
        self.assertGreater(pairs["absent"], nxt + 27)
        self.assertIn("It holds 107 of the 177 miss pairs (0.605).", self.text)
        self.assertIn("the five causes hold 184 entries for 177 pairs", self.text)
        self.assertIn("107 is 0.40 or more of 204 and is above 33 plus 27. So the first block is stable.", self.text)
        self.assertEqual(text("d0/d0-ordering.txt").split(":")[0], "absent holds 107 of 177 miss pairs (0.605)")

    def test_the_fixed_words_are_the_summary_s(self):
        summary = text("summary.txt").split("\n")
        heads = ("D0 ordering", "order of the input group", "coverage:", "ruling 6:", "removed, as counted")
        fixed = [line for line in summary if line.startswith(heads)]
        self.assertEqual(len(fixed), 5)
        for line in fixed:
            self.assertEqual(self.lines.count(line), 1, line)
        self.assertIn("order of the input group: C3, then C4 and C5; then C7", fixed[1])
        self.assertTrue(fixed[3].endswith("the first block is stable"))
        gates = [line for line in text("gates.txt").split("\n") if line]
        self.assertEqual(len(gates), 5)
        for line in gates:
            self.assertIn(f"gate passed: {line}", summary)
        for name in ("summary.txt", "gates.txt"):
            self.assertIn("the label count's datasets.tsv that are not M", text(name))

    def test_each_learner_s_pairs_and_the_datasets_file_agree(self):
        learners = {r["learner"]: r for r in rows("page-d0-learners.tsv")}
        pairs = rows("page-d0-datasets.tsv")
        self.assertEqual(len(pairs), 177)
        self.assertEqual(len({r["dataset"] for r in pairs}), 65)
        self.assertEqual(len({r["species"] for r in pairs}), 16)
        self.assertEqual(len({(r["dataset"], r["learner"]) for r in pairs}), 177)
        want = {"scm": 23, "cart": 28, "rf": 52, "gbt": 41, "lasso": 33}
        for name in LEARNERS:
            mine = [r for r in pairs if r["learner"] == name]
            self.assertEqual((len(mine), int(learners[name]["miss_pairs"])), (want[name], want[name]))
            for cause in page_tables.CAUSES:
                held = sum(1 for r in mine if cause in causes_of(r))
                self.assertEqual(held, int(learners[name][f"datasets_{cause}"]), (name, cause))
            self.assertEqual(sum(1 for r in mine if "+" in r["cause"]), int(learners[name]["ties"]))
            largest = max(page_tables.CAUSES, key=lambda c: int(learners[name][f"datasets_{c}"]))
            self.assertEqual(largest, "absent", name)
        self.assertIn("Absent is the largest cause for each of the five learners.", self.text)
        self.assertIn(
            "They give 177 miss pairs: 23 for `scm`, 28 for `cart`, 52 for `rf`, 41 for `gbt` and 33 for `lasso`.",
            self.text,
        )

    def test_what_was_seen_after_is_what_the_datasets_file_holds(self):
        pairs = rows("page-d0-datasets.tsv")
        alone = collections.Counter(r["cause"] for r in pairs if "+" not in r["cause"])
        self.assertEqual(
            dict(alone),
            {"absent": 104, "unused": 27, "chosen": 23, "removed": 14, "lineage": 1, "no missed fold": 1},
        )
        self.assertEqual(sum(1 for r in pairs if "+" in r["cause"]), 7)
        self.assertIn(
            "Absent is the one cause of 104 of the 177 pairs. Unused is the one cause of 27 pairs, "
            "Chosen of 23 pairs, Removed of 14 pairs and Lineage of 1 pair. Seven pairs tie and 1 pair has no missed fold:",
            self.text,
        )
        absent = collections.Counter(r["species"] for r in pairs if "absent" in causes_of(r))
        self.assertEqual(sum(absent.values()), 107)
        self.assertEqual((absent["Salmonella enterica"], absent["Escherichia coli"]), (40, 20))
        five = ["Salmonella enterica", "Escherichia coli", "Staphylococcus aureus", "Enterobacter cloacae",
                "Klebsiella pneumoniae"]
        self.assertEqual([name for name, _ in absent.most_common(5)], five)
        self.assertEqual(sum(absent[name] for name in five), 96)
        self.assertEqual(len({r["species"] for r in pairs}) - len(absent), 7)
        removed = collections.Counter(r["learner"] for r in pairs if "removed" in causes_of(r))
        self.assertEqual((sum(removed.values()), sorted(removed)), (15, ["gbt", "lasso", "rf"]))
        lineage = [r for r in pairs if "lineage" in causes_of(r)]
        self.assertEqual(
            [(r["species"], r["antibiotic"], r["learner"]) for r in lineage], [("Shigella sonnei", "colistin", "rf")]
        )
        chosen = [r for r in pairs if "chosen" in causes_of(r)]
        self.assertEqual((len(chosen), sum(1 for r in chosen if r["learner"] == "rf")), (28, 17))
        coli = [r for r in pairs if r["species"] == "Campylobacter coli"]
        self.assertEqual((len(coli), sum(1 for r in coli if "chosen" in causes_of(r))), (7, 7))
        unused = collections.Counter(r["learner"] for r in pairs if "unused" in causes_of(r))
        self.assertEqual((unused["scm"], unused["cart"], unused["rf"]), (8, 10, 3))

    def test_the_overlaps_of_the_reading_are_the_summary_table_s(self):
        summary = {(r["scheme"], r["learner"]): r for r in rows("page-c8-summary.tsv")}
        self.assertEqual(len(summary), 12)

        def medians(scheme):
            return [page_tables.d3(summary[(scheme, name)]["median"]) for name in LEARNERS]

        self.assertEqual(medians("lineage"), ["0.050", "0.053", "0.183", "0.075", "0.066"])
        self.assertEqual(medians("random"), ["0.155", "0.125", "0.250", "0.123", "0.126"])
        for name in LEARNERS:
            self.assertLess(float(summary[("lineage", name)]["median"]), float(summary[("random", name)]["median"]))
        self.assertIn(
            "The median overlap is 0.050 for `scm`, 0.053 for `cart`, 0.183 for `rf`, 0.075 for `gbt` and 0.066 for `lasso`. "
            "On the random split each is higher: 0.155, 0.125, 0.250, 0.123 and 0.126.",
            self.text,
        )
        shares = {name: page_tables.d3(summary[("lineage", name)]["catalog_family_share_of_carded"]) for name in LEARNERS}
        self.assertEqual((min(shares, key=shares.get), max(shares, key=shares.get)), ("lasso", "scm"))
        self.assertEqual((shares["lasso"], shares["scm"]), ("0.101", "0.259"))
        self.assertIn("between 0.101 (`lasso`) and 0.259 (`scm`)", self.text)
        five = summary[("lineage", "over_five")]
        self.assertEqual((five["with_card_summed"], five["distinct_summed"]), ("20580", "105166"))
        self.assertIn("20,580 of the 105,166 distinct fold-chosen features have a card", self.text)
        for scheme in ("lineage", "random"):
            for name in ("scm", "cart"):
                self.assertEqual(page_tables.d3(summary[(scheme, name)]["q1"]), "0.000")

    def test_the_units_named_on_the_page_are_the_units_file_s(self):
        units = rows("page-c8-datasets.tsv")
        self.assertEqual(len(units), 238 * 5 * 2)
        self.assertEqual(len({(r["dataset"], r["learner"], r["scheme"]) for r in units}), len(units))
        undefined = sorted({(r["species"], r["antibiotic"]) for r in units if r["defined"] != "yes"})
        self.assertEqual(
            undefined,
            [("Enterococcus faecium", "chloramphenicol"), ("Salmonella enterica", "temocillin"),
             ("Staphylococcus aureus", "ceftaroline"), ("Staphylococcus aureus", "levofloxacin"),
             ("Staphylococcus aureus", "moxifloxacin")],
        )
        self.assertEqual({r["scheme"] for r in units if r["defined"] != "yes"}, {"lineage"})
        self.assertEqual(sum(1 for r in units if r["defined"] != "yes"), 25)
        under = sorted(
            (r["learner"], r["species"], r["antibiotic"]) for r in units
            if r["defined"] == "yes" and int(r["folds_choosing"].split("-")[0]) < 2
        )
        self.assertEqual(
            under,
            [("gbt", "Acinetobacter baumannii", "cefepime"), ("gbt", "Escherichia coli", "amikacin"),
             ("gbt", "Neisseria gonorrhoeae", "ceftriaxone"), ("lasso", "Klebsiella pneumoniae", "cefotetan"),
             ("lasso", "Klebsiella pneumoniae", "doripenem")],
        )
        self.assertEqual({r["scheme"] for r in units if r["defined"] == "yes"
                          and int(r["folds_choosing"].split("-")[0]) < 2}, {"lineage"})
        for words in (
            "*Acinetobacter baumannii* with cefepime, *Escherichia coli* with amikacin and *Neisseria gonorrhoeae* "
            "with ceftriaxone for `gbt`, and *Klebsiella pneumoniae* with cefotetan and with doripenem for `lasso`",
            "*Enterococcus faecium* with chloramphenicol, *Salmonella enterica* with temocillin, and "
            "*Staphylococcus aureus* with ceftaroline, with levofloxacin and with moxifloxacin",
        ):
            self.assertIn(words, self.text)

    def test_the_species_table_adds_up_to_the_page_s_row(self):
        species = rows("page-species.tsv")
        self.assertEqual(len(species), 23)
        cols = ["datasets", "d0_datasets", "check_dataset_rows", "check_learner_rows",
                "breakpoints_set_to_the_release", "feature_lists_fetched",
                "samples_with_other_descriptive_fields", "cards_with_other_locus"]
        sums = [sum(int(r[c]) for r in species) for c in cols]
        self.assertEqual(sums, [238, 65, 238, 2790, 52, 991, 171103, 4])
        self.assertEqual(self.lines.count("| All | 238 | 65 | 238 | 2,790 | 52 | 991 | 171,103 | 4 |"), 1)
        self.assertEqual(sum(int(r["feature_lists_equal_to_master"]) for r in species), 0)
        self.assertEqual(sum(int(r["datasets_with_other_locus_cards"]) for r in species), 3)
        self.assertEqual(sum(1 for r in species if int(r["d0_datasets"]) > 0), 16)
        for r in species:
            self.assertEqual(len(r["antibiotics"].split("; ")), int(r["datasets"]), r["species"])
            d0 = [] if r["d0_antibiotics"] == "none" else r["d0_antibiotics"].split("; ")
            self.assertEqual(len(d0), int(r["d0_datasets"]), r["species"])
            self.assertLessEqual(set(d0), set(r["antibiotics"].split("; ")))
        pairs = rows("page-d0-datasets.tsv")
        for r in species:
            mine = sorted({p["antibiotic"] for p in pairs if p["species"] == r["species"]})
            d0 = [] if r["d0_antibiotics"] == "none" else r["d0_antibiotics"].split("; ")
            self.assertEqual(mine, sorted(d0), r["species"])

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in text("SHA256SUMS").splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        found = []
        for folder, _, names in os.walk(RECORD):
            for name in names:
                if name != "SHA256SUMS":
                    found.append(os.path.relpath(os.path.join(folder, name), RECORD))
        self.assertEqual(sorted(listed), sorted(found))
        for name, digest in listed.items():
            with open(os.path.join(RECORD, name), "rb") as handle:
                self.assertEqual(hashlib.sha256(handle.read()).hexdigest(), digest, name)


if __name__ == "__main__":
    unittest.main()
