"""docs/mechanism-split-result.md against the committed tables of reference/mechanism-split/.

This is the one test file of the suite that reads committed result tables:
the aggregates the table script wrote from the one reading, by dataset,
species, cell and learner. It fits nothing and reads no genome.
"""

import csv
import hashlib
import re
import unittest

from tests.mechanism import common

page_tables = common.load("page_tables")
page_faults = common.load("page_faults")

RECORD = common.ROOT / "reference" / "mechanism-split"
PAGE = common.ROOT / "docs" / "mechanism-split-result.md"
LISTS = common.ROOT / "scripts" / "c3-inputs" / "lists.tsv"
SECTIONS = ["Reading", "Data", "M1", "M2", "Guards", "Seen after", "Limits", "Method", "Datasets"]
TABLES = {"Reading": ["m1_reading"], "M1": ["m1_beside", "m1_neither", "m1_strata"], "M2": ["m2_cells"],
          "Guards": ["m2_summary", "m2_secondary"], "Datasets": ["species", "datasets"]}
LEARNERS = page_tables.LEARNERS
d3, of, pval = page_tables.d3, page_tables.of, page_tables.pval


def rows(name, path=None):
    with open(path or RECORD / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def blocks(text):
    return [b for b in re.split(r"\n\s*\n", text) if b.strip()]


def count(n):
    return f"{int(n):,}"


def span(values):
    values = sorted(values)
    return values[0], values[-1]


def m1(split="lineage", set_="primary"):
    return {r["row"]: r for r in rows("combined/m1-reading.tsv") if (r["split"], r["set"]) == (split, set_)}


def side(x):
    x = float(x)
    return "gene" if x > 0.02 else "point" if x < -0.02 else "same"


def reading_datasets(set_="primary"):
    """The datasets of the M1 reading: lineage split, both strata at 10 or more, per row."""
    out = {}
    for r in rows("page-m1-datasets.tsv"):
        if r["split"] == "lineage" and r["set"] == set_ and r["gene_less_point"] != "":
            out.setdefault((r["species"], r["antibiotic"], r["drug_class"]), {})[r["row"]] = r
    return out


def deciding_cells(periods="primary"):
    return [r for r in rows("page-m2-cells.tsv")
            if (r["cell"], r["periods"], r["guard"]) == ("species_class", periods, "bioproject")]


def summary(cell, periods, guard, statistic):
    found = [r for r in rows("combined/m2-summary.tsv")
             if (r["cell"], r["periods"], r["guard"], r["statistic"]) == (cell, periods, guard, statistic)]
    assert len(found) == 1, (cell, periods, guard, statistic)
    return found[0]


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.lines = self.text.split("\n")
        self.sections = sections(self.text)

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_every_printed_table_stands_whole_in_its_section(self):
        tables = page_tables.tables(RECORD)
        self.assertEqual(sorted(tables), sorted(t for names in TABLES.values() for t in names))
        for section, names in TABLES.items():
            for name in names:
                printed = tables[name]
                self.assertGreater(len(printed), 2, name)
                self.assertIn("\n".join(printed) + "\n", self.sections[section], name)
                for line in printed[:1] + printed[2:]:
                    self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_m1_reading_stands_in_plain_words(self):
        read = m1()
        learners = [read[x] for x in LEARNERS]
        self.assertEqual({(r["datasets"], r["species"]) for r in learners}, {("19", "5")})
        self.assertEqual({(r["lo"], r["hi"]) for r in learners}, {("", "")})
        self.assertEqual({r["reading"] for r in learners}, {"no difference shown"})
        self.assertEqual({pval(r["holm_p"]) for r in learners}, {"0.318"})
        gene, point = span(int(r["gene_ahead"]) for r in learners), span(int(r["point_ahead"]) for r in learners)
        lo, hi = span(float(r["median"]) for r in learners)
        reading = self.sections["Reading"]
        self.assertIn("On the lineage split of the primary set, 19 of the 204 datasets, in 5 species, hold 10 or "
                      "more resistant genomes in both.", reading)
        self.assertIn(f"each learner's recall is higher on gene only, by more than 0.02, on {gene[0]} or {gene[1]} "
                      f"datasets, and higher on point only on {point[0]} or {point[1]}, with a median difference of "
                      f"{d3(lo)} to {d3(hi)}.", reading)
        self.assertIn("With 5 species there is no species interval, since the method gives none under 6.", reading)
        self.assertIn("The sign test after Holm over the five learners gives 0.318 for each, so each learner reads "
                      "\"no difference shown\".", reading)
        (datasets, species), = {(r["datasets"], r["species"]) for r in learners}
        self.assertIn(f"**No difference shown.** On M1, \"no difference shown\" is the rule's reading on {datasets} "
                      f"datasets in {species} species and is no statement that recall is equal.", self.sections["Limits"])

    def test_the_two_callers_beside_stand_as_their_rows_give(self):
        read, ds = m1(), reading_datasets()
        cat, prog = read["catalog"], read["program"]
        self.assertEqual((cat["gene_ahead"], cat["point_ahead"], cat["same"], d3(cat["median"])), ("17", "0", "2", "0.984"))
        self.assertEqual((prog["gene_ahead"], prog["point_ahead"], prog["same"], d3(prog["median"])), ("0", "0", "19", "0.000"))
        cats = [d["catalog"] for d in ds.values()]
        self.assertEqual(len(cats), 19)
        self.assertEqual(sum(float(r["gene_only_recall"]) == 1 for r in cats), 17)
        self.assertEqual(d3(min(float(r["gene_only_recall"]) for r in cats)), "0.984")
        self.assertEqual(sum(float(r["point_only_recall"]) == 0 for r in cats), 10)
        progs = [d["program"] for d in ds.values()]
        self.assertTrue(all(float(r["gene_only_recall"]) == float(r["point_only_recall"]) == 1 for r in progs))
        reading = self.sections["Reading"]
        self.assertIn("Its recall is higher on gene only on 17 of the 19 datasets and the same on 2, with a median "
                      "difference of 0.984. Its recall on gene only is 1.000 on 17 of the 19 and 0.984 or more on "
                      "all 19. On point only it is 0.000 on 10 of the 19.", reading)
        self.assertIn("Its recall is 1.000 on both strata in each of the 19 datasets, so its difference is 0.000 on "
                      "19 of 19.", reading)

    def test_the_m1_datasets_are_named_in_full(self):
        def named(set_):
            by = {}
            for sp, ab, _ in sorted(reading_datasets(set_)):
                by.setdefault(sp, []).append(ab)
            return by
        primary, in_m = named("primary"), named("M")
        self.assertEqual(sum(map(len, primary.values())), 19)
        self.assertEqual(sum(map(len, in_m.values())), 5)
        data = self.sections["Data"]

        def words(names):
            return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]
        for sp, abs_ in primary.items():
            self.assertIn(f"*{sp}* with {words(abs_)}", data)
        for sp, abs_ in in_m.items():
            self.assertIn(f"*{sp}* with {words(abs_)}", data)
        most = max(primary, key=lambda sp: len(primary[sp]))
        sizes = sorted((len(v) for v in primary.values()), reverse=True)
        self.assertGreater(sizes[0], sizes[1])
        self.assertIn(f"The M1 reading stands on {sum(map(len, primary.values()))} datasets in {len(primary)} species, "
                      f"{len(primary[most])} of them in *{most}*.", self.sections["Limits"])

    def test_beside_neither_and_strata_stand_as_their_rows_give(self):
        text = self.sections["M1"]
        for (split, set_), n in ((("lineage", "M"), 5), (("random", "primary"), 19)):
            r = [m1(split, set_)[x] for x in LEARNERS]
            self.assertEqual({x["datasets"] for x in r}, {str(n)})
            lo, hi = span(float(x["median"]) for x in r)
            plo, phi = span(float(x["sign_p"]) for x in r)
            glo, ghi = span(int(x["gene_ahead"]) for x in r)
            if split == "lineage":
                self.assertEqual((glo, ghi, pval(plo), pval(phi)), (5, 5, "0.062", "0.062"))
                self.assertIn(f"each learner is gene ahead on 5 of 5, with a median of {d3(lo)} to {d3(hi)} and a "
                              f"sign p of 0.062.", text)
            else:
                self.assertIn(f"each learner is gene ahead on {glo} or {ghi} of the 19 primary datasets, with a "
                              f"median of {d3(lo)} to {d3(hi)} and a sign p of {pval(plo)} to {pval(phi)}.", text)
        neither = {r["row"]: r for r in rows("combined/m1-neither.tsv")
                   if (r["split"], r["set"]) == ("lineage", "primary")}
        lo, hi = span(float(neither[x]["median_recall"]) for x in LEARNERS)
        self.assertEqual({neither[x]["datasets"] for x in neither}, {"92"})
        self.assertIn(f"Over the 92 primary datasets with 10 or more resistant genomes in neither stratum", text)
        self.assertIn(f"the learners' median recall on the lineage split is {d3(lo)} to {d3(hi)}. The median recall "
                      f"of GRM's known-gene rule is {d3(neither['catalog']['median_recall'])}, and that of the "
                      f"AMRFinderPlus program {d3(neither['program']['median_recall'])}.", text)
        st = {r["stratum"]: r for r in rows("page-m1-strata.tsv")
              if (r["set"], r["split"], r["row"]) == ("primary", "lineage", "program")}
        self.assertEqual(sum(int(st[s]["resistant"]) for s in page_tables.STRATA), int(st["all"]["resistant"]))
        self.assertIn(f"Of the {count(st['all']['resistant'])} resistant genomes of the primary set, "
                      f"{count(st['gene_only']['resistant'])} are gene only, {count(st['point_only']['resistant'])} "
                      f"point only, {count(st['both']['resistant'])} both and {count(st['neither']['resistant'])} "
                      f"neither. Point only holds 10 or more resistant genomes in {st['point_only']['ten_or_more']} of "
                      f"the 204 datasets and gene only in {st['gene_only']['ten_or_more']}.", text)

    def test_the_m2_reading_stands_as_the_deciding_cells_give(self):
        cells = deciding_cells()
        passing = [r for r in cells if r["passes"] == "1"]
        point, tvd = summary("species_class", "primary", "bioproject", "point"), \
            summary("species_class", "primary", "bioproject", "tvd")
        self.assertEqual((len(cells), len(passing), point["pass"], point["move"], tvd["move"]),
                         (285, 110, "110", "0", "0"))
        mixed = [r for r in passing if r["p_tvd"] != ""]
        self.assertEqual(len(mixed), int(tvd["pass"]))
        species = sorted({r["species"] for r in passing})
        self.assertEqual(len(species), 6)
        low = d3(min(float(r["holm_p_point"]) for r in passing))
        self.assertEqual({float(r["holm_p_tvd"]) for r in mixed}, {1.0})
        sec = summary("species_class", "secondary", "bioproject", "point")
        reading = self.sections["Reading"]
        names = ", ".join(f"*{s}*" for s in species[:-1]) + f" and *{species[-1]}*"
        self.assertIn(f"110 of the 285 species and class cells pass the vote, all in 6 of the 21 species: {names}. "
                      f"No cell moves. On the point share 0 of the 110 move, and the smallest Holm p is {low}. On the "
                      f"family mix 0 of the {len(mixed)} passing cells with a mix in both periods move, every Holm p "
                      f"1.000. On 2000 to 2009 against 2010 to 2019, {sec['pass']} of the 285 pass and none moves.",
                      reading)
        text = self.sections["M2"]
        per = {}
        for r in cells:
            per.setdefault(r["species"], [0, 0])
            per[r["species"]][0] += 1
            per[r["species"]][1] += r["passes"] == "1"
        words = ", ".join(f"*{s}* {per[s][1]} of {per[s][0]}" for s in species[:-1])
        self.assertIn(f"On the primary periods 110 of the 285 cells pass, in 6 species: {words} and *{species[-1]}* "
                      f"{per[species[-1]][1]} of {per[species[-1]][0]}. No cell of the other 15 species passes.", text)
        self.assertTrue(all(per[s][0] == per[s][1] for s in species))
        sec_cells = [r for r in deciding_cells("secondary") if r["passes"] == "1"]
        sec_species = sorted({r["species"] for r in sec_cells})
        self.assertEqual(len(sec_cells), int(sec["pass"]))
        self.assertIn(f"{sec['pass']} of the 285 cells pass, in {len(sec_species)} species: "
                      + ", ".join(f"*{s}*" for s in sec_species[:-1]) + f" and *{sec_species[-1]}*. "
                      f"{summary('species_class', 'secondary', 'bioproject', 'tvd')['pass']} of the {sec['pass']} "
                      "have a family mix in both periods, and none moves on either statistic.", text)

    def test_no_guard_moves_a_cell_and_the_guards_sentences_stand(self):
        all_rows = rows("combined/m2-summary.tsv")
        self.assertEqual({r["move"] for r in all_rows}, {"0"})
        voted = [r for r in all_rows if r["cell"] == "species_class" and r["periods"] == "primary"
                 and r["guard"] not in ("bioproject", "pooled")]
        low, high = min(voted, key=lambda r: int(r["pass"])), max(voted, key=lambda r: int(r["pass"]))
        self.assertEqual((low["pass"], low["guard"], low["statistic"]), ("34", "flagged_out", "tvd"))
        self.assertEqual((high["pass"], high["guard"], high["statistic"]), ("145", "lineage", "point"))
        self.assertEqual({r["pass"] for r in all_rows if r["guard"] == "pooled"}, {"0"})
        guards = self.sections["Guards"]
        self.assertIn("that runs from 34 cells with a mix, with flagged BioProjects left out, to 145 cells on the "
                      "point share, with one vote per lineage cluster.", guards)
        ds = [summary("dataset_resistant", p, "bioproject", s)["pass"]
              for p in ("primary", "secondary") for s in ("point", "tvd")]
        self.assertEqual(summary("dataset_resistant", "primary", "bioproject", "point")["cells"], "204")
        self.assertIn(f"On the primary periods {ds[0]} dataset cells pass on the point share and {ds[1]} on the mix, "
                      f"on the secondary periods {ds[2]} and {ds[3]}, and none moves.", guards)

    def test_the_time_split_link_stands_beside_recall_and_log_loss(self):
        link = {(r["learner"], r["measure"]): r for r in rows("combined/timesplit-link.tsv")}
        self.assertEqual({r["datasets"] for r in link.values()}, {"24"})
        guards = self.sections["Guards"]
        for x in LEARNERS:
            line = (f"| {page_tables.LEARNER_WORDS[x]} | 24 | {d3(link[(x, 'd_recall')]['spearman'])} | "
                    f"{d3(link[(x, 'd_log_loss')]['spearman'])} |")
            self.assertEqual(self.lines.count(line), 1, line)
            self.assertIn(line, guards)
        for measure, words in (("d_recall", "On recall"), ("d_log_loss", "on log loss")):
            low = min(LEARNERS, key=lambda x: float(link[(x, measure)]["spearman"]))
            high = max(LEARNERS, key=lambda x: float(link[(x, measure)]["spearman"]))
            self.assertIn(f"{words} it runs from {d3(link[(low, measure)]['spearman'])} (`{low}`) to "
                          f"{d3(link[(high, measure)]['spearman'])} (`{high}`)", guards)

    def test_the_seen_after_figures_follow_from_the_tables(self):
        seen = self.sections["Seen after"]
        ds = reading_datasets()
        quin = sorted(k for k in ds if k[2] == "quinolone")
        rest = [k for k in ds if k[2] != "quinolone"]
        self.assertEqual((len(quin), len(rest)), (5, 14))
        for k in quin:
            self.assertEqual({side(ds[k][x]["gene_less_point"]) for x in LEARNERS}, {"point"})
            self.assertEqual((float(ds[k]["catalog"]["gene_only_recall"]), float(ds[k]["catalog"]["point_only_recall"])),
                             (1.0, 0.0))
        sides = {x: [side(ds[k][x]["gene_less_point"]) for k in rest] for x in LEARNERS}
        for x in ("scm", "gbt", "lasso"):
            self.assertEqual(sides[x].count("gene"), 14, x)
        odd = [k for k in rest if side(ds[k]["cart"]["gene_less_point"]) != "gene"]
        self.assertEqual(odd, [k for k in rest if side(ds[k]["rf"]["gene_less_point"]) != "gene"])
        self.assertEqual([k[:2] for k in odd], [("Klebsiella pneumoniae", "doripenem")])
        self.assertEqual((sides["cart"].count("point"), sides["rf"].count("same")), (1, 1))
        dor = ds[odd[0]]
        self.assertIn(f"*Klebsiella pneumoniae* with doripenem ({d3(dor['cart']['gene_only_recall'])} against "
                      f"{d3(dor['cart']['point_only_recall'])}). `rf` is gene ahead on 13 and the same on that one "
                      f"({d3(dor['rf']['gene_only_recall'])} against {d3(dor['rf']['point_only_recall'])}).", seen)
        self.assertIn("*Klebsiella pneumoniae* with ciprofloxacin and levofloxacin, *Salmonella enterica* with "
                      "nalidixic acid, and *Shigella sonnei* with ciprofloxacin and nalidixic acid", seen)
        self.assertIn(f"Of the {len(ds)} datasets of the M1 reading, {len(quin)} have the drug class quinolone", seen)
        self.assertIn(f"On the {len(rest)} others,", seen)
        self.assertIn(f"`scm`, `gbt` and `lasso` are gene ahead on {sides['scm'].count('gene')} of {len(rest)}. `cart` "
                      f"is gene ahead on {sides['cart'].count('gene')} and point ahead on {sides['cart'].count('point')},",
                      seen)
        self.assertEqual(sorted({k[2] for k in rest}), sorted(["cephalosporin (1st-2nd gen)", "cephalosporin (3rd gen)",
                                                               "carbapenem", "penicillin + inhibitor", "aminoglycoside"]))

    def test_the_two_part_rule_stands_as_the_program_rows_give(self):
        rules = rows(None, common.KIT / "drug_rules.tsv")
        two = [r for r in rules if "&" in r["rule"]]
        self.assertEqual({r["antibiotic"] for r in two}, {"sulfamethoxazole/trimethoprim"})
        self.assertEqual({r["rule"] for r in two}, {"C:TRIMETHOPRIM & C:SULFONAMIDE"})
        prog = [r for r in rows("page-m1-datasets.tsv") if r["split"] == "lineage" and r["row"] == "program"]
        strata = ("gene_only", "point_only", "both")

        def short(r):
            return any(r[s + "_resistant"] not in ("", "0") and r[s + "_resistant"] != r[s + "_called"] for s in strata)
        smx = [r for r in prog if r["antibiotic"] == "sulfamethoxazole/trimethoprim"]
        self.assertEqual((len(smx), sum(map(short, smx))), (len(two), 9))
        self.assertEqual(sum(short(r) for r in prog if r["antibiotic"] != "sulfamethoxazole/trimethoprim"), 0)
        st = {r["stratum"]: r for r in rows("page-m1-strata.tsv")
              if (r["set"], r["split"], r["row"]) == ("primary", "lineage", "program")}
        gap = {s: sum(int(r[s + "_resistant"]) - int(r[s + "_called"]) for r in smx if r["set"] == "primary")
               for s in ("gene_only", "point_only")}
        self.assertEqual(gap["gene_only"], int(st["gene_only"]["resistant"]) - int(st["gene_only"]["resistant_called"]))
        self.assertEqual(gap["point_only"], int(st["point_only"]["resistant"]) - int(st["point_only"]["resistant_called"]))
        seen = self.sections["Seen after"]
        self.assertIn(f"So on 9 of the {len(smx)} sulfamethoxazole/trimethoprim datasets", seen)
        self.assertIn(f"it calls {count(st['gene_only']['resistant_called'])} of the "
                      f"{count(st['gene_only']['resistant'])} gene-only resistant genomes and "
                      f"{count(st['point_only']['resistant_called'])} of the {count(st['point_only']['resistant'])} "
                      f"point-only ones, and the {count(gap['gene_only'])} and {gap['point_only']} it does not call",
                      seen)
        keys = set(reading_datasets()) | set(reading_datasets("M"))
        self.assertFalse([k for k in keys if k[1] == "sulfamethoxazole/trimethoprim"])
        for r in smx:
            line = (f"| *{r['species']}* | {r['set']} | {of(r['gene_only_called'], r['gene_only_resistant'])} | "
                    f"{of(r['point_only_called'], r['point_only_resistant'])} |")
            self.assertEqual(self.lines.count(line), 1, line)

    def test_fusidic_acid_and_the_p_before_holm_stand_as_the_cells_give(self):
        seen = self.sections["Seen after"]
        names = sorted({r["class"] for r in rows("page-m2-cells.tsv") if r["cell"] == "species_class"})
        self.assertEqual(len(names), 32)
        folded = {}
        for n in names:
            folded.setdefault(n.replace("_", " "), []).append(n)
        self.assertEqual([v for v in folded.values() if len(v) > 1], [["FUSIDIC ACID", "FUSIDIC_ACID"]])
        fus = [r for r in rows("page-m2-cells.tsv") if r["guard"] == "bioproject" and r["cell"] == "species_class"
               and r["class"] in ("FUSIDIC ACID", "FUSIDIC_ACID")]
        self.assertEqual({r["species"] for r in fus}, {"Staphylococcus aureus"})
        for r in fus:
            share = ("earlier_point", "later_point") if r["class"] == "FUSIDIC ACID" else ("earlier_gene", "later_gene")
            self.assertEqual({float(r[c]) for c in share}, {0.0}, r["periods"])
        self.assertIn("No other two of the 32 class names differ only by a space or an underscore.", seen)
        passing = [r for r in deciding_cells() if r["passes"] == "1"]
        mixed = [r for r in passing if r["p_tvd"] != ""]
        low = min(passing, key=lambda r: float(r["p_point"]))
        self.assertIn(f"On the primary periods, {sum(float(r['p_point']) < 0.05 for r in passing)} of the 110 point "
                      f"p values and {sum(float(r['p_tvd']) < 0.05 for r in mixed)} of the {len(mixed)} mix p values "
                      f"are under 0.05 before Holm. The smallest point p is {pval(low['p_point'])}, *{low['species']}* "
                      f"with `{low['class']}` (a point share of {d3(low['earlier_point'])} against "
                      f"{d3(low['later_point'])}, {low['earlier_units']} and {low['later_units']} voting BioProjects), "
                      f"and its Holm p is {pval(low['holm_p_point'])}.", seen)

    def test_the_data_section_stands_as_the_record_gives(self):
        data = self.sections["Data"]
        sp = rows("page-species.tsv")

        def total(k):
            return sum(int(r[k]) for r in sp)
        self.assertEqual((len(sp), total("datasets"), total("in_m"), total("entered"), total("hard_check_pass")),
                         (21, 237, 33, 237, 237))
        self.assertIn(f"2,844 of 2,844 learner and catalog rows", data)
        self.assertIn(f"A: the {total('datasets')} datasets of the {len(sp)} species with the program's tables "
                      f"(`scripts/c3-inputs/lists.tsv`), {total('datasets') - total('in_m')} in the primary set and "
                      f"{total('in_m')} in M, which wait for the maintainer's label rule. All {total('entered')} "
                      f"entered.", data)
        self.assertEqual((total("learner_rows_equal"), total("learner_rows")), (2844, 2844))
        ds = rows("page-datasets.tsv")
        for col in ("rule_check", "hard_check", "program_check", "record_check", "entered"):
            self.assertEqual(sum(r[col] == "1" for r in ds), 237, col)
        years = {}
        for r in rows("page-m2-years.tsv"):
            years[r["years"]] = years.get(r["years"], 0) + int(r["genomes"])
        self.assertEqual(total("m2_genomes"), sum(v for k, v in years.items() if k != "no table"))
        no_table = [r["species"] for r in sp if r["m2_no_table"] != "0"]
        self.assertEqual(no_table, ["Staphylococcus aureus"])
        self.assertIn(f"{count(total('m2_genomes'))} genomes with a program table appear in the release's reports of "
                      f"the 21 species, and {years['no table']} genomes of *Staphylococcus aureus* in a report have none. "
                      f"Of the {count(total('m2_genomes'))}, {count(years['2010 to 2014'])} were collected in 2010 to "
                      f"2014 and {count(years['2015 to 2019'])} in 2015 to 2019, {count(years['2000 to 2009'])} in 2000 "
                      f"to 2009, {count(years['before 2000'])} before 2000 and {count(years['2020 on'])} from 2020 on, "
                      f"and {count(years['no year'])} have no year.", data)
        lin = rows("page-m2-lineage.tsv")
        self.assertEqual(({r["given"] for r in lin}, {r["differing_genome_pairs"] for r in lin}, len(lin)),
                         ({"1"}, {"0"}, 21))
        self.assertIn("with 0 genome pairs differing, so the lineage guard is given in 21 of 21.", data)

    def test_the_datasets_are_those_of_a_and_named_in_full(self):
        lrows = rows(None, LISTS)
        a = [r for r in lrows if r["in_A"] == "1"]
        self.assertEqual((len(a), sum(r["in_M"] == "1" for r in a)), (237, 33))
        named = {(r["species"], r["antibiotic"], r["set"]) for r in rows("page-datasets.tsv")}
        self.assertEqual(named, {(r["species"], r["antibiotic"], "M" if r["in_M"] == "1" else "primary") for r in a})
        datasets = self.sections["Datasets"]
        for sp, ab, set_ in named:
            self.assertIn(f"| *{sp}* | {ab} | {set_} |", datasets)

    def test_the_checks_and_the_reading_stand_as_their_files_give(self):
        checks = (RECORD / "checks.txt").read_text(encoding="utf-8")
        method = self.sections["Method"]
        for in_file, on_page in (
            ("checks: all pass", "All six checks passed."),
            ("species copies: 21; species of A in lists.tsv: 21; equal: 1", "| species | 21 copies, the 21 species of A |"),
            ("holds in 21 of 21 species, 210 of 210 lines, every file listed and none more",
             "holds in 21 of 21 species, 210 of 210 lines, every file listed and none more"),
            ("21 of 21 name their species, commit 11c299674c68, 10 results and the sha256 of their SHA256SUMS; FAILED "
             "files: 0", "21 of 21 name their species, commit `11c299674c68`, 10 results and the sha256 of their "
             "`SHA256SUMS`; no `FAILED` file"),
            ("19 first-run copies with DONE (commit 253d60e12759); 19 of 19 SHA256SUMS equal the second run's byte for "
             "byte; not in the first run: streptococcus_agalactiae, streptococcus_pneumoniae",
             "19 first-run copies with `DONE`, at commit `253d60e12759`; 19 of 19 `SHA256SUMS` equal the second run's "
             "byte for byte; *Streptococcus agalactiae* and *Streptococcus pneumoniae* not in the first run"),
            ("against the species' own ids: 0 files with a hit in 609 files (189 tables) over 21 species; scans "
             "failing: 0", "0 files with a hit in 609 files (189 tables) over 21 species"),
            ("against every species' ids: 0 files with a hit in 609 files (189 tables)",
             "against the 23,757 genome ids of the release's 23 species files"),
            ("23757 genome ids in the release's 23 species files", "the 23,757 genome ids of the release's 23"),
            ("237 of 237 datasets pass, 2844 of 2844 rows equal, 237 entered",
             "| hard check | 237 of 237 datasets pass, 2,844 of 2,844 rows equal, 237 entered |"),
            ("21 species read the program's tables made before and 0 ran the program again; 21 logged that nothing was "
             "installed (TABLES and DB given), 0 installed the program; 21 of 21 tables logs hold every table equal to "
             "its recorded sha256, 21625 tables in all; 21 of 21 database logs give 43 files equal to the kit's list",
             "21 species read the program's tables made before and 0 ran it again, and none installed it; 21,625 "
             "tables, each equal to its recorded sha256; 21 of 21 database logs give 43 files equal to the kit's list"),
        ):
            self.assertIn(in_file, checks)
            self.assertIn(on_page, method)
        reading = (RECORD / "reading.txt").read_text(encoding="utf-8").splitlines()
        printed = [line[len("# printed: "):] for line in reading if line.startswith("# printed: ")]
        self.assertEqual(printed, ["combine: 21 species runs; m1 rows 28, m2 rows 10752, time-split link rows 10"])
        self.assertIn(f"`combine.py` printed `{printed[0]}`.", method)
        self.assertIn("# runs: the 21 species of A at commit 11c299674c681581bbdaea0717058ed47034cd68", reading[1])
        self.assertIn("The reading ran once on 23 Nov 2026", self.text)
        self.assertEqual(reading[0], "# The mechanism split's one reading, 2026-11-23")
        for line in reading:
            if line.startswith("python3 -I scripts/mechanism/combine.py"):
                self.assertIn(line, method)
        gates = (RECORD / "gates.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(gates), len(page_tables.GATES))
        for n, line in enumerate(gates, 1):
            self.assertTrue(line.startswith(f"{n} "))
            name, words = line.split(" ", 1)[1].split(": ", 1)
            self.assertEqual(self.lines.count(f"| {n}. {name} | {words} |"), 1, line)
        self.assertEqual(len(gates), 6)
        self.assertIn("checked the reading against six gates and wrote `reference/mechanism-split/`. All six passed:",
                      method)

    def test_the_fault_table_names_each_planted_faults_test(self):
        method = self.sections["Method"]
        self.assertEqual([f[0] for f in page_faults.FAULTS], [f"P{n}" for n in range(1, 8)])
        for fault, words, test, _ in page_faults.FAULTS:
            self.assertEqual(self.lines.count(f"| {fault}. {words[0].upper()}{words[1:]} | "
                                              f"`{test.split('test_page_tables.')[1]}` |"), 1, fault)
        self.assertIn(f"All {len(page_faults.FAULTS)} were caught.", method)

    def test_the_page_gives_the_suites_own_test_counts(self):
        loader = unittest.TestLoader()
        here = loader.loadTestsFromName("tests.mechanism.test_page").countTestCases()
        tables = loader.loadTestsFromName("tests.mechanism.test_page_tables").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_page_tables.py` holds {tables} SYNTHETIC tests", method)
        self.assertIn(f"`tests/mechanism/test_page.py` holds {here} tests of this page", method)

    def test_the_words_keep_the_limits(self):
        for sentence in re.split(r"(?<=[.:])\s", self.text):
            if "external validation" in sentence:
                self.assertRegex(sentence, r"\b(never|not)\b", sentence)
            if "prevalence" in sentence:
                self.assertRegex(sentence, r"\b(no|never|not)\b", sentence)
        for block in blocks(self.text):
            if "log loss" in block:
                self.assertIn("recall", block, block[:80])
        self.assertNotIn("—", self.text)
        self.assertNotIn("<!--", self.text)
        self.assertNotRegex(self.text, r"\b(GCA|GCF)_\d|\bPRJ[A-Z]{2}\d")

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        found = sorted(
            str(p.relative_to(RECORD)) for p in RECORD.rglob("*") if p.is_file() and p.name != "SHA256SUMS"
        )
        self.assertEqual(sorted(listed), found)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)


if __name__ == "__main__":
    unittest.main()
