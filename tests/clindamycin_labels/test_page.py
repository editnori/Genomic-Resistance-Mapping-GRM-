"""The bundled clindamycin kinds against the page docs/clindamycin-labels.md.

Every count of the bundled file is read back from the page's own table cell
or sentence for that dataset, so the page and the site cannot drift apart.
The kind of each dataset is read from the page three times (its opening
list, its table and its reading section), and the writer's kind rule,
applied to the page's own cells, must give it. No count is typed here.
"""

import re
import unittest

from . import common

S = common.load_site_data()

KIND_WORDS = {
    "measured dependence": "measured_dependence",
    "one study": "one_study",
    "class balance, dependence not measured": "class_balance",
    "no dependence under the rule": "no_dependence",
}

# The integer fields of a dataset in the bundle; each one is checked against the page below.
COUNTS = (
    "samples",
    "study_samples",
    "disagreements",
    "study_disagreements",
    "no_pmid_samples",
    "left_resistant",
    "left_susceptible",
    "rows_passing",
    "rows_tested",
    "distinct_passing",
    "distinct_tested",
)


def section(text, title):
    """The lines of one `## title` section of the page."""
    m = re.search(rf"^## {re.escape(title)}\n(.*?)(?=^## |\Z)", text, re.M | re.S)
    if m is None:
        raise AssertionError(f"the page has no section {title!r}")
    return m.group(1)


def table(text, first_header):
    """The rows of the Markdown table whose header row starts with `first_header`, as dicts."""
    lines = text.split("\n")
    start = next(i for i, line in enumerate(lines) if line.startswith(first_header))
    head = [c.strip() for c in lines[start].strip("|").split("|")]
    rows = []
    for line in lines[start + 2 :]:
        if not line.startswith("|"):
            break
        rows.append(dict(zip(head, (c.strip() for c in line.strip("|").split("|")))))
    return rows


def dataset_of(cell):
    m = re.fullmatch(r"`([a-z_]+__clindamycin)`", cell)
    if m is None:
        raise AssertionError(f"{cell!r} names no dataset")
    return m.group(1)


def of(cell):
    """'182 of 219, 83.1% (...)' -> (182, 219)."""
    m = re.match(r"(\d+) of (\d+)\b", cell)
    if m is None:
        raise AssertionError(f"{cell!r} is not 'x of n'")
    return int(m.group(1)), int(m.group(2))


def read_page():
    """What the page says of each dataset, from its cells and sentences."""
    text = common.PAGE.read_text(encoding="utf-8")
    six = section(text, "The six datasets")
    page = {}
    for row in table(six, "| Dataset | PMID the plan picks |"):
        d = dataset_of(row["Dataset"])
        study_disagreements, disagreements = of(row["Its disagreements"])
        study_samples, samples = of(row["Its samples"])
        page[d] = {
            "pmid": row["PMID the plan picks"],
            "samples": samples,
            "study_samples": study_samples,
            "disagreements": disagreements,
            "study_disagreements": study_disagreements,
            "no_pmid_samples": int(row["Samples with no PMID"]),
        }
    for row in table(six, "| Dataset | Left without it, R / S |"):
        d = dataset_of(row["Dataset"])
        left = re.fullmatch(r"(\d+) / (\d+)", row["Left without it, R / S"])
        verdict = row["Plan's verdict"]
        tested = re.search(r"(\d+) of (\d+) rows pass \((\d+) of (\d+) distinct\)", verdict)
        existence = re.search(r"existence: (?:(\d+) reported (resistant|susceptible) left, under (\d+)|the study is the whole dataset)", verdict)
        page[d].update(
            {
                "left_resistant": int(left.group(1)),
                "left_susceptible": int(left.group(2)),
                "scored_without": {"yes": True, "no": False}[row["Scored without it"]],
                "rows_passing": None if tested is None else int(tested.group(1)),
                "rows_tested": None if tested is None else int(tested.group(2)),
                "distinct_passing": None if tested is None else int(tested.group(3)),
                "distinct_tested": None if tested is None else int(tested.group(4)),
                "kind_in_table": KIND_WORDS[row["Kind"].lower()],
                "existence": None if existence is None else existence.groups(),
            }
        )
    reading = section(text, "What a released clindamycin figure can be read as")
    for m in re.finditer(r"^- \*\*`([a-z_]+__clindamycin)`: ([^*]+)\.\*\* (.*)$", reading, re.M):
        moves = set(re.findall(r"\b(rise|fall) by (0\.\d+) or more", m.group(3)))
        page[m.group(1)].update({"kind_in_reading": KIND_WORDS[m.group(2).lower()], "moves": moves})
    opening = text.split("\n## ", 1)[0]
    page_kinds = {}
    for m in re.finditer(r"^- \*\*([^*]+):\*\* (.*)$", opening, re.M):
        for d in re.findall(r"`([a-z_]+__clindamycin)`", m.group(2).split(".")[0]):
            page_kinds[d] = KIND_WORDS[m.group(1).lower()]
    for d, kind in page_kinds.items():
        page[d]["kind_in_opening"] = kind
    return text, page


class PageAndBundle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text, cls.page = read_page()
        cls.data = common.bundled()

    def test_the_page_and_the_bundle_hold_the_same_six_datasets(self):
        self.assertEqual(len(self.page), 6)
        self.assertEqual(sorted(self.page), sorted(self.data["datasets"]))

    def test_the_page_names_each_kind_the_same_way_three_times(self):
        for d, p in self.page.items():
            self.assertEqual({p["kind_in_opening"], p["kind_in_table"], p["kind_in_reading"]}, {p["kind_in_table"]}, d)

    def test_the_writers_rule_gives_the_pages_kind_from_the_pages_own_cells(self):
        for d, p in self.page.items():
            self.assertEqual(S.kind_of(p["scored_without"], p["left_resistant"], p["left_susceptible"], p["distinct_passing"]), p["kind_in_table"], d)

    def test_each_kind_in_the_bundle_is_the_pages(self):
        for d, p in self.page.items():
            self.assertEqual(self.data["datasets"][d]["kind"], p["kind_in_table"], d)

    def test_every_count_of_the_bundle_is_the_pages_cell_for_that_dataset(self):
        for d, entry in self.data["datasets"].items():
            ints = {k for k, v in entry.items() if isinstance(v, int) and not isinstance(v, bool)}
            self.assertLessEqual(ints, set(COUNTS), f"{d}: a count the page check does not read")
            for field in COUNTS:
                self.assertEqual(entry[field], self.page[d][field], f"{d} {field}")
            self.assertEqual(entry["pmid"], self.page[d]["pmid"], d)
            self.assertEqual(entry["scored_without"], self.page[d]["scored_without"], d)

    def test_a_class_left_under_the_minimum_is_the_pages_class_and_count(self):
        for d, p in self.page.items():
            entry = self.data["datasets"][d]
            if p["existence"] is None:
                continue
            left, which, minimum = p["existence"]
            if left is None:
                # "the study is the whole dataset"
                self.assertEqual((entry["study_samples"], entry["left_resistant"], entry["left_susceptible"]), (entry["samples"], 0, 0), d)
                continue
            self.assertEqual(entry[f"left_{which}"], int(left), d)
            self.assertEqual(self.data["min_class"], int(minimum), d)
            self.assertLess(entry[f"left_{which}"], self.data["min_class"], d)

    def test_the_passing_changes_rise_or_fall_as_the_reading_says(self):
        words = {"rise": "rises", "fall": "falls"}
        for d, p in self.page.items():
            entry = self.data["datasets"][d]
            if entry["kind"] == "measured_dependence":
                self.assertEqual(len(p["moves"]), 1, d)
                ((move, threshold),) = p["moves"]
                self.assertEqual(entry["direction"], words[move], d)
                self.assertEqual(self.data["threshold"], float(threshold), d)
            else:
                self.assertIsNone(entry["direction"], d)

    def test_the_caveat_rules_threshold_and_minimum_are_the_pages(self):
        rule = section(self.text, "The plan and the release read")
        self.assertIn(f"leaves a class under {self.data['min_class']} samples", rule)
        self.assertIn(f"is at least {self.data['threshold']} in size", rule)

    def test_the_page_counts_the_datasets_of_each_kind_as_the_bundle_does(self):
        counted = {k: 0 for k in S.KINDS}
        for entry in self.data["datasets"].values():
            counted[entry["kind"]] += 1
        from_page = {k: 0 for k in S.KINDS}
        for p in self.page.values():
            from_page[p["kind_in_opening"]] += 1
        self.assertEqual(counted, from_page)

    def test_the_release_is_the_pages(self):
        self.assertIn(f"The release `{self.data['release']}` (evaluation revision {self.data['evaluation_revision']})", self.text)
        self.assertEqual(self.data["page"], "docs/clindamycin-labels.md")
        self.assertTrue((common.ROOT / self.data["reference"]).is_dir())


if __name__ == "__main__":
    unittest.main()
