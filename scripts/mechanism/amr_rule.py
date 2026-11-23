"""The AMRFinderPlus comparison's rule, as that comparison applied it, and the
database's own family names.

The table reader, the variants and the partial methods are copied from the
specification's `amrfinder_adapter.py` (sha256 ADAPTER_SHA256), and
`token_match` and `explain` from its `grm_eval.py` (sha256 GRM_EVAL_SHA256),
which port `catalog::token_match` and `catalog::explain` (src/catalog.rs).
`run_species.sh` runs that adapter itself, unchanged, and `strata.py` stops
unless the calls made here equal the adapter's on every genome.

This reads the AMRFinderPlus program's output. GRM's known-gene rule is
GRM's own reader of the AMRFinderPlus catalog; its calls are in the reports
and never pass through here.

Python 3.12, standard library only; run scripts with -I.
"""

import csv
import os

ADAPTER_SHA256 = "4fda485048373a645385fa0a487f7d9a4897ab89557813dc65a8b1338183c833"
GRM_EVAL_SHA256 = "855002e74e1f596b7cb561caf70dfc92581134188b8f4f4ac88ac66321695003"

# Output column names by tool generation. 4.x renamed three of these.
COLUMNS = {
    "symbol": ("Element symbol", "Gene symbol"),
    "scope": ("Scope",),
    "type": ("Type", "Element type"),
    "subtype": ("Subtype", "Element subtype"),
    "class": ("Class",),
    "subclass": ("Subclass",),
    "method": ("Method",),
}
PARTIAL_METHODS = frozenset(("INTERNAL_STOP", "PARTIALX", "PARTIALP", "PARTIAL_CONTIG_ENDX", "PARTIAL_CONTIG_ENDP"))
GENE, POINT = frozenset(("AMR",)), frozenset(("AMR", "POINT", "POINT_DISRUPT"))
CORE, BOTH = frozenset(("core",)), frozenset(("core", "plus"))
# variant -> (subtypes counted, scopes counted, methods dropped)
VARIANTS = {
    "genes_all_scopes": (GENE, BOTH, frozenset()),
    "genes_core": (GENE, CORE, frozenset()),
    "tool_default": (POINT, CORE, frozenset()),
    "tool_default_intact": (POINT, CORE, PARTIAL_METHODS),
}
NEVER_COUNTED = frozenset(("AMR-SUSCEPTIBLE",))
# The subtypes read with point (the lead's ruling on question 1, 04:13Z on 23
# Nov 2026): a POINT row, and a POINT_DISRUPT row, a change in a core gene.
# Subtype AMR is gene. M1, M2 and the time-split link read them alike.
POINT_SUBTYPES = frozenset(("POINT", "POINT_DISRUPT"))
DISRUPT_SUBTYPE = "POINT_DISRUPT"


class Refusal(Exception):
    """An input breaks the contract. Nothing is scored after one."""


def refuse(message):
    raise Refusal(message)


def read_table(path):
    """Type AMR rows of one AMRFinderPlus table, as dicts with the adapter's
    field names. Refuses a file with no recognisable header."""
    if not os.path.isfile(path):
        # The message names no genome: logs leave the box.
        refuse("a genome of the cohort has no tool output; it was not run, or its table was not copied")
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        try:
            header = next(reader)
        except StopIteration:
            refuse("empty table; a genome with no hit still has the header line")
        index = {}
        for field, names in COLUMNS.items():
            found = [header.index(n) for n in names if n in header]
            if len(found) != 1:
                refuse(f"expected exactly one of the columns {names}; the header starts {header[:6]}")
            index[field] = found[0]
        rows = []
        for n, line in enumerate(reader, start=2):
            if len(line) != len(header):
                refuse(f"table line {n}: {len(line)} fields against {len(header)} in the header")
            row = {field: line[i] for field, i in index.items()}
            if row["type"] != "AMR" or row["subtype"] in NEVER_COUNTED:
                continue
            if row["subtype"] not in POINT:
                refuse(f"table line {n}: Type AMR row with subtype {row['subtype']!r}, which this adapter does not know")
            if row["scope"] not in BOTH:
                refuse(f"table line {n}: scope {row['scope']!r} is neither core nor plus")
            rows.append(row)
    return rows


def counted(rows, variant="tool_default"):
    """The rows a variant counts, whole."""
    subtypes, scopes, dropped = VARIANTS[variant]
    return [r for r in rows if r["subtype"] in subtypes and r["scope"] in scopes and r["method"] not in dropped]


def hits_for(rows, variant):
    return [{"group": r["symbol"], "class": r["class"], "subclass": r["subclass"]} for r in counted(rows, variant)]


def token_match(token, cls, subclass):
    kind, want = token[:2], token[2:]
    if kind not in ("C:", "S:"):
        refuse(f"rule token {token!r} is neither C: nor S:")
    field = cls if kind == "C:" else subclass
    return any(t.strip() == want for t in field.split("/"))


def explain(rule, hits):
    """(every rule group met, names of the hits used). `hits` are dicts
    with `group`, `class` and `subclass`."""
    used, met = set(), True
    for group in rule:
        matching = [h for h in hits if any(token_match(t, h["class"], h["subclass"]) for t in group)]
        if not matching:
            met = False
        used.update(h["group"] for h in matching)
    return met, sorted(used)


def program_call(rows, rule, variant="tool_default"):
    """The program's call for a dataset, as the comparison read it."""
    return explain(rule, hits_for(rows, variant))[0]


def rule_rows(rows, rule, variant="tool_default"):
    """The counted rows that meet any token of any group of the rule: the
    rows `explain` uses, whole."""
    tokens = [t for group in rule for t in group]
    return [r for r in counted(rows, variant) if any(token_match(t, r["class"], r["subclass"]) for t in tokens)]


def classes_of(row):
    """The class names of a row, as its Class field gives them split on '/'."""
    return [t.strip() for t in row["class"].split("/") if t.strip()]


def is_point(row):
    return row["subtype"] in POINT_SUBTYPES


def is_disrupt(row):
    return row["subtype"] == DISRUPT_SUBTYPE


def is_partial(row):
    return row["method"] in PARTIAL_METHODS


def parse_rule(text):
    """A rule as drug_rules.tsv writes it: groups joined by ' & ', tokens of a
    group by ' | '."""
    groups = [[t.strip() for t in g.split(" | ")] for g in text.split(" & ")]
    for g in groups:
        for t in g:
            if t[:2] not in ("C:", "S:") or len(t) < 3:
                refuse(f"rule token {t!r} is neither C: nor S:")
    return groups


def format_rule(rule):
    return " & ".join(" | ".join(g) for g in rule)


class Database:
    """The family names of the AMRFinderPlus database, read from its own
    ReferenceGeneCatalog.txt and ReferenceGeneHierarchy.txt."""

    def __init__(self, directory):
        cat = os.path.join(directory, "ReferenceGeneCatalog.txt")
        hier = os.path.join(directory, "ReferenceGeneHierarchy.txt")
        self.allele_family, self.families = {}, set()
        self.subclass_classes = {}
        with open(cat, encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                fam = (row.get("gene_family") or "").strip()
                if not fam:
                    continue
                self.families.add(fam)
                allele = (row.get("allele") or "").strip()
                if allele:
                    self.allele_family[allele] = fam
                classes = [t.strip() for t in (row.get("class") or "").split("/") if t.strip()]
                for s in (row.get("subclass") or "").split("/"):
                    s = s.strip()
                    if s:
                        self.subclass_classes.setdefault(s, set()).update(classes)
        self.parent = {}
        with open(hier, encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                node = (row.get("node_id") or "").strip()
                if node:
                    self.parent[node] = (row.get("parent_node_id") or "").strip()
        if not self.families:
            refuse("the database's catalogue names no gene_family")

    def family(self, symbol, subtype="AMR"):
        """The family of an element symbol, as the database names it: the
        allele's gene_family; else the symbol where it is a gene_family; else
        its first ancestor in the hierarchy that is a gene_family. Only a row
        of Subtype AMR enters a family mix: a POINT or POINT_DISRUPT row is
        read with point and has no family here (the lead's rulings on
        questions 1 and 10, 23 Nov 2026)."""
        if subtype != "AMR":
            refuse(f"a row of Subtype {subtype} is read with point and enters no family mix")
        if symbol in self.allele_family:
            return self.allele_family[symbol]
        if symbol in self.families:
            return symbol
        seen, node = set(), self.parent.get(symbol)
        while node and node not in seen:
            if node in self.families:
                return node
            seen.add(node)
            node = self.parent.get(node)
        # The message names no symbol: a mutation's name holds a position.
        refuse(f"an element symbol of Subtype {subtype}: neither the catalogue nor the hierarchy gives it a family")

    def rule_names_class(self, rule, cls):
        """Whether a published rule names a drug class: a C: token equal to it,
        or an S: token that the catalogue lists under it."""
        for group in rule:
            for t in group:
                if t[:2] == "C:" and t[2:] == cls:
                    return True
                if t[:2] == "S:" and cls in self.subclass_classes.get(t[2:], ()):
                    return True
        return False
