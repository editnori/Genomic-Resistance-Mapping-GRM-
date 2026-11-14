"""export.py [--check FILE]: GRM's entries for Givens, written from committed files.

Givens is a shared browser of the definitions and assumptions of four projects. This script writes GRM's part of it:
evaluation terms, the assumptions a report figure rests on, and aggregate counts. docs/givens-export.md describes the
file.

  python3 -I scripts/givens-export/export.py > reference/givens-export/grm.json
  python3 -I scripts/givens-export/export.py --check reference/givens-export/grm.json

With --check it prints nothing and exits 1 unless FILE is exactly the text it would write.

The words come from spec.json beside this script. Every number comes from one of two places:

  basis "file"      read here from a committed summary file. The export lists that file's sha256.
  basis "document"  quoted from a sentence of the repository's text. The script requires the sentence to be in the
                    document and the number to be in the sentence. No committed file rebuilds such a number unless the
                    line names the suite that does.

An entry GRM has no count for carries "unmeasured" with the reason. The script fills in nothing.

It is an allow-list. It stops with a message when the spec holds a key it does not know, reads a file or a field
outside the lists below, quotes a sentence its document does not hold, or writes a string that names a dataset, a
species, an antibiotic or a drug class. Standard library only.
"""

import hashlib
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SPEC = pathlib.Path(__file__).resolve().with_name("spec.json")
FORMAT = "grm-givens-export/1"

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLITS = ("random", "lineage")
KINDS = {"concept", "state", "role", "rule", "measure", "artifact"}
TRUTH_STATES = {"asserted", "negated", "false", "unknown", "absent", "missing", "unchecked", "unsupported", "inferred",
                "unresolved"}
EDGE_TYPES = {"is_a", "part_of", "value_of", "input_of", "output_of", "guards", "depends_on", "measured_by",
              "contrasts_with", "compared_with"}

# The committed files a number may be read from, and the fields of each. Nothing outside these is opened or read.
_L, _S = "|".join(LEARNERS), "|".join(SPLITS)
FIELDS = {
    "index-counts.json": re.compile(r"/(datasets|files|rows|learner_rows|flagged_rows|unflagged_with_no_call|"
                                    r"rates_compared|mismatches)"),
    "summary.json": re.compile(rf"/datasets|/per_learner|/per_learner/({_L})/(median_lineage|n_both)|/best/n_both|"
                               r"/no_lineage"),
    "guidance.json": re.compile(rf"/datasets|/rows/({_L}):({_S})/(n|fewer_than_every_genome)"),
}
SUMMARY_DIR = re.compile(r"reference/release-summary/rev(\d+)-[0-9a-z-]+/")
OUTSIDE = "web/src/data/outside-study.json"
DOCUMENT = re.compile(r"README\.md|docs/[a-z0-9-]+\.md")

MODULE_KEYS = ["id", "namespace", "name", "job", "repo", "release", "binding", "binding_note", "files", "nodes",
               "edges", "assumptions", "import_rule", "left_out"]
NODE_KEYS = {"id", "kind", "label", "definition", "truth_state", "quotes", "measured", "unmeasured"}
EDGE_KEYS = {"from", "type", "to"}
GIVEN_KEYS = {"id", "statement", "under", "hidden_in", "when_false", "quotes", "measured", "unmeasured"}
LINE_KEYS = {"what", "revision", "count", "of", "note", "document", "quote", "rebuilt_by", "values"}
OUTSIDE_ROWS = {
    "read": "that are read",
    "holds": "that hold",
    "sensitivity above": "whose sensitivity reads above",
    "sensitivity below": "whose sensitivity reads below",
}
OUTSIDE_TOTALS = {"datasets": "datasets", "species": "species", "frozen": "genomes frozen", "scored": "genomes scored"}
NUMBER = re.compile(r"\d[\d,]*(?:\.\d+)?")


class Fault(Exception):
    """The spec or a source file breaks one of the export's rules."""


def _keys(at, got, allowed, required=()):
    extra, missing = sorted(set(got) - set(allowed)), sorted(set(required) - set(got))
    if extra:
        raise Fault(f"{at}: the export has no field {', '.join(extra)}")
    if missing:
        raise Fault(f"{at}: {', '.join(missing)} is required")


class Sources:
    """Reads committed files under the allow-list and remembers which ones were read."""

    def __init__(self, root, aliases):
        self.root, self.aliases = pathlib.Path(root), aliases
        self.files, self.documents, self.used, self._json, self._text = {}, set(), set(), {}, {}

    def _bytes(self, path):
        full = self.root / path
        if not full.is_file():
            raise Fault(f"{path}: no such committed file")
        return full.read_bytes()

    def data(self, alias):
        if alias not in self.aliases:
            raise Fault(f"file {alias!r} is not named under files")
        path = self.aliases[alias]
        name = path.rsplit("/", 1)[-1]
        if not (path == OUTSIDE or (SUMMARY_DIR.fullmatch(path[: -len(name)]) and name in FIELDS)):
            raise Fault(f"{path}: the export may not read this file")
        if path not in self._json:
            raw = self._bytes(path)
            self.files[path] = hashlib.sha256(raw).hexdigest()
            self._json[path] = json.loads(raw.decode("utf-8"))
        self.used.add(alias)
        return path, self._json[path]

    def revision(self, alias):
        """The evaluation revision a data file names for itself."""
        path, doc = self.data(alias)
        value = doc.get("evaluation_revision", doc.get("revision"))
        if not isinstance(value, int):
            raise Fault(f"{path}: the file names no evaluation revision")
        return value

    def document(self, path):
        if not DOCUMENT.fullmatch(path):
            raise Fault(f"{path}: a quote comes from README.md or a page under docs/")
        if path not in self._text:
            self._text[path] = self._bytes(path).decode("utf-8")
        self.documents.add(path)
        return self._text[path]

    def quote(self, at, path, quote):
        if not isinstance(quote, str) or len(quote.split()) < 4:
            raise Fault(f"{at}: a quote is at least four words")
        if quote not in self.document(path):
            raise Fault(f"{at}: {path} does not hold the sentence {quote!r}")

    # A number, with the text that says where it was read. Every form returns (value, text, revisions).
    def number(self, at, expr):
        if not isinstance(expr, dict):
            raise Fault(f"{at}: a number is read from a file, never written in the spec")
        if "file" in expr:
            _keys(at, expr, {"file", "at", "length"}, {"file", "at"})
            path, value = self.data(expr["file"])
            pointer = "".join(f"/{k}" for k in expr["at"])
            if path == OUTSIDE or not FIELDS[path.rsplit("/", 1)[-1]].fullmatch(pointer):
                raise Fault(f"{at}: {path}#{pointer} is not a field the export may read")
            for key in expr["at"]:
                if not isinstance(value, dict) or key not in value:
                    raise Fault(f"{at}: {path}#{pointer} is not in the file")
                value = value[key]
            text = f"{path}#{pointer}"
            if expr.get("length"):
                if not isinstance(value, (list, dict)):
                    raise Fault(f"{at}: {text} has no length")
                value, text = len(value), f"length of {text}"
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise Fault(f"{at}: {text} is not a number")
            return value, text, {self.revision(expr["file"])}
        if "outside" in expr:
            return self._outside(at, expr)
        for op, sign in (("sum", " + "), ("minus", " - "), ("product", " * ")):
            if op in expr:
                _keys(at, expr, {op})
                parts = [self.number(at, x) for x in expr[op]]
                if len(parts) < 2 or (op == "minus" and len(parts) != 2):
                    raise Fault(f"{at}: {op} takes two numbers" + ("" if op == "minus" else " or more"))
                value = parts[0][0]
                for p in parts[1:]:
                    value = value + p[0] if op == "sum" else value - p[0] if op == "minus" else value * p[0]
                return value, sign.join(f"({p[1]})" if " " in p[1] else p[1] for p in parts), set().union(
                    *(p[2] for p in parts))
        raise Fault(f"{at}: unknown number form {sorted(expr)}")

    def _outside(self, at, expr):
        if self.aliases.get("outside") != OUTSIDE:
            raise Fault(f"{at}: files.outside must be {OUTSIDE}")
        path, doc = self.data("outside")
        population = next((p for p in doc["populations"] if p["id"] == expr["outside"]), None)
        if population is None:
            raise Fault(f"{at}: no population {expr['outside']!r}")
        revisions = {self.revision("outside")}
        if "total" in expr:
            _keys(at, expr, {"outside", "total"})
            kind = expr["total"]
            if kind not in OUTSIDE_TOTALS:
                raise Fault(f"{at}: unknown total {kind!r}")
            value = (len(population[kind]) if kind in ("datasets", "species")
                     else sum(s[kind] for s in population["species"]))
            return value, f"{path}: {OUTSIDE_TOTALS[kind]} of population {population['id']}", revisions
        _keys(at, expr, {"outside", "learner", "rows"}, {"learner", "rows"})
        learner, kind = expr["learner"], expr["rows"]
        if learner not in LEARNERS or kind not in OUTSIDE_ROWS:
            raise Fault(f"{at}: unknown learner {learner!r} or rows {kind!r}")
        count = 0
        for row in population["rows"]:
            if row["learner"] != learner:
                continue
            unread = [row[m].startswith("not read") for m in ("sensitivity", "specificity")]
            if unread[0] != unread[1]:
                raise Fault(f"{at}: a row of {path} is read on one measure only")
            if unread[0]:
                continue
            if (kind == "read" or (kind == "holds" and "below" not in (row["sensitivity"], row["specificity"]))
                    or (kind.startswith("sensitivity ") and row["sensitivity"] == kind.split()[1])):
                count += 1
        return count, f"{path}: rows of population {population['id']}, learner {learner}, {OUTSIDE_ROWS[kind]}", revisions


def _numbers(text):
    return {m.group().replace(",", "").rstrip(".") for m in NUMBER.finditer(text)}


def _count_of(at, count, of, note):
    for name, value in (("count", count), ("of", of)):
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 0):
            raise Fault(f"{at}: {name} {value!r} is not a count")
    if count is not None and of is not None and count > of:
        raise Fault(f"{at}: count {count} of {of} is not a count within its denominator")
    if of is None and not note:
        raise Fault(f"{at}: a line with no denominator says why in note")


def _line(at, spec, src):
    """One measured line of a node or of an assumption."""
    _keys(at, spec, LINE_KEYS, {"what"})
    out = {"count": None, "of": None, "what": spec["what"]}
    if "document" in spec:
        _keys(at, spec, {"what", "revision", "count", "of", "note", "document", "quote", "rebuilt_by"},
              {"document", "quote", "count"})
        src.quote(at, spec["document"], spec["quote"])
        terms = spec["count"] if isinstance(spec["count"], list) else [spec["count"]]
        stated = _numbers(spec["quote"])
        for n in terms + ([spec["of"]] if spec.get("of") is not None else []):
            if isinstance(n, bool) or not isinstance(n, int) or str(n) not in stated:
                raise Fault(f"{at}: {n!r} is not a number of the quoted sentence")
        out.update(count=sum(terms), of=spec.get("of"), source=spec["document"], basis="document")
        if "revision" in spec:
            out["evaluation_revision"] = spec["revision"]
        out["quote"] = spec["quote"]
        if "rebuilt_by" in spec:
            if not (src.root / spec["rebuilt_by"]).is_dir():
                raise Fault(f"{at}: no suite {spec['rebuilt_by']}")
            out["rebuilt_by"] = spec["rebuilt_by"]
    elif "values" in spec:
        _keys(at, spec, {"what", "of", "note", "values"}, {"values", "of", "note"})
        of, of_text, _ = src.number(at, spec["of"])
        values, paths = [], []
        for v in spec["values"]:
            _keys(at, v, {"revision", "value"}, {"revision", "value"})
            value, text, revisions = src.number(at, v["value"])
            if revisions != {v["revision"]}:
                raise Fault(f"{at}: {text} is of revision {sorted(revisions)}, and the line says {v['revision']}")
            values.append({"evaluation_revision": v["revision"], "value": value, "read": text})
            paths += re.findall(r"[\w./-]+\.json", text)
        paths += re.findall(r"[\w./-]+\.json", of_text)
        out.update(of=of, source="; ".join(dict.fromkeys(paths)), basis="file", values=values, read={"of": of_text})
    else:
        _keys(at, spec, {"what", "revision", "count", "of", "note"}, {"revision", "count"})
        count, count_text, revisions = src.number(at, spec["count"])
        read, paths = {"count": count_text}, re.findall(r"[\w./-]+\.json", count_text)
        if spec.get("of") is not None:
            out["of"], of_text, more = src.number(at, spec["of"])
            revisions, read["of"] = revisions | more, of_text
            paths += re.findall(r"[\w./-]+\.json", of_text)
        if revisions != {spec["revision"]}:
            raise Fault(f"{at}: the files read are of revision {sorted(revisions)}, and the line says {spec['revision']}")
        out.update(count=count, source="; ".join(dict.fromkeys(paths)), basis="file",
                   evaluation_revision=spec["revision"], read=read)
    if spec.get("note"):
        out["note"] = spec["note"]
    _count_of(at, out["count"], out["of"], out.get("note"))
    return out


def _entry(at, spec, src, out):
    """The parts a node and an assumption share: quotes, measured lines and what is unmeasured."""
    quotes = spec.get("quotes") or []
    if not quotes:
        raise Fault(f"{at}: an entry quotes at least one sentence of the repository")
    for q in quotes:
        _keys(at, q, {"path", "quote"}, {"path", "quote"})
        src.quote(at, q["path"], q["quote"])
    out["source"] = "; ".join(dict.fromkeys(q["path"] for q in quotes))
    out["quotes"] = [{"path": q["path"], "quote": q["quote"]} for q in quotes]
    out["measured"] = [_line(f"{at} line {i + 1}", x, src) for i, x in enumerate(spec.get("measured") or [])]
    if "unmeasured" in spec:
        if not isinstance(spec["unmeasured"], str) or not spec["unmeasured"].strip():
            raise Fault(f"{at}: unmeasured says what GRM has not counted")
        out["unmeasured"] = spec["unmeasured"]
    if not out["measured"] and "unmeasured" not in out and spec.get("counted", True):
        raise Fault(f"{at}: an entry with no measured line says so in unmeasured")
    return out


def _barred(root, src):
    """The names the export must never hold: datasets, species, antibiotics and drug classes of the files read."""
    names = set()
    for path, doc in src._json.items():
        if path == OUTSIDE:
            for p in doc["populations"]:
                names |= {s["slug"] for s in p["species"]}
                for d in p["datasets"]:
                    names |= {d["id"], d["species"], d["antibiotic"]}
        elif path.endswith("summary.json"):
            names |= set(doc.get("no_lineage") or [])
            names |= {s["species"] for s in doc.get("per_species") or []}
            names |= {c["class"] for c in doc.get("per_drug_class_min5") or []}
    barred = set()
    for name in names:
        low = name.lower()
        barred |= {low, low.replace("_", " ")}
        words = low.replace("_", " ").split()
        if len(words) == 2 and "/" not in low and "__" not in low and all(w.isalpha() for w in words):
            barred.add(f"{words[0][0]}. {words[1]}")  # the short form of a species name
    return {b for b in barred if len(b) > 3}


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for k, v in value.items():
            if k not in ("read", "path", "source", "files", "rebuilt_by"):
                yield from _strings(v)
    elif isinstance(value, list):
        for v in value:
            yield from _strings(v)


def screen(module, barred):
    """Fails when a string of the export names what the export must leave out."""
    for text in _strings(module):
        low = text.lower()
        hit = next((b for b in sorted(barred) if re.search(rf"(?<![a-z]){re.escape(b)}(?![a-z])", low)), None)
        if hit:
            raise Fault(f"the export may not name {hit!r}: {text!r}")
        if re.search(r"[ACGTacgt]{12,}", text):
            raise Fault(f"the export holds no sequence: {text!r}")


def build_module(spec, root=ROOT):
    """The export as a dict, from a spec dict. Raises Fault when a rule is broken."""
    _keys("spec", spec, MODULE_KEYS, MODULE_KEYS)
    if spec["id"] != "grm" or spec["namespace"] != "grm":
        raise Fault("the export is GRM's: id and namespace are grm")
    if spec["binding"] != "typed_export":
        raise Fault("the binding of a generated export is typed_export")
    src = Sources(root, spec["files"])
    module = {"format": FORMAT}
    for key in ("id", "namespace", "name", "job", "repo", "release", "binding", "binding_note"):
        if not isinstance(spec[key], str) or not spec[key].strip():
            raise Fault(f"{key} is required")
        module[key] = spec[key]
    nodes, ids = [], set()
    for n in spec["nodes"]:
        at = n.get("id", "a node")
        _keys(at, n, NODE_KEYS, {"id", "kind", "label", "definition", "quotes"})
        if not n["id"].startswith("grm:") or n["id"].startswith("grm:given/") or n["id"] in ids:
            raise Fault(f"{at}: an id starts with grm: and is defined once")
        if n["kind"] not in KINDS:
            raise Fault(f"{at}: unknown kind {n['kind']}")
        if "truth_state" in n and n["truth_state"] not in TRUTH_STATES:
            raise Fault(f"{at}: unknown truth state {n['truth_state']}")
        out = {"id": n["id"], "kind": n["kind"], "label": n["label"], "definition": n["definition"]}
        if "truth_state" in n:
            out["truth_state"] = n["truth_state"]
        ids.add(n["id"])
        nodes.append(_entry(at, dict(n, counted=False), src, out))
    edges = []
    for e in spec["edges"]:
        at = f"{e.get('from')} {e.get('type')} {e.get('to')}"
        _keys(at, e, EDGE_KEYS, EDGE_KEYS)
        if e["type"] not in EDGE_TYPES:
            raise Fault(f"{at}: unknown edge type")
        if e["from"] not in ids or e["to"] not in ids:
            raise Fault(f"{at}: both ends of an edge are nodes of this export")
        if e in edges:
            raise Fault(f"{at}: written twice")
        edges.append({"from": e["from"], "type": e["type"], "to": e["to"]})
    givens, given_ids = [], set()
    for g in spec["assumptions"]:
        at = g.get("id", "an assumption")
        _keys(at, g, GIVEN_KEYS, {"id", "statement", "under", "hidden_in", "when_false", "quotes"})
        if not g["id"].startswith("grm:given/") or g["id"] in given_ids:
            raise Fault(f"{at}: an assumption's id starts with grm:given/ and is defined once")
        if not g["under"] or any(u not in ids for u in g["under"]):
            raise Fault(f"{at}: an assumption sits under at least one node of this export")
        given_ids.add(g["id"])
        out = {"id": g["id"], "statement": g["statement"], "under": list(g["under"]), "hidden_in": g["hidden_in"],
               "when_false": g["when_false"]}
        givens.append(_entry(at, g, src, out))
    for key in ("import_rule", "left_out"):
        if not spec[key] or not all(isinstance(x, str) and x.strip() for x in spec[key]):
            raise Fault(f"{key} is a list of sentences")
    module["sources"] = {
        "program": "scripts/givens-export/export.py",
        "spec": "scripts/givens-export/spec.json",
        "files": [{"path": p, "sha256": h} for p, h in sorted(src.files.items())],
        "documents": sorted(src.documents),
    }
    unused = sorted(set(spec["files"]) - src.used)
    if unused:
        raise Fault(f"files names {', '.join(unused)}, and no line reads it")
    module.update(nodes=nodes, edges=edges, assumptions=givens, import_rule=list(spec["import_rule"]),
                  left_out=list(spec["left_out"]))
    screen(module, _barred(root, src))
    return module


def build(root=ROOT, spec_path=SPEC):
    """The export's text."""
    spec = json.loads(pathlib.Path(spec_path).read_text(encoding="utf-8"))
    return json.dumps(build_module(spec, root), indent=2, ensure_ascii=False) + "\n"


def main(argv):
    args = argv[1:]
    if args and (len(args) != 2 or args[0] != "--check"):
        raise SystemExit(__doc__)
    try:
        text = build()
    except Fault as fault:
        print(f"givens export: {fault}", file=sys.stderr)
        return 1
    if args:
        if pathlib.Path(args[1]).read_text(encoding="utf-8") == text:
            return 0
        print("the file differs from what the committed sources give", file=sys.stderr)
        return 1
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
