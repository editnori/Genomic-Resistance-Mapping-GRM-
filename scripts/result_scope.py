"""Select current, run-scoped dataset results before either checker reads them."""
import fnmatch
import json
import re
import sys
from pathlib import Path

CLASS_MINIMUM = "too few labeled genomes after download, skipped"
SKIP = re.compile(r"\s+(\S+__\S+): " + CLASS_MINIMUM + r"\s*$")
DONE = re.compile(r"\s+(\S+__\S+) done in \S+ \(.*\)\s*$")
RECOMPUTE = re.compile(r"\s+(\S+__\S+): .+; recomputing report\s*$")


def class_minimum_skips(log):
    """Only an exact skip line, with no later done/recompute, excuses absence."""
    skipped = set()
    if log is not None:
        with open(log, errors="replace") as stream:
            for line in stream:
                if match := SKIP.fullmatch(line):
                    skipped.add(match.group(1))
                elif match := DONE.fullmatch(line) or RECOMPUTE.fullmatch(line):
                    skipped.discard(match.group(1))
    return skipped


def select_results(root, pattern="*", log=None):
    """Return result paths and explicit exclusions, skips, and missing results.

    Existing results keep their species in scope. Otherwise both run inputs,
    genomes.txt and counted.json, are required. A specific species glob also
    selects matching curated species even before either run input exists.
    """
    root = Path(root)
    curated = json.loads((root / "curated/datasets.json").read_text())
    expected = {
        root / "species" / row["species_slug"] / "results" / (row["id"] + ".json"): row
        for row in curated if fnmatch.fnmatchcase(row["species_slug"], pattern)
    }
    files = sorted((root / "species").glob(pattern + "/results/*.json"))
    species = {row["species_slug"] for row in expected.values()}
    judged = {path.parent.parent.name for path in files}
    judged.update(slug for slug in species
                  if (root / "species" / slug / "genomes.txt").is_file()
                  and (root / "species" / slug / "counted.json").is_file())
    if pattern != "*":
        judged.update(species)

    scope = {
        "not_in_current_curation": [
            {"file": str(path), "reason": "not in current curation"}
            for path in files if path not in expected
        ],
        "not_run": [{"species": slug, "reason": "not run"} for slug in sorted(species - judged)],
        "class_minimum_skipped": [],
        "missing_curated_results": [],
    }
    skipped = class_minimum_skips(log)
    present = set(files)
    for path, row in sorted(expected.items()):
        if row["species_slug"] not in judged or path in present:
            continue
        if row["id"] in skipped:
            scope["class_minimum_skipped"].append({"file": str(path), "reason": CLASS_MINIMUM})
        else:
            scope["missing_curated_results"].append(
                {"file": str(path), "error": "missing results for curated dataset"})
    for key in ("not_in_current_curation", "not_run", "class_minimum_skipped"):
        for entry in scope[key]:
            print("%s: %s" % (entry.get("file", entry.get("species")), entry["reason"]), file=sys.stderr)
    return [path for path in files if path in expected], scope
