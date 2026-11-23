"""The record of an arm's study fit that grm learn refused for no candidate (docs/input-arms.md, ruling 16).

  python3 -I refused.py <fit log> <dataset> <arm> <out file>

grm learn --fold-file refuses the whole run of a dataset when an outer
fold's training genomes give no candidate feature (eval::evaluate_assigned);
the drawn splits never refuse for it. run_species.sh calls this on the log
of an arm's failed study fit. When the log holds grm's sentence for that
refusal at the study-whole scheme, and the arm is neither the default nor
the twin, it writes <out file>: the format, the dataset, arm, scheme and
seed of the run, the reason, the fold's number, its training genomes and the
support floor, and grm's sentence, which holds those counts and nothing
else. score.py scores the record as the constant model on the default's
folds. Otherwise it writes nothing and exits 1, and the run stops.
"""

import json
import re
import sys
from pathlib import Path

FORMAT = "grm-c3-refused-1"
REASON = "no candidate"
SCHEME = "study-whole"
SEED = 42
DEFAULT, TWIN = "k31", "twin"
# grm's sentence, src/learn/eval.rs, evaluate_assigned.
SENTENCE = re.compile(
    r"(?P<scheme>[\w-]+) outer fold (?P<fold>\d+)'s (?P<train>\d+) training samples give no candidate feature: "
    r"none is present in at least (?P<floor>\d+) of them and absent from at least (?P=floor)")


def record(log, dataset, arm):
    """The refusal record for an arm's study fit whose log is `log`, or None."""
    if arm in (DEFAULT, TWIN):
        return None
    m = SENTENCE.search(log)
    if m is None or m["scheme"] != SCHEME:
        return None
    return {"format": FORMAT, "dataset": dataset, "arm": arm, "scheme": SCHEME, "seed": SEED, "reason": REASON,
            "fold": int(m["fold"]), "training_genomes": int(m["train"]), "support_floor": int(m["floor"]),
            "detail": m.group(0)}


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    log, dataset, arm, out = argv
    rec = record(Path(log).read_text(errors="replace"), dataset, arm)
    if rec is None:
        return 1
    out = Path(out)
    part = out.with_name(out.name + ".part")
    part.write_text(json.dumps(rec, indent=1) + "\n")
    part.replace(out)
    print(f"refused for no candidate: outer fold {rec['fold']}, {rec['training_genomes']} training genomes")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
