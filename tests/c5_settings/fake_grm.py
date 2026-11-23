"""A SYNTHETIC stand-in for `grm learn --fold-file F --fold-out O --exp-settings S --force`.

  python3 -I fake_grm.py --data <root> learn --species <slug> --fold-file F --fold-out O --exp-settings S --force

It reads <root>/fake.json (the made-up genomes' labels and how good each
learner's trial is), refuses a fold file that does not list every genome
exactly once, appends one line to <root>/calls.jsonl (the scheme, the
genomes the fold file tests and trains, and whether the driver's pick.json
stood beside the fold-out folder when it ran), and writes one trial result
per trial in the form of src/learn/settings.rs. A genome's score is a
made-up function of its label, the learner and the trial; no fit runs and
no real data is read.
"""

import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from hyperopt import protocol as hyper_protocol  # noqa: E402
from hyperopt import search as hyper_search  # noqa: E402

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
PROTOCOL = {"families": hyper_protocol.FAMILIES, "seeds": {"sampler": 20261030}}
CLIP = 1e-4
NOISE = 1.0


def unit(*parts):
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(h[:8], "big") / 2**64


def ident(learner, setting):
    """The trial a setting is: 0 for the release's, else the lowest of the protocol's draws 1 to 30 it equals.
    A genome's score follows the setting, not the id the settings file gives it."""
    if setting == "release":
        return 0
    for t in range(1, 31):
        if hyper_search.sample_config(PROTOCOL, learner, t) == setting:
            return t
    return -1


def quality(world, learner, trial, outer):
    key = f"{learner}:{trial}"
    if outer and key in world.get("outer_quality", {}):
        return world["outer_quality"][key]
    if key in world.get("quality", {}):
        return world["quality"][key]
    return 0.75 if trial == 0 else 0.55 + 0.2 * unit("q", learner, trial)


def score(world, genome, learner, trial, outer):
    """A made-up score: near q for a resistant genome, near 1 - q for a susceptible one."""
    q = quality(world, learner, trial, outer)
    noise_trial = world.get("alias", {}).get(f"{learner}:{trial}", trial)
    y = world["labels"][genome]
    s = (q if y else 1 - q) + (unit("n", genome, learner, noise_trial) - 0.5) * NOISE
    if not outer:
        s += world.get("inner_shift", 0.0)
    if learner == "scm":
        s = 1.0 if s >= 0.5 else 0.0
    return min(max(s, 0.0), 1.0)


def calls(world, genomes, learner, trial, cutoff, outer):
    c = {"n": 0, "tp": 0, "fp": 0, "tn": 0, "fn": 0, "loss_sum": 0.0}
    for g in genomes:
        s = score(world, g, learner, trial, outer)
        y = world["labels"][g]
        call = s >= cutoff
        c["n"] += 1
        c["tp" if y and call else "fn" if y else "fp" if call else "tn"] += 1
        p = min(max(s, CLIP), 1 - CLIP)
        c["loss_sum"] -= math.log(p if y else 1 - p)
    return c


def balanced(c):
    pos, neg = c["tp"] + c["fn"], c["tn"] + c["fp"]
    return None if not pos or not neg else (c["tp"] / pos + c["tn"] / neg) / 2


def main(argv):
    if argv[0] != "--data" or argv[2] != "learn":
        return 2
    root = Path(argv[1])
    opts = {}
    rest = argv[3:]
    i = 0
    while i < len(rest):
        if rest[i] == "--force":
            i += 1
            continue
        opts[rest[i]] = rest[i + 1]
        i += 2
    world = json.loads((root / "fake.json").read_text())
    fold_file = json.loads(Path(opts["--fold-file"]).read_text())
    settings = json.loads(Path(opts["--exp-settings"]).read_text())
    out = Path(opts["--fold-out"])
    listed = [g for fold in fold_file["test"] for g in fold] + fold_file["train"] + fold_file["unused"]
    if sorted(listed) != sorted(world["labels"]) or len(set(listed)) != len(listed):
        print("fake grm: the fold file does not list every genome once", file=sys.stderr)
        return 2
    outer = fold_file["scheme"].startswith("c5o-")
    with open(root / "calls.jsonl", "a") as f:
        f.write(json.dumps({"scheme": fold_file["scheme"],
                            "tested": sorted(g for fold in fold_file["test"] for g in fold),
                            "trained": sorted(fold_file["train"]),
                            "pick_written": (out.parent / "pick.json").is_file(),
                            "trials": [t["id"] for t in settings["trials"]]}) + "\n")
    out.mkdir(parents=True, exist_ok=True)
    for t in settings["trials"]:
        rows = []
        for learner in t["learners"]:
            tid = ident(learner, t["learners"][learner])
            cut = t.get("cutoffs", {}).get(learner)
            cutoffs = sorted({0.5, *settings["cutoffs"], *([cut] if cut is not None else [])})
            folds = fold_file["test"]
            rows.append({
                "learner": learner, "setting": t["learners"][learner], "status": "ok", "reason": None,
                "cutoff": cut,
                "fold_balanced_accuracy": [balanced(calls(world, fold, learner, tid, 0.5, outer)) for fold in folds],
                "fold_complexity": [1.0 for _ in folds], "no_call_resistant": 0, "no_call_susceptible": 0,
                "calls": [{"cutoff": c, "folds": [calls(world, fold, learner, tid, c, outer) for fold in folds]}
                          for c in cutoffs],
                "seconds": 0.001})
        doc = {"format": "grm-c5-trial-1", "evaluation_revision": 14, "dataset": fold_file["dataset"],
               "scheme": fold_file["scheme"], "fold_file_sha256": "0" * 64, "settings_sha256": "0" * 64,
               "seed": settings["seed"], "tuning": settings["tuning"], "trial": t["id"], "min_support": 3,
               "folds": [{"fold": k, "test": len(fold)} for k, fold in enumerate(fold_file["test"])],
               "learners": rows}
        (out / f"{fold_file['dataset']}.{fold_file['scheme']}.t{t['id']}.json").write_text(json.dumps(doc))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
