"""revlasso: recall beside log loss for the 250 test samples of lassoscan's
one-fit measurement (brief point 5.4), from grm's refit scores in the dump.

usage: python3 -I testrecall.py DUMP.json OUT_JSON
"""
import json
import math
import sys
from pathlib import Path


def main(dump, outp):
    d = json.loads(Path(dump).read_text())
    y = set(d["y"])
    tr = set(d["train"])
    test = [i for i in range(d["n"]) if i not in tr]
    s = d["scores"]
    tp = sum(1 for i in test if i in y and s[i] >= 0.5)
    fn = sum(1 for i in test if i in y and s[i] < 0.5)
    fp = sum(1 for i in test if i not in y and s[i] >= 0.5)
    tn = sum(1 for i in test if i not in y and s[i] < 0.5)
    terms = []
    for i in test:
        q = min(max(s[i], 1e-15), 1 - 1e-15)
        terms.append(-math.log(q) if i in y else -math.log(1 - q))
    out = {"test": len(test), "resistant": tp + fn, "tp": tp, "fn": fn, "fp": fp, "tn": tn,
           "recall": tp / (tp + fn), "specificity": tn / (tn + fp),
           "log_loss_clip_1e-15": math.fsum(terms) / len(test)}
    print(json.dumps(out))
    Path(outp).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:3])
