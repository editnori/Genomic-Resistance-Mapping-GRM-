"""Solves some of item 5's path points with the independent solver, exactly as
item5.py does (same problem, same clarabel_solve), and saves each solution in
OUTDIR/ind/pointKK.npz, where item5.py picks it up. Lets the 30 independent
points run in parallel.

usage: venv/bin/python -I item5-solve.py DUMP.json OUTDIR K [K ...]
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from item5 import Prob, clarabel_solve, dense  # noqa: E402


def main(argv):
    dump, outdir, ks = Path(argv[0]), Path(argv[1]), [int(a) for a in argv[2:]]
    d = json.loads(dump.read_text())
    Xall, yall = dense(d)
    P = Prob(Xall, yall, sorted(d["fit"]))
    (outdir / "ind").mkdir(parents=True, exist_ok=True)
    for k in ks:
        saved = outdir / "ind" / f"point{k:02d}.npz"
        if saved.exists():
            continue
        lam = d["path"][k]["model"]["lambda"]
        ib0, ib, info = clarabel_solve(P, lam)
        tmp = outdir / "ind" / f"point{k:02d}.tmp.npz"
        np.savez(tmp, b0=ib0, b=ib, info=json.dumps(info))
        tmp.rename(saved)
        print(f"point {k}: {info['status']} {info['seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
