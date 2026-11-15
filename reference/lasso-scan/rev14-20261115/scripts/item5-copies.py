"""Item 5, exploratory check (decided at 15:31Z after reading item5.py's
validation log losses): do the two solves' different validation losses come
from exact copies? For each path point, the independent solution's weight on
each group of columns identical on the fitting samples is moved to the group's
lowest-index column (what grm's collapse does). The fitting objective does not
change; the validation log loss is recomputed and set beside grm's.

usage: venv/bin/python -I item5-copies.py DUMP.json ITEM5_DIR OUT.json
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from item5 import Prob, balanced_log_loss, dense, grm_coef, sigmoid  # noqa: E402


def main(argv):
    dump, idir, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    d = json.loads(dump.read_text())
    Xall, yall = dense(d)
    fit, val = sorted(d["fit"]), np.array(sorted(d["val"]))
    P = Prob(Xall, yall, fit)
    Xf = Xall[fit]
    # Lowest-index identical column on the fitting samples, by exact column bytes.
    first, rep = {}, np.arange(Xf.shape[1])
    for k in range(Xf.shape[1]):
        key = Xf[:, k].tobytes()
        rep[k] = first.setdefault(key, k)
    it = json.loads((idir / "item5.json").read_text())
    rows = []
    for k, f in enumerate(d["path"]):
        lam = f["model"]["lambda"]
        z = np.load(idir / "ind" / f"point{k:02d}.npz")
        ib0, ib = float(z["b0"]), z["b"]
        cb = np.zeros_like(ib)
        np.add.at(cb, rep, ib)
        gb0, gb = grm_coef(f["model"], d["features"])
        XA = Xall.astype(np.float64)
        v_ind = balanced_log_loss(sigmoid(ib0 + XA @ ib), yall, val)
        v_col = balanced_log_loss(sigmoid(ib0 + XA @ cb), yall, val)
        v_grm = balanced_log_loss(sigmoid(gb0 + XA @ gb), yall, val)
        rows.append({"point": k, "lambda": lam, "status": it["points"][k]["clarabel"]["status"],
                     "P_ind": P.objective(lam, ib0, ib), "P_ind_collapsed": P.objective(lam, ib0, cb),
                     "val_ind": v_ind, "val_ind_collapsed": v_col, "val_grm": v_grm,
                     "val_ind_minus_grm": v_ind - v_grm, "val_collapsed_minus_grm": v_col - v_grm})
        print(f"point {k}: val ind-grm {v_ind - v_grm:+.2e} collapsed-grm {v_col - v_grm:+.2e} "
              f"P change {P.objective(lam, ib0, cb) - P.objective(lam, ib0, ib):+.1e}", flush=True)
    out.write_text(json.dumps({"label": "exploratory, decided after reading item5.py's validation log losses",
                               "groups_with_copies": int(len(set(rep[rep != np.arange(len(rep))]))),
                               "redundant_copies": int((rep != np.arange(len(rep))).sum()), "points": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
