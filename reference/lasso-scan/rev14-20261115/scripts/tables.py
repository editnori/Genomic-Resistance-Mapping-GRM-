"""Markdown tables for REPORT.md from the scan summary, the comparison and item 5.

usage: python3 -I tables.py SUMMARY.json COMPARE.json [ITEM5.json] [ITEM5_POLISH.json]
"""
import json
import sys
from pathlib import Path


def short(s):
    return s.replace("mycobacterium_tuberculosis", "MTB").replace("klebsiella_pneumoniae", "KP").replace(
        "salmonella_enterica", "SE").replace("escherichia_coli", "EC").replace("staphylococcus_aureus", "SA")


def main(argv):
    summ = json.loads(Path(argv[0]).read_text())
    comp = {r["id"]: r for r in json.loads(Path(argv[1]).read_text())} if len(argv) > 1 and Path(argv[1]).exists() else {}
    print("| Dataset | Fingerprint = release | Fits | Converged (scan rule) | Path points | Largest IRLS steps | Points with a short solve | Lasso s (10 outer fits) | Final fit s | Release lasso s | Rev 7 lasso s |")
    print("|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|")
    for d in summ["datasets"]:
        pts = sum(f["path_points"] for f in d["fits"])
        print(f"| `{short(d['id'])}` | {'yes' if d['fingerprint_match'] else 'NO'} | {len(d['fits'])} | {'yes' if d['converged'] else 'NO: ' + '; '.join(d['not_converged_reasons'])} "
              f"| {pts} | {d['irls_max']} | {d['capped_points']} of {d['points']} | {d['lasso_seconds_outer']:.2f} | {d['final_fit_seconds']:.2f} "
              f"| {d['release_lasso_seconds']:.2f} | {d['rev7_lasso_seconds']:.1f} |")
    if comp:
        print()
        print("| Dataset | Converged (scan rule) | Fits | Path + refit points | IRLS steps, total (largest at a point) | Calls changed, random + lineage | Recall random | Log loss random | Recall lineage | Log loss lineage | Largest change in a released measure |")
        print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for d in summ["datasets"]:
            c = comp.get(d["id"])
            if not c:
                continue
            sp = {s["scheme"]: s for s in c["splits"]}
            big = max((abs(v) for s in c["splits"] for v in s["delta"].values() if v is not None), default=None)
            cell = lambda sch, k: f"{sp[sch]['scan'][k]:.4f}" if sch in sp else "n/a"
            print(f"| `{short(d['id'])}` | {'yes' if d['converged'] else 'NO'} | {len(d['fits'])} | {d['points']} "
                  f"| {sum(f['irls_total'] for f in d['fits'])} ({d['irls_max']}) | {sum(s['calls_changed'] for s in c['splits'])} "
                  f"| {cell('random', 'recall')} | {cell('random', 'log_loss')} | {cell('lineage', 'recall')} | {cell('lineage', 'log_loss')} "
                  f"| {big if big is None else format(big, 'g')}{'' if not c['differences'] else ' (differences: see compare.json)'} |")
        print()
        print("| Dataset | Split | Scores equal to release | Largest score change | Calls changed | Recall release | Recall scan | Log loss release | Log loss scan | Audits, weights, final model |")
        print("|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for d in summ["datasets"]:
            c = comp.get(d["id"])
            if not c:
                continue
            same = "equal" if not c["differences"] else "; ".join(c["differences"])[:200]
            for s in c["splits"]:
                print(f"| `{short(d['id'])}` | {s['scheme']} | {s['scores_identical']} of {s['n']} | {s['max_abs_score_diff']:.1e} | {s['calls_changed']} "
                      f"| {s['release']['recall']:.4f} | {s['scan']['recall']:.4f} | {s['release']['log_loss']:.4f} | {s['scan']['log_loss']:.4f} | {same} |")
    if len(argv) > 2 and Path(argv[2]).exists():
        it = json.loads(Path(argv[2]).read_text())
        pol = {}
        if len(argv) > 3 and Path(argv[3]).exists():
            pol = {r["point"]: r for r in json.loads(Path(argv[3]).read_text())["points"]}
        e = lambda v, f=".1e": "none" if v is None else format(v, f)
        print()
        print("| Point | lambda | grm active | grm IRLS steps | grm KKT | P_grm - P_ind | Independent KKT | Independent active | Clarabel status, s "
              "| Bound for grm, planned dual | Bound for grm, exploratory | Active at certified optimum |")
        print("|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|")
        for r in it.get("points", []):
            q = pol.get(r["point"], {})
            print(f"| {r['point']} | {r['lambda']:.3e} | {r['grm_active']} | {r['grm_irls_steps']} | {e(r['kkt_grm'])} | {r['P_grm_minus_P_ind']:+.1e} "
                  f"| {e(r['kkt_ind_thresholded'])} | {r['ind_active_thresholded']} | {r['clarabel']['status']}, {r['clarabel']['seconds']:.1f} "
                  f"| {e(r['certified_bound_grm'])} | {e(q.get('certified_bound_grm'))} | {q.get('polished_active', 'n/a')} |")

if __name__ == "__main__":
    main(sys.argv[1:])
