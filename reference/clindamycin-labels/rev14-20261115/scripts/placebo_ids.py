#!/usr/bin/env python3
"""Placebo removals (PLAN.md section 5): for seed s in 1..10, remove from the samples outside the chosen study as
many reported-resistant and as many reported-susceptible samples as the study holds, drawn with
random.Random(20261114 * 100 + s) over the sorted ids (resistant draw first, then susceptible). No placebo when the
rest holds too few of a class, or when the removal would leave a class under min_class.

Run: python3 -I scripts/placebo_ids.py <membership tsv> <min_class> <out dir>
Writes <out dir>/placebo-<s>.ids and prints one line per seed.
"""
import os
import random
import sys


def main():
    member, min_class, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    os.makedirs(out, exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in open(member)][1:]
    # columns: genome_id reported known_gene_rule study_key source_field lineage n_label_records in_chosen_study
    in_r = sum(1 for r in rows if r[7] == "1" and r[1] == "R")
    in_s = sum(1 for r in rows if r[7] == "1" and r[1] == "S")
    pool_r = sorted(r[0] for r in rows if r[7] == "0" and r[1] == "R")
    pool_s = sorted(r[0] for r in rows if r[7] == "0" and r[1] == "S")
    tot_r = in_r + len(pool_r)
    tot_s = in_s + len(pool_s)
    if len(pool_r) < in_r or len(pool_s) < in_s:
        print(f"no placebo: the rest holds R {len(pool_r)} S {len(pool_s)}, the study R {in_r} S {in_s}")
        return
    if tot_r - in_r < min_class or tot_s - in_s < min_class:
        print(f"no placebo: removal leaves R {tot_r - in_r} S {tot_s - in_s}, under {min_class}")
        return
    for s in range(1, 11):
        rng = random.Random(20261114 * 100 + s)
        pick = rng.sample(pool_r, in_r) + rng.sample(pool_s, in_s)
        with open(os.path.join(out, f"placebo-{s}.ids"), "w") as fh:
            fh.write("\n".join(sorted(pick)) + "\n")
        print(f"placebo {s}: removed R {in_r} S {in_s} (pool R {len(pool_r)} S {len(pool_s)})")


if __name__ == "__main__":
    main()
