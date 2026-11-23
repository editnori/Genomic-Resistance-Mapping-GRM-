"""Time and memory per step and fit of one species' C2 run.

  python3 -I timing.py <logs dir>

Reads every <name>.time that GNU time wrote in run_species.sh's form
("<name> <wall> s wall, <peak> KB peak, <user> s user, <sys> s sys") and
prints step, wall seconds, CPU seconds (user plus system) and peak memory
in KB, one row per file, sorted by name.
"""

import re
import sys
from pathlib import Path

LINE = re.compile(r"^(\S+) ([0-9.]+) s wall, ([0-9]+) KB peak, ([0-9.]+) s user, ([0-9.]+) s sys$")


def rows(logs):
    out = []
    for p in sorted(Path(logs).glob("*.time")):
        text = [t for t in p.read_text(encoding="utf-8").splitlines() if LINE.match(t)]
        if not text:
            out.append((p.stem, "-", "-", "-"))
            continue
        _, wall, peak, user, sys_ = LINE.match(text[-1]).groups()
        out.append((p.stem, wall, f"{float(user) + float(sys_):.2f}", peak))
    return out


def main(argv):
    if len(argv) != 2:
        raise SystemExit(__doc__)
    print("step\twall_s\tcpu_s\tpeak_kb")
    for r in rows(argv[1]):
        print("\t".join(r))


if __name__ == "__main__":
    main(sys.argv)
