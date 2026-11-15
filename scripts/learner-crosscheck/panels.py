"""panels.py FIXTURES

The saved test panels of the learner cross-check, decoded without the Rust
code: tests/data/lasso-stopping-fixtures.txt. It needs no library.

A panel is a dict: `name`, `n` samples, `y` (resistant), `fit` (the training
mask) and `rows`, one presence set per feature. A set of samples is an int
used as a bitset: bit i is sample i. The file writes a set as hex nibbles, low
bit first, so nibble k holds samples 4k to 4k + 3.

Run as a script it prints one line per panel: its name and counts.
"""

import sys

LINE_TYPES = ("fixture", "y", "fit", "r", "point")


def decode_bits(n, text):
    """The sample set of one `y`, `fit` or `r` field."""
    if len(text) != -(-n // 4):
        raise ValueError(f"a set of {n} samples has {len(text)} nibbles")
    bits = int(text[::-1], 16) if text else 0
    if bits >> n:
        raise ValueError(f"a bit is set past sample {n}")
    return bits


def encode_bits(n, bits):
    """The inverse of decode_bits, as the Rust harness writes sets."""
    return "".join("0123456789abcdef"[(bits >> (4 * k)) & 15] for k in range(-(-n // 4)))


def members(bits):
    """The samples of a set, ascending."""
    out = []
    while bits:
        low = bits & -bits
        out.append(low.bit_length() - 1)
        bits ^= low
    return out


def count(bits):
    return bin(bits).count("1")


def read(path):
    """Every panel of the fixture file, in file order."""
    panels = []
    with open(path, encoding="ascii") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split(" ")
            kind = fields[0]
            where = f"{path} line {number}"
            if kind not in LINE_TYPES:
                raise ValueError(f"{where}: unknown line type {kind!r}")
            if kind == "fixture":
                if len(fields) != 6:
                    raise ValueError(f"{where}: a fixture line has 6 fields")
                panels.append({"name": fields[1], "n": int(fields[2]), "declared_features": int(fields[3]), "rows": []})
                continue
            if not panels:
                raise ValueError(f"{where}: no fixture yet")
            panel = panels[-1]
            if kind == "point":
                continue  # the saved lasso paths, which this check does not use
            if len(fields) != 2:
                raise ValueError(f"{where}: a set line has 2 fields")
            try:
                bits = decode_bits(panel["n"], fields[1])
            except ValueError as error:
                raise ValueError(f"{where}: {error}") from None
            if kind == "r":
                panel["rows"].append(bits)
            elif kind in panel:
                raise ValueError(f"{where}: a second {kind} line")
            else:
                panel[kind] = bits
    for panel in panels:
        if "y" not in panel or "fit" not in panel:
            raise ValueError(f"panel {panel['name']} lacks y or fit")
        if len(panel["rows"]) != panel.pop("declared_features"):
            raise ValueError(f"panel {panel['name']} has {len(panel['rows'])} feature rows, not the count it declares")
    names = [panel["name"] for panel in panels]
    if len(set(names)) != len(names):
        raise ValueError("a panel name is used twice")
    return panels


def counts(panel):
    """The counts the Rust harness writes for a panel, for the decode check."""
    return {
        "samples": panel["n"],
        "features": len(panel["rows"]),
        "resistant": count(panel["y"]),
        "fit": count(panel["fit"]),
        "fit_resistant": count(panel["fit"] & panel["y"]),
        "present": [count(row) for row in panel["rows"]],
    }


def main(argv):
    if len(argv) != 2:
        sys.exit(__doc__.splitlines()[0])
    for panel in read(argv[1]):
        c = counts(panel)
        print(panel["name"], c["samples"], c["features"], c["resistant"], c["fit"], c["fit_resistant"], sep="\t")


if __name__ == "__main__":
    main(sys.argv)
