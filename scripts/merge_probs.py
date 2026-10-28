#!/usr/bin/env python3
"""Merge probability tables (header 'dt mu3 mu8 p01..p48') into one master table.

    scripts/merge_probs.py out.txt master_in.txt other.txt ...

Rows are keyed on (dt, mu3, mu8). When the same key appears in several inputs the row from the
LAST input that has it is kept (give the more precise table last) and the largest difference is
reported. Number text is copied verbatim.
"""
import sys

N = 48
HEADER = ["dt", "mu3", "mu8"] + [f"p{i:02d}" for i in range(1, N + 1)]


def read(path):
    rows = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            tok = line.split("#")[0].split()
            if not tok or tok == HEADER:
                continue
            if len(tok) != len(HEADER):
                sys.exit(f"{path}:{n}: expected {len(HEADER)} values, found {len(tok)}")
            key = tuple(float(t) for t in tok[:3])
            if key in rows:
                sys.exit(f"{path}:{n}: duplicate key {key}")
            rows[key] = tok
    return rows


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    out, inputs = argv[1], argv[2:]
    merged, source = {}, {}
    for path in inputs:
        for key, tok in read(path).items():
            if key in merged:
                diff = max(abs(float(a) - float(b)) for a, b in zip(merged[key], tok))
                print(f"key dt={key[0]} mu3={key[1]} mu8={key[2]}: in {source[key]} and {path}; "
                      f"max difference {diff:.3g}; keeping {path}")
            merged[key], source[key] = tok, path
    with open(out, "w") as f:
        f.write("# Master update-probability table. dt = 4*beta/Lt; one row per (dt, mu3, mu8).\n")
        f.write("# The simulator looks rows up by the run's parameters; keys must be unique.\n")
        f.write(" ".join(HEADER) + "\n")
        for key in sorted(merged):
            f.write(" ".join(merged[key]) + "\n")
    print(f"wrote {len(merged)} rows to {out}")


main(sys.argv)
