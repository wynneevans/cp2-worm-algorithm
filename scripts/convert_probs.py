#!/usr/bin/env python3
"""Convert an old headerless probabilities file (51 numbers per line:
dt mu3 mu8 p1..p48) into the labelled format the simulator now reads.

    scripts/convert_probs.py old.txt new.txt

Number text is copied verbatim, so values are bit-for-bit identical.
"""
import sys

N_WEIGHTS = 48


def convert(src, dst):
    header = ["dt", "mu3", "mu8"] + [f"p{i:02d}" for i in range(1, N_WEIGHTS + 1)]
    rows = []
    with open(src) as f:
        for n, line in enumerate(f, 1):
            tok = line.split()
            if not tok:
                continue
            if len(tok) != len(header):
                sys.exit(f"{src}:{n}: expected {len(header)} numbers, found {len(tok)}")
            for t in tok:
                float(t)  # raises on anything that is not a number
            rows.append(tok)
    if not rows:
        sys.exit(f"{src}: no rows found")
    with open(dst, "w") as out:
        out.write("# Update-probability weights. dt = 4*beta/Lt; one row per (dt, mu3, mu8).\n")
        out.write(" ".join(header) + "\n")
        for r in rows:
            out.write(" ".join(r) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    convert(sys.argv[1], sys.argv[2])
