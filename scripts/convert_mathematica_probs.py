#!/usr/bin/env python3
"""Convert Mathematica notebook output of update probabilities into probability-table rows.

    scripts/convert_mathematica_probs.py OUTPUT.txt --dt 0.05 --vary mu8 [--fixed 0] > rows.txt

OUTPUT.txt is the copied cell output: a line with the chemical potential (e.g. "0.01"),
then 12 lines of SequenceForm[...] holding the 48 weights p01..p48 (six lines of 6, six of 2),
repeated for each value. --vary names the potential on those lines (mu3 or mu8); the other
one is --fixed. Prints "dt mu3 mu8 p01 ... p48" rows for configs/probabilities.txt, with
values within 1e-12 of 0 or 1 snapped to 0 or 1 (round-off in the notebook) and the rest
written with 15 significant digits, as in the existing table. Merge them with
scripts/merge_probs.py.
"""

import argparse
import re
import sys

N_WEIGHTS = 48
SEQ = re.compile(r"SequenceForm\[(.*?)\]")
NUM = re.compile(r"[-+]?\d+\.?\d*(?:`[\d.]*)?(?:\*\^[-+]?\d+)?")


def mathematica_number(tok: str) -> float:
    """'2.138220058923884`*^-15' -> 2.138220058923884e-15; '1`' -> 1.0."""
    mantissa, _, exponent = tok.partition("*^")
    return float(mantissa.split("`")[0]) * 10.0 ** int(exponent or 0)


def parse(text: str) -> list[tuple[str, list[float]]]:
    """(potential as written, 48 weights) per block."""
    blocks: list[tuple[str, list[float]]] = []
    for n, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        if re.fullmatch(r"\d+(\.\d+)?", line):
            blocks.append((line, []))
            continue
        m = SEQ.search(line)
        if not m or not blocks:
            sys.exit(f"line {n}: expected a potential or a SequenceForm[...] line")
        body = re.sub(r'"[^"]*"', " ", m.group(1))  # drop the " " separators
        blocks[-1][1].extend(mathematica_number(t) for t in NUM.findall(body))
    for mu, w in blocks:
        if len(w) != N_WEIGHTS:
            sys.exit(f"value {mu}: expected {N_WEIGHTS} weights, found {len(w)}")
        if any(not -1e-12 <= x <= 1 + 1e-12 for x in w):
            sys.exit(f"value {mu}: a weight is outside [0, 1]")
    return blocks


def fmt(x: float) -> str:
    if abs(x) < 1e-12:
        return "0"
    if abs(x - 1) < 1e-12:
        return "1"
    return f"{x:.15g}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("output", help="the copied Mathematica output")
    p.add_argument("--dt", required=True, help="dt = 4*beta/Lt the weights were computed for")
    p.add_argument("--vary", required=True, choices=["mu3", "mu8"], help="potential given on each block")
    p.add_argument("--fixed", default="0", help="value of the other potential (default 0)")
    args = p.parse_args()
    with open(args.output, encoding="utf-8") as f:
        blocks = parse(f.read())
    for mu, w in blocks:
        mu3, mu8 = (mu, args.fixed) if args.vary == "mu3" else (args.fixed, mu)
        print(" ".join([args.dt, fmt(float(mu3)), fmt(float(mu8))] + [fmt(x) for x in w]))


if __name__ == "__main__":
    main()
