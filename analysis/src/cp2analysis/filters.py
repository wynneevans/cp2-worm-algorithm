"""Selecting runs by parameter: --beta 56.3, --mu3 0,0.02, --mu3 0.01:0.03, --lx 100, --lattice 4x4x200."""

import argparse
from collections.abc import Callable

from .runs import RunInfo

TOL = 1e-9

# option -> (run attribute, kind)
FILTERS = {
    "bc": ("bc", "text"),
    "measurements": ("measurements", "number"),
    "beta": ("beta", "number"),
    "mu3": ("mu3", "number"),
    "mu8": ("mu8", "number"),
    "lx": ("lx", "number"),
    "ly": ("ly", "number"),
    "lt": ("lt", "number"),
    "lattice": ("size", "lattice"),
}


def _number_term(term: str) -> Callable[[float], bool]:
    if ":" in term:
        lo, hi = term.split(":", 1)
        lo_v = float(lo) if lo else float("-inf")
        hi_v = float(hi) if hi else float("inf")
        return lambda v: lo_v - TOL <= v <= hi_v + TOL
    x = float(term)
    return lambda v: abs(v - x) <= TOL * max(1.0, abs(x))


def _lattice_term(term: str) -> Callable[[str], bool]:
    # "LxxLy" matches any Lt; "LxxLyxLt" matches exactly.
    parts = term.lower().split("x")
    if len(parts) not in (2, 3) or not all(p.isdigit() for p in parts):
        raise ValueError(f"lattice '{term}' should be LxxLy or LxxLyxLt")
    return lambda size: size.split("x")[: len(parts)] == parts


def parse_filter(text: str, kind: str) -> Callable:
    """A predicate for one option value: a comma list of values (numbers may be ranges a:b, a:, :b)."""
    make = {"number": _number_term, "lattice": _lattice_term, "text": lambda t: t.__eq__}[kind]
    try:
        terms = [make(t.strip()) for t in text.split(",") if t.strip()]
    except ValueError as e:
        raise SystemExit(f"error: bad filter '{text}': {e}") from None
    if not terms:
        raise SystemExit(f"error: empty filter '{text}'")
    return lambda v: any(t(v) for t in terms)


def add_filter_arguments(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("select runs by parameter (values, comma lists, ranges a:b)")
    for opt, (_, kind) in FILTERS.items():
        metavar = {"number": "V", "text": "NAME", "lattice": "LxxLy[xLt]"}[kind]
        g.add_argument(f"--{opt}", metavar=metavar)


def active_filters(args: argparse.Namespace) -> dict[str, str]:
    return {opt: getattr(args, opt) for opt in FILTERS if getattr(args, opt, None) is not None}


def apply_filters(runs: list[RunInfo], filters: dict[str, str]) -> list[RunInfo]:
    """The runs matching every given filter."""
    preds = [(FILTERS[opt][0], parse_filter(text, FILTERS[opt][1])) for opt, text in filters.items()]
    return [r for r in runs if all(pred(getattr(r, attr)) for attr, pred in preds)]
