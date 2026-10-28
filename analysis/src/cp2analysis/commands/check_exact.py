"""Pass/fail check of M and chi against the exact results (used by the exact_* ctest tests).

Each parameter folder is pooled as in `summary` (the first --skip fraction of each run dropped). The
charge with the chemical potential is compared, or both charges when mu3 = mu8 = 0 (their exact values
are then equal). Exits with status 1 if any value is more than --max-pull errors from the exact one, or
if a selected parameter set has no exact value.
"""

import argparse

import polars as pl

from ..exact import QUANTITIES, exact_value, load_exact
from ..filters import add_filter_arguments
from ..runs import RunInfo, read_m
from . import _common
from .summary import compute as summarise_groups
from .summary import group_runs


def add_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "inputs",
        nargs="*",
        help=f"parameter folders, M_<time>.dat files or folders with runs (default: {_common.DEFAULT_INPUT})",
    )
    add_filter_arguments(p)
    _common.add_exact_argument(p)
    _common.add_skip_argument(p)
    _common.add_nblocks_argument(p)
    p.add_argument("--max-pull", type=float, default=4.0, help="largest |MC - exact| / error that passes")
    p.add_argument(
        "--typed",
        action="store_true",
        help="also use exact values typed into the notebook (see exact.py)",
    )


def compare(summary: pl.DataFrame, runs: list[RunInfo], exact: pl.DataFrame, max_pull: float) -> pl.DataFrame:
    """One row per compared value: label, quantity, value, error, exact, pull, ok.

    `runs` holds one run per summary row, giving its parameters.
    """
    rows = []
    for s, run in zip(summary.iter_rows(named=True), runs, strict=True):
        charges = [c for c, mu in (("3", run.mu3), ("8", run.mu8)) if mu != 0] or ["3", "8"]
        for c in charges:
            for key, quantity in QUANTITIES.items():
                value, error = s[f"{key}{c}"], s[f"d{key}{c}"]
                ex = exact_value(exact, run, quantity)
                if ex is None:
                    pull, ok = None, False
                elif error > 0:
                    pull = (value - ex) / error
                    ok = abs(pull) <= max_pull
                else:  # constant data: only an exact match passes
                    pull, ok = None, abs(value - ex) < 1e-12
                rows.append(
                    {
                        "label": s["label"],
                        "quantity": f"{key}{c}",
                        "value": value,
                        "error": error,
                        "exact": ex,
                        "pull": pull,
                        "ok": ok,
                    }
                )
    return pl.DataFrame(rows, schema_overrides={"exact": pl.Float64, "pull": pl.Float64})


def run(args: argparse.Namespace) -> None:
    runs = _common.select(args)
    exact = load_exact(args.exact, typed=args.typed)
    groups = [[(r, read_m(r.path)) for r in g] for g in group_runs(runs, per_run=False)]
    summary = summarise_groups(groups, exact, args.skip, args.nblocks)
    result = compare(summary, [g[0][0] for g in groups], exact, args.max_pull)
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=200, float_precision=8):
        print(result)
    bad = result.filter(~pl.col("ok"))
    if bad.height:
        raise SystemExit(
            f"error: {bad.height} of {result.height} values disagree with the exact results "
            f"(more than {args.max_pull} errors, or no exact value)"
        )
    print(f"all {result.height} values within {args.max_pull} errors of the exact results")
