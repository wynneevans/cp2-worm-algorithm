"""Defaults and options shared by the subcommands (paths are relative to the working directory)."""

import argparse
import os
from pathlib import Path

from ..filters import active_filters, add_filter_arguments, apply_filters
from ..runs import RunInfo, find_runs

DEFAULT_INPUT = "data"
DEFAULT_EXACT = "tests/exact/exactCP2_results.csv"


def add_io_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "inputs",
        nargs="*",
        help=f"parameter folders, M_<time>.dat files or folders searched for runs (default: {DEFAULT_INPUT})",
    )
    p.add_argument(
        "--out",
        help="output directory (default: <parameter folder>/plots for one parameter set, "
        "otherwise plots/ in the folder that holds them)",
    )
    add_filter_arguments(p)


def add_exact_argument(p: argparse.ArgumentParser) -> None:
    p.add_argument("--exact", default=DEFAULT_EXACT, help="exact results CSV (default: %(default)s)")


def add_skip_argument(p: argparse.ArgumentParser) -> None:
    p.add_argument("--skip", type=float, default=0.01, help="fraction of each run's leading rows to discard")


def add_nblocks_argument(p: argparse.ArgumentParser) -> None:
    p.add_argument("--nblocks", type=int, default=50, help="jackknife blocks")


def default_out(runs: list[RunInfo]) -> Path:
    folders = {r.path.parent for r in runs}
    if len(folders) == 1:
        return folders.pop() / "plots"
    return Path(os.path.commonpath([f.parent for f in folders])) / "plots"


def select(args: argparse.Namespace) -> list[RunInfo]:
    """The runs chosen by the inputs and filters; exits with an error if there are none."""
    inputs = args.inputs or [DEFAULT_INPUT]
    runs = find_runs(inputs)
    filters = active_filters(args)
    runs = apply_filters(runs, filters)
    if not runs:
        what = " ".join(f"--{k} {v}" for k, v in filters.items())
        raise SystemExit(
            f"error: no finished runs in {' '.join(inputs)}" + (f" matching {what}" if what else "")
        )
    return runs


def resolve(args: argparse.Namespace) -> tuple[list[RunInfo], Path]:
    """The selected runs and the (created) output directory."""
    runs = select(args)
    out = Path(args.out) if args.out else default_out(runs)
    out.mkdir(parents=True, exist_ok=True)
    return runs, out
