"""Defaults and options shared by the subcommands (paths are relative to the working directory)."""

import argparse
from pathlib import Path

from ..runs import find_m_files

DEFAULT_INPUT = "data"
DEFAULT_EXACT = "tests/exact/exactCP2_results.csv"


def add_io_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "inputs", nargs="*", help=f"M files or folders searched for finished runs (default: {DEFAULT_INPUT})"
    )
    p.add_argument("--out", help="output directory (default: <first input folder>/plots)")


def add_exact_argument(p: argparse.ArgumentParser) -> None:
    p.add_argument("--exact", default=DEFAULT_EXACT, help="exact results CSV (default: %(default)s)")


def add_skip_argument(p: argparse.ArgumentParser) -> None:
    p.add_argument("--skip", type=float, default=0.01, help="fraction of leading rows to discard")


def add_nblocks_argument(p: argparse.ArgumentParser) -> None:
    p.add_argument("--nblocks", type=int, default=50, help="jackknife blocks")


def resolve(args: argparse.Namespace) -> tuple[list[Path], Path]:
    """The M_ files to process and the (created) output directory."""
    inputs = args.inputs or [DEFAULT_INPUT]
    files = find_m_files(inputs)
    if not files:
        raise SystemExit(
            "error: no M files found (folders are searched for M_<time>.dat with a run_<time>.toml)"
        )
    first = Path(inputs[0])
    out = Path(args.out) if args.out else (first if first.is_dir() else first.parent) / "plots"
    out.mkdir(parents=True, exist_ok=True)
    return files, out
