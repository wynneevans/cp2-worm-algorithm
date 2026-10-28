"""How many initial measurements to discard: skip scan of chi3 and M3, and MSER-5 per column."""

import argparse
from pathlib import Path

import polars as pl

from .. import plotting
from ..exact import QUANTITIES, exact_value, load_exact
from ..observables import skip_scan
from ..runs import SERIES, RunInfo, read_m, run_info
from ..stats import mser
from . import _common


def add_arguments(p: argparse.ArgumentParser) -> None:
    _common.add_io_arguments(p)
    _common.add_exact_argument(p)
    _common.add_nblocks_argument(p)


def plot_scan(run: RunInfo, scan: pl.DataFrame, exact: pl.DataFrame, outdir: Path) -> None:
    plotting.apply_theme()
    fig, axes = plotting.new_figure(1, 2, figsize=(11, 4))
    for ax, key, label in zip(axes, ("chi3", "M3"), ("susceptibility", "magnetisation"), strict=True):
        plotting.errorbar(
            ax,
            scan["skip"].to_numpy(),
            scan[key].to_numpy(),
            scan["d" + key].to_numpy(),
            fmt="o-",
            ms=5,
            lw=1,
            capsize=3,
            label="Monte Carlo",
        )
        value = exact_value(exact, run, QUANTITIES[key.removesuffix("3")])
        if value is not None:
            ax.axhline(value, color=plotting.ORANGE, lw=1.5, ls="--", label="exact")
        ax.set_xscale("symlog", linthresh=10)
        ax.set_xlim(left=0)
        ax.set_xlabel("measurements skipped")
        ax.set_ylabel(label)
        ax.legend(frameon=False)
    fig.suptitle(
        f"{run.name} ({run.lattice} {run.bc}, beta={run.beta}): "
        "estimate vs number of initial measurements skipped"
    )
    plotting.save(fig, outdir / f"skipscan_{run.name}.png")


def compute_mser(run: RunInfo, df: pl.DataFrame) -> pl.DataFrame:
    rows = [
        {
            "name": run.name,
            "lattice": run.lattice,
            "boundary": run.bc,
            "column": c,
            "mser5_skip": mser(df[c].to_numpy()),
        }
        for c in SERIES
    ]
    return pl.DataFrame(rows, schema_overrides={"mser5_skip": pl.Int64})


def run(args: argparse.Namespace) -> None:
    files, out = _common.resolve(args)
    exact = load_exact(args.exact)
    tables = []
    for f in files:
        info = run_info(f)
        df = read_m(f)
        plot_scan(info, skip_scan(info, df, args.nblocks), exact, out)
        tables.append(compute_mser(info, df))
    table = pl.concat(tables)
    table.write_csv(out / "mser.csv")
    with pl.Config(tbl_rows=-1):
        print(table)
    print(f"wrote {out}/mser.csv and skipscan_*.png")
