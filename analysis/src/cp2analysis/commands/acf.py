"""Autocorrelation function and integrated autocorrelation time of the M_*.dat columns."""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
import seaborn as sns

from .. import plotting
from ..runs import SERIES, RunInfo, read_m, run_info
from ..stats import SOKAL_C, acf, tau_int
from . import _common


def add_arguments(p: argparse.ArgumentParser) -> None:
    _common.add_io_arguments(p)
    _common.add_skip_argument(p)
    p.add_argument("--maxlag", type=int, default=50, help="largest lag drawn")


def compute(run: RunInfo, df: pl.DataFrame) -> tuple[pl.DataFrame, dict[str, np.ndarray | None]]:
    """tau_int table (constant columns are null) and the ACF of each column."""
    rows, rhos = [], {}
    for col in SERIES:
        rho = acf(df[col].to_numpy())
        rhos[col] = rho
        if rho is None:
            rows.append({"column": col, "tau_int": None, "tau_err": None, "window": None})
        else:
            t, e, w = tau_int(rho, n=df.height)
            rows.append({"column": col, "tau_int": t, "tau_err": e, "window": w})
    table = pl.DataFrame(
        rows, schema_overrides={"tau_int": pl.Float64, "tau_err": pl.Float64, "window": pl.Int64}
    )
    table = table.with_columns(
        pl.lit(run.name).alias("name"), pl.lit(run.lattice).alias("lattice"), pl.lit(run.bc).alias("boundary")
    )
    return table.select("name", "lattice", "boundary", "column", "tau_int", "tau_err", "window"), rhos


def plot(run: RunInfo, table: pl.DataFrame, rhos: dict, maxlag: int, outdir: Path) -> None:
    plotting.apply_theme()
    fig, axes = plotting.new_figure(3, 3, figsize=(12, 8), sharex=True)
    for ax, row in zip(axes.flat, table.iter_rows(named=True), strict=False):
        col, rho = row["column"], rhos[row["column"]]
        ax.set_title(col, fontsize=10)
        if rho is None:
            ax.text(0.5, 0.5, "constant", ha="center", va="center", transform=ax.transAxes, color="0.4")
            continue
        ax.axhline(0, color="0.6", lw=0.8)
        sns.lineplot(x=np.arange(maxlag + 1), y=rho[: maxlag + 1], ax=ax, color=plotting.BLUE, lw=1.2)
        ax.axvline(min(row["window"], maxlag), color="0.55", lw=0.8, ls="--")
        ax.text(
            0.97,
            0.95,
            f"tau_int = {row['tau_int']:.3g} ± {row['tau_err']:.1g}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=9,
        )
    for ax in axes[-1]:
        ax.set_xlabel("lag (measurements)")
    for ax in axes[:, 0]:
        ax.set_ylabel("autocorrelation")
    fig.suptitle(
        f"{run.name} ({run.lattice} {run.bc}, beta={run.beta}): "
        f"autocorrelation (dashed: Sokal window, C={SOKAL_C:g})"
    )
    plotting.save(fig, outdir / f"acf_{run.name}.png")


def run(args: argparse.Namespace) -> None:
    files, out = _common.resolve(args)
    tables = []
    for f in files:
        info = run_info(f)
        df = read_m(f)
        df = df.slice(int(args.skip * df.height))
        table, rhos = compute(info, df)
        plot(info, table, rhos, args.maxlag, out)
        tables.append(table)
    table = pl.concat(tables)
    table.write_csv(out / "tau_int.csv")
    with pl.Config(tbl_rows=-1):
        print(table)
    print(f"wrote {out}/tau_int.csv and acf_*.png")
