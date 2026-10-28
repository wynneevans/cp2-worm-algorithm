"""Block-averaged raw M_*.dat columns, for inspecting thermalisation and autocorrelation by eye."""

import argparse
from pathlib import Path

import polars as pl
import seaborn as sns

from .. import plotting
from ..runs import SERIES, RunInfo, read_m
from . import _common


def add_arguments(p: argparse.ArgumentParser) -> None:
    _common.add_io_arguments(p)
    p.add_argument("--block", type=int, default=100, help="measurements averaged per plotted point")


def block_average(df: pl.DataFrame, block: int) -> pl.DataFrame:
    """One row per full block: x = measurement number at the block centre, then each column's mean."""
    nfull = df.height // block
    return (
        df.head(nfull * block)
        .with_row_index("i")
        .group_by(pl.col("i") // block, maintain_order=True)
        .agg([pl.col("i").mean().alias("x")] + [pl.col(c).mean() for c in SERIES])
        .drop("i")
    )


def plot_file(run: RunInfo, block: int, outdir: Path) -> tuple[Path, int]:
    df = read_m(run.path)
    avg = block_average(df, block)
    x = avg["x"].to_numpy()
    plotting.apply_theme()
    fig, axes = plotting.new_figure(len(SERIES), 1, figsize=(11, 1.5 * len(SERIES)), sharex=True)
    for ax, col in zip(axes, SERIES, strict=True):
        ax.axhline(df[col].mean(), color="0.55", lw=0.8, ls="--")
        sns.lineplot(x=x, y=avg[col].to_numpy(), ax=ax, color=plotting.BLUE, lw=0.8)
        ax.set_ylabel(col, rotation=0, ha="right", va="center")
    axes[-1].set_xlabel("measurement number")
    fig.suptitle(
        f"{run.label} ({run.size} {run.bc}, beta={run.beta}, mu3={run.mu3}): "
        f"mean of every {block} measurements (dashed: whole-run mean)"
    )
    out = outdir / f"raw_{run.label}.png"
    plotting.save(fig, out)
    return out, avg.height


def run(args: argparse.Namespace) -> None:
    runs, out = _common.resolve(args)
    for r in runs:
        png, n = plot_file(r, args.block, out)
        print(f"{png} ({n} points per column)")
