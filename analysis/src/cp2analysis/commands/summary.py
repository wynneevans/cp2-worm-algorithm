"""Magnetisation and susceptibility from M_*.dat, plotted against the exact results."""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
import seaborn as sns

from .. import plotting
from ..exact import QUANTITIES, exact_value, load_exact
from ..observables import summarise
from ..runs import read_m, run_info
from . import _common

PLOTS = (  # value, error, exact column, axis label, file name
    ("chi3", "dchi3", "exact_chi", "susceptibility", "susceptibility.png"),
    ("M3", "dM3", "exact_M", "magnetisation", "magnetisation.png"),
)


def add_arguments(p: argparse.ArgumentParser) -> None:
    _common.add_io_arguments(p)
    _common.add_exact_argument(p)
    _common.add_skip_argument(p)
    _common.add_nblocks_argument(p)


def compute(files: list[Path], exact: pl.DataFrame, skip: float, nblocks: int) -> pl.DataFrame:
    rows = []
    for f in files:
        run = run_info(f)
        df = read_m(f)
        row = summarise(run, df.slice(int(skip * df.height)), nblocks)
        for key, quantity in QUANTITIES.items():
            row[f"exact_{key}"] = exact_value(exact, run, quantity)
        rows.append(row)
    return pl.DataFrame(rows)


def plot(summary: pl.DataFrame, outdir: Path) -> None:
    cases = summary["name"].to_list()
    x = np.arange(len(cases))
    plotting.apply_theme("talk")
    for val, err, exact, ylabel, fname in PLOTS:
        fig, ax = plotting.new_figure(figsize=(7, 4.5))
        plotting.errorbar(
            ax,
            x,
            summary[val].to_numpy(),
            summary[err].to_numpy(),
            fmt="o",
            ms=9,
            capsize=4,
            lw=2,
            label="Monte Carlo",
        )
        sns.scatterplot(
            x=x,
            y=summary[exact].to_numpy(),
            marker="D",
            s=90,
            color=plotting.ORANGE,
            ax=ax,
            label="exact",
            zorder=5,
        )
        ax.set_xticks(x, cases)
        ax.set_xlim(-0.5, len(cases) - 0.5)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        if val == "M3":
            ax.axhline(0, color="0.6", lw=1)
        ax.legend(frameon=False)
        plotting.save(fig, outdir / fname, dpi=150)


def run(args: argparse.Namespace) -> None:
    files, out = _common.resolve(args)
    summary = compute(files, load_exact(args.exact), args.skip, args.nblocks)
    summary.write_csv(out / "m_summary.csv")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250):
        print(summary)
    plot(summary, out)
    print(f"wrote {out}/m_summary.csv, susceptibility.png, magnetisation.png")
