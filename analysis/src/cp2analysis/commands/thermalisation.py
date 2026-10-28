"""How many initial measurements to discard, and whether the run is stationary.

Top row: chi3 and M3 with the first `skip` measurements dropped (cumulative: neighbouring points share
most of their data). Bottom row: chi3 and M3 on --chunks non-overlapping parts of the run, with
chi^2/dof of the chunks about their mean (about 1 if they agree within errors). MSER-5 per column.
"""

import argparse
from pathlib import Path

import polars as pl

from .. import plotting
from ..exact import QUANTITIES, exact_value, load_exact
from ..observables import chunk_scan, consistency, skip_scan
from ..runs import SERIES, RunInfo, read_m
from ..stats import mser
from . import _common


def add_arguments(p: argparse.ArgumentParser) -> None:
    _common.add_io_arguments(p)
    _common.add_exact_argument(p)
    _common.add_nblocks_argument(p)
    p.add_argument("--chunks", type=int, default=8, help="non-overlapping parts for the bottom row")


PANELS = (("chi3", "susceptibility"), ("M3", "magnetisation"))


def plot_scan(
    run: RunInfo, scan: pl.DataFrame, chunks: pl.DataFrame, exact: pl.DataFrame, outdir: Path
) -> None:
    plotting.apply_theme()
    fig, axes = plotting.new_figure(2, 2, figsize=(11, 7.5))
    for col, (key, label) in enumerate(PANELS):
        value = exact_value(exact, run, QUANTITIES[key.removesuffix("3")])

        ax = axes[0, col]
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
        if value is not None:
            ax.axhline(value, color=plotting.ORANGE, lw=1.5, ls="--", label="exact")
        ax.set_xscale("symlog", linthresh=10)
        ax.set_xlim(left=0)
        ax.set_xlabel("measurements skipped")
        ax.set_ylabel(label)
        ax.set_title("cumulative: all measurements after the skip", fontsize=10)
        ax.legend(frameon=False)

        ax = axes[1, col]
        start, end = chunks["start"].to_numpy(), chunks["end"].to_numpy()
        v, e = chunks[key].to_numpy(), chunks["d" + key].to_numpy()
        ax.axhline(scan[key][0], color="0.55", lw=0.8, ls="--", label="whole run")
        plotting.errorbar(
            ax, (start + end) / 2, v, e, xerr=(end - start) / 2, fmt="o", ms=5, lw=1, capsize=3, label="chunk"
        )
        if value is not None:
            ax.axhline(value, color=plotting.ORANGE, lw=1.5, ls="--", label="exact")
        c2 = consistency(v, e)
        if c2 is not None:
            ax.text(
                0.98, 0.95, f"χ²/dof = {c2:.2f}", ha="right", va="top", transform=ax.transAxes, fontsize=10
            )
        ax.set_xlabel("measurement number")
        ax.set_ylabel(label)
        ax.set_title(f"{chunks.height} separate chunks (should agree within errors)", fontsize=10)
        ax.legend(frameon=False, loc="lower left", fontsize="small")
    fig.suptitle(f"{run.label} ({run.size} {run.bc}, beta={run.beta}, mu3={run.mu3})")
    plotting.save(fig, outdir / f"skipscan_{run.label}.png")


def chunk_table(run: RunInfo, chunks: pl.DataFrame) -> pl.DataFrame:
    """The chunks with the run's label and chi^2/dof per quantity."""
    extra = {
        f"chi2dof_{key}": consistency(chunks[key].to_numpy(), chunks["d" + key].to_numpy())
        for key, _ in PANELS
    }
    return chunks.select(
        pl.lit(run.label).alias("label"),
        pl.lit(run.folder).alias("folder"),
        pl.all(),
        *[pl.lit(v, dtype=pl.Float64).alias(k) for k, v in extra.items()],
    )


def compute_mser(run: RunInfo, df: pl.DataFrame) -> pl.DataFrame:
    rows = [
        {
            "label": run.label,
            "folder": run.folder,
            "column": c,
            "mser5_skip": mser(df[c].to_numpy()),
        }
        for c in SERIES
    ]
    return pl.DataFrame(rows, schema_overrides={"mser5_skip": pl.Int64})


def run(args: argparse.Namespace) -> None:
    runs, out = _common.resolve(args)
    exact = load_exact(args.exact)
    tables, chunk_tables = [], []
    for info in runs:
        df = read_m(info.path)
        chunks = chunk_scan(info, df, args.chunks, args.nblocks)
        plot_scan(info, skip_scan(info, df, args.nblocks), chunks, exact, out)
        tables.append(compute_mser(info, df))
        chunk_tables.append(chunk_table(info, chunks))
    table = pl.concat(tables)
    table.write_csv(out / "mser.csv")
    chunk_all = pl.concat(chunk_tables)
    chunk_all.write_csv(out / "chunks.csv")
    with pl.Config(tbl_rows=-1):
        print(table)
        print(
            chunk_all.group_by("label", maintain_order=True).agg(pl.col("chi2dof_chi3", "chi2dof_M3").first())
        )
    print(f"wrote {out}/mser.csv, chunks.csv and skipscan_*.png")
    print("chi2/dof: about 1 if the chunks agree; well above 1 means drift or too-small errors")
