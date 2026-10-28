"""Magnetisation and susceptibility from M_*.dat, plotted against a parameter or per result.

By default the runs in each parameter folder are pooled: each run's first --skip fraction is dropped,
the rest of all runs is joined and block-jackknifed. --per-run gives one result per run instead.
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
import seaborn as sns

from .. import plotting
from ..exact import QUANTITIES, exact_value, load_exact
from ..observables import summarise
from ..runs import RunInfo, read_m
from . import _common

PLOTS = (  # value, error, exact column, axis label, file name
    ("chi3", "dchi3", "exact_chi", "susceptibility", "susceptibility.png"),
    ("M3", "dM3", "exact_M", "magnetisation", "magnetisation.png"),
)
# Parameters that identify a parameter set, in the order used for series names.
PARAMS = ("bc", "lx", "ly", "lt", "measurements", "beta", "mu3", "mu8")
# Default x when several parameters vary: the first of these that varies; the others become series.
X_PREFERENCE = ("mu3", "mu8", "beta", "ly", "lx", "lt", "measurements")
X_CHOICES = ("beta", "mu3", "mu8", "lx", "ly", "lt", "measurements", "run")


def add_arguments(p: argparse.ArgumentParser) -> None:
    _common.add_io_arguments(p)
    _common.add_exact_argument(p)
    _common.add_skip_argument(p)
    _common.add_nblocks_argument(p)
    p.add_argument(
        "--per-run", action="store_true", help="one result per run instead of pooling each parameter folder"
    )
    p.add_argument(
        "--x",
        choices=X_CHOICES,
        help="x axis (default: first varying of mu3, mu8, beta, ly, lx, lt, measurements; else per result)",
    )


def group_runs(runs: list[RunInfo], per_run: bool) -> list[list[RunInfo]]:
    """Runs analysed together: one group per parameter folder, or one per run."""
    if per_run:
        return [[r] for r in runs]
    groups: dict[Path, list[RunInfo]] = {}
    for r in runs:
        groups.setdefault(r.path.parent, []).append(r)
    return list(groups.values())


def compute(
    groups: list[list[tuple[RunInfo, pl.DataFrame]]], exact: pl.DataFrame, skip: float, nblocks: int
) -> pl.DataFrame:
    """One row per group: the first `skip` fraction of each run dropped, the rest pooled."""
    rows = []
    for group in groups:
        run = group[0][0]
        used = pl.concat([df.slice(int(skip * df.height)) for _, df in group])
        row = {
            "label": group[0][0].label
            if len(group) == 1
            else "+".join(dict.fromkeys(r.name for r, _ in group)),
            "folder": run.folder,
            "n_runs": len(group),
            "bc": run.bc,
            "lattice": run.size,
            "lx": run.lx,
            "ly": run.ly,
            "lt": run.lt,
            "measurements": run.measurements,
            "beta": run.beta,
            "mu3": run.mu3,
            "mu8": run.mu8,
            **summarise(run, used, nblocks),
        }
        for key, quantity in QUANTITIES.items():
            row[f"exact_{key}"] = exact_value(exact, run, quantity)
        rows.append(row)
    return pl.DataFrame(rows, schema_overrides={"exact_chi": pl.Float64, "exact_M": pl.Float64})


def varying(summary: pl.DataFrame) -> list[str]:
    return [c for c in PARAMS if summary[c].n_unique() > 1]


def choose_x(summary: pl.DataFrame, x: str | None) -> str:
    """The requested x, else the first varying parameter in X_PREFERENCE, else "run"."""
    if x:
        return x
    v = varying(summary)
    return next((p for p in X_PREFERENCE if p in v), "run")


def dodge(x: np.ndarray) -> np.ndarray:
    """Spread points sharing an x value a little, so their error bars stay apart."""
    span = float(np.ptp(x)) or 1.0
    out = x.astype(float).copy()
    for v in np.unique(x):
        idx = np.flatnonzero(x == v)
        out[idx] += (np.arange(len(idx)) - (len(idx) - 1) / 2) * 0.015 * span
    return out


def series(summary: pl.DataFrame, x: str) -> list[tuple[str, pl.DataFrame]]:
    """(legend label, rows) per combination of the varying parameters other than x."""
    others = [c for c in varying(summary) if c != x] if x != "run" else []
    if not others:
        return [("Monte Carlo", summary)]
    parts = summary.partition_by(others, maintain_order=True)
    return [(", ".join(f"{c}={p[c][0]}" for c in others), p) for p in parts]


def plot(summary: pl.DataFrame, outdir: Path, x: str) -> None:
    plotting.apply_theme("talk")
    groups = series(summary, x)
    colours = plotting.palette(len(groups))
    for val, err, exact, ylabel, fname in PLOTS:
        fig, ax = plotting.new_figure(figsize=(13 if len(groups) > 1 else 8, 5))
        for (name, rows), colour in zip(groups, colours, strict=True):
            xs = np.arange(rows.height) if x == "run" else dodge(rows[x].to_numpy())
            plotting.errorbar(
                ax,
                xs,
                rows[val].to_numpy(),
                rows[err].to_numpy(),
                fmt="o",
                ms=8,
                capsize=4,
                lw=2,
                color=colour,
                label=name,
            )
        has_exact = summary.filter(pl.col(exact).is_not_null())
        if has_exact.height:
            ex = (
                np.flatnonzero(summary[exact].is_not_null().to_numpy())
                if x == "run"
                else has_exact[x].to_numpy()
            )
            sns.scatterplot(
                x=ex,
                y=has_exact[exact].to_numpy(),
                marker="D",
                s=80,
                color=plotting.ORANGE,
                ax=ax,
                label="exact",
                zorder=5,
            )
        if x == "run":
            ax.set_xticks(np.arange(summary.height), summary["label"].to_list(), rotation=30, ha="right")
            ax.set_xlim(-0.5, summary.height - 0.5)
            ax.set_xlabel("")
        else:
            ax.set_xlabel(x)
        ax.set_ylabel(ylabel)
        if val == "M3":
            ax.axhline(0, color="0.6", lw=1)
        if len(groups) > 1:  # long series names: legend beside the axes, not over the data
            ax.legend(frameon=False, fontsize="x-small", loc="upper left", bbox_to_anchor=(1.01, 1))
        else:
            ax.legend(frameon=False, fontsize="small")
        plotting.save(fig, outdir / fname, dpi=150)


def run(args: argparse.Namespace) -> None:
    runs, out = _common.resolve(args)
    groups = [[(r, read_m(r.path)) for r in g] for g in group_runs(runs, args.per_run)]
    summary = compute(groups, load_exact(args.exact), args.skip, args.nblocks)
    summary.write_csv(out / "m_summary.csv")
    shown = summary.drop("folder", "lx", "ly", "lt", "exact_M")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250):
        print(shown)
    x = choose_x(summary, args.x)
    plot(summary, out, x)
    print(f"{len(runs)} runs, {summary.height} results (x = {x}); wrote {out}/m_summary.csv, ")
    print("  susceptibility.png, magnetisation.png")
