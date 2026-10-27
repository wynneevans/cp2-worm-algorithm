"""Plot block-averaged raw M_*.dat columns to inspect thermalisation and autocorrelation.

Every `--block` consecutive measurements are averaged into one point. The `tail` column is a
string label, not a series, so it is not plotted.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns

from analyse_m import COLUMNS, parse_name, read_m

SERIES = COLUMNS[:9]
BLUE = "#0072B2"


def block_average(df, block):
    """One row per full block: x = measurement number at the block centre, then each column's mean."""
    nfull = df.height // block
    return (
        df.head(nfull * block)
        .with_row_index("i")
        .group_by(pl.col("i") // block, maintain_order=True)
        .agg([pl.col("i").mean().alias("x")] + [pl.col(c).mean() for c in SERIES])
        .drop("i")
    )


def plot_file(path, block, outdir):
    p = parse_name(path)
    df = read_m(path)
    avg = block_average(df, block)
    x = avg["x"].to_numpy()
    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(len(SERIES), 1, figsize=(11, 1.5 * len(SERIES)), sharex=True)
    for ax, col in zip(axes, SERIES):
        ax.axhline(df[col].mean(), color="0.55", lw=0.8, ls="--")
        sns.lineplot(x=x, y=avg[col].to_numpy(), ax=ax, color=BLUE, lw=0.8)
        ax.set_ylabel(col, rotation=0, ha="right", va="center")
    axes[-1].set_xlabel("measurement number")
    fig.suptitle(f"{p['lx']}x{p['ly']} {p['bc']}, beta={p['beta']}: mean of every {block} measurements "
                 f"(dashed: whole-run mean)")
    sns.despine(fig)
    fig.tight_layout()
    out = outdir / f"raw_{p['lx']}x{p['ly']}_{p['bc']}.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out, avg.height


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="*", help="M_*.dat files or directories (default: ../data/exact)")
    ap.add_argument("--block", type=int, default=100, help="measurements averaged per plotted point")
    ap.add_argument("--out", default="../data/exact/plots")
    a = ap.parse_args(argv)

    files = []
    for i in a.inputs or ["../data/exact"]:
        p = Path(i)
        files += sorted(p.glob("M_*.dat")) if p.is_dir() else [p]
    if not files:
        sys.exit("error: no M_*.dat files found")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for f in files:
        png, n = plot_file(f, a.block, out)
        print(f"{png} ({n} points per column)")


if __name__ == "__main__":
    main()
