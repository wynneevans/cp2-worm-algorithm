"""Estimate how many initial measurements to discard, from M_*.dat.

1. Skip scan: chi3 and M3 recomputed with the first `s` rows dropped, for increasing `s`, plotted with
   the exact value (if tests/exact has one). The estimate should reach a plateau once s is past the transient.
2. MSER-5: for each column, the truncation d (a multiple of 5) minimising the variance of the mean of the
   remaining 5-measurement block means, searched over the first half of the run. Constant columns are skipped.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns

from analyse_m import CHI_NORM, M_NORM, block_stats, load_exact, parse_name, read_m
from plot_m_raw import BLUE, SERIES

MSER_BATCH = 5
ORANGE = "#D55E00"


def mser(x, batch=MSER_BATCH):
    """MSER-`batch` truncation point (in measurements), or None for a constant series."""
    x = np.asarray(x, dtype=float)
    if not np.any(x - x[0]):
        return None
    m = len(x) // batch
    y = x[: m * batch].reshape(m, batch).mean(axis=1)
    k = np.arange(m, 0, -1)                       # points kept when dropping d = 0..m-1
    s1 = np.cumsum(y[::-1])[::-1]                 # suffix sums
    s2 = np.cumsum((y ** 2)[::-1])[::-1]
    stat = (s2 - s1 ** 2 / k) / k ** 2
    d = int(np.argmin(stat[: m // 2]))
    return d * batch


def skip_values(n, nblocks):
    """0 and log-spaced skips, up to half the run (and at least 100 rows per block)."""
    top = min(n // 2, n - 100 * nblocks)
    s = [0] + [int(v) for v in np.unique(np.round(np.logspace(1, np.log10(top), 14)))]
    return [v for v in s if v <= top]


def skip_scan(df, p, nblocks):
    q3 = df["q3"].to_numpy()
    n = p["lx"] * p["ly"]
    rows = []
    for s in skip_values(len(q3), nblocks):
        m, dm = block_stats(q3[s:], nblocks)
        c, dc = block_stats(q3[s:] ** 2, nblocks)
        sm, sc = 1 / (M_NORM * p["lx"]), p["beta"] / (CHI_NORM * n)
        rows.append(dict(skip=s, M3=m * sm, dM3=dm * sm, chi3=c * sc, dchi3=dc * sc))
    return pl.DataFrame(rows)


def analyse_file(path, exact, nblocks, outdir):
    p = parse_name(path)
    df = read_m(path)
    scan = skip_scan(df, p, nblocks)
    lat = f"{p['lx']}x{p['ly']}"
    ex = {}
    for q, key in (("susceptibility", "chi3"), ("magnetization", "M3")):
        v = exact.filter((pl.col("lattice") == lat) & (pl.col("boundary") == p["bc"]) & (pl.col("quantity") == q)
                         & ((pl.col("beta_label") - p["beta"]).abs() < 1e-9)
                         & ((pl.col("mu3") - p["mu3"]).abs() < 1e-9) & ((pl.col("mu8") - p["mu8"]).abs() < 1e-9))
        ex[key] = v["value"][0] if v.height else None

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, key, label in zip(axes, ("chi3", "M3"), ("susceptibility", "magnetisation")):
        ax.errorbar(scan["skip"].to_numpy(), scan[key].to_numpy(), yerr=scan["d" + key].to_numpy(), fmt="o-",
                    ms=5, lw=1, capsize=3, color=BLUE, label="Monte Carlo")
        if ex[key] is not None:
            ax.axhline(ex[key], color=ORANGE, lw=1.5, ls="--", label="exact")
        ax.set_xscale("symlog", linthresh=10)
        ax.set_xlim(left=0)
        ax.set_xlabel("measurements skipped")
        ax.set_ylabel(label)
        ax.legend(frameon=False)
    fig.suptitle(f"{lat} {p['bc']}, beta={p['beta']}: estimate vs number of initial measurements skipped")
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(outdir / f"skipscan_{lat}_{p['bc']}.png", dpi=130)
    plt.close(fig)

    rows = [dict(lattice=lat, boundary=p["bc"], column=c, mser5_skip=mser(df[c].to_numpy())) for c in SERIES]
    return pl.DataFrame(rows, schema_overrides={"mser5_skip": pl.Int64})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="*", help="M_*.dat files or directories (default: ../data/exact)")
    ap.add_argument("--exact", default="../tests/exact/exactCP2_results.csv")
    ap.add_argument("--nblocks", type=int, default=50)
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
    exact = load_exact(a.exact)
    table = pl.concat([analyse_file(f, exact, a.nblocks, out) for f in files])
    table.write_csv(out / "mser.csv")
    with pl.Config(tbl_rows=-1):
        print(table)
    print(f"wrote {out}/mser.csv and skipscan_*.png")


if __name__ == "__main__":
    main()
