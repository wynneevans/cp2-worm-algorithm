"""Autocorrelation function and integrated autocorrelation time of the M_*.dat columns.

tau_int = 1/2 + sum_{t=1}^{W} rho(t), with the window W chosen automatically (Sokal): the smallest W
with W >= C * tau_int(W). Error: tau * sqrt((4W + 2) / N) (Madras-Sokal). Constant columns are skipped.
Intended use: compare tau_int with the block size used for errors (nblocks should give blocks >> 2*tau_int).
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

from analyse_m import parse_name, read_m
from plot_m_raw import BLUE, SERIES

WINDOW_C = 6.0


def acf(x):
    """Normalised autocorrelation function rho(0..N-1) via FFT, or None for a constant series."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    if not np.any(x):
        return None
    f = np.fft.rfft(x, 2 * n)
    c = np.fft.irfft(f * np.conj(f))[:n]
    return c / c[0]


def tau_int(rho, c=WINDOW_C, n=None):
    """(tau_int, error, window) with Sokal's automatic windowing."""
    n = n or len(rho)
    tau = 0.5 + np.cumsum(rho[1:])
    tau = np.concatenate([[0.5], tau])
    w = np.arange(len(tau))
    ok = np.nonzero(w >= c * tau)[0]
    win = int(ok[0]) if len(ok) else len(tau) - 1
    t = tau[win]
    return t, t * np.sqrt((4 * win + 2) / n), win


def analyse_file(path, skip, maxlag, outdir):
    p = parse_name(path)
    df = read_m(path)
    df = df.slice(int(skip * df.height))
    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(3, 3, figsize=(12, 8), sharex=True)
    rows = []
    for ax, col in zip(axes.flat, SERIES):
        rho = acf(df[col].to_numpy())
        ax.set_title(col, fontsize=10)
        if rho is None:
            ax.text(0.5, 0.5, "constant", ha="center", va="center", transform=ax.transAxes, color="0.4")
            rows.append(dict(column=col, tau_int=None, tau_err=None, window=None))
            continue
        t, e, w = tau_int(rho, n=df.height)
        lags = np.arange(maxlag + 1)
        ax.axhline(0, color="0.6", lw=0.8)
        sns.lineplot(x=lags, y=rho[: maxlag + 1], ax=ax, color=BLUE, lw=1.2)
        ax.axvline(min(w, maxlag), color="0.55", lw=0.8, ls="--")
        ax.text(0.97, 0.95, f"tau_int = {t:.3g} ± {e:.1g}", ha="right", va="top", transform=ax.transAxes,
                fontsize=9)
        rows.append(dict(column=col, tau_int=t, tau_err=e, window=w))
    for ax in axes[-1]:
        ax.set_xlabel("lag (measurements)")
    for ax in axes[:, 0]:
        ax.set_ylabel("autocorrelation")
    fig.suptitle(f"{p['lx']}x{p['ly']} {p['bc']}, beta={p['beta']}: autocorrelation "
                 f"(dashed: Sokal window, C={WINDOW_C:g})")
    sns.despine(fig)
    fig.tight_layout()
    stem = f"acf_{p['lx']}x{p['ly']}_{p['bc']}"
    fig.savefig(outdir / f"{stem}.png", dpi=130)
    plt.close(fig)
    table = pl.DataFrame(rows).with_columns(pl.lit(f"{p['lx']}x{p['ly']}").alias("lattice"),
                                            pl.lit(p["bc"]).alias("boundary"))
    return table.select("lattice", "boundary", "column", "tau_int", "tau_err", "window")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="*", help="M_*.dat files or directories (default: ../data/exact)")
    ap.add_argument("--skip", type=float, default=0.01, help="fraction of leading rows to discard")
    ap.add_argument("--maxlag", type=int, default=50, help="largest lag drawn")
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
    table = pl.concat([analyse_file(f, a.skip, a.maxlag, out) for f in files])
    table.write_csv(out / "tau_int.csv")
    with pl.Config(tbl_rows=-1):
        print(table)
    print(f"wrote {out}/tau_int.csv and acf_*.png")


if __name__ == "__main__":
    main()
