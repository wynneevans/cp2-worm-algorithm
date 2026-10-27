"""Magnetisation and susceptibility from M_*.dat, plotted against exact results.

M_*.dat columns: w3x, w8x, w3y, w8y, q3, q8, l, lp, lm, tail. q3 and q8 are the
temporal winding charges in doubled units (q3 = 2*Q3, q8 = (2/sqrt(3))*Q8).
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns

COLUMNS = ["w3x", "w8x", "w3y", "w8y", "q3", "q8", "l", "lp", "lm", "tail"]

# Normalisations (see README): M = <q>/(2 Lx), chi = beta <q^2> / (4 Lx Ly).
M_NORM = 2.0
CHI_NORM = 4.0
Q8_M_FACTOR = np.sqrt(3.0)   # w9 = sqrt(3) <q8>
Q8_CHI_FACTOR = 3.0          # w6 = 3 <q8^2>

NAME_RE = re.compile(
    r"BC_(?P<bc>[^_]+)_nMeas_(?P<nmeas>\d+)_beta_(?P<beta>[\d.]+)_XYT_(?P<lx>\d+)_(?P<ly>\d+)_(?P<lt>\d+)"
    r"_mu3_(?P<mu3>-?[\d.]+)_mu8_(?P<mu8>-?[\d.]+)_nSim_(?P<sim>\d+)"
)
OKABE_BLUE, OKABE_ORANGE = "#0072B2", "#D55E00"


def parse_name(path):
    m = NAME_RE.search(Path(path).name)
    if not m:
        raise ValueError(f"cannot parse run parameters from '{path}'")
    d = m.groupdict()
    return dict(bc=d["bc"], beta=float(d["beta"]), lx=int(d["lx"]), ly=int(d["ly"]), lt=int(d["lt"]),
                mu3=float(d["mu3"]), mu8=float(d["mu8"]))


def read_m(path):
    df = pl.read_csv(path, has_header=False, new_columns=COLUMNS, infer_schema=False)
    return df.with_columns([pl.col(c).str.strip_chars().cast(pl.Float64) for c in COLUMNS[:9]])


def block_stats(x, nblocks):
    """Mean and standard error from non-overlapping block jackknife."""
    x = np.asarray(x, dtype=float)
    size = len(x) // nblocks
    blocks = x[: size * nblocks].reshape(nblocks, size).mean(axis=1)
    mean = blocks.mean()
    jack = (blocks.sum() - blocks) / (nblocks - 1)
    err = np.sqrt((nblocks - 1) / nblocks * ((jack - jack.mean()) ** 2).sum())
    return mean, err


def summarise(path, skip_frac=0.01, nblocks=50):
    p = parse_name(path)
    df = read_m(path)
    df = df.slice(int(skip_frac * df.height))
    n = p["lx"] * p["ly"]
    out = dict(lattice=f"{p['lx']}x{p['ly']}", boundary=p["bc"], beta=p["beta"], mu3=p["mu3"], mu8=p["mu8"],
               n_used=df.height)
    q3, q8 = df["q3"].to_numpy(), df["q8"].to_numpy()
    for name, q, mfac, cfac in (("3", q3, 1.0, 1.0), ("8", q8, Q8_M_FACTOR, Q8_CHI_FACTOR)):
        m, dm = block_stats(q, nblocks)
        c, dc = block_stats(q ** 2, nblocks)
        sm = mfac / (M_NORM * p["lx"])
        sc = cfac * p["beta"] / (CHI_NORM * n)
        out.update({f"M{name}": m * sm, f"dM{name}": dm * sm, f"chi{name}": c * sc, f"dchi{name}": dc * sc})
    return out


def load_exact(path):
    df = pl.read_csv(path, schema_overrides={"beta_label": pl.Float64, "mu3": pl.Float64, "mu8": pl.Float64})
    df = df.filter(pl.col("status") == "labelled file")
    df = df.unique(subset=["lattice", "boundary", "beta_label", "mu3", "mu8", "quantity"], keep="first")
    return df.select("lattice", "boundary", "beta_label", "mu3", "mu8", "quantity", "value")


def attach_exact(summary, exact):
    for q, col in (("susceptibility", "exact_chi"), ("magnetization", "exact_M")):
        e = exact.filter(pl.col("quantity") == q).select(
            pl.col("lattice"), pl.col("boundary"), pl.col("beta_label").alias("beta"), "mu3", "mu8",
            pl.col("value").alias(col))
        summary = summary.join(e, on=["lattice", "boundary", "beta", "mu3", "mu8"], how="left")
    return summary


def plot(summary, outdir):
    df = summary.with_columns((pl.col("lattice") + " " + pl.col("boundary")).alias("case"))
    cases = df["case"].to_list()
    x = np.arange(len(cases))
    sns.set_theme(style="whitegrid", context="talk")
    for val, err, exact, ylabel, fname in (
        ("chi3", "dchi3", "exact_chi", "susceptibility", "susceptibility.png"),
        ("M3", "dM3", "exact_M", "magnetisation", "magnetisation.png"),
    ):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.errorbar(x, df[val].to_numpy(), yerr=df[err].to_numpy(), fmt="o", ms=9, capsize=4, lw=2,
                    color=OKABE_BLUE, label="Monte Carlo")
        sns.scatterplot(x=x, y=df[exact].to_numpy(), marker="D", s=90, color=OKABE_ORANGE, ax=ax,
                        label="exact", zorder=5)
        ax.set_xticks(x, cases)
        ax.set_xlim(-0.5, len(cases) - 0.5)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        if val == "M3":
            ax.axhline(0, color="0.6", lw=1)
        ax.legend(frameon=False)
        sns.despine(fig)
        fig.tight_layout()
        outdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir / fname, dpi=150)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("inputs", nargs="*", help="M_*.dat files or directories (default: ../data/exact)")
    ap.add_argument("--exact", default="../tests/exact/exactCP2_results.csv")
    ap.add_argument("--skip", type=float, default=0.01, help="fraction of leading rows to discard")
    ap.add_argument("--nblocks", type=int, default=50)
    ap.add_argument("--out", default="../data/exact/plots")
    a = ap.parse_args(argv)

    files = []
    for i in a.inputs or ["../data/exact"]:
        p = Path(i)
        files += sorted(p.glob("M_*.dat")) if p.is_dir() else [p]
    if not files:
        sys.exit("error: no M_*.dat files found")

    summary = pl.DataFrame([summarise(f, a.skip, a.nblocks) for f in files])
    summary = attach_exact(summary, load_exact(a.exact))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    summary.write_csv(out / "m_summary.csv")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250):
        print(summary)
    plot(summary, out)
    print(f"wrote {out}/m_summary.csv, susceptibility.png, magnetisation.png")


if __name__ == "__main__":
    main()
