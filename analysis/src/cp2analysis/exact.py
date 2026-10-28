"""Exact results (tests/exact/exactCP2_results.csv)."""

from pathlib import Path

import polars as pl

from .runs import RunInfo

QUANTITIES = {"chi": "susceptibility", "M": "magnetization"}


def load_exact(path: str | Path) -> pl.DataFrame:
    """Labelled-file rows, de-duplicated (the 'plotted data' rows repeat them)."""
    df = pl.read_csv(path, schema_overrides={"beta_label": pl.Float64, "mu3": pl.Float64, "mu8": pl.Float64})
    df = df.filter(pl.col("status") == "labelled file")
    df = df.unique(subset=["lattice", "boundary", "beta_label", "mu3", "mu8", "quantity"], keep="first")
    return df.select("lattice", "boundary", "beta_label", "mu3", "mu8", "quantity", "value")


def exact_value(exact: pl.DataFrame, run: RunInfo, quantity: str) -> float | None:
    """The exact `quantity` ('susceptibility' or 'magnetization') for this run, or None."""
    v = exact.filter(
        (pl.col("lattice") == run.lattice)
        & (pl.col("boundary") == run.bc)
        & (pl.col("quantity") == quantity)
        & ((pl.col("beta_label") - run.beta).abs() < 1e-9)
        & ((pl.col("mu3") - run.mu3).abs() < 1e-9)
        & ((pl.col("mu8") - run.mu8).abs() < 1e-9)
    )
    return float(v["value"][0]) if v.height else None
