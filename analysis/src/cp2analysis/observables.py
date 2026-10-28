"""Magnetisation and susceptibility from the temporal winding charges q3, q8.

q3 = 2*Q3 and q8 = (2/sqrt(3))*Q8 (doubled units). With L = Lx and N = Lx*Ly:
M = <q>/(2L), chi = beta <q^2>/(4N); for q8 the factors sqrt(3) and 3 restore the Q8 normalisation.
"""

import numpy as np
import polars as pl

from .runs import RunInfo
from .stats import block_stats

M_NORM = 2.0
CHI_NORM = 4.0
Q8_M_FACTOR = float(np.sqrt(3.0))
Q8_CHI_FACTOR = 3.0

# charge column -> (M factor, chi factor)
CHARGES = {"3": ("q3", 1.0, 1.0), "8": ("q8", Q8_M_FACTOR, Q8_CHI_FACTOR)}


def scales(run: RunInfo, m_factor: float = 1.0, chi_factor: float = 1.0) -> tuple[float, float]:
    """Factors turning <q> into M and <q^2> into chi."""
    return m_factor / (M_NORM * run.lx), chi_factor * run.beta / (CHI_NORM * run.n_sites)


def summarise(run: RunInfo, df: pl.DataFrame, nblocks: int) -> dict:
    """M and chi (with errors) for both charges, one summary row."""
    out = {
        "name": run.name,
        "lattice": run.lattice,
        "boundary": run.bc,
        "beta": run.beta,
        "mu3": run.mu3,
        "mu8": run.mu8,
        "n_used": df.height,
    }
    for name, (col, mfac, cfac) in CHARGES.items():
        q = df[col].to_numpy()
        m, dm = block_stats(q, nblocks)
        c, dc = block_stats(q**2, nblocks)
        sm, sc = scales(run, mfac, cfac)
        out.update({f"M{name}": m * sm, f"dM{name}": dm * sm, f"chi{name}": c * sc, f"dchi{name}": dc * sc})
    return out


def skip_values(n: int, nblocks: int) -> list[int]:
    """0 and log-spaced skips, up to half the run (and at least 100 rows per block)."""
    top = min(n // 2, n - 100 * nblocks)
    s = [0] + [int(v) for v in np.unique(np.round(np.logspace(1, np.log10(top), 14)))]
    return [v for v in s if v <= top]


def skip_scan(run: RunInfo, df: pl.DataFrame, nblocks: int) -> pl.DataFrame:
    """M3 and chi3 recomputed with the first `skip` rows dropped."""
    q3 = df["q3"].to_numpy()
    sm, sc = scales(run)
    rows = []
    for s in skip_values(len(q3), nblocks):
        m, dm = block_stats(q3[s:], nblocks)
        c, dc = block_stats(q3[s:] ** 2, nblocks)
        rows.append({"skip": s, "M3": m * sm, "dM3": dm * sm, "chi3": c * sc, "dchi3": dc * sc})
    return pl.DataFrame(rows)
