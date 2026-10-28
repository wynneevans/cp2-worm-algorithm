import numpy as np
import polars as pl

from cp2analysis.observables import chunk_scan, consistency, scales, skip_scan, skip_values, summarise
from cp2analysis.runs import RunInfo

RUN = RunInfo(bc="closed", beta=0.05, lx=2, ly=2, lt=4, mu3=0.0, mu8=0.0)


def test_scales():
    m, chi = scales(RUN)
    assert m == 1 / (2 * 2)
    assert chi == 0.05 / (4 * 4)


def test_summarise_of_constant_charge():
    df = pl.DataFrame({"q3": [1.0] * 1000, "q8": [-1.0] * 1000})
    s = summarise(RUN, df, 10)
    assert s["n_used"] == 1000
    assert s["M3"] == 1 / 4 and s["dM3"] == 0
    assert abs(s["chi3"] - 0.05 / 16) < 1e-15


def test_skip_scan_starts_at_zero_and_stays_in_first_half():
    values = skip_values(100_000, 50)
    assert values[0] == 0 and max(values) <= 50_000
    df = pl.DataFrame({"q3": [1.0, -1.0] * 50_000})
    scan = skip_scan(RUN, df, 50)
    assert scan["skip"].to_list() == values


def test_chunk_scan_covers_the_run_without_overlap():
    df = pl.DataFrame({"q3": [1.0] * 1000})
    chunks = chunk_scan(RUN, df, 4, 10)
    assert chunks["start"].to_list() == [0, 250, 500, 750] and chunks["end"].to_list()[-1] == 1000
    assert chunks["chi3"].to_list() == [0.05 / 16] * 4 and set(chunks["dchi3"]) == {0.0}
    assert consistency(chunks["chi3"].to_numpy(), chunks["dchi3"].to_numpy()) is None  # no errors to compare


def test_consistency_is_near_one_for_independent_data_and_large_for_a_drift():
    rng = np.random.default_rng(0)
    q = rng.choice([-1.0, 1.0], size=80_000) * rng.integers(1, 4, size=80_000)
    chunks = chunk_scan(RUN, pl.DataFrame({"q3": q}), 8, 20)
    assert 0.2 < consistency(chunks["chi3"].to_numpy(), chunks["dchi3"].to_numpy()) < 3
    drift = q * np.linspace(1, 2, len(q))
    chunks = chunk_scan(RUN, pl.DataFrame({"q3": drift}), 8, 20)
    assert consistency(chunks["chi3"].to_numpy(), chunks["dchi3"].to_numpy()) > 10
