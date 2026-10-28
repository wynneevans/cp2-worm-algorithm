import polars as pl

from cp2analysis.observables import scales, skip_scan, skip_values, summarise
from cp2analysis.runs import RunInfo

RUN = RunInfo(bc="closed", beta=0.05, lx=2, ly=2, lt=4, mu3=0.0, mu8=0.0)


def test_scales():
    m, chi = scales(RUN)
    assert m == 1 / (2 * 2)
    assert chi == 0.05 / (4 * 4)


def test_summarise_of_constant_charge():
    df = pl.DataFrame({"q3": [1.0] * 1000, "q8": [-1.0] * 1000})
    s = summarise(RUN, df, 10)
    assert s["lattice"] == "2x2"
    assert s["M3"] == 1 / 4 and s["dM3"] == 0
    assert abs(s["chi3"] - 0.05 / 16) < 1e-15


def test_skip_scan_starts_at_zero_and_stays_in_first_half():
    values = skip_values(100_000, 50)
    assert values[0] == 0 and max(values) <= 50_000
    df = pl.DataFrame({"q3": [1.0, -1.0] * 50_000})
    scan = skip_scan(RUN, df, 50)
    assert scan["skip"].to_list() == values
