"""Run with: uv run python test_analyse_m.py"""
import numpy as np

from analyse_m import block_stats, parse_name

p = parse_name("M_wormCP2AFMmu2dplaq_BC_y-open_nMeas_1000000_beta_0.05_XYT_2_3_4_mu3_0.000_mu8_0.000_nSim_1.dat")
assert p == dict(bc="y-open", beta=0.05, lx=2, ly=3, lt=4, mu3=0.0, mu8=0.0), p

rng = np.random.default_rng(1)
x = rng.normal(3.0, 2.0, 100_000)
mean, err = block_stats(x, 50)
assert abs(mean - 3.0) < 4 * err, (mean, err)
expected = 2.0 / np.sqrt(len(x))
assert abs(err - expected) < 0.3 * expected, (err, expected)
print("ok")
