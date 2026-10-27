"""Run with: uv run python test_thermalisation.py"""
import numpy as np

from thermalisation import mser

rng = np.random.default_rng(3)
n = 100_000
stationary = rng.normal(size=n)
assert mser(stationary) < 0.1 * n
assert mser(np.ones(100)) is None

# exponential transient of time scale 500 starting 20 sigma away
transient = rng.normal(size=n) + 20 * np.exp(-np.arange(n) / 500)
d = mser(transient)
assert 500 <= d <= 8000, d
print("ok", d)
