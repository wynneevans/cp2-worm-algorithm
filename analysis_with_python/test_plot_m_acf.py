"""Run with: uv run python test_plot_m_acf.py"""
import numpy as np

from plot_m_acf import acf, tau_int

rng = np.random.default_rng(2)
n, phi = 400_000, 0.9          # AR(1): tau_int = 1/2 + phi/(1-phi) -> (1+phi)/(2(1-phi)) = 9.5
x = np.empty(n)
x[0] = 0
e = rng.normal(size=n)
for i in range(1, n):
    x[i] = phi * x[i - 1] + e[i]
t, err, w = tau_int(acf(x), n=n)
assert abs(t - 9.5) < 4 * err + 0.3, (t, err)
white = acf(rng.normal(size=n))
assert abs(tau_int(white, n=n)[0] - 0.5) < 0.05
assert acf(np.ones(10)) is None
print("ok", t, err)
