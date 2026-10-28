import numpy as np

from cp2analysis.stats import acf, block_stats, mser, tau_int


def test_block_stats_recovers_mean_and_error():
    x = np.random.default_rng(1).normal(3.0, 2.0, 100_000)
    mean, err = block_stats(x, 50)
    expected = 2.0 / np.sqrt(len(x))
    assert abs(mean - 3.0) < 4 * err
    assert abs(err - expected) < 0.3 * expected


def test_tau_int_of_ar1():
    rng = np.random.default_rng(2)
    n, phi = 400_000, 0.9  # tau_int = (1 + phi) / (2 (1 - phi)) = 9.5
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = 0
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    t, err, _ = tau_int(acf(x), n=n)
    assert abs(t - 9.5) < 4 * err + 0.3


def test_white_noise_and_constant_series():
    rng = np.random.default_rng(3)
    white = rng.normal(size=100_000)
    assert abs(tau_int(acf(white), n=len(white))[0] - 0.5) < 0.05
    assert acf(np.ones(10)) is None
    assert mser(np.ones(100)) is None


def test_mser_stationary_and_transient():
    rng = np.random.default_rng(3)
    n = 100_000
    assert mser(rng.normal(size=n)) < 0.1 * n
    transient = rng.normal(size=n) + 20 * np.exp(-np.arange(n) / 500)
    assert 500 <= mser(transient) <= 8000
