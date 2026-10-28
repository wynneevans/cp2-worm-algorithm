"""Statistics for Monte Carlo time series (numpy only, no I/O)."""

import numpy as np

SOKAL_C = 6.0  # automatic-window constant
MSER_BATCH = 5


def block_stats(x, nblocks: int) -> tuple[float, float]:
    """Mean and standard error from a non-overlapping block jackknife."""
    x = np.asarray(x, dtype=float)
    size = len(x) // nblocks
    blocks = x[: size * nblocks].reshape(nblocks, size).mean(axis=1)
    mean = blocks.mean()
    jack = (blocks.sum() - blocks) / (nblocks - 1)
    err = np.sqrt((nblocks - 1) / nblocks * ((jack - jack.mean()) ** 2).sum())
    return mean, err


def acf(x) -> np.ndarray | None:
    """Normalised autocorrelation function rho(0..N-1) via FFT, or None for a constant series."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    if not np.any(x):
        return None
    f = np.fft.rfft(x, 2 * n)
    c = np.fft.irfft(f * np.conj(f))[:n]
    return c / c[0]


def tau_int(rho, c: float = SOKAL_C, n: int | None = None) -> tuple[float, float, int]:
    """Integrated autocorrelation time, its error and the window, using Sokal's automatic windowing.

    tau_int = 1/2 + sum_{t=1}^{W} rho(t), with the smallest W >= c * tau_int(W);
    error = tau * sqrt((4W + 2) / N) (Madras-Sokal).
    """
    n = n or len(rho)
    tau = np.concatenate([[0.5], 0.5 + np.cumsum(rho[1:])])
    w = np.arange(len(tau))
    ok = np.nonzero(w >= c * tau)[0]
    win = int(ok[0]) if len(ok) else len(tau) - 1
    t = tau[win]
    return t, t * np.sqrt((4 * win + 2) / n), win


def mser(x, batch: int = MSER_BATCH) -> int | None:
    """MSER-`batch` truncation point in measurements (searched in the first half), None if constant."""
    x = np.asarray(x, dtype=float)
    if not np.any(x - x[0]):
        return None
    m = len(x) // batch
    y = x[: m * batch].reshape(m, batch).mean(axis=1)
    k = np.arange(m, 0, -1)  # points kept when dropping d = 0..m-1
    s1 = np.cumsum(y[::-1])[::-1]  # suffix sums
    s2 = np.cumsum((y**2)[::-1])[::-1]
    stat = (s2 - s1**2 / k) / k**2
    return int(np.argmin(stat[: m // 2])) * batch
