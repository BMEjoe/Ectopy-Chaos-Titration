"""titration11 at the default reading reproduces titration10.titrate exactly (20 seeded series)."""
import numpy as np
import pytest

from experiments.phase10_final import titration10 as T
from experiments.phase11_robustness import titration11 as T11


def _series(i):
    rng = np.random.default_rng(11000 + i)
    n = 400 + 10 * i
    kind = i % 4
    if kind == 0:                                   # logistic map + noise (positive)
        x = rng.uniform(0.2, 0.8)
        y = np.empty(n)
        for k in range(500 + n):
            x = 3.9 * x * (1 - x)
            if k >= 500:
                y[k - 500] = x
        y = y + 0.05 * rng.standard_normal(n)
    elif kind == 1:                                 # AR(2) (negative)
        e = rng.standard_normal(n + 200)
        y = np.zeros(n + 200)
        for k in range(2, n + 200):
            y[k] = 1.5 * y[k - 1] - 0.75 * y[k - 2] + e[k]
        y = y[200:]
    elif kind == 2:                                 # RR-like with couplets (ectopy-driven positive)
        y = 0.8 + 0.03 * rng.standard_normal(n)
        for j in range(30, n - 3, 60):
            y[j], y[j + 1] = 0.5, 0.5
            y[j + 2] = 1.4
    else:                                           # Henon x + noise
        a, b, u, v = 1.4, 0.3, 0.1, 0.1
        y = np.empty(n)
        for k in range(500 + n):
            u, v = 1 - a * u * u + v, b * u
            if k >= 500:
                y[k - 500] = u
        y = y + 0.05 * rng.standard_normal(n)
    mask = np.zeros(n, bool)
    mask[rng.integers(0, n, 2)] = True
    return y, mask


@pytest.mark.parametrize("i", range(20))
def test_default_reading_reproduces_titration10(i):
    y, mask = _series(i)
    for vr in (None, T.rows_from_mask(mask)):
        a = T.titrate(y, np.random.default_rng(i), valid_rows=vr)
        b = T11.titrate11(y, np.random.default_rng(i), valid_rows=vr)
        c = T11.titrate_reading(y, np.random.default_rng(i), "R0", mask=None if vr is None else mask)
        for o in (b, c):
            assert o["NL_i"] == a["NL_i"]
            assert o["positive"] == a["positive"]
            assert o["base"] == a["base"]
