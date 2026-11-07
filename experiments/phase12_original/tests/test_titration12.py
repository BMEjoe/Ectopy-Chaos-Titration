"""titration12 at the default reading R0 reproduces titration10.vwk_test and titration10.titrate exactly
(20 seeded series, raw and masked); surrogate generators preserve what they should."""
import numpy as np
import pytest

from experiments.phase10_final import titration10 as T
from experiments.phase12_original import titration12 as T12


def _series(i):
    rng = np.random.default_rng(12000 + i)
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


KEYS = ("analysable", "n_rows", "C_lin", "C_nl", "r_lin", "r_nl", "kappa_lin", "kappa_nl", "d_nl", "eps_lin",
        "eps_nl", "p", "nonlinear")


@pytest.mark.parametrize("i", range(20))
def test_r0_reproduces_titration10(i):
    y, mask = _series(i)
    for vr in (None, T.rows_from_mask(mask)):
        a = T.vwk_test(y, valid_rows=vr)
        b = T12.vwk12(y, valid_rows=vr, name="R0")
        for k in KEYS:
            if k in a:
                assert a[k] == b[k], k
        ta = T.titrate(y, np.random.default_rng(i), valid_rows=vr)
        tb = T12.titrate12(y, np.random.default_rng(i), valid_rows=vr, name="R0")
        assert ta["NL_i"] == tb["NL_i"]
        assert ta["positive"] == tb["positive"]
        assert ta["NL"] == tb["NL"]


def test_surrogates_preserve_spectrum_and_values():
    rng = np.random.default_rng(12099)
    y = np.cumsum(rng.standard_normal(512)) * 0.1 + rng.standard_normal(512)
    s = T12.ft_surrogate(y, np.random.default_rng(1))
    assert np.allclose(np.abs(np.fft.rfft(s)), np.abs(np.fft.rfft(y)))
    assert not np.allclose(s, y)
    a = T12.aaft_surrogate(y, np.random.default_rng(2))
    assert np.allclose(np.sort(a), np.sort(y))


def test_chi2_gate_behaviour():
    rng = np.random.default_rng(12098)
    assert T12.chi2_normal_p(rng.standard_normal(900)) > 0.01
    assert T12.chi2_normal_p(rng.exponential(size=900)) < 0.01
