"""
Phase 12 opt-in titration readings (FIDELITY.md section 4; PREREGISTRATION.md section 3).  Reuses the frozen Phase 10
implementation (experiments.phase10_final.titration10: design, leading_rss, _residuals, constants) and changes
nothing in it.  At the default reading R0, `vwk12` reproduces `titration10.vwk_test` and `titrate12` reproduces
`titration10.titrate` exactly (tests/test_titration12.py, 20 seeded series, raw and masked).

Reading options:
  f_stat     'df'    (R0) F = (RSS_lin/df_lin)/(RSS_nl/df_nl)
             'plain' (BP96 Fig. 1 legend) F_calc = eps_lin^2 / eps_nl^2 = RSS_lin / RSS_nl
  n_mode     'rows'  (R0) N in C(r) = log eps + r/N and in df = N - r is the number of regression rows
             'length' (BP96 Fig. 2 df, PM97 Fig. 2 legend) N = number of points in the series
  gate       False   (R0) always the F-test
             True    (BP96 p. 216) chi-square normality test (99 % cutoff) of both best models' residuals; if either
                     rejects, the one-sided Mann-Whitney test on squared residuals (alpha) replaces the F-test
  surrogate  None    (R0) | 'FT' | 'AAFT': R0 decision AND, for each of K = 10 surrogates of the analysed series,
                     the original's best nonlinear error is significantly smaller (the reading's F convention, alpha)
                     than min(best linear, best nonlinear) error of the surrogate (worst case over surrogates)
Decision: nonlinear iff C_nl < C_lin (d_opt > 1) and the test rejects at alpha (and the surrogate requirement holds).
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats

from experiments.phase10_final import titration10 as T

READINGS = {
    "R0": {},
    "R0-BP96F": {"f_stat": "plain", "n_mode": "length"},
    "R0-BP96": {"f_stat": "plain", "n_mode": "length", "gate": True},
    "R0-S-FT": {"surrogate": "FT"},
    "R0-S-AAFT": {"surrogate": "AAFT"},
}
K_SURR = 10
CHI2_ALPHA = 0.01


def reading(name):
    r = {"f_stat": "df", "n_mode": "rows", "gate": False, "surrogate": None}
    r.update(READINGS[name])
    return r


# ------------------------------------------------------------------------------------------------ surrogates
def ft_surrogate(y, rng):
    """Phase-randomised (FT) surrogate: same amplitude spectrum (hence autocorrelation and variance)."""
    y = np.asarray(y, float)
    n = len(y)
    Y = np.fft.rfft(y)
    ph = rng.uniform(0.0, 2 * np.pi, len(Y))
    ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    return np.fft.irfft(Y * np.exp(1j * ph), n)          # DC (and Nyquist, n even) kept; other phases uniform


def aaft_surrogate(y, rng):
    """Amplitude-adjusted FT surrogate (Theiler et al. 1992): same values as y, approximately the same spectrum."""
    y = np.asarray(y, float)
    n = len(y)
    ry = np.argsort(np.argsort(y))
    g = np.sort(rng.standard_normal(n))[ry]
    s = ft_surrogate(g, rng)
    return np.sort(y)[np.argsort(np.argsort(s))]


def chi2_normal_p(e):
    """Chi-square goodness of fit of standardised residuals to N(0,1): k = ceil(2 n^0.4) equiprobable bins,
    df = k - 3 (mean and SD estimated)."""
    e = np.asarray(e, float)
    n = len(e)
    sd = e.std()
    if n < 20 or sd <= 0:
        return 0.0
    z = (e - e.mean()) / sd
    k = int(math.ceil(2 * n ** 0.4))
    edges = stats.norm.ppf(np.arange(1, k) / k)
    obs = np.bincount(np.searchsorted(edges, z), minlength=k)
    exp = n / k
    stat = float(np.sum((obs - exp) ** 2) / exp)
    return float(stats.chi2.sf(stat, k - 3))


# ------------------------------------------------------------------------------------------------ fits
def fits12(y, rows, kappa_max, d_max, kappa_lin_max, n_pen):
    """titration10._fits with the penalty denominator N as an argument (identical code otherwise)."""
    t = y[rows]
    sst = float(np.sum((t - t.mean()) ** 2))
    N = n_pen
    Xl = T.design(y, kappa_lin_max, 1, rows)
    rss_l, _ = T.leading_rss(Xl, t)
    r = np.arange(1, Xl.shape[1] + 1)
    C_l = T.EPS_POWER * np.log(np.maximum(rss_l / sst, 1e-300)) + r / N
    il = int(np.argmin(C_l))
    best_lin = {"C": float(C_l[il]), "r": int(il + 1), "kappa": int(il), "rss": float(rss_l[il])}
    best_nl = None
    for kappa in range(1, kappa_max + 1):
        X = T.design(y, kappa, d_max, rows)
        rss, _ = T.leading_rss(X, t)
        rr = np.arange(1, X.shape[1] + 1)
        C = T.EPS_POWER * np.log(np.maximum(rss / sst, 1e-300)) + rr / N
        lo = kappa + 1
        if lo >= X.shape[1]:
            continue
        j = lo + int(np.argmin(C[lo:]))
        cand = {"C": float(C[j]), "r": int(j + 1), "kappa": kappa,
                "d": int(next(dd for dd in range(2, d_max + 1) if j + 1 <= T.n_terms(kappa, dd))),
                "rss": float(rss[j])}
        if best_nl is None or cand["C"] < best_nl["C"]:
            best_nl = cand
    return best_lin, best_nl, len(rows), sst


def _std(y):
    y = np.asarray(y, float)
    s = np.std(y)
    return (y - y.mean()) / (s if s > 0 else 1.0)


def _resid(y, rows, lin, nl, d_max):
    rl = T._residuals(y, rows, lin["r"], max(lin["kappa"], 1), 1) if lin["r"] > 1 else y[rows] - y[rows].mean()
    rn = T._residuals(y, rows, nl["r"], nl["kappa"], d_max)
    return rl, rn


def vwk12(y, valid_rows=None, name="R0", rng=None, surr_source=None, kappa_max=T.KAPPA_MAX, d_max=T.D_MAX,
          alpha=T.ALPHA, kappa_lin_max=None, details=False, K=K_SURR):
    """Barahona-Poon nonlinearity test under reading `name`.  `rng` is needed only for surrogate readings;
    `surr_source` is the series the surrogates are generated from (default: y itself; masked arm: the interpolated
    series), fitted on the same rows as the original."""
    R = reading(name)
    kl = T.KAPPA_LIN_MAX if kappa_lin_max is None else kappa_lin_max
    y = _std(y)
    n = len(y)
    rows = np.arange(max(kappa_max, kl), n)
    if valid_rows is not None:
        rows = rows[np.asarray(valid_rows, bool)[rows]]
    if len(rows) < max(T.MIN_ROWS if valid_rows is not None else 0, 3 * T.n_terms(kappa_max, d_max)):
        return {"nonlinear": False, "analysable": False, "n_rows": int(len(rows))}
    n_pen = len(rows) if R["n_mode"] == "rows" else n
    lin, nl, Nrows, sst = fits12(y, rows, kappa_max, d_max, kl, n_pen)
    out = {"analysable": True, "n_rows": int(Nrows), "C_lin": lin["C"], "C_nl": nl["C"], "r_lin": lin["r"],
           "r_nl": nl["r"], "kappa_lin": lin["kappa"], "kappa_nl": nl["kappa"], "d_nl": nl["d"],
           "eps_lin": float(np.sqrt(lin["rss"] / sst)), "eps_nl": float(np.sqrt(nl["rss"] / sst))}
    better = nl["C"] < lin["C"]
    df1, df2 = n_pen - lin["r"], n_pen - nl["r"]
    F = (lin["rss"] / df1) / max(nl["rss"] / df2, 1e-300) if R["f_stat"] == "df" else lin["rss"] / max(nl["rss"], 1e-300)
    test = "F"
    rl = rn = None
    if R["gate"] or details:
        rl, rn = _resid(y, rows, lin, nl, d_max)
        out["chi2_p_lin"], out["chi2_p_nl"] = chi2_normal_p(rl), chi2_normal_p(rn)
        if R["gate"] and min(out["chi2_p_lin"], out["chi2_p_nl"]) < CHI2_ALPHA:
            test = "MWW"
    if test == "F":
        p = float(stats.f.sf(F, df1, df2))
    else:
        p = float(stats.mannwhitneyu(rl ** 2, rn ** 2, alternative="greater").pvalue)
    out.update(F=float(F), df1=int(df1), df2=int(df2), test=test, p=p)
    out["nonlinear"] = bool(better and p < alpha)
    if details:
        out["F_calc_plain"] = lin["rss"] / max(nl["rss"], 1e-300)
        out["F_crit"] = float(stats.f.isf(alpha, df1, df2))
        u = stats.mannwhitneyu(rn ** 2, rl ** 2, alternative="less")
        n1 = n2 = len(rn)
        out["MW_Z"] = float((u.statistic - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12))
        out["d_opt_gt1"] = bool(better)
    if R["surrogate"] and out["nonlinear"]:
        src = _std(y if surr_source is None else surr_source)
        gen = ft_surrogate if R["surrogate"] == "FT" else aaft_surrogate
        e2o = nl["rss"] / sst
        worst = 0.0
        for _ in range(K):
            s = _std(gen(src, rng))
            ls, ns, _, sst_s = fits12(s, rows, kappa_max, d_max, kl, n_pen)
            best = ls if ls["rss"] <= ns["rss"] else ns
            e2s = best["rss"] / sst_s
            d_s, d_o = n_pen - best["r"], n_pen - nl["r"]
            Fs = (e2s / d_s) / max(e2o / d_o, 1e-300) if R["f_stat"] == "df" else e2s / max(e2o, 1e-300)
            ps = float(stats.f.sf(Fs, d_s, d_o))
            worst = max(worst, ps)
        out["surr_p_worst"] = worst
        out["nonlinear"] = bool(worst < alpha)
    return out


def titrate12(y, rng, valid_rows=None, name="R0", n_real=T.N_REAL, kappa_max=T.KAPPA_MAX, d_max=T.D_MAX,
              kappa_lin_max=None):
    """titration10.titrate with the reading's test inside the bisection (surrogate readings: not titrated;
    see `nl_surrogate_reading`)."""
    if reading(name)["surrogate"]:
        raise ValueError("surrogate readings are not titrated (PREREGISTRATION.md 3.3)")
    y = np.asarray(y, float)
    sy = float(np.std(y))
    kw = dict(valid_rows=valid_rows, name=name, kappa_max=kappa_max, d_max=d_max, kappa_lin_max=kappa_lin_max)
    base = vwk12(y, **kw)
    out = {"base": base, "analysable": base.get("analysable", False)}
    if not base.get("nonlinear"):
        out.update(NL=0.0, NL_i=[0.0] * n_real, positive=False, NL_power=0.0)
        return out
    nls = []
    for _ in range(n_real):
        xi = rng.standard_normal(len(y))
        xi = (xi - xi.mean()) / xi.std()
        det = lambda a: vwk12(y + a * sy * xi, **kw)["nonlinear"]  # noqa: E731
        lo, hi = 0.0, T.HI0
        while det(hi) and hi < T.HI_MAX:
            lo, hi = hi, 2 * hi
        if lo >= T.HI_MAX:
            nls.append(T.HI_MAX)
            continue
        for _ in range(T.N_BISECT):
            mid = 0.5 * (lo + hi)
            if det(mid):
                lo = mid
            else:
                hi = mid
        nls.append(lo)
    NL = float(np.mean(nls))
    out.update(NL=NL, NL_i=[float(v) for v in nls], positive=bool(NL > 0),
               NL_power=float(np.mean(np.square(nls))))
    return out


def rows_masked(mask_intervals):
    return T.rows_from_mask(mask_intervals)
