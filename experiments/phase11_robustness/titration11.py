"""
Phase 11 titration VARIANTS (A5 readings, A2 fixed-seed titration).  Reuses the frozen Phase 10
implementation (experiments.phase10_final.titration10: design, leading_rss, rows_from_mask, constants)
and changes nothing in it.  At the default reading (R0) `titrate11` reproduces `titration10.titrate`
exactly (tests/test_titration11.py, 20 seeded series, raw and masked).

Reading options (PREREGISTRATION.md A5):
  kappa_lin_max   linear memory limit (R0 83; R1 6; R2 20; R3 50)
  test            'F_var' (R0), 'MWW' (R4, Mann-Whitney on squared residuals), 'F_nested' (R5)
  eps_power       C(r) = eps_power * log(RSS/SST) + r/N: 0.5 = eps as RMS (R0); 1.0 = eps as variance (R6)
  scheme          'bisection' (R0) or 'increment' (R7, noise added in steps of 1 % of sigma_y until the
                  first non-detection; titration10.titrate scheme='increment')
The masked-arm row rule uses the longest lag of any candidate model, max(KAPPA_MAX, kappa_lin_max)
(R0: 83, as titration10.rows_from_mask's default).
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from experiments.phase10_final import titration10 as T

READINGS = {
    "R0": {},
    "R1": {"kappa_lin_max": 6},
    "R2": {"kappa_lin_max": 20},
    "R3": {"kappa_lin_max": 50},
    "R4": {"test": "MWW"},
    "R5": {"test": "F_nested"},
    "R6": {"eps_power": 1.0},
    "R7": {"scheme": "increment"},
}


def reading(name):
    r = {"kappa_lin_max": T.KAPPA_LIN_MAX, "test": "F_var", "eps_power": T.EPS_POWER, "scheme": "bisection"}
    r.update(READINGS[name])
    return r


def _fits11(y, rows, kappa_max, d_max, kappa_lin_max, eps_power):
    """titration10._fits with EPS_POWER as an argument (identical code otherwise)."""
    t = y[rows]
    sst = float(np.sum((t - t.mean()) ** 2))
    N = len(rows)
    Xl = T.design(y, kappa_lin_max, 1, rows)
    rss_l, _ = T.leading_rss(Xl, t)
    r = np.arange(1, Xl.shape[1] + 1)
    C_l = eps_power * np.log(np.maximum(rss_l / sst, 1e-300)) + r / N
    il = int(np.argmin(C_l))
    best_lin = {"C": float(C_l[il]), "r": int(il + 1), "kappa": int(il), "rss": float(rss_l[il])}
    best_nl = None
    for kappa in range(1, kappa_max + 1):
        X = T.design(y, kappa, d_max, rows)
        rss, _ = T.leading_rss(X, t)
        rr = np.arange(1, X.shape[1] + 1)
        C = eps_power * np.log(np.maximum(rss / sst, 1e-300)) + rr / N
        lo = kappa + 1
        if lo >= X.shape[1]:
            continue
        j = lo + int(np.argmin(C[lo:]))
        cand = {"C": float(C[j]), "r": int(j + 1), "kappa": kappa,
                "d": int(next(dd for dd in range(2, d_max + 1) if j + 1 <= T.n_terms(kappa, dd))),
                "rss": float(rss[j])}
        if best_nl is None or cand["C"] < best_nl["C"]:
            best_nl = cand
    return best_lin, best_nl, N, sst


def vwk_test11(y, valid_rows=None, kappa_max=T.KAPPA_MAX, d_max=T.D_MAX, alpha=T.ALPHA, test="F_var",
               kappa_lin_max=None, eps_power=T.EPS_POWER):
    """titration10.vwk_test with eps_power as an argument (identical code otherwise)."""
    kl = T.KAPPA_LIN_MAX if kappa_lin_max is None else kappa_lin_max
    y = np.asarray(y, float)
    s = np.std(y)
    y = (y - y.mean()) / (s if s > 0 else 1.0)
    n = len(y)
    rows = np.arange(max(kappa_max, kl), n)
    if valid_rows is not None:
        rows = rows[np.asarray(valid_rows, bool)[rows]]
    if len(rows) < max(T.MIN_ROWS if valid_rows is not None else 0, 3 * T.n_terms(kappa_max, d_max)):
        return {"nonlinear": False, "analysable": False, "n_rows": int(len(rows))}
    lin, nl, N, sst = _fits11(y, rows, kappa_max, d_max, kl, eps_power)
    out = {"analysable": True, "n_rows": int(N), "C_lin": lin["C"], "C_nl": nl["C"], "r_lin": lin["r"],
           "r_nl": nl["r"], "kappa_lin": lin["kappa"], "kappa_nl": nl["kappa"], "d_nl": nl["d"],
           "eps_lin": float(np.sqrt(lin["rss"] / sst)), "eps_nl": float(np.sqrt(nl["rss"] / sst))}
    better = nl["C"] < lin["C"]
    if test == "F_var":
        df1, df2 = N - lin["r"], N - nl["r"]
        F = (lin["rss"] / df1) / max(nl["rss"] / df2, 1e-300)
        p = float(stats.f.sf(F, df1, df2))
    elif test == "F_nested":
        dr = nl["r"] - lin["r"]
        if dr <= 0:
            p = 0.0 if nl["rss"] < lin["rss"] else 1.0
        else:
            F = ((lin["rss"] - nl["rss"]) / dr) / max(nl["rss"] / (N - nl["r"]), 1e-300)
            p = float(stats.f.sf(F, dr, N - nl["r"]))
    elif test == "MWW":
        rl = T._residuals(y, rows, lin["r"], max(lin["kappa"], 1), 1) if lin["r"] > 1 else y[rows] - y[rows].mean()
        rn = T._residuals(y, rows, nl["r"], nl["kappa"], d_max)
        p = float(stats.mannwhitneyu(rl ** 2, rn ** 2, alternative="greater").pvalue)
    else:
        raise ValueError(test)
    out["p"] = p
    out["nonlinear"] = bool(better and p < alpha)
    return out


def titrate11(y, rng, valid_rows=None, n_real=T.N_REAL, test="F_var", kappa_max=T.KAPPA_MAX, d_max=T.D_MAX,
              scheme="bisection", step=0.01, kappa_lin_max=None, eps_power=T.EPS_POWER):
    """titration10.titrate with eps_power as an argument (identical code otherwise)."""
    y = np.asarray(y, float)
    sy = float(np.std(y))
    kw = dict(valid_rows=valid_rows, test=test, kappa_max=kappa_max, d_max=d_max, kappa_lin_max=kappa_lin_max,
              eps_power=eps_power)
    base = vwk_test11(y, **kw)
    out = {"base": base, "analysable": base.get("analysable", False)}
    if not base.get("nonlinear"):
        out.update(NL=0.0, NL_i=[0.0] * n_real, positive=False, NL_power=0.0)
        return out
    nls = []
    for _ in range(n_real):
        xi = rng.standard_normal(len(y))
        xi = (xi - xi.mean()) / xi.std()
        det = lambda a: vwk_test11(y + a * sy * xi, **kw)["nonlinear"]  # noqa: E731
        if scheme == "bisection":
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
        else:
            a = 0.0
            while a + step <= T.HI_MAX and det(a + step):
                a += step
            nls.append(a)
    NL = float(np.mean(nls))
    out.update(NL=NL, NL_i=[float(v) for v in nls], positive=bool(NL > 0),
               NL_power=float(np.mean(np.square(nls))))
    return out


def rows_from_mask11(mask_intervals, kappa_lin_max=None):
    kl = T.KAPPA_LIN_MAX if kappa_lin_max is None else kappa_lin_max
    return T.rows_from_mask(mask_intervals, kappa_max=max(T.KAPPA_MAX, kl))


def titrate_reading(y, rng, name="R0", mask=None):
    """Titration under reading `name`; mask (interval mask) -> masked arm with the reading's row rule."""
    r = reading(name)
    vr = None if mask is None else rows_from_mask11(mask, r["kappa_lin_max"])
    return titrate11(y, rng, valid_rows=vr, test=r["test"], scheme=r["scheme"], kappa_lin_max=r["kappa_lin_max"],
                     eps_power=r["eps_power"])
