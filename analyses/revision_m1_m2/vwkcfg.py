"""
Revision M1/M2 (EXPLORATORY): the Barahona-Poon VWK nonlinearity test with explicit switches (PLAN.md Section 3.2).

Reuses the frozen Phase 10 / 12 primitives (titration10.design, leading_rss, constants; titration12 surrogates and the
chi-square test) and changes nothing in them.  With every switch at its default, `vwk_cfg` reproduces
`titration10.vwk_test` and `titrate_cfg` reproduces `titration10.titrate` exactly; with the Phase 12 switch sets it
reproduces `titration12.vwk12` exactly (tests/test_regression.py).

Switches (Cfg fields; defaults = R0):
  n_mode        'rows' | 'length'      S1: N in C(r) and df = series length
  f_stat        'df' | 'plain'         S2: F = RSS_lin / RSS_nl
  gate          False | True           S3: chi-square (99 %) gate -> one-sided Mann-Whitney on squared residuals
  test          'F' | 'MW'             S3b: Mann-Whitney always
  surrogate     None | 'FT' | 'AAFT'   S4: BP96 confirmatory surrogate step (10 surrogates, worst case)
  kappa_max, d_max, m_cap              S5a: d_max 8 with m_cap 84; S5b: kappa_max 7 (M 120)
  kappa_lin_max None -> (largest nonlinear model size) - 1   (matched linear search)
  rows          'common' | 'zeropad'   S6: every point a target, pre-sample lags = 0
  decision      'both' | 'test'        S7: positive iff the test rejects (no C_nl < C_lin requirement)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy import stats

from experiments.phase10_final import titration10 as T
from experiments.phase12_original import titration12 as T12


@dataclass(frozen=True)
class Cfg:
    n_mode: str = "rows"
    f_stat: str = "df"
    gate: bool = False
    test: str = "F"
    surrogate: str | None = None
    kappa_max: int = T.KAPPA_MAX
    d_max: int = T.D_MAX
    m_cap: int | None = None
    kappa_lin_max: int | None = None
    rows: str = "common"
    decision: str = "both"
    alpha: float = T.ALPHA
    K_surr: int = T12.K_SURR
    switches: tuple = field(default=(), compare=False)

    # ---------------------------------------------------------------- derived
    def d_for(self, kappa):
        """Degree used for embedding kappa: d_max, or the largest d <= d_max with M <= m_cap (0 if d < 2)."""
        if self.m_cap is None:
            return self.d_max
        d = 0
        for dd in range(1, self.d_max + 1):
            if T.n_terms(kappa, dd) <= self.m_cap:
                d = dd
        return d if d >= 2 else 0

    @property
    def m_max(self):
        return max(T.n_terms(k, self.d_for(k)) for k in range(1, self.kappa_max + 1) if self.d_for(k) >= 2)

    @property
    def kl(self):
        return self.m_max - 1 if self.kappa_lin_max is None else self.kappa_lin_max

    @property
    def K(self):
        """Longest lag of any candidate model (masked-row span; common first target)."""
        return max(self.kappa_max, self.kl)


R0 = Cfg()


def _std(y):
    y = np.asarray(y, float)
    s = np.std(y)
    return (y - y.mean()) / (s if s > 0 else 1.0)


def _prep(y, cfg, valid_rows):
    """Standardised series (padded for 'zeropad'), target rows (indices into the returned array), offset."""
    y = _std(y)
    n = len(y)
    K = cfg.K
    if cfg.rows == "common":
        rows = np.arange(K, n)
        if valid_rows is not None:
            rows = rows[np.asarray(valid_rows, bool)[rows]]
        return y, rows, 0
    if cfg.rows == "zeropad":
        yy = np.concatenate([np.zeros(K), y])
        rows = np.arange(K, K + n)
        if valid_rows is not None:
            rows = rows[np.asarray(valid_rows, bool)[rows - K]]
        return yy, rows, K
    raise ValueError(cfg.rows)


def _fits(yy, rows, cfg, n_pen, curves=False):
    """titration10._fits / titration12.fits12 generalised to per-kappa degrees (identical arithmetic at defaults)."""
    t = yy[rows]
    sst = float(np.sum((t - t.mean()) ** 2))
    N = n_pen
    kl = cfg.kl
    Xl = T.design(yy, kl, 1, rows)
    rss_l, _ = T.leading_rss(Xl, t)
    r = np.arange(1, Xl.shape[1] + 1)
    C_l = T.EPS_POWER * np.log(np.maximum(rss_l / sst, 1e-300)) + r / N
    il = int(np.argmin(C_l))
    best_lin = {"C": float(C_l[il]), "r": int(il + 1), "kappa": int(il), "rss": float(rss_l[il])}
    best_nl = None
    cur = {"lin": C_l} if curves else None
    for kappa in range(1, cfg.kappa_max + 1):
        dk = cfg.d_for(kappa)
        if dk < 2:
            continue
        X = T.design(yy, kappa, dk, rows)
        rss, _ = T.leading_rss(X, t)
        rr = np.arange(1, X.shape[1] + 1)
        C = T.EPS_POWER * np.log(np.maximum(rss / sst, 1e-300)) + rr / N
        if curves:
            cur[kappa] = C
        lo = kappa + 1
        if lo >= X.shape[1]:
            continue
        j = lo + int(np.argmin(C[lo:]))
        cand = {"C": float(C[j]), "r": int(j + 1), "kappa": kappa,
                "d": int(next(dd for dd in range(2, dk + 1) if j + 1 <= T.n_terms(kappa, dd))),
                "dk": dk, "rss": float(rss[j])}
        if best_nl is None or cand["C"] < best_nl["C"]:
            best_nl = cand
    return best_lin, best_nl, sst, cur


def _resid(yy, rows, lin, nl):
    rl = T._residuals(yy, rows, lin["r"], max(lin["kappa"], 1), 1) if lin["r"] > 1 else yy[rows] - yy[rows].mean()
    rn = T._residuals(yy, rows, nl["r"], nl["kappa"], nl["dk"])
    return rl, rn


def min_rows(cfg, masked):
    return max(T.MIN_ROWS if masked else 0, 3 * cfg.m_max)


def vwk_cfg(y, cfg=R0, valid_rows=None, rng=None, surr_source=None, details=False, internals=False):
    """Barahona-Poon test under configuration `cfg`.  Returns the titration10 keys (+ titration12 keys).
    internals=True adds '_int': rows (original indices), residual vectors, design matrices of both best models and the
    C(r) curves (for the M2 diagnostics; not JSON-serialisable)."""
    yy, rows, off = _prep(y, cfg, valid_rows)
    n = len(yy) - off
    if len(rows) < min_rows(cfg, valid_rows is not None):
        return {"nonlinear": False, "analysable": False, "n_rows": int(len(rows))}
    n_pen = len(rows) if cfg.n_mode == "rows" else n
    lin, nl, sst, cur = _fits(yy, rows, cfg, n_pen, curves=internals)
    out = {"analysable": True, "n_rows": int(len(rows)), "C_lin": lin["C"], "C_nl": nl["C"], "r_lin": lin["r"],
           "r_nl": nl["r"], "kappa_lin": lin["kappa"], "kappa_nl": nl["kappa"], "d_nl": nl["d"],
           "eps_lin": float(np.sqrt(lin["rss"] / sst)), "eps_nl": float(np.sqrt(nl["rss"] / sst))}
    better = nl["C"] < lin["C"]
    df1, df2 = n_pen - lin["r"], n_pen - nl["r"]
    if cfg.f_stat == "df":
        F = (lin["rss"] / df1) / max(nl["rss"] / df2, 1e-300)
    else:
        F = lin["rss"] / max(nl["rss"], 1e-300)
    test = "MW" if cfg.test == "MW" else "F"
    rl = rn = None
    if cfg.gate or details or internals or cfg.test == "MW":
        rl, rn = _resid(yy, rows, lin, nl)
        if cfg.gate or details:
            out["chi2_p_lin"], out["chi2_p_nl"] = T12.chi2_normal_p(rl), T12.chi2_normal_p(rn)
            if cfg.gate and cfg.test == "F" and min(out["chi2_p_lin"], out["chi2_p_nl"]) < T12.CHI2_ALPHA:
                test = "MWW"
    if test == "F":
        p = float(stats.f.sf(F, df1, df2))
    else:
        p = float(stats.mannwhitneyu(rl ** 2, rn ** 2, alternative="greater").pvalue)
    out.update(F=float(F), df1=int(df1), df2=int(df2), test=test, p=p)
    out["d_opt_gt1"] = bool(better)
    out["nonlinear"] = bool((better or cfg.decision == "test") and p < cfg.alpha)
    if details:
        out["F_calc_plain"] = lin["rss"] / max(nl["rss"], 1e-300)
        out["F_crit"] = float(stats.f.isf(cfg.alpha, df1, df2))
        u = stats.mannwhitneyu(rn ** 2, rl ** 2, alternative="less")
        n1 = n2 = len(rn)
        out["MW_Z"] = float((u.statistic - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12))
    if cfg.surrogate and out["nonlinear"]:
        src = T12._std(yy[off:] if surr_source is None else surr_source)
        gen = T12.ft_surrogate if cfg.surrogate == "FT" else T12.aaft_surrogate
        e2o = nl["rss"] / sst
        worst = 0.0
        for _ in range(cfg.K_surr):
            s = T12._std(gen(src, rng))
            ss = np.concatenate([np.zeros(off), s]) if off else s
            ls, ns, sst_s, _ = _fits(ss, rows, cfg, n_pen)
            best = ls if ls["rss"] <= ns["rss"] else ns
            e2s = best["rss"] / sst_s
            d_s, d_o = n_pen - best["r"], n_pen - nl["r"]
            Fs = (e2s / d_s) / max(e2o / d_o, 1e-300) if cfg.f_stat == "df" else e2s / max(e2o, 1e-300)
            ps = float(stats.f.sf(Fs, d_s, d_o))
            worst = max(worst, ps)
        out["surr_p_worst"] = worst
        out["nonlinear"] = bool(worst < cfg.alpha)
    if internals:
        Xn = T.design(yy, nl["kappa"], nl["dk"], rows)[:, :nl["r"]]
        Xl = T.design(yy, max(lin["kappa"], 1), 1, rows)[:, :lin["r"]] if lin["r"] > 1 else np.ones((len(rows), 1))
        out["_int"] = {"rows": rows - off, "e_lin": rl, "e_nl": rn, "X_nl": Xn, "X_lin": Xl, "sst": sst,
                       "rss_lin": lin["rss"], "rss_nl": nl["rss"], "curves": cur, "t": yy[rows]}
    return out


def titrate_cfg(y, rng, cfg=R0, valid_rows=None, n_real=T.N_REAL):
    """titration10.titrate with the configuration's test inside the bisection.  Surrogate configurations are not
    titrated (PLAN.md 3.2: NL taken from R0 when the surrogate requirement holds noise-free)."""
    if cfg.surrogate:
        raise ValueError("surrogate configurations are not titrated")
    y = np.asarray(y, float)
    sy = float(np.std(y))
    base = vwk_cfg(y, cfg, valid_rows=valid_rows)
    out = {"base": base, "analysable": base.get("analysable", False)}
    if not base.get("nonlinear"):
        out.update(NL=0.0, NL_i=[0.0] * n_real, positive=False, NL_power=0.0)
        return out
    nls = []
    for _ in range(n_real):
        xi = rng.standard_normal(len(y))
        xi = (xi - xi.mean()) / xi.std()
        det = lambda a: vwk_cfg(y + a * sy * xi, cfg, valid_rows=valid_rows)["nonlinear"]  # noqa: E731
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
    out.update(NL=NL, NL_i=[float(v) for v in nls], positive=bool(NL > 0), NL_power=float(np.mean(np.square(nls))))
    return out


def rows_mask(mask_intervals, cfg=R0):
    """Masked-arm rows: target and every lag up to cfg.K free of ectopy-related intervals (R0: titration10 rule)."""
    m = np.asarray(mask_intervals, bool)
    if cfg.rows == "common":
        return T.rows_from_mask(m, kappa_max=cfg.K)
    n = len(m)
    c = np.concatenate([[0], np.cumsum(m.astype(int))])
    idx = np.arange(n)
    return (c[idx + 1] - c[np.maximum(idx - cfg.K, 0)]) == 0


def with_switches(base, **kw):
    return replace(base, **kw)
