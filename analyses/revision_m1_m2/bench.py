"""
Revision M1/M2 (EXPLORATORY): benchmark generators and per-unit evaluation helpers (PLAN.md Sections 3-4).
Every design reuses the Phase 10-12 generators and seeding, so that R0 reproduces the stored results.
"""
from __future__ import annotations

import os
import zlib
from pathlib import Path

import numpy as np

from analyses.revision_m1_m2 import vwkcfg as V
from experiments.phase10_final import data10 as D
from experiments.phase10_final import harness10 as H
from experiments.phase10_final import verify_titration as VT
from experiments.phase12_original import run12 as R12
from experiments.phase12_original import systems12 as S

ENTROPY = 20261117
ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "results" / "revision_m1_m2" / "_cache" / "series"
T1_SYSTEMS = ("logistic", "henon", "ikeda", "ecological", "gopy_sna", "rossler", "duffing", "lorenz",
              "mackey_glass_100")
T1_SECONDARY = ("mackey_glass_17",)
FLOWS = ("rossler", "duffing", "lorenz", "mackey_glass_100", "mackey_glass_17")
P1_M = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0)
NOISE_ECO = tuple(range(0, 65, 5))
B2_HI, B2_STEPS = R12.B2_HI, R12.B2_STEPS
PART = {"t1": 1, "eco": 2, "log_x0": 3, "log_noise": 4, "pm": 5, "surr": 7, "win": 8, "seg": 9, "tr": 10, "synth": 11,
        "pm97": 12, "m2": 13}


def crc(s):
    return zlib.crc32(str(s).encode())


def ssq(part, *key):
    return np.random.SeedSequence([ENTROPY, PART[part], *[int(k) for k in key]])


def surr_rng(cfg, *key):
    """RNG for the surrogate step (only used by surrogate configurations)."""
    return np.random.default_rng(ssq("surr", *key)) if cfg.surrogate else None


# ------------------------------------------------------------------------------------------- Table 1 series
def _flow_series(system, seed, m):
    """systems12 flows with the sampling step multiplied by m (P1); m = 1 reproduces systems12 exactly."""
    n = S.N_BENCH
    if system == "rossler":
        a, b, c = 0.2, 0.2, 5.7
        f = lambda t, s: np.array([-s[1] - s[2], s[0] + a * s[1], b + s[2] * (s[0] - c)])  # noqa: E731
        return S._retry(lambda rng: S._flow(f, rng.uniform(-1, 1, 3) + [1, 1, 0], n, 0.5 * m, 0.01, 500.0), 6, seed)[0]
    if system == "duffing":
        delta, gamma, omega = 0.25, 0.3, 1.0
        T0 = 2 * np.pi / omega
        f = lambda t, s: np.array([s[1], -delta * s[1] + s[0] - s[0] ** 3 + gamma * np.cos(omega * t)])  # noqa: E731
        return S._retry(lambda rng: S._flow(f, rng.uniform(-1, 1, 2), n, T0 / 12 * m, T0 / 600, 100 * T0), 7, seed)[0]
    if system == "lorenz":
        sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0
        f = lambda t, s: np.array([sigma * (s[1] - s[0]), s[0] * (rho - s[2]) - s[1], s[0] * s[1] - beta * s[2]])  # noqa: E731
        return S._retry(lambda rng: S._flow(f, rng.uniform([-10, -10, 10], [10, 10, 40]), n, 0.075 * m, 0.0025, 100.0),
                        8, seed)[0]
    if system.startswith("mackey_glass"):
        tau = 100.0 if system.endswith("100") else 17.0
        return np.asarray(VT.mackey_glass(tau, n=n, ts=12.5 * m, seed=int(seed)), float)
    raise ValueError(system)


def t1_series(system, seed, m=1.0, cache=True):
    if system not in FLOWS:
        return S.SYSTEMS[system](seed)
    p = CACHE / f"{system}_{seed}_{m:g}.npy"
    if cache and p.exists():
        return np.load(p)
    y = S.SYSTEMS[system](seed) if m == 1.0 else _flow_series(system, seed, m)
    if cache:
        CACHE.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(f".{os.getpid()}.tmp.npy")
        np.save(tmp, y)
        os.replace(tmp, p)
    return y


def t1_noise(system, y, seed):
    """Phase 12 B2 noise shape (one FT surrogate of the clean series per seed)."""
    return S.correlated_noise(y, np.random.default_rng(R12.ss("b2", R12.crc(system), seed, 1)))


def t1_max_noise(system, seed, cfg, m=1.0):
    """Largest correlated-noise % still detected (Phase 12 B2 statistic; bisection on [0, 200] %, 8 halvings)."""
    y = t1_series(system, seed, m)
    xi = t1_noise(system, y, seed)
    key = (crc(system), seed, int(round(1000 * m)))

    def det(p):
        return V.vwk_cfg(y + p / 100.0 * y.std() * xi, cfg, rng=surr_rng(cfg, *key, int(round(100 * p))))["nonlinear"]
    base = V.vwk_cfg(y, cfg, details=True, rng=surr_rng(cfg, *key, 0))
    out = {"base": {k: base.get(k) for k in ("nonlinear", "kappa_nl", "d_nl", "r_lin", "r_nl", "kappa_lin",
                                             "F_calc_plain", "p", "eps_lin", "eps_nl", "test", "d_opt_gt1")}}
    if not base["nonlinear"]:
        out.update(max_noise=0.0, censored=False)
        return out
    if det(B2_HI):
        out.update(max_noise=B2_HI, censored=True)
        return out
    lo, hi = 0.0, B2_HI
    for _ in range(B2_STEPS):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if det(mid) else (lo, mid)
    out.update(max_noise=lo, censored=False)
    return out


def p1_score(system, seed, m):
    """P1 selection score: eps_lin - eps_nl, noise-free, R0 (BP96 p. 216)."""
    y = t1_series(system, seed, m)
    o = V.vwk_cfg(y, V.R0)
    return {"eps_lin": o["eps_lin"], "eps_nl": o["eps_nl"], "diff": o["eps_lin"] - o["eps_nl"], "nonlinear": o["nonlinear"]}


# ------------------------------------------------------------------------------------------- ecological (Fig. 1d)
def eco_levels(seed, cfg):
    """Phase 12 B1 design: correlated noise p = 0..60 % (5 % steps), one noise shape per seed."""
    y = S.ecological(seed)
    xi = S.correlated_noise(y, np.random.default_rng(R12.ss("b1", seed, 1)))
    out = {}
    for p in NOISE_ECO:
        o = V.vwk_cfg(y + p / 100.0 * y.std() * xi, cfg, details=True, rng=surr_rng(cfg, 1, seed, p))
        out[str(p)] = {k: o.get(k) for k in ("nonlinear", "d_opt_gt1", "F", "F_calc_plain", "F_crit", "p", "test",
                                             "r_lin", "r_nl", "kappa_nl", "d_nl", "MW_Z", "C_lin", "C_nl")}
    neg = [p for p in NOISE_ECO if not out[str(p)]["nonlinear"]]
    return {"levels": out, "loss": float(neg[0]) if neg else 65.0, "censored": not neg}


# ------------------------------------------------------------------------------------------- logistic NL (PB01)
def log_x0(seed):
    return float(np.random.default_rng(ssq("log_x0", seed)).uniform(0.1, 0.9))


def log_nl(seed, cfg, r=3.7):
    """NL of the logistic map at r from initial condition x0(seed) (1,000 transient iterations, 1,000 points)."""
    y = VT.logistic(r, x0=log_x0(seed))
    if cfg.surrogate:
        o0 = V.titrate_cfg(y, np.random.default_rng(ssq("log_noise", seed)), V.R0)
        g = V.vwk_cfg(y, cfg, rng=surr_rng(cfg, 3, seed))
        pos = bool(o0["positive"] and g["nonlinear"])
        return {"NL": o0["NL"] if pos else 0.0, "positive": pos, "R0_NL": o0["NL"], "x0": log_x0(seed)}
    o = V.titrate_cfg(y, np.random.default_rng(ssq("log_noise", seed)), cfg)
    return {"NL": o["NL"], "NL_i": o["NL_i"], "positive": o["positive"], "x0": log_x0(seed),
            "kappa_nl": o["base"].get("kappa_nl"), "d_nl": o["base"].get("d_nl"), "r_nl": o["base"].get("r_nl")}


def scan_nl(r, cfg):
    """Manuscript Fig. 1 / Phase 10 V4 series and noise seeds."""
    y = VT.logistic(r)
    rng = np.random.default_rng(np.random.SeedSequence([20261003, zlib.crc32(f"logistic|{r}".encode())]))
    if cfg.surrogate:
        o0 = V.titrate_cfg(y, rng, V.R0)
        g = V.vwk_cfg(y, cfg, rng=surr_rng(cfg, 4, crc(r)))
        pos = bool(o0["positive"] and g["nonlinear"])
        return {"NL": o0["NL"] if pos else 0.0, "positive": pos}
    o = V.titrate_cfg(y, rng, cfg)
    return {"NL": o["NL"], "NL_i": o["NL_i"], "positive": o["positive"]}


# ------------------------------------------------------------------------------------------- arms (Phase 10 rules)
def arm(rr, mask, arm_name, cfg, rng=None, keep=("nonlinear", "analysable", "n_rows", "p", "r_lin", "r_nl",
                                                     "kappa_lin", "kappa_nl", "d_nl", "test")):
    """Noise-free decision of one arm: raw | masked (rows free of ectopy over the configuration's longest lag) |
    edited (Phase 6 edit, NN fraction >= 0.80) | wu (masked intervals deleted, >= 400 left)."""
    x = np.asarray(rr, float)
    m = None if mask is None else np.asarray(mask, bool)
    if arm_name == "raw":
        o = V.vwk_cfg(x, cfg, rng=rng)
    elif arm_name == "masked":
        xe, _ = H.edit(x, imask=m)
        o = V.vwk_cfg(x, cfg, valid_rows=V.rows_mask(m, cfg), rng=rng, surr_source=xe if xe is not None else x)
    elif arm_name == "edited":
        xe, nnf = H.edit(x, imask=m)
        if xe is None or nnf < H.NN_EDIT_MIN:
            return {"eligible": False, "nonlinear": False, "nn_fraction": float(nnf)}
        o = V.vwk_cfg(xe, cfg, rng=rng)
    elif arm_name == "wu":
        nn = H.wu_nn(x, m)
        if len(nn) < H.SEG_MIN_N:
            return {"eligible": False, "nonlinear": False, "n": int(len(nn))}
        o = V.vwk_cfg(nn, cfg, rng=rng)
    else:
        raise ValueError(arm_name)
    out = {k: o.get(k) for k in keep if k in o}
    out["eligible"] = True
    if "surr_p_worst" in o:
        out["surr_p_worst"] = o["surr_p_worst"]
    return out


# ------------------------------------------------------------------------------------------- PM97 (Phase 12 protocol)
_PMCTX = {}


def pm_ctx(pilot=False):
    if pilot not in _PMCTX:
        _PMCTX[pilot] = R12.Ctx(pilot)
    return _PMCTX[pilot]


def pm_subject(db, rec, cfg, arms=("raw",), lengths=(500, 2000), pilot=False):
    """Phase 12 PM97 block and segments; decisions of the requested arms (raw, rawf, masked, edited, wu)."""
    ctx = pm_ctx(pilot)
    blk = R12.pm97_block(ctx, db, rec)
    if "excluded" in blk:
        return {"excluded": blk["excluded"]}
    key = crc(f"{db}|{rec}")
    out = {"group": "NSR" if db == "nsrdb" else "CHF", "rule": blk["rule"]}
    for L in lengths:
        res = []
        for sg in R12.pm97_segments(blk, L, ctx):
            r = {"k": sg["k"], "ectopy": sg["ectopy"], "n_out": sg["n_out_of_range"]}
            for a in arms:
                rng = surr_rng(cfg, key, L, sg["k"], crc(a))
                if a == "raw":
                    r[a] = arm(sg["raw"], None, "raw", cfg, rng)
                elif a == "rawf":
                    r[a] = arm(sg["rawf"], None, "raw", cfg, rng)
                else:
                    r[a] = arm(sg["rawf"], sg["m"], a, cfg, rng)
            res.append(r)
        out[f"seg{L}"] = res
    return out


# ------------------------------------------------------------------------------------------- confirmatory units
def subject_of(db, rec):
    return f"{db}:{rec}"


def records(db, root):
    return D.records(db, root)
