"""
Revision M2 (EXPLORATORY): k inserted beats of five types in N1 linear RR windows, and per-window diagnostics of the
VWK test (PLAN.md Section 6).  The ventricular / random design equals Phase 11 explore11 for seeds 11000-11099.
"""
from __future__ import annotations

from itertools import combinations_with_replacement

import numpy as np

from analyses.revision_m1_m2 import vwkcfg as V
from experiments.phase10_final import harness10 as H

ENTROPY = 20261117
TYPES = ("ventricular", "atrial", "spike", "missed", "extra")
WIDTH = {"none": 2, "ventricular": 2, "atrial": 2, "spike": 1, "missed": 1, "extra": 2}
TYPE_CODE = {"none": 9, "ventricular": 0, "atrial": 1, "spike": 2, "missed": 3, "extra": 4}
PLACE_CODE = {"random": 0, "periodic": 1}
KS = (0, 1, 2, 4, 8, 16, 26, 51)
K_PERIODIC = (16, 26, 51)
SEEDS = range(11000, 11200)
SPACING = 3                     # Phase 5 ECTOPIC_MIN_SPACING
KSTAR = 6                       # nonlinear memory: affected rows = target or lags 1..6 on an affected interval
COOK_M = (1, 2, 4, 8, 16)
REP_SEEDS = (11000, 11001, 11002)
FS = 360.0


def base_window(seed):
    from experiments.phase11_robustness import run11 as R11
    x, _ = R11.synth_series("null", "N1_linear_rr", seed)
    return np.asarray(x, float)


def positions(n, k, seed, placement):
    """Positions and prematurity factors c.  random: explore11.insert_isolated draws (positions, then c in sorted
    order); periodic: spacing floor((n-2)/k), random offset."""
    if k == 0:
        return [], []
    if placement == "random":
        rng = np.random.default_rng(np.random.SeedSequence([20261103, 99, k, seed, 1]))
        chosen = []
        while len(chosen) < k:
            i = int(rng.integers(1, n - 1))
            if all(abs(i - j) >= SPACING for j in chosen):
                chosen.append(i)
        chosen.sort()
    else:
        rng = np.random.default_rng(np.random.SeedSequence([ENTROPY, 13, 1, k, seed]))
        s = (n - 2) // k
        off = int(rng.integers(1, s + 1))
        chosen = [off + j * s for j in range(k)]
    cs = [float(rng.uniform(0.6, 0.8)) for _ in chosen]
    return chosen, cs


def insert(x, k, seed, btype, placement):
    """Returns (series, affected-interval boolean array, positions).  Re-quantized to 1/360 s for k > 0."""
    x = np.array(x, dtype=float)
    n = len(x)
    pos, cs = positions(n, k, seed, placement)
    if k == 0:
        return x, np.zeros(n, bool), pos
    rt = np.random.default_rng(np.random.SeedSequence([ENTROPY, 60 + TYPE_CODE[btype], PLACE_CODE[placement], k, seed]))
    if btype == "ventricular":
        for i, c in zip(pos, cs):
            short = c * x[i]
            x[i + 1] += x[i] - short
            x[i] = short
        aff = np.zeros(n, bool)
        aff[pos] = True
        aff[np.asarray(pos) + 1] = True
        y = x
    elif btype == "atrial":
        for i, c in zip(pos, cs):
            x[i] = c * x[i]
            x[i + 1] = float(rt.uniform(1.0, 1.1)) * x[i + 1]
        aff = np.zeros(n, bool)
        aff[pos] = True
        aff[np.asarray(pos) + 1] = True
        y = x
    elif btype == "spike":
        for i, c in zip(pos, cs):
            s = 1.0 if rt.random() < 0.5 else -1.0
            x[i] = x[i] + s * (1.0 - c) * x[i]
        aff = np.zeros(n, bool)
        aff[pos] = True
        y = x
    elif btype in ("missed", "extra"):
        P = set(pos)
        out, a = [], []
        i = 0
        while i < n:
            if i in P and btype == "missed":
                out.append(x[i] + x[i + 1])
                a.append(True)
                i += 2
                continue
            if i in P and btype == "extra":
                u = float(rt.uniform(0.3, 0.7))
                out += [u * x[i], (1.0 - u) * x[i]]
                a += [True, True]
                i += 1
                continue
            out.append(x[i])
            a.append(False)
            i += 1
        y, aff = np.asarray(out, float), np.asarray(a, bool)
    else:
        raise ValueError(btype)
    return H.quantize_times(0.0, y, FS), aff, pos


# ------------------------------------------------------------------------------------------- diagnostics
def affected_rows(aff, idx, kstar=KSTAR):
    """Row (target index) n is affected iff any of intervals n-kstar..n is affected."""
    c = np.concatenate([[0], np.cumsum(aff.astype(int))])
    idx = np.asarray(idx)
    return (c[idx + 1] - c[np.maximum(idx - kstar, 0)]) > 0


def term_labels(kappa, d):
    lab = ["1"]
    for deg in range(1, d + 1):
        for c in combinations_with_replacement(range(1, kappa + 1), deg):
            lab.append("y" + "*y".join(str(j) for j in c))
    return lab


def term_lags(label):
    return [] if label == "1" else [int(t) for t in label[1:].split("*y")]


def term_kind(label):
    L = term_lags(label)
    if len(L) <= 1:
        return "constant" if not L else "linear"
    if len(set(L)) == 1:
        return "power"
    s = sorted(set(L))
    return "adjacent_product" if any(b - a == 1 for a, b in zip(s, s[1:])) else "product"


def leverage_cooks(X, e):
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    keep = s > s[0] * 1e-10 if len(s) else np.zeros(0, bool)
    p = int(keep.sum())
    h = np.sum(U[:, keep] ** 2, axis=1)
    n = len(e)
    s2 = float(np.sum(e ** 2) / max(n - p, 1))
    D = (e ** 2 / max(p * s2, 1e-300)) * h / np.maximum((1 - h) ** 2, 1e-12)
    return h, D


def sham_positions(n, aff, k, width, rng, gap=KSTAR + 2, tries=400):
    """k pseudo-positions whose would-be affected intervals [j, j + width) keep >= gap intervals from every truly
    affected interval and from each other; None if unavailable."""
    if k == 0:
        return []
    bad = np.zeros(n, bool)
    for i in np.nonzero(aff)[0]:
        bad[max(0, i - gap):min(n, i + gap + 1)] = True
    chosen = []
    for _ in range(tries):
        j = int(rng.integers(1, n - width))
        if bad[j:j + width].any():
            continue
        chosen.append(j)
        bad[max(0, j - gap):min(n, j + width + gap)] = True
        if len(chosen) == k:
            return sorted(chosen)
    return None


def diagnose(y, aff, cfg, rng_sham, rep=False, kstar=KSTAR, k=0, width=2):
    """Full M2 diagnostics of one window under configuration cfg."""
    o = V.vwk_cfg(y, cfg, internals=True)
    if not o.get("analysable"):
        return {"analysable": False}
    I = o.pop("_int")
    rows, el, en, t = I["rows"], I["e_lin"], I["e_nl"], I["t"]
    z = (y - y.mean()) / y.std()
    a6 = affected_rows(aff, rows, kstar)
    ak = affected_rows(aff, rows, max(o["kappa_nl"], 1))
    dl = el ** 2 - en ** 2
    tot = float(dl.sum())
    tm = t - t.mean()
    res = {k_: o[k_] for k_ in ("nonlinear", "p", "F", "r_lin", "r_nl", "kappa_lin", "kappa_nl", "d_nl", "eps_lin",
                                 "eps_nl", "C_lin", "C_nl", "n_rows", "d_opt_gt1", "test")}
    res["n"] = int(len(y))
    res["relimp"] = float(1 - I["rss_nl"] / I["rss_lin"])
    res["lag1"] = float(np.corrcoef(z[:-1], z[1:])[0, 1])
    res["lag2"] = float(np.corrcoef(z[:-2], z[2:])[0, 1])
    res["frac_aff6"] = float(a6.mean())
    res["frac_affk"] = float(ak.mean())
    res["delta_rss"] = tot
    res["share_aff6"] = float(dl[a6].sum() / tot) if tot != 0 else None
    res["share_affk"] = float(dl[ak].sum() / tot) if tot != 0 else None
    for g, msk in (("aff", a6), ("oth", ~a6)):
        ss = float(np.sum(tm[msk] ** 2))
        res[f"nrss_lin_{g}"] = float(np.sum(el[msk] ** 2) / ss) if msk.any() and ss > 0 else None
        res[f"nrss_nl_{g}"] = float(np.sum(en[msk] ** 2) / ss) if msk.any() and ss > 0 else None
    res["max_abs_z_aff"] = float(np.max(np.abs(z[aff]))) if aff.any() else None
    res["max_abs_z_oth"] = float(np.max(np.abs(z[~aff])))
    # model anatomy: nonlinear terms of the best nonlinear model, structural kind and ectopic concentration
    labs = term_labels(o["kappa_nl"], cfg.d_for(o["kappa_nl"]))[:o["r_nl"]]
    Xn = I["X_nl"]
    terms = []
    for j, lb in enumerate(labs):
        L = term_lags(lb)
        if len(L) < 2:
            continue
        col = Xn[:, j] ** 2
        hit = np.zeros(len(rows), bool)
        for lag in set(L):
            ix = rows - lag
            ok = ix >= 0
            hit[ok] |= aff[ix[ok]]
        terms.append([lb, float(col[hit].sum() / col.sum()) if col.sum() > 0 else 0.0])
    res["terms"] = terms
    # leverage and Cook's distance
    hn, Dn = leverage_cooks(Xn, en)
    hl, Dl = leverage_cooks(I["X_lin"], el)
    order = np.argsort(-Dn)
    res["lev_nl_aff_median"] = float(np.median(hn[a6])) if a6.any() else None
    res["lev_nl_oth_median"] = float(np.median(hn[~a6]))
    res["cook_nl_aff_max"] = float(Dn[a6].max()) if a6.any() else None
    res["cook_nl_oth_max"] = float(Dn[~a6].max())
    res["cook_lin_aff_max"] = float(Dl[a6].max()) if a6.any() else None
    res["cook_lin_oth_max"] = float(Dl[~a6].max())
    res["top8_cook_aff_share"] = float(a6[order[:8]].mean())
    # influence (a): remove the affected rows (all windows)
    n = len(y)
    full_aff = affected_rows(aff, np.arange(n), kstar)
    oa = V.vwk_cfg(y, cfg, valid_rows=~full_aff)
    res["rm_aff"] = {"nonlinear": bool(oa["nonlinear"]), "analysable": bool(oa.get("analysable")),
                     "n_removed": int(full_aff[rows].sum()), "n_rows": oa.get("n_rows")}
    # influence (b): sham removal of the same structure at pseudo-positions
    sp = sham_positions(n, aff, k, width, rng_sham)
    if sp is None or k == 0:
        res["rm_sham"] = {"available": bool(k == 0), "nonlinear": res["nonlinear"] if k == 0 else None}
    else:
        sa = np.zeros(n, bool)
        for j in sp:
            sa[j:j + width] = True
        sfull = affected_rows(sa, np.arange(n), kstar)
        os_ = V.vwk_cfg(y, cfg, valid_rows=~sfull)
        res["rm_sham"] = {"available": True, "nonlinear": bool(os_["nonlinear"]), "analysable": bool(os_.get("analysable")),
                          "n_removed": int(sfull[rows].sum())}
    # influence (c): top-m Cook's rows removed (raw-positive windows)
    if res["nonlinear"]:
        res["rm_cook"] = {}
        for m in COOK_M:
            vr = np.ones(n, bool)
            vr[rows[order[:m]]] = False
            oc = V.vwk_cfg(y, cfg, valid_rows=vr)
            res["rm_cook"][str(m)] = {"nonlinear": bool(oc["nonlinear"]), "analysable": bool(oc.get("analysable")),
                                      "aff_removed": int(a6[order[:m]].sum())}
    # k = 1: the R0 masked rule (rows free of ectopy over the configuration's longest lag)
    if k == 1:
        om = V.vwk_cfg(y, cfg, valid_rows=V.rows_mask(aff, cfg))
        res["rm_masked_rule"] = {"nonlinear": bool(om["nonlinear"]), "analysable": bool(om.get("analysable"))}
    if rep:
        cur = I["curves"]
        res["curves"] = {"lin": [float(v) for v in cur["lin"]],
                         "nl": [float(v) for v in cur[o["kappa_nl"]]], "kappa_nl": o["kappa_nl"]}
    return res
