"""
Phase 11 A2: heart rate turbulence (HRT) measurement and an opt-in HRT model for synthetic ectopy.
These are DESCRIPTIVE measurements and generator options; no detector parameter depends on them.

Definitions (METHODS.md Phase 11 Section 2 lists the sources and what could not be verified):
  Tachogram of a usable isolated ventricular premature beat (VPB) at beat v (symbol 'V'):
      RR_-5 .. RR_-1   the 5 sinus intervals before the coupling interval (intervals v-6 .. v-2)
      CI               coupling interval (v-1 -> v),  CP compensatory pause (v -> v+1)
      RR_1 .. RR_20    the 20 sinus intervals after the compensatory pause (intervals v+1 .. v+20)
  Usable (filter rules of the Schmidt group, as restated by Cygankiewicz et al. 2003 and the HRT
  consensus, Bauer et al. 2008): beats v-6 .. v-1 and v+1 .. v+21 normal (label 'N'); every sinus
  interval in [300, 2000] ms; consecutive sinus intervals differ by <= 200 ms; every sinus interval within
  20 % of the reference (mean of RR_-5 .. RR_-1); CI <= 0.8 * reference (prematurity >= 20 %);
  CP >= 1.2 * reference.
  TO (per VPB, %) = 100 * [(RR_1 + RR_2) - (RR_-2 + RR_-1)] / (RR_-2 + RR_-1); subject TO = mean over VPBs.
  TS (ms per RR interval) = maximum regression slope over any 5 consecutive intervals within RR_1 .. RR_15
  of the subject's AVERAGED tachogram (consensus); TS20 (secondary) within RR_1 .. RR_20 (Schmidt 1999).
  Normalised tachogram of a VPB: h_j = RR_j / mean(RR_-2, RR_-1) - 1, j = 1..20.

HRT model (opt-in; A2 generators): after the compensatory pause the next 20 sinus intervals are multiplied
by (1 + h_j).  Profiles: 'empirical' (a group's subject-weighted mean normalised tachogram), 'normal'
(parametric two-phase profile calibrated so that TO = -2.5 % and TS = 13 ms/RR on a constant baseline equal
to the pre-VPB interval), 'blunted' (TO = 0 %, TS = 2.5 ms/RR).
"""
from __future__ import annotations

import numpy as np

N_PRE, N_POST = 5, 20
RR_LO_MS, RR_HI_MS = 300.0, 2000.0
MAX_STEP_MS = 200.0
MAX_DEV = 0.20
PREMATURITY = 0.80
COMP_MIN = 1.20
TS_RANGE = 15
TARGETS = {"normal": (-2.5, 13.0), "blunted": (0.0, 2.5)}   # (TO %, TS ms/RR): Grimm et al. 2003 medians /
#                                                             Schmidt et al. 1999 abnormal cut-offs (as given)


# ----------------------------------------------------------------------------------- measurement
def max_slope(tach, n_range=TS_RANGE, width=5):
    """Maximum least-squares slope over `width` consecutive values among the first n_range of tach."""
    tach = np.asarray(tach, float)
    k = np.arange(width, dtype=float)
    kc = k - k.mean()
    best = -np.inf
    for s in range(0, n_range - width + 1):
        y = tach[s:s + width]
        best = max(best, float(np.sum(kc * (y - y.mean())) / np.sum(kc * kc)))
    return best


def turbulence_onset(pre2, pre1, post1, post2):
    return 100.0 * ((post1 + post2) - (pre2 + pre1)) / (pre2 + pre1)


def tachograms(t, sym, lab):
    """Usable VPB tachograms of one recording.  t: beat times (s); sym: symbols; lab: 'N'/'X' labels.
    Returns a list of dicts with the intervals in ms and the per-VPB TO."""
    rr = np.diff(np.asarray(t, float)) * 1000.0
    sym = np.asarray(sym).astype(str)
    lab = np.asarray(lab).astype(str)
    nb = len(t)
    out = []
    for v in np.nonzero(sym == "V")[0]:
        if v - 6 < 0 or v + 21 >= nb:
            continue
        if not (np.all(lab[v - 6:v] == "N") and np.all(lab[v + 1:v + 22] == "N")):
            continue
        pre = rr[v - 6:v - 1]                       # RR_-5 .. RR_-1
        ci, cp = rr[v - 1], rr[v]
        post = rr[v + 1:v + 21]                     # RR_1 .. RR_20
        sinus = np.concatenate([pre, post])
        ref = float(pre.mean())
        if np.any(sinus < RR_LO_MS) or np.any(sinus > RR_HI_MS):
            continue
        if np.any(np.abs(np.diff(pre)) > MAX_STEP_MS) or np.any(np.abs(np.diff(post)) > MAX_STEP_MS):
            continue
        if np.any(np.abs(sinus - ref) > MAX_DEV * ref):
            continue
        if ci > PREMATURITY * ref or cp < COMP_MIN * ref:
            continue
        base = 0.5 * (pre[-2] + pre[-1])
        out.append({"beat": int(v), "pre": pre.tolist(), "ci": float(ci), "cp": float(cp), "post": post.tolist(),
                    "ref": ref, "TO": turbulence_onset(pre[-2], pre[-1], post[0], post[1]),
                    "h": (post / base - 1.0).tolist(), "ci_ratio": float(ci / base), "cp_ratio": float(cp / base)})
    return out


def subject_hrt(tg):
    """Subject TO (mean of per-VPB TO), TS and TS20 (on the averaged tachogram), mean normalised profile."""
    if not tg:
        return {"n_vpb": 0}
    post = np.mean([g["post"] for g in tg], axis=0)
    return {"n_vpb": len(tg), "TO": float(np.mean([g["TO"] for g in tg])), "TS": max_slope(post, TS_RANGE),
            "TS20": max_slope(post, N_POST), "h": np.mean([g["h"] for g in tg], axis=0).tolist(),
            "mean_tachogram_ms": post.tolist(), "ref_ms": float(np.mean([g["ref"] for g in tg])),
            "ci_ratio": float(np.mean([g["ci_ratio"] for g in tg])),
            "cp_ratio": float(np.mean([g["cp_ratio"] for g in tg]))}


# ----------------------------------------------------------------------------------- profiles
def _shapes(n=N_POST):
    j = np.arange(1, n + 1, dtype=float)
    s = np.exp(-(j - 1.0) / 1.5)                    # early shortening (acceleration), decays within ~5 beats
    d = np.exp(-((j - 7.0) / 4.0) ** 2)             # deceleration peaking at beat 7, back to ~0 by beat 20
    return s, d


def profile_metrics(h, rr0_s):
    """TO (%) and TS (ms/RR) of the profile applied to a constant baseline rr0 (pre-VPB interval)."""
    h = np.asarray(h, float)
    tach = 1000.0 * rr0_s * (1.0 + h)
    return 100.0 * (h[0] + h[1]) / 2.0, max_slope(tach, TS_RANGE)


def parametric_profile(kind, rr0_s):
    """h = -a s + b d with (a, b) >= 0 solving TO, TS = TARGETS[kind] at baseline rr0_s."""
    from scipy.optimize import least_squares
    s, d = _shapes()
    to_t, ts_t = TARGETS[kind]

    def resid(p):
        h = -p[0] * s + p[1] * d
        to, ts = profile_metrics(h, rr0_s)
        return [to - to_t, (ts - ts_t) / 5.0]
    best = None
    for a0, b0 in ((0.02, 0.05), (0.05, 0.02), (0.01, 0.01), (0.03, 0.10)):
        r = least_squares(resid, [a0, b0], bounds=([0, 0], [1, 1]), xtol=1e-12, ftol=1e-12)
        if best is None or r.cost < best.cost:
            best = r
    return (-best.x[0] * s + best.x[1] * d).tolist()


# ----------------------------------------------------------------------------------- insertion
def insert_event(x, j, variant, draws, h_emp=None):
    """One ectopic event in RR series x (seconds) at premature interval index j (ectopic beat j+1).
    variant: V0 none | V1 ventricular | V2 V1 + HRT-empirical | V3 V1 + HRT-normal | V4 V1 + HRT-blunted |
             V5 atrial-type | V6 couplet | V7 HRT-empirical without a premature beat.
    draws: dict with c, c2 (coupling factors U(0.6, 0.8)) and d (atrial post-ectopic factor U(1.0, 1.1)),
    shared by every variant of a window (paired design).
    Ventricular (Phase 6 model): x_j -> c x_j; compensatory x_{j+1} -> x_j + x_{j+1} - c x_j (sinus node not
    reset).  Couplet: x_j -> c x_j, x_{j+1} -> c2 x_j, x_{j+2} -> x_j + x_{j+1} + x_{j+2} - (c + c2) x_j.
    Atrial (Phase 6 E5): x_j -> c x_j, x_{j+1} -> d x_{j+1}.  HRT multiplies the 20 intervals after the
    compensatory interval by (1 + h); its pre-VPB baseline is rr0 = mean(x_{j-2}, x_{j-1})."""
    x = np.asarray(x, float).copy()
    if variant == "V0":
        return x
    rr0 = 0.5 * (x[j - 2] + x[j - 1])
    if variant in ("V1", "V2", "V3", "V4"):
        b0, b1 = x[j], x[j + 1]
        x[j] = draws["c"] * b0
        x[j + 1] = b0 + b1 - x[j]
        k = j + 2
    elif variant == "V5":
        x[j] = draws["c"] * x[j]
        x[j + 1] = draws["d"] * x[j + 1]
        k = None
    elif variant == "V6":
        b0, b1, b2 = x[j], x[j + 1], x[j + 2]
        x[j], x[j + 1] = draws["c"] * b0, draws["c2"] * b0
        x[j + 2] = b0 + b1 + b2 - x[j] - x[j + 1]
        k = None
    elif variant == "V7":
        k = j + 2
    else:
        raise ValueError(variant)
    prof = {"V2": "empirical", "V3": "normal", "V4": "blunted", "V7": "empirical"}.get(variant)
    if prof is not None:
        h = np.asarray(h_emp if prof == "empirical" else parametric_profile(prof, rr0), float)
        x[k:k + N_POST] = x[k:k + N_POST] * (1.0 + h)
    return x


VARIANTS = ("V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7")
EDGE = 30                     # event (premature .. last HRT interval) >= 30 intervals from both ends
SPAN = 2 + N_POST + 1         # intervals j .. j + 22 may change (couplet: j .. j + 2; HRT: j + 2 .. j + 21)


def event_position(n, rng):
    return int(rng.integers(EDGE, n - EDGE - SPAN + 1))


def draws_for(rng):
    return {"c": float(rng.uniform(0.6, 0.8)), "c2": float(rng.uniform(0.6, 0.8)), "d": float(rng.uniform(1.0, 1.1))}
