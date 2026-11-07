"""
Phase 11 shared harness: sham masks (A1), extended masks (A2.5), transplant events (A2.3-4), window /
segment covariates (A3), K3 with a degeneracy diagnostic (B2), variance-preserving spike-ins (B3).
Everything reuses the frozen Phase 10 harness (experiments.phase10_final.harness10) and changes nothing
in it.  RNG streams: SeedSequence([ENTROPY11, part code, ...task identity]).
"""
from __future__ import annotations

import zlib

import numpy as np

from experiments.phase10_final import harness10 as H

ENTROPY11 = 20261103           # Phase 11 RNG entropy (distinct from every earlier phase's)
N_SHAM_WIN, N_SHAM_SEG = 5, 2
SHAM_TRIES = 200
HRT_GUARD = 20                 # intervals after an abnormal beat that a sham may not touch (HRT)
PART_CODE = {"sham_win": 1, "sham_seg": 2, "transplant": 3, "synth1": 4, "donor_win": 5, "donor_seg": 6,
             "b3": 7, "b4": 8}


def ss(part, *key):
    return np.random.SeedSequence([ENTROPY11, PART_CODE[part], *[int(k) for k in key]])


def key_of(*parts):
    """Stable integer key of strings / ints (CRC32)."""
    return zlib.crc32("|".join(str(p) for p in parts).encode())


# ---------------------------------------------------------------------------------------- masks
def abnormal_beats(labels):
    return np.nonzero(np.asarray(labels).astype(str) != "N")[0]


def forbidden(mask, labels, guard=HRT_GUARD):
    """True ectopy-related intervals plus the `guard` intervals after every abnormal beat."""
    f = np.asarray(mask, bool).copy()
    n = len(f)
    for b in abnormal_beats(labels):
        f[min(b, n):min(b + guard + 1, n)] = True        # intervals b .. b + guard (interval b starts at beat b)
    return f


def sham_masks(mask, fb, n_sham, rng, tries=SHAM_TRIES):
    """Circular shifts of the true interval mask that touch no forbidden interval (fb = forbidden(...));
    distinct offsets.  Returns a list of (offset, sham_mask); fewer than n_sham entries if a sham is
    unavailable in `tries` draws ("sham unavailable")."""
    m = np.asarray(mask, bool)
    fb = np.asarray(fb, bool)
    n = len(m)
    if not m.any():
        return []
    out, used = [], set()
    for _ in range(n_sham):
        for _ in range(tries):
            o = int(rng.integers(1, n))
            if o in used:
                continue
            s = np.roll(m, o)
            if not np.any(s & fb):
                out.append((o, s))
                used.add(o)
                break
        else:
            break
    return out


def borrow_mask(donor_mask, n):
    """A donor's interval mask applied at the same positions to a series of n intervals (truncated or padded
    with unmasked intervals)."""
    d = np.asarray(donor_mask, bool)
    out = np.zeros(n, bool)
    k = min(n, len(d))
    out[:k] = d[:k]
    return out


def extended_mask(mask, labels, guard=HRT_GUARD):
    """A2.5: the K3 mask plus the `guard` intervals after every abnormal beat."""
    return forbidden(mask, labels, guard)


# ---------------------------------------------------------------------------------------- covariates
def nn_covariates(rr, mask, t_start):
    """Per window / segment, from normal-to-normal intervals only (intervals not in the K3 mask):
    mean NN, SDNN, RMSSD (successive differences of adjacent NN intervals only), the half-mean shift
    |mean(2nd half) - mean(1st half)| / SD (halves by interval index), the K1 detrend ratio
    (final_pipeline.linear_trend_ratio on the NN series), elapsed hours at the start."""
    import final_pipeline as fp
    rr = np.asarray(rr, float)
    m = np.asarray(mask, bool)
    nn = rr[~m]
    n = len(rr)
    adj = (~m[:-1]) & (~m[1:])
    d = np.diff(rr)[adj]
    half = np.arange(n) < n // 2
    a, b = rr[(~m) & half], rr[(~m) & ~half]
    sd = float(np.std(nn)) if len(nn) > 1 else float("nan")
    return {"mean_nn": float(np.mean(nn)) if len(nn) else float("nan"), "sdnn_nn": sd,
            "rmssd_nn": float(np.sqrt(np.mean(d ** 2))) if len(d) else float("nan"),
            "half_shift": float(abs(b.mean() - a.mean()) / sd) if len(a) and len(b) and sd > 0 else float("nan"),
            "trend_ratio": float(fp.linear_trend_ratio(nn)) if len(nn) > 2 else float("nan"),
            "elapsed_h": float(t_start / 3600.0), "n_nn": int(len(nn))}


# ---------------------------------------------------------------------------------------- K3 diagnostic
def k3_diag(rr, labels, config=None):
    """final_pipeline.masked_growth_chaos_test re-run with the surrogate prediction-error mean and SD per m
    exposed (identical computation; the decision, z_m and G equal the frozen function's: test in
    tests/test_harness11.py).  Degenerate (B2) iff the SD at the m attaining z_max is < 1e-6."""
    import final_pipeline as fp
    cfg = fp.CFG if config is None else config
    x = fp._as_1d(rr)
    lab = np.asarray(labels).astype(str)
    bad = lab != "N"
    mask = bad[:-1] | bad[1:]
    xi = x.copy()
    if mask.any() and (~mask).sum() >= 2:
        idx = np.arange(len(x))
        xi[mask] = np.interp(idx[mask], idx[~mask], x[~mask])
    h = zlib.crc32(np.ascontiguousarray(np.asarray(x, dtype=float)).tobytes())
    rng = np.random.default_rng(np.random.SeedSequence([cfg.masked_growth_seed_entropy, int(h), 1]))
    ms = (2, 3, 4, 5)
    mv = cfg.masked_growth_min_valid
    obs = {m: fp._masked_nlp_errors(x, mask, m, (1, 2, 3, 4, 5), min_valid=mv) for m in ms}
    out = {"n_masked": int(mask.sum())}
    if obs[2] is None:
        out.update({"analysable": False, "detected": False})
        return out
    sur = {m: [] for m in ms}
    for _ in range(cfg.masked_growth_surrogates):
        s = fp.iaaft_surrogate(xi, rng)
        for m in ms:
            e = fp._masked_nlp_errors(s, mask, m, (1,), min_valid=mv)
            sur[m].append(e[1] if e is not None else np.nan)
    out["analysable"] = all(obs[m] is not None for m in ms)
    for m in ms:
        if obs[m] is None:
            out[f"z_m{m}"] = float("nan")
            continue
        sv = np.array(sur[m], float)
        out[f"sur_mean_m{m}"] = float(np.nanmean(sv))
        out[f"sur_sd_m{m}"] = float(np.nanstd(sv, ddof=1))
        out[f"obs_m{m}"] = float(obs[m][1])
        out[f"z_m{m}"] = float((np.nanmean(sv) - obs[m][1]) / max(np.nanstd(sv, ddof=1), 1e-12))
    zs = [out[f"z_m{m}"] for m in ms]
    out["zmax"] = float(np.nanmax(zs))
    marg = ms[int(np.nanargmax(zs))]
    out["m_argmax"] = int(marg)
    out["sd_at_argmax"] = out.get(f"sur_sd_m{marg}", float("nan"))
    out["degenerate"] = bool(np.isfinite(out["sd_at_argmax"]) and out["sd_at_argmax"] < 1e-6)
    out["G"] = float(obs[3][5] - obs[3][2]) if obs[3] is not None else float("nan")
    out["detected"] = bool(out["analysable"] and out["zmax"] >= cfg.masked_growth_z and out["G"] >= cfg.masked_growth_g)
    xs = x[~mask] if (~mask).any() else x
    out["mad_nn"] = float(np.median(np.abs(xs - np.median(xs))))
    out["sdnn_nn"] = float(np.std(xs))
    return out


# ---------------------------------------------------------------------------------------- B3
def spike_varpres(rr, labels, c, f, fs):
    """Variance-preserving design x'_k = mu_NN + sqrt(1 - f) (x_k - mu_NN) + sqrt(f) sigma_NN c_k for every
    interval k (mu_NN, sigma_NN over the window's non-masked intervals), then quantized to 1/fs."""
    rr = np.asarray(rr, float)
    m = H.k3_mask(labels)
    nn = rr[~m] if (~m).sum() >= 2 else rr
    mu, sd = float(nn.mean()), float(nn.std())
    x = mu + np.sqrt(1.0 - f) * (rr - mu) + np.sqrt(f) * sd * np.asarray(c, float)
    q = H.quantize_times(0.0, x, fs)
    return q, {"mu_nn": mu, "sd_nn": sd, "n_out_of_range": int(np.sum((q < H.LO) | (q > H.HI))),
               "n_nonpositive": int(np.sum(q <= 0)), "min_rr": float(q.min())}
