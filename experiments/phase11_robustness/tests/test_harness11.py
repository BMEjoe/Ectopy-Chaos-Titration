"""harness11: K3 diagnostic equals the frozen K3; sham masks respect the forbidden zone."""
import numpy as np

import final_pipeline as fp
from experiments.phase10_final import harness10 as H
from experiments.phase11_robustness import harness11 as H11


def _series(seed, n=512):
    rng = np.random.default_rng(11000 + seed)
    x = rng.uniform(0.2, 0.8)
    out = np.empty(n)
    for k in range(300 + n):
        x = 3.9 * x * (1 - x)
        if k >= 300:
            out[k - 300] = x
    rr = 0.8 + 0.05 * (out - out.mean()) / out.std() + 0.01 * rng.standard_normal(n)
    lab = np.array(["N"] * (n + 1))
    lab[rng.integers(1, n, 3)] = "V"
    return rr, lab


def test_k3_diag_matches_frozen():
    for s in range(3):
        rr, lab = _series(s)
        a = fp.masked_growth_chaos_test(rr, lab)
        b = H11.k3_diag(rr, lab)
        assert a["detected"] == b["detected"] and a["analysable"] == b["analysable"]
        for k in ("zmax", "G", "z_m2", "z_m3", "z_m4", "z_m5"):
            assert a[k] == b[k]


def test_sham_masks_avoid_forbidden():
    rng = np.random.default_rng(11000)
    lab = np.array(["N"] * 513)
    lab[[50, 300]] = "V"
    m = H.k3_mask(lab)
    fb = H11.forbidden(m, lab)
    assert fb[49] and fb[50] and fb[70] and not fb[71]
    sh = H11.sham_masks(m, fb, 5, rng)
    assert len(sh) == 5 and len({o for o, _ in sh}) == 5
    for _, s in sh:
        assert s.sum() == m.sum() and not np.any(s & fb)
