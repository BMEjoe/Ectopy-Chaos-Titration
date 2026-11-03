"""HRT generator profiles hit their TO / TS targets when measured by hrt11's own TO / TS code."""
import numpy as np
import pytest

from experiments.phase11_robustness import hrt11 as H


@pytest.mark.parametrize("kind", ["normal", "blunted"])
@pytest.mark.parametrize("rr0", [0.6, 0.8, 1.0, 1.2])
def test_profile_targets(kind, rr0):
    h = H.parametric_profile(kind, rr0)
    to, ts = H.profile_metrics(h, rr0)
    to_t, ts_t = H.TARGETS[kind]
    assert abs(to - to_t) < 0.05 and abs(ts - ts_t) < 0.05
    assert abs(h[-1]) < 0.01                       # recovery by beat 20


@pytest.mark.parametrize("kind", ["normal", "blunted"])
def test_generated_tachogram_measured(kind):
    """Insert V3/V4 into a constant 0.8 s series, re-measure with tachograms(): TO, TS on target."""
    n = 200
    x = np.full(n, 0.8)
    v = "V3" if kind == "normal" else "V4"
    xe = H.insert_event(x, 100, v, {"c": 0.7, "c2": 0.7, "d": 1.05})
    t = np.concatenate([[0.0], np.cumsum(xe)])
    sym = np.array(["N"] * (n + 1)); sym[101] = "V"
    lab = np.where(sym == "N", "N", "X")
    tg = H.tachograms(t, sym, lab)
    assert len(tg) == 1
    s = H.subject_hrt(tg)
    to_t, ts_t = H.TARGETS[kind]
    assert abs(s["TO"] - to_t) < 0.05 and abs(s["TS"] - ts_t) < 0.05


def test_empirical_profile_roundtrip():
    """V2 with a given empirical profile reproduces that profile in the measured normalised tachogram."""
    rng = np.random.default_rng(11000)
    h = list(0.02 * np.sin(np.arange(20) / 3.0))
    x = np.full(200, 0.9)
    xe = H.insert_event(x, 80, "V2", {"c": 0.65, "c2": 0.7, "d": 1.0}, h_emp=h)
    t = np.concatenate([[0.0], np.cumsum(xe)])
    sym = np.array(["N"] * 201); sym[81] = "V"
    tg = H.tachograms(t, sym, np.where(sym == "N", "N", "X"))
    assert np.allclose(tg[0]["h"], h, atol=1e-12)
    assert rng is not None


def test_compensatory_pause_preserves_timing():
    x = np.linspace(0.7, 0.9, 100)
    for v in ("V1", "V6"):
        xe = H.insert_event(x, 40, v, {"c": 0.7, "c2": 0.75, "d": 1.05})
        assert np.isclose(np.sum(xe), np.sum(x))
