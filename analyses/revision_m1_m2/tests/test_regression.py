"""
Revision M1/M2 regression check (PLAN.md Section 1): with every switch at its default, the configurable test
reproduces R0 exactly, and with the Phase 12 switch sets it reproduces the Phase 12 readings exactly, both against the
frozen code and against stored Phase 10-12 outputs.  Must pass before any new run.

    PYTHONPATH=. python -m pytest -q analyses/revision_m1_m2/tests
"""
from __future__ import annotations

import json
import zlib
from pathlib import Path

import numpy as np
import pytest

from analyses.revision_m1_m2 import bench as B
from analyses.revision_m1_m2 import configs as CF
from analyses.revision_m1_m2 import vwkcfg as V
from experiments.phase10_final import data10 as D
from experiments.phase10_final import harness10 as H
from experiments.phase10_final import run10 as R10
from experiments.phase10_final import titration10 as T
from experiments.phase10_final import verify_titration as VT
from experiments.phase12_original import run12 as R12
from experiments.phase12_original import systems12 as S
from experiments.phase12_original import titration12 as T12

ROOT = Path(__file__).resolve().parents[3]
P10 = ROOT / "experiments" / "phase10_final" / "results"
P11 = ROOT / "experiments" / "phase11_robustness" / "results" / "conf"
P12 = ROOT / "experiments" / "phase12_original" / "results" / "conf"
T10_KEYS = ("analysable", "n_rows", "C_lin", "C_nl", "r_lin", "r_nl", "kappa_lin", "kappa_nl", "d_nl", "eps_lin",
            "eps_nl", "p", "nonlinear")
CONF_OK = (D.CONF_DIR / "nsr2db" / "RECORDS").exists() and (D.CONF_DIR / "chf2db" / "RECORDS").exists()
DEV_OK = (D.DEV_DIR / "nsrdb" / "RECORDS").exists() and (D.DEV_DIR / "chfdb" / "RECORDS").exists()


def jl(path):
    return [json.loads(line) for line in open(path) if line.strip()]


def same(a, b):
    """Exact equality after a JSON round trip (floats bit for bit)."""
    return json.loads(json.dumps(a)) == json.loads(json.dumps(b))


def series_bank():
    out = []
    for s in range(30):
        rng = np.random.default_rng(1000 + s)
        n = (512, 1000, 800)[s % 3]
        kind = s % 5
        if kind == 0:
            y = rng.standard_normal(n)
        elif kind == 1:
            y = np.zeros(n)
            e = rng.standard_normal(n)
            for i in range(1, n):
                y[i] = 0.8 * y[i - 1] + e[i]
        elif kind == 2:
            y = VT.logistic(3.6 + 0.01 * s, x0=rng.uniform(0.2, 0.8), n=n)
        elif kind == 3:
            y = 0.8 + 0.05 * rng.standard_normal(n)
            for i in rng.choice(np.arange(5, n - 5), 3, replace=False):
                y[i] *= 0.7
                y[i + 1] += 0.3 * y[i] / 0.7
        else:
            y = np.sin(0.3 * np.arange(n)) + 0.3 * rng.standard_normal(n)
        mask = np.zeros(n, bool)
        mask[rng.choice(n, 2, replace=False)] = True
        out.append((y, mask))
    return out


# ------------------------------------------------------------------------------------------- frozen code
def test_defaults_equal_titration10_raw_and_masked():
    for y, mask in series_bank():
        a = T.vwk_test(y)
        b = V.vwk_cfg(y, CF.get("R0"))
        assert {k: a[k] for k in a} == {k: b[k] for k in a}
        vr = T.rows_from_mask(mask)
        a = T.vwk_test(y, valid_rows=vr)
        b = V.vwk_cfg(y, CF.get("R0"), valid_rows=V.rows_mask(mask, CF.get("R0")))
        assert {k: a[k] for k in a} == {k: b[k] for k in a}


def test_rows_mask_default_equals_titration10():
    for _, mask in series_bank():
        assert np.array_equal(T.rows_from_mask(mask), V.rows_mask(mask, V.R0))


def test_phase12_readings_equal_titration12():
    for i, (y, mask) in enumerate(series_bank()[:18]):
        for name in ("R0", "R0-BP96F", "R0-BP96", "R0-S-FT"):
            cfg = CF.get({"R0-S-FT": "S4"}.get(name, name))
            a = T12.vwk12(y, name=name, rng=np.random.default_rng(i), details=True)
            b = V.vwk_cfg(y, cfg, rng=np.random.default_rng(i), details=True)
            assert {k: a[k] for k in a} == {k: b[k] for k in a}, name


def test_titrate_defaults_equal_titration10():
    for i, (y, mask) in enumerate(series_bank()[:12]):
        a = T.titrate(y, np.random.default_rng(i))
        b = V.titrate_cfg(y, np.random.default_rng(i), V.R0)
        assert a["NL_i"] == b["NL_i"] and a["positive"] == b["positive"]
        assert {k: a["base"][k] for k in a["base"]} == {k: b["base"][k] for k in a["base"]}


def test_flow_wrappers_reproduce_systems12():
    for sy in ("rossler", "duffing", "lorenz"):
        assert np.array_equal(B._flow_series(sy, 12000, 1.0), S.SYSTEMS[sy](12000))


def test_switch_configs_build():
    for s in CF.SINGLE + tuple(CF.COMPOSITES):
        cfg = CF.get(s)
        assert cfg.kl == cfg.m_max - 1
    assert CF.get("S5a").m_max == 84 and [CF.get("S5a").d_for(k) for k in range(1, 7)] == [8, 8, 6, 4, 3, 3]
    assert CF.get("S5b").m_max == 120 and CF.get("S5b").K == 119


# ------------------------------------------------------------------------------------------- stored outputs
def test_stored_phase12_b1_details():
    R = {r["seed"]: r for r in jl(P12 / "b1.jsonl") if not r.get("error")}
    for s in (12000, 12001, 12002):
        y = S.ecological(s)
        xi = S.correlated_noise(y, np.random.default_rng(R12.ss("b1", s, 1)))
        for p in R12.NOISE_B1:
            yp = y + p / 100.0 * y.std() * xi
            for rd in ("R0", "R0-BP96F", "R0-BP96"):
                st = R[s]["levels"][str(p)][rd]
                o = V.vwk_cfg(yp, CF.get(rd), details=True)
                for k in ("nonlinear", "F", "p", "r_lin", "r_nl", "C_lin", "C_nl", "F_calc_plain", "test"):
                    assert same(st[k], o[k]), (s, p, rd, k)


def test_stored_phase12_b2_max_noise():
    R = {(r["system"], r["seed"], r["reading"]): r for r in jl(P12 / "b2.jsonl") if not r.get("error")}
    for sy in ("logistic", "ecological", "gopy_sna"):
        for s in (12000, 12001):
            for rd in ("R0", "R0-BP96F", "R0-BP96"):
                o = B.t1_max_noise(sy, s, CF.get(rd))
                assert o["max_noise"] == R[(sy, s, rd)]["max_noise"], (sy, s, rd)


def test_stored_phase10_verification_nl():
    R = json.load(open(P10 / "verification" / "titration_verification.json"))["records"]
    get = {(r["system"], r["arg"]): r for r in R}
    for r in (3.7, 3.575, 3.565, 3.6, 3.8, 3.85, 3.95, 4.0):
        o = B.scan_nl(r, V.R0)
        assert o["NL_i"] == get[("logistic", r)]["NL_i"], r


def test_stored_phase10_synthetic_windows():
    R = [r for r in jl(P10 / "conf" / "synth.jsonl") if not r.get("error") and r["n"] == 512]
    R = sorted(R, key=lambda r: r["id"])[::max(1, len(R) // 20)][:20]
    for r in R:
        cond = tuple(r["cond"]) if r["kind"] == "vdp" else r["cond"]
        rr, types = H.synthetic_window(r["kind"], cond, r["seed"], r["n"])
        m = H.k3_mask(types)
        for a in ("raw", "masked", "edited"):
            st = r[f"TIT_{a}"]
            rng = np.random.default_rng(np.random.SeedSequence([H.ENTROPY, H.crc(rr), 77, 0]))
            if a == "raw":
                o = V.titrate_cfg(rr, rng, V.R0)
            elif a == "masked":
                o = V.titrate_cfg(rr, rng, V.R0, valid_rows=V.rows_mask(m, V.R0))
            else:
                xe, nnf = H.edit(rr, imask=m)
                if xe is None or nnf < H.NN_EDIT_MIN:
                    assert st.get("eligible") is False
                    continue
                o = V.titrate_cfg(xe, rng, V.R0)
            assert o["positive"] == st["positive"] and o["NL"] == st["NL"], (r["id"], a)
            if o["base"].get("analysable"):
                assert o["base"]["p"] == st["p"] and o["base"]["r_nl"] == st["r_nl"], (r["id"], a)
            d = B.arm(rr, m, a, V.R0)
            assert d["nonlinear"] == st["positive"], (r["id"], a)


def test_stored_phase11_explore_dose():
    from experiments.phase11_robustness import explore11 as E11
    from experiments.phase11_robustness import run11 as R11
    R = {r["id"]: r for r in jl(P11 / "explore_dose.jsonl") if not r.get("error")}
    for k in (1, 26):
        for s in range(11000, 11010):
            x, _ = R11.synth_series("null", "N1_linear_rr", s)
            xk, idx = E11.insert_isolated(x, np.random.default_rng(np.random.SeedSequence([20261103, 99, k, s, 1])), k)
            xk = H.quantize_times(0.0, xk, 360.0)
            assert idx == R[f"{k}|{s}"]["idx"]
            assert V.vwk_cfg(xk, V.R0)["nonlinear"] == R[f"{k}|{s}"]["positive"], (k, s)


@pytest.mark.skipif(not CONF_OK, reason="confirmatory data not downloaded")
def test_stored_phase10_real_windows_and_phase12_b4w():
    R = {r["id"]: r for r in jl(P10 / "conf" / "real.jsonl") if not r.get("error")}
    W = R12.window_list(R12.Ctx(False))
    W = W[::max(1, len(W) // 40)][:40]
    b4w = {r["id"].split("|", 1)[1]: r for r in jl(P12 / "b4w.jsonl") if not r.get("error")}
    for w in W:
        d = R12.win_data(R12.Ctx(False), w)
        rr, lab = d["rr"], d["lab"]
        m = H.k3_mask(lab)
        st = R[w["wid"]]
        for a in ("raw", "masked", "edited"):
            s = st[f"TIT_{a}"]
            rng = np.random.default_rng(np.random.SeedSequence([H.ENTROPY, H.crc(rr), 77, 0]))
            if a == "raw":
                o = V.titrate_cfg(rr, rng, V.R0)
            elif a == "masked":
                o = V.titrate_cfg(rr, rng, V.R0, valid_rows=V.rows_mask(m, V.R0))
            else:
                xe, nnf = H.edit(rr, imask=m)
                if xe is None or nnf < H.NN_EDIT_MIN:
                    assert s.get("eligible") is False
                    continue
                o = V.titrate_cfg(xe, rng, V.R0)
            assert o["positive"] == s["positive"] and o["NL"] == s["NL"], (w["wid"], a)
            if o["base"].get("analysable"):
                assert o["base"]["p"] == s["p"] and o["base"]["n_rows"] == s["n_rows"], (w["wid"], a)
            assert B.arm(rr, m, a, V.R0)["nonlinear"] == s["positive"]
        if w["wid"] in b4w:
            assert B.arm(rr, m, "raw", CF.get("R0-BP96"))["nonlinear"] == b4w[w["wid"]]["raw"]["R0-BP96"]["nonlinear"]


@pytest.mark.skipif(not DEV_OK, reason="development data not downloaded")
def test_stored_phase12_pm97():
    R = {(r["db"], r["rec"]): r for r in jl(P12 / "pm97.jsonl") if not r.get("error") and "excluded" not in r}
    for db, rec in (("nsrdb", D.records("nsrdb", D.DEV_DIR)[0]), ("chfdb", D.records("chfdb", D.DEV_DIR)[1])):
        st = R[(db, rec)]
        for rd in ("R0", "R0-BP96F", "R0-BP96"):
            o = B.pm_subject(db, rec, CF.get(rd), arms=("raw", "masked", "edited", "wu"), lengths=(500,))
            for sg, ss in zip(o["seg500"], st["seg500"]):
                assert sg["raw"]["nonlinear"] == ss["raw"][rd]["nonlinear"], (db, rec, rd, sg["k"])
                assert sg["masked"]["nonlinear"] == ss["masked"][rd]["nonlinear"], (db, rec, rd, sg["k"])
                if ss.get("edited_eligible"):
                    assert sg["edited"]["nonlinear"] == ss["edited"][rd]["nonlinear"]
                if ss.get("wu_eligible"):
                    assert sg["wu"]["nonlinear"] == ss["wu"][rd]["nonlinear"]
