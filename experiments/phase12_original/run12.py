"""
Phase 12 runner (PREREGISTRATION.md).  Resumable: finished task ids are skipped.  Confirmatory runs refuse to start
unless PREREGISTRATION.md is committed, unmodified and pushed (Phase 3 guard).  --pilot runs a small subset on
synthetic pilot seeds (12900+) and mitdb only (no nsr2db / chf2db / nsrdb / chfdb) to time each part.

    PYTHONPATH=. python -m experiments.phase12_original.run12 --part tier1b --workers 4
    PYTHONPATH=. python -m experiments.phase12_original.run12 --part pm97 --workers 4
    PYTHONPATH=. python -m experiments.phase12_original.run12 --pilot --part all --workers 4

Parts:
  b1      BP96 Fig. 1 ecological model: 20 seeds x correlated noise 0-60 % (5 % steps) x readings (details)
  b2      BP96 Table 1 systems: 20 seeds x readings, largest correlated-noise % still detected (bisection)
  b4a     static transforms of Gaussian AR(1): phi 0.5 / 0.9 x n 1,000 / 512 x 100 seeds x 6 transforms x readings
  b4v     verification V1-V6 under R0-BP96F, R0-BP96 (titrated) and R0-S-FT, R0-S-AAFT (gated R0 NL)
  b4s     Phase 10 synthetic conditions (N1-N6, E3 couplets, E4 runs, S2 isolated 5 %) x 100 seeds x readings
  b4t     Phase 11 transplant variants V0, V1, V5, V6 in the 439 clean confirmatory windows x readings
  b4w     the 988 confirmatory windows, raw / masked / edited x readings
  pm97    PM97 replication on nsrdb / chfdb: 20,000 beats, 40 x 500 (+ 10 x 2,000) segments, arms x readings
  pm97nl  titration NL (R0) of the raw 500-beat segments (secondary)
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from experiments.phase10_final import data10 as D
from experiments.phase10_final import harness10 as H
from experiments.phase10_final import run10 as R10
from experiments.phase10_final import titration10 as T
from experiments.phase12_original import systems12 as S
from experiments.phase12_original import titration12 as T12

HERE = Path(__file__).resolve().parent
PREREG = HERE / "PREREGISTRATION.md"
METHODS_LISTED = ("TIT12",)
ENTROPY12 = S.ENTROPY12
PART_CODE = {"b1": 1, "b2": 2, "b4a": 3, "b4v": 4, "b4s": 5, "b4t": 6, "b4w": 7, "pm97": 8, "pm97nl": 9}
READ_ALL = ("R0", "R0-BP96F", "R0-BP96", "R0-S-FT", "R0-S-AAFT")
READ_BP = ("R0", "R0-BP96F", "R0-BP96")
TIER1B = ("b1", "b2", "b4a", "b4v", "b4s", "b4t", "b4w")
TIER1C = ("pm97", "pm97nl")
NOISE_B1 = tuple(range(0, 65, 5))
B2_HI, B2_STEPS = 200.0, 8
PM97_BEATS, PM97_STEP = 20000, 1000
SEG_LENS = (500, 2000)
COEF_KAPPA, COEF_D, COEF_N = 7, 3, 80
B4S_CONDS = ("N1_linear_rr", "N2_power_law", "N3_linear_rr_trend", "N4_linear_rr_step", "N5_linear_rr_warped",
             "N6_noisy_rsa", "E3_couplets_5pct", "E4_runs", "S2_ectopic_5pct")
TRANSPLANT_VARIANTS = ("V0", "V1", "V5", "V6")


def crc(s):
    return zlib.crc32(str(s).encode())


def ss(part, *key):
    return np.random.SeedSequence([ENTROPY12, PART_CODE[part], *[int(k) for k in key]])


def guard():
    from experiments.phase3_lle import run_phase3 as R3
    R3.PREREG = PREREG
    R3.check_preregistration(list(METHODS_LISTED))


class Ctx:
    def __init__(self, pilot):
        self.pilot = pilot
        self.out = HERE / "results" / ("pilot" if pilot else "conf")
        self.win_root, self.win_dbs = (D.DEV_DIR, ("mitdb",)) if pilot else (D.CONF_DIR, D.CONF_DBS)
        self.pm_root, self.pm_dbs = (D.DEV_DIR, ("mitdb",)) if pilot else (D.DEV_DIR, ("nsrdb", "chfdb"))
        self.seeds = range(12900, 12902) if pilot else range(12000, 12020)
        self.seeds100 = range(12900, 12903) if pilot else range(12000, 12100)


# ------------------------------------------------------------------------------------------------ helpers
def decide(y, name, rng=None, valid_rows=None, surr_source=None, details=False):
    o = T12.vwk12(y, valid_rows=valid_rows, name=name, rng=rng, surr_source=surr_source, details=details)
    keep = ("nonlinear", "analysable", "n_rows", "r_lin", "r_nl", "kappa_lin", "kappa_nl", "d_nl", "eps_lin", "eps_nl",
            "F", "df1", "df2", "test", "p", "surr_p_worst", "F_calc_plain", "F_crit", "MW_Z", "d_opt_gt1", "chi2_p_lin",
            "chi2_p_nl", "C_lin", "C_nl")
    return {k: o[k] for k in keep if k in o}


def readings_block(y, part, key, valid_rows=None, surr_source=None, names=READ_ALL):
    out = {}
    for ri, nm in enumerate(names):
        rng = np.random.default_rng(ss(part, *key, 100 + ri)) if "S-" in nm else None
        out[nm] = decide(y, nm, rng=rng, valid_rows=valid_rows, surr_source=surr_source)
    return out


def window_list(ctx):
    W = []
    for db in ctx.win_dbs:
        for si, rec in enumerate(D.records(db, ctx.win_root)):
            t = R10.series(ctx.win_root, db, rec)[0]
            for j, i0 in H.window_starts(t, n_per=R10.n_per(db)).items():
                W.append({"db": db, "rec": rec, "si": si, "j": j, "i0": i0, "wid": f"real|{db}|{rec}|{j}"})
    return W


def win_data(ctx, w):
    return R10._window(ctx.win_root, {"db": w["db"], "rec": w["rec"], "i0": w["i0"]}, H.L_PRIMARY)


# ------------------------------------------------------------------------------------------------ tasks
def tasks(ctx, part):
    if part == "b1":
        return [{"id": f"b1|{s}", "part": part, "seed": s} for s in ctx.seeds]
    if part == "b2":
        systems = list(S.SYSTEMS) if not ctx.pilot else ["logistic", "ecological", "lorenz"]
        return [{"id": f"b2|{sy}|{s}|{rd}", "part": part, "system": sy, "seed": s, "reading": rd}
                for sy in systems for s in ctx.seeds for rd in READ_BP]
    if part == "b4a":
        return [{"id": f"b4a|{phi}|{n}|{s}", "part": part, "phi": phi, "n": n, "seed": s}
                for phi in (0.5, 0.9) for n in (1000, 512) for s in ctx.seeds100]
    if part == "b4v":
        from experiments.phase10_final import verify_titration as VT
        VTT = VT.tasks()
        if ctx.pilot:
            VTT = VTT[::25]
        return [{"id": f"b4v|{rd}|{nm}|{a}", "part": part, "reading": rd, "system": nm, "arg": a}
                for rd in ("R0-BP96F", "R0-BP96", "R0-S-FT", "R0-S-AAFT") for nm, a in VTT]
    if part == "b4s":
        return [{"id": f"b4s|{c}|{s}", "part": part, "cond": c, "seed": s} for c in B4S_CONDS for s in ctx.seeds100]
    if part == "b4t":
        if ctx.pilot:
            W = [w for w in window_list(ctx)][:8]
            return [{"id": f"b4t|{w['wid']}", "part": part, "w": w, "group": "pilot", "subject": w["rec"]} for w in W]
        P11 = HERE.parent / "phase11_robustness" / "results" / "conf" / "transplant.jsonl"
        ids = {}
        for line in open(P11):
            r = json.loads(line)
            if not r.get("error"):
                ids[r["id"].split("|", 1)[1]] = r
        W = {w["wid"]: w for w in window_list(ctx)}
        return [{"id": f"b4t|{wid}", "part": part, "w": W[wid], "group": r["group"], "subject": r["subject"]}
                for wid, r in sorted(ids.items())]
    if part == "b4w":
        W = window_list(ctx)
        if ctx.pilot:
            W = W[::12]
        return [{"id": f"b4w|{w['wid']}", "part": part, "w": w} for w in W]
    if part in ("pm97", "pm97nl"):
        out = []
        for db in ctx.pm_dbs:
            recs = D.records(db, ctx.pm_root)
            if ctx.pilot:
                recs = recs[:2]
            out += [{"id": f"{part}|{db}|{rec}", "part": part, "db": db, "rec": rec} for rec in recs]
        return out
    raise ValueError(part)


# ------------------------------------------------------------------------------------------------ B1 / B2
def run_b1(ctx, task):
    s = task["seed"]
    y = S.ecological(s)
    xi = S.correlated_noise(y, np.random.default_rng(ss("b1", s, 1)))
    out = {"seed": s, "levels": {}}
    for p in NOISE_B1:
        yp = y + p / 100.0 * y.std() * xi
        out["levels"][str(p)] = {rd: decide(yp, rd, details=True) for rd in READ_BP}
    return out


def run_b2(ctx, task):
    y = S.SYSTEMS[task["system"]](task["seed"])
    xi = S.correlated_noise(y, np.random.default_rng(ss("b2", crc(task["system"]), task["seed"], 1)))
    det = lambda p: T12.vwk12(y + p / 100.0 * y.std() * xi, name=task["reading"])["nonlinear"]  # noqa: E731
    base = T12.vwk12(y, name=task["reading"], details=True)
    out = {"system": task["system"], "seed": task["seed"], "reading": task["reading"],
           "base": {k: base.get(k) for k in ("nonlinear", "kappa_nl", "d_nl", "r_lin", "r_nl", "F_calc_plain", "p")}}
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


# ------------------------------------------------------------------------------------------------ B4
def run_b4a(ctx, task):
    z = S.ar1(task["seed"], task["n"], task["phi"])
    out = {"phi": task["phi"], "n": task["n"], "seed": task["seed"], "transforms": {}}
    for ti, (nm, f) in enumerate(S.TRANSFORMS.items()):
        out["transforms"][nm] = readings_block(f(z), "b4a", (int(100 * task["phi"]), task["n"], task["seed"], ti))
    return out


def run_b4v(ctx, task):
    from experiments.phase10_final import verify_titration as VT
    name, arg, rd = task["system"], task["arg"], task["reading"]
    if name == "logistic":
        y, extra = VT.logistic(arg), {"LE": VT.logistic_le(arg)}
    elif name == "henon":
        y, extra = VT.henon(), {}
    elif name == "lorenz":
        y, extra = VT.lorenz(arg), {}
    elif name == "mg":
        y, extra = VT.mackey_glass(arg), {}
    else:
        y, extra = VT.controls(name, arg), {}
    if "S-" not in rd:
        rng = np.random.default_rng(np.random.SeedSequence([20261003, zlib.crc32(f"{name}|{arg}".encode())]))
        o = T12.titrate12(y, rng, name=rd)
        return {"reading": rd, "system": name, "arg": arg, "NL": o["NL"], "positive": o["positive"], **extra}
    # surrogate readings: NL = Phase 11 stored R0 NL (same series, same seeds) if the surrogate requirement holds
    r0 = R0_VERIFY().get(f"{name}|{arg}")
    g = decide(y, rd, rng=np.random.default_rng(ss("b4v", crc(f"{name}|{arg}"), crc(rd))))
    pos = bool(r0 and r0["positive"] and g["nonlinear"])
    return {"reading": rd, "system": name, "arg": arg, "NL": float(r0["NL"]) if pos else 0.0, "positive": pos,
            "R0_NL": r0["NL"] if r0 else None, "surr_p_worst": g.get("surr_p_worst"), **extra}


_R0V = {}


def R0_VERIFY():
    if not _R0V:
        p = HERE.parent / "phase11_robustness" / "results" / "conf" / "rdg_verify.jsonl"
        for line in open(p):
            r = json.loads(line)
            if r.get("reading") == "R0" and not r.get("error"):
                _R0V[f"{r['system']}|{r['arg']}"] = r
    return _R0V


def run_b4s(ctx, task):
    from experiments.phase11_robustness import run11 as R11
    x, types = R11.synth_series("null", task["cond"], task["seed"])
    return {"cond": task["cond"], "seed": task["seed"], "n_ectopic": int(np.sum(np.asarray(types) != "N")),
            "raw": readings_block(x, "b4s", (crc(task["cond"]), task["seed"]))}


def run_b4t(ctx, task):
    from experiments.phase11_robustness import harness11 as H11
    from experiments.phase11_robustness import hrt11 as HRT
    w = task["w"]
    d = win_data(ctx, w)
    x = d["rr"]
    key = H11.key_of(w["wid"])
    rng_ev = np.random.default_rng(H11.ss("transplant", key, 1))
    j = HRT.event_position(len(x), rng_ev)
    dr = HRT.draws_for(rng_ev)
    out = {"subject": task["subject"], "group": task["group"], "pos": j, "variants": {}}
    for vi, v in enumerate(TRANSPLANT_VARIANTS):
        xv = HRT.insert_event(x, j, v, dr, h_emp=None)
        if v != "V0":
            xv = H.quantize_times(0.0, xv, d["fs"])
        out["variants"][v] = readings_block(xv, "b4t", (key, vi))
    return out


def run_b4w(ctx, task):
    w = task["w"]
    d = win_data(ctx, w)
    rr, lab = d["rr"], d["lab"]
    m = H.k3_mask(lab)
    key = H11key(w["wid"])
    out = {"subject": f"{w['db']}:{w['rec']}", "db": w["db"], "j": w["j"], "n_masked": int(m.sum())}
    out["raw"] = readings_block(rr, "b4w", (key, 0))
    xe, nnf = H.edit(rr, imask=m)
    vr = T.rows_from_mask(m)
    out["masked"] = readings_block(rr, "b4w", (key, 1), valid_rows=vr, surr_source=xe if xe is not None else rr)
    if xe is not None and nnf >= H.NN_EDIT_MIN:
        out["edited"] = readings_block(xe, "b4w", (key, 2))
    out["edited_eligible"] = bool(xe is not None and nnf >= H.NN_EDIT_MIN)
    out["nn_fraction"] = float(nnf)
    return out


def H11key(wid):
    return crc(wid)


# ------------------------------------------------------------------------------------------------ PM97
def pm97_block(ctx, db, rec):
    """20,000 consecutive beat intervals: from midnight if the header has a clock start time and >= 20,000 beats
    follow midnight; otherwise the block (starts every 1,000 beats) with the lowest mean heart rate, computed from its
    intervals within [0.25, 3.0] s.  Pilot (mitdb): 2,000 beats from the first beat."""
    t, lab, sym, fs, info = R10.series(ctx.pm_root, db, rec)
    rr = np.diff(t)
    nb = 2000 if ctx.pilot else PM97_BEATS
    if len(rr) < nb:
        return {"excluded": f"fewer than {nb} beat intervals ({len(rr)})"}
    rule, start = None, None
    if info.get("base_time") and not ctx.pilot:
        h, mi, se = (float(v) for v in str(info["base_time"]).split(":"))
        t0 = h * 3600 + mi * 60 + se
        to_midnight = (24 * 3600 - t0) % (24 * 3600)
        i_mid = int(np.searchsorted(t, to_midnight))
        if i_mid + nb <= len(rr):
            rule, start = "midnight", i_mid
    if rule is None:
        if ctx.pilot:
            rule, start = "pilot_first", 0
        else:
            best = None
            for i0 in range(0, len(rr) - nb + 1, PM97_STEP):
                x = rr[i0:i0 + nb]
                xin = x[(x >= H.LO) & (x <= H.HI)]
                hr = 60.0 / xin.mean() if len(xin) else np.inf
                if best is None or hr < best[0]:
                    best = (hr, i0)
            rule, start = "lowest_mean_hr", best[1]
    i0, i1 = start, start + nb
    return {"rule": rule, "start_beat": int(i0), "end_beat": int(i1), "t_start": float(t[i0]), "t_end": float(t[i1]),
            "base_time": info.get("base_time"), "fs": fs, "rr": rr[i0:i1], "lab": lab[i0:i1 + 1], "sym": sym[i0:i1 + 1]}


def detrend(x):
    i = np.arange(len(x), dtype=float)
    b = np.polyfit(i, x, 1)
    return x - (b[0] * i + b[1]) + x.mean()


def coef80(y):
    """Fig. 3 reading: full VWK model kappa = 7, d = 3 (M = 120) by least squares on the standardised segment (rows
    n >= 83, as R0); coefficients of the first 80 nonlinear terms in Eq. 1 order."""
    y = T12._std(y)
    rows = np.arange(max(T.KAPPA_MAX, T.KAPPA_LIN_MAX), len(y))
    X = T.design(y, COEF_KAPPA, COEF_D, rows)
    c, *_ = np.linalg.lstsq(X, y[rows], rcond=None)
    return [float(v) for v in c[COEF_KAPPA + 1:COEF_KAPPA + 1 + COEF_N]]


def pm97_segments(blk, L, ctx):
    rr_raw = blk["rr"]
    x = detrend(rr_raw)
    lab, sym = blk["lab"], blk["sym"]
    segs = []
    for k in range(len(x) // L):
        a, b = k * L, (k + 1) * L
        xs, rs, lb, sy = x[a:b], rr_raw[a:b], lab[a:b + 1], sym[a:b + 1]
        keep = (rs >= H.LO) & (rs <= H.HI)
        m = H.k3_mask(lb)[keep]
        segs.append({"k": k, "raw": xs, "rawf": xs[keep], "m": m, "n_out_of_range": int((~keep).sum()),
                     "ectopy": int(np.sum(~np.isin(sy[1:], list(D.NORMAL)))), "mean_rr": float(rs[keep].mean())})
    return segs


def run_pm97(ctx, task):
    blk = pm97_block(ctx, task["db"], task["rec"])
    if "excluded" in blk:
        return {"db": task["db"], "rec": task["rec"], **blk}
    key = crc(f"{task['db']}|{task['rec']}")
    out = {k: blk[k] for k in ("rule", "start_beat", "end_beat", "t_start", "t_end", "base_time", "fs")}
    out.update(db=task["db"], rec=task["rec"], group="NSR" if task["db"] in ("nsrdb",) else "CHF",
               n_nonnormal_block=int(np.sum(~np.isin(blk["sym"][1:], list(D.NORMAL)))),
               n_out_of_range_block=int(np.sum((blk["rr"] < H.LO) | (blk["rr"] > H.HI))))
    for L in SEG_LENS:
        res = []
        for sg in pm97_segments(blk, L, ctx):
            r = {"k": sg["k"], "ectopy": sg["ectopy"], "n_out_of_range": sg["n_out_of_range"], "mean_rr": sg["mean_rr"]}
            r["raw"] = readings_block(sg["raw"], "pm97", (key, L, sg["k"], 0), names=READ_BP)
            r["rawf"] = readings_block(sg["rawf"], "pm97", (key, L, sg["k"], 1), names=READ_BP)
            m = sg["m"]
            r["masked"] = readings_block(sg["rawf"], "pm97", (key, L, sg["k"], 2), valid_rows=T.rows_from_mask(m),
                                         names=READ_BP)
            xe, nnf = H.edit(sg["rawf"], imask=m)
            r["edited_eligible"] = bool(xe is not None and nnf >= H.NN_EDIT_MIN)
            if r["edited_eligible"]:
                r["edited"] = readings_block(xe, "pm97", (key, L, sg["k"], 3), names=READ_BP)
            nn = H.wu_nn(sg["rawf"], m)
            r["wu_n"] = int(len(nn))
            r["wu_eligible"] = bool(len(nn) >= H.SEG_MIN_N)
            if r["wu_eligible"]:
                r["wu"] = readings_block(nn, "pm97", (key, L, sg["k"], 4), names=READ_BP)
            if L == 500 and r["raw"]["R0"]["nonlinear"]:
                r["coef80"] = coef80(sg["raw"])
            res.append(r)
        out[f"seg{L}"] = res
    return out


def run_pm97nl(ctx, task):
    blk = pm97_block(ctx, task["db"], task["rec"])
    if "excluded" in blk:
        return {"db": task["db"], "rec": task["rec"], **blk}
    key = crc(f"{task['db']}|{task['rec']}")
    res = []
    for sg in pm97_segments(blk, 500, ctx):
        o = T.titrate(sg["raw"], np.random.default_rng(ss("pm97nl", key, sg["k"])))
        res.append({"k": sg["k"], "NL": float(o["NL"]), "positive": bool(o["positive"])})
    return {"db": task["db"], "rec": task["rec"], "segments": res}


RUN = {"b1": run_b1, "b2": run_b2, "b4a": run_b4a, "b4v": run_b4v, "b4s": run_b4s, "b4t": run_b4t, "b4w": run_b4w,
       "pm97": run_pm97, "pm97nl": run_pm97nl}


def execute(args):
    pilot, task = args
    ctx = Ctx(pilot)
    t0, c0 = time.time(), time.process_time()
    out = {"id": task["id"], "part": task["part"]}
    try:
        out.update(RUN[task["part"]](ctx, task))
    except Exception as exc:                                    # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
        out["traceback"] = traceback.format_exc(limit=4)
    out["t_wall"], out["t_cpu"] = time.time() - t0, time.process_time() - c0
    return out


def _done(path):
    done = set()
    if path.exists():
        for line in open(path):
            try:
                r = json.loads(line)
                if not r.get("error"):
                    done.add(r["id"])
            except Exception:                                   # noqa: BLE001
                pass
    return done


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="tier1b")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--pilot", action="store_true")
    a = ap.parse_args(argv)
    if not a.pilot:
        guard()
    ctx = Ctx(a.pilot)
    ctx.out.mkdir(parents=True, exist_ok=True)
    parts = {"tier1b": TIER1B, "tier1c": TIER1C, "all": TIER1B + TIER1C}.get(a.part, tuple(a.part.split(",")))
    for part in parts:
        path = ctx.out / f"{part}.jsonl"
        done = _done(path)
        todo = [t for t in tasks(ctx, part) if t["id"] not in done]
        todo.sort(key=lambda t: crc(t["id"]))
        print(f"{part}: {len(todo)} tasks to run ({len(done)} done)", flush=True)
        t0 = time.time()
        with ProcessPoolExecutor(a.workers) as ex, open(path, "a") as f:
            for i, r in enumerate(ex.map(execute, [(a.pilot, t) for t in todo], chunksize=1), 1):
                f.write(json.dumps(r, default=float) + "\n")
                f.flush()
                if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                    el = (time.time() - t0) / 60
                    print(f"  {part} {i}/{len(todo)} {el:.1f} min; eta {el / i * (len(todo) - i):.1f} min", flush=True)
    print("RUN12 DONE", parts, flush=True)


if __name__ == "__main__":
    main()
