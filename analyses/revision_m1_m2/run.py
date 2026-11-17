"""
Revision M1/M2 runner (EXPLORATORY; PLAN.md).  Resumable: one JSON line per task, finished ids skipped.  Refuses to run
(except --pilot) unless analyses/revision_m1_m2/PLAN.md is committed, unmodified and pushed (Phase 3 guard).

    PYTHONPATH=. python -m analyses.revision_m1_m2.run --pilot --part all --workers 4      # timing only
    PYTHONPATH=. python -m analyses.revision_m1_m2.run --part p1sel,abl --workers 4
    PYTHONPATH=. python -m analyses.revision_m1_m2.run --part abl2 --workers 4              # R0-ORIG + leave-one-out
    PYTHONPATH=. python -m analyses.revision_m1_m2.run --part head,m2 --workers 4

Parts
  p1sel   P1 sampling selection for the flows: m x 5 design seeds, noise-free, R0
  gen     generate and cache the Table 1 flow series (seeds 12000-12099, m = 1 and the chosen m)
  abl     benchmarks (t1, t1p1, eco, log, scan, pm) for R0, S1-S7, R0-BP96F, R0-BP96, R0-TEXT
  abl2    the same benchmarks for R0-ORIG and its leave-one-out variants
  ver     Phase 10 verification controls V5-V6 (V1-V4 come from the logistic scan) for every configuration
  head    headline units (synth, win, tr, seg, pm97) for R0, R0-BP96, R0-ORIG
  m2      M2 beat-type x k x placement windows: R0 and R0-ORIG diagnostics, R0-BP96 decision
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

from analyses.revision_m1_m2 import bench as B  # noqa: E402
from analyses.revision_m1_m2 import configs as CF  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PLAN = HERE / "PLAN.md"
RES = ROOT / "results" / "revision_m1_m2"
ABL_CONFIGS = ("R0",) + CF.SINGLE + tuple(CF.COMPOSITES)
HEAD_CONFIGS = ("R0", "R0-BP96", "R0-ORIG")
T1_SEEDS = range(12000, 12100)
P1_SEEDS = range(12900, 12905)
LOG_SEEDS = range(0, 100)
SCAN_R = tuple(round(3.5 + 0.005 * i, 3) for i in range(101))
SYNTH_SEEDS = range(10100, 10200)


def guard():
    from experiments.phase3_lle import run_phase3 as R3
    R3.PREREG = PLAN
    R3.check_preregistration(["REVM1M2"])


def out_dir(pilot):
    return RES / ("pilot" if pilot else "raw")


# ------------------------------------------------------------------------------------------- P1 selection
def p1_choice():
    p = RES / "p1_choice.json"
    return json.load(open(p)) if p.exists() else None


def write_p1_choice():
    R = [json.loads(line) for line in open(out_dir(False) / "p1sel.jsonl") if line.strip()]
    R = [r for r in R if not r.get("error")]
    out = {"rule": "PLAN.md 3.2 P1: m maximising the median (eps_lin - eps_nl) over design seeds 12900-12904, "
                   "noise-free, R0; ties -> smaller m", "systems": {}}
    for sy in B.FLOWS:
        meds = {}
        for m in B.P1_M:
            v = [r["diff"] for r in R if r["system"] == sy and r["m"] == m]
            if len(v) == len(P1_SEEDS):
                meds[m] = float(np.median(v))
        best = max(sorted(meds), key=lambda m: (meds[m], -m))
        out["systems"][sy] = {"m_opt": best, "median_diff": {str(k): v for k, v in meds.items()}}
    (RES / "p1_choice.json").write_text(json.dumps(out, indent=1))
    return out


# ------------------------------------------------------------------------------------------- tasks
def tasks(part, pilot):
    T = []
    if part == "p1sel":
        systems = ("rossler",) if pilot else B.FLOWS
        ms = (0.25, 4.0) if pilot else B.P1_M
        seeds = (12950,) if pilot else P1_SEEDS
        return [{"id": f"p1sel|{sy}|{m:g}|{s}", "part": part, "system": sy, "m": m, "seed": s}
                for sy in systems for m in ms for s in seeds]
    if part == "gen":
        p1 = p1_choice() if not pilot else {"systems": {"rossler": {"m_opt": 2.0}}}
        systems = ("rossler", "mackey_glass_100") if pilot else B.FLOWS
        seeds = (12950, 12951) if pilot else T1_SEEDS
        for sy in systems:
            ms = {1.0} | ({p1["systems"][sy]["m_opt"]} if p1 and sy in p1["systems"] else set())
            T += [{"id": f"gen|{sy}|{m:g}|{s}", "part": part, "system": sy, "m": m, "seed": s} for m in sorted(ms)
                  for s in seeds]
        return T
    if part in ("abl", "abl2"):
        if part == "abl":
            cfgs = ABL_CONFIGS
        else:
            sw = CF.r0orig_switches()
            cfgs = ("R0-ORIG",) + tuple(f"R0-ORIG-minus-{s}" for s in sw if len(sw) > 1)
        t1sys = ("logistic", "rossler", "mackey_glass_100") if pilot else B.T1_SYSTEMS + B.T1_SECONDARY
        seeds = (12950, 12951) if pilot else T1_SEEDS
        logs = (950, 951) if pilot else LOG_SEEDS
        scan = (3.7, 3.9) if pilot else SCAN_R
        pmdb = (("mitdb", ("100", "101")),) if pilot else (("nsrdb", None), ("chfdb", None))
        p1 = p1_choice() if not pilot else {"systems": {"rossler": {"m_opt": 2.0}}}
        for c in cfgs:
            T += [{"id": f"t1|{c}|{sy}|{s}", "part": "t1", "cfg": c, "system": sy, "seed": s, "m": 1.0}
                  for sy in t1sys for s in seeds]
            if p1 is not None:
                for sy, v in p1["systems"].items():
                    if v["m_opt"] != 1.0 and (sy in t1sys or not pilot):
                        T += [{"id": f"t1p1|{c}|{sy}|{s}", "part": "t1p1", "cfg": c, "system": sy, "seed": s,
                               "m": v["m_opt"]} for s in seeds]
            T += [{"id": f"eco|{c}|{s}", "part": "eco", "cfg": c, "seed": s} for s in seeds]
            T += [{"id": f"log|{c}|{s}", "part": "log", "cfg": c, "seed": s} for s in logs]
            T += [{"id": f"scan|{c}|{r}", "part": "scan", "cfg": c, "r": r} for r in scan]
            for db, recs in pmdb:
                recs = recs or B.records(db, B.D.DEV_DIR)
                T += [{"id": f"pm|{c}|{db}|{rec}", "part": "pm", "cfg": c, "db": db, "rec": rec} for rec in recs]
        return T
    if part == "ver":
        from experiments.phase10_final import verify_titration as VT
        cfgs = ABL_CONFIGS + ("R0-ORIG",)
        VTT = [(n, a) for n, a in VT.tasks() if n != "logistic"]
        if pilot:
            VTT = VTT[::20]
        return [{"id": f"ver|{c}|{n}|{a}", "part": "ver", "cfg": c, "system": n, "arg": a} for c in cfgs for n, a in VTT]
    if part == "head":
        return head_tasks(pilot)
    if part == "m2":
        return m2_tasks(pilot)
    raise ValueError(part)


def head_tasks(pilot):
    from experiments.phase10_final import run10 as R10
    from experiments.phase12_original import run12 as R12
    T = []
    cfgs = ("R0", "R0-BP96") if pilot else HEAD_CONFIGS
    ctx = R12.Ctx(pilot)
    if pilot:
        conds = [("null", "N1_linear_rr"), ("null", "E3_couplets_5pct"), ("vdp", [26, "E3_10"])]
        syn = [{"kind": k, "cond": c, "seed": 12950} for k, c in conds]
    else:
        syn = [{"kind": t["kind"], "cond": t["cond"], "seed": t["seed"]} for t in R10.tasks("conf", "synth")
               if t["n"] == 512]
    W = R12.window_list(ctx)
    if pilot:
        W = W[:6]
    TR = R12.tasks(ctx, "b4t")
    if pilot:
        TR = TR[:4]
    if pilot:
        segsub = [("mitdb", r) for r in ("100", "101")]
        pmsub = [("mitdb", r) for r in ("100", "101")]
    else:
        segsub = [(db, rec) for db in ("nsr2db", "chf2db") for rec in B.records(db, B.D.CONF_DIR)]
        pmsub = [(db, rec) for db in ("nsrdb", "chfdb") for rec in B.records(db, B.D.DEV_DIR)]
    for c in cfgs:
        for s in syn:
            cid = s["cond"] if s["kind"] == "null" else f"vdp{s['cond'][0]}_{s['cond'][1]}"
            T.append({"id": f"synth|{c}|{cid}|{s['seed']}", "part": "synth", "cfg": c, **s})
        T += [{"id": f"win|{c}|{w['wid']}", "part": "win", "cfg": c, "w": w} for w in W]
        T += [{"id": f"tr|{c}|{t['w']['wid']}", "part": "tr", "cfg": c, "w": t["w"], "group": t["group"],
               "subject": t["subject"]} for t in TR]
        T += [{"id": f"seg|{c}|{db}|{rec}", "part": "seg", "cfg": c, "db": db, "rec": rec} for db, rec in segsub]
        T += [{"id": f"pm97|{c}|{db}|{rec}", "part": "pm97", "cfg": c, "db": db, "rec": rec} for db, rec in pmsub]
    return T


def m2_tasks(pilot):
    from analyses.revision_m1_m2 import m2gen as M
    seeds = (12950, 12951) if pilot else M.SEEDS
    ks = (1, 51) if pilot else M.KS
    T = []
    for s in seeds:
        for k in ks:
            if k == 0:
                T.append({"id": f"m2|none|random|0|{s}", "part": "m2", "type": "none", "placement": "random", "k": 0,
                          "seed": s})
                continue
            for ty in M.TYPES:
                T.append({"id": f"m2|{ty}|random|{k}|{s}", "part": "m2", "type": ty, "placement": "random", "k": k,
                          "seed": s})
                if k in M.K_PERIODIC:
                    T.append({"id": f"m2|{ty}|periodic|{k}|{s}", "part": "m2", "type": ty, "placement": "periodic",
                              "k": k, "seed": s})
    return T


# ------------------------------------------------------------------------------------------- execution
def run_task(task, pilot):
    p = task["part"]
    if p == "p1sel":
        return B.p1_score(task["system"], task["seed"], task["m"])
    if p == "gen":
        y = B.t1_series(task["system"], task["seed"], task["m"])
        return {"n": int(len(y)), "sd": float(np.std(y))}
    cfg = CF.get(task["cfg"]) if "cfg" in task else None
    if p in ("t1", "t1p1"):
        return B.t1_max_noise(task["system"], task["seed"], cfg, task["m"])
    if p == "eco":
        return B.eco_levels(task["seed"], cfg)
    if p == "log":
        return B.log_nl(task["seed"], cfg)
    if p == "scan":
        return B.scan_nl(task["r"], cfg)
    if p == "pm":
        return B.pm_subject(task["db"], task["rec"], cfg, arms=("raw",), lengths=(500, 2000), pilot=pilot)
    if p == "ver":
        return B.verify_control(task["system"], task["arg"], cfg)
    if p == "synth":
        return run_synth(task, cfg)
    if p == "win":
        return run_win(task, cfg, pilot)
    if p == "tr":
        return run_tr(task, cfg, pilot)
    if p == "seg":
        return run_seg(task, cfg, pilot)
    if p == "pm97":
        return B.pm_subject(task["db"], task["rec"], cfg, arms=("raw", "rawf", "masked", "edited", "wu"),
                            lengths=(500, 2000), pilot=pilot)
    if p == "m2":
        return run_m2(task, pilot)
    raise ValueError(p)


def run_synth(task, cfg):
    from experiments.phase10_final import harness10 as H
    cond = tuple(task["cond"]) if task["kind"] == "vdp" else task["cond"]
    rr, types = synth_window_cached(task["kind"], cond, task["seed"])
    m = H.k3_mask(types)
    key = (B.PART["synth"], B.crc(task["id"].split("|", 2)[2]))
    return {"n_ectopic": int(np.sum(np.asarray(types) != "N")),
            **{a: B.arm(rr, m, a, cfg, B.surr_rng(cfg, *key, i)) for i, a in enumerate(("raw", "masked", "edited"))}}


def synth_window_cached(kind, cond, seed, n=512):
    """harness10.synthetic_window, cached on disk (the van der Pol generator takes ~5 s per window)."""
    from experiments.phase10_final import harness10 as H
    cid = cond if kind == "null" else f"vdp{cond[0]}_{cond[1]}"
    p = B.CACHE.parent / "synth" / f"{cid}_{seed}_{n}.npz"
    if p.exists():
        z = np.load(p)
        return z["rr"], z["types"]
    rr, types = H.synthetic_window(kind, cond, seed, n)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f"{p.stem}.{os.getpid()}.tmp.npz")
    np.savez(tmp, rr=np.asarray(rr, float), types=np.asarray(types))
    os.replace(tmp, p)
    return np.asarray(rr, float), np.asarray(types)


def run_win(task, cfg, pilot):
    from experiments.phase10_final import harness10 as H
    from experiments.phase12_original import run12 as R12
    w = task["w"]
    d = R12.win_data(R12.Ctx(pilot), w)
    rr, lab = d["rr"], d["lab"]
    m = H.k3_mask(lab)
    key = (B.PART["win"], B.crc(w["wid"]))
    return {"n_masked": int(m.sum()),
            **{a: B.arm(rr, m, a, cfg, B.surr_rng(cfg, *key, i)) for i, a in enumerate(("raw", "masked", "edited", "wu"))}}


def run_tr(task, cfg, pilot):
    from experiments.phase10_final import harness10 as H
    from experiments.phase11_robustness import harness11 as H11
    from experiments.phase11_robustness import hrt11 as HRT
    from experiments.phase12_original import run12 as R12
    w = task["w"]
    d = R12.win_data(R12.Ctx(pilot), w)
    x = d["rr"]
    key = H11.key_of(w["wid"])
    rng_ev = np.random.default_rng(H11.ss("transplant", key, 1))
    j = HRT.event_position(len(x), rng_ev)
    dr = HRT.draws_for(rng_ev)
    out = {"subject": task["subject"], "group": task["group"], "pos": j}
    for vi, v in enumerate(R12.TRANSPLANT_VARIANTS):
        xv = HRT.insert_event(x, j, v, dr, h_emp=None)
        if v != "V0":
            xv = H.quantize_times(0.0, xv, d["fs"])
        out[v] = B.arm(xv, None, "raw", cfg, B.surr_rng(cfg, B.PART["tr"], key, vi))
    return out


def run_seg(task, cfg, pilot):
    from experiments.phase10_final import data10 as D
    from experiments.phase10_final import harness10 as H
    from experiments.phase10_final import run10 as R10
    root = D.DEV_DIR if pilot else D.CONF_DIR
    t, lab, sym, fs, info = R10.series(root, task["db"], task["rec"])
    segs = []
    for s in H.segments(t, lab):
        rec = {"seg": s["seg"], "gap": s["gap"]}
        if not s["gap"]:
            rr_all = np.diff(t[s["i0"]:s["i1"] + 1])
            lb = lab[s["i0"]:s["i1"] + 1]
            sy = sym[s["i0"]:s["i1"] + 1]
            keep = (rr_all >= H.LO) & (rr_all <= H.HI)
            rr, m = rr_all[keep], H.k3_mask(lb)[keep]
            rec["burden"] = int(np.sum(~np.isin(sy, list(D.NORMAL))))
            rec["n_masked"] = int(m.sum())
            key = (B.PART["seg"], B.crc(f"{task['db']}|{task['rec']}"), s["seg"])
            for i, a in enumerate(("raw", "masked", "edited", "wu")):
                rec[a] = B.arm(rr, m, a, cfg, B.surr_rng(cfg, *key, i), keep=("nonlinear", "analysable"))
        segs.append(rec)
    return {"subject": f"{task['db']}:{task['rec']}", "group": "NSR" if "nsr" in task["db"] else "CHF",
            "segments": segs}


def run_m2(task, pilot):
    from analyses.revision_m1_m2 import m2gen as M
    x = M.base_window(task["seed"]) if not pilot else M.base_window(task["seed"])
    ty = task["type"]
    y, aff, pos = M.insert(x, task["k"], task["seed"], "ventricular" if ty == "none" else ty, task["placement"])
    out = {"type": ty, "placement": task["placement"], "k": task["k"], "seed": task["seed"], "pos": pos,
           "n": int(len(y)), "n_aff_intervals": int(aff.sum())}
    rep = task["seed"] in M.REP_SEEDS and ty in ("none", "ventricular")
    names = ("R0", "R0-ORIG") if not pilot else ("R0", "S2")
    for ci, nm in enumerate(names):
        cfg = CF.get(nm)
        rs = np.random.default_rng(np.random.SeedSequence([M.ENTROPY, 13, 2, M.TYPE_CODE[ty],
                                                           M.PLACE_CODE[task["placement"]], task["k"], task["seed"]]))
        out[nm] = M.diagnose(y, aff, cfg, rs, rep=rep, k=task["k"], width=M.WIDTH[ty])
    out["R0-BP96"] = {"nonlinear": bool(B.V.vwk_cfg(y, CF.get("R0-BP96"))["nonlinear"])}
    return out


def execute(args):
    pilot, task = args
    t0, c0 = time.time(), time.process_time()
    out = {"id": task["id"], "part": task["part"]}
    for k in ("cfg", "system", "seed", "m", "r", "db", "rec", "kind", "cond", "type", "placement", "k"):
        if k in task:
            out[k] = task[k]
    try:
        out.update(run_task(task, pilot))
    except Exception as exc:                                    # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
        out["traceback"] = traceback.format_exc(limit=6)
    out["t_wall"], out["t_cpu"] = time.time() - t0, time.process_time() - c0
    return out


def _default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


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


def run_part(part, pilot, workers, limit=None):
    od = out_dir(pilot)
    od.mkdir(parents=True, exist_ok=True)
    path = od / f"{part}.jsonl"
    done = _done(path)
    todo = [t for t in tasks(part, pilot) if t["id"] not in done]
    if limit:
        todo = todo[:limit]
    todo.sort(key=lambda t: B.crc(t["id"]))
    print(f"{part}: {len(todo)} tasks to run ({len(done)} done)", flush=True)
    if not todo:
        return
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex, open(path, "a") as f:
        for i, r in enumerate(ex.map(execute, [(pilot, t) for t in todo], chunksize=1), 1):
            f.write(json.dumps(r, default=_default) + "\n")
            f.flush()
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                el = (time.time() - t0) / 60
                print(f"  {part} {i}/{len(todo)} {el:.1f} min; eta {el / i * (len(todo) - i):.1f} min", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="abl")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    if not a.pilot:
        guard()
    parts = ("p1sel", "gen", "abl", "head", "m2") if a.part == "all" else tuple(a.part.split(","))
    for part in parts:
        run_part(part, a.pilot, a.workers, a.limit)
        if part == "p1sel" and not a.pilot:
            print(json.dumps(write_p1_choice(), indent=1), flush=True)
    print("RUN DONE", parts, flush=True)


if __name__ == "__main__":
    main()
