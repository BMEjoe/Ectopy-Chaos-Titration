"""Phase 11 EXPLORATORY analysis (PREREGISTRATION.md, Amendment 1): ectopy dose-response on N1.

Not confirmatory: no prediction, no verdict. k isolated ventricular-type beats (Phase 5 add_ectopics rule:
c ~ U(0.6, 0.8), full compensatory pause, indices 1..N-2, spacing >= 3) in N1 linear RR (512 intervals,
seeds 11000-11099), re-quantized to 1/360 s, frozen titration (R0, raw arm).  The titration noise seed is shared
by every k of a series (paired across doses).

    PYTHONPATH=. python -m experiments.phase11_robustness.explore11 --workers 4       # run (resumable)
    PYTHONPATH=. python -m experiments.phase11_robustness.explore11 --analyse         # summary JSON
"""
import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "conf" / "explore_dose.jsonl"
SUMMARY = HERE / "results" / "conf" / "analysis" / "explore_dose.json"
KS = (0, 1, 2, 4, 8, 16, 26, 51)
SEEDS = range(11000, 11100)
ENTROPY11 = 20261103
CODE = 99
FS_SYN = 360.0


def insert_isolated(x, rng, k):
    """Phase 5 add_ectopics with a fixed count k (same logic: prematurity, compensatory pause, spacing)."""
    from experiments.phase5_rr import config as C
    x = np.array(x, dtype=float)
    chosen = []
    while len(chosen) < k:
        i = int(rng.integers(1, len(x) - 1))
        if all(abs(i - j) >= C.ECTOPIC_MIN_SPACING for j in chosen):
            chosen.append(i)
    chosen.sort()
    for i in chosen:
        c = rng.uniform(*C.ECTOPIC_PREMATURITY_RANGE)
        short = c * x[i]
        x[i + 1] += x[i] - short
        x[i] = short
    return x, chosen


def task(args):
    k, seed = args
    from experiments.phase10_final import harness10 as H
    from experiments.phase10_final import titration10 as T
    from experiments.phase11_robustness import run11 as R
    t0 = time.time()
    out = {"id": f"{k}|{seed}", "k": k, "seed": seed}
    try:
        x, _ = R.synth_series("null", "N1_linear_rr", seed)
        rng_ev = np.random.default_rng(np.random.SeedSequence([ENTROPY11, CODE, k, seed, 1]))
        xk, idx = insert_isolated(x, rng_ev, k)
        if k > 0:
            xk = H.quantize_times(0.0, xk, FS_SYN)
        o = T.titrate(xk, np.random.default_rng(np.random.SeedSequence([ENTROPY11, CODE, seed, 2])))
        out.update(positive=bool(o["positive"]), NL=float(o["NL"]), idx=[int(i) for i in idx])
    except Exception as exc:                                     # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    out["t"] = time.time() - t0
    return out


def run(workers):
    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["id"])
    todo = [(k, s) for k in KS for s in SEEDS if f"{k}|{s}" not in done]
    print(f"explore_dose: {len(todo)} tasks to run ({len(done)} done)", flush=True)
    with ProcessPoolExecutor(workers) as ex, OUT.open("a") as f:
        for n, r in enumerate(ex.map(task, todo, chunksize=2), 1):
            f.write(json.dumps(r) + "\n")
            f.flush()
            if n % 50 == 0:
                print(f"  {n}/{len(todo)}", flush=True)
    print("EXPLORE DONE", flush=True)


def analyse():
    from experiments.phase10_final import analysis10 as A10
    rows = [json.loads(line) for line in OUT.read_text().splitlines() if line.strip()]
    res = {"label": "EXPLORATORY (PREREGISTRATION.md Amendment 1); no prediction, no verdict",
           "n_errors": sum(1 for r in rows if r.get("error")), "by_k": {}}
    for k in KS:
        rk = [r for r in rows if r["k"] == k and not r.get("error")]
        pos = [r for r in rk if r["positive"]]
        res["by_k"][str(k)] = {"n": len(rk), "k_pos": len(pos), "rate": len(pos) / len(rk) if rk else None,
                               "wilson": list(A10.wilson(len(pos), len(rk))) if rk else None,
                               "mean_NL_pos": float(np.mean([r["NL"] for r in pos])) if pos else None}
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(res, indent=1))
    for k, v in res["by_k"].items():
        print(k, v["k_pos"], "/", v["n"], v["wilson"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    analyse() if a.analyse else run(a.workers)
