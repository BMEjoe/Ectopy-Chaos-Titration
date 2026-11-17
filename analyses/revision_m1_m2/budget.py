"""
Revision M1/M2: runtime projection from the pilot (PLAN.md Section 7).  Mean CPU per task (by part and configuration)
x the number of full tasks / workers.  Unknown at pilot time: P1's chosen sampling (worst case assumed: every flow
resampled) and R0-ORIG (assumed as costly as the most costly single switch).

    PYTHONPATH=. python -m analyses.revision_m1_m2.budget
"""
from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from analyses.revision_m1_m2 import bench as B
from analyses.revision_m1_m2 import m2gen as M
from analyses.revision_m1_m2 import run as RN


def main():
    R = []
    for part in ("p1sel", "gen", "abl", "head", "m2"):
        p = RN.out_dir(True) / f"{part}.jsonl"
        if p.exists():
            R += [json.loads(line) for line in open(p) if line.strip()]
    err = [r for r in R if r.get("error")]
    cpu = defaultdict(list)
    for r in R:
        if not r.get("error"):
            cpu[(r["part"], r.get("cfg"))].append(r["t_cpu"])
    mean = {k: float(np.mean(v)) for k, v in cpu.items()}
    lines = [f"pilot tasks {len(R)}, errors {len(err)}"]
    for e in err[:5]:
        lines.append(f"  ERROR {e['id']}: {e['error']}")
    n_t1 = len(B.T1_SYSTEMS + B.T1_SECONDARY) * len(RN.T1_SEEDS)
    n_t1p1 = len(B.FLOWS) * len(RN.T1_SEEDS)
    full = {"t1": n_t1, "t1p1": n_t1p1, "eco": len(RN.T1_SEEDS), "log": len(RN.LOG_SEEDS), "scan": len(RN.SCAN_R),
            "pm": 33}
    pilot_t1 = {"logistic": 0, "rossler": 0, "mackey_glass_100": 0}
    tot = 0.0
    rows = []
    max_switch = defaultdict(float)
    for c in RN.ABL_CONFIGS:
        sec = 0.0
        for p_, n in full.items():
            v = mean.get((p_, c))
            if v is None:
                v = mean.get(("t1", c), 0.0) if p_ == "t1p1" else 0.0
            # pm pilot used a 2,000-beat mitdb block (1 + 4 segments); full: 20,000 beats (40 + 10 segments)
            sec += v * n * (10.0 if p_ == "pm" else 1.0)
            if c != "R0":
                max_switch[p_] = max(max_switch[p_], v / max(mean.get((p_, "R0"), 1e-9), 1e-9))
        rows.append((c, sec))
        tot += sec
    abl2 = 4 * max(s for _, s in rows if not CF_is_surr(_))
    p1 = mean.get(("p1sel", None), 3.0) * len(B.FLOWS) * len(B.P1_M) * len(RN.P1_SEEDS) * 2.0
    gen = mean.get(("gen", None), 3.0) * len(B.FLOWS) * len(RN.T1_SEEDS) * 2
    ratio = max(max_switch.values()) if max_switch else 2.0
    head_full = {"synth": 4100, "win": 988, "tr": 439, "seg": 83, "pm97": 33}
    head_scale = {"synth": 1.0, "win": 1.0, "tr": 1.0, "seg": 60.0, "pm97": 10.0}   # pilot seg/pm97 on short mitdb records
    head = 0.0
    for p_, n in head_full.items():
        r0 = mean.get((p_, "R0"), 0.0)
        bp = mean.get((p_, "R0-BP96"), r0)
        head += (r0 + bp + min(ratio, 11.0) * r0) * n * head_scale[p_]
    n_m2 = len(RN.m2_tasks(False))
    m2 = float(np.mean(cpu.get(("m2", None), [0.5]))) * n_m2 * 1.5
    W = 4
    lines.append(f"ablation (13 configurations): {tot / 3600 / W:.2f} h wall")
    for c, s in rows:
        lines.append(f"  {c:9s} {s / 3600 / W:.2f} h")
    lines.append(f"abl2 (R0-ORIG + up to 3 leave-one-out, non-surrogate worst case): {abl2 / 3600 / W:.2f} h")
    lines.append(f"P1 selection + series generation: {(p1 + gen) / 3600 / W:.2f} h")
    lines.append(f"headline (R0, R0-BP96, R0-ORIG at {min(ratio, 11.0):.1f} x R0): {head / 3600 / W:.2f} h")
    lines.append(f"M2 ({n_m2} windows): {m2 / 3600 / W:.2f} h")
    total = (tot + abl2 + p1 + gen + head + m2) / 3600 / W
    lines.append(f"TOTAL projected: {total:.2f} h wall on {W} workers (budget 12 h)")
    txt = "\n".join(lines)
    (RN.RES / "pilot" / "budget.md").write_text("# Pilot runtime projection\n\n```\n" + txt + "\n```\n")
    print(txt)


def CF_is_surr(name):
    from analyses.revision_m1_m2 import configs as CF
    return CF.is_surrogate(name)


if __name__ == "__main__":
    main()
