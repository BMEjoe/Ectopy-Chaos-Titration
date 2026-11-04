"""
Phase 11 runtime budget from the pilot (results/pilot/*.jsonl): mean CPU seconds per task (measured with 4 workers
busy) x confirmatory task count / workers.

    python -m experiments.phase11_robustness.budget11 [--workers 4]
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent


def main(argv=None):
    from experiments.phase11_robustness import run11 as R
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    conf = R.Ctx(False)
    rows = []
    for part in R.TIER1 + R.TIER2:
        p = HERE / "results" / "pilot" / f"{part}.jsonl"
        if not p.exists():
            rows.append({"part": part, "pilot_tasks": 0})
            continue
        t = [json.loads(line)["t_total"] for line in open(p)]
        n_conf = len(R.tasks(conf, part)) if part != "transplant" else 439
        mean = float(np.mean(t))
        rows.append({"part": part, "pilot_tasks": len(t), "mean_s": mean, "conf_tasks": n_conf,
                     "wall_h": mean * n_conf / a.workers / 3600.0})
    tot1 = sum(r.get("wall_h", 0) for r in rows if r["part"] in R.TIER1)
    tot2 = sum(r.get("wall_h", 0) for r in rows if r["part"] in R.TIER2)
    out = {"rows": rows, "tier1_h": tot1, "tier2_h": tot2, "total_h": tot1 + tot2, "workers": a.workers}
    (HERE / "results" / "pilot" / "budget.json").write_text(json.dumps(out, indent=1))
    lines = ["| part | tier | pilot tasks | mean CPU s / task (4 workers busy) | confirmatory tasks | projected wall (4 workers) |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        tier = 1 if r["part"] in R.TIER1 else 2
        if r["pilot_tasks"]:
            lines.append(f"| {r['part']} | {tier} | {r['pilot_tasks']} | {r['mean_s']:.1f} | {r['conf_tasks']} | {r['wall_h']:.2f} h |")
        else:
            lines.append(f"| {r['part']} | {tier} | 0 | – | – | – |")
    lines.append(f"| **total** | | | | | **Tier 1 {tot1:.1f} h + Tier 2 {tot2:.1f} h = {tot1 + tot2:.1f} h** |")
    md = "\n".join(lines)
    (HERE / "results" / "pilot" / "budget.md").write_text(md + "\n")
    print(md)


if __name__ == "__main__":
    main()
