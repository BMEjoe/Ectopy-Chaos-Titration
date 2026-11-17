"""
Revision M1/M2: one command that regenerates every result, table and figure (EXPLORATORY; PLAN.md).

    python analyses/revision_m1_m2/run_all.py              # resumes: finished tasks in results/revision_m1_m2/raw are kept
    python analyses/revision_m1_m2/run_all.py --fresh      # moves raw/ aside and recomputes everything (~6 h, 4 workers)
    python analyses/revision_m1_m2/run_all.py --analysis   # tables, figures and summary from the stored raw results only

Steps: regression test (must pass) -> data check (downloads only if the cache is missing, SHA-256 checked) -> P1
selection, series generation, ablation -> R0-ORIG by the PLAN.md 3.5 rule (must equal the committed r0orig.json) ->
R0-ORIG benchmarks, headline units, M2 -> analysis, deviations, tables, figures -> environment record.
Run from the repository root with the pinned environment (requirements-lock.txt) and OMP_NUM_THREADS=1.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
RES = ROOT / "results" / "revision_m1_m2"
PY = sys.executable
ENV = {**os.environ, "PYTHONPATH": str(ROOT), "OMP_NUM_THREADS": "1"}


def sh(*args, check=True):
    print("+", " ".join(args), flush=True)
    return subprocess.run(list(args), cwd=ROOT, env=ENV, check=check)


def data_ok():
    from experiments.phase10_final import data10 as D
    return all((D.DEV_DIR / db / "RECORDS").exists() for db in ("nsrdb", "chfdb", "mitdb")) and \
        all((D.CONF_DIR / db / "RECORDS").exists() for db in ("nsr2db", "chf2db"))


def env_record(t0):
    freeze = subprocess.run([PY, "-m", "pip", "freeze"], capture_output=True, text=True, env=ENV).stdout
    if not freeze.strip():
        freeze = subprocess.run(["uv", "pip", "freeze", "-p", PY], capture_output=True, text=True).stdout
    git = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    cpu = {}
    for p in sorted((RES / "raw").glob("*.jsonl")):
        tot = 0.0
        for line in open(p):
            try:
                tot += json.loads(line).get("t_cpu", 0.0)
            except Exception:                                   # noqa: BLE001
                pass
        cpu[p.stem] = round(tot / 3600, 3)
    rec = {"python": platform.python_version(), "platform": platform.platform(), "git_commit": git("rev-parse", "HEAD"),
           "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"), "plan_commit": git("log", "-1", "--format=%H %cI", "--",
                                                                                     "analyses/revision_m1_m2/PLAN.md"),
           "cpu_hours_by_part": cpu, "cpu_hours_total": round(sum(cpu.values()), 2),
           "this_invocation_wall_minutes": round((time.time() - t0) / 60, 1), "workers": 4, "OMP_NUM_THREADS": 1,
           "seeds": {"new streams": "SeedSequence([20261117, part code, ...])",
                     "Table 1 / ecological": "12000-12099 (Phase 12 generators and noise seeding)",
                     "P1 design": "12900-12904", "logistic r = 3.7": "x0 and noise from seeds 0-99",
                     "logistic scan": "Phase 10 verification seeds", "headline": "Phase 10-12 designs (synthetic 10100-10199)",
                     "M2": "N1 seeds 11000-11199; ventricular random = explore11 draws", "bootstrap": "20261117-20261118"},
           "packages": [ln for ln in freeze.splitlines() if ln.strip()]}
    (RES / "environment.json").write_text(json.dumps(rec, indent=1))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--analysis", action="store_true")
    ap.add_argument("--workers", default="4")
    a = ap.parse_args()
    t0 = time.time()
    if not a.analysis:
        sh(PY, "-m", "pytest", "-q", "analyses/revision_m1_m2/tests")
        if not data_ok():
            sh(PY, "-m", "experiments.phase10_final.data10", "--download", "dev")
            sh(PY, "-m", "experiments.phase10_final.run10", "--phase", "conf", "--download")
        if a.fresh and (RES / "raw").exists():
            shutil.move(str(RES / "raw"), str(RES / f"raw_old_{int(time.time())}"))
        committed = json.load(open(RES / "r0orig.json")) if (RES / "r0orig.json").exists() else None
        sh(PY, "-m", "analyses.revision_m1_m2.run", "--part", "p1sel,gen,abl", "--workers", a.workers)
        sh(PY, "-m", "analyses.revision_m1_m2.analyse", "--only", "abl", "--write-r0orig")
        new = json.load(open(RES / "r0orig.json"))
        if committed is not None and committed["switches"] != new["switches"]:
            raise SystemExit(f"R0-ORIG rule gave {new['switches']}, committed {committed['switches']}: stop and review")
        sh(PY, "-m", "analyses.revision_m1_m2.run", "--part", "abl2,head,m2", "--workers", a.workers)
    sh(PY, "-m", "analyses.revision_m1_m2.analyse")
    sh(PY, "-m", "analyses.revision_m1_m2.deviations")
    sh(PY, "-m", "analyses.revision_m1_m2.tables")
    sh(PY, "-m", "analyses.revision_m1_m2.figures")
    rec = env_record(t0)
    print(json.dumps({k: rec[k] for k in ("git_commit", "cpu_hours_total", "this_invocation_wall_minutes")}, indent=1))
    print("RUN_ALL DONE")


if __name__ == "__main__":
    main()
