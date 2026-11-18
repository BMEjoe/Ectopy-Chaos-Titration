"""
Revision M1/M2: numbers ledger (EXPLORATORY).  Every number quoted in REPORT.md and manuscript_text.md, with its source
(CSV row or summary.json key), written to results/revision_m1_m2/numbers.md.  `--check` scans REPORT.md and
manuscript_text.md for decimal numbers and lists those that do not appear in the ledger or the CSV tables.

    PYTHONPATH=. python -m analyses.revision_m1_m2.numbers [--check]
"""
from __future__ import annotations

import argparse
import json
import re

import pandas as pd

from analyses.revision_m1_m2 import analyse as AN

R = AN.RES


def f1(x, d=1):
    return "n/a" if x is None else f"{x:.{d}f}"


def ci(v, lo, hi, s=1.0, d=1):
    return f"{s * v:.{d}f} [{s * lo:.{d}f}, {s * hi:.{d}f}]"


def ledger():
    S = json.load(open(AN.SUMMARY))
    L = []
    ab = pd.read_csv(R / "ablation_benchmarks.csv")
    for r in ab.itertuples():
        d = 3 if r.benchmark in ("LOG", "composite") else 1
        sc = 100.0 if r.benchmark == "PM" else 1.0
        L.append((f"abl:{r.config}:{r.benchmark}{':P1' if 'P1' in str(r.switches) else ''}",
                  ci(r.value, r.ci_lo, r.ci_hi, sc, d) + f"; closure {r.closure:+.2f} [{r.closure_ci_lo:+.2f}, "
                  f"{r.closure_ci_hi:+.2f}]; {r.extra if isinstance(r.extra, str) else ''}", "ablation_benchmarks.csv"))
    t1 = pd.read_csv(R / "ablation_table1_by_system.csv")
    for r in t1.itertuples():
        L.append((f"t1:{r.config}:{r.sampling}:{r.system}", f"{r.median_max_noise:.1f} [{r.ci_lo:.1f}, {r.ci_hi:.1f}] "
                  f"(Table 1: {r.table1:.0f})", "ablation_table1_by_system.csv"))
    hb = pd.read_csv(R / "headline_by_reading.csv")
    for r in hb.itertuples():
        sc = 100.0 if r.unit == "rate" else 1.0
        d = 2 if r.unit in ("OR", "ratio") else (2 if r.unit == "fraction" else 1)
        if r.unit == "fraction":
            sc = 100.0
            d = 1
        v = f"{sc * r.value:.{d}f}" if pd.isna(r.ci_lo) else ci(r.value, r.ci_lo, r.ci_hi, sc, d)
        L.append((f"head:{r.reading}:{r.analysis}:{r.quantity}", f"{v} {r.note if isinstance(r.note, str) else ''}",
                  "headline_by_reading.csv"))
    for f in ("m2_influence.csv", "m2_ssr_share.csv", "m2_downslope.csv"):
        if (R / f).exists():
            df = pd.read_csv(R / f)
            for r in df.itertuples():
                key = f"{f[:-4]}:{r.reading}:{r.type}:{r.placement}:k{r.k}" + (f":{r.subset}" if "subset" in df else "")
                vals = {c: getattr(r, c) for c in df.columns if c not in ("reading", "type", "placement", "k", "subset")}
                L.append((key, "; ".join(f"{c}={v:.3g}" if isinstance(v, float) else f"{c}={v}" for c, v in vals.items()
                                         if not (isinstance(v, float) and pd.isna(v))), f))
    m2 = S.get("m2", {}).get("rules", {})
    L.append(("m2:rules", json.dumps(m2, default=float), "summary.json m2.rules"))
    h = S.get("headline", {})
    L.append(("head:agreement", json.dumps(h.get("agreement", {})), "summary.json headline.agreement"))
    L.append(("m2:explore11_agreement", json.dumps(S.get("m2", {}).get("explore11_agreement")), "summary.json"))
    lines = ["# Numbers ledger (EXPLORATORY): every value quoted in REPORT.md, with its source", "",
             "| key | value | source |", "|---|---|---|"]
    lines += [f"| {k} | {v} | {s} |" for k, v, s in L]
    (R / "numbers.md").write_text("\n".join(lines) + "\n")
    return L


def check():
    led = (R / "numbers.md").read_text()
    pool = led
    for f in R.glob("*.csv"):
        pool += f.read_text()
    pool += AN.SUMMARY.read_text()
    missing = {}
    for doc in ("REPORT.md", "manuscript_text.md"):
        p = R / doc
        if not p.exists():
            continue
        txt = p.read_text()
        nums = set(re.findall(r"(?<![\w.])[-+−]?\d+\.\d+", txt))
        miss = []
        for n in sorted(nums):
            v = n.replace("−", "-").lstrip("+")
            cand = {v, v.lstrip("-")}
            try:
                x = float(v)
                for d in (1, 2, 3):
                    cand |= {f"{x:.{d}f}", f"{abs(x):.{d}f}", f"{x / 100:.{d + 2}f}", f"{abs(x) / 100:.{d + 2}f}"}
            except ValueError:
                pass
            if not any(c in pool for c in cand):
                miss.append(n)
        missing[doc] = miss
    print(json.dumps(missing, indent=1))
    return missing


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        check()
    else:
        ledger()
