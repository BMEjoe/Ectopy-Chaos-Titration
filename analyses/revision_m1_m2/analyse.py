"""
Revision M1/M2 analysis (EXPLORATORY; PLAN.md).  Reads results/revision_m1_m2/raw/*.jsonl and the stored Phase 10-12
results; writes results/revision_m1_m2/summary.json and the CSV tables.

    PYTHONPATH=. python -m analyses.revision_m1_m2.analyse --only abl          # ablation + R0-ORIG rule
    PYTHONPATH=. python -m analyses.revision_m1_m2.analyse                     # everything available
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyses.revision_m1_m2 import bench as B
from analyses.revision_m1_m2 import configs as CF
from experiments.phase10_final import analysis10 as A10

HERE = Path(__file__).resolve().parent
RES = HERE.parents[1] / "results" / "revision_m1_m2"
RAW = RES / "raw"
SUMMARY = RES / "summary.json"
P10 = HERE.parents[1] / "experiments" / "phase10_final" / "results"
SEED = 20261117
NB = 10000
TARGET = {"T1": 0.0, "ECO": 50.0, "LOG": 0.75, "PM": 1.0}
TABLE1 = {k: float(v) for k, v in B.S.TABLE1.items()}
BENCH = ("T1", "ECO", "LOG", "PM")


def load(part):
    p = RAW / f"{part}.jsonl"
    out, seen = [], set()
    if not p.exists():
        return out
    for line in open(p):
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("error") or r["id"] in seen:
            continue
        seen.add(r["id"])
        out.append(r)
    return out


def n_errors(part):
    p = RAW / f"{part}.jsonl"
    return sum(1 for line in open(p) if line.strip() and json.loads(line).get("error")) if p.exists() else None


def wil(k, n):
    lo, hi = A10.wilson(int(k), int(n))
    return {"k": int(k), "n": int(n), "rate": k / n if n else None, "ci": [lo, hi]}


def pct(v, lo=2.5, hi=97.5):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return [float(np.percentile(v, lo)), float(np.percentile(v, hi))] if len(v) else [None, None]


# ============================================================================================ M1.2 ablation
def abl_frames():
    R = load("abl") + load("abl2")
    t1 = pd.DataFrame([{"cfg": r["cfg"], "system": r["system"], "seed": r["seed"], "m": r["m"], "p1": r["part"] == "t1p1",
                        "max_noise": r["max_noise"], "base_pos": r["base"]["nonlinear"],
                        "kappa_nl": r["base"].get("kappa_nl"), "d_nl": r["base"].get("d_nl")}
                       for r in R if r["part"] in ("t1", "t1p1")])
    eco = pd.DataFrame([{"cfg": r["cfg"], "seed": r["seed"], "loss": r["loss"], "censored": r["censored"],
                         **{f"pos_{p}": r["levels"][str(p)]["nonlinear"] for p in B.NOISE_ECO},
                         "dopt0": r["levels"]["0"]["d_opt_gt1"], "Fcalc0": r["levels"]["0"]["F_calc_plain"],
                         **{f"dopt_{p}": r["levels"][str(p)]["d_opt_gt1"] for p in B.NOISE_ECO},
                         **{f"Fcalc_{p}": r["levels"][str(p)]["F_calc_plain"] for p in B.NOISE_ECO},
                         **{f"Fcrit_{p}": r["levels"][str(p)]["F_crit"] for p in B.NOISE_ECO},
                         **{f"test_{p}": r["levels"][str(p)]["test"] for p in B.NOISE_ECO}}
                        for r in R if r["part"] == "eco"])
    log = pd.DataFrame([{"cfg": r["cfg"], "seed": r["seed"], "NL": r["NL"], "pos": r["positive"],
                         "kappa_nl": r.get("kappa_nl"), "d_nl": r.get("d_nl"), "r_nl": r.get("r_nl")}
                        for r in R if r["part"] == "log"])
    scan = pd.DataFrame([{"cfg": r["cfg"], "r": r["r"], "NL": r["NL"], "pos": r["positive"]}
                         for r in R if r["part"] == "scan"])
    pm = []
    for r in R:
        if r["part"] != "pm" or "excluded" in r:
            continue
        for L in (500, 2000):
            for s in r[f"seg{L}"]:
                pm.append({"cfg": r["cfg"], "subject": f"{r['db']}:{r['rec']}", "group": r["group"], "L": L, "k": s["k"],
                           "pos": bool(s["raw"]["nonlinear"]), "ectopy": s["ectopy"]})
    return t1, eco, log, scan, pd.DataFrame(pm)


def _le_map():
    R = json.load(open(P10 / "verification" / "titration_verification.json"))["records"]
    return {r["arg"]: r["LE"] for r in R if r["system"] == "logistic"}


def boot_indices(n, seed, B_=NB):
    return np.random.default_rng(seed).integers(0, n, size=(B_, n))


def abl_section():
    t1, eco, log, scan, pm = abl_frames()
    cfgs = [c for c in CF_ORDER() if c in set(t1.cfg)]
    p1 = json.load(open(RES / "p1_choice.json")) if (RES / "p1_choice.json").exists() else None
    seeds = sorted(t1.seed.unique())
    ix_seed = boot_indices(len(seeds), SEED + 1)
    ix_log = boot_indices(100, SEED + 2)
    nsr = [f"nsrdb:{r}" for r in B.records("nsrdb", B.D.DEV_DIR)]
    ix_pm = boot_indices(len(nsr), SEED + 3)
    LE = _le_map()
    out = {"configs": {}, "p1_choice": p1, "n_errors": {"abl": n_errors("abl"), "abl2": n_errors("abl2")}}
    boot = {}

    def t1_matrix(c, p1sampling=False):
        d = t1[(t1.cfg == c) & (~t1.p1)]
        M = {}
        for sy in B.T1_SYSTEMS + B.T1_SECONDARY:
            v = d[d.system == sy].set_index("seed").max_noise.reindex(seeds).values
            if p1sampling and p1 and p1["systems"].get(sy, {}).get("m_opt", 1.0) != 1.0:
                dp = t1[(t1.cfg == c) & t1.p1 & (t1.system == sy)].set_index("seed").max_noise.reindex(seeds).values
                v = dp
            M[sy] = v
        return M

    for c in cfgs:
        res = {"switches": list(CF.get(c).switches)}
        bvals = {}
        # ---------------- T1
        for tag, p1s in (("T1", False), ("T1_P1", True)):
            M = t1_matrix(c, p1s)
            if any(np.isnan(M[sy]).any() for sy in B.T1_SYSTEMS):
                continue
            sysres = {}
            meds_b = []
            for sy in B.T1_SYSTEMS + B.T1_SECONDARY:
                v = M[sy]
                if np.isnan(v).any():
                    continue
                mb = np.median(v[ix_seed], axis=1)
                sysres[sy] = {"median": float(np.median(v)), "ci": pct(mb), "q25": float(np.percentile(v, 25)),
                              "q75": float(np.percentile(v, 75)), "table1": TABLE1[sy],
                              "within15": bool(abs(np.median(v) - TABLE1[sy]) <= 15),
                              "frac_detected_noise_free": float(np.mean(v > 0))}
                if sy in B.T1_SYSTEMS:
                    meds_b.append(np.abs(mb - TABLE1[sy]))
            mad = float(np.mean([abs(sysres[sy]["median"] - TABLE1[sy]) for sy in B.T1_SYSTEMS]))
            mad_b = np.mean(meds_b, axis=0)
            res[tag] = {"systems": sysres, "MAD": mad, "MAD_ci": pct(mad_b),
                        "n_within15": int(sum(sysres[sy]["within15"] for sy in B.T1_SYSTEMS))}
            bvals[tag] = (mad, mad_b)
        # ---------------- ECO
        d = eco[eco.cfg == c].set_index("seed").reindex(seeds)
        if d.loss.notna().all():
            v = d.loss.values
            mb = np.median(v[ix_seed], axis=1)
            curve = {str(p): wil(int(d[f"pos_{p}"].sum()), len(d)) for p in B.NOISE_ECO}
            res["ECO"] = {"median_loss": float(np.median(v)), "ci": pct(mb), "q25": float(np.percentile(v, 25)),
                          "q75": float(np.percentile(v, 75)), "n_censored": int(d.censored.sum()), "curve": curve,
                          "dopt_curve": {str(p): float(d[f"dopt_{p}"].mean()) for p in B.NOISE_ECO},
                          "invF_median": {str(p): float(np.median(1 / d[f"Fcalc_{p}"].values)) for p in B.NOISE_ECO},
                          "Fcrit_median": {str(p): float(np.median(d[f"Fcrit_{p}"].values)) for p in B.NOISE_ECO},
                          "mw_fraction": {str(p): float((d[f"test_{p}"] == "MWW").mean()) for p in B.NOISE_ECO},
                          "pass_rule": {"dopt0_count": int(d.dopt0.sum()), "Fcalc0_median": float(d.Fcalc0.median()),
                                        "pass": bool(d.dopt0.sum() >= 95 and d.Fcalc0.median() >= 100
                                                     and 35 <= np.median(v) <= 65)}}
            bvals["ECO"] = (abs(np.median(v) - 50.0), np.abs(mb - 50.0))
        # ---------------- LOG
        d = log[log.cfg == c].set_index("seed").reindex(range(100))
        if d.NL.notna().all():
            v = d.NL.values
            mb = v[ix_log].mean(axis=1)
            res["LOG"] = {"mean_NL": float(v.mean()), "ci": pct(mb), "sd": float(v.std(ddof=1)),
                          "n_positive": int(d.pos.sum()),
                          "model": d[["kappa_nl", "d_nl"]].astype(str).agg("/".join, axis=1).value_counts().to_dict()
                          if d.kappa_nl.notna().any() else None}
            bvals["LOG"] = (abs(v.mean() - 0.75), np.abs(mb - 0.75))
        # ---------------- scan (secondary)
        d = scan[scan.cfg == c]
        if len(d) == 101:
            le = np.array([LE[r] for r in d.r])
            ch, pe = le > 0.02, le < -0.02
            from scipy import stats
            res["SCAN"] = {"r": list(map(float, d.r)), "NL": list(map(float, d.NL)), "LE": list(map(float, le)),
                           "frac_pos_chaotic": float(np.mean(d.NL.values[ch] > 0)), "n_chaotic": int(ch.sum()),
                           "frac_zero_periodic": float(np.mean(d.NL.values[pe] == 0)), "n_periodic": int(pe.sum()),
                           "spearman": float(stats.spearmanr(d.NL.values[ch], le[ch])[0]),
                           "NL_3.7": float(d[d.r == 3.7].NL.iloc[0])}
        # ---------------- PM
        d = pm[(pm.cfg == c)]
        if len(d):
            d5 = d[d.L == 500]
            per = d5.groupby(["subject", "group"]).pos.mean().reset_index()
            v = per[per.group == "NSR"].set_index("subject").pos.reindex(nsr).values
            if np.isnan(v).any() or d.subject.nunique() < 33:
                v = None
        if len(d) and v is not None:
            mb = v[ix_pm].mean(axis=1)
            gr = A10.group_rates(d5.assign(posf=d5.pos.astype(float)), "pos")
            d2 = d[d.L == 2000]
            gr2 = A10.group_rates(d2, "pos")
            res["PM"] = {"nsr_rate_500": float(v.mean()), "ci": pct(mb), "group_rates_500": gr, "group_rates_2000": gr2}
            bvals["PM"] = (abs(v.mean() - 1.0), np.abs(mb - 1.0))
        out["configs"][c] = res
        boot[c] = bvals
    # ---------------- closures (paired bootstrap; R0 reference)
    if "R0" in boot:
        r0 = boot["R0"]
        for c, bv in boot.items():
            cl = {}
            comp = []
            for b in BENCH:
                if b not in bv or b not in r0:
                    continue
                g0, g0b = r0[b]
                g, gb = bv[b]
                est = (g0 - g) / g0 if g0 > 0 else None
                cb = (g0b - gb) / np.where(g0b > 0, g0b, np.nan)
                cl[b] = {"gap": float(g), "closure": est, "ci": pct(cb)}
                comp.append(cb)
            if len(comp) == 4:
                cbm = np.nanmean(np.vstack(comp), axis=0)
                cl["composite"] = {"closure": float(np.mean([cl[b]["closure"] for b in BENCH])), "ci": pct(cbm)}
            if "T1_P1" in bv:
                g0, g0b = r0["T1"]
                g, gb = bv["T1_P1"]
                cl["T1_with_P1"] = {"gap": float(g), "closure": (g0 - g) / g0, "ci": pct((g0b - gb) / g0b)}
            out["configs"][c]["closure"] = cl
            gc = {}
            res = out["configs"][c]
            if "T1" in res:
                gc["T1"] = res["T1"]["n_within15"] >= 8
            if "ECO" in res:
                gc["ECO"] = 45 <= res["ECO"]["median_loss"] <= 55
            if "LOG" in res:
                gc["LOG"] = abs(res["LOG"]["mean_NL"] - 0.75) <= 0.05
            if "PM" in res:
                gc["PM"] = res["PM"]["nsr_rate_500"] >= 0.90
            out["configs"][c]["gap_closed"] = gc
        # P1 as a protocol switch on R0 (affects T1 only)
        if "T1_P1" in boot["R0"]:
            g0 = boot["R0"]["T1"][0]
            g = boot["R0"]["T1_P1"][0]
            out["P1_on_R0"] = {"T1_MAD_fixed": g0, "T1_MAD_P1": g, "closure_T1": (g0 - g) / g0,
                               "composite_closure": (g0 - g) / g0 / 4}
    out["r0orig_rule"] = r0orig_rule(out)
    return out


def CF_ORDER():
    sw = []
    try:
        s = CF.r0orig_switches()
        sw = ["R0-ORIG"] + [f"R0-ORIG-minus-{x}" for x in s]
    except FileNotFoundError:
        pass
    return ["R0"] + list(CF.SINGLE) + list(CF.COMPOSITES) + sw


def r0orig_rule(abl):
    """PLAN.md 3.5: switches with composite closure > 0; best option in two-option groups."""
    comp = {s: abl["configs"].get(s, {}).get("closure", {}).get("composite", {}).get("closure") for s in CF.SINGLE}
    if any(v is None for v in comp.values()):
        return {"complete": False, "composite": comp}
    elig = [s for s, v in comp.items() if v > 0]
    chosen = []
    for s in elig:
        g = CF.GROUPS.get(s)
        if g is None:
            chosen.append(s)
            continue
        rivals = [x for x in elig if CF.GROUPS.get(x) == g]
        if max(rivals, key=lambda x: comp[x]) == s:
            chosen.append(s)
    order = list(CF.SINGLE)
    chosen = sorted(chosen, key=order.index)
    return {"complete": True, "composite": comp, "eligible": elig, "switches": chosen,
            "rule": "PLAN.md 3.5: every single switch with composite closure > 0; in the groups {S3, S3b} and "
                    "{S5a, S5b} only the option with the larger positive composite closure"}


def write_r0orig(abl):
    rule = abl["r0orig_rule"]
    if not rule.get("complete"):
        raise SystemExit("ablation incomplete: R0-ORIG cannot be built")
    (RES / "r0orig.json").write_text(json.dumps({"switches": rule["switches"], "composite": rule["composite"],
                                                 "eligible": rule["eligible"], "rule": rule["rule"]}, indent=1))
    return rule


def abl_csv(abl):
    rows = []
    for c, res in abl["configs"].items():
        sw = "+".join(res["switches"]) or "-"
        cl = res.get("closure", {})
        vals = {"T1": (res.get("T1", {}).get("MAD"), res.get("T1", {}).get("MAD_ci"), "MAD vs Table 1 (pp)"),
                "ECO": (res.get("ECO", {}).get("median_loss"), res.get("ECO", {}).get("ci"), "median loss point (% noise)"),
                "LOG": (res.get("LOG", {}).get("mean_NL"), res.get("LOG", {}).get("ci"), "mean NL at r = 3.7"),
                "PM": (res.get("PM", {}).get("nsr_rate_500"), res.get("PM", {}).get("ci"), "nsrdb raw 500-beat positive rate")}
        for b in BENCH:
            v, ci, lab = vals[b]
            if v is None:
                continue
            rows.append({"config": c, "switches": sw, "benchmark": b, "metric": lab, "value": v,
                         "ci_lo": ci[0], "ci_hi": ci[1], "target": TARGET[b],
                         "gap": cl.get(b, {}).get("gap"), "closure": cl.get(b, {}).get("closure"),
                         "closure_ci_lo": (cl.get(b, {}).get("ci") or [None, None])[0],
                         "closure_ci_hi": (cl.get(b, {}).get("ci") or [None, None])[1],
                         "gap_closed": res.get("gap_closed", {}).get(b),
                         "extra": (f"{res['T1']['n_within15']}/9 within 15 pp" if b == "T1" else
                                   (f"pass rule {res['ECO']['pass_rule']['pass']}" if b == "ECO" else ""))})
        if "composite" in cl:
            rows.append({"config": c, "switches": sw, "benchmark": "composite", "metric": "mean closure of 4 benchmarks",
                         "value": cl["composite"]["closure"], "ci_lo": cl["composite"]["ci"][0],
                         "ci_hi": cl["composite"]["ci"][1], "target": 1.0, "gap": None,
                         "closure": cl["composite"]["closure"], "closure_ci_lo": cl["composite"]["ci"][0],
                         "closure_ci_hi": cl["composite"]["ci"][1], "gap_closed": None, "extra": ""})
        if "T1_with_P1" in cl:
            rows.append({"config": c, "switches": sw + "+P1", "benchmark": "T1", "metric": "MAD vs Table 1 (pp), P1 sampling",
                         "value": res["T1_P1"]["MAD"], "ci_lo": res["T1_P1"]["MAD_ci"][0], "ci_hi": res["T1_P1"]["MAD_ci"][1],
                         "target": 0.0, "gap": cl["T1_with_P1"]["gap"], "closure": cl["T1_with_P1"]["closure"],
                         "closure_ci_lo": cl["T1_with_P1"]["ci"][0], "closure_ci_hi": cl["T1_with_P1"]["ci"][1],
                         "gap_closed": res["T1_P1"]["n_within15"] >= 8,
                         "extra": f"{res['T1_P1']['n_within15']}/9 within 15 pp"})
    pd.DataFrame(rows).to_csv(RES / "ablation_benchmarks.csv", index=False, float_format="%.4g")
    rows = []
    for c, res in abl["configs"].items():
        for tag in ("T1", "T1_P1"):
            for sy, v in res.get(tag, {}).get("systems", {}).items():
                rows.append({"config": c, "sampling": "Phase 12" if tag == "T1" else "P1", "system": sy,
                             "median_max_noise": v["median"], "ci_lo": v["ci"][0], "ci_hi": v["ci"][1],
                             "q25": v["q25"], "q75": v["q75"], "table1": v["table1"], "within15": v["within15"],
                             "frac_detected_noise_free": v["frac_detected_noise_free"]})
    pd.DataFrame(rows).to_csv(RES / "ablation_table1_by_system.csv", index=False, float_format="%.4g")


# ============================================================================================ main
def save(sections):
    old = json.load(open(SUMMARY)) if SUMMARY.exists() else {}
    old.update(sections)
    SUMMARY.write_text(json.dumps(old, indent=1, default=float))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="all")
    ap.add_argument("--write-r0orig", action="store_true")
    a = ap.parse_args(argv)
    want = None if a.only == "all" else set(a.only.split(","))
    sec = {}
    if want is None or "abl" in want:
        sec["ablation"] = abl_section()
        abl_csv(sec["ablation"])
        if a.write_r0orig:
            print(json.dumps(write_r0orig(sec["ablation"]), indent=1))
    try:
        from analyses.revision_m1_m2 import analyse_more as AM
        sec.update(AM.sections(want))
    except ImportError:
        pass
    save(sec)
    print("ANALYSE DONE", list(sec))


if __name__ == "__main__":
    main()
