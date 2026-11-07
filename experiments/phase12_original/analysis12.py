"""
Phase 12 analysis (PREREGISTRATION.md).  Reads results/conf/*.jsonl (+ Phase 10/11 stored results) and writes
results/conf/analysis/analysis12.json; every table and figure of Phase 12 is built from that file only.

    PYTHONPATH=. python -m experiments.phase12_original.analysis12 [--only B1,B2,...]
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from experiments.phase10_final import analysis10 as A10
from experiments.phase11_robustness import analysis11 as A11
from experiments.phase12_original import run12 as R12
from experiments.phase12_original import systems12 as S

HERE = Path(__file__).resolve().parent
RES = HERE / "results" / "conf"
OUT = RES / "analysis" / "analysis12.json"
GROUPS = ("NSR", "CHF")
SEED = 20261107
READ_ALL = R12.READ_ALL
READ_BP = R12.READ_BP
BP96_REF = {"F_calc": 9387.0, "F_crit": 1.17, "df": [915, 915], "MW_Z": -37.13, "loss_pct": 50.0}

FIDELITY_ROWS = [
    # id, point, source, verdict, decision-relevant
    ("F1", "Model form and term order", "BP96 p. 215, Eq. 1", "matches", False),
    ("F2", "Truncated models (leading r terms)", "BP96 p. 216, Eq. 3; Fig. 1b", "matches", False),
    ("F3", "Constant term counted in r", "BP96 p. 215, Eq. 1; Figs. 1b, 2a", "matches", False),
    ("F4", "Least-squares (Korenberg Gram-Schmidt) estimation", "BP96 p. 215", "matches", False),
    ("F5", "eps^2 = normalised residual variance (eps RMS-like)", "BP96 p. 215, Eq. 2; p. 216", "matches", False),
    ("F6", "C(r) = ln eps(r) + r/N (natural log)", "BP96 p. 216, Eq. 3; PM97 p. 494, Fig. 2 legend", "matches", False),
    ("F7", "N = series length (code: number of regression rows)", "BP96 p. 217, Fig. 2 df; PM97 p. 494", "deviates", True),
    ("F8", "Linear search over the same term range (kappa_lin <= 83)", "BP96 p. 216; Figs. 1b, 2a", "matches", False),
    ("F9", "kappa, d ranges (code kappa <= 6, d <= 3)", "BP96 Figs. 1b, 1c, 2a (d up to 5 and 8)", "not specified", False),
    ("F10", "Nonlinear family d > 1", "BP96 p. 216", "matches", False),
    ("F11", "Decision: d_opt > 1 corroborated by the test", "BP96 p. 216; Figs. 1-2", "matches", False),
    ("F12", "F_calc = eps_lin^2/eps_nl^2 (code: df-adjusted ratio)", "BP96 p. 216, Fig. 1 legend", "deviates", True),
    ("F13", "df = N - r (code: N = rows)", "BP96 p. 217, Fig. 2a/2b; p. 216, Fig. 1", "deviates (via F7)", True),
    ("F14", "1 % level", "BP96 Figs. 1-2 legends; Table 1 legend", "matches", False),
    ("F15", "One-sided F-test", "BP96 p. 216, Fig. 1 legend", "matches", False),
    ("F16", "chi-square Gaussianity check (99 %), Mann-Whitney alternative", "BP96 p. 216", "deviates (absent)", True),
    ("F17", "Mann-Whitney one-sided; ranked quantity", "BP96 p. 216; Figs. 1, 2a", "sidedness matches; quantity not specified", False),
    ("F18", "Surrogate step (FT-type, confirmatory)", "BP96 p. 216", "deviates (absent)", True),
    ("F19", "No preprocessing (code standardises: decision-invariant)", "BP96 p. 217", "matches in effect", False),
    ("F20", "Handling of the first kappa points", "BP96 p. 215, Eq. 2", "not specified", False),
    ("F21", "Minimum series length", "-", "not specified", False),
    ("F22", "Noise titration (NL)", "absent from BP96 and PM97 (Poon & Barahona 2001)", "outside BP96/PM97", False),
]


def load(part):
    p = RES / f"{part}.jsonl"
    out, seen = [], set()
    if not p.exists():
        return out
    for line in open(p):
        r = json.loads(line)
        if r["id"] in seen or r.get("error"):
            continue
        seen.add(r["id"])
        out.append(r)
    return out


def n_errors(part):
    p = RES / f"{part}.jsonl"
    return sum(1 for line in open(p) if json.loads(line).get("error")) if p.exists() else None


def wil(k, n):
    lo, hi = A10.wilson(int(k), int(n))
    return {"k": int(k), "n": int(n), "rate": (k / n) if n else None, "wilson": [float(lo), float(hi)]}


def q(v):
    v = np.asarray(v, float)
    return {"median": float(np.median(v)), "q25": float(np.percentile(v, 25)), "q75": float(np.percentile(v, 75)),
            "n": int(len(v))}


# ------------------------------------------------------------------------------------------------ sections
def sec_sanity():
    p = "experiments/phase11_robustness/results/conf/analysis/analysis11.json"
    root = HERE.parent.parent
    head = subprocess.run(["git", "rev-parse", f"HEAD:{p}"], cwd=root, capture_output=True, text=True).stdout.strip()
    disk = subprocess.run(["git", "hash-object", p], cwd=root, capture_output=True, text=True).stdout.strip()
    return {"analysis11_rerun_identical_to_committed": bool(head and head == disk), "blob_committed": head,
            "blob_rerun": disk}


def sec_fidelity():
    return {"rows": [dict(zip(("id", "point", "source", "verdict", "decision_relevant"), r)) for r in FIDELITY_ROWS],
            "n_decision_relevant_deviations": sum(1 for r in FIDELITY_ROWS if r[4])}


def sec_b1():
    R = load("b1")
    out = {"n_seeds": len(R), "BP96": BP96_REF, "readings": {}}
    for rd in READ_BP:
        seeds, losses, curve = [], [], {}
        for r in R:
            lv = r["levels"]
            p0 = lv["0"][rd]
            seeds.append({"seed": r["seed"], "nonlinear": p0["nonlinear"], "d_opt_gt1": p0["d_opt_gt1"],
                          "F_calc": p0["F_calc_plain"], "F_reading": p0["F"], "F_crit": p0["F_crit"],
                          "df": [p0["df1"], p0["df2"]], "MW_Z": p0["MW_Z"], "kappa_nl": p0["kappa_nl"],
                          "d_nl": p0["d_nl"], "r_lin": p0["r_lin"], "r_nl": p0["r_nl"], "test": p0.get("test")})
            loss = next((int(p) for p in sorted(lv, key=int) if not lv[p][rd]["nonlinear"]), 65)
            losses.append(loss)
            for p in lv:
                curve.setdefault(int(p), []).append((1 / lv[p][rd]["F_calc_plain"], lv[p][rd]["nonlinear"]))
        n_dopt = sum(s["d_opt_gt1"] for s in seeds)
        medF = float(np.median([s["F_calc"] for s in seeds]))
        medL = float(np.median(losses))
        res = {"per_seed": seeds, "n_d_opt_gt1": int(n_dopt), "n_detected_p0": int(sum(s["nonlinear"] for s in seeds)),
               "F_calc": q([s["F_calc"] for s in seeds]), "MW_Z": q([s["MW_Z"] for s in seeds]),
               "F_crit": q([s["F_crit"] for s in seeds]), "loss_point": {**q(losses), "per_seed": losses,
                                                                          "n_censored_gt60": int(sum(x == 65 for x in losses))},
               "kappa_d": {f"{s['kappa_nl']},{s['d_nl']}": sum(1 for t in seeds if (t['kappa_nl'], t['d_nl']) == (s['kappa_nl'], s['d_nl'])) for s in seeds},
               "curve": {str(p): {"inv_F_median": float(np.median([a for a, _ in v])),
                                  "inv_F_q25": float(np.percentile([a for a, _ in v], 25)),
                                  "inv_F_q75": float(np.percentile([a for a, _ in v], 75)),
                                  "detected": int(sum(b for _, b in v)), "n": len(v)} for p, v in sorted(curve.items())}}
        res["rule"] = {"a_dopt_ge19": bool(n_dopt >= 19), "b_medF_ge100": bool(medF >= 100),
                       "c_loss_35_65": bool(35 <= medL <= 65)}
        res["pass"] = bool(all(res["rule"].values()))
        out["readings"][rd] = res
    return out


def sec_b2():
    R = load("b2")
    out = {"systems": {}, "misses": {}}
    for rd in READ_BP:
        miss = []
        for sy in S.SYSTEMS:
            v = [r["max_noise"] for r in R if r["system"] == sy and r["reading"] == rd]
            if not v:
                continue
            qq = q(v)
            t1 = S.TABLE1[sy]
            res = {**qq, "table1": t1, "n_censored_200": int(sum(r["censored"] for r in R if r["system"] == sy and r["reading"] == rd)),
                   "n_not_detected_noise_free": int(sum(x == 0 for x in v)),
                   "consistent": bool(abs(qq["median"] - t1) <= 15), "per_seed": v,
                   "secondary": sy == "mackey_glass_17"}
            out["systems"].setdefault(sy, {})[rd] = res
            if not res["consistent"] and sy != "mackey_glass_17":
                miss.append(sy)
        out["misses"][rd] = miss
    out["note"] = "Series D (Santa Fe D) is a data set and is not generated (PREREGISTRATION 3, B2)"
    return out


def sec_b3():
    p = RES / "b3_status.json"
    return json.load(open(p)) if p.exists() else {"status": "not attempted"}


def sec_b4_static():
    R = load("b4a")
    cells = {}
    for r in R:
        for tr, blk in r["transforms"].items():
            for rd, o in blk.items():
                cells.setdefault((r["phi"], r["n"], tr, rd), []).append(bool(o["nonlinear"]))
    out = {"cells": {}}
    for (phi, n, tr, rd), v in sorted(cells.items()):
        out["cells"][f"{phi}|{n}|{tr}|{rd}"] = wil(sum(v), len(v))
    r0_flags = [k for k, c in out["cells"].items() if k.endswith("|R0") and "|identity|" not in k and c["wilson"][0] > 0.05]
    out["R0_flags_static_transforms"] = bool(r0_flags)
    out["R0_flagged_cells"] = r0_flags
    rem = {}
    for rd in ("R0-S-FT", "R0-S-AAFT", "R0-BP96", "R0-BP96F"):
        for k, c in out["cells"].items():
            if not k.endswith("|R0"):
                continue
            base = k[:-3]
            o = out["cells"].get(f"{base}|{rd}")
            if o and c["k"] >= 50:
                rem.setdefault(rd, {})[base] = bool(o["k"] <= 5)
    out["removes"] = rem
    return out


def sec_b4_verify():
    from experiments.phase10_final import verify_titration as VT
    R = load("b4v")
    out = {}
    a11 = json.load(open(HERE.parent / "phase11_robustness" / "results" / "conf" / "analysis" / "analysis11.json"))
    v0 = a11["A5"]["verification"]["R0"]
    out["R0"] = {k: v0[k] for k in ("V1", "V2", "V3", "V4", "V5", "V6", "VERIFIED")}
    out["R0"]["source"] = "Phase 11 analysis11.json A5.verification.R0"
    for rd in ("R0-BP96F", "R0-BP96", "R0-S-FT", "R0-S-AAFT"):
        rr = [r for r in R if r["reading"] == rd]
        try:
            v = VT.evaluate(rr)
            out[rd] = {k: (v[k]["pass"] if isinstance(v[k], dict) else v[k]) for k in ("V1", "V2", "V3", "V4", "V5", "V6", "VERIFIED")}
            out[rd]["detail"] = {k: {kk: vv for kk, vv in v[k].items() if kk != "scan"} for k in ("V1", "V2", "V3", "V4", "V5", "V6")}
        except Exception as exc:                                # noqa: BLE001
            out[rd] = {"error": f"{type(exc).__name__}: {exc} ({len(rr)} records)"}
    return out


def sec_b4_synth():
    R = load("b4s")
    out = {}
    for c in R12.B4S_CONDS:
        rr = [r for r in R if r["cond"] == c]
        out[c] = {rd: wil(sum(r["raw"][rd]["nonlinear"] for r in rr), len(rr)) for rd in READ_ALL}
    return out


def sec_b4_transplant():
    R = load("b4t")
    rows = []
    for r in R:
        wid = r["id"].split("|", 1)[1]
        row = {"wid": wid, "subject": r["subject"], "group": r["group"]}
        for v, blk in r["variants"].items():
            for rd, o in blk.items():
                row[f"{v}|{rd}"] = bool(o["nonlinear"])
        rows.append(row)
    df = pd.DataFrame(rows)
    out = {"n_windows": len(df), "n_subjects": int(df.subject.nunique()) if len(df) else 0, "readings": {}}
    for rd in READ_ALL:
        d = df.rename(columns={f"{v}|{rd}": v for v in R12.TRANSPLANT_VARIANTS})
        res = {"rates": {v: float(d[v].mean()) for v in R12.TRANSPLANT_VARIANTS}}
        for v in ("V1", "V5", "V6"):
            res[f"{v}_minus_V0"] = A11.paired_delta(d, "V0", v)
        out["readings"][rd] = res
    # agreement of the recomputed R0 with the Phase 11 stored transplant decisions
    P11 = {json.loads(line)["id"].split("|", 1)[1]: json.loads(line)
           for line in open(HERE.parent / "phase11_robustness" / "results" / "conf" / "transplant.jsonl")}
    agree = tot = 0
    for _, row in df.iterrows():
        s = P11.get(row["wid"])
        if not s:
            continue
        for v in R12.TRANSPLANT_VARIANTS:
            tot += 1
            agree += int(bool(s[v]["positive"]) == bool(row[f"{v}|R0"]))
    out["R0_agreement_with_phase11"] = {"agree": agree, "n": tot}
    return out


def window_frames():
    df10 = A11.p10_real("conf")
    base = df10.set_index("wid")[["subject", "group", "burden", "x", "chf"]]
    frames = {"R0": df10[["wid", "subject", "group", "burden", "x", "chf", "TIT", "TITm", "TITm_analysable", "TITe",
                          "TITe_eligible"]].rename(columns={"TIT": "raw", "TITm": "masked", "TITm_analysable": "masked_ok",
                                                            "TITe": "edited", "TITe_eligible": "edited_ok"})}
    R = load("b4w")
    recomputed = {}
    for rd in READ_ALL:
        rows = []
        for r in R:
            wid = r["id"].split("|", 1)[1]
            e = r.get("edited", {}).get(rd, {})
            rows.append({"wid": wid, "raw": bool(r["raw"][rd]["nonlinear"]), "masked": bool(r["masked"][rd]["nonlinear"]),
                         "masked_ok": bool(r["masked"][rd].get("analysable")),
                         "edited": bool(e.get("nonlinear", False)),
                         "edited_ok": bool(r["edited_eligible"] and e.get("analysable", False))})
        f = pd.DataFrame(rows).set_index("wid").join(base, how="inner").reset_index()
        if rd == "R0":
            recomputed["R0"] = f
        else:
            frames[rd] = f
    return frames, recomputed.get("R0")


def claims(f):
    m = {"n_windows": len(f)}
    g = A10.gee(f.assign(TITr=f.raw.astype(bool)), "TITr", ["x", "chf"])
    m["i_gee"] = {k: g.get(k) for k in ("estimated", "working_correlation", "n_pos")}
    if g.get("estimated"):
        c = g["coef"]["x"]
        m["i_gee"].update(OR=c["OR"], ci95=c["OR_ci95"])
        m["i"], m["i_reversed"] = bool(c["OR_ci95"][0] > 1), bool(c["OR_ci95"][1] < 1)
    else:
        m["i"] = None
    gr = A10.group_rates(f, "raw")
    m["ii_windows"] = gr.get("CHF_minus_NSR")
    m["ii"] = bool(m["ii_windows"]["ci95"][0] > 0) if m["ii_windows"] else None
    post = {arm: A10.group_rates(f, arm, f"{arm}_ok").get("CHF_minus_NSR") for arm in ("masked", "edited")}
    m["iii_windows"] = post
    m["iii"] = bool(all(post[a] and post[a]["ci95"][0] <= 0 <= post[a]["ci95"][1] for a in post))
    m["iii_reversed"] = bool(any(post[a] and not (post[a]["ci95"][0] <= 0 <= post[a]["ci95"][1]) for a in post))
    m["rates"] = {"raw": float(f.raw.mean()), "masked_analysable": float(f[f.masked_ok].masked.mean()),
                  "edited_eligible": float(f[f.edited_ok].edited.mean()), "masked_analysable_frac": float(f.masked_ok.mean())}
    return m


def sec_b4_matrix(sections):
    frames, r0re = window_frames()
    out = {"readings": {}}
    if r0re is not None:
        st = frames["R0"].set_index("wid")
        j = r0re.set_index("wid")[["raw", "masked", "edited", "masked_ok", "edited_ok"]].join(
            st[["raw", "masked", "edited"]], rsuffix="_p10", how="inner")
        out["R0_recomputed_vs_phase10"] = {"n": len(j), "raw_agree": int((j.raw == j.raw_p10).sum()),
                                          "masked_agree": int((j.masked == j.masked_p10).sum()),
                                          "edited_agree": int((j.edited == j.edited_p10).sum())}
    ver = sections.get("B4_verify", {})
    syn = sections.get("B4_synth", {})
    tp = sections.get("B4_transplant", {}).get("readings", {})
    st = sections.get("B4_static", {})
    for rd in READ_ALL:
        m = claims(frames[rd])
        e3, n1 = syn.get("E3_couplets_5pct", {}).get(rd), syn.get("N1_linear_rr", {}).get(rd)
        m["iv"] = bool(e3 and n1 and e3["k"] >= 50 and n1["k"] <= 5) if (e3 and n1) else None
        d = tp.get(rd, {}).get("V1_minus_V0", {}).get("pooled")
        m["v"] = bool(d and d["ci95"][0] > 0) if d else None
        m["v_delta"] = d
        m["verified"] = ver.get(rd, {}).get("VERIFIED")
        m["against_H"] = bool(m["verified"] and (m.get("i") is False or m.get("iii") is False or m.get("v") is False))
        n5 = syn.get("N5_linear_rr_warped", {}).get(rd)
        m["N5_synthetic"] = n5
        rem = st.get("removes", {}).get(rd, {})
        m["static_cells_removed"] = {k: v for k, v in rem.items()}
        m["removes_N5"] = bool(n5 and n5["k"] <= 5 and all(v for k, v in rem.items() if "|N5" in k)) if rd != "R0" else False
        m["removes_ectopy"] = bool((d and d["ci95"][0] <= 0) or (e3 and e3["k"] < 50)) if rd != "R0" else False
        out["readings"][rd] = m
    # registered N5 sentence (PREREGISTRATION B4.3 (3))
    rr = out["readings"]
    ft = rr["R0-S-FT"]
    aa = rr["R0-S-AAFT"]
    if not ((ft["removes_N5"] and ft["verified"]) or (aa["removes_N5"] and aa["verified"])):
        n5 = "R0, like titration as published in 2001, flags static transforms, and so does BP96's surrogate step (FT and AAFT)."
        case = "no verified surrogate reading removes N5"
    elif aa["removes_N5"] and aa["verified"] and not ft["removes_N5"]:
        n5 = ("R0 flags static transforms; amplitude-adjusted surrogates would remove them, phase-randomised surrogates "
              "(BP96's description) would not.")
        case = "AAFT removes, FT does not"
    else:
        n5 = None
        case = "other pattern: described as observed"
    out["N5_sentence"] = {"case": case, "text": n5}
    out["any_headline_change_verified"] = {rd: bool(rr[rd]["verified"] and any(rr[rd].get(k) != rr["R0"].get(k)
                                                                             for k in ("i", "ii", "iii", "iv", "v")))
                                           for rd in READ_ALL if rd != "R0"}
    return out


# ------------------------------------------------------------------------------------------------ PM97
def pm97_frames():
    R = load("pm97")
    subj, frames = [], {500: [], 2000: []}
    for r in R:
        sid = f"{r['db']}:{r['rec']}"
        if "excluded" in r:
            subj.append({"subject": sid, "db": r["db"], "rec": r["rec"], "excluded": r["excluded"]})
            continue
        s500 = r["seg500"]
        subj.append({"subject": sid, "db": r["db"], "rec": r["rec"], "group": r["group"], "rule": r["rule"],
                     "start_beat": r["start_beat"], "end_beat": r["end_beat"], "t_start_h": r["t_start"] / 3600,
                     "t_end_h": r["t_end"] / 3600, "base_time": r["base_time"], "fs": r["fs"],
                     "n_nonnormal": r["n_nonnormal_block"], "n_out_of_range": r["n_out_of_range_block"],
                     "sd_segment_means": float(np.std([s["mean_rr"] for s in s500], ddof=1))})
        for L in (500, 2000):
            for s in r[f"seg{L}"]:
                row = {"subject": sid, "group": r["group"], "k": s["k"], "ectopy": s["ectopy"],
                       "n_out": s["n_out_of_range"], "coef80": s.get("coef80")}
                for arm in ("raw", "rawf", "masked", "edited", "wu"):
                    blk = s.get(arm)
                    elig = {"edited": s.get("edited_eligible"), "wu": s.get("wu_eligible")}.get(arm, True)
                    for rd in READ_BP:
                        o = (blk or {}).get(rd, {})
                        row[f"{arm}|{rd}"] = bool(o.get("nonlinear", False))
                        row[f"{arm}|{rd}|ok"] = bool(elig and o.get("analysable", False))
                frames[L].append(row)
    return pd.DataFrame(subj), {L: pd.DataFrame(v) for L, v in frames.items()}


def pm97_stats(df, col, ok):
    d = df[df[ok]]
    per = d.groupby(["subject", "group"])[col].mean().reset_index()
    gr = A10.group_rates(d, col)
    out = {"group_rates": gr}
    v = {g: per[per.group == g][col].values for g in GROUPS}
    for g in GROUPS:
        if len(v[g]):
            out[g] = {"n_subjects": int(len(v[g])), "mean": float(v[g].mean()), "median": float(np.median(v[g]))}
    if len(v["NSR"]) and len(v["CHF"]):
        out["mw_p"] = float(stats.mannwhitneyu(v["CHF"], v["NSR"], alternative="two-sided").pvalue)
        hm = float(v["NSR"].mean())
        drop = hm - v["CHF"]
        out["chf_20_80pp_below_nsr_mean"] = int(np.sum((drop >= 0.20) & (drop <= 0.80)))
        out["chf_relative_reduction_20_80pct"] = int(np.sum(((drop / hm) >= 0.20) & ((drop / hm) <= 0.80))) if hm > 0 else None
        out["n_chf"] = int(len(v["CHF"]))
        cmn = gr.get("CHF_minus_NSR")
        out["pattern"] = {"a_nsr_mean_ge_0.90": bool(hm >= 0.90),
                          "b_diff_ci_below_0": bool(cmn and cmn["ci95"][1] < 0),
                          "c_ge_2_3_chf_20_80pp": bool(out["chf_20_80pp_below_nsr_mean"] >= 2 / 3 * len(v["CHF"]))}
        out["pattern"]["reproduced"] = bool(all(out["pattern"].values()))
    out["per_subject"] = {s: float(x) for s, x in zip(per.subject, per[col])}
    return out


def sec_pm97():
    subj, F = pm97_frames()
    out = {"subjects": subj.to_dict("records"), "n_excluded": int(subj.get("excluded", pd.Series(dtype=object)).notna().sum())
           if "excluded" in subj else 0}
    inc = subj[subj.get("excluded").isna()] if "excluded" in subj else subj
    # PM97-like selection (registered rank rule)
    sel = []
    for db, k in (("nsrdb", 8), ("chfdb", 11)):
        d = inc[inc.db == db].copy()
        d["r1"] = d.n_nonnormal.rank(method="average")
        d["r2"] = d.sd_segment_means.rank(method="average")
        d["rs"] = d.r1 + d.r2
        d = d.sort_values(["rs", "rec"])
        sel += list(d.subject.iloc[:k])
    out["pm97_like_selection"] = sel
    out["cells"] = {}
    for L in (500, 2000):
        df = F[L]
        for subset, dd in (("all", df), ("selected", df[df.subject.isin(sel)])):
            for arm in ("raw", "rawf", "masked", "edited", "wu"):
                for rd in READ_BP:
                    key = f"{L}|{subset}|{arm}|{rd}"
                    out["cells"][key] = pm97_stats(dd, f"{arm}|{rd}", f"{arm}|{rd}|ok")
    prim = out["cells"]["500|all|raw|R0"]
    out["primary"] = {"pattern": prim.get("pattern"), "NSR": prim.get("NSR"), "CHF": prim.get("CHF"),
                      "diff": prim["group_rates"].get("CHF_minus_NSR"), "mw_p": prim.get("mw_p")}
    after = {arm: out["cells"][f"500|all|{arm}|R0"]["group_rates"].get("CHF_minus_NSR") for arm in ("masked", "edited", "wu")}
    out["after_handling"] = after
    out["against_claims_decrease_survives_handling"] = bool(all(a and a["ci95"][1] < 0 for a in after.values()))
    # within-subject ectopy association (conditional logit), 500-beat raw R0, all subjects
    d5 = F[500]
    cl = d5.assign(pos=d5["raw|R0"].astype(float), x=np.log2(1 + d5.ectopy.astype(float)))
    out["clogit_ectopy"] = A11.clogit(cl, "pos", ["x"])
    out["ectopy_by_group"] = {g: q(d5[d5.group == g].groupby("subject").ectopy.sum().values) for g in GROUPS}
    # DESCRIPTIVE / POST HOC (not registered; added after the PM97 outputs were seen, labelled as such)
    out["descriptive_handling"] = {
        "label": "descriptive, not registered (added after the outputs were seen)",
        "analysable_fraction_500": {g: {arm: float(d5[d5.group == g][f"{arm}|R0|ok"].mean())
                                        for arm in ("raw", "rawf", "masked", "edited", "wu")} for g in GROUPS},
        "segments_with_ectopy_fraction_500": {g: float((d5[d5.group == g].ectopy > 0).mean()) for g in GROUPS},
        "ectopic_beats_per_segment_500": {g: q(d5[d5.group == g].ectopy.values) for g in GROUPS}}
    clean = d5[(d5.ectopy == 0) & (d5.n_out == 0)]
    out["exploratory_clean_segments"] = {
        "label": "EXPLORATORY / POST HOC (not registered): 500-beat segments with no non-normal beat and no out-of-range "
                 "interval, so raw = every handled arm; raw decision",
        "n_segments": {g: int((clean.group == g).sum()) for g in GROUPS},
        "n_subjects": {g: int(clean[clean.group == g].subject.nunique()) for g in GROUPS},
        **{rd: A10.group_rates(clean, f"raw|{rd}") for rd in ("R0", "R0-BP96")}}
    # Fig. 3 reading: mean SD of the first 80 nonlinear coefficients across positive raw 500-beat segments
    vals = {}
    for sid, g in d5[d5["raw|R0"] & d5.coef80.notna()].groupby("subject"):
        C = np.array(list(g.coef80.values), float)
        if len(C) >= 2:
            vals[sid] = (g.group.iloc[0], float(np.mean(np.std(C, axis=0, ddof=1))))
    fig3 = {"per_subject": {k: v[1] for k, v in vals.items()}}
    for gname in GROUPS:
        x = np.array([v[1] for v in vals.values() if v[0] == gname])
        if len(x):
            fig3[gname] = {"n": int(len(x)), "mean": float(x.mean()), "se": float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else None}
    xs = {gname: np.array([v[1] for v in vals.values() if v[0] == gname]) for gname in GROUPS}
    if all(len(x) > 1 for x in xs.values()):
        fig3["CHF_minus_NSR"] = {"diff": float(xs["CHF"].mean() - xs["NSR"].mean()),
                                 "ci95": A10.ci(A10.boot_diff(xs["CHF"], xs["NSR"])),
                                 "mw_p": float(stats.mannwhitneyu(xs["CHF"], xs["NSR"], alternative="two-sided").pvalue)}
    fig3["reading"] = "full kappa = 7, d = 3 model, first 80 nonlinear coefficients, positive raw R0 500-beat segments"
    out["fig3"] = fig3
    # secondary: titration NL (R0, raw 500-beat)
    NL = load("pm97nl")
    rows = []
    for r in NL:
        if "excluded" in r:
            continue
        sid = f"{r['db']}:{r['rec']}"
        g = "NSR" if r["db"] == "nsrdb" else "CHF"
        for s in r["segments"]:
            rows.append({"subject": sid, "group": g, "NL": s["NL"], "pos": s["positive"]})
    nd = pd.DataFrame(rows)
    if len(nd):
        out["titration"] = {"positive_rate": A10.group_rates(nd, "pos"),
                            "mean_NL_positive": {g: q(nd[(nd.group == g) & nd.pos].NL.values) if ((nd.group == g) & nd.pos).any() else None for g in GROUPS},
                            "agreement_with_noise_free_decision": None}
        dec = d5.set_index(["subject", "k"])["raw|R0"]
        agr = tot = 0
        for r in NL:
            if "excluded" in r:
                continue
            sid = f"{r['db']}:{r['rec']}"
            for s in r["segments"]:
                if (sid, s["k"]) in dec.index:
                    tot += 1
                    agr += int(bool(dec[(sid, s["k"])]) == bool(s["positive"]))
        out["titration"]["agreement_with_noise_free_decision"] = {"agree": agr, "n": tot}
    return out


SECTIONS = {"sanity": sec_sanity, "fidelity": sec_fidelity, "B1": sec_b1, "B2": sec_b2, "B3": sec_b3,
            "B4_static": sec_b4_static, "B4_verify": sec_b4_verify, "B4_synth": sec_b4_synth,
            "B4_transplant": sec_b4_transplant, "B4_matrix": None, "PM97": sec_pm97}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    a = ap.parse_args(argv)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res = json.load(open(OUT)) if OUT.exists() else {}
    names = list(SECTIONS) if a.only is None else a.only.split(",")
    for n in names:
        try:
            res[n] = sec_b4_matrix(res) if n == "B4_matrix" else SECTIONS[n]()
            print(n, "ok", flush=True)
        except Exception as exc:                                # noqa: BLE001
            import traceback
            res[n] = {"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc(limit=6)}
            print(n, "ERROR", exc, flush=True)
    res["errors_per_part"] = {p: n_errors(p) for p in R12.PART_CODE}
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
