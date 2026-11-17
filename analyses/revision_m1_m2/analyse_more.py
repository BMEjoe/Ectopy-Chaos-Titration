"""
Revision M1/M2 analysis, part 2 (EXPLORATORY): headline results by reading (PLAN.md 4), M2 diagnostics (6) and the
claims check (5).  Called from analyse.py.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analyses.revision_m1_m2 import analyse as AN
from analyses.revision_m1_m2 import m2gen as M
from experiments.phase10_final import analysis10 as A10
from experiments.phase11_robustness import analysis11 as A11
from experiments.phase12_original import analysis12 as A12

ROOT = Path(__file__).resolve().parents[2]
P10 = ROOT / "experiments" / "phase10_final" / "results" / "conf"
P11 = ROOT / "experiments" / "phase11_robustness" / "results" / "conf"
P12 = ROOT / "experiments" / "phase12_original" / "results" / "conf"
GROUPS = ("NSR", "CHF")
HEADS = ("R0", "R0-BP96", "R0-ORIG")
SEED = 20261118


def jl(path):
    return [json.loads(line) for line in open(path) if line.strip()]


def ok_flag(o, arm):
    if arm == "raw":
        return True
    if arm == "masked":
        return bool(o.get("analysable", False))
    return bool(o.get("eligible", False) and o.get("analysable", False))


# ============================================================================================ headline frames
def head_records():
    return AN.load("head")


def win_frame(R):
    base = A11.p10_real("conf").set_index("wid")[["subject", "group", "burden", "x", "chf", "n_masked"]]
    rows = []
    for r in R:
        if r["part"] != "win":
            continue
        wid = r["id"].split("|", 2)[2]
        row = {"cfg": r["cfg"], "wid": wid}
        for a in ("raw", "masked", "edited", "wu"):
            row[a] = bool(r[a]["nonlinear"])
            row[f"{a}_ok"] = ok_flag(r[a], a)
        rows.append(row)
    return pd.DataFrame(rows).join(base, on="wid", how="inner")


def tr_frame(R):
    rows = []
    for r in R:
        if r["part"] != "tr":
            continue
        rows.append({"cfg": r["cfg"], "wid": r["id"].split("|", 2)[2], "subject": r["subject"], "group": r["group"],
                     **{v: bool(r[v]["nonlinear"]) for v in ("V0", "V1", "V5", "V6")}})
    return pd.DataFrame(rows)


def seg_frame(R):
    rows = []
    for r in R:
        if r["part"] != "seg":
            continue
        for s in r["segments"]:
            if s["gap"]:
                continue
            row = {"cfg": r["cfg"], "subject": r["subject"], "group": r["group"], "seg": s["seg"], "burden": s["burden"],
                   "n_masked": s["n_masked"]}
            for a in ("raw", "masked", "edited", "wu"):
                row[a] = bool(s[a]["nonlinear"])
                row[f"{a}_ok"] = ok_flag(s[a], a)
            rows.append(row)
    return pd.DataFrame(rows)


def synth_frame(R):
    rows = []
    for r in R:
        if r["part"] != "synth":
            continue
        cid = r["id"].split("|")[2]
        rows.append({"cfg": r["cfg"], "cond": cid, "seed": r["seed"], "kind": r["kind"],
                     **{a: bool(r[a]["nonlinear"]) for a in ("raw", "masked", "edited")},
                     "masked_ok": ok_flag(r["masked"], "masked"), "edited_ok": ok_flag(r["edited"], "edited")})
    return pd.DataFrame(rows)


def pm97_frame(R):
    rows = []
    for r in R:
        if r["part"] != "pm97" or "excluded" in r:
            continue
        sid = f"{r['db']}:{r['rec']}"
        for L in (500, 2000):
            for s in r[f"seg{L}"]:
                row = {"cfg": r["cfg"], "subject": sid, "group": r["group"], "L": L, "k": s["k"], "ectopy": s["ectopy"],
                       "n_out": s["n_out"]}
                for a in ("raw", "rawf", "masked", "edited", "wu"):
                    row[a] = bool(s[a]["nonlinear"])
                    row[f"{a}_ok"] = ok_flag(s[a], a) if a not in ("raw", "rawf") else True
                rows.append(row)
    return pd.DataFrame(rows)


# ============================================================================================ agreement with stored
def agreement(W, TR, SG, SY, PM):
    out = {}
    st = {r["id"]: r for r in jl(P10 / "real.jsonl") if not r.get("error")}
    w = W[W.cfg == "R0"]
    out["windows_R0_vs_phase10"] = {a: {"agree": int(sum(bool(st[i][f"TIT_{a}"]["positive"]) == v
                                                         for i, v in zip(w.wid, w[a]))), "n": len(w)}
                                    for a in ("raw", "masked", "edited")}
    b4w = {r["id"].split("|", 1)[1]: r for r in jl(P12 / "b4w.jsonl") if not r.get("error")}
    w = W[W.cfg == "R0-BP96"]
    out["windows_R0BP96_vs_phase12"] = {"raw": {"agree": int(sum(bool(b4w[i]["raw"]["R0-BP96"]["nonlinear"]) == v
                                                                 for i, v in zip(w.wid, w.raw) if i in b4w)),
                                                "n": int(sum(i in b4w for i in w.wid))}}
    st = {r["id"].split("|", 1)[1]: r for r in jl(P11 / "transplant.jsonl") if not r.get("error")}
    t = TR[TR.cfg == "R0"]
    out["transplant_R0_vs_phase11"] = {"agree": int(sum(bool(st[i][v]["positive"]) == bool(row[v])
                                                        for i, row in zip(t.wid, t.to_dict("records"))
                                                        for v in ("V0", "V1", "V5", "V6"))), "n": 4 * len(t)}
    s10 = {}
    for r in jl(P10 / "seg.jsonl"):
        if r.get("error"):
            continue
        for s in r["segments"]:
            if not s["gap"]:
                s10[(r["subject"], s["seg"])] = s
    g = SG[SG.cfg == "R0"]
    out["segments_R0_vs_phase10"] = {a: {"agree": int(sum(bool(s10[(sj, sg)][f"TIT_{a}"].get("positive", False)) == v
                                                          for sj, sg, v in zip(g.subject, g.seg, g[a]))), "n": len(g)}
                                     for a in ("raw", "masked", "edited", "wu")}
    s10 = {r["id"]: r for r in jl(P10 / "synth.jsonl") if not r.get("error")}
    y = SY[SY.cfg == "R0"]
    out["synth_R0_vs_phase10"] = {a: {"agree": int(sum(bool(s10[f"synth|{c}|512|{sd}"][f"TIT_{a}"]["positive"]) == v
                                                       for c, sd, v in zip(y.cond, y.seed, y[a]))), "n": len(y)}
                                  for a in ("raw", "masked", "edited")}
    p12 = {}
    for r in jl(P12 / "pm97.jsonl"):
        if r.get("error") or "excluded" in r:
            continue
        for L in (500, 2000):
            for s in r[f"seg{L}"]:
                p12[(f"{r['db']}:{r['rec']}", L, s["k"])] = s
    for rd in ("R0", "R0-BP96"):
        p = PM[PM.cfg == rd]
        out[f"pm97_{rd}_vs_phase12"] = {a: {"agree": int(sum(bool(p12[(sj, L, k)][a][rd]["nonlinear"]) == v
                                                             for sj, L, k, v in zip(p.subject, p.L, p.k, p[a]))),
                                            "n": len(p)} for a in ("raw", "rawf", "masked")}
    return out


# ============================================================================================ headline statistics
def burden_rows(d, col="raw"):
    out = []
    for lo, hi in A10.BURDEN_BINS:
        x = d[(d.burden >= lo) & (d.burden <= hi)]
        lab = f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 10 ** 9 else f">={lo}")
        out.append({"burden": lab, "n": len(x), "n_subjects": int(x.subject.nunique()), **AN.wil(int(x[col].sum()), len(x))})
    return out


def af_emp(d):
    dd = d.assign(_y=d.raw.astype(float), _n=1.0, _y0=(d.raw & (d.burden == 0)).astype(float),
                  _n0=(d.burden == 0).astype(float))
    est, ci = A11.ratio_boot(dd, ["_y", "_n", "_y0", "_n0"],
                             lambda s: 1 - (s[:, 2] / s[:, 3]) / (s[:, 0] / s[:, 1]), by_group=True, seed=SEED + 11)
    return {"AF": est, "ci": ci, "rate_all": float(d.raw.mean()), "rate_burden0": float(d[d.burden == 0].raw.mean())}


def af_gee(d, B_=1000):
    import statsmodels.api as sm

    def fit(dd):
        X = sm.add_constant(dd[["x", "chf"]].astype(float).values)
        y = dd.raw.astype(float).values
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = sm.GEE(y, X, groups=dd.subject.values, family=sm.families.Binomial(),
                       cov_struct=sm.cov_struct.Exchangeable()).fit()
            if not (np.all(np.isfinite(m.params)) and np.all(np.isfinite(m.bse))):
                m = sm.GEE(y, X, groups=dd.subject.values, family=sm.families.Binomial(),
                           cov_struct=sm.cov_struct.Independence()).fit()
        b = m.params
        p0 = 1 / (1 + np.exp(-(b[0] + b[2] * dd.chf.values)))
        return 1 - p0.mean() / y.mean()
    est = float(fit(d))
    rng = np.random.default_rng(SEED + 12)
    subs = {g: d[d.group == g].subject.unique() for g in GROUPS}
    by = {s: d[d.subject == s] for s in d.subject.unique()}
    bs = []
    for b in range(B_):
        parts = []
        for g in GROUPS:
            pick = rng.choice(subs[g], len(subs[g]), replace=True)
            parts += [by[s].assign(subject=f"{s}#{i}") for i, s in enumerate(pick)]
        try:
            bs.append(fit(pd.concat(parts)))
        except Exception:                                       # noqa: BLE001
            pass
    return {"AF": est, "ci": AN.pct(bs), "n_boot_ok": len(bs)}


def clogit_block(d):
    cw, _ = A11.covar_frames("conf")
    df = d.merge(cw.drop(columns=["subject", "group", "j"]), on="wid", how="inner")
    dd, nd = A11._prep_cov(df.assign(TITr=df.raw.astype(bool)))
    return {"n": len(dd), "n_dropped": nd, "cov": A11.clogit(dd, "TITr", ["x"] + A11.COVS),
            "x_only": A11.clogit(dd, "TITr", ["x"]),
            "gee_cov_chf": A10.gee(dd, "TITr", ["x", "chf"] + A11.COVS)}


def head_section(want_af_gee=True):
    R = head_records()
    W, TR, SG, SY, PM = win_frame(R), tr_frame(R), seg_frame(R), synth_frame(R), pm97_frame(R)
    out = {"n_errors": AN.n_errors("head"), "agreement": agreement(W, TR, SG, SY, PM), "readings": {}}
    for rd in HEADS:
        res = {}
        w = W[W.cfg == rd]
        if len(w):
            res["dose_burden"] = burden_rows(w)
            res["gee"] = A10.gee(w.assign(TIT=w.raw), "TIT", ["x", "chf"])
            res["gee_chf_unadjusted"] = A10.gee(w.assign(TIT=w.raw), "TIT", ["chf"])
            res["clogit"] = clogit_block(w)
            res["AF_emp"] = af_emp(w)
            if want_af_gee:
                res["AF_gee"] = af_gee(w)
            res["ratio_1_vs_0"] = (float(w[w.burden == 1].raw.mean() / w[w.burden == 0].raw.mean()))
            res["win_groups"] = {a: A10.group_rates(w, a, None if a == "raw" else f"{a}_ok")
                                 for a in ("raw", "masked", "edited", "wu")}
            res["win_paired_change"] = {a: A10.paired_change(w, "raw", a, f"{a}_ok") for a in ("masked", "edited", "wu")}
            res["win_ok_frac"] = {a: {g: float(w[w.group == g][f"{a}_ok"].mean()) for g in GROUPS}
                                  for a in ("masked", "edited", "wu")}
        t = TR[TR.cfg == rd]
        if len(t):
            res["transplant"] = {"rates": {v: float(t[v].mean()) for v in ("V0", "V1", "V5", "V6")},
                                 **{f"{v}_minus_V0": A11.paired_delta(t, "V0", v) for v in ("V1", "V5", "V6")}}
        s = SG[SG.cfg == rd]
        if len(s):
            res["seg_groups"] = {a: A10.group_rates(s, a, None if a == "raw" else f"{a}_ok")
                                 for a in ("raw", "masked", "edited", "wu")}
            res["seg_paired_change"] = {a: A10.paired_change(s, "raw", a, f"{a}_ok") for a in ("masked", "edited", "wu")}
            res["seg_gee"] = A10.gee(s.assign(TIT=s.raw), "TIT", ["x", "chf"])
            res["n_segments"] = {g: int((s.group == g).sum()) for g in GROUPS}
        y = SY[SY.cfg == rd]
        if len(y):
            res["synth"] = {}
            for c, g in y.groupby("cond"):
                res["synth"][c] = {a: AN.wil(int(g[a].sum()), len(g)) for a in ("raw",)}
                res["synth"][c]["edited"] = AN.wil(int(g[g.edited_ok].edited.sum()), int(g.edited_ok.sum()))
                res["synth"][c]["masked"] = AN.wil(int(g[g.masked_ok].masked.sum()), int(g.masked_ok.sum()))
            for e in ("none", "S2_5", "E3_10", "E1"):
                g = y[y.cond.str.startswith("vdp") & y.cond.str.endswith(f"_{e}")]
                if len(g):
                    res["synth"][f"vdp_pooled_{e}"] = {"raw": AN.wil(int(g.raw.sum()), len(g)),
                                                       "edited": AN.wil(int(g[g.edited_ok].edited.sum()),
                                                                        int(g.edited_ok.sum()))}
        p = PM[PM.cfg == rd]
        if len(p):
            res["pm97"] = {}
            for L in (500, 2000):
                d = p[p.L == L].rename(columns={a: f"{a}|x" for a in ("raw", "rawf", "masked", "edited", "wu")})
                d = d.rename(columns={f"{a}_ok": f"{a}|x|ok" for a in ("raw", "rawf", "masked", "edited", "wu")})
                for a in ("raw", "rawf", "masked", "edited", "wu"):
                    st = A12.pm97_stats(d, f"{a}|x", f"{a}|x|ok")
                    st.pop("per_subject", None)
                    res["pm97"][f"{L}|{a}"] = st
            d5 = p[p.L == 500]
            res["pm97"]["clogit_ectopy"] = A11.clogit(d5.assign(pos=d5.raw.astype(float),
                                                                x=np.log2(1 + d5.ectopy.astype(float))), "pos", ["x"])
            clean = d5[(d5.ectopy == 0) & (d5.n_out == 0)]
            res["pm97"]["clean_segments"] = A10.group_rates(clean, "raw")
        out["readings"][rd] = res
    return out


# ============================================================================================ stored other readings
def stored_readings():
    """Phase 12 stored results for the claims check (R0-BP96F, R0-S-FT, R0-S-AAFT)."""
    A = json.load(open(P12 / "analysis" / "analysis12.json"))
    out = {}
    for rd in ("R0-BP96F", "R0-S-FT", "R0-S-AAFT", "R0", "R0-BP96"):
        m = A.get("B4_matrix", {}).get("readings", {}).get(rd, {})
        tp = A.get("B4_transplant", {}).get("readings", {}).get(rd, {})
        sy = A.get("B4_synth", {})
        out[rd] = {"gee_OR": m.get("i_gee", {}).get("OR"), "gee_ci": m.get("i_gee", {}).get("ci95"),
                   "win_raw_diff": (m.get("ii_windows") or {}).get("diff"), "win_raw_ci": (m.get("ii_windows") or {}).get("ci95"),
                   "win_masked": (m.get("iii_windows", {}).get("masked") or {}),
                   "win_edited": (m.get("iii_windows", {}).get("edited") or {}),
                   "transplant_V1": tp.get("V1_minus_V0", {}).get("pooled"), "transplant_rates": tp.get("rates"),
                   "couplets5": sy.get("E3_couplets_5pct", {}).get(rd), "runs": sy.get("E4_runs", {}).get(rd),
                   "isolated5": sy.get("S2_ectopic_5pct", {}).get(rd), "N1": sy.get("N1_linear_rr", {}).get(rd),
                   "N5": sy.get("N5_linear_rr_warped", {}).get(rd), "raw_rate": m.get("rates", {}).get("raw")}
    return out


# ============================================================================================ M2
def m2_frame():
    R = AN.load("m2")
    rows, terms, curves = [], [], []
    for r in R:
        for rd in ("R0", "R0-ORIG"):
            o = r.get(rd)
            if not o or not o.get("analysable", True) or "nonlinear" not in o:
                continue
            row = {"reading": rd, "type": r["type"], "placement": r["placement"], "k": r["k"], "seed": r["seed"]}
            for k_ in ("nonlinear", "p", "F", "r_lin", "r_nl", "kappa_lin", "kappa_nl", "d_nl", "eps_lin", "eps_nl",
                       "relimp", "lag1", "lag2", "frac_aff6", "frac_affk", "share_aff6", "share_affk", "nrss_lin_aff",
                       "nrss_nl_aff", "nrss_lin_oth", "nrss_nl_oth", "max_abs_z_aff", "max_abs_z_oth",
                       "lev_nl_aff_median", "lev_nl_oth_median", "cook_nl_aff_max", "cook_nl_oth_max",
                       "top8_cook_aff_share", "n", "test"):
                row[k_] = o.get(k_)
            ra, rs = o["rm_aff"], o["rm_sham"]
            row["rm_aff_ok"], row["rm_aff_pos"] = ra["analysable"], ra["nonlinear"]
            row["rm_aff_n"] = ra["n_removed"]
            row["sham_ok"] = bool(rs.get("available") and rs.get("analysable", True))
            row["sham_pos"] = rs.get("nonlinear")
            for m in M.COOK_M:
                c = o.get("rm_cook", {}).get(str(m))
                row[f"cook{m}_pos"] = None if c is None else (c["nonlinear"] if c["analysable"] else None)
            mr = o.get("rm_masked_rule")
            row["masked_rule_ok"] = None if mr is None else mr["analysable"]
            row["masked_rule_pos"] = None if mr is None else mr["nonlinear"]
            rows.append(row)
            if o["nonlinear"]:
                for lb, conc in o.get("terms", []):
                    terms.append({"reading": rd, "type": r["type"], "placement": r["placement"], "k": r["k"],
                                  "seed": r["seed"], "term": lb, "conc": conc})
            if "curves" in o:
                curves.append({"reading": rd, "type": r["type"], "placement": r["placement"], "k": r["k"],
                               "seed": r["seed"], **o["curves"], "nonlinear": o["nonlinear"]})
        rows_bp = {"reading": "R0-BP96", "type": r["type"], "placement": r["placement"], "k": r["k"], "seed": r["seed"],
                   "nonlinear": r["R0-BP96"]["nonlinear"]}
        rows.append(rows_bp)
    df = pd.DataFrame(rows)
    # k = 0 windows (type 'none') serve every type
    z = df[df.type == "none"]
    df = pd.concat([df[df.type != "none"]] + [z.assign(type=t) for t in M.TYPES], ignore_index=True)
    tm = pd.DataFrame(terms)
    if len(tm):
        z = tm[tm.type == "none"]
        tm = pd.concat([tm[tm.type != "none"]] + [z.assign(type=t) for t in M.TYPES], ignore_index=True)
    return df, tm, curves


def mean_ci(v, seed=SEED):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if not len(v):
        return None, [None, None]
    rng = np.random.default_rng(seed)
    b = v[rng.integers(0, len(v), (2000, len(v)))].mean(1)
    return float(v.mean()), AN.pct(b)


def newcombe(k1, n1, k2, n2):
    l1, u1 = A10.wilson(k1, n1)
    l2, u2 = A10.wilson(k2, n2)
    p1, p2 = k1 / n1, k2 / n2
    d = p1 - p2
    return d, [d - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), d + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)]


def m2_tables(df, tm):
    inf, ssr, down = [], [], []
    for (rd, ty, pl, k), g in df.groupby(["reading", "type", "placement", "k"]):
        n = len(g)
        npos = int(g.nonlinear.sum())
        base = {"reading": rd, "type": ty, "placement": pl, "k": k, "n": n, "n_pos": npos,
                "pos_rate": npos / n, "pos_ci_lo": A10.wilson(npos, n)[0], "pos_ci_hi": A10.wilson(npos, n)[1]}
        if rd == "R0-BP96":
            inf.append(base)
            continue
        ok = g[g.rm_aff_ok == True]                                 # noqa: E712
        row = dict(base)
        row.update(rm_aff_n_ok=len(ok), rm_aff_rate=float(ok.rm_aff_pos.mean()) if len(ok) else None,
                   rm_aff_ci_lo=A10.wilson(int(ok.rm_aff_pos.sum()), len(ok))[0] if len(ok) else None,
                   rm_aff_ci_hi=A10.wilson(int(ok.rm_aff_pos.sum()), len(ok))[1] if len(ok) else None,
                   mean_rows_removed=float(ok.rm_aff_n.mean()) if len(ok) else None)
        sh = g[g.sham_ok == True]                                   # noqa: E712
        row.update(sham_n_ok=len(sh), sham_rate=float(sh.sham_pos.astype(bool).mean()) if len(sh) else None)
        P = g[g.nonlinear]
        Pa = P[P.rm_aff_ok == True]                                 # noqa: E712
        Ps = P[P.sham_ok == True]                                   # noqa: E712
        row.update(pos_flip_aff=float((~Pa.rm_aff_pos.astype(bool)).mean()) if len(Pa) else None, pos_flip_aff_n=len(Pa),
                   pos_flip_sham=float((~Ps.sham_pos.astype(bool)).mean()) if len(Ps) else None, pos_flip_sham_n=len(Ps))
        for m in M.COOK_M:
            c = P[f"cook{m}_pos"].dropna()
            row[f"pos_flip_cook{m}"] = float((~c.astype(bool)).mean()) if len(c) else None
        if k == 1:
            mr = P[P.masked_rule_ok == True]                        # noqa: E712
            row["pos_flip_masked_rule"] = float((~mr.masked_rule_pos.astype(bool)).mean()) if len(mr) else None
        row["top8_cook_aff_share"] = float(P.top8_cook_aff_share.mean()) if len(P) else None
        row["cook_nl_aff_max_median"] = float(P.cook_nl_aff_max.median()) if len(P) and k else None
        row["cook_nl_oth_max_median"] = float(P.cook_nl_oth_max.median()) if len(P) else None
        row["lev_aff_median"] = float(P.lev_nl_aff_median.median()) if len(P) and k else None
        row["lev_oth_median"] = float(P.lev_nl_oth_median.median()) if len(P) else None
        inf.append(row)
        if k > 0:
            for subset, G in (("positive", P), ("all", g)):
                sm, sci = mean_ci(G.share_aff6)
                fm = float(G.frac_aff6.mean()) if len(G) else None
                km, kci = mean_ci(G.share_affk)
                ssr.append({"reading": rd, "type": ty, "placement": pl, "k": k, "subset": subset, "n": len(G),
                            "share_aff_rows": sm, "share_ci_lo": sci[0], "share_ci_hi": sci[1],
                            "frac_aff_rows": fm, "enrichment": (sm / fm) if (sm is not None and fm) else None,
                            "share_aff_rows_kappa_nl": km, "share_kappa_nl_ci_lo": kci[0], "share_kappa_nl_ci_hi": kci[1]})
        drow = dict(base)
        for c in ("lag1", "lag2", "eps_lin", "eps_nl", "relimp", "nrss_lin_aff", "nrss_nl_aff", "nrss_lin_oth",
                  "nrss_nl_oth", "max_abs_z_aff", "max_abs_z_oth"):
            v = g[c].astype(float)
            if c in ("eps_lin", "eps_nl"):
                v = v ** 2
                c = c + "_sq"
            mval, ci = mean_ci(v)
            drow[c] = mval
            drow[f"{c}_ci_lo"], drow[f"{c}_ci_hi"] = ci
        for c in ("kappa_lin", "r_nl", "kappa_nl"):
            drow[f"{c}_median"] = float(g[c].median())
            drow[f"{c}_q25"], drow[f"{c}_q75"] = float(g[c].quantile(0.25)), float(g[c].quantile(0.75))
        drow["d_nl_mode"] = int(g.d_nl.mode().iloc[0])
        drow["F_median"] = float(g.F.median())
        down.append(drow)
    # term selection among positive windows
    trows = []
    if len(tm):
        npos = df[df.nonlinear & (df.reading != "R0-BP96")].groupby(["reading", "type", "placement", "k"]).size()
        for (rd, ty, pl, k, term), g in tm.groupby(["reading", "type", "placement", "k", "term"]):
            N = int(npos.get((rd, ty, pl, k), 0))
            trows.append({"reading": rd, "type": ty, "placement": pl, "k": k, "term": term, "kind": M.term_kind(term),
                          "lags": term[1:].replace("*y", ","), "n_pos_windows": N, "n_selected": len(g),
                          "freq": len(g) / N if N else None, "mean_ectopic_concentration": float(g.conc.mean()),
                          "flag_power_or_adjacent": M.term_kind(term) in ("power", "adjacent_product")})
    return pd.DataFrame(inf), pd.DataFrame(ssr), pd.DataFrame(down), pd.DataFrame(trows)


def paired_diff(df, rd, col, k1, k2, ty="ventricular", pl="random", seed=SEED + 21):
    a = df[(df.reading == rd) & (df.type == ty) & (df.placement == pl) & (df.k == k1)].set_index("seed")[col]
    b = df[(df.reading == rd) & (df.type == ty) & (df.placement == pl) & (df.k == k2)].set_index("seed")[col]
    j = pd.concat([a, b], axis=1, keys=["a", "b"]).dropna()
    d = (j.b - j.a).astype(float).values
    rng = np.random.default_rng(seed)
    bs = d[rng.integers(0, len(d), (5000, len(d)))].mean(1)
    return {"k1": k1, "k2": k2, "mean_k1": float(j.a.mean()), "mean_k2": float(j.b.mean()), "diff": float(d.mean()),
            "ci": AN.pct(bs), "n": len(d)}


def m2_rules(df, inf):
    out = {}
    for rd in ("R0", "R0-ORIG"):
        res = {}
        for ty in M.TYPES:
            g = df[(df.reading == rd) & (df.type == ty) & (df.placement == "random") & df.k.isin([1, 2, 4, 8]) & df.nonlinear]
            ga = g[g.rm_aff_ok == True]                             # noqa: E712
            gs = g[g.sham_ok == True]                               # noqa: E712
            fa = float((~ga.rm_aff_pos.astype(bool)).mean()) if len(ga) else None
            fs = float((~gs.sham_pos.astype(bool)).mean()) if len(gs) else None
            verdict = None
            if fa is not None and fs is not None:
                verdict = ("driven by a few high-leverage equations" if fa >= 0.8 and fs <= 0.2 else
                           ("largely driven by a few equations" if fa >= 0.5 else "not driven by a few equations"))
            cook = {}
            for k in (1, 2, 4, 8, 16):
                r = inf[(inf.reading == rd) & (inf.type == ty) & (inf.placement == "random") & (inf.k == k)]
                if len(r):
                    r = r.iloc[0]
                    mm = [m for m in M.COOK_M if r.get(f"pos_flip_cook{m}") is not None and r[f"pos_flip_cook{m}"] >= 0.8]
                    cook[str(k)] = min(mm) if mm else None
            res[ty] = {"n_pos_analysable": len(ga), "flip_after_removal": fa, "n_pos_sham": len(gs),
                       "flip_after_sham": fs, "verdict": verdict, "smallest_cook_m_80pct": cook}
        out[rd] = {"few_equations": res}
        # H2b downslope (ventricular, random)
        d = df[(df.reading == rd)]
        lag = {k: float(d[(d.type == "ventricular") & (d.placement == "random") & (d.k == k)].lag1.mean())
               for k in (8, 16, 26, 51)}
        a = paired_diff(df, rd, "lag1", 8, 51)
        b = paired_diff(df, rd, "nrss_lin_aff", 8, 51)
        c = paired_diff(df, rd, "relimp", 8, 51)
        e = paired_diff(df, rd, "nrss_nl_aff", 8, 51)
        f = paired_diff(df, rd, "max_abs_z_aff", 8, 51)
        g = paired_diff(df, rd, "r_nl", 8, 51)
        mono = all(lag[x] > lag[y] for x, y in ((8, 16), (16, 26), (26, 51)))
        ha = bool(a["ci"][1] < 0 and mono)
        hb = bool(b["ci"][1] < 0)
        hc = bool(c["ci"][1] < 0)
        per = {}
        for k in (16, 26, 51):
            gr = df[(df.reading == rd) & (df.type == "ventricular") & (df.placement == "random") & (df.k == k)]
            gp = df[(df.reading == rd) & (df.type == "ventricular") & (df.placement == "periodic") & (df.k == k)]
            dd, ci = newcombe(int(gp.nonlinear.sum()), len(gp), int(gr.nonlinear.sum()), len(gr))
            per[str(k)] = {"periodic": float(gp.nonlinear.mean()), "random": float(gr.nonlinear.mean()),
                           "diff_periodic_minus_random": dd, "ci": [float(ci[0]), float(ci[1])]}
        out[rd]["downslope"] = {"lag1_by_k": lag, "a_lag1": a, "a_monotone": mono, "b_nrss_lin_aff": b,
                                "c_relimp": c, "alt_nrss_nl_aff": e, "alt_max_abs_z_aff": f, "alt_r_nl": g,
                                "H2b_supported": bool(ha and hb and hc), "parts": {"a": ha, "b": hb, "c": hc},
                                "periodic_vs_random": per,
                                "periodic_lower_at_16": bool(per["16"]["ci"][1] < 0)}
    return out


def m2_section():
    df, tm, curves = m2_frame()
    if df.empty:
        return {}
    inf, ssr, down, trows = m2_tables(df, tm)
    R = AN.RES
    inf.to_csv(R / "m2_influence.csv", index=False, float_format="%.4g")
    ssr.to_csv(R / "m2_ssr_share.csv", index=False, float_format="%.4g")
    down.to_csv(R / "m2_downslope.csv", index=False, float_format="%.4g")
    trows.to_csv(R / "m2_term_selection.csv", index=False, float_format="%.4g")
    (R / "m2_curves.json").write_text(json.dumps(curves))
    return {"n_errors": AN.n_errors("m2"), "n_windows": int(len(df[df.reading == "R0"])), "rules": m2_rules(df, inf),
            "explore11_agreement": explore_agreement(df)}


def explore_agreement(df):
    st = {r["id"]: r["positive"] for r in jl(P11 / "explore_dose.jsonl") if not r.get("error")}
    d = df[(df.reading == "R0") & (df.type == "ventricular") & (df.placement == "random") & (df.seed < 11100)]
    ag = sum(bool(st[f"{k}|{s}"]) == bool(v) for k, s, v in zip(d.k, d.seed, d.nonlinear) if f"{k}|{s}" in st)
    return {"agree": int(ag), "n": int(sum(f"{k}|{s}" in st for k, s in zip(d.k, d.seed)))}


# ============================================================================================ entry
def sections(want):
    out = {}
    if (want is None or "head" in want) and (AN.RAW / "head.jsonl").exists():
        out["headline"] = head_section()
        out["stored_readings"] = stored_readings()
    if (want is None or "m2" in want) and (AN.RAW / "m2.jsonl").exists():
        out["m2"] = m2_section()
    return out
