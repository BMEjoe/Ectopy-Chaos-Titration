"""
Phase 11 analyses (PREREGISTRATION.md Sections 2-7) on run11's JSONL output and the Phase 10 results.

    python -m experiments.phase11_robustness.analysis11 --phase conf [--only A1,A2,...]
    python -m experiments.phase11_robustness.analysis11 --phase pilot      # development pilot (debugging)

Writes results/<phase>/analysis/analysis11.json (merged: sections already present are kept unless recomputed).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from experiments.phase10_final import analysis10 as A10

HERE = pathlib.Path(__file__).resolve().parent
P10 = HERE.parent / "phase10_final" / "results"
B = 10000
SEED = 20261104
GROUPS = ("NSR", "CHF")


# ----------------------------------------------------------------------------------------------- io
def p10_phase(phase):
    return "dev" if phase == "pilot" else "conf"


def load(phase, part):
    p = HERE / "results" / phase / f"{part}.jsonl"
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


def n_errors(phase, part):
    p = HERE / "results" / phase / f"{part}.jsonl"
    return sum(1 for line in open(p) if json.loads(line).get("error")) if p.exists() else None


def p10_real(phase):
    df = A10.real_frame(A10.load(p10_phase(phase), "real"))
    df = df[df.group.isin(GROUPS)].copy()
    df["wid"] = df["id"]
    return df


def p10_seg(phase):
    sg = A10.seg_frame(A10.load(p10_phase(phase), "seg"))
    return sg[~sg.gap].copy()


# ----------------------------------------------------------------------------------------------- statistics
def cluster_sums(df, cols, subj="subject"):
    g = df.groupby(subj)[list(cols)].sum()
    return g.values.astype(float), list(g.index)


def boot_idx(n, B_=B, seed=SEED):
    return np.random.default_rng(seed).integers(0, n, size=(B_, n))


def ratio_boot(df, cols, fn, by_group=False, seed=SEED, B_=B):
    """Subject-cluster bootstrap of a statistic fn(sums) of column sums (vectorised: sums has shape (B, k)).
    by_group: resample subjects within NSR and CHF separately (stratified)."""
    if df.empty:
        return float("nan"), [float("nan"), float("nan")]
    est = float(fn(df[list(cols)].sum().values[None, :].astype(float))[0])
    if by_group:
        tot = 0
        for i, g in enumerate(GROUPS):
            d = df[df.group == g]
            if d.empty:
                continue
            S, _ = cluster_sums(d, cols)
            tot = tot + S[boot_idx(len(S), B_, seed + 7 * i)].sum(1)
        bs = fn(tot)
    else:
        S, _ = cluster_sums(df, cols)
        bs = fn(S[boot_idx(len(S), B_, seed)].sum(1))
    bs = np.asarray(bs, float)
    bs = bs[np.isfinite(bs)]
    return est, ([float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if len(bs) else [float("nan")] * 2)


def mcnemar(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    n10, n01 = int((a & ~b).sum()), int((~a & b).sum())
    p = float(stats.binomtest(min(n10, n01), n10 + n01, 0.5).pvalue) if n10 + n01 else 1.0
    return {"a_only": n10, "b_only": n01, "p_exact": p}


def wilson(k, n):
    return A10.wilson(int(k), int(n))


def paired_delta(df, a, b, by_group_list=GROUPS, seed=SEED):
    """rate(b) - rate(a) over units (window level, paired), cluster bootstrap, exact McNemar; by group and pooled."""
    out = {}
    for g in tuple(by_group_list) + ("pooled",):
        d = df if g == "pooled" else df[df.group == g]
        if d.empty:
            continue
        dd = d.assign(_a=d[a].astype(float), _b=d[b].astype(float), _n=1.0)
        est, ci = ratio_boot(dd, ["_a", "_b", "_n"], lambda s: (s[:, 1] - s[:, 0]) / s[:, 2], seed=seed)
        out[g] = {"n_units": len(d), "n_subjects": int(d.subject.nunique()), "rate_a": float(d[a].mean()),
                  "rate_b": float(d[b].mean()), "delta": est, "ci95": ci, "mcnemar": mcnemar(d[a], d[b])}
    return out


# ----------------------------------------------------------------------------------------------- A1
def _a1_units(units, arm):
    """units: DataFrame rows with true_pos (bool), true_ok, raw, kind, shams (list of dicts)."""
    rows = []
    for u in units:
        sh = [s for s in u["shams"] if s["ok"]]
        rows.append({**{k: u[k] for k in ("subject", "group", "uid", "kind", "raw", "true_ok", "true_pos")},
                     "n_sham_ok": len(sh), "n_sham": u["n_sham"],
                     "S": float(np.mean([not s["pos"] for s in sh])) if sh else np.nan,
                     "S_created": float(np.mean([s["pos"] for s in sh])) if sh else np.nan})
    return pd.DataFrame(rows)


def a1_cell(df):
    """Primary A1 statistic on ectopic raw-positive units usable in both arms."""
    d = df[(df.kind == "own") & df.raw & df.true_ok & (df.n_sham_ok >= 1)].copy()
    out = {"n_units": len(d), "n_subjects": int(d.subject.nunique()) if len(d) else 0}
    if d.empty:
        return out
    d["R"] = (~d.true_pos.astype(bool)).astype(float)
    d["_n"] = 1.0
    res = {}
    for g in GROUPS + ("pooled",):
        dd = d if g == "pooled" else d[d.group == g]
        if dd.empty:
            continue
        D, Dci = ratio_boot(dd, ["R", "S", "_n"], lambda s: (s[:, 0] - s[:, 1]) / s[:, 2])
        Q, Qci = ratio_boot(dd, ["R", "S", "_n"], lambda s: s[:, 1] / np.where(s[:, 0] > 0, s[:, 0], np.nan))
        ft, fs = float(dd.R.mean()), float(dd.S.mean())
        verdict = "supports" if (Dci[0] > 0 and Q < 0.5) else "against"
        res[g] = {"n_units": len(dd), "n_subjects": int(dd.subject.nunique()), "F_true": ft, "F_sham": fs,
                  "D": D, "D_ci95": Dci, "Q": Q, "Q_ci95": Qci, "verdict": verdict}
    out.update(res)
    # secondary: clean units (donor masks) and created positives
    c = df[(df.kind == "donor") & df.raw & (df.n_sham_ok >= 1)]
    out["clean_raw_positive_units"] = len(c)
    if len(c):
        cc = c.assign(_n=1.0)
        est, ci = ratio_boot(cc, ["S", "_n"], lambda s: s[:, 0] / s[:, 1])
        out["clean_removed_by_donor_mask"] = {"fraction": est, "ci95": ci, "n_subjects": int(c.subject.nunique())}
    neg = df[(df.kind == "own") & ~df.raw & df.true_ok & (df.n_sham_ok >= 1)]
    if len(neg):
        out["created_among_raw_negative"] = {"n": len(neg), "true": float(neg.true_pos.mean()),
                                             "sham": float(neg.S_created.mean())}
    out["sham_unavailable_units"] = int((df.kind == "own").sum() - ((df.kind == "own") & (df.n_sham >= 1)).sum())
    out["units_with_fewer_shams"] = int(((df.kind == "own") & (df.n_sham < df.n_sham.max())).sum())
    return out


def a1(phase):
    res = {}
    df10 = p10_real(phase).set_index("wid")
    W = load(phase, "sham_win")
    for arm, key, okf in (("masked", "TITm", "TITm_analysable"), ("edited", "TITe", "TITe_eligible")):
        units = []
        for r in W:
            wid = r["id"].split("|", 1)[1]
            if wid not in df10.index:
                continue
            w = df10.loc[wid]
            sh = []
            for s in r["shams"]:
                o = s[arm]
                ok = bool(o.get("analysable")) if arm == "masked" else bool(o.get("eligible") and o.get("analysable", True))
                sh.append({"ok": ok, "pos": bool(o.get("positive"))})
            units.append({"subject": r["subject"], "group": r["group"], "uid": wid, "kind": r["kind"],
                          "raw": bool(w.TIT), "true_ok": bool(w[okf]), "true_pos": bool(w[key]), "shams": sh,
                          "n_sham": r["n_sham"]})
        res[f"windows_{arm}"] = a1_cell(_a1_units(units, arm))
    S = load(phase, "sham_seg")
    if S:
        sg = p10_seg(phase).set_index(["subject", "seg"])
        for arm, key, okf in (("masked", "masked", "masked_ok"), ("edited", "edited", "edited_ok"), ("wu", "wu", "wu_ok")):
            units = []
            for r in S:
                for s in r["segments"]:
                    k = (r["subject"], s["seg"])
                    if k not in sg.index:
                        continue
                    w = sg.loc[k]
                    sh = []
                    for x in s["shams"]:
                        o = x[arm]
                        if arm == "masked":
                            ok = bool(o.get("analysable"))
                        elif arm == "edited":
                            ok = bool(o.get("eligible") and o.get("analysable", True))
                        else:
                            ok = bool(o.get("eligible", True) and o.get("analysable", True))
                        sh.append({"ok": ok, "pos": bool(o.get("positive"))})
                    units.append({"subject": r["subject"], "group": r["group"], "uid": f"{k[0]}|{k[1]}",
                                  "kind": s["kind"], "raw": bool(w.raw), "true_ok": bool(w[okf]),
                                  "true_pos": bool(w[key]), "shams": sh, "n_sham": s["n_sham"]})
            res[f"segments_{arm}"] = a1_cell(_a1_units(units, arm))
    cells = [k for k in ("windows_masked", "windows_edited", "segments_masked", "segments_edited", "segments_wu")
             if k in res and "pooled" in res[k]]
    res["summary_verdict"] = ("supports" if cells and all(res[k]["pooled"]["verdict"] == "supports" for k in cells)
                              else "weakens") if len(cells) == 5 else f"incomplete ({len(cells)}/5 cells)"
    return res


# ----------------------------------------------------------------------------------------------- A2
def a2_hrt(phase):
    R = load(phase, "hrt")
    out = {}
    for g in GROUPS:
        rr = [r for r in R if r["group"] == g]
        s1 = [r["summary"] for r in rr if r["summary"]["n_vpb"] >= 1]
        s5 = [s for s in s1 if s["n_vpb"] >= 5]

        def desc(L, k):
            v = np.array([s[k] for s in L], float)
            return {"median": float(np.median(v)), "q25": float(np.percentile(v, 25)), "q75": float(np.percentile(v, 75)),
                    "n": len(v)} if len(v) else None
        out[g] = {"n_subjects": len(rr), "n_subjects_with_V": int(sum(r["n_V"] > 0 for r in rr)),
                  "n_subjects_ge1_usable": len(s1), "n_subjects_ge5_usable": len(s5),
                  "n_usable_vpb": int(sum(s["n_vpb"] for s in s1)), "n_V_total": int(sum(r["n_V"] for r in rr)),
                  "TO_ge5": desc(s5, "TO"), "TS_ge5": desc(s5, "TS"), "TS20_ge5": desc(s5, "TS20"),
                  "TO_ge1": desc(s1, "TO"), "TS_ge1": desc(s1, "TS"),
                  "frac_TO_abnormal_ge5": float(np.mean([s["TO"] >= 0 for s in s5])) if s5 else None,
                  "frac_TS_abnormal_ge5": float(np.mean([s["TS"] <= 2.5 for s in s5])) if s5 else None,
                  "frac_both_abnormal_ge5": float(np.mean([(s["TO"] >= 0) and (s["TS"] <= 2.5) for s in s5])) if s5 else None,
                  "profile_h": np.mean([s["h"] for s in s1], axis=0).tolist() if s1 else None,
                  "ci_ratio_mean": float(np.mean([s["ci_ratio"] for s in s1])) if s1 else None,
                  "cp_ratio_mean": float(np.mean([s["cp_ratio"] for s in s1])) if s1 else None,
                  "per_subject": [{"subject": r["subject"], "n_V": r["n_V"], "n_vpb": r["summary"]["n_vpb"],
                                   "TO": r["summary"].get("TO"), "TS": r["summary"].get("TS")} for r in rr]}
    from experiments.phase11_robustness import hrt11 as HRT
    for g in GROUPS:
        h = out[g]["profile_h"]
        if h:
            out[g]["profile_TO_TS_at_0.8s"] = HRT.profile_metrics(h, 0.8)
    out["reference_values"] = {"healthy_median_TS_ms_per_RR": 13.0, "healthy_median_TO_pct": -2.5,
                               "abnormal_TS_le": 2.5, "abnormal_TO_ge": 0.0,
                               "source": "Grimm et al. 2003 (values as given in the task; not verified here)"}
    return out


def _variant_frame(R, keys=("subject", "group")):
    rows = []
    for r in R:
        base = {k: r.get(k) for k in keys}
        for v in ("V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7"):
            if "positive" in r.get(v, {}):
                base[v] = bool(r[v]["positive"])
        rows.append(base)
    return pd.DataFrame(rows)


def a2_transplant(phase):
    R = load(phase, "transplant")
    df = _variant_frame(R)
    out = {"n_windows": len(df), "n_subjects": int(df.subject.nunique()) if len(df) else 0}
    if df.empty:
        return out
    out["rates"] = {v: {g: {"k": int(df[df.group == g][v].sum()), "n": int((df.group == g).sum()),
                            "wilson": wilson(df[df.group == g][v].sum(), (df.group == g).sum())}
                        for g in GROUPS if (df.group == g).any()} for v in df.columns if v.startswith("V")}
    out["delta_vs_V0"] = {v: paired_delta(df.dropna(subset=[v]), "V0", v) for v in df.columns
                          if v.startswith("V") and v != "V0"}
    return out


def a2_synth(phase):
    R = load(phase, "synth1")
    out = {}
    rows = []
    for r in R:
        c = r["cond"] if r["kind"] == "null" else f"vdp{r['cond']}"
        rows.append({"cond": c, **{v: bool(r[v]["positive"]) for v in ("V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7")
                                   if "positive" in r.get(v, {})}})
    df = pd.DataFrame(rows)
    rng = np.random.default_rng(SEED + 5)
    for c, g in df.groupby("cond"):
        n = len(g)
        idx = rng.integers(0, n, size=(B, n))
        res = {"n": n}
        for v in [x for x in g.columns if x.startswith("V")]:
            k = int(g[v].sum())
            res[v] = {"k": k, "wilson": wilson(k, n)}
            if v != "V0":
                a, b = g["V0"].values.astype(float), g[v].values.astype(float)
                bs = (b[idx] - a[idx]).mean(1)
                res[v]["delta"] = float(b.mean() - a.mean())
                res[v]["ci95"] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
                res[v]["mcnemar"] = mcnemar(g["V0"], g[v])
        out[c] = res
    return out


def a2_extmask(phase):
    R = load(phase, "extmask")
    df10 = p10_real(phase).set_index("wid")
    rows = []
    for r in R:
        wid = r["id"].split("|", 1)[1]
        rows.append({"subject": r["subject"], "group": r["group"], "raw": bool(df10.loc[wid].TIT) if wid in df10.index else None,
                     "std_ok": bool(r["masked_std"].get("analysable")), "ext_ok": bool(r["masked_ext"].get("analysable")),
                     "std": bool(r["masked_std"].get("positive")), "ext": bool(r["masked_ext"].get("positive")),
                     "n_ext_masked": r["n_ext_masked"], "n_masked": r["n_masked"]})
    df = pd.DataFrame(rows)
    if df.empty:
        return {}
    both = df[df.std_ok & df.ext_ok]
    return {"n_windows": len(df), "n_analysable_both": len(both), "n_std_only": int((df.std_ok & ~df.ext_ok).sum()),
            "all_units": paired_delta(both, "std", "ext"),
            "raw_positive_units": paired_delta(both[both.raw == True], "std", "ext"),   # noqa: E712
            "mean_masked_std": float(df.n_masked.mean()), "mean_masked_ext": float(df.n_ext_masked.mean())}


SV = set("AaJS")


def covar_frames(phase):
    R = load(phase, "covar")
    W, S = [], []
    for r in R:
        for w in r["windows"]:
            W.append({"subject": r["subject"], "group": r["group"], **w})
        for s in r["segments"]:
            S.append({"subject": r["subject"], "group": r["group"], **s})
    return pd.DataFrame(W), pd.DataFrame(S)


def a2_beat_type(phase):
    cw, _ = covar_frames(phase)
    df = p10_real(phase).merge(cw[["wid", "abnormal_symbols", "n_rr_flag_only"]], on="wid", how="left")

    def cat(r):
        if r.burden == 1:
            s = r.abnormal_symbols[0]
            return "V" if s == "V" else ("supraventricular" if s in SV else "other_unknown")
        if r.burden == 0 and r.n_masked > 0:
            return "rr_rule_only"
        if r.burden == 0:
            return "clean"
        return None
    df["cat"] = [cat(r) for r in df.itertuples()]
    out = {}
    for c in ("clean", "rr_rule_only", "V", "supraventricular", "other_unknown"):
        d = df[df.cat == c]
        k, n = int(d.TIT.sum()), len(d)
        out[c] = {"n": n, "n_subjects": int(d.subject.nunique()), "k": k, "rate": k / n if n else None,
                  "wilson": wilson(k, n),
                  "symbols": d.abnormal_symbols.apply(lambda s: s[0] if len(s) else "").value_counts().to_dict()
                  if c in ("other_unknown", "supraventricular") else None}
    return out


# ----------------------------------------------------------------------------------------------- A3
CONT = ("mean_nn", "log_sdnn", "log_rmssd", "half_shift", "log1p_trend")


def _prep_cov(d):
    d = d.copy()
    d["log_sdnn"] = np.log(d["sdnn_nn"])
    d["log_rmssd"] = np.log(d["rmssd_nn"])
    d["log1p_trend"] = np.log1p(d["trend_ratio"])
    d["eh_bin"] = np.clip(np.floor(d["elapsed_h"] / 6.0), 0, 3).astype(int)
    for b in (1, 2, 3):
        d[f"eh{b}"] = (d["eh_bin"] == b).astype(float)
    ok = np.all(np.isfinite(d[list(CONT)].values), axis=1)
    n_drop = int((~ok).sum())
    d = d[ok].copy()
    for c in CONT:
        sd = d[c].std()
        d[c] = (d[c] - d[c].mean()) / (sd if sd > 0 else 1.0)
    if "x" not in d:
        d["x"] = np.log2(1 + d["burden"])
    d["chf"] = (d["group"] == "CHF").astype(int)
    return d, n_drop


COVS = list(CONT) + ["eh1", "eh2", "eh3"]


def clogit(d, outcome, covars, boot=0, seed=SEED + 31):
    from statsmodels.discrete.conditional_models import ConditionalLogit
    y = d[outcome].astype(float).values
    out = {"outcome": outcome, "covariates": covars, "n": len(d), "n_pos": int(y.sum()),
           "n_subjects": int(d.subject.nunique())}
    g = d.groupby("subject")[outcome].agg(["min", "max"])
    inf = g[g["min"] != g["max"]].index
    out["n_informative_subjects"] = int(len(inf))
    di = d[d.subject.isin(inf)]
    out["n_informative_units"] = len(di)
    if y.sum() < A10.MIN_POS_GEE or len(d) - y.sum() < A10.MIN_POS_GEE or len(inf) < 2:
        out["estimated"] = False
        return out
    cv = [c for c in covars if di[c].std() > 0]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = ConditionalLogit(di[outcome].astype(float).values, di[cv].astype(float).values,
                                 groups=di.subject.values).fit(disp=0)
        b, se = np.asarray(m.params), np.asarray(m.bse)
        out["estimated"] = bool(np.all(np.isfinite(b)) and np.all(np.isfinite(se)))
        out["coef"] = {c: {"beta": float(bb), "se": float(s), "OR": float(np.exp(bb)),
                           "OR_ci95": [float(np.exp(bb - 1.96 * s)), float(np.exp(bb + 1.96 * s))],
                           "p": float(2 * stats.norm.sf(abs(bb / s)))} for c, bb, s in zip(cv, b, se)}
        if boot and out["estimated"]:
            rng = np.random.default_rng(seed)
            subs = np.array(sorted(inf))
            by = {s: di[di.subject == s] for s in subs}
            bx = []
            for _ in range(boot):
                pick = rng.choice(subs, len(subs), replace=True)
                dd = pd.concat([by[s].assign(subject=f"{s}#{i}") for i, s in enumerate(pick)])
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        mm = ConditionalLogit(dd[outcome].astype(float).values, dd[cv].astype(float).values,
                                              groups=dd.subject.values).fit(disp=0)
                    bx.append(float(np.asarray(mm.params)[cv.index("x")]))
                except Exception:                               # noqa: BLE001
                    pass
            bx = np.asarray(bx)
            bx = bx[np.isfinite(bx)]
            out["x_boot"] = {"n_ok": len(bx), "OR_ci95": [float(np.exp(np.percentile(bx, 2.5))),
                                                          float(np.exp(np.percentile(bx, 97.5)))] if len(bx) else None}
    except Exception as exc:                                    # noqa: BLE001
        out["estimated"] = False
        out["reason"] = f"{type(exc).__name__}: {exc}"
    return out


def a3_models(d, outcome, boot):
    res = {"gee_cov": A10.gee(d, outcome, ["x", "chf"] + COVS),
           "gee_prereg": A10.gee(d, outcome, ["x", "chf"]),
           "clogit_cov": clogit(d, outcome, ["x"] + COVS, boot=boot),
           "clogit_x": clogit(d, outcome, ["x"])}
    for g in GROUPS:
        dg = d[d.group == g]
        res[f"{g}_gee_cov"] = A10.gee(dg, outcome, ["x"] + COVS)
        res[f"{g}_clogit_cov"] = clogit(dg, outcome, ["x"] + COVS)
        res[f"{g}_clogit_x"] = clogit(dg, outcome, ["x"])
    return res


def a3(phase, boot=2000):
    cw, cs = covar_frames(phase)
    out = {"clock_time": "unavailable (no header start time); elapsed hours since the first beat used"}
    df = p10_real(phase).merge(cw.drop(columns=["subject", "group", "j"]), on="wid", how="inner")
    d, nd = _prep_cov(df)
    out["windows"] = {"n": len(d), "n_dropped_undefined_covariate": nd,
                      "TIT": a3_models(d, "TIT", boot), "LLE": a3_models(d, "LLE", 0), "K1": a3_models(d, "K1", 0)}
    sg = p10_seg(phase)
    if len(cs) and len(sg):
        ds = sg.merge(cs.drop(columns=["group"]), on=["subject", "seg"], how="inner")
        ds, nds = _prep_cov(ds.assign(TIT=ds.raw.astype(bool)))
        out["segments"] = {"n": len(ds), "n_dropped_undefined_covariate": nds, "TIT": a3_models(ds, "TIT", 0)}
    c = out["windows"]["TIT"]["clogit_cov"]
    if c.get("estimated"):
        lo, hi = c["coef"]["x"]["OR_ci95"]
        out["verdict"] = "supports" if lo > 1 else ("against" if lo <= 1 <= hi else "against (OR < 1)")
    return out


# ----------------------------------------------------------------------------------------------- A4
def a4(phase):
    sg = p10_seg(phase)
    out = {"pm97_pdf": "not available; Poon & Merrill 1997 metric not extracted"}
    arms = {"raw": (None, "raw", "raw_NL"), "masked": ("masked_ok", "masked", None), "edited": ("edited_ok", "edited", None),
            "wu": ("wu_ok", "wu", "wu_NL")}
    sg["block3h"] = np.clip(np.floor(sg["t_start"].astype(float) / 10800.0), 0, 7).fillna(-1).astype(int)
    for arm, (ok, col, nl) in arms.items():
        d = sg if ok is None else sg[sg[ok]]
        res = {}
        for g in GROUPS:
            dg = d[d.group == g]
            fr = dg.groupby("subject")[col].mean()
            n = len(fr)
            res[g] = {"n_subjects": n}
            for th, name in ((0, "ge1_positive"), (0.10, "ge10pct"), (0.25, "ge25pct"), (0.50, "ge50pct")):
                k = int((fr > 0).sum()) if th == 0 else int((fr >= th).sum())
                res[g][name] = {"k": k, "frac": k / n if n else None, "wilson": wilson(k, n)}
            blocks = {}
            for b_ in range(8):
                db = dg[dg.block3h == b_]
                if db.empty:
                    continue
                v, _ = A10.subject_rates(db, col)
                blocks[f"{3 * b_}-{3 * b_ + 3}h"] = {"n_subjects": len(v), "DR": float(np.mean(v)),
                                                      "ci95": A10.ci(A10.boot_mean(v))}
            res[g]["DR_per_3h_block"] = blocks
            if nl is not None:
                pos = dg[dg[col]]
                res[g]["mean_NL_positive"] = float(pos[nl].mean()) if len(pos) else None
        out[arm] = res
    return out


# ----------------------------------------------------------------------------------------------- A5
def verify_reading(R):
    from experiments.phase10_final import verify_titration as VT
    return VT.evaluate(R)


def a5(phase):
    out = {"A5_1": "Barahona & Poon 1996 PDF not available: definitions not extracted; R8 not run"}
    V = load(phase, "rdg_verify")
    ver = {}
    for rd in sorted({r["reading"] for r in V}):
        rr = [r for r in V if r["reading"] == rd]
        try:
            v = verify_reading(rr)
            ver[rd] = {k: (v[k]["pass"] if isinstance(v[k], dict) else v[k]) for k in ("V1", "V2", "V3", "V4", "V5", "V6", "VERIFIED")}
            ver[rd]["detail"] = {k: {kk: vv for kk, vv in v[k].items() if kk != "scan"} for k in ("V1", "V2", "V3", "V4", "V5", "V6")}
        except Exception as exc:                                # noqa: BLE001
            ver[rd] = {"error": f"{type(exc).__name__}: {exc} (incomplete: {len(rr)} records)"}
    out["verification"] = ver
    df10 = p10_real(phase)
    frames = {"R0": df10[["wid", "subject", "group", "burden", "x", "chf", "TIT", "TITm", "TITm_analysable", "TITe",
                          "TITe_eligible"]].rename(columns={"TIT": "raw", "TITm": "masked", "TITm_analysable": "masked_ok",
                                                            "TITe": "edited", "TITe_eligible": "edited_ok"})}
    base = df10.set_index("wid")[["subject", "group", "burden", "x", "chf"]]
    W = load(phase, "rdg_win")
    for rd in sorted({r["reading"] for r in W}):
        rows = []
        for r in W:
            if r["reading"] != rd:
                continue
            wid = r["id"].split("|", 1)[1].rsplit("|", 1)[0]
            rows.append({"wid": wid, "raw": bool(r["raw"]["positive"]), "masked": bool(r["masked"]["positive"]),
                         "masked_ok": bool(r["masked"].get("analysable")), "edited": bool(r["edited"]["positive"]),
                         "edited_ok": bool(r["edited"].get("eligible") and r["edited"].get("analysable", True))})
        f = pd.DataFrame(rows).set_index("wid").join(base, how="inner").reset_index()
        frames[rd] = f
    Sy = load(phase, "rdg_synth")
    syn = {}
    for r in Sy:
        c = r["cond"] if r["kind"] == "null" else f"vdp{r['cond']}"
        syn.setdefault(r["reading"], {}).setdefault(c, []).append(bool(r["raw"]["positive"]))
    segs = {"R0": p10_seg(phase).rename(columns={})}
    Sg = load(phase, "rdg_seg")
    for rd in sorted({r["reading"] for r in Sg}):
        rows = []
        for r in Sg:
            if r["reading"] != rd:
                continue
            for s in r["segments"]:
                rows.append({"subject": r["subject"], "group": r["group"], "seg": s["seg"],
                             "raw": bool(s["raw"]["positive"]), "masked": bool(s["masked"]["positive"]),
                             "masked_ok": bool(s["masked"].get("analysable")), "wu": bool(s["wu"]["positive"]),
                             "wu_ok": bool(s["wu"].get("eligible", True) and s["wu"].get("analysable", True))})
        segs[rd] = pd.DataFrame(rows)
    matrix = {}
    for rd, f in frames.items():
        m = {"n_windows": len(f)}
        g = A10.gee(f.assign(TITr=f.raw.astype(bool)), "TITr", ["x", "chf"])
        m["i_gee"] = {k: g.get(k) for k in ("estimated", "working_correlation", "n_pos")}
        if g.get("estimated"):
            c = g["coef"]["x"]
            m["i_gee"].update(OR=c["OR"], ci95=c["OR_ci95"])
            m["i"] = bool(c["OR_ci95"][0] > 1)
            m["i_reversed"] = bool(c["OR_ci95"][1] < 1)
        else:
            m["i"] = None
        gr = A10.group_rates(f, "raw")
        m["ii_windows"] = gr.get("CHF_minus_NSR")
        m["ii"] = bool(m["ii_windows"]["ci95"][0] > 0) if m["ii_windows"] else None
        post = {}
        for arm in ("masked", "edited"):
            ga = A10.group_rates(f, arm, f"{arm}_ok").get("CHF_minus_NSR")
            post[arm] = ga
        m["iii_windows"] = post
        m["iii"] = bool(all(post[a] and post[a]["ci95"][0] <= 0 <= post[a]["ci95"][1] for a in post))
        m["iii_reversed"] = bool(any(post[a] and not (post[a]["ci95"][0] <= 0 <= post[a]["ci95"][1]) for a in post))
        m["rates"] = {"raw": float(f.raw.mean()), "masked_analysable": float(f[f.masked_ok].masked.mean()),
                      "edited_eligible": float(f[f.edited_ok].edited.mean()),
                      "masked_analysable_frac": float(f.masked_ok.mean())}
        sy = syn.get(rd, {})
        if sy:
            m["synthetic"] = {c: {"k": int(sum(v)), "n": len(v)} for c, v in sy.items()}
            cp, n1 = sy.get("E3_couplets_5pct"), sy.get("N1_linear_rr")
            m["iv"] = bool(cp and n1 and sum(cp) >= 50 and sum(n1) <= 5) if (cp and n1) else None
        if rd in segs and len(segs[rd]):
            s = segs[rd]
            sg = {"ii": A10.group_rates(s, "raw").get("CHF_minus_NSR"),
                  "iii_masked": A10.group_rates(s, "masked", "masked_ok").get("CHF_minus_NSR"),
                  "iii_wu": A10.group_rates(s, "wu", "wu_ok").get("CHF_minus_NSR")}
            sg["ii_holds"] = bool(sg["ii"] and sg["ii"]["ci95"][0] > 0)
            sg["iii_holds"] = bool(all(sg[k] and sg[k]["ci95"][0] <= 0 <= sg[k]["ci95"][1] for k in ("iii_masked", "iii_wu")))
            sg["raw_DR"] = {g: float(A10.subject_rates(s[s.group == g], "raw")[0].mean()) for g in GROUPS}
            m["segments"] = sg
        vv = ver.get(rd, {})
        m["verified"] = vv.get("VERIFIED")
        m["against_H"] = bool(m["verified"] and (m.get("i") is False or m.get("iii") is False))
        matrix[rd] = m
    out["matrix"] = matrix
    return out


# ----------------------------------------------------------------------------------------------- B1
def headline_windows(f):
    g = A10.gee(f, "TIT", ["x", "chf"])
    bt = A10.burden_table(f, ("TIT",))
    return {"n_pos": int(f.TIT.sum()), "OR": g["coef"]["x"]["OR"] if g.get("estimated") else None,
            "OR_ci95": g["coef"]["x"]["OR_ci95"] if g.get("estimated") else None,
            "CHF_OR": g["coef"]["chf"]["OR"] if g.get("estimated") else None,
            "CHF_OR_ci95": g["coef"]["chf"]["OR_ci95"] if g.get("estimated") else None,
            "burden_rates": {b["burden"]: b["TIT"]["rate"] for b in bt}}


def headline_segments(s):
    out = {}
    for arm, ok in (("raw", None), ("masked", "masked_ok"), ("wu", "wu_ok")):
        gr = A10.group_rates(s, arm, ok)
        out[arm] = {"NSR": gr["NSR"]["mean_subject_rate"], "CHF": gr["CHF"]["mean_subject_rate"],
                    "diff": gr["CHF_minus_NSR"]["diff"], "ci95": gr["CHF_minus_NSR"]["ci95"]}
    pc = A10.paired_change(s, "raw", "masked", "masked_ok")
    out["masked_change"] = {g: {"change": pc[g]["change"], "ci95": pc[g]["ci95"]} for g in pc}
    return out


def b1(phase):
    out = {}
    df10 = p10_real(phase)
    out["salt0_windows"] = headline_windows(df10)
    sg10 = p10_seg(phase)
    out["salt0_segments"] = headline_segments(sg10)
    W = load(phase, "seed_win")
    base = df10.set_index("wid")
    for s in sorted({r["salt"] for r in W}):
        rows = []
        for r in W:
            if r["salt"] != s:
                continue
            wid = r["id"].split("|", 1)[1].rsplit("|", 1)[0]
            rows.append({"wid": wid, "TIT_s": bool(r["raw"]["positive"]), "TITm_s": bool(r["masked"]["positive"]),
                         "TITe_s": bool(r["edited"]["positive"])})
        f = pd.DataFrame(rows).set_index("wid").join(base, how="inner")
        res = {"n": len(f), "flips_raw": int((f.TIT_s != f.TIT).sum()), "flips_masked": int((f.TITm_s != f.TITm).sum()),
               "flips_edited": int((f.TITe_s != f.TITe).sum())}
        res["headline"] = headline_windows(f.assign(TIT=f.TIT_s).reset_index())
        nl = {}
        for r in W:
            if r["salt"] == s and r["raw"]["positive"]:
                nl[r["id"].split("|", 1)[1].rsplit("|", 1)[0]] = r["raw"]["NL"]
        dnl = np.array([nl[w] - base.loc[w, "TIT_NL"] for w in nl if w in base.index and "TIT_NL" in base.columns])
        if len(dnl):
            res["raw_positive_NL_change"] = {"n": int(len(dnl)), "n_identical": int((np.abs(dnl) < 1e-12).sum()),
                                             "mean_abs": float(np.mean(np.abs(dnl))), "max_abs": float(np.max(np.abs(dnl)))}
        out[f"salt{s}_windows"] = res
    S = load(phase, "seed_seg")
    sgi = sg10.set_index(["subject", "seg"])
    for s in sorted({r["salt"] for r in S}):
        rows = []
        for r in S:
            if r["salt"] != s:
                continue
            for x in r["segments"]:
                rows.append({"subject": r["subject"], "group": r["group"], "seg": x["seg"],
                             "raw": bool(x["raw"]["positive"]), "masked": bool(x["masked"]["positive"]),
                             "masked_ok": bool(x["masked"].get("analysable")), "wu": bool(x["wu"]["positive"]),
                             "wu_ok": bool(x["wu"].get("eligible", True) and x["wu"].get("analysable", True))})
        f = pd.DataFrame(rows)
        j = f.set_index(["subject", "seg"]).join(sgi[["raw", "masked", "wu"]], rsuffix="_0", how="inner")
        out[f"salt{s}_segments"] = {"n": len(f), "flips_raw": int((j.raw != j.raw_0).sum()),
                                    "flips_masked": int((j.masked != j.masked_0).sum()),
                                    "flips_wu": int((j.wu != j.wu_0).sum()), "headline": headline_segments(f)}
    K3 = load(phase, "k3seed")
    for e in sorted({r["entropy"] for r in K3}):
        rr = [r for r in K3 if r["entropy"] == e]
        det = {r["id"].split("|", 1)[1].rsplit("|", 1)[0]: bool(r["K3"]["detected"]) for r in rr}
        f = base.join(pd.Series(det, name="K3_s"), how="inner")
        out[f"K3_entropy{e}"] = {"n": len(f), "detected": int(f.K3_s.sum()), "detected_seed0": int(f.K3.sum()),
                                 "flips": int((f.K3_s != f.K3).sum()),
                                 "windows": sorted(f.index[f.K3_s].tolist())}
    K1 = load(phase, "k1seed")
    for sd in sorted({tuple(r["seeds"]) for r in K1}):
        rr = [r for r in K1 if tuple(r["seeds"]) == sd]
        rows = {r["id"].split("|", 1)[1].rsplit("|", 1)[0]: (bool(r["K1"]["and_detected"]), bool(r["K1"]["lle_detected"]),
                                                             bool(r["K1"]["upo_detected"])) for r in rr}
        f = base.join(pd.DataFrame.from_dict(rows, orient="index", columns=["K1_s", "LLE_s", "UPO_s"]), how="inner")
        out[f"K1_seeds{sd[0]}_{sd[1]}"] = {"n": len(f), **{f"{k}_detected": int(f[f"{k}_s"].sum()) for k in ("K1", "LLE", "UPO")},
                                           **{f"{k}_seed0": int(f[k].sum()) for k in ("K1", "LLE", "UPO")},
                                           **{f"{k}_flips": int((f[f"{k}_s"] != f[k]).sum()) for k in ("K1", "LLE", "UPO")},
                                           "LLE_gee": headline_lle(f.assign(LLE=f.LLE_s).reset_index())}
    out["conclusions"] = b1_conclusions(out, headline_lle(df10))
    return out


def _side(ci, null):
    """'above' / 'includes' / 'below' the null value (the registered B1 'conclusion' of a CI)."""
    if ci is None or any(c is None or not np.isfinite(c) for c in ci):
        return None
    return "above" if ci[0] > null else ("below" if ci[1] < null else "includes")


def b1_conclusions(out, lle0):
    """B1 registered rule (PREREGISTRATION 7): a headline conclusion 'changes' when the CI's side of its null differs
    from salt 0 (Phase 10).  Written before the seed results were inspected."""
    items = []

    def add(name, ci0, null, others):
        s0 = _side(ci0, null)
        for lab, ci in others:
            si = _side(ci, null)
            items.append({"headline": name, "seed": lab, "salt0": s0, "seed_side": si, "changed": si != s0})

    w0 = out["salt0_windows"]
    wk = [k for k in out if k.startswith("salt") and k.endswith("_windows") and k != "salt0_windows"]
    add("P3 burden OR (windows)", w0["OR_ci95"], 1.0, [(k, out[k]["headline"]["OR_ci95"]) for k in wk])
    add("P3 CHF OR (windows)", w0["CHF_OR_ci95"], 1.0, [(k, out[k]["headline"]["CHF_OR_ci95"]) for k in wk])
    s0 = out["salt0_segments"]
    sk = [k for k in out if k.startswith("salt") and k.endswith("_segments") and k != "salt0_segments"]
    for arm in ("raw", "masked", "wu"):
        add(f"P2 CHF-NSR {arm} (segments)", s0[arm]["ci95"], 0.0, [(k, out[k]["headline"][arm]["ci95"]) for k in sk])
    for g in s0["masked_change"]:
        add(f"P2 masked change {g} (segments)", s0["masked_change"][g]["ci95"], 0.0,
            [(k, out[k]["headline"]["masked_change"].get(g, {}).get("ci95")) for k in sk])
    if lle0:
        kk = [k for k in out if k.startswith("K1_seeds")]
        add("LLE burden OR (windows)", lle0["OR_ci95"], 1.0, [(k, (out[k]["LLE_gee"] or {}).get("OR_ci95")) for k in kk])
    n_changed = sum(1 for i in items if i["changed"])
    return {"items": items, "n_changed": n_changed,
            "verdict": "weakens (a headline conclusion changes under a new seed)" if n_changed else
                       "supports (no headline conclusion changes)"}


def headline_lle(f):
    g = A10.gee(f, "LLE", ["x", "chf"])
    return {"OR": g["coef"]["x"]["OR"], "OR_ci95": g["coef"]["x"]["OR_ci95"]} if g.get("estimated") else None


# ----------------------------------------------------------------------------------------------- B2
def b2(phase):
    R = load(phase, "b2")
    rows = []
    for r in R:
        rows.append({"group": r["group"], "cond": f"{r['name']}|{r['ectopy']}|{r['level']}", "label": r["label"],
                     "family": r.get("family"), "sdnn128": r["K3_128"].get("sdnn_nn", np.nan),
                     "ana128": bool(r["K3_128"].get("analysable", True)),
                     "d360": bool(r["K3_360"]["detected"]), "d128": bool(r["K3_128"]["detected"]),
                     "deg128": bool(r["K3_128"].get("degenerate")), "deg360": bool(r["K3_360"].get("degenerate"))})
    df = pd.DataFrame(rows)
    out = {"n_windows": len(df)}
    if df.empty:
        return out
    fp = df[df.group.isin(["null_input", "null_real", "noncha_input", "noncha_real"])]
    cond = fp.groupby(["group", "cond"]).agg(n=("d128", "size"), k128=("d128", "sum"), k360=("d360", "sum"),
                                             deg128=("deg128", "sum")).reset_index()
    cond["rate128"] = cond.k128 / cond.n
    cond["rate360"] = cond.k360 / cond.n
    out["fp_conditions"] = int(len(cond))
    out["fp_fail_128"] = cond[cond.rate128 > 0.07][["group", "cond", "n", "k128"]].to_dict("records")
    out["fp_fail_360"] = cond[cond.rate360 > 0.07][["group", "cond", "n", "k360"]].to_dict("records")
    out["fp_pooled"] = {"n": len(fp), "k128": int(fp.d128.sum()), "k360": int(fp.d360.sum()),
                        "deg128": int(fp.deg128.sum()), "deg360": int(fp.deg360.sum())}
    out["fp_by_group"] = {g: {"n": int(len(d)), "k128": int(d.d128.sum()), "k360": int(d.d360.sum())}
                          for g, d in fp.groupby("group")}
    out["pass_128"] = bool(len(out["fp_fail_128"]) == 0)
    out["pass_360"] = bool(len(out["fp_fail_360"]) == 0)
    out["max_condition_k128"] = int(cond.k128.max()) if len(cond) else None
    pw = df[df.group.isin(["chaos_primary", "chaos_second"])]
    out["power"] = {g: {"n": int(len(d)), "k128": int(d.d128.sum()), "k360": int(d.d360.sum()),
                        "rate128": float(d.d128.mean()), "rate360": float(d.d360.mean())} for g, d in pw.groupby("group")}
    # descriptive breakdown (where the 1/128 s detections come from)
    det = fp[fp.d128]
    out["fp_detail"] = {"n_errors_generation": n_errors(phase, "b2"),
                        "det128_by_family": {k: int(v) for k, v in det.groupby("family").size().items()},
                        "det128_degenerate": int(det.deg128.sum()),
                        "non_ktz": {"n": int((fp.family != "ktz").sum()), "k128": int(fp[fp.family != "ktz"].d128.sum()),
                                          "k360": int(fp[fp.family != "ktz"].d360.sum())},
                        "ktz_sdnn_nn_ms_median_128": float(np.nanmedian(fp[fp.family == "ktz"].sdnn128)) * 1000,
                        "not_analysable_128": int((~fp.ana128).sum())}
    out["power_by_family"] = {f"{g}|{f}": {"n": int(len(d)), "k128": int(d.d128.sum()), "k360": int(d.d360.sum())}
                              for (g, f), d in pw.groupby(["group", "family"])}
    out["power_det128_degenerate"] = int(pw[pw.d128].deg128.sum())
    return out


def real_k3_detections(phase):
    """B2: degeneracy diagnostic on the Phase 10 real K3 detections."""
    from experiments.phase11_robustness import harness11 as H11
    from experiments.phase11_robustness import run11 as R11
    ctx = R11.Ctx(phase == "pilot")
    df = p10_real(phase)
    det = df[df.K3]
    out = []
    W = {w["wid"]: w for w in R11.windows(ctx)}
    for wid in det.wid:
        w = W[wid]
        d = R11.win_data(ctx, w)
        r = H11.k3_diag(d["rr"], d["lab"])
        out.append({"wid": wid, "sdnn_ms": 1000 * float(np.std(d["rr"])), "sdnn_nn_ms": 1000 * r["sdnn_nn"],
                    "mad_nn_ms": 1000 * r["mad_nn"], "zmax": r["zmax"], "G": r["G"], "detected": r["detected"],
                    "m_argmax": r["m_argmax"], "sd_at_argmax": r["sd_at_argmax"], "degenerate": r["degenerate"]})
    return out


# ----------------------------------------------------------------------------------------------- B3 / B4
def spike_rows(R, prefix):
    rows = []
    for r in R:
        for k, v in r["designs"].items():
            if k.startswith(prefix):
                rows.append({"subject": r["subject"], "group": r["group"], "family": r["family"], "param": r["param"],
                             "lam": r["lam"], "f": v["f"], "K3": bool(v["K3"]["detected"]),
                             "n_out": v.get("n_out_of_range", 0)})
    return pd.DataFrame(rows)


def exclusion_table(sp, U):
    """Phase 10 primary model-free exclusion bound (analysis10.exclusion 'design' and 'design_wilson_lower')."""
    out = {}
    for (fam, par), g in sp.groupby(["family", "param"]):
        fs = sorted(g.f.unique())
        res = {"family": fam, "param": par, "lam": float(g.lam.iloc[0]), "f": [float(x) for x in fs], "k": [], "n": []}
        pe, pw = [], []
        for f in fs:
            gf = g[g.f == f]
            k, n = int(gf.K3.sum()), len(gf)
            res["k"].append(k)
            res["n"].append(n)
            pe.append(k / n if n else 0.0)
            pw.append(A10.wilson(k, n)[0])
        for name, p in (("design", pe), ("wilson_lower", pw)):
            pi = [float(min(1.0, U / q)) if q > 0 else 1.0 for q in p]
            fm = {}
            for th in (0.05, 0.20):
                ok = [i for i in range(len(fs)) if all(v < th for v in pi[i:])]
                fm[f"{th:g}"] = float(fs[ok[0]]) if ok else None
            res[name] = {"pi_upper": pi, "f_min": fm}
        out[f"{fam}|{par}"] = res
    return out


def b3_b4(phase):
    from experiments.phase10_final import analysis10 as A
    an10 = json.load(open(P10 / p10_phase(phase) / "analysis" / "analysis.json"))
    U = an10["Q1"]["detectors"]["K3"]["real"]["U"]
    sp10 = A.spike_frame(A.load(p10_phase(phase), "spike"))
    sp10 = sp10[sp10.design == "add"]
    out = {"U_K3": U, "primary_phase10_recomputed": exclusion_table(sp10, U)}
    R3 = load(phase, "b3")
    if R3:
        sp = spike_rows(R3, "varpres_")
        out["B3_varpres"] = exclusion_table(sp, U)
        out["B3_out_of_range"] = {"n_series": len(sp), "n_with_out_of_range": int((sp.n_out > 0).sum()),
                                  "total_out_of_range": int(sp.n_out.sum())}
    R4 = load(phase, "b4")
    if R4:
        sp4 = spike_rows(R4, "add_")
        out["B4_second_window"] = exclusion_table(sp4, U)
        out["B4_pooled_166"] = exclusion_table(pd.concat([sp10[["subject", "group", "family", "param", "lam", "f", "K3"]],
                                                         sp4[["subject", "group", "family", "param", "lam", "f", "K3"]]]), U)
    return out


def b1_k3_bounds(phase):
    """NOT REGISTERED (descriptive addition, added after the B1 results were seen): Phase 10's real K3 bound U and the
    primary-design f_min recomputed with the K3 detections of each new K3 seed, because the K3 count changes."""
    from experiments.phase10_final import analysis10 as A
    df = p10_real(phase)
    K3 = load(phase, "k3seed")
    sp10 = A.spike_frame(A.load(p10_phase(phase), "spike"))
    sp10 = sp10[sp10.design == "add"]
    out = {"label": "descriptive, not registered"}
    for e in [0] + sorted({r["entropy"] for r in K3}):
        d = df.copy()
        if e:
            det = {r["id"].split("|", 1)[1].rsplit("|", 1)[0]: bool(r["K3"]["detected"]) for r in K3 if r["entropy"] == e}
            d["K3"] = d.wid.map(det).fillna(False).astype(bool)
        up = A.real_upper(d, "K3")
        ex = exclusion_table(sp10, up["U"])
        out[f"seed{e}"] = {"k": up["k"], "U": up["U"],
                           "f_min": {k: v["design"]["f_min"] for k, v in ex.items()}}
    return out


# ----------------------------------------------------------------------------------------------- main
SECTIONS = {"A1": a1, "A2_hrt": a2_hrt, "A2_transplant": a2_transplant, "A2_synth": a2_synth,
            "A2_extmask": a2_extmask, "A2_beat_type": a2_beat_type, "A3": a3, "A4": a4, "A5": a5, "B1": b1, "B2": b2,
            "B2_real_K3": real_k3_detections, "B3_B4": b3_b4,
            "B1_K3_bounds": b1_k3_bounds}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["conf", "pilot"], default="conf")
    ap.add_argument("--only", default=None)
    a = ap.parse_args(argv)
    outd = HERE / "results" / a.phase / "analysis"
    outd.mkdir(parents=True, exist_ok=True)
    p = outd / "analysis11.json"
    res = json.load(open(p)) if p.exists() else {}
    names = list(SECTIONS) if a.only is None else a.only.split(",")
    for n in names:
        try:
            res[n] = SECTIONS[n](a.phase)
            print(n, "ok", flush=True)
        except Exception as exc:                                # noqa: BLE001
            import traceback
            res[n] = {"error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc(limit=5)}
            print(n, "ERROR", exc, flush=True)
    res["errors_per_part"] = {pt: n_errors(a.phase, pt) for pt in ("sham_win", "sham_seg", "transplant", "synth1", "extmask",
                                                                   "rdg_verify", "rdg_win", "rdg_synth", "rdg_seg", "seed_win",
                                                                   "seed_seg", "k3seed", "k1seed", "b2", "b3", "b4", "hrt", "covar")}
    json.dump(res, open(p, "w"), indent=1, default=A10._json_default)



# ----------------------------------------------------------------------------------------------- C1 support
def c1_runs(phase):
    """C1: source search for the draft statement 'a single run of three premature beats was enough in 7 of 8
    windows': windows whose only abnormal beats form exactly one run of k >= 3 consecutive abnormal beats
    (symbol-based), with the Phase 10 raw titration outcome, by run length."""
    from experiments.phase10_final import data10 as D
    from experiments.phase11_robustness import run11 as R11
    ctx = R11.Ctx(phase == "pilot")
    df = p10_real(phase).set_index("wid")
    rows = []
    for w in R11.windows(ctx):
        if w["wid"] not in df.index:
            continue
        sym = R11.win_data(ctx, w)["sym"]
        ab = ~np.isin(sym, list(D.NORMAL))
        if not ab.any():
            continue
        runs, k, i = [], 0, 0
        for v in ab:
            if v:
                k += 1
            elif k:
                runs.append(k)
                k = 0
        if k:
            runs.append(k)
        if len(runs) == 1 and runs[0] >= 3:
            r = df.loc[w["wid"]]
            rows.append({"wid": w["wid"], "run_length": runs[0], "symbols": "".join(sym[ab]), "TIT": bool(r.TIT),
                         "TITm": bool(r.TITm), "TITm_ok": bool(r.TITm_analysable), "TITe": bool(r.TITe)})
    out = {"windows": rows, "n": len(rows), "n_positive_raw": int(sum(r["TIT"] for r in rows))}
    for L in sorted({r["run_length"] for r in rows}):
        rr = [r for r in rows if r["run_length"] == L]
        out[f"run_{L}"] = {"n": len(rr), "positive_raw": int(sum(r["TIT"] for r in rr))}
    return out


SECTIONS["C1_runs"] = c1_runs


if __name__ == "__main__":
    main()
