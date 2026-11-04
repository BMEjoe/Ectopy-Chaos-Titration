"""
Phase 11 C1: manuscript tables built ONLY from the analysis JSON (Phase 10 analysis.json, Phase 11 analysis11.json,
titration verification JSON) and from the frozen code constants (Table I).  Nothing is typed by hand.

    python -m experiments.phase11_robustness.paper_tables

Writes results/conf/analysis/paper_tables/<table>.md and .csv, and all_tables.md.
Table mapping (Phase 10 report order; the manuscript draft was not available to confirm the numbering):
  I detectors and parameters; II titration verification; III synthetic titration false positives (Q2e);
  IV P2 segment detection rates by arm; V Q3 dose-response; VI Q1 K3 exclusion; VII Q1 K1 exclusion;
  S-A1 .. S-B4 Phase 11 supplementary tables.
"""
from __future__ import annotations

import csv
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
P10 = HERE.parent / "phase10_final" / "results"
OUT = HERE / "results" / "conf" / "analysis" / "paper_tables"
FAM = {"henon": "Hénon", "logistic": "logistic", "coupled_vdp": "coupled vdP", "phase_reset": "phase-reset"}


def pct(x, d=1):
    return "–" if x is None else f"{100 * x:.{d}f}"


def cis(c, d=1, scale=100):
    return "–" if not c or c[0] is None else f"[{scale * c[0]:.{d}f}, {scale * c[1]:.{d}f}]"


def orci(o, c):
    return "–" if o is None else f"{o:.2f} [{c[0]:.2f}, {c[1]:.2f}]"


def write(name, title, header, rows, note=""):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"{name}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    md = [f"### {title}", "", "| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    md += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    if note:
        md += ["", note]
    (OUT / f"{name}.md").write_text("\n".join(md) + "\n")
    return "\n".join(md)


# --------------------------------------------------------------------------------------------- Phase 10 tables
def table_I():
    import final_pipeline as fp
    from experiments.phase10_final import titration10 as T
    from experiments.phase7_mitbih import detector as DT
    c, k = fp.CFG, DT.DETECTOR
    rows = [
        ["K1", "LLE test: IAAFT surrogates", k.lle_chaos_test_surrogates, "final_pipeline / phase7 detector.DETECTOR"],
        ["K1", "LLE test: alpha", k.lle_chaos_test_alpha, ""],
        ["K1", "LLE test: embedding m", k.lle_chaos_test_m, ""],
        ["K1", "UPO test: surrogates", k.so_surrogate_count, ""],
        ["K1", "keep_upo_on_short_lle_embedding", k.keep_upo_on_short_lle_embedding, ""],
        ["K1", "detrend method (D2)", fp.PHASE6_DETREND["rr_detrend"], "PHASE6_DETREND"],
        ["K1", "detrend gate: trend change / residual SD >=", fp.PHASE6_DETREND["rr_detrend_min_trend_sd"], "PHASE6_DETREND"],
        ["K3", "z threshold (masked_growth_z)", c.masked_growth_z, "PipelineConfig"],
        ["K3", "growth gate G (masked_growth_g)", c.masked_growth_g, ""],
        ["K3", "IAAFT surrogates", c.masked_growth_surrogates, ""],
        ["K3", "minimum valid delay vectors", c.masked_growth_min_valid, ""],
        ["TIT", "nonlinear memory kappa <=", T.KAPPA_MAX, "titration10"],
        ["TIT", "degree d <=", T.D_MAX, ""],
        ["TIT", "linear memory <=", T.KAPPA_LIN_MAX, ""],
        ["TIT", "test level", T.ALPHA, "variance-ratio F-test"],
        ["TIT", "noise realisations", T.N_REAL, "bisection, 10 halvings"],
        ["TIT", "eps power in C(r)", T.EPS_POWER, "0.5 = RMS"],
    ]
    return write("table_I", "Table I. Detectors and frozen parameters (read from the code)", ["method", "parameter", "value", "source"], rows)


def table_II():
    v = json.load(open(P10 / "verification" / "titration_verification.json"))["verification"]
    rows = [["V1 logistic r = 3.7", "~0.75", f"{v['V1']['NL_3.7']:.3f}", v["V1"]["pass"]],
            ["V2 logistic r = 3.575", "~0.09", f"{v['V2']['NL_3.575']:.3f}", v["V2"]["pass"]],
            ["V3 logistic r = 3.565 (periodic)", "~0", f"{v['V3']['NL_3.565']:.3f}", v["V3"]["pass"]],
            ["V4 bifurcation scan", "NL > 0 iff LE > 0",
             f"NL>0 {v['V4']['frac_NL_pos_chaotic'] * v['V4']['n_chaotic']:.0f}/{v['V4']['n_chaotic']}; NL=0 "
             f"{v['V4']['frac_NL_zero_periodic'] * v['V4']['n_periodic']:.0f}/{v['V4']['n_periodic']}; Spearman {v['V4']['spearman_NL_LE']:.2f}",
             v["V4"]["pass"]],
            ["V5 controls (20 each)", "NL = 0", "; ".join(f"{k} {c['positive']}/{c['n']} positive" for k, c in v["V5"]["controls"].items()),
             v["V5"]["pass"]],
            ["V6 flows + Hénon", "chaotic > 0, periodic 0",
             f"Lorenz28 {v['V6']['lorenz28']:.2f}, Lorenz160 {v['V6']['lorenz160']:.2f}, Hénon {v['V6']['henon']:.2f}, "
             f"MG chaotic {', '.join(f'{x:.2f}' for x in v['V6']['mg_chaotic_17_23_30'])}, MG periodic "
             f"{', '.join(f'{x:.2f}' for x in v['V6']['mg_periodic_14_15_16'])}", v["V6"]["pass"]]]
    return write("table_II", f"Table II. Titration verification (VERIFIED = {v['VERIFIED']})",
                 ["check", "published", "this implementation", "pass"], rows)


def table_III(A):
    Q = A["Q2e"]
    rows = []
    for key in sorted(Q, key=lambda s: (s.split("|")[0].startswith("vdp"), s)):
        c, n = key.split("|")
        v = Q[key]
        rows.append([c, n, v["n_windows"], v["raw"], v.get("raw_q128"), f"{v['masked']} ({v['masked_analysable']})",
                     f"{v['edited']} ({v['edited_eligible']})"])
    return write("table_III", "Table III. Titration false positives on validated non-chaotic synthetic windows (Q2e)",
                 ["condition", "intervals", "windows", "raw", "raw 1/128 s", "masked (analysable)", "edited (eligible)"], rows)


def table_IV(A):
    Q = A["Q2"]
    rows = []
    for key, name in (("P2_raw_by_group", "raw (P2)"), ("masked_by_group", "masked (analysable)"),
                      ("edited_by_group", "edited (eligible)"), ("wu_by_group", "Wu preprocessing")):
        g = Q[key]
        d = g["CHF_minus_NSR"]
        rows.append([name, f"{pct(g['NSR']['mean_subject_rate'])} {cis(g['NSR']['ci95'])}",
                     f"{pct(g['CHF']['mean_subject_rate'])} {cis(g['CHF']['ci95'])}", f"{pct(d['diff'])} {cis(d['ci95'])}",
                     f"{d['perm_p']:.4f}"])
    for key, name in (("P2_change_masked", "change masked"), ("change_edited", "change edited"), ("change_wu", "change Wu")):
        for g, v in Q[key].items():
            rows.append([f"{name} ({g}, {v['n_units']} segments)", f"raw {pct(v['rate_a'])}", f"arm {pct(v['rate_b'])}",
                         f"{pct(v['change'])} {cis(v['ci95'])}", f"removed {pct(v['a_pos_removed_frac'])} %"])
    return write("table_IV", "Table IV. Titration detection rate (%) on 12-min segments by arm (P2)",
                 ["arm", "NSR mean subject DR [95 % CI]", "CHF", "CHF − NSR", "permutation p / note"], rows)


def table_V(A):
    Q = A["Q3"]
    rows = []
    for c, nm in (("TIT", "titration (P3)"), ("K3", "K3 (P3)"), ("LLE", "LLE alone"), ("UPO", "UPO alone"), ("K1", "K1")):
        g = Q["P3_gee"][c]
        if g.get("estimated"):
            x, h = g["coef"]["x"], g["coef"]["chf"]
            rows.append([nm, g["n"], g["n_pos"], orci(x["OR"], x["OR_ci95"]), f"{x['p']:.2g}", orci(h["OR"], h["OR_ci95"])])
        else:
            rows.append([nm, g["n"], g["n_pos"], "not estimable", "", ""])
    t1 = write("table_V", "Table V. Ectopy dose-response (GEE, OR per doubling of 1 + burden)",
               ["outcome", "windows", "positive", "OR [95 % CI]", "p", "CHF OR adjusted"], rows)
    rows = []
    for b in Q["burden_table"]:
        rows.append([b["burden"], b["n"], b["n_subjects"]] + [f"{b[c]['k']} ({pct(b[c]['rate'], 0)} %)" for c in
                                                              ("TIT", "LLE", "UPO", "K1", "K3")])
    t2 = write("table_V_burden", "Table V (b). Positive windows by ectopy burden",
               ["burden", "windows", "subjects", "titration", "LLE", "UPO", "K1", "K3"], rows)
    return t1 + "\n\n" + t2


def _q1(A, det, name, title):
    D = A["Q1"]["detectors"][det]
    keys = sorted(D["curves"], key=lambda k: (list(FAM).index(D["curves"][k]["family"]), D["curves"][k]["lam"]))
    rows = []
    for k in keys:
        c, e = D["curves"][k], D["exclusion_L0"][k]
        fs = sorted(c["empirical"], key=float)
        rows.append([FAM[c["family"]], c["param"], f"{c['lam']:.3g}", " ".join(str(c["empirical"][f]["k"]) for f in fs),
                     c["replacement"]["k"], e["design"]["f_min"]["0.05"] or "none", e["design"]["f_min"]["0.2"] or "none",
                     e["design_wilson_lower"]["f_min"]["0.2"] or "none"])
    r = D["real"]
    note = (f"Real detections {r['k']}/{r['n_windows']} windows; U = {r['U']:.4f}. Detected counts are out of "
            f"{D['curves'][keys[0]]['empirical'][fs[0]]['n']} base windows at f = {', '.join(fs)}.")
    return write(name, title, ["family", "parameter", "λ", "detected at f grid", "replacement", "f_min π<0.05",
                               "f_min π<0.20", "f_min π<0.20 (Wilson lower)"], rows, note)


# --------------------------------------------------------------------------------------------- Phase 11 tables
def s_a1(B):
    rows = []
    A1 = B.get("A1", {})
    for cell in ("windows_masked", "windows_edited", "segments_masked", "segments_edited", "segments_wu"):
        c = A1.get(cell, {})
        for g in ("NSR", "CHF", "pooled"):
            v = c.get(g)
            if not v:
                continue
            rows.append([cell, g, v["n_units"], v["n_subjects"], pct(v["F_true"]), pct(v["F_sham"]),
                         f"{pct(v['D'])} {cis(v['D_ci95'])}", f"{v['Q']:.2f} {cis(v['Q_ci95'], 2, 1)}", v["verdict"]])
        cl = c.get("clean_removed_by_donor_mask")
        if cl:
            rows.append([cell, "clean units (donor masks)", c["clean_raw_positive_units"], cl["n_subjects"], "–",
                         pct(cl["fraction"]), "–", "–", f"CI {cis(cl['ci95'])}"])
    return write("S_A1", f"S-A1. True vs sham handling (raw-positive ectopic units; summary verdict: {A1.get('summary_verdict')})",
                 ["level / arm", "group", "units", "subjects", "removed by true (%)", "removed by sham (%)", "D [95 % CI] (points)",
                  "Q = sham/true", "verdict"], rows)


def s_a2(B):
    out = []
    T = B.get("A2_transplant", {})
    rows = []
    for v, d in T.get("delta_vs_V0", {}).items():
        for g in ("NSR", "CHF", "pooled"):
            x = d.get(g)
            if x:
                rows.append([v, g, x["n_units"], pct(x["rate_a"]), pct(x["rate_b"]), f"{pct(x['delta'])} {cis(x['ci95'])}",
                             f"{x['mcnemar']['a_only']}/{x['mcnemar']['b_only']}", f"{x['mcnemar']['p_exact']:.3g}"])
    out.append(write("S_A2_transplant", "S-A2a. One ectopic event transplanted into clean confirmatory windows (raw titration)",
                     ["variant", "group", "windows", "V0 rate (%)", "variant rate (%)", "Δ [95 % CI] (points)",
                      "discordant V0+/v+", "McNemar p"], rows))
    S = B.get("A2_synth", {})
    rows = []
    for c, d in S.items():
        rows.append([c, d["n"]] + [f"{d[v]['k']}" + (f" (Δ {pct(d[v]['delta'], 0)})" if v != "V0" else "")
                                   for v in ("V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7") if v in d])
    out.append(write("S_A2_synth", "S-A2b. Single events in synthetic non-chaotic windows (positives / 100)",
                     ["condition", "n", "V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7"], rows))
    H = B.get("A2_hrt", {})
    rows = []
    for g in ("NSR", "CHF"):
        h = H.get(g)
        if not h:
            continue
        to, ts = h.get("TO_ge5"), h.get("TS_ge5")
        rows.append([g, h["n_subjects"], h["n_subjects_with_V"], h["n_subjects_ge1_usable"], h["n_subjects_ge5_usable"],
                     h["n_usable_vpb"], "–" if not to else f"{to['median']:.2f} [{to['q25']:.2f}, {to['q75']:.2f}]",
                     "–" if not ts else f"{ts['median']:.2f} [{ts['q25']:.2f}, {ts['q75']:.2f}]",
                     pct(h.get("frac_TO_abnormal_ge5"), 0), pct(h.get("frac_TS_abnormal_ge5"), 0)])
    out.append(write("S_A2_hrt", "S-A2c. Heart rate turbulence (subjects with >= 5 usable VPBs: median [IQR])",
                     ["group", "subjects", "with V beats", ">= 1 usable VPB", ">= 5 usable", "usable VPBs", "TO (%)",
                      "TS (ms/RR)", "TO >= 0 (%)", "TS <= 2.5 (%)"], rows,
                     "Reference: healthy medians TS ~ 13 ms/RR, TO ~ -2.5 % (Grimm et al. 2003, as given; not verified)."))
    E = B.get("A2_extmask", {})
    rows = []
    for sub in ("all_units", "raw_positive_units"):
        for g, x in E.get(sub, {}).items():
            rows.append([sub, g, x["n_units"], pct(x["rate_a"]), pct(x["rate_b"]), f"{pct(x['delta'])} {cis(x['ci95'])}",
                         f"{x['mcnemar']['p_exact']:.3g}"])
    out.append(write("S_A2_extmask", "S-A2d. Standard vs extended (+20 post-ectopic intervals) masking, burden 1-4 windows",
                     ["units", "group", "windows", "standard (%)", "extended (%)", "Δ [95 % CI]", "McNemar p"], rows))
    BT = B.get("A2_beat_type", {})
    rows = [[c, v["n"], v["n_subjects"], v["k"], pct(v["rate"]), cis(v["wilson"])] for c, v in BT.items()]
    out.append(write("S_A2_beat_type", "S-A2e. Raw titration positive rate by the single abnormal beat's type",
                     ["category", "windows", "subjects", "positive", "rate (%)", "Wilson 95 % CI"], rows))
    return "\n\n".join(out)


def _coef_row(name, m):
    if not m or not m.get("estimated"):
        return [name, m.get("n") if m else "", m.get("n_pos") if m else "", "not estimable", "", ""]
    x = m["coef"]["x"]
    chf = m["coef"].get("chf")
    extra = (m.get("x_boot") or {}).get("OR_ci95")
    return [name, m["n"], m["n_pos"], orci(x["OR"], x["OR_ci95"]), f"{x['p']:.2g}",
            (orci(chf["OR"], chf["OR_ci95"]) if chf else "–") + (f"; boot CI [{extra[0]:.2f}, {extra[1]:.2f}]" if extra else "")]


def s_a3(B):
    A3 = B.get("A3", {})
    rows = []
    for lvl in ("windows", "segments"):
        L = A3.get(lvl, {})
        for oc in ("TIT", "LLE", "K1"):
            M = L.get(oc)
            if not M:
                continue
            for k, m in M.items():
                rows.append([f"{lvl} {oc}"] + _coef_row(k, m))
    return write("S_A3", f"S-A3. Confounding: burden OR per doubling (within-subject verdict: {A3.get('verdict')})",
                 ["level / outcome", "model", "units", "positive", "burden OR [95 % CI]", "p", "CHF OR / bootstrap CI"], rows)


def s_a4(B):
    A4 = B.get("A4", {})
    rows = []
    for arm in ("raw", "masked", "edited", "wu"):
        for g in ("NSR", "CHF"):
            v = A4.get(arm, {}).get(g)
            if not v:
                continue
            rows.append([arm, g, v["n_subjects"]] + [f"{v[k]['k']} ({pct(v[k]['frac'], 0)} %)" for k in
                                                       ("ge1_positive", "ge10pct", "ge25pct", "ge50pct")]
                        + ["–" if v.get("mean_NL_positive") is None else f"{v['mean_NL_positive']:.3f}"])
    return write("S_A4", "S-A4. Per-subject metrics (subjects with >= 1 / >= 10 / 25 / 50 % positive segments)",
                 ["arm", "group", "subjects", ">= 1 positive", ">= 10 %", ">= 25 %", ">= 50 %", "mean NL (positive)"], rows,
                 A4.get("pm97_pdf", ""))


def s_a5(B):
    A5 = B.get("A5", {})
    rows = []
    for rd, m in A5.get("matrix", {}).items():
        g = m.get("i_gee", {})
        ii = m.get("ii_windows") or {}
        po = m.get("iii_windows", {})
        sg = m.get("segments", {})
        rows.append([rd, m.get("verified"), orci(g.get("OR"), g.get("ci95")) if g.get("OR") else "–", m.get("i"),
                     f"{pct(ii.get('diff'))} {cis(ii.get('ci95'))}", m.get("ii"),
                     " / ".join(f"{a} {pct(po[a]['diff'])} {cis(po[a]['ci95'])}" for a in po if po[a]), m.get("iii"),
                     " ".join(f"{c}:{v['k']}" for c, v in m.get("synthetic", {}).items()), m.get("iv"),
                     (f"ii {sg.get('ii_holds')}, iii {sg.get('iii_holds')}" if sg else "–"), m.get("against_H")])
    return write("S_A5", "S-A5. Titration readings x headline claims", ["reading", "V1-V6 verified", "(i) burden OR", "(i) holds",
                 "(ii) raw CHF − NSR (windows)", "(ii) holds", "(iii) post-handling CHF − NSR (windows)", "(iii) holds",
                 "synthetic positives /100", "(iv) holds", "segments", "against H"], rows, A5.get("A5_1", ""))


def s_b(B):
    out = []
    B1 = B.get("B1", {})
    rows = []
    for k, v in B1.items():
        if k.startswith("salt") and k.endswith("_windows") and "headline" in v:
            h = v["headline"]
            rows.append([k, v["n"], h["n_pos"], orci(h["OR"], h["OR_ci95"]), orci(h["CHF_OR"], h["CHF_OR_ci95"]),
                         f"{v['flips_raw']}/{v['flips_masked']}/{v['flips_edited']}"])
        elif k == "salt0_windows":
            rows.append([k, "", v["n_pos"], orci(v["OR"], v["OR_ci95"]), orci(v["CHF_OR"], v["CHF_OR_ci95"]), "–"])
    out.append(write("S_B1_windows", "S-B1a. Titration seed robustness (windows)", ["seed", "windows", "raw positive", "burden OR",
                                                                                 "CHF OR adjusted", "flips raw/masked/edited"], rows))
    rows = []
    for k, v in B1.items():
        if k.endswith("_segments"):
            h = v.get("headline", v)
            rows.append([k, f"{pct(h['raw']['diff'])} {cis(h['raw']['ci95'])}", f"{pct(h['masked']['diff'])} {cis(h['masked']['ci95'])}",
                         f"{pct(h['wu']['diff'])} {cis(h['wu']['ci95'])}",
                         " ".join(f"{g} {pct(x['change'])}" for g, x in h["masked_change"].items()),
                         "–" if "flips_raw" not in v else f"{v['flips_raw']}/{v['flips_masked']}/{v['flips_wu']} of {v['n']}"])
    out.append(write("S_B1_segments", "S-B1b. Titration seed robustness (segments; CHF − NSR)",
                     ["seed", "raw", "masked", "Wu", "masked change", "flips raw/masked/Wu"], rows))
    rows = [[k, v["n"], v.get("detected", v.get("K1_detected")), v.get("detected_seed0", v.get("K1_seed0")),
             v.get("flips", v.get("K1_flips")),
             (f"LLE {v['LLE_detected']} (seed0 {v['LLE_seed0']}, flips {v['LLE_flips']}); UPO {v['UPO_detected']} "
              f"(seed0 {v['UPO_seed0']}, flips {v['UPO_flips']})") if "LLE_detected" in v else ", ".join(v.get("windows", []))]
            for k, v in B1.items() if k.startswith(("K3_", "K1_"))]
    out.append(write("S_B1_detectors", "S-B1c. K3 / K1 seed robustness", ["seed", "windows", "detected", "detected (Phase 10 seed)",
                                                                         "flips", "detail"], rows))
    B2 = B.get("B2", {})
    if B2.get("fp_by_group"):
        rows = [[g, v["n"], v["k360"], v["k128"]] for g, v in B2["fp_by_group"].items()]
        rows += [[f"power {g}", v["n"], f"{v['k360']} ({pct(v['rate360'])} %)", f"{v['k128']} ({pct(v['rate128'])} %)"]
                 for g, v in B2.get("power", {}).items()]
        out.append(write("S_B2", f"S-B2. K3 at 1/360 s vs re-quantized 1/128 s (pass 1/128: {B2.get('pass_128')}; "
                                 f"pass 1/360: {B2.get('pass_360')}; max detections in one FP condition at 1/128: "
                                 f"{B2.get('max_condition_k128')})", ["group", "windows", "K3 1/360 s", "K3 1/128 s"], rows))
    RK = B.get("B2_real_K3", [])
    if isinstance(RK, list) and RK:
        rows = [[r["wid"], f"{r['sdnn_ms']:.1f}", f"{r['sdnn_nn_ms']:.1f}", f"{r['zmax']:.3g}", f"{r['G']:.2f}", r["m_argmax"],
                 f"{r['sd_at_argmax']:.3g}", r["degenerate"]] for r in RK]
        out.append(write("S_B2_real", "S-B2b. Degeneracy diagnostic of the real K3 detections",
                         ["window", "SDNN (ms)", "SDNN of NN (ms)", "z_max", "G", "m at z_max", "surrogate SD at m", "degenerate"], rows))
    BB = B.get("B3_B4", {})
    rows = []
    for key in BB.get("primary_phase10_recomputed", {}):
        p = BB["primary_phase10_recomputed"][key]
        r = [FAM[p["family"]], p["param"], " ".join(map(str, p["k"])), p["design"]["f_min"]["0.05"] or "none",
             p["design"]["f_min"]["0.2"] or "none"]
        for sec in ("B3_varpres", "B4_pooled_166"):
            q = BB.get(sec, {}).get(key)
            r += ([" ".join(map(str, q["k"])), q["design"]["f_min"]["0.05"] or "none", q["design"]["f_min"]["0.2"] or "none",
                   q["wilson_lower"]["f_min"]["0.2"] or "none"] if q else ["–"] * 4)
        rows.append(r)
    out.append(write("S_B3_B4", "S-B3/B4. K3 exclusion: additive (Phase 10) vs variance-preserving (B3) vs 166 base windows (B4)",
                     ["family", "parameter", "add: detected", "add f_min .05", "add f_min .20", "varpres: detected", "vp f_min .05",
                      "vp f_min .20", "vp Wilson .20", "166: detected", "166 f_min .05", "166 f_min .20", "166 Wilson .20"], rows))
    return "\n\n".join(out)


def main():
    A = json.load(open(P10 / "conf" / "analysis" / "analysis.json"))
    p11 = HERE / "results" / "conf" / "analysis" / "analysis11.json"
    B = json.load(open(p11)) if p11.exists() else {}
    parts = ["# Manuscript tables (generated by paper_tables.py from the analysis JSON)", "", table_I(), table_II(),
             table_III(A), table_IV(A), table_V(A),
             _q1(A, "K3", "table_VI", "Table VI. K3 detection of additive chaos and exclusion bounds (Q1, P1)"),
             _q1(A, "K1", "table_VII", "Table VII. K1 detection of additive chaos and exclusion bounds (Q1, P1)"),
             "## Phase 11 supplementary tables", s_a1(B), s_a2(B), s_a3(B), s_a4(B), s_a5(B), s_b(B)]
    (OUT / "all_tables.md").write_text("\n\n".join(parts) + "\n")
    print("wrote", OUT / "all_tables.md")


if __name__ == "__main__":
    main()
