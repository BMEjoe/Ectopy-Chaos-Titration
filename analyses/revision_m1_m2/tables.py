"""
Revision M1/M2 (EXPLORATORY): headline_by_reading.csv and claims_check.csv / .md from summary.json only
(PLAN.md Sections 4-5).

    PYTHONPATH=. python -m analyses.revision_m1_m2.tables
"""
from __future__ import annotations

import json

import pandas as pd

from analyses.revision_m1_m2 import analyse as AN

HEADS = ("R0", "R0-BP96", "R0-ORIG")
TABLE3 = [("N1_linear_rr", "Linear (N1), none"), ("S2_ectopic_5pct", "Isolated, 5 %"),
          ("S2_ectopic_10pct", "Isolated, 10 %"), ("E6_ectopic10_trend", "Isolated, 10 %, with trend"),
          ("E5_atrial_5pct", "Atrial-type, 5 %"), ("E5_atrial_10pct", "Atrial-type, 10 %"),
          ("E3_couplets_5pct", "Couplets, 5 %"), ("E3_couplets_10pct", "Couplets, 10 %"), ("E4_runs", "Runs"),
          ("E2_trigeminy", "Trigeminy"), ("E1_bigeminy", "Bigeminy"), ("N5_linear_rr_warped", "N5 static warp"),
          ("vdp_pooled_none", "Van der Pol (6 regimes), none"), ("vdp_pooled_S2_5", "Van der Pol, isolated 5 %"),
          ("vdp_pooled_E3_10", "Van der Pol, couplets 10 %"), ("vdp_pooled_E1", "Van der Pol, bigeminy")]


def S():
    return json.load(open(AN.SUMMARY))


def row(analysis, quantity, unit, rd, value, ci=None, n=None, note=""):
    ci = ci or [None, None]
    return {"analysis": analysis, "quantity": quantity, "unit": unit, "reading": rd, "value": value,
            "ci_lo": ci[0], "ci_hi": ci[1], "n": n, "note": note}


def headline_rows(s):
    H = s["headline"]["readings"]
    out = []
    for rd in HEADS:
        h = H.get(rd, {})
        if not h:
            continue
        for key, lab in TABLE3:
            v = h.get("synth", {}).get(key)
            if v:
                out.append(row("H-a synthetic (Table III)", f"{lab}: raw positive", "rate", rd, v["raw"]["rate"],
                               v["raw"]["ci"], v["raw"]["n"], f"{v['raw']['k']}/{v['raw']['n']}"))
                if v.get("edited", {}).get("n"):
                    out.append(row("H-a synthetic (Table III)", f"{lab}: edited positive", "rate", rd,
                                   v["edited"]["rate"], v["edited"]["ci"], v["edited"]["n"],
                                   f"{v['edited']['k']}/{v['edited']['n']}"))
        for b in h.get("dose_burden", []):
            out.append(row("H-b dose-response (988 windows)", f"burden {b['burden']}: raw positive", "rate", rd,
                           b["rate"], b["ci"], b["n"], f"{b['k']}/{b['n']} windows, {b['n_subjects']} subjects"))
        g = h.get("gee", {})
        if g.get("estimated"):
            c = g["coef"]["x"]
            out.append(row("H-b dose-response (988 windows)", "GEE OR per doubling of (1 + burden), adjusted for CHF",
                           "OR", rd, c["OR"], c["OR_ci95"], g["n"], g.get("working_correlation", "")))
            c = g["coef"]["chf"]
            out.append(row("H-b dose-response (988 windows)", "GEE CHF OR adjusted for burden", "OR", rd, c["OR"],
                           c["OR_ci95"], g["n"]))
        g = h.get("gee_chf_unadjusted", {})
        if g.get("estimated"):
            c = g["coef"]["chf"]
            out.append(row("H-b dose-response (988 windows)", "GEE CHF OR, unadjusted", "OR", rd, c["OR"], c["OR_ci95"],
                           g["n"], g.get("working_correlation", "")))
        cl = h.get("clogit", {})
        for key, lab in (("cov", "within-subject OR per doubling (conditional logit, 5 NN covariates + elapsed time)"),
                         ("x_only", "within-subject OR per doubling (conditional logit, burden only)")):
            c = cl.get(key, {})
            if c.get("estimated"):
                out.append(row("H-b dose-response (988 windows)", lab, "OR", rd, c["coef"]["x"]["OR"],
                               c["coef"]["x"]["OR_ci95"], c["n_informative_units"],
                               f"{c['n_informative_subjects']} informative subjects"))
        g = h.get("seg_gee", {})
        if g.get("estimated"):
            c = g["coef"]["x"]
            out.append(row("H-b dose-response (segments)", "segment GEE OR per doubling of (1 + burden)", "OR", rd,
                           c["OR"], c["OR_ci95"], g["n"]))
        a = h.get("AF_emp")
        if a:
            out.append(row("H-b attributable fraction", "AF = 1 - rate(burden 0)/rate(all), windows", "fraction", rd,
                           a["AF"], a["ci"], None, f"rate all {a['rate_all']:.3f}, burden 0 {a['rate_burden0']:.3f}"))
        a = h.get("AF_gee")
        if a:
            out.append(row("H-b attributable fraction", "GEE-based AF (predicted at burden 0)", "fraction", rd, a["AF"],
                           a["ci"], None, f"{a['n_boot_ok']} bootstrap refits"))
        if "ratio_1_vs_0" in h:
            out.append(row("H-b dose-response (988 windows)", "rate ratio, 1 ectopic beat vs 0", "ratio", rd,
                           h["ratio_1_vs_0"]))
        t = h.get("transplant", {})
        for v in ("V0", "V1", "V5", "V6"):
            if v in t.get("rates", {}):
                out.append(row("H-c single-beat transplant (439 windows)", f"{v} positive rate", "rate", rd,
                               t["rates"][v], None, 439))
        for v, lab in (("V1", "ventricular"), ("V5", "atrial-type"), ("V6", "couplet")):
            d = t.get(f"{v}_minus_V0", {}).get("pooled")
            if d:
                out.append(row("H-c single-beat transplant (439 windows)", f"{v} - V0 ({lab})", "pp", rd,
                               100 * d["delta"], [100 * d["ci95"][0], 100 * d["ci95"][1]], d["n_units"]))
        for lvl, key in (("windows", "win_groups"), ("segments", "seg_groups")):
            for arm, gr in h.get(key, {}).items():
                if "CHF_minus_NSR" not in gr:
                    continue
                c = gr["CHF_minus_NSR"]
                out.append(row(f"H-d group comparison ({lvl})", f"{arm}: CHF - NSR", "pp", rd, 100 * c["diff"],
                               [100 * c["ci95"][0], 100 * c["ci95"][1]], None,
                               f"NSR {100 * gr['NSR']['mean_subject_rate']:.1f} %, CHF {100 * gr['CHF']['mean_subject_rate']:.1f} %, "
                               f"perm p {c['perm_p']:.4f}; subjects {gr['NSR']['n_subjects']}/{gr['CHF']['n_subjects']}"))
        for key, st in h.get("pm97", {}).items():
            if key in ("clogit_ectopy", "clean_segments"):
                continue
            L, arm = key.split("|")
            gr = st["group_rates"]
            if "CHF_minus_NSR" in gr:
                c = gr["CHF_minus_NSR"]
                out.append(row(f"H-e PM97 replication ({L} beats)", f"{arm}: NSR rate", "rate", rd,
                               gr["NSR"]["mean_subject_rate"], gr["NSR"]["ci95"], gr["NSR"]["n_subjects"]))
                out.append(row(f"H-e PM97 replication ({L} beats)", f"{arm}: CHF rate", "rate", rd,
                               gr["CHF"]["mean_subject_rate"], gr["CHF"]["ci95"], gr["CHF"]["n_subjects"]))
                out.append(row(f"H-e PM97 replication ({L} beats)", f"{arm}: CHF - NSR", "pp", rd, 100 * c["diff"],
                               [100 * c["ci95"][0], 100 * c["ci95"][1]], None,
                               f"MW p {st.get('mw_p', float('nan')):.3f}; pattern reproduced "
                               f"{st.get('pattern', {}).get('reproduced')}"))
        c = h.get("pm97", {}).get("clogit_ectopy", {})
        if c.get("estimated"):
            out.append(row("H-e PM97 replication (500 beats)", "within-subject OR per doubling of segment ectopy", "OR",
                           rd, c["coef"]["x"]["OR"], c["coef"]["x"]["OR_ci95"], c["n_informative_units"]))
    m2 = s.get("m2")
    if m2:
        inf = pd.read_csv(AN.RES / "m2_influence.csv")
        d = inf[(inf.type == "ventricular") & (inf.placement == "random")]
        for r in d.itertuples():
            out.append(row("H-a synthetic isolated-beat dose (M2, 200 N1 windows)", f"k = {r.k} ventricular beats",
                           "rate", r.reading, r.pos_rate, [r.pos_ci_lo, r.pos_ci_hi], r.n, f"{r.n_pos}/{r.n}"))
    df = pd.DataFrame(out)
    df.to_csv(AN.RES / "headline_by_reading.csv", index=False, float_format="%.4g")
    return df


# ============================================================================================ claims
def fmt_ci(v, ci, scale=1.0, d=1):
    if v is None:
        return "n/a"
    if ci is None or ci[0] is None:
        return f"{scale * v:.{d}f}"
    return f"{scale * v:.{d}f} [{scale * ci[0]:.{d}f}, {scale * ci[1]:.{d}f}]"


def claims(s):
    H = s["headline"]["readings"]
    ST = s.get("stored_readings", {})
    rows = []

    def add(cid, where, text, vals, holds, wording, criterion):
        rows.append({"id": cid, "where": where, "claim": text, **{f"value_{k}": v for k, v in vals.items()},
                     "criterion": criterion, "holds_under_every_reading": holds, "suggested_wording": wording})

    def g(rd, *path):
        x = H.get(rd, {})
        for p in path:
            if isinstance(x, list):
                x = x[p]
            else:
                x = (x or {}).get(p)
        return x

    rds = [r for r in HEADS if r in H]
    # A0 / C8: "mainly"
    af = {rd: g(rd, "AF_emp") for rd in rds}
    vals = {rd: fmt_ci(a["AF"], a["ci"], 100, 0) + " %" for rd, a in af.items() if a}
    minaf = min(a["AF"] for a in af.values() if a)
    holds = all(a["AF"] > 0.5 for a in af.values() if a)
    lo_af = min(a["AF"] for a in af.values() if a)
    hi_af = max(a["AF"] for a in af.values() if a)
    add("A0", "abstract; conclusions (C8 'most often caused by ectopic beats')",
        "In heart rate, titration mainly detects ectopic beats", {**vals}, holds,
        None if holds else (f"In heart rate, ectopic beats are a major source of titration positives: they account for "
                            f"about {100 * lo_af:.0f}-{100 * hi_af:.0f} % of positive windows, depending on the test's "
                            f"conventions, and burden, not group, predicts a positive window"),
        "attributable fraction AF > 0.50 under every reading (PLAN.md 5)")
    # L1: nearly tripled
    vals = {rd: f"{g(rd, 'ratio_1_vs_0'):.2f}x" for rd in rds if g(rd, "ratio_1_vs_0")}
    rr = [g(rd, "ratio_1_vs_0") for rd in rds if g(rd, "ratio_1_vs_0")]
    add("L1", "lead paragraph", "One ectopic beat nearly tripled the chance of titration reporting chaos", vals,
        all(x >= 2.5 for x in rr), f"One ectopic beat raised the chance of titration reporting chaos {min(rr):.1f}- to "
        f"{max(rr):.1f}-fold, depending on the test's conventions", "rate ratio (burden 1 vs 0) >= 2.5 under every reading")
    # burden, not group (Sec. V B)
    vals, ok = {}, True
    for rd in rds:
        ge = g(rd, "gee")
        if ge and ge.get("estimated"):
            c = ge["coef"]["chf"]
            vals[rd] = fmt_ci(c["OR"], c["OR_ci95"], 1, 2)
            ok &= c["OR_ci95"][0] <= 1 <= c["OR_ci95"][1]
    add("V-B", "Sec. V B; Discussion", "With burden in the model there was no separate effect of CHF (OR 0.77); ectopy "
        "accounts for the group difference at the window level", vals, ok, None, "adjusted CHF OR CI includes 1 under every "
        "reading")
    # isolated ectopy rarely triggers titration (Sec. IV A, Table III)
    vals, ok = {}, True
    for rd in rds:
        v5, v10 = g(rd, "synth", "S2_ectopic_5pct", "raw"), g(rd, "synth", "S2_ectopic_10pct", "raw")
        a5 = g(rd, "synth", "E5_atrial_5pct", "raw")
        if v5:
            vals[rd] = f"isolated 5 %: {v5['k']}/100; 10 %: {v10['k']}/100; atrial 5 %: {a5['k']}/100"
            ok &= v5["k"] <= 10 and a5["k"] <= 10
    add("IV-A", "Sec. IV A; Table III", "Isolated beats (5 %, 10 %) and atrial-type beats rarely triggered titration",
        vals, ok, "Under my implementation, isolated beats at 5-10 % density rarely triggered titration; under the original "
        "test's conventions, 5 % isolated and atrial-type beats made most windows positive", "isolated 5 % and atrial 5 % "
        "<= 10/100 under every reading")
    # A1: couplets 100/100 every version
    vals = {}
    allh = True
    for rd in rds:
        v = g(rd, "synth", "E3_couplets_5pct", "raw")
        vals[rd] = f"{v['k']}/{v['n']}" if v else "n/a"
        allh &= bool(v and v["k"] == v["n"])
    for rd in ("R0-BP96F", "R0-S-FT", "R0-S-AAFT"):
        v = ST.get(rd, {}).get("couplets5")
        if v:
            vals[rd + " (Phase 12, other seeds)"] = f"{v['k']}/{v['n']}"
            allh &= v["k"] == v["n"]
    add("A1", "abstract", "For all 100 linear stochastic RR series with ectopic couplets, every titration version tested "
        "positive", vals, allh, None, "couplets 5 %: 100/100 under every reading")
    # A2: 22 % -> 63 %
    vals, ok = {}, True
    for rd in rds:
        b = g(rd, "dose_burden")
        if b:
            vals[rd] = f"{100 * b[0]['rate']:.0f} % -> {100 * b[1]['rate']:.0f} %"
            ok &= b[1]["ci"][0] > b[0]["ci"][1]
    add("A2", "abstract; lead", "In 988 real windows, positive rates increased from 22 % (no ectopic beats) to 63 % "
        "(1 ectopic beat)", vals, ok, "positive rates rose from 22-37 % (no ectopic beats) to 50-63 % (1 ectopic beat), "
        "depending on the test's conventions", "rate at burden 1 above burden 0 (Wilson CIs disjoint); the quoted "
        "numbers are R0's")
    # A3 / C4: single beat 28-64 pp
    vals, lo, hi = {}, [], []
    for rd in rds:
        d = g(rd, "transplant", "V1_minus_V0", "pooled")
        if d:
            vals[rd] = fmt_ci(d["delta"], d["ci95"], 100) + " pp"
            lo.append(d["ci95"][0])
            hi.append(d["delta"])
    for rd in ("R0-BP96F", "R0-S-FT", "R0-S-AAFT"):
        d = ST.get(rd, {}).get("transplant_V1")
        if d:
            vals[rd + " (Phase 12)"] = fmt_ci(d["delta"], d["ci95"], 100) + " pp"
            lo.append(d["ci95"][0])
            hi.append(d["delta"])
    ests = [float(v.split()[0]) for v in vals.values()]
    holds = all(x > 0 for x in lo) and min(ests) >= 27.5 and max(ests) <= 64.5
    add("A3", "abstract; conclusions", "One premature beat inserted into clean windows added 28-64 percentage points, "
        "depending on conventions", vals, holds, f"added {min(ests):.0f}-{max(ests):.0f} percentage points, depending on "
        "the test's conventions" if not holds else None, "V1 - V0 CI > 0 under every reading and the range covers every "
        "reading's estimate")
    # A4 / C5 / L2: segments raw +22.4, after handling -4.9 to +5.4 n.s.
    vals, ok_raw, ok_post = {}, True, True
    ranges = []
    for rd in rds:
        sg = g(rd, "seg_groups")
        if not sg:
            continue
        raw = sg["raw"]["CHF_minus_NSR"]
        post = [sg[a]["CHF_minus_NSR"] for a in ("masked", "edited", "wu")]
        vals[rd] = (f"raw {fmt_ci(raw['diff'], raw['ci95'], 100)}; after handling "
                    + ", ".join(fmt_ci(p["diff"], p["ci95"], 100) for p in post))
        ok_raw &= raw["ci95"][0] > 0
        ok_post &= all(p["ci95"][0] <= 0 <= p["ci95"][1] for p in post)
        ranges += [100 * p["diff"] for p in post]
    add("A4", "abstract; lead; conclusions", "Raw CHF segments were positive 22.4 points more often than healthy ones, but "
        "-4.9 to +5.4 points (not significant) after ectopy handling", vals, bool(ok_raw and ok_post),
        f"holds under every reading; across readings the handled differences span {min(ranges):+.1f} to {max(ranges):+.1f} "
        "points (all not significant)", "raw CI > 0 and every handled-arm CI includes 0, under every reading")
    # windows raw excess (Sec. V A) incl. stored AAFT
    vals, ok = {}, True
    for rd in rds:
        c = g(rd, "win_groups", "raw", "CHF_minus_NSR")
        if c:
            vals[rd] = fmt_ci(c["diff"], c["ci95"], 100) + " pp"
            ok &= c["ci95"][0] > 0
    for rd in ("R0-BP96F", "R0-S-FT", "R0-S-AAFT"):
        d = ST.get(rd, {})
        if d.get("win_raw_diff") is not None:
            vals[rd + " (Phase 12)"] = fmt_ci(d["win_raw_diff"], d["win_raw_ci"], 100) + " pp"
            ok &= d["win_raw_ci"][0] > 0
    add("V-A", "Sec. V A; Discussion", "Raw window-level CHF excess +23.5 (R0), +14.2 (original conventions), +4.2 "
        "(AAFT surrogates)", vals, ok, "keep the manuscript's wording: the raw excess is reading-dependent and disappears "
        "with amplitude-adjusted surrogates", "raw window CHF - NSR CI > 0 under every reading (known to fail for AAFT)")
    # A5 / L3 / C6: PM97
    vals, okr, okh = {}, True, True
    pm_def = []
    for rd in rds:
        p = g(rd, "pm97")
        if not p:
            continue
        raw = p["500|raw"]["group_rates"]["CHF_minus_NSR"]
        post = [p[f"500|{a}"]["group_rates"]["CHF_minus_NSR"] for a in ("masked", "edited", "wu")]
        vals[rd] = (f"NSR raw {100 * p['500|raw']['group_rates']['NSR']['mean_subject_rate']:.1f} %; raw diff "
                    f"{fmt_ci(raw['diff'], raw['ci95'], 100)}; handled " + ", ".join(fmt_ci(x["diff"], x["ci95"], 100) for x in post))
        okr &= not p["500|raw"].get("pattern", {}).get("reproduced", False)
        pm_def += [-100 * x["diff"] for x in post]
        okh &= all(x["ci95"][1] < 0 for x in post)
    add("A5", "abstract; lead; conclusions", "On 2 older databases, Poon and Merrill's procedure did not reproduce their "
        "lower positive rate in CHF on raw intervals; it appeared only after ectopy handling (26-46 points lower)",
        vals, bool(okr and okh), f"holds under every reading; after handling CHF was {min(pm_def):.0f}-{max(pm_def):.0f} "
        "points lower across readings", "raw PM97 pattern not reproduced and every handled-arm CI < 0, under every reading")
    # A7: static transforms (N5)
    vals, ok = {}, True
    for rd in rds:
        v = g(rd, "synth", "N5_linear_rr_warped", "raw")
        if v:
            vals[rd] = f"{v['k']}/{v['n']}"
            ok &= v["k"] >= 50
    for rd in ("R0-S-FT", "R0-S-AAFT"):
        v = ST.get(rd, {}).get("N5")
        if v:
            vals[rd + " (Phase 12)"] = f"{v['k']}/{v['n']}"
            ok &= v["k"] >= 50
    add("A7", "abstract; Discussion", "Titration also detects some static transforms of linear noise (N5 100/100)",
        vals, ok, None, "N5 >= 50/100 under every reading")
    # C3: burden steep; within subjects
    vals, ok = {}, True
    for rd in rds:
        ge, cl = g(rd, "gee"), g(rd, "clogit", "cov")
        if ge and ge.get("estimated") and cl and cl.get("estimated"):
            vals[rd] = (f"GEE {fmt_ci(ge['coef']['x']['OR'], ge['coef']['x']['OR_ci95'], 1, 2)}; within-subject "
                        f"{fmt_ci(cl['coef']['x']['OR'], cl['coef']['x']['OR_ci95'], 1, 2)}")
            ok &= ge["coef"]["x"]["OR_ci95"][0] > 1 and cl["coef"]["x"]["OR_ci95"][0] > 1
    add("C3", "conclusions; Sec. V B", "Positives increased steeply with ectopic burden (GEE OR 2.31); the same within "
        "subjects (OR 3.40)", vals, ok, None, "both OR CIs > 1 under every reading")
    # C2: few isolated beats under R0
    m2 = s.get("m2")
    if m2:
        inf = pd.read_csv(AN.RES / "m2_influence.csv")
        d = inf[(inf.type == "ventricular") & (inf.placement == "random")]
        vals = {}
        for rd in ("R0", "R0-BP96", "R0-ORIG"):
            x = d[d.reading == rd].set_index("k")
            if len(x):
                vals[rd] = ", ".join(f"k={k}: {100 * x.loc[k, 'pos_rate']:.0f} %" for k in (0, 1, 2, 4, 8, 16, 26, 51))
        add("C2", "conclusions; Sec. IV A", "1-8 isolated premature beats made 91-99 % of linear windows positive; 26 and "
            "51 beats 5 % and 0 %", vals, None, None, "descriptive (exploratory dose-response)")
    df = pd.DataFrame(rows)
    df.to_csv(AN.RES / "claims_check.csv", index=False)
    vcols = [c for c in df.columns if c.startswith("value_")]
    md = ["# Claims check (EXPLORATORY, post hoc)", "",
          "One row per quantitative claim in the abstract, lead paragraph and conclusions (plus the Results statements "
          "they rest on). Values per reading; 'holds' applies the stated criterion to every reading listed. Phase 12 "
          "values are stored results on that phase's seeds.", ""]
    for _, r in df.iterrows():
        md.append(f"## {r['id']} ({r['where']})")
        md.append(f"**Claim.** {r['claim']}")
        md.append("")
        for c in vcols:
            if isinstance(r[c], str):
                md.append(f"- {c[6:]}: {r[c]}")
        hv = r["holds_under_every_reading"]
        md.append(f"- **Criterion:** {r['criterion']}. **Holds under every reading:** "
                  f"{'yes' if hv is True else ('no' if hv is False else 'descriptive')}")
        if isinstance(r["suggested_wording"], str):
            md.append(f"- **Suggested wording:** {r['suggested_wording']}")
        md.append("")
    (AN.RES / "claims_check.md").write_text("\n".join(md) + "\n")
    return df


def main():
    s = S()
    headline_rows(s)
    claims(s)


if __name__ == "__main__":
    main()
