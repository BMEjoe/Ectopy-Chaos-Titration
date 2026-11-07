"""
Phase 12 Tier 2: manuscript tables built from results/conf/analysis/analysis12.json ONLY (markdown + CSV).

    PYTHONPATH=. python -m experiments.phase12_original.paper_tables12
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
A = json.load(open(HERE / "results" / "conf" / "analysis" / "analysis12.json"))
OUT = HERE / "results" / "conf" / "analysis" / "paper_tables12"
READ_ALL = ("R0", "R0-BP96F", "R0-BP96", "R0-S-FT", "R0-S-AAFT")
NAMES = {"ecological": "Ecological model", "logistic": "Logistic (r = 4)", "henon": "Hénon (1.4, 0.3)", "ikeda": "Ikeda (u = 0.9)",
         "gopy_sna": "Non-chaotic fractal (GOPY)", "rossler": "Rössler", "duffing": "Duffing", "lorenz": "Lorenz",
         "mackey_glass_100": "Mackey–Glass τ = 100", "mackey_glass_17": "Mackey–Glass τ = 17 (secondary)"}


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


def pp(d, key="diff"):
    return f"{100 * d[key]:+.1f} [{100 * d['ci95'][0]:+.1f}, {100 * d['ci95'][1]:+.1f}]" if d else "–"


def t_fidelity():
    rows = [[r["id"], r["point"], r["source"], r["verdict"], "yes" if r["decision_relevant"] else ""] for r in A["fidelity"]["rows"]]
    return write("T12_1_fidelity", "Table 12-1. Fidelity of the implementation (R0) to BP96 / PM97",
                 ["#", "point", "source location", "verdict", "decision-relevant"], rows,
                 "Full row text: experiments/phase12_original/FIDELITY.md.")


def t_b1():
    rows = []
    for rd, r in A["B1"]["readings"].items():
        rows.append([rd, f"{r['n_d_opt_gt1']}/20", f"{r['F_calc']['median']:.0f} [{r['F_calc']['q25']:.0f}, {r['F_calc']['q75']:.0f}]",
                     f"{r['F_crit']['median']:.3f}", f"{r['MW_Z']['median']:.2f}",
                     f"{r['loss_point']['median']:.0f} [{r['loss_point']['q25']:.0f}, {r['loss_point']['q75']:.0f}]",
                     "pass" if r["pass"] else "FAIL"])
    rows.append(["BP96 (Fig. 1)", "d_opt > 1", "9,387", "1.17", "−37.13", "~50", "–"])
    sa = A["B3"].get("series_A", {})
    note = ""
    if sa:
        r0 = sa["R0"]
        note = (f"Santa Fe series A (1,000 points): R0 F_calc {r0['F_calc_plain']:.2f} (BP96 10.23), F crit {r0['F_crit']:.2f} "
                f"(1.16), Mann–Whitney Z {r0['MW_Z']:.2f} (−24.60). Series E: not obtainable.")
    return write("T12_2_B1", "Table 12-2. BP96 Fig. 1 ecological model, 20 seeds (median [IQR])",
                 ["reading", "d_opt > 1", "F_calc", "F crit", "Mann–Whitney Z", "detection lost at noise (%)", "registered rule"],
                 rows, note)


def t_b2():
    rows = []
    for sy, v in A["B2"]["systems"].items():
        row = [NAMES.get(sy, sy), v["R0"]["table1"]]
        for rd in ("R0", "R0-BP96F", "R0-BP96"):
            x = v[rd]
            row.append(f"{x['median']:.0f} [{x['q25']:.0f}, {x['q75']:.0f}]" + ("" if x["consistent"] else " ✗"))
        row.append(f"{v['R0']['n_not_detected_noise_free']}/20")
        rows.append(row)
    return write("T12_3_B2", "Table 12-3. Maximum correlated noise (%) still detected, 1,000 points, 20 seeds (median [IQR]); "
                 "✗ = outside ±15 pp of BP96 Table 1 (approximate benchmark)",
                 ["system", "BP96 Table 1", "R0", "R0-BP96F", "R0-BP96", "R0 not detected noise-free"], rows,
                 "Mackey–Glass τ = 100 meets the ±15 pp rule only because Table 1 gives 15 %: it is never detected.")


def t_static():
    rows = []
    cells = A["B4_static"]["cells"]
    for phi in ("0.5", "0.9"):
        for n in ("1000", "512"):
            for tr in ("identity", "cube", "exp", "sqrt", "tanh", "N5"):
                rows.append([phi, n, tr] + [cells[f"{phi}|{n}|{tr}|{rd}"]["k"] for rd in READ_ALL])
    return write("T12_4_static", "Table 12-4. Static transforms of Gaussian AR(1): positive / 100", ["φ", "n", "transform", *READ_ALL], rows)


def t_synth_matrix():
    rows = [[c] + [v[rd]["k"] for rd in READ_ALL] for c, v in A["B4_synth"].items()]
    t1 = write("T12_5_synthetic", "Table 12-5. Phase 10 synthetic conditions (512 intervals, seeds 12000–12099): positive / 100",
               ["condition", *READ_ALL], rows)
    rows = []
    for rd, m in A["B4_matrix"]["readings"].items():
        iii = m["iii_windows"]
        rows.append([rd, "yes" if m["verified"] else "NO", f"{m['i_gee']['OR']:.2f} [{m['i_gee']['ci95'][0]:.2f}, {m['i_gee']['ci95'][1]:.2f}]",
                     pp(m["ii_windows"]), pp(iii["masked"]), pp(iii["edited"]), "yes" if m["iv"] else "no", pp(m["v_delta"], "delta"),
                     "AGAINST" if m["against_H"] else "–"])
    t2 = write("T12_6_claims", "Table 12-6. Claims matrix by reading (988 windows; transplant 439 windows)",
               ["reading", "V1–V6", "(i) burden OR", "(ii) raw CHF−NSR (pp)", "(iii) masked (pp)", "(iii) edited (pp)",
                "(iv) couplets vs N1", "(v) V1−V0 (pp)", "against H"], rows,
               f"N5 sentence (registered rule): {A['B4_matrix']['N5_sentence']['text']}")
    tp = A["B4_transplant"]["readings"]
    rows = [[rd] + [f"{100 * tp[rd]['rates'][v]:.1f}" for v in ("V0", "V1", "V5", "V6")] +
            [pp(tp[rd]["V1_minus_V0"]["pooled"], "delta"), pp(tp[rd]["V6_minus_V0"]["pooled"], "delta")] for rd in READ_ALL]
    t3 = write("T12_7_transplant", "Table 12-7. Transplant variants by reading (positive %, 439 clean windows)",
               ["reading", "V0", "V1 ventricular", "V5 atrial", "V6 couplet", "V1 − V0 (pp)", "V6 − V0 (pp)"], rows)
    return "\n\n".join([t1, t2, t3])


def t_pm97():
    P = A["PM97"]
    rows = [[s["subject"], s.get("group"), s.get("rule"), s.get("start_beat"), s.get("end_beat"), s.get("base_time"),
             f"{s.get('t_start_h', 0):.2f}", s.get("n_nonnormal"), s.get("n_out_of_range"),
             "yes" if s["subject"] in P["pm97_like_selection"] else ""] for s in sorted(P["subjects"], key=lambda s: s["subject"])]
    t1 = write("T12_8_pm97_blocks", "Table 12-8. PM97 replication: the 20,000-beat block of each subject",
               ["subject", "group", "rule", "start beat", "end beat", "header start time", "block start (h after record start)",
                "non-normal beats", "intervals outside 0.25–3 s", "PM97-like selection"], rows)
    rows = []
    for key, c in P["cells"].items():
        L, sub, arm, rd = key.split("|")
        if rd == "R0-BP96F":
            continue
        cmn = c["group_rates"].get("CHF_minus_NSR")
        rows.append([L, sub, arm, rd, f"{100 * c['NSR']['mean']:.1f} ({c['NSR']['n_subjects']})" if "NSR" in c else "–",
                     f"{100 * c['CHF']['mean']:.1f} ({c['CHF']['n_subjects']})" if "CHF" in c else "–", pp(cmn),
                     f"{c.get('mw_p', float('nan')):.3g}", f"{c.get('chf_20_80pp_below_nsr_mean')}/{c.get('n_chf')}",
                     "yes" if c.get("pattern", {}).get("reproduced") else "no"])
    t2 = write("T12_9_pm97", "Table 12-9. PM97 replication: detection rate (% of segments, mean of subjects) by arm and reading",
               ["segment", "subjects", "arm", "reading", "NSR % (n)", "CHF % (n)", "CHF − NSR (pp)", "Mann–Whitney p",
                "CHF 20–80 pp below NSR", "PM97 pattern"], rows)
    per = P["cells"]["500|all|raw|R0"]["per_subject"]
    perm = P["cells"]["500|all|masked|R0"]["per_subject"]
    perb = P["cells"]["500|all|raw|R0-BP96"]["per_subject"]
    rows = [[s, f"{100 * per[s]:.0f}", f"{100 * perb[s]:.0f}" if s in perb else "–", f"{100 * perm[s]:.0f}" if s in perm else "–"]
            for s in sorted(per)]
    t3 = write("T12_10_pm97_subjects", "Table 12-10. Per-subject detection rates (%; 40 × 500-beat segments; cf. PM97 Fig. 2b)",
               ["subject", "raw R0", "raw R0-BP96", "masked R0"], rows)
    cl, f3, cs = P["clogit_ectopy"], P["fig3"], P["exploratory_clean_segments"]
    rows = [["within-subject OR of ectopy (per doubling), conditional logit", f"{cl['coef']['x']['OR']:.2f} [{cl['coef']['x']['OR_ci95'][0]:.2f}, {cl['coef']['x']['OR_ci95'][1]:.2f}]", f"{cl['n']} segments, {cl['n_informative_subjects']} informative subjects"],
            ["Fig. 3 reading: mean SD of the first 80 nonlinear coefficients, NSR", f"{f3['NSR']['mean']:.2f} ± {f3['NSR']['se']:.2f} (s.e.)", f"{f3['NSR']['n']} subjects"],
            ["same, CHF", f"{f3['CHF']['mean']:.2f} ± {f3['CHF']['se']:.2f}", f"{f3['CHF']['n']} subjects"],
            ["same, CHF − NSR", f"{f3['CHF_minus_NSR']['diff']:.2f} [{f3['CHF_minus_NSR']['ci95'][0]:.2f}, {f3['CHF_minus_NSR']['ci95'][1]:.2f}], MW p {f3['CHF_minus_NSR']['mw_p']:.2f}", "descriptive"],
            ["EXPLORATORY: ectopy-free segments, CHF − NSR (R0)", pp(cs["R0"]["CHF_minus_NSR"]) + " pp", f"{cs['n_subjects']['NSR']} NSR / {cs['n_subjects']['CHF']} CHF subjects"]]
    t4 = write("T12_11_pm97_other", "Table 12-11. PM97 replication: ectopy association, coefficient variability, ectopy-free segments",
               ["quantity", "value", "n"], rows)
    return "\n\n".join([t1, t2, t3, t4])


def main():
    parts = ["# Phase 12 tables (generated by paper_tables12.py from analysis12.json)", "", t_fidelity(), t_b1(), t_b2(),
             t_static(), t_synth_matrix(), t_pm97()]
    (OUT / "all_tables12.md").write_text("\n\n".join(parts) + "\n")
    print("wrote", OUT / "all_tables12.md")


if __name__ == "__main__":
    main()
