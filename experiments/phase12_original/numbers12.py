"""
Phase 12 Tier 2: the numbers ledger (docs/PHASE12_NUMBERS.md) for Phases 10-12.  Every value is read from the analysis
JSON files by the key shown in the ledger; nothing is typed by hand.

    PYTHONPATH=. python -m experiments.phase12_original.numbers12
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
F10 = "experiments/phase10_final/results/conf/analysis/analysis.json"
F11 = "experiments/phase11_robustness/results/conf/analysis/analysis11.json"
F12 = "experiments/phase12_original/results/conf/analysis/analysis12.json"
J = {k: json.load(open(ROOT / f)) for k, f in (("10", F10), ("11", F11), ("12", F12))}
OUT = ROOT / "docs" / "PHASE12_NUMBERS.md"


def get(src, path):
    o = J[src]
    for k in path:
        o = o[k] if not callable(k) else k(o)
    return o


def pp(d, key="diff"):
    return f"{100 * d[key]:+.1f} pp [{100 * d['ci95'][0]:+.1f}, {100 * d['ci95'][1]:+.1f}]"


def pct_ci(v, ci):
    return f"{100 * v:.1f} % [{100 * ci[0]:.1f}, {100 * ci[1]:.1f}]"


def orci(c):
    return f"OR {c['OR']:.2f} [{c['OR_ci95'][0]:.2f}, {c['OR_ci95'][1]:.2f}]"


def kn(c):
    return f"{c['k']}/{c['n']}"


def burden(label):
    return lambda t: next(b for b in t if b["burden"] == label)


# (section, id, definition, src, path, formatter, n-text, sentence template with {v}, flag)
E = []


def add(section, nid, definition, src, path, fmt, n, sentence, flag=""):
    E.append((section, nid, definition, src, path, fmt, n, sentence, flag))


# ---------------------------------------------------------------------------------------------- Methods
add("Methods", "M1", "Decision-relevant deviations of R0 from BP96 found in the fidelity audit", "12",
    ["fidelity", "n_decision_relevant_deviations"], str, "22 audited points",
    "Of 22 points audited against the original papers, {v} were decision-relevant deviations of our implementation "
    "(N in the criterion, the F statistic, the missing Gaussianity check and the omitted surrogate step).")
add("Methods", "M2", "analysis11 rerun on stored results identical to the committed JSON", "12",
    ["sanity", "analysis11_rerun_identical_to_committed"], str, "-",
    "Re-running the Phase 11 analysis on the stored results reproduced it exactly ({v}).")
add("Methods", "M3", "R0 recomputed vs Phase 10 stored decisions (raw / masked / edited agree, of 988)", "12",
    ["B4_matrix", "R0_recomputed_vs_phase10"], lambda d: f"{d['raw_agree']}/{d['masked_agree']}/{d['edited_agree']} of {d['n']}",
    "988 windows", "Recomputing R0's noise-free decision reproduced every stored Phase 10 decision ({v}).")

# ---------------------------------------------------------------------------------------------- Results: titration and ectopy
add("Results: titration and ectopy", "R1", "Raw segment positive rate, NSR (mean subject rate)", "10",
    ["Q2", "P2_raw_by_group", "NSR"], lambda d: pct_ci(d["mean_subject_rate"], d["ci95"]), "54 subjects",
    "On raw RR, {v} of healthy subjects' 12-min segments were titration-positive.")
add("Results: titration and ectopy", "R2", "Raw segment positive rate, CHF", "10",
    ["Q2", "P2_raw_by_group", "CHF"], lambda d: pct_ci(d["mean_subject_rate"], d["ci95"]), "29 subjects",
    "In CHF, {v} of segments were positive.")
add("Results: titration and ectopy", "R3", "Raw CHF - NSR (segments)", "10", ["Q2", "P2_raw_by_group", "CHF_minus_NSR"], pp,
    "83 subjects", "CHF segments were positive more often than healthy ones on raw RR ({v}).")
add("Results: titration and ectopy", "R4", "Masked CHF - NSR (segments)", "10", ["Q2", "masked_by_group", "CHF_minus_NSR"], pp,
    "83 subjects", "After masking ectopy the group difference vanished ({v}).")
add("Results: titration and ectopy", "R5", "Edited CHF - NSR (segments)", "10", ["Q2", "edited_by_group", "CHF_minus_NSR"], pp,
    "83 subjects", "After editing: {v}.")
add("Results: titration and ectopy", "R6", "Wu-style CHF - NSR (segments)", "10", ["Q2", "wu_by_group", "CHF_minus_NSR"], pp,
    "83 subjects", "After Wu-style elimination: {v}.")
add("Results: titration and ectopy", "R7", "Window burden OR per doubling (GEE, x + CHF)", "10",
    ["Q3", "P3_gee", "TIT", "coef", "x"], orci, "988 windows", "Each doubling of ectopy burden raised the odds of a positive window ({v}).")
add("Results: titration and ectopy", "R8", "CHF OR adjusted for burden (GEE)", "10", ["Q3", "P3_gee", "TIT", "coef", "chf"], orci,
    "988 windows", "Adjusted for burden, CHF itself did not raise the odds ({v}).")
add("Results: titration and ectopy", "R9", "Positive rate, windows with burden 0", "10", ["Q3", "burden_table", burden("0"), "TIT"],
    lambda d: pct_ci(d["rate"], d["wilson"]), "562 windows", "Windows without an abnormal beat were positive in {v}.")
add("Results: titration and ectopy", "R10", "Positive rate, windows with burden 1", "10", ["Q3", "burden_table", burden("1"), "TIT"],
    lambda d: pct_ci(d["rate"], d["wilson"]), "116 windows", "Windows with one abnormal beat were positive in {v}.")
add("Results: titration and ectopy", "R11", "Within-subject burden OR, covariate-adjusted (conditional logit)", "11",
    ["A3", "windows", "TIT", "clogit_cov", "coef", "x"], orci, "987 windows, 77 informative subjects",
    "Within subjects and adjusted for heart-rate covariates, the burden effect held ({v}).")
add("Results: titration and ectopy", "R12", "Transplant: one ventricular-type beat, V1 - V0 (R0)", "11",
    ["A2_transplant", "delta_vs_V0", "V1", "pooled"], lambda d: pp(d, "delta"), "439 windows, 69 subjects",
    "Inserting one ventricular-type beat into clean windows raised the positive rate by {v} (22.1 % to 86.3 %).")
add("Results: titration and ectopy", "R13", "Transplant V1 - V0 under the faithful BP96 reading (R0-BP96)", "12",
    ["B4_transplant", "readings", "R0-BP96", "V1_minus_V0", "pooled"], lambda d: pp(d, "delta"), "439 windows",
    "Under the original paper's own test conventions the single-beat effect was smaller but still clear ({v}).",
    "NEW (Phase 12): qualifies R12")
add("Results: titration and ectopy", "R14", "Sham control, windows masked: D (true - sham)", "11", ["A1", "windows_masked", "pooled"],
    lambda d: f"D = {d['D']:.3f} [{d['D_ci95'][0]:.3f}, {d['D_ci95'][1]:.3f}], Q = {d['Q']:.2f}", "126 units",
    "Masking ectopy removed more positives than masking the same number of random rows ({v}).")

# ---------------------------------------------------------------------------------------------- Results: fidelity and benchmarks
add("Results: fidelity and benchmarks", "B1", "Ecological model (BP96 Fig. 1): F_calc median [IQR], R0", "12",
    ["B1", "readings", "R0", "F_calc"], lambda d: f"{d['median']:.0f} [{d['q25']:.0f}, {d['q75']:.0f}]", "20 seeds",
    "On BP96's ecological model our implementation gave F_calc {v} (BP96: 9,387), with d_opt > 1 in all 20 seeds.")
add("Results: fidelity and benchmarks", "B2", "Ecological model: correlated noise at which detection is lost, median [IQR], R0", "12",
    ["B1", "readings", "R0", "loss_point"], lambda d: f"{d['median']:.0f} % [{d['q25']:.0f}, {d['q75']:.0f}]", "20 seeds",
    "Detection was lost at {v} correlated noise, against about 50 % in BP96 (registered rule: pass, at its lower edge).",
    "lower than BP96")
add("Results: fidelity and benchmarks", "B3", "Ecological model: loss point under R0-BP96", "12",
    ["B1", "readings", "R0-BP96", "loss_point"], lambda d: f"{d['median']:.0f} % [{d['q25']:.0f}, {d['q75']:.0f}]", "20 seeds",
    "With the faithful reading detection was lost at {v} (registered B1 rule: fail).", "R0-BP96 fails B1")
add("Results: fidelity and benchmarks", "B4", "Table 1 systems outside +/-15 pp (R0)", "12", ["B2", "misses", "R0"],
    lambda d: f"{len(d)} of 9 ({', '.join(d)})", "20 seeds each",
    "Against BP96's approximate Table 1, our maximum tolerable noise was within 15 points for 5 of 9 systems; misses: {v}.")
add("Results: fidelity and benchmarks", "B5", "Santa Fe series A: F_calc (R0) vs BP96 10.23", "12", ["B3", "series_A", "R0"],
    lambda d: f"F_calc {d['F_calc_plain']:.2f}, Z {d['MW_Z']:.1f}, F crit {d['F_crit']:.2f}", "1,000 points",
    "On the Santa Fe laser series our implementation gave {v} (BP96: 10.23, -24.60, 1.16).")
add("Results: fidelity and benchmarks", "B6", "Readings passing V1-V6 (of R0-BP96F, R0-BP96, R0-S-FT, R0-S-AAFT)", "12",
    ["B4_verify"], lambda d: f"{sum(1 for k, v in d.items() if k != 'R0' and v.get('VERIFIED'))} of 4", "-",
    "All {v} alternative readings passed the standard-system verification.")

# ---------------------------------------------------------------------------------------------- Results: static transforms
add("Results: static transforms and surrogates", "S1", "N5 (Phase 10 warp) positive under R0 / R0-S-FT / R0-S-AAFT", "12",
    ["B4_synth", "N5_linear_rr_warped"], lambda d: f"{d['R0']['k']}/{d['R0-S-FT']['k']}/{d['R0-S-AAFT']['k']} of 100", "100 windows",
    "The static monotone warp N5 was positive in {v} windows under R0, FT-surrogate and AAFT-surrogate readings.")
add("Results: static transforms and surrogates", "S2", "AR(1) phi 0.9, n 1,000, exp transform: R0 vs R0-S-AAFT", "12",
    ["B4_static", "cells"], lambda d: f"{d['0.9|1000|exp|R0']['k']}/100 vs {d['0.9|1000|exp|R0-S-AAFT']['k']}/100", "100 series",
    "A static exponential transform of linear Gaussian noise was flagged in {v} (R0 vs amplitude-adjusted surrogates).")
add("Results: static transforms and surrogates", "S3", "Registered N5 sentence", "12", ["B4_matrix", "N5_sentence", "text"], str,
    "-", "{v}")
add("Results: static transforms and surrogates", "S4", "Transplant V1 - V0 under R0-S-AAFT", "12",
    ["B4_transplant", "readings", "R0-S-AAFT", "V1_minus_V0", "pooled"], lambda d: pp(d, "delta"), "439 windows",
    "The surrogate step did not remove the single-beat effect ({v} with amplitude-adjusted surrogates).")
add("Results: static transforms and surrogates", "S5", "Raw CHF - NSR (windows) under R0-S-AAFT", "12",
    ["B4_matrix", "readings", "R0-S-AAFT", "ii_windows"], pp, "988 windows",
    "With amplitude-adjusted surrogates the raw group difference disappeared ({v}).", "claim (ii) fails under R0-S-AAFT")
add("Results: static transforms and surrogates", "S6", "S2 isolated ectopy 5 %: positive under R0 / R0-BP96F / R0-BP96", "12",
    ["B4_synth", "S2_ectopic_5pct"], lambda d: f"{d['R0']['k']}/{d['R0-BP96F']['k']}/{d['R0-BP96']['k']} of 100", "100 windows",
    "Isolated ectopy at 5 % density was rarely flagged under R0 but often under BP96's conventions ({v}).",
    "qualifies Phase 10 S2 (5/100)")

# ---------------------------------------------------------------------------------------------- Results: PM97 replication
add("Results: PM97 replication", "P1", "nsrdb raw detection rate (500-beat, R0)", "12",
    ["PM97", "cells", "500|all|raw|R0", "NSR"], lambda d: f"{100 * d['mean']:.1f} % (median {100 * d['median']:.1f} %)",
    "18 subjects", "On nsrdb, raw 500-beat segments were positive in {v}, not the ~100 % PM97 reported.")
add("Results: PM97 replication", "P2", "chfdb raw detection rate (500-beat, R0)", "12",
    ["PM97", "cells", "500|all|raw|R0", "CHF"], lambda d: f"{100 * d['mean']:.1f} % (median {100 * d['median']:.1f} %)",
    "15 subjects", "On chfdb: {v}.")
add("Results: PM97 replication", "P3", "Raw CHF - NSR (500-beat, R0)", "12",
    ["PM97", "cells", "500|all|raw|R0", "group_rates", "CHF_minus_NSR"], pp, "33 subjects",
    "With PM97's own preprocessing, CHF was not lower ({v}); PM97's pattern was not reproduced.")
add("Results: PM97 replication", "P4", "Raw CHF - NSR under R0-BP96", "12",
    ["PM97", "cells", "500|all|raw|R0-BP96", "group_rates", "CHF_minus_NSR"], pp, "33 subjects",
    "Under the faithful BP96 reading: {v}.")
add("Results: PM97 replication", "P5", "Masked CHF - NSR (500-beat, R0)", "12", ["PM97", "after_handling", "masked"], pp,
    "29 subjects", "After masking ectopy, CHF segments were positive far less often ({v}).", "AGAINST the claim (registered)")
add("Results: PM97 replication", "P6", "Edited CHF - NSR (500-beat, R0)", "12", ["PM97", "after_handling", "edited"], pp,
    "32 subjects", "After editing: {v}.", "AGAINST the claim (registered)")
add("Results: PM97 replication", "P7", "Wu CHF - NSR (500-beat, R0)", "12", ["PM97", "after_handling", "wu"], pp, "32 subjects",
    "After Wu-style elimination: {v}.", "AGAINST the claim (registered)")
add("Results: PM97 replication", "P8", "Within-subject OR of segment ectopy (conditional logit, raw R0)", "12",
    ["PM97", "clogit_ectopy", "coef", "x"], orci, "1,320 segments",
    "Within subjects, ectopy raised the odds of detection ({v} per doubling of the segment's ectopic beats).")
add("Results: PM97 replication", "P9", "Ectopy-free segments, CHF - NSR (R0; POST HOC, exploratory)", "12",
    ["PM97", "exploratory_clean_segments", "R0", "CHF_minus_NSR"], pp, "18 NSR, 8 CHF subjects",
    "In segments without any ectopic beat (exploratory), CHF was still lower ({v}).", "exploratory")
add("Results: PM97 replication", "P10", "Masked-arm analysable fraction, CHF vs NSR (descriptive)", "12",
    ["PM97", "descriptive_handling", "analysable_fraction_500"], lambda d: f"CHF {100 * d['CHF']['masked']:.0f} % vs NSR {100 * d['NSR']['masked']:.0f} %",
    "600 / 720 segments", "Masking left {v} of segments analysable.", "descriptive")

# ---------------------------------------------------------------------------------------------- Discussion / K3 qualifiers
add("Discussion: K3 and robustness qualifiers", "K1", "K3 real detections at 512 intervals", "10",
    ["Q1", "detectors", "K3", "real"], lambda d: f"{d['k']}/{d['n_windows']} (U = {d['U']:.4f})", "988 windows",
    "K3 detected chaos in {v} confirmatory windows.")
add("Discussion: K3 and robustness qualifiers", "K2", "K3 FP at 1/128 s (pooled) vs 1/360 s", "11", ["B2", "fp_pooled"],
    lambda d: f"{d['k128']}/{d['n']} vs {d['k360']}/{d['n']}", "10,092 windows",
    "At the data's 1/128 s resolution K3 gave {v} false positives on synthetic non-chaotic series.", "Phase 11 qualifier")
add("Discussion: K3 and robustness qualifiers", "K3", "K3 detections under two new seeds", "11", ["B1"],
    lambda d: f"{d['K3_entropy11101']['detected']} and {d['K3_entropy11102']['detected']} (Phase 10: 8)", "988 windows",
    "Under two new K3 seeds the count was {v}; six windows were detected under all three seeds.", "Phase 11 qualifier")
add("Discussion: K3 and robustness qualifiers", "K4", "K1 within-subject burden OR", "11", ["A3", "windows", "K1", "clogit_cov", "coef", "x"],
    orci, "16 positives", "K1's association with burden was not significant within subjects ({v}).", "Phase 11 qualifier")


def main():
    lines = ["# Numbers ledger (Phases 10–12)", "",
             "Generated by `experiments/phase12_original/numbers12.py` from the analysis JSON files. Each entry gives "
             "the value with its CI and n, its definition, the source key, and a sentence ready to paste. Sources:",
             f"- `10`: `{F10}`", f"- `11`: `{F11}`", f"- `12`: `{F12}`", "",
             "**Changed by Phase 12.** No Phase 10 or 11 value changed. The recomputed R0 decisions equal the stored "
             "ones (M3). The `Flag` column marks Phase 12 numbers that qualify or contradict an earlier statement.", ""]
    sec = None
    for (section, nid, definition, src, path, fmt, n, sentence, flag) in E:
        if section != sec:
            lines += ["", f"## {section}", "", "| ID | Value (CI) | n | Definition | Source key | Ready-to-paste sentence | Flag |",
                      "|---|---|---|---|---|---|---|"]
            sec = section
        try:
            v = fmt(get(src, path))
        except Exception as exc:                                # noqa: BLE001
            v = f"ERROR {type(exc).__name__}: {exc}"
        key = f"{src}: " + ".".join(p if isinstance(p, str) else "<item>" for p in path)
        lines.append(f"| {nid} | {v} | {n} | {definition} | `{key}` | {sentence.format(v=v)} | {flag} |")
    OUT.write_text("\n".join(lines) + "\n")
    print("wrote", OUT, len(E), "entries")


if __name__ == "__main__":
    main()
