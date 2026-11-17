"""
Revision M1.1 (EXPLORATORY): every way R0 differs, or might differ, from BP96, PM97 and PB01 (PLAN.md Section 2).
Sources: BP96 (Nature 381:215), PM97 (Nature 389:492), PB01 (PNAS 98:7107, PMC34630), FIDELITY.md, Phase 10
METHODS.md, titration10.py.  Quotes are kept to a few words; page numbers refer to the journal pages.  The switch
column links to the ablation (ablation_benchmarks.csv); the composite closure is filled in from summary.json.
"""
from __future__ import annotations

import json

import pandas as pd

from analyses.revision_m1_m2 import analyse as AN

COLS = ("id", "aspect", "source", "original", "R0", "alternative_closest_to_original", "confidence",
        "decision_relevant", "switch")

ROWS = [
    ("D1", "N in C(r) = log eps(r) + r/N", "BP96 p. 216 Eq. 3, Eq. 2 (sums over n = 1..N); BP96 Fig. 2 legend df; "
     "PM97 p. 494 Fig. 2 legend",
     "N is the length of the series (PM97: N = series length; BP96 Fig. 2a df 968 and 944 for 1,000 points)",
     "N = number of regression rows (917 for 1,000 points, 429 for 512)", "N = series length", "high", "yes", "S1"),
    ("D2", "F statistic", "BP96 p. 216 Fig. 1 legend", "F_calc = (eps_lin)^2/(eps_nl)^2, the plain ratio of normalised "
     "error variances ('9,387 >> F(0.01, 915, 915)')", "(RSS_lin/df_lin)/(RSS_nl/df_nl)", "plain ratio RSS_lin/RSS_nl",
     "high", "yes", "S2"),
    ("D3", "F-test degrees of freedom", "BP96 p. 217 Fig. 2a/2b legends", "df = N - r with N = series length "
     "(F(0.01, 968, 944); F(0.06, 1855, 1783))", "df = rows - r", "df = length - r", "high", "yes (via D1)", "S1"),
    ("D4", "Test sidedness", "BP96 p. 216 Fig. 1 legend", "one-sided: nonlinear 'significantly more predictive'; "
     "F_calc >> critical value; Z << -2.33", "one-sided upper tail", "same", "high", "matches", "none"),
    ("D5", "Significance level", "BP96 Figs. 1-2 legends, Table 1 ('99 % cutoff'); PB01 Fig. 1 legend ('1% "
     "significance level')", "1 %", "1 %", "same", "high", "matches", "none"),
    ("D6", "Gaussianity check", "BP96 p. 216", "F-test 'for models with gaussian residuals'; the assumption 'was "
     "verified ... by using a chi2-test with a 99% cutoff'; what happens on rejection is not stated", "no check",
     "chi-square (99 %) on both residual sets; on rejection the Mann-Whitney test", "medium (binning and fallback "
     "unspecified)", "yes", "S3"),
    ("D7", "Mann-Whitney test", "BP96 p. 216; Fig. 1 (Z = -37.13), Fig. 2a (Z = -24.60)", "'Alternatively, the results "
     "are confirmed by' the Mann-Whitney rank-sum statistic; the ranked quantity is not stated", "not used",
     "one-sided Mann-Whitney on squared residuals, always", "low-medium", "yes", "S3b"),
    ("D8", "Surrogate step", "BP96 p. 216", "best models also fitted to surrogates 'with the same autocorrelation (and "
     "power spectrum)'; the original's nonlinear error must beat all surrogate errors; surrogates 'play only a "
     "confirmatory role'; number and type not stated", "none",
     "10 phase-randomised (FT) surrogates, worst case decides (AAFT stored from Phase 12; IAAFT postdates BP96 and is "
     "not original)", "high that a step existed; medium on its form", "yes (can only remove positives)", "S4"),
    ("D9", "Linear search range", "BP96 p. 216; Figs. 1b, 2a", "best linear model 'by searching for kappa_lin which "
     "minimizes C(r) with d = 1'; linear curves span the same term range as the nonlinear (about 85 and 56 terms)",
     "kappa_lin <= 83 (84 terms)", "same (kept matched in S5a/S5b)", "medium-high", "matches", "none"),
    ("D10", "Matched maximum term counts", "BP96 Figs. 1b, 2a; PM97 p. 494 Fig. 2 legend", "linear and nonlinear "
     "curves over the same r range; PM97 compares with the best linear model 'with a similar number of' terms",
     "84 linear vs 84 nonlinear", "same", "medium-high", "matches", "none"),
    ("D11", "Polynomial degree", "BP96 p. 216 ('increasing kappa and d > 1'); Fig. 2a (kappa = 3, d up to 5); Fig. 1c "
     "(kappa = 2, d = 8 model 'obtained with the algorithm')", "no fixed degree limit; d up to 5 and 8 in the "
     "examples", "d <= 3 (Wysocki 2006 routine default)", "d up to 8 with M <= 84 (kappa 1-2: d 8; 3: 6; 4: 4; 5-6: 3)",
     "medium", "possibly", "S5a"),
    ("D12", "Embedding dimension kappa", "BP96 Figs. 1b, 1c, 2a; PM97 p. 494 Fig. 3", "no fixed limit (kappa = 6, 3, 2 "
     "in the examples); PM97 reports 'the first 80 nonlinear coefficient estimates', more than R0's largest model "
     "has (77)", "kappa <= 6", "kappa <= 7, d <= 3 (M = 120), linear memory 119", "low-medium", "possibly", "S5b"),
    ("D13", "Definition of eps", "BP96 p. 215 Eq. 2; p. 216", "eps^2 = normalised residual variance; eps called the "
     "'error standard deviations'; C uses log eps", "eps = sqrt(RSS/SST)", "same (eps as variance: Phase 11 R6, "
     "stored)", "high", "matches", "none"),
    ("D14", "Centring of SST", "BP96 Eq. 2", "y-bar = mean over the series", "mean over the target rows",
     "same", "high", "no by construction (constant shift of C; cancels in F)", "none"),
    ("D15", "Standardisation / preprocessing", "BP96 p. 217", "no preprocessing", "series standardised",
     "none", "high", "no by construction (affine invariance)", "none"),
    ("D16", "First kappa points (target rows)", "BP96 p. 215 Eq. 2; coefficient estimation from autocorrelations",
     "sums over n = 1..N; pre-sample lags not described", "common rows n >= 83 for every model",
     "every point a target, pre-sample lags 0 (zero padding); N = length", "low", "possibly", "S6"),
    ("D17", "Decision rule", "BP96 p. 216; Fig. 1d; Table 1 legend; PB01 Fig. 1 legend", "'nonlinear determinism is "
     "indicated if d_opt > 1', corroborated by the test; Fig. 1d and Table 1 measure loss of detection with 'a 99% "
     "cutoff'", "C_nl < C_lin and F-test at 1 %", "test at 1 % alone (as Fig. 1d / Table 1)", "medium", "yes", "S7"),
    ("D18", "Coefficient estimation", "BP96 p. 215", "recursive Gram-Schmidt from linear and nonlinear "
     "autocorrelations", "modified Gram-Schmidt least squares on the target rows", "same in substance (edges: D16)",
     "medium", "no (edges via D16)", "none"),
    ("D19", "Titration noise", "PB01 p. 7107, Fig. 1 legend", "'white (or linearly correlated) noise'; 'a white noise "
     "series ... of equal length and unit variance'", "Gaussian white, standardised sample", "same", "high",
     "NL only", "none"),
    ("D20", "NL definition", "PB01 Fig. 1 legend", "'NL = alpha/sigma_y'", "alpha/sigma_y", "same", "high", "NL only",
     "none"),
    ("D21", "Titration search and step", "PB01 Fig. 1 legend; Wu 2009", "'a bisection scheme'; step not stated "
     "(Wu 2009: increments)", "bisection, 10 halvings (0.2 % of sigma_y)", "same (increments: Phase 11 R7, stored)",
     "high", "NL only", "none"),
    ("D22", "Noise realisations", "PB01 Fig. 1 legend", "'repeated for 10 different examples of white noise'", "10",
     "same", "high", "NL only", "none"),
    ("D23", "Series length of benchmarks", "BP96 Table 1; PB01 Fig. 1", "1,000 points", "1,000 (benchmarks)",
     "same", "high", "matches", "none"),
    ("D24", "Initial conditions and transients of simulated systems", "BP96 (not stated)", "not stated",
     "random initial condition per seed, 1,000 iterations discarded (Phase 12)", "not identifiable; 100 seeds "
     "spread the initial conditions", "low", "benchmark only", "none"),
    ("D25", "Sampling of continuous flows", "BP96 p. 216", "'the optimal time delay tau_opt may be chosen so as to "
     "maximize the difference between' eps_lin and eps_nl", "fixed sampling (Phase 12)", "sampling step chosen to "
     "maximise eps_lin - eps_nl (P1)", "medium", "benchmark only", "P1"),
    ("D26", "Table 1 system parameters", "BP96 Table 1", "not given", "standard parameters (Phase 12)",
     "not identifiable", "low", "benchmark only", "none"),
    ("D27", "Correlated noise in benchmarks", "BP96 p. 217; Fig. 1d legend", "'coloured measurement noises ... with the "
     "same autocorrelation as the original series'; % = ratio of SDs", "FT surrogate of the clean series, scaled",
     "same", "high", "matches", "none"),
    ("D28", "Minimum length guard", "not stated", "-", ">= 3 M rows (>= 300 when masked)", "-", "-", "rare", "none"),
    ("D29", "Titration of surrogate readings", "PB01 (no surrogates in titration)", "-", "surrogate readings take NL "
     "from R0 when the requirement holds noise-free", "-", "-", "NL only", "none"),
    ("D30", "PM97 subject selection", "PM97 p. 493", "subjects selected for a stable mean heart rate and few ectopic "
     "and undetected beats", "all subjects (primary); rank-based PM97-like subset (Phase 12)", "-", "medium",
     "protocol", "none"),
    ("D31", "PM97 kappa and d", "PM97", "not stated", "R0 defaults", "see D11, D12", "low", "possibly", "S5a/S5b"),
    ("D32", "Constant term counted in r", "BP96 Eq. 1; Figs. 1b, 2a", "a_0 is a term; C(1) = 0", "same", "same",
     "high", "matches", "none"),
]


def build():
    df = pd.DataFrame(ROWS, columns=COLS)
    S = json.load(open(AN.SUMMARY)) if AN.SUMMARY.exists() else {}
    cfgs = S.get("ablation", {}).get("configs", {})
    comp = []
    for sw in df.switch:
        parts = [s for s in sw.split("/") if s in cfgs]
        if not parts:
            comp.append(None if sw != "P1" else S.get("ablation", {}).get("P1_on_R0", {}).get("composite_closure"))
            continue
        comp.append("; ".join(f"{s}: {cfgs[s]['closure']['composite']['closure']:+.2f}" for s in parts))
    df["composite_closure_of_switch"] = comp
    df.to_csv(AN.RES / "deviations.csv", index=False)
    lines = ["# Deviations of R0 from BP96, PM97 and PB01 (EXPLORATORY, post hoc)", "",
             "One row per way R0 differs, or might differ, from the original papers. 'Switch' links to the ablation "
             "(`ablation_benchmarks.csv`); the last column is that switch's composite closure of the benchmark gap "
             "(1 = gap closed, 0 = no change, negative = moves away).", "",
             "| id | aspect | original (source) | R0 | closest alternative | confidence | decision-relevant | switch | "
             "closure |", "|---|---|---|---|---|---|---|---|---|"]
    for r in df.itertuples():
        lines.append(f"| {r.id} | {r.aspect} | {r.original} ({r.source}) | {r.R0} | {r.alternative_closest_to_original} | "
                     f"{r.confidence} | {r.decision_relevant} | {r.switch} | {r.composite_closure_of_switch or ''} |")
    (AN.RES / "deviations.md").write_text("\n".join(lines) + "\n")
    return df


if __name__ == "__main__":
    build()
