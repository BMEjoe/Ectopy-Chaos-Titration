"""
Phase 11 C2: publication figures F1-F7 (vector PDF + 300-dpi PNG; AIP widths 3.37 in single, 6.69 in double column;
colour-blind-safe categorical palette validated with the dataviz validator: blue, orange, aqua, yellow; every series
also has its own marker / line style).

    python -m experiments.phase11_robustness.figures11

Inputs: Phase 10 analysis.json, real.jsonl, seg.jsonl, verification JSON; Phase 11 analysis11.json.
Output: results/conf/analysis/figures/F*.pdf / .png
"""
from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
P10 = HERE.parent / "phase10_final" / "results"
OUT = HERE / "results" / "conf" / "analysis" / "figures"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # categorical slots 1-4 (validated)
MK = ["o", "s", "^", "D"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
W1, W2 = 3.37, 6.69
plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5, "legend.fontsize": 6.5,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.2})


def _grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0, c - h), min(1, c + h)


# ------------------------------------------------------------------------------------------------ F1
def f1(A):
    Q = A["Q2e"]
    keys = sorted([k for k in Q if k.endswith("|512")], key=lambda s: (s.startswith("vdp"), s))
    arms = [("raw", "raw"), ("raw_q128", "raw, 1/128 s"), ("edited", "edited"), ("masked", "masked")]
    M = np.full((len(keys), len(arms)), np.nan)
    for i, k in enumerate(keys):
        v = Q[k]
        for j, (a, _) in enumerate(arms):
            den = v["n_windows"] if a in ("raw", "raw_q128") else v[f"{a}_eligible" if a == "edited" else "masked_analysable"]
            M[i, j] = 100 * v[a] / den if den else np.nan
    fig, ax = plt.subplots(figsize=(W1, 0.15 * len(keys) + 0.8))
    cm = matplotlib.colormaps["Blues"].copy()
    cm.set_bad("#f0efec")
    im = ax.imshow(M, cmap=cm, vmin=0, vmax=100, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            t = "n/a" if np.isnan(M[i, j]) else f"{M[i, j]:.0f}"
            ax.text(j, i, t, ha="center", va="center", fontsize=5,
                    color="white" if (not np.isnan(M[i, j]) and M[i, j] > 60) else INK)
    ax.set_yticks(range(len(keys)), [k.split("|")[0].replace("_", " ") for k in keys], fontsize=5)
    ax.set_xticks(range(len(arms)), [a for _, a in arms])
    ax.xaxis.tick_top()
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("titration-positive windows (%)")
    ax.set_title("False positives on non-chaotic synthetic RR (512 intervals; masked/edited: % of usable)", fontsize=6, pad=18)
    save(fig, "F1_false_positive_heatmap")


# ------------------------------------------------------------------------------------------------ F2
def f2(A):
    R = [json.loads(line) for line in open(P10 / "conf" / "real.jsonl")]
    b = np.array([r["burden"] for r in R])
    y = np.array([r["TIT_raw"]["positive"] for r in R])
    chf = np.array([r["group"] == "CHF" for r in R])
    bins = [(0, 0, "0"), (1, 1, "1"), (2, 4, "2–4"), (5, 15, "5–15"), (16, 10 ** 9, "≥16")]
    g = A["Q3"]["P3_gee"]["TIT"]["coef"]
    xs = np.linspace(0, np.log2(1 + 200), 200)
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    for gi, (grp, mask) in enumerate((("NSR", ~chf), ("CHF", chf))):
        px, py, lo, hi = [], [], [], []
        for a, c, _ in bins:
            m = mask & (b >= a) & (b <= c)
            if m.sum() == 0:
                continue
            k, n = int(y[m].sum()), int(m.sum())
            px.append(np.mean(np.log2(1 + b[m])) + (gi - 0.5) * 0.08)
            py.append(k / n)
            w = wilson(k, n)
            lo.append(k / n - w[0])
            hi.append(w[1] - k / n)
        ax.errorbar(px, py, yerr=[lo, hi], fmt=MK[gi], color=S[gi], ms=4, capsize=1.5, lw=0.8, label=f"{grp} (bins, Wilson CI)")
        eta = g["const"]["beta"] + g["x"]["beta"] * xs + g["chf"]["beta"] * gi
        ax.plot(xs, 1 / (1 + np.exp(-eta)), color=S[gi], ls="-" if gi == 0 else "--", lw=1.2, label=f"{grp} GEE fit")
    ax.set_xticks([np.log2(1 + v) for v in (0, 1, 3, 7, 15, 31, 63, 127)], ["0", "1", "3", "7", "15", "31", "63", "127"])
    ax.set_xlabel("ectopic beats in the 512-interval window (log2(1 + burden) scale)")
    ax.set_ylabel("titration-positive windows")
    ax.set_ylim(0, 1)
    c = g["x"]
    ax.set_title(f"Dose-response: OR {c['OR']:.2f} [{c['OR_ci95'][0]:.2f}, {c['OR_ci95'][1]:.2f}] per doubling", loc="left")
    ax.legend(frameon=False, loc="lower right")
    _grid(ax)
    save(fig, "F2_dose_response")


# ------------------------------------------------------------------------------------------------ F3
def f3(A):
    R = [json.loads(line) for line in open(P10 / "conf" / "seg.jsonl")]
    arms = [("raw", "TIT_raw", "P2_raw_by_group"), ("masked", "TIT_masked", "masked_by_group"),
            ("edited", "TIT_edited", "edited_by_group"), ("Wu", "TIT_wu", "wu_by_group")]
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    rng = np.random.default_rng(11000)
    for ai, (name, key, akey) in enumerate(arms):
        for gi, grp in enumerate(("NSR", "CHF")):
            vals = []
            for r in R:
                if r["group"] != grp:
                    continue
                ss = [s[key] for s in r["segments"] if not s["gap"]]
                if key == "TIT_masked":
                    ss = [x for x in ss if x.get("analysable")]
                elif key in ("TIT_edited", "TIT_wu"):
                    ss = [x for x in ss if x.get("eligible", True) and x.get("analysable", True)]
                if ss:
                    vals.append(np.mean([x["positive"] for x in ss]))
            x0 = ai + (gi - 0.5) * 0.36
            ax.scatter(x0 + rng.uniform(-0.09, 0.09, len(vals)), vals, s=3, color=S[gi], alpha=0.45, lw=0)
            gr = A["Q2"][akey][grp]
            m, (lo, hi) = gr["mean_subject_rate"], gr["ci95"]
            ax.errorbar([x0], [m], yerr=[[m - lo], [hi - m]], fmt=MK[gi], color=INK, mfc=S[gi], ms=5, capsize=2, lw=1,
                        label=grp if ai == 0 else None, zorder=5)
    ax.set_xticks(range(4), [a for a, _, _ in arms])
    ax.set_ylabel("per-subject detection rate (12-min segments)")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="upper right", title="mean [95 % CI]; dots = subjects", title_fontsize=6)
    _grid(ax)
    save(fig, "F3_group_rate_by_arm")


# ------------------------------------------------------------------------------------------------ F4
def f4(A):
    D = A["Q1"]["detectors"]["K3"]
    fams = [("henon", "Hénon"), ("logistic", "logistic"), ("coupled_vdp", "coupled vdP"), ("phase_reset", "phase-reset")]
    fig, axs = plt.subplots(1, 4, figsize=(W2, 1.9), sharey=True)
    for ax, (fam, title) in zip(axs, fams):
        keys = sorted([k for k in D["curves"] if D["curves"][k]["family"] == fam], key=lambda k: D["curves"][k]["lam"])
        for i, k in enumerate(keys):
            e = D["exclusion_L0"][k]
            f = e["design"]["f"]
            ax.plot(f, e["design"]["pi_upper"], marker=MK[i], ms=3, color=S[i], lw=1.1, label=f"{D['curves'][k]['param']}")
            ax.plot(f, e["design_wilson_lower"]["pi_upper"], color=S[i], lw=0.8, ls="--")
        for th in (0.05, 0.20):
            ax.axhline(th, color=INK2, lw=0.6, ls=":")
        ax.set_yscale("log")
        ax.set_ylim(0.008, 1.2)
        ax.set_title(title, loc="left")
        ax.set_xlabel("chaos fraction f")
        ax.legend(frameon=False, fontsize=5, loc="lower left")
        _grid(ax, "both")
    axs[0].set_ylabel("exclusion bound π_upper")
    fig.suptitle("K3: largest prevalence of additive chaos compatible with the real detection rate "
                 "(solid: empirical p; dashed: Wilson-lower p; dotted: 0.05, 0.20)", fontsize=6.5, x=0.01, ha="left", y=1.04)
    save(fig, "F4_K3_exclusion_bound")


# ------------------------------------------------------------------------------------------------ F5
def f5(A, B):
    v = json.load(open(P10 / "verification" / "titration_verification.json"))["verification"]["V4"]["scan"]
    r = np.array([a for a, _, _ in v])
    le = np.array([b_ for _, b_, _ in v])
    nl = np.array([c for _, _, c in v])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 2.0), gridspec_kw={"width_ratios": [1.4, 1]})
    a1.plot(r, nl, "-o", color=S[0], ms=2, lw=1, label="noise limit NL")
    a1b = a1.inset_axes([0.0, -0.62, 1.0, 0.42], sharex=a1)
    a1b.plot(r, le, color=S[1], lw=1)
    a1b.axhline(0, color=INK2, lw=0.6)
    a1b.set_ylabel("Lyapunov exp.")
    a1b.set_xlabel("logistic r")
    a1.set_ylabel("NL (fraction of SD)")
    a1.set_title("Titration verification: logistic bifurcation scan", loc="left")
    plt.setp(a1.get_xticklabels(), visible=False)
    _grid(a1)
    ver = (B.get("A5") or {}).get("verification", {})
    rds = [k for k in ("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7") if k in ver and "V1" in ver[k]]
    checks = ("V1", "V2", "V3", "V4", "V5", "V6")
    for i, rd in enumerate(rds):
        for j, c in enumerate(checks):
            ok = ver[rd][c]
            a2.add_patch(plt.Rectangle((j, i), 0.92, 0.92, color=S[2] if ok else S[1], alpha=0.85 if ok else 0.9))
            a2.text(j + 0.46, i + 0.46, "pass" if ok else "FAIL", ha="center", va="center", fontsize=5, color=INK)
    a2.set_xlim(0, len(checks))
    a2.set_ylim(0, max(len(rds), 1))
    a2.set_xticks(np.arange(len(checks)) + 0.46, checks)
    a2.set_yticks(np.arange(len(rds)) + 0.46, rds)
    a2.invert_yaxis()
    for s in a2.spines.values():
        s.set_visible(False)
    a2.tick_params(length=0)
    a2.set_title("Verification V1–V6 under readings R0–R7", loc="left")
    save(fig, "F5_validation_summary")


# ------------------------------------------------------------------------------------------------ F6
def f6(B):
    A1 = B.get("A1", {})
    cells = [("windows_masked", "windows\nmasked"), ("windows_edited", "windows\nedited"), ("segments_masked", "segments\nmasked"),
             ("segments_edited", "segments\nedited"), ("segments_wu", "segments\nWu")]
    cells = [(k, n) for k, n in cells if A1.get(k, {}).get("pooled")]
    if not cells:
        return
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 2.2), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.62})
    x = np.arange(len(cells))
    ft = [A1[k]["pooled"]["F_true"] for k, _ in cells]
    fs = [A1[k]["pooled"]["F_sham"] for k, _ in cells]
    a1.bar(x - 0.19, ft, 0.36, color=S[0], label="true ectopy mask")
    a1.bar(x + 0.19, fs, 0.36, color=S[3], label="sham mask (same size, shifted)", hatch="///", edgecolor="white", lw=0)
    a1.set_xticks(x, [n for _, n in cells], fontsize=6)
    a1.set_ylabel("fraction of raw positives removed")
    a1.set_ylim(0, 1)
    a1.legend(frameon=False, loc="upper left")
    a1.set_title("True vs sham handling (raw-positive ectopic units)", loc="left")
    _grid(a1)
    for gi, g in enumerate(("NSR", "CHF", "pooled")):
        for i, (k, _) in enumerate(cells):
            v = A1[k].get(g)
            if not v:
                continue
            y = i + (gi - 1) * 0.22
            a2.errorbar([v["D"]], [y], xerr=[[v["D"] - v["D_ci95"][0]], [v["D_ci95"][1] - v["D"]]], fmt=MK[gi],
                        color=[S[0], S[1], INK][gi], ms=3.5, capsize=1.5, lw=0.8, label=g if i == 0 else None)
    a2.axvline(0, color=INK2, lw=0.7)
    a2.set_yticks(x, [n.replace("\n", " ") for _, n in cells])
    a2.invert_yaxis()
    a2.set_xlabel("D = removed (true) − removed (sham), 95 % cluster CI")
    a2.set_xlim(0, None)
    a2.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handletextpad=0.3)
    _grid(a2, "x")
    save(fig, "F6_true_vs_sham")


# ------------------------------------------------------------------------------------------------ F7
def f7(B):
    T = (B.get("A2_transplant") or {}).get("delta_vs_V0", {})
    if not T:
        return
    lab = {"V1": "V1 ventricular beat", "V2": "V2 V1 + HRT (empirical)", "V3": "V3 V1 + HRT (normal)",
           "V4": "V4 V1 + HRT (blunted)", "V5": "V5 atrial beat", "V6": "V6 couplet", "V7": "V7 HRT only (no beat)"}
    vs = [v for v in lab if v in T]
    fig, ax = plt.subplots(figsize=(W1, 2.4))
    for gi, g in enumerate(("NSR", "CHF", "pooled")):
        for i, v in enumerate(vs):
            d = T[v].get(g)
            if not d:
                continue
            y = i + (gi - 1) * 0.22
            ax.errorbar([d["delta"]], [y], xerr=[[d["delta"] - d["ci95"][0]], [d["ci95"][1] - d["delta"]]], fmt=MK[gi],
                        color=[S[0], S[1], INK][gi], ms=3.5, capsize=1.5, lw=0.8, label=g if i == 0 else None)
    ax.axvline(0, color=INK2, lw=0.7)
    ax.set_yticks(range(len(vs)), [lab[v] for v in vs])
    ax.invert_yaxis()
    p = T[vs[0]].get("pooled", {})
    ax.set_xlabel(f"Δ titration-positive rate vs untouched window (V0 = {100 * p.get('rate_a', np.nan):.0f} %)")
    ax.set_title("One ectopic event transplanted into clean real windows", loc="left")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.98), handletextpad=0.3)
    _grid(ax, "x")
    save(fig, "F7_transplant_variants")


def main():
    A = json.load(open(P10 / "conf" / "analysis" / "analysis.json"))
    p = HERE / "results" / "conf" / "analysis" / "analysis11.json"
    B = json.load(open(p)) if p.exists() else {}
    for fn, args in ((f1, (A,)), (f2, (A,)), (f3, (A,)), (f4, (A,)), (f5, (A, B)), (f6, (B,)), (f7, (B,))):
        try:
            fn(*args)
            print(fn.__name__, "ok")
        except Exception as exc:                                # noqa: BLE001
            print(fn.__name__, "ERROR", type(exc).__name__, exc)


if __name__ == "__main__":
    main()
