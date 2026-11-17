"""
Revision M1/M2 figures (EXPLORATORY): vector PDF + 300-dpi PNG, Phase 11 style (palette, markers, AIP widths).
Built from summary.json and the CSV tables only.

    PYTHONPATH=. python -m analyses.revision_m1_m2.figures
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from analyses.revision_m1_m2 import analyse as AN
from analyses.revision_m1_m2 import configs as CF
from analyses.revision_m1_m2 import m2gen as M
from experiments.phase11_robustness import figures11 as F11

plt = F11.plt
S, MK, INK, INK2, GRID, W1, W2 = F11.S, F11.MK, F11.INK, F11.INK2, F11.GRID, F11.W1, F11.W2
OUT = AN.RES / "figures"
RCOL = {"R0": S[0], "R0-BP96": S[1], "R0-ORIG": S[2], "R0-BP96F": S[3], "R0-ORIG-minus-S3b": S[3]}
RMK = {"R0": "o", "R0-BP96": "s", "R0-ORIG": "^", "R0-BP96F": "D", "R0-ORIG-minus-S3b": "D"}
TYPE_LABEL = {"ventricular": "ventricular", "atrial": "atrial-type", "spike": "single spike", "missed": "missed beat",
              "extra": "extra beat"}


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def load():
    return json.load(open(AN.SUMMARY))


# ------------------------------------------------------------------------------------------- M1 ablation
def fig_ablation(A):
    ab = A["ablation"]["configs"]
    order = [c for c in AN.CF_ORDER() if c in ab and "closure" in ab[c] and c != "R0"]
    lab = {c: (c if c not in CF.SWITCHES else f"{c}: {CF.SWITCH_LABEL[c]}") for c in order}
    panels = [("T1", "Table 1 (MAD)"), ("ECO", "Ecological loss point"), ("LOG", "Logistic NL(3.7)"),
              ("PM", "nsrdb raw rate"), ("composite", "Mean of four")]
    fig, axes = plt.subplots(1, 5, figsize=(W2, 0.22 * len(order) + 1.0), sharey=True,
                             gridspec_kw={"wspace": 0.12})
    y = np.arange(len(order))[::-1]
    for ax, (b, title) in zip(axes, panels):
        for yi, c in zip(y, order):
            cl = ab[c]["closure"].get(b)
            if not cl or cl.get("closure") is None:
                continue
            col = S[2] if c.startswith("R0-ORIG") else (S[1] if c in CF.COMPOSITES else S[0])
            ax.errorbar(cl["closure"], yi, xerr=[[cl["closure"] - cl["ci"][0]], [cl["ci"][1] - cl["closure"]]],
                        fmt="o", ms=3, color=col, ecolor=col, elinewidth=0.8, capsize=0)
        ax.axvline(0, color=INK2, lw=0.6)
        ax.axvline(1, color=INK2, lw=0.6, ls=":")
        ax.set_title(title, loc="left")
        ax.set_xlim(-1.6, 1.6)
        F11._grid(ax, "x")
        ax.set_xlabel("gap closed")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([lab[c] for c in order])
    fig.suptitle("How far each configuration moves R0 toward the published benchmark values (1 = gap closed)",
                 x=0.0, ha="left", fontsize=7.5, y=1.0)
    save(fig, "M1_ablation_closure")


def fig_logistic(A):
    ab = A["ablation"]["configs"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(W1, 3.4), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    for i, rd in enumerate(("R0", "R0-BP96", "R0-ORIG")):
        sc = ab.get(rd, {}).get("SCAN")
        if not sc:
            continue
        a1.plot(sc["r"], sc["NL"], color=RCOL[rd], lw=0.9, marker=RMK[rd], ms=1.8, label=rd)
        if rd == "R0":
            a2.plot(sc["r"], sc["LE"], color=INK2, lw=0.8)
    a1.plot([3.7, 3.575, 3.565], [0.75, 0.09, 0.0], "x", color=INK, ms=5, label="PB01 Fig. 2 (read)")
    for i, rd in enumerate(("R0", "R0-BP96", "R0-ORIG")):
        lg = ab.get(rd, {}).get("LOG")
        if lg:
            a1.errorbar(3.705 + 0.006 * i, lg["mean_NL"], yerr=[[lg["mean_NL"] - lg["ci"][0]], [lg["ci"][1] - lg["mean_NL"]]],
                        fmt=RMK[rd], ms=3, color=RCOL[rd], mfc="white", elinewidth=0.8)
    a1.set_ylabel("noise limit NL")
    a1.set_title("Logistic map: NL vs r (one series per r; open markers: mean of 100 x0 at r = 3.7)", loc="left")
    a1.legend(frameon=False, loc="upper left")
    F11._grid(a1)
    a2.axhline(0, color=INK2, lw=0.5)
    a2.set_ylabel("Lyapunov exp.")
    a2.set_xlabel("r")
    F11._grid(a2)
    save(fig, "M1_logistic_NL_vs_r")


def fig_dose(A):
    H = A["headline"]["readings"]
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    for i, rd in enumerate(("R0", "R0-BP96", "R0-ORIG", "R0-ORIG-minus-S3b")):
        b = H.get(rd, {}).get("dose_burden")
        if not b:
            continue
        x = np.arange(len(b)) + (i - 1.5) * 0.15
        r = np.array([v["rate"] for v in b]) * 100
        lo = np.array([v["ci"][0] for v in b]) * 100
        hi = np.array([v["ci"][1] for v in b]) * 100
        ge = H[rd]["gee"]["coef"]["x"]
        lab = "S1+S2+S5a (Amendment 3)" if rd == "R0-ORIG-minus-S3b" else rd
        ax.errorbar(x, r, yerr=[r - lo, hi - r], fmt=RMK[rd], ms=3.5, color=RCOL[rd], elinewidth=0.8,
                    label=f"{lab} (GEE OR {ge['OR']:.2f})")
        ax.plot(x, r, color=RCOL[rd], lw=0.7)
    ax.set_xticks(range(5))
    ax.set_xticklabels([f"{v['burden']}\n({v['n']})" for v in H["R0"]["dose_burden"]])
    ax.set_xlabel("annotated ectopic beats in the window (windows)")
    ax.set_ylabel("titration-positive windows (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Positive rate by ectopic burden, 988 windows, by reading", loc="left")
    ax.legend(frameon=False, loc="lower right")
    F11._grid(ax)
    save(fig, "M1_dose_response_by_reading")


# ------------------------------------------------------------------------------------------- M2
def fig_k_types(A):
    inf = pd.read_csv(AN.RES / "m2_influence.csv")
    fig, axes = plt.subplots(2, 5, figsize=(W2, 3.4), sharex=True, sharey=True, gridspec_kw={"hspace": 0.35, "wspace": 0.1})
    xs = {k: i for i, k in enumerate(M.KS)}
    for row, rd in enumerate(("R0", "R0-ORIG")):
        for col, ty in enumerate(M.TYPES):
            ax = axes[row, col]
            d = inf[(inf.reading == rd) & (inf.type == ty) & (inf.placement == "random")].sort_values("k")
            x = [xs[k] for k in d.k]
            ax.fill_between(x, 100 * d.pos_ci_lo, 100 * d.pos_ci_hi, color=RCOL[rd], alpha=0.15, lw=0)
            ax.plot(x, 100 * d.pos_rate, color=RCOL[rd], marker="o", ms=2.5, label="raw")
            dd = d[d.rm_aff_n_ok >= 20]
            ax.plot([xs[k] for k in dd.k], 100 * dd.rm_aff_rate, color=INK, ls="--", marker="v", ms=2.5,
                    label="affected rows removed")
            ds = d[d.sham_n_ok >= 20]
            ax.plot([xs[k] for k in ds.k], 100 * ds.sham_rate, color=INK2, ls=":", marker="x", ms=2.5,
                    label="sham rows removed")
            p = inf[(inf.reading == rd) & (inf.type == ty) & (inf.placement == "periodic")].sort_values("k")
            ax.plot([xs[k] for k in p.k], 100 * p.pos_rate, ls="none", marker="s", mfc="white", ms=3, color=RCOL[rd],
                    label="periodic placement")
            if rd == "R0":
                bp = inf[(inf.reading == "R0-BP96") & (inf.type == ty) & (inf.placement == "random")].sort_values("k")
                ax.plot([xs[k] for k in bp.k], 100 * bp.pos_rate, color=S[1], lw=0.7, alpha=0.8, label="R0-BP96")
            ax.set_title(f"{rd}: {TYPE_LABEL[ty]}", loc="left", fontsize=6.5)
            F11._grid(ax)
            ax.set_ylim(-3, 103)
    for ax in axes[1]:
        ax.set_xticks(range(len(M.KS)))
        ax.set_xticklabels([str(k) for k in M.KS])
        ax.set_xlabel("inserted beats k")
    for ax in axes[:, 0]:
        ax.set_ylabel("positive windows (%)")
    h, l = axes[0, 0].get_legend_handles_labels()
    h2, l2 = axes[1, 0].get_legend_handles_labels()
    fig.legend(h + [h2[0]], ["R0, raw", "affected rows removed", "sham rows removed", "periodic placement", "R0-BP96, raw",
                             "R0-ORIG, raw"], loc="lower center", ncol=6, frameon=False, fontsize=6,
               bbox_to_anchor=(0.5, -0.11))
    save(fig, "M2_positive_rate_vs_k_by_type")


def fig_heatmap(A):
    tm = pd.read_csv(AN.RES / "m2_term_selection.csv")
    fig, axes = plt.subplots(1, 2, figsize=(W2, 3.6), gridspec_kw={"wspace": 0.55})
    for ax, rd in zip(axes, ("R0", "R0-ORIG")):
        d = tm[(tm.reading == rd) & (tm.type == "ventricular") & (tm.placement == "random")]
        if d.empty:
            continue
        top = d.groupby("term").n_selected.sum().sort_values(ascending=False).index[:24]
        ks = [k for k in M.KS if k in set(d.k)]
        Z = np.full((len(top), len(ks)), np.nan)
        for i, t in enumerate(top):
            for j, k in enumerate(ks):
                v = d[(d.term == t) & (d.k == k)]
                Z[i, j] = float(v.freq.iloc[0]) if len(v) else 0.0
        im = ax.imshow(Z, aspect="auto", cmap="Blues", vmin=0, vmax=1)
        ax.set_xticks(range(len(ks)))
        ax.set_xticklabels([str(k) for k in ks])
        ax.set_yticks(range(len(top)))
        labs = []
        for t in top:
            kind = M.term_kind(t)
            labs.append(t.replace("*", "·") + (" *" if kind in ("power", "adjacent_product") else ""))
        ax.set_yticklabels(labs, fontsize=5.5)
        ax.set_xlabel("inserted ventricular beats k (random placement)")
        ax.set_title(f"{rd}: share of positive windows whose best\nnonlinear model contains the term", loc="left")
    fig.colorbar(im, ax=axes, shrink=0.6, label="frequency")
    fig.text(0.0, -0.02, "* power of one lag or product of adjacent lags (can see a premature interval and its pause "
                         "together); yj = interval at lag j", fontsize=5.5)
    save(fig, "M2_term_selection_heatmap")


def fig_ssr(A):
    ss = pd.read_csv(AN.RES / "m2_ssr_share.csv")
    fig, axes = plt.subplots(1, 2, figsize=(W2, 2.3), sharey=True)
    xs = {k: i for i, k in enumerate(M.KS[1:])}
    for ax, rd in zip(axes, ("R0", "R0-ORIG")):
        for i, ty in enumerate(M.TYPES):
            d = ss[(ss.reading == rd) & (ss.type == ty) & (ss.placement == "random") & (ss.subset == "positive")]
            d = d[d.n >= 10].sort_values("k")
            x = np.array([xs[k] for k in d.k]) + (i - 2) * 0.07
            ax.errorbar(x, 100 * d.share_median, yerr=[100 * (d.share_median - d.share_median_ci_lo),
                                                       100 * (d.share_median_ci_hi - d.share_median)],
                        fmt="-o", ms=2.5, lw=0.8, color=(S + ["#7d4cdb"])[i], elinewidth=0.6, label=TYPE_LABEL[ty])
        d = ss[(ss.reading == rd) & (ss.type == "ventricular") & (ss.placement == "random") & (ss.subset == "positive")]
        d = d.sort_values("k")
        ax.plot([xs[k] for k in d.k], 100 * d.frac_aff_rows, color=INK2, ls=":", label="share of rows affected (ventr.)")
        ax.set_xticks(range(len(xs)))
        ax.set_xticklabels([str(k) for k in M.KS[1:]])
        ax.set_xlabel("inserted beats k (random)")
        ax.set_title(f"{rd}: share of RSS_lin - RSS_nl from affected rows", loc="left")
        ax.set_ylim(0, 105)
        F11._grid(ax)
    axes[0].set_ylabel("% of the improvement\n(median, positive windows)")
    axes[0].legend(frameon=False, fontsize=5.5)
    save(fig, "M2_ssr_share_vs_k")


def fig_downslope(A):
    dn = pd.read_csv(AN.RES / "m2_downslope.csv")
    curves = json.load(open(AN.RES / "m2_curves.json"))
    fig, axes = plt.subplots(2, 3, figsize=(W2, 4.0), gridspec_kw={"hspace": 0.75, "wspace": 0.75})
    xs = {k: i for i, k in enumerate(M.KS)}
    d = dn[(dn.reading == "R0") & (dn.type == "ventricular")]
    for pl, mk, ls in (("random", "o", "-"), ("periodic", "s", "none")):
        g = d[d.placement == pl].sort_values("k")
        x = [xs[k] for k in g.k]
        ax = axes[0, 0]
        ax.errorbar(x, g["lag1"], yerr=[g["lag1"] - g["lag1_ci_lo"], g["lag1_ci_hi"] - g["lag1"]], fmt=mk, ls=ls, ms=2.8,
                    lw=0.8, color=S[0], mfc="white" if pl == "periodic" else S[0], elinewidth=0.6, label=pl)
        ax.set_title("lag-1 autocorrelation", loc="left")
        if pl == "random":
            gg = g[g.k > 0]
            axes[0, 1].plot([xs[k] for k in gg.k], gg.relimp_median, marker="o", ms=2.8, color=S[0], label="random (median)")
        else:
            axes[0, 1].text(0.02, 0.22, "periodic k = 16, 26, 51:\n" + ", ".join(f"{v:.1f}" for v in g.relimp_median) +
                            " (off scale)", transform=axes[0, 1].transAxes, fontsize=5.5, color=INK2)
        axes[0, 1].set_title("1 - RSS_nl/RSS_lin (median)", loc="left")
        axes[0, 1].axhline(0, color=INK2, lw=0.5)
        axes[0, 2].plot(x, g.pos_rate * 100, marker=mk, ls=ls, color=S[0], mfc="white" if pl == "periodic" else S[0],
                        ms=2.8, label=pl)
    g = d[d.placement == "random"].sort_values("k")
    x = [xs[k] for k in g.k]
    ax = axes[1, 0]
    for c, col, lab in (("nrss_lin_aff", S[1], "linear, affected rows"), ("nrss_nl_aff", S[2], "nonlinear, affected rows"),
                        ("nrss_lin_oth", S[1], "linear, other rows"), ("nrss_nl_oth", S[2], "nonlinear, other rows")):
        gg = g[g[f"{c}_median"].notna()]
        ax.plot([xs[k] for k in gg.k], gg[f"{c}_median"], color=col, ls="-" if "aff" in c else ":",
                marker="o" if "aff" in c else None, ms=2.5, label=lab)
    ax.set_ylim(0, 1.4)
    ax.set_title("normalised RSS (median)", loc="left")
    ax.legend(frameon=False, fontsize=5, loc="upper left")
    ax = axes[1, 1]
    ax.plot(x, g.kappa_lin_median, color=S[1], marker="o", ms=2.5, label="best linear order (median)")
    ax2 = ax.twinx()
    ax2.plot(x, g.max_abs_z_aff, color=S[3], marker="^", ms=2.5, label="max |z| of affected intervals")
    ax.set_title("linear order; beat size", loc="left")
    ax.set_ylabel("best linear order (orange)", fontsize=6)
    ax2.set_ylabel("max |z| of beats (yellow)", fontsize=6)
    ax2.spines["right"].set_visible(True)
    ax = axes[1, 2]
    for kk, col in ((1, S[0]), (8, S[2]), (51, S[1])):
        cv = [c for c in curves if c["reading"] == "R0" and c["type"] in ("ventricular", "none") and c["k"] == kk
              and c["placement"] == "random" and c["seed"] == 11000]
        if cv:
            c = cv[0]
            ax.plot(np.arange(1, len(c["lin"]) + 1), c["lin"], color=col, ls=":", lw=0.9)
            ax.plot(np.arange(1, len(c["nl"]) + 1), c["nl"], color=col, lw=0.9, label=f"k = {kk} (kappa_nl {c['kappa_nl']})")
    ax.set_xlabel("number of terms r")
    ax.set_title("C(r) curves, seed 11000", loc="left")
    ax.text(0.98, 0.98, "dotted: linear\nsolid: nonlinear", transform=ax.transAxes, ha="right", va="top", fontsize=5)
    ax.legend(frameon=False, fontsize=5)
    for a in list(axes[0]) + list(axes[1, :2]):
        a.set_xticks(range(len(M.KS)))
        a.set_xticklabels([str(k) for k in M.KS])
        a.set_xlabel("inserted ventricular beats k")
        F11._grid(a)
    axes[0, 2].set_title("positive windows (%)", loc="left")
    axes[0, 0].legend(frameon=False, fontsize=5.5)
    save(fig, "M2_downslope_diagnostics")


def main():
    A = load()
    if "ablation" in A:
        fig_ablation(A)
        fig_logistic(A)
    if "headline" in A:
        fig_dose(A)
    if "m2" in A:
        fig_k_types(A)
        fig_heatmap(A)
        fig_ssr(A)
        fig_downslope(A)
    print("FIGURES DONE")


if __name__ == "__main__":
    main()
