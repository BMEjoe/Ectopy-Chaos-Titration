"""
Phase 12 Tier 2 figures (from analysis12.json only; vector PDF + 300-dpi PNG; AIP widths; Phase 11 palette and style).
  F12-1  BP96 Fig. 1d analogue: 1/F_calc vs % added correlated noise, ecological model (R0, R0-BP96)
  F12-2  static transforms and the surrogate step; transplant V1 - V0 by reading
  F12-3  PM97 replication: per-subject detection rates by group and arm (R0)

    PYTHONPATH=. python -m experiments.phase12_original.figures12
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from experiments.phase11_robustness import figures11 as F11

plt = F11.plt
S, MK, INK, INK2, W1, W2 = F11.S, F11.MK, F11.INK, F11.INK2, F11.W1, F11.W2
HERE = Path(__file__).resolve().parent
A = json.load(open(HERE / "results" / "conf" / "analysis" / "analysis12.json"))
OUT = HERE / "results" / "conf" / "analysis" / "figures12"


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def f1():
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    for i, rd in enumerate(("R0", "R0-BP96")):
        c = A["B1"]["readings"][rd]["curve"]
        p = np.array(sorted(int(k) for k in c))
        med = np.array([c[str(k)]["inv_F_median"] for k in p])
        lo = np.array([c[str(k)]["inv_F_q25"] for k in p])
        hi = np.array([c[str(k)]["inv_F_q75"] for k in p])
        ax.fill_between(p + 0.6 * i, lo, hi, color=S[i], alpha=0.18, lw=0)
        ax.plot(p + 0.6 * i, med, marker=MK[i], ms=3, color=S[i], label=f"{rd} (median, IQR; 20 seeds)")
    fc = A["B1"]["readings"]["R0"]["F_crit"]["median"]
    ax.axhline(1 / fc, color=INK2, ls="--", lw=0.8)
    ax.text(1, 1 / fc + 0.04, "1 % critical value", fontsize=6, color=INK2)
    ax.axvline(50, color=S[3], lw=0.8, ls=":")
    ax.text(50.8, 0.05, "BP96: lost at ~50 %", fontsize=6, color=INK2)
    ax.set_xlabel("added correlated noise (% of the series SD)")
    ax.set_ylabel(r"$1/F_{\mathrm{calc}}$")
    ax.set_ylim(0, 1.6)
    ax.set_title("BP96 Fig. 1d analogue: ecological model, 1,000 points", loc="left")
    ax.legend(frameon=False, loc="upper left")
    F11._grid(ax)
    save(fig, "F12_1_B1_fig1d")


def f2():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 2.3), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.45})
    cells = A["B4_static"]["cells"]
    trs = ["cube", "exp", "N5", "tanh"]
    rds = ["R0", "R0-BP96", "R0-S-FT", "R0-S-AAFT"]
    x = np.arange(len(trs) + 1)
    w = 0.2
    for i, rd in enumerate(rds):
        vals = [cells[f"0.9|1000|{t}|{rd}"]["k"] for t in trs] + [A["B4_synth"]["N5_linear_rr_warped"][rd]["k"]]
        a1.bar(x + (i - 1.5) * w, vals, w, color=S[i], label=rd)
    a1.set_xticks(x, [f"AR(1)\n{t}" for t in trs] + ["Phase 10\nN5 (RR)"], fontsize=6)
    a1.set_ylabel("positive / 100")
    a1.set_title("Static transforms of linear Gaussian data (φ = 0.9, n = 1,000)", loc="left")
    a1.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.22))
    F11._grid(a1)
    tp = A["B4_transplant"]["readings"]
    order = ["R0", "R0-BP96F", "R0-BP96", "R0-S-FT", "R0-S-AAFT"]
    for i, rd in enumerate(order):
        d = tp[rd]["V1_minus_V0"]["pooled"]
        a2.errorbar([100 * d["delta"]], [i], xerr=[[100 * (d["delta"] - d["ci95"][0])], [100 * (d["ci95"][1] - d["delta"])]],
                    fmt="o", color=INK, ms=3.5, capsize=1.5, lw=0.8)
    a2.axvline(0, color=INK2, lw=0.7)
    a2.set_yticks(range(len(order)), order)
    a2.invert_yaxis()
    a2.set_xlim(0, None)
    a2.set_xlabel("one ventricular-type beat: V1 − V0 (pp, 95 % CI)")
    a2.set_title("Transplant effect by reading", loc="left")
    F11._grid(a2, "x")
    save(fig, "F12_2_B4_static_surrogates")


def f3():
    P = A["PM97"]
    arms = [("raw", "raw (PM97)"), ("masked", "masked"), ("edited", "edited"), ("wu", "Wu")]
    fig, ax = plt.subplots(figsize=(W2, 2.4))
    rng = np.random.default_rng(12000)
    grp = {s["subject"]: s.get("group") for s in P["subjects"]}
    for ai, (arm, lab) in enumerate(arms):
        per = P["cells"][f"500|all|{arm}|R0"]["per_subject"]
        for gi, g in enumerate(("NSR", "CHF")):
            v = np.array([100 * x for s, x in per.items() if grp.get(s) == g])
            xpos = ai * 3 + gi
            ax.scatter(xpos + rng.uniform(-0.18, 0.18, len(v)), v, s=7, color=S[gi], marker=MK[gi], alpha=0.75, lw=0,
                       label=(f"{'nsrdb (healthy)' if g == 'NSR' else 'chfdb (CHF)'}") if ai == 0 else None)
            ax.plot([xpos - 0.3, xpos + 0.3], [v.mean()] * 2, color=INK, lw=1.2)
    ax.set_xticks([ai * 3 + 0.5 for ai in range(len(arms))], [lab for _, lab in arms])
    ax.set_ylabel("segments detected (%)")
    ax.set_ylim(-3, 103)
    ax.set_title("PM97 replication on nsrdb / chfdb: 40 × 500-beat segments per subject (R0; bar = group mean)", loc="left")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
    F11._grid(ax)
    save(fig, "F12_3_pm97")


def main():
    for fn in (f1, f2, f3):
        fn()
        print(fn.__name__, "ok")


if __name__ == "__main__":
    main()
