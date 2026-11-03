"""
Phase 11 runner (robustness analyses; resumable: one JSON line per task, finished ids skipped).

    python -m experiments.phase11_robustness.run11 --part prep,hrt,covar --workers 4      # guarded
    python -m experiments.phase11_robustness.run11 --part tier1 --workers 4               # guarded
    python -m experiments.phase11_robustness.run11 --part tier2 --workers 4               # guarded
    python -m experiments.phase11_robustness.run11 --pilot --part all --workers 4         # dev data / pilot seeds

The confirmatory run refuses to start unless experiments/phase11_robustness/PREREGISTRATION.md is committed,
unmodified and pushed (Phase 3 guard).  --pilot runs ~2 % of each part's tasks on DEVELOPMENT data (nsrdb,
chfdb) and pilot seeds 11900+ (never used by the confirmatory synthetic parts) for the runtime budget.

Parts (PREREGISTRATION.md sections in brackets)
  prep        masks of every window / segment (donor pool for A1 clean units)
  hrt         [A2.1] usable VPB tachograms and subject TO / TS (full recordings)
  covar       [A2.6, A3] window / segment NN covariates, abnormal-beat symbols
  sham_win    [A1] 5 sham masks per 512-interval window; masked + edited titration
  sham_seg    [A1] 2 sham masks per 12-min segment; masked + edited + Wu titration
  transplant  [A2.3] V0-V7 in clean confirmatory windows; raw titration (fixed seed per window)
  synth1      [A2.4] V0-V7 in N1 and 6 non-chaotic coupled-vdP regimes, seeds 11000-11099
  extmask     [A2.5] masked arm with the mask extended by 20 intervals after every abnormal beat
  rdg_verify  [A5.3] verification V1-V6 under readings R0-R7
  rdg_win     [A5.4] readings R1-R7: titration raw / masked / edited on all windows
  rdg_synth   [A5.4] readings R0-R7 on the key synthetic conditions, seeds 11000-11099
  rdg_seg     [A5.4] readings R1, R4: segments raw / masked / Wu
  seed_win    [B1] titration salts 1, 2 on all windows (raw / masked / edited)
  seed_seg    [B1] titration salts 1, 2 on all segments (raw / masked / Wu)
  k3seed      [B1] K3 with 2 new masked_growth_seed_entropy values
  k1seed      [B1] K1 with 2 new (lle_chaos_test_seed_entropy, random_seed) pairs
  b2          [B2] Phase 9 TEST groups (new seeds), K3 at 1/360 s and re-quantized at 1/128 s
  b3          [B3] variance-preserving spike-ins (K3)
  b4          [B4] additive spike-ins in a second base window per subject (K3)
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import time
import traceback
import zlib
from dataclasses import replace
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

from experiments.phase10_final import data10 as D  # noqa: E402
from experiments.phase10_final import harness10 as H  # noqa: E402
from experiments.phase10_final import run10 as R10  # noqa: E402
from experiments.phase11_robustness import harness11 as H11  # noqa: E402
from experiments.phase11_robustness import hrt11 as HRT  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
PREREG = HERE / "PREREGISTRATION.md"
METHODS_LISTED = ("TIT", "K3", "K1", "HRT")
TIER1 = ("prep", "hrt", "covar", "sham_win", "transplant", "synth1", "extmask", "rdg_verify", "rdg_synth",
         "rdg_win", "sham_seg", "rdg_seg")
TIER2 = ("seed_win", "k3seed", "b2", "b3", "seed_seg", "k1seed", "b4")
SEEDS = range(11000, 11100)
PILOT_SEEDS = range(11900, 11902)
READINGS_ALL = ("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7")
READINGS_SEG = ("R1", "R4")
SALTS = (1, 2)
K3_ENTROPIES = (11101, 11102)
K1_SEEDS = ((11201, 11211), (11202, 11212))
SYNTH1_CONDS = (("null", "N1_linear_rr"),) + tuple(("vdp", r) for r in H.VDP_NONCHAOTIC)
RDG_SYNTH = (("null", "N1_linear_rr"), ("null", "E3_couplets_5pct"), ("null", "E4_runs"),
             ("null", "N5_linear_rr_warped"), ("vdp", 26))
N_SYN = 512
FS_SYN = 360.0


def guard():
    from experiments.phase3_lle import run_phase3 as R3
    R3.PREREG = PREREG
    R3.check_preregistration(list(METHODS_LISTED))


class Ctx:
    def __init__(self, pilot):
        self.pilot = pilot
        self.phase = "dev" if pilot else "conf"
        self.root = D.DEV_DIR if pilot else D.CONF_DIR
        self.dbs = ("nsrdb", "chfdb") if pilot else D.CONF_DBS
        self.out = HERE / "results" / ("pilot" if pilot else "conf")
        self.seeds = PILOT_SEEDS if pilot else SEEDS


def subjects(ctx):
    out = []
    for db in ctx.dbs:
        for si, rec in enumerate(D.records(db, ctx.root)):
            out.append((db, rec, si))
    return out


def windows(ctx):
    W = []
    for db, rec, si in subjects(ctx):
        t = R10.series(ctx.root, db, rec)[0]
        for j, i0 in H.window_starts(t).items():
            W.append({"db": db, "rec": rec, "si": si, "j": j, "i0": i0, "wid": f"real|{db}|{rec}|{j}"})
    return W


def win_data(ctx, w, L=H.L_PRIMARY):
    return R10._window(ctx.root, {"db": w["db"], "rec": w["rec"], "i0": w["i0"]}, L)


def seg_units(ctx, db, rec):
    """Phase 10 run_seg preprocessing of every analysable 12-min segment: (seg, rr, mask, forbidden, t_start)."""
    t, lab, sym, fs, info = R10.series(ctx.root, db, rec)
    out = []
    for s in H.segments(t, lab):
        if s["gap"]:
            continue
        rr_all = np.diff(t[s["i0"]:s["i1"] + 1])
        lb = lab[s["i0"]:s["i1"] + 1]
        sy = sym[s["i0"]:s["i1"] + 1]
        keep = (rr_all >= H.LO) & (rr_all <= H.HI)
        m_all = H.k3_mask(lb)
        out.append({"seg": s["seg"], "rr": rr_all[keep], "mask": m_all[keep], "fb": H11.forbidden(m_all, lb)[keep],
                    "t_start": float(t[s["i0"]]), "burden": int(np.sum(~np.isin(sy, list(D.NORMAL))))})
    return out


def load_jsonl(path):
    R = {}
    if path.exists():
        for line in open(path):
            try:
                r = json.loads(line)
                R[r["id"]] = r
            except Exception:                                   # noqa: BLE001
                pass
    return R


def hrt_profiles(ctx):
    """Group HRT-empirical profile: subject-weighted mean normalised tachogram over subjects with >= 1 usable
    VPB (results/<phase>/hrt.jsonl)."""
    R = load_jsonl(ctx.out / "hrt.jsonl")
    prof = {}
    for g in ("NSR", "CHF"):
        hs = [r["summary"]["h"] for r in R.values() if r.get("group") == g and r.get("summary", {}).get("n_vpb", 0) >= 1]
        if hs:
            prof[g] = {"h": np.mean(hs, axis=0).tolist(), "n_subjects": len(hs)}
    return prof


# ------------------------------------------------------------------------------------------- tasks
def tasks(ctx, part):
    T = []
    pick = (lambda L: L[::50]) if ctx.pilot else (lambda L: L)
    if part in ("prep",):
        return [{"id": "prep", "part": "prep"}]
    if part in ("hrt", "covar", "sham_seg", "seed_seg", "rdg_seg"):
        subs = subjects(ctx)
        if ctx.pilot and part in ("sham_seg", "seed_seg", "rdg_seg"):
            subs = subs[::16]
        base = [{"db": db, "rec": rec, "si": si, "subject": R10.subject_of(db, rec), "group": R10.GROUP[db]}
                for db, rec, si in subs]
        if part == "seed_seg":
            return [{**b, "id": f"{part}|{b['db']}|{b['rec']}|s{s}", "part": part, "salt": s} for s in SALTS for b in base]
        if part == "rdg_seg":
            return [{**b, "id": f"{part}|{b['db']}|{b['rec']}|{rd}", "part": part, "reading": rd}
                    for rd in READINGS_SEG for b in base]
        return [{**b, "id": f"{part}|{b['db']}|{b['rec']}", "part": part} for b in base]
    if part in ("sham_win", "transplant", "extmask", "rdg_win", "seed_win", "k3seed", "k1seed"):
        W = windows(ctx)
        ref = load_jsonl(HERE.parent / "phase10_final" / "results" / ctx.phase / "real.jsonl")
        for w in W:
            w["subject"], w["group"] = R10.subject_of(w["db"], w["rec"]), R10.GROUP[w["db"]]
            r = ref.get(w["wid"], {})
            w["burden"], w["n_masked"] = r.get("burden"), r.get("n_masked")
        if part == "sham_win":
            T = [{**w, "id": f"{part}|{w['wid']}", "part": part} for w in pick(W)]
        elif part == "transplant":
            prof = hrt_profiles(ctx)
            T = [{**w, "id": f"{part}|{w['wid']}", "part": part, "h_emp": prof.get(w["group"], {}).get("h")}
                 for w in pick([w for w in W if w["burden"] == 0 and w["n_masked"] == 0])]
        elif part == "extmask":
            T = [{**w, "id": f"{part}|{w['wid']}", "part": part}
                 for w in pick([w for w in W if w["burden"] is not None and 1 <= w["burden"] <= 4])]
        elif part == "rdg_win":
            T = [{**w, "id": f"{part}|{w['wid']}|{rd}", "part": part, "reading": rd}
                 for rd in READINGS_ALL[1:] for w in pick(W)]
        elif part == "seed_win":
            T = [{**w, "id": f"{part}|{w['wid']}|s{s}", "part": part, "salt": s} for s in SALTS for w in pick(W)]
        elif part == "k3seed":
            T = [{**w, "id": f"{part}|{w['wid']}|e{e}", "part": part, "entropy": e} for e in K3_ENTROPIES for w in pick(W)]
        elif part == "k1seed":
            T = [{**w, "id": f"{part}|{w['wid']}|e{a}_{b}", "part": part, "lle_entropy": a, "random_seed": b}
                 for a, b in K1_SEEDS for w in pick(W)]
        return T
    if part == "synth1":
        prof = hrt_profiles(ctx).get("NSR", {}).get("h")
        for kind, c in SYNTH1_CONDS:
            for s in ctx.seeds:
                T.append({"id": f"{part}|{kind}|{c}|{s}", "part": part, "kind": kind, "cond": c, "seed": s, "h_emp": prof})
        return T
    if part == "rdg_synth":
        for rd in READINGS_ALL:
            for kind, c in RDG_SYNTH:
                for s in ctx.seeds:
                    T.append({"id": f"{part}|{rd}|{kind}|{c}|{s}", "part": part, "reading": rd, "kind": kind,
                              "cond": c, "seed": s})
        return T
    if part == "rdg_verify":
        from experiments.phase10_final import verify_titration as VT
        VTT = VT.tasks()
        if ctx.pilot:
            VTT = VTT[::25]
        return [{"id": f"{part}|{rd}|{n}|{a}", "part": part, "reading": rd, "system": n, "arg": a}
                for rd in READINGS_ALL for n, a in VTT]
    if part == "b2":
        from experiments.phase9_waveform import run_test9 as RT
        RT.S100 = ctx.seeds if not ctx.pilot else range(11900, 11901)
        RT.S30 = range(11000, 11030) if not ctx.pilot else range(11900, 11901)
        RT.S10 = range(11000, 11010) if not ctx.pilot else range(11900, 11901)
        TT = RT.tasks()
        if ctx.pilot:
            TT = TT[::3]
        return [{**t, "id": f"b2|{t['task_id']}", "part": "b2"} for t in TT]
    if part in ("b3", "b4"):
        for db, rec, si in subjects(ctx):
            t = R10.series(ctx.root, db, rec)[0]
            ws = H.window_starts(t)
            if not ws:
                continue
            jb = R10.base_window(si, ws)
            if part == "b4":
                jb = second_base_window(si, ws, jb)
                if jb is None:
                    continue
            for fam, lst in H.FAMILIES.items():
                for p in range(len(lst)):
                    T.append({"id": f"{part}|{db}|{rec}|{jb}|{fam}|{p}", "part": part, "db": db, "rec": rec, "j": jb,
                              "i0": ws[jb], "family": fam, "pidx": p, "si": si})
        return pick(T) if not ctx.pilot else T[::50]
    raise ValueError(part)


def second_base_window(si, ws, jb):
    """B4: j' = 1 + ((si + 6) mod 12), the next available window cyclically, never the Phase 10 base window."""
    for k in range(H.N_PER):
        j = 1 + (si + 6 + k) % H.N_PER
        if j in ws and j != jb:
            return j
    return None


# ------------------------------------------------------------------------------------------- execution
def tit(rr, arm, mask, salt=0):
    """Compact Phase 10 titration (harness10.titration) with an explicit interval mask."""
    o = H.titration(rr, None, arm, salt=salt, imask=mask)
    return {k: o.get(k) for k in ("eligible", "analysable", "positive", "NL", "n_rows", "nn_fraction")}


def tit_wu(rr, mask, salt=0):
    nn = H.wu_nn(rr, mask)
    if len(nn) < H.SEG_MIN_N:
        return {"eligible": False, "n": int(len(nn)), "positive": False, "NL": 0.0}
    o = tit(nn, "raw", None, salt)
    o["n"] = int(len(nn))
    return o


def tit_reading(rr, arm, mask, reading):
    """harness10.titration's arms and seeding (crc of the data, salt 0) under titration11 reading."""
    from experiments.phase11_robustness import titration11 as T11
    x = np.asarray(rr, float)
    rng = np.random.default_rng(np.random.SeedSequence([H.ENTROPY, H.crc(x), 77, 0]))
    if arm == "raw":
        o = T11.titrate_reading(x, rng, reading)
    elif arm == "masked":
        o = T11.titrate_reading(x, rng, reading, mask=np.asarray(mask, bool))
    elif arm == "edited":
        xe, nnf = H.edit(x, imask=mask)
        if xe is None or nnf < H.NN_EDIT_MIN:
            return {"eligible": False, "positive": False, "NL": 0.0, "nn_fraction": nnf}
        o = T11.titrate_reading(xe, rng, reading)
    elif arm == "wu":
        nn = H.wu_nn(x, mask)
        if len(nn) < H.SEG_MIN_N:
            return {"eligible": False, "positive": False, "NL": 0.0, "n": int(len(nn))}
        rng = np.random.default_rng(np.random.SeedSequence([H.ENTROPY, H.crc(nn), 77, 0]))
        o = T11.titrate_reading(nn, rng, reading)
    else:
        raise ValueError(arm)
    b = o["base"]
    return {"eligible": True, "analysable": bool(b.get("analysable", False)), "positive": bool(o["positive"]),
            "NL": float(o["NL"]), "n_rows": b.get("n_rows")}


def run_prep(ctx, task):
    win, seg = {}, {}
    for w in windows(ctx):
        d = win_data(ctx, w)
        m = H.k3_mask(d["lab"])
        win[w["wid"]] = {"group": R10.GROUP[w["db"]], "n": len(m), "idx": np.nonzero(m)[0].tolist()}
    for db, rec, si in subjects(ctx):
        for u in seg_units(ctx, db, rec):
            seg[f"{db}|{rec}|{u['seg']}"] = {"group": R10.GROUP[db], "n": len(u["mask"]),
                                             "idx": np.nonzero(u["mask"])[0].tolist()}
    p = ctx.out / "prep_masks.json"
    json.dump({"win": win, "seg": seg}, open(p, "w"))
    return {"n_win": len(win), "n_seg": len(seg)}


_PREP = {}


def prep(ctx):
    if ctx.phase not in _PREP:
        P = json.load(open(ctx.out / "prep_masks.json"))
        pools = {}
        for unit in ("win", "seg"):
            for g in ("NSR", "CHF"):
                pools[(unit, g)] = sorted(k for k, v in P[unit].items() if v["group"] == g and v["idx"])
        _PREP[ctx.phase] = (P, pools)
    return _PREP[ctx.phase]


def donor_masks(ctx, unit, group, n, k, rng):
    P, pools = prep(ctx)
    pool = pools[(unit, group)]
    pick = rng.choice(len(pool), size=min(k, len(pool)), replace=False)
    out = []
    for i in pick:
        d = P[unit][pool[int(i)]]
        m = np.zeros(d["n"], bool)
        m[d["idx"]] = True
        out.append((pool[int(i)], H11.borrow_mask(m, n)))
    return out


def run_hrt(ctx, task):
    t, lab, sym, fs, info = R10.series(ctx.root, task["db"], task["rec"])
    tg = HRT.tachograms(t, sym, lab)
    s = HRT.subject_hrt(tg)
    return {"subject": task["subject"], "group": task["group"], "n_V": int(np.sum(np.asarray(sym) == "V")),
            "summary": s, "TO_per_vpb": [g["TO"] for g in tg]}


def run_covar(ctx, task):
    t, lab, sym, fs, info = R10.series(ctx.root, task["db"], task["rec"])
    win = []
    for j, i0 in H.window_starts(t).items():
        rr = np.diff(t[i0:i0 + 513])
        lb, sy = lab[i0:i0 + 513], sym[i0:i0 + 513]
        m = H.k3_mask(lb)
        ab = [str(s) for s in sy if s not in D.NORMAL]
        win.append({"wid": f"real|{task['db']}|{task['rec']}|{j}", "j": j,
                    **H11.nn_covariates(rr, m, float(t[i0] - t[0])), "abnormal_symbols": ab,
                    "n_rr_flag_only": int(np.sum((lb != "N") & np.isin(sy, list(D.NORMAL))))})
    seg = []
    for u in seg_units(ctx, task["db"], task["rec"]):
        seg.append({"seg": u["seg"], **H11.nn_covariates(u["rr"], u["mask"], u["t_start"] - float(t[0]))})
    return {"subject": task["subject"], "group": task["group"], "windows": win, "segments": seg}


def run_sham_win(ctx, task):
    d = win_data(ctx, task)
    rr, lab = d["rr"], d["lab"]
    m = H.k3_mask(lab)
    rng = np.random.default_rng(H11.ss("sham_win", H11.key_of(task["wid"])))
    out = {"subject": task["subject"], "group": task["group"], "j": task["j"], "n_masked": int(m.sum()),
           "burden": int(np.sum(~np.isin(d["sym"], list(D.NORMAL)))), "shams": []}
    if m.any():
        out["kind"] = "own"
        shams = [(o, s) for o, s in H11.sham_masks(m, H11.forbidden(m, lab), H11.N_SHAM_WIN, rng)]
    else:
        out["kind"] = "donor"
        shams = donor_masks(ctx, "win", task["group"], len(rr), H11.N_SHAM_WIN, rng)
    out["n_sham"] = len(shams)
    for o, s in shams:
        out["shams"].append({"shift_or_donor": o, "n_sham_masked": int(s.sum()), "masked": tit(rr, "masked", s),
                             "edited": tit(rr, "edited", s)})
    return out


def run_sham_seg(ctx, task):
    out = {"subject": task["subject"], "group": task["group"], "segments": []}
    for u in seg_units(ctx, task["db"], task["rec"]):
        rr, m = u["rr"], u["mask"]
        rng = np.random.default_rng(H11.ss("sham_seg", H11.key_of(task["db"], task["rec"], u["seg"])))
        if m.any():
            kind, shams = "own", H11.sham_masks(m, u["fb"], H11.N_SHAM_SEG, rng)
        else:
            kind, shams = "donor", donor_masks(ctx, "seg", task["group"], len(rr), H11.N_SHAM_SEG, rng)
        rec = {"seg": u["seg"], "kind": kind, "n_masked": int(m.sum()), "n_sham": len(shams), "shams": []}
        for o, s in shams:
            rec["shams"].append({"shift_or_donor": o, "n_sham_masked": int(s.sum()), "masked": tit(rr, "masked", s),
                                 "edited": tit(rr, "edited", s), "wu": tit_wu(rr, s)})
        out["segments"].append(rec)
    return out


def _transplant(x, h_emp, rng_ev, rng_tit_seed, fs, quantize=True):
    """V0-V7 in series x at one seeded position; raw titration with the same rng seed for every variant."""
    from experiments.phase10_final import titration10 as T
    j = HRT.event_position(len(x), rng_ev)
    dr = HRT.draws_for(rng_ev)
    res = {"pos": j, "draws": dr}
    for v in HRT.VARIANTS:
        if v in ("V2", "V7") and h_emp is None:
            res[v] = {"error": "no HRT-empirical profile"}
            continue
        xv = HRT.insert_event(x, j, v, dr, h_emp=h_emp)
        if v != "V0" and quantize:
            xv = H.quantize_times(0.0, xv, fs)
        o = T.titrate(xv, np.random.default_rng(rng_tit_seed))
        res[v] = {"positive": bool(o["positive"]), "NL": float(o["NL"]), "min_rr": float(np.min(xv)),
                  "max_rr": float(np.max(xv))}
    return res


def run_transplant(ctx, task):
    d = win_data(ctx, task)
    key = H11.key_of(task["wid"])
    rng_ev = np.random.default_rng(H11.ss("transplant", key, 1))
    seed = H11.ss("transplant", key, 2)
    return {"subject": task["subject"], "group": task["group"], "j": task["j"],
            **_transplant(d["rr"], task["h_emp"], rng_ev, seed, d["fs"])}


def synth_series(kind, cond, seed, n=N_SYN):
    """Phase 10 synthetic window (harness10.synthetic_window): Phase 5-6 null or non-chaotic coupled vdP."""
    if kind == "vdp":
        return H.synthetic_window("vdp", (cond, "none"), seed, n)
    return H.synthetic_window("null", cond, seed, n)


def run_synth1(ctx, task):
    x, _ = synth_series(task["kind"], task["cond"], task["seed"])
    cc = 0 if task["kind"] == "null" else int(task["cond"])
    rng_ev = np.random.default_rng(H11.ss("synth1", cc, task["seed"], 1))
    seed = H11.ss("synth1", cc, task["seed"], 2)
    return {"kind": task["kind"], "cond": task["cond"], "seed": task["seed"],
            **_transplant(x, task["h_emp"], rng_ev, seed, FS_SYN)}


def run_extmask(ctx, task):
    d = win_data(ctx, task)
    rr, lab = d["rr"], d["lab"]
    m = H.k3_mask(lab)
    ext = H11.extended_mask(m, lab)
    return {"subject": task["subject"], "group": task["group"], "j": task["j"], "n_masked": int(m.sum()),
            "n_ext_masked": int(ext.sum()), "masked_std": tit(rr, "masked", m), "masked_ext": tit(rr, "masked", ext)}


def run_rdg_verify(ctx, task):
    from experiments.phase10_final import verify_titration as VT
    from experiments.phase11_robustness import titration11 as T11
    name, arg = task["system"], task["arg"]
    if name == "logistic":
        y, extra = VT.logistic(arg), {"LE": VT.logistic_le(arg)}
    elif name == "henon":
        y, extra = VT.henon(), {}
    elif name == "lorenz":
        y, extra = VT.lorenz(arg), {}
    elif name == "mg":
        y, extra = VT.mackey_glass(arg), {}
    else:
        y, extra = VT.controls(name, arg), {}
    rng = np.random.default_rng(np.random.SeedSequence([20261003, zlib.crc32(f"{name}|{arg}".encode())]))
    o = T11.titrate_reading(y, rng, task["reading"])
    return {"reading": task["reading"], "system": name, "arg": arg, "NL": o["NL"], "positive": o["positive"], **extra}


def run_rdg_win(ctx, task):
    d = win_data(ctx, task)
    m = H.k3_mask(d["lab"])
    return {"subject": task["subject"], "group": task["group"], "j": task["j"], "reading": task["reading"],
            **{a: tit_reading(d["rr"], a, m, task["reading"]) for a in ("raw", "masked", "edited")}}


def run_rdg_synth(ctx, task):
    x, _ = synth_series(task["kind"], task["cond"], task["seed"])
    return {"reading": task["reading"], "kind": task["kind"], "cond": task["cond"], "seed": task["seed"],
            "raw": tit_reading(x, "raw", None, task["reading"])}


def run_rdg_seg(ctx, task):
    out = {"subject": task["subject"], "group": task["group"], "reading": task["reading"], "segments": []}
    for u in seg_units(ctx, task["db"], task["rec"]):
        out["segments"].append({"seg": u["seg"], **{a: tit_reading(u["rr"], a, u["mask"], task["reading"])
                                                    for a in ("raw", "masked", "wu")}})
    return out


def run_seed_win(ctx, task):
    d = win_data(ctx, task)
    m = H.k3_mask(d["lab"])
    return {"subject": task["subject"], "group": task["group"], "j": task["j"], "salt": task["salt"],
            **{a: tit(d["rr"], a, m, task["salt"]) for a in ("raw", "masked", "edited")}}


def run_seed_seg(ctx, task):
    out = {"subject": task["subject"], "group": task["group"], "salt": task["salt"], "segments": []}
    for u in seg_units(ctx, task["db"], task["rec"]):
        out["segments"].append({"seg": u["seg"], "raw": tit(u["rr"], "raw", None, task["salt"]),
                                "masked": tit(u["rr"], "masked", u["mask"], task["salt"]),
                                "wu": tit_wu(u["rr"], u["mask"], task["salt"])})
    return out


def run_k3seed(ctx, task):
    import final_pipeline as fp
    d = win_data(ctx, task)
    cfg = replace(fp.CFG, masked_growth_seed_entropy=int(task["entropy"]))
    r = H11.k3_diag(d["rr"], d["lab"], cfg)
    return {"subject": task["subject"], "group": task["group"], "j": task["j"], "entropy": task["entropy"],
            "K3": {k: r.get(k) for k in ("detected", "analysable", "zmax", "G", "degenerate", "sd_at_argmax")}}


def run_k1seed(ctx, task):
    from experiments.phase7_mitbih import detector as DT
    d = win_data(ctx, task)
    cfg = replace(DT.DETECTOR, lle_chaos_test_seed_entropy=int(task["lle_entropy"]), random_seed=int(task["random_seed"]))
    r = DT.evaluate(np.asarray(d["rr"], float), cfg)
    return {"subject": task["subject"], "group": task["group"], "j": task["j"],
            "seeds": [task["lle_entropy"], task["random_seed"]], "K1": {k: r.get(k) for k in H.K1_KEYS}}


def run_b2(ctx, task):
    from experiments.phase9_waveform import run_test9 as RT
    w = RT.window(task)
    rr = np.asarray(w["rr"], float)
    types = np.asarray(w["beat_types"])
    q = H.quantize_times(0.0, rr, 128.0)
    keep = ("detected", "analysable", "zmax", "G", "degenerate", "sd_at_argmax", "mad_nn", "sdnn_nn")
    out = {k: task[k] for k in ("group", "source", "name", "label", "ectopy", "level", "seed") if k in task}
    out["family"] = task.get("family")
    out["n_nonN"] = int(np.sum(types != "N"))
    out["K3_360"] = {k: v for k, v in H11.k3_diag(rr, types).items() if k in keep}
    out["K3_128"] = {k: v for k, v in H11.k3_diag(q, types).items() if k in keep}
    return out


def run_spike11(ctx, task):
    w = win_data(ctx, task)
    rr, lab, sym, fs = w["rr"], w["lab"], w["sym"], w["fs"]
    fam, p = task["family"], task["pidx"]
    dbc = {"nsrdb": 1, "chfdb": 2, "nsr2db": 3, "chf2db": 4}[task["db"]]
    c = H.chaos_component(fam, p, len(rr), [dbc, task["si"], task["j"]])
    out = {"subject": R10.subject_of(task["db"], task["rec"]), "group": R10.GROUP[task["db"]], "j": task["j"],
           "family": fam, "param": H.FAMILIES[fam][p][0], "lam": H.FAMILIES[fam][p][2], "fs": fs,
           "burden": int(np.sum(~np.isin(sym, list(D.NORMAL)))), "designs": {}}
    for f in H.F_GRID:
        if task["part"] == "b3":
            q, meta = H11.spike_varpres(rr, lab, c, f, fs)
            name = f"varpres_{f}"
        else:
            q, meta = H.spike_additive(rr, lab, c, f, fs)
            name = f"add_{f}"
        out["designs"][name] = {"f": f, **meta, "K3": H.k3(q, lab)}
    return out


RUN = {"prep": run_prep, "hrt": run_hrt, "covar": run_covar, "sham_win": run_sham_win, "sham_seg": run_sham_seg,
       "transplant": run_transplant, "synth1": run_synth1, "extmask": run_extmask, "rdg_verify": run_rdg_verify,
       "rdg_win": run_rdg_win, "rdg_synth": run_rdg_synth, "rdg_seg": run_rdg_seg, "seed_win": run_seed_win,
       "seed_seg": run_seed_seg, "k3seed": run_k3seed, "k1seed": run_k1seed, "b2": run_b2, "b3": run_spike11,
       "b4": run_spike11}


def execute(args):
    pilot, task = args
    ctx = Ctx(pilot)
    t0 = time.time()
    out = {"id": task["id"], "part": task["part"]}
    try:
        out.update(RUN[task["part"]](ctx, task))
    except Exception as exc:                                    # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
        out["traceback"] = traceback.format_exc(limit=6)
    out["t_total"] = time.time() - t0
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="tier1")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    if not a.pilot:
        guard()
    ctx = Ctx(a.pilot)
    ctx.out.mkdir(parents=True, exist_ok=True)
    parts = {"tier1": TIER1, "tier2": TIER2, "all": TIER1 + TIER2}.get(a.part, tuple(a.part.split(",")))
    for part in parts:
        path = ctx.out / f"{part}.jsonl"
        done = set(load_jsonl(path))
        T = [t for t in tasks(ctx, part) if t["id"] not in done]
        if a.limit:
            T = T[:a.limit]
        print(f"{part}: {len(T)} tasks to run ({len(done)} done)", flush=True)
        if not T:
            continue
        t0 = time.time()
        T.sort(key=lambda t: zlib.crc32(t["id"].encode()))      # interleave subjects / conditions
        with Pool(a.workers, maxtasksperchild=200) as pool, open(path, "a") as fh:
            for i, r in enumerate(pool.imap_unordered(execute, [(a.pilot, t) for t in T], chunksize=1)):
                fh.write(json.dumps(r, default=R10._default) + "\n")
                fh.flush()
                if (i + 1) % 50 == 0 or i + 1 == len(T):
                    el = time.time() - t0
                    print(f"  {part} {i + 1}/{len(T)} {el / 60:.1f} min; eta {el / (i + 1) * (len(T) - i - 1) / 60:.1f} min",
                          flush=True)
    print("RUN11 DONE", parts, flush=True)


if __name__ == "__main__":
    main()
