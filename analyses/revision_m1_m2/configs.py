"""
Revision M1/M2 (EXPLORATORY): named configurations (PLAN.md Section 3.2).  R0-ORIG is built by the fixed rule of
PLAN.md 3.5 from the single-switch benchmark results and stored in results/revision_m1_m2/r0orig.json.
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from analyses.revision_m1_m2.vwkcfg import Cfg

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "revision_m1_m2"
R0ORIG_FILE = RES / "r0orig.json"

# switch id -> Cfg field changes (one deviation each)
SWITCHES = {
    "S1": {"n_mode": "length"},
    "S2": {"f_stat": "plain"},
    "S3": {"gate": True},
    "S3b": {"test": "MW"},
    "S4": {"surrogate": "FT"},
    "S5a": {"d_max": 8, "m_cap": 84},
    "S5b": {"kappa_max": 7},
    "S6": {"rows": "zeropad"},
    "S7": {"decision": "test"},
}
SWITCH_LABEL = {
    "S1": "N = series length in C(r) and df",
    "S2": "plain F = RSS_lin/RSS_nl",
    "S3": "chi-square gate -> Mann-Whitney",
    "S3b": "Mann-Whitney only",
    "S4": "FT surrogate step",
    "S5a": "degree up to 8 (M <= 84)",
    "S5b": "kappa <= 7 (M = 120)",
    "S6": "zero-padded start (all points targets)",
    "S7": "test-only decision (no d_opt > 1)",
}
GROUPS = {"S3": "test", "S3b": "test", "S5a": "range", "S5b": "range"}
SINGLE = tuple(SWITCHES)
COMPOSITES = {"R0-BP96F": ("S1", "S2"), "R0-BP96": ("S1", "S2", "S3"), "R0-TEXT": ("S1", "S2", "S3", "S4")}


def build(switches):
    kw = {}
    for s in switches:
        kw.update(SWITCHES[s])
    return replace(Cfg(**kw), switches=tuple(switches))


def r0orig_switches():
    if not R0ORIG_FILE.exists():
        raise FileNotFoundError("R0-ORIG not built yet (run the ablation analysis first)")
    return tuple(json.load(open(R0ORIG_FILE))["switches"])


def get(name):
    """R0 | S1..S7 | R0-BP96F | R0-BP96 | R0-TEXT | R0-ORIG | R0-ORIG-minus-<switch>."""
    if name == "R0":
        return Cfg()
    if name in SWITCHES:
        return build((name,))
    if name in COMPOSITES:
        return build(COMPOSITES[name])
    if name == "R0-ORIG":
        return build(r0orig_switches())
    if name.startswith("R0-ORIG-minus-"):
        drop = name.split("R0-ORIG-minus-", 1)[1]
        return build(tuple(s for s in r0orig_switches() if s != drop))
    raise KeyError(name)


def is_surrogate(name):
    return bool(get(name).surrogate)
