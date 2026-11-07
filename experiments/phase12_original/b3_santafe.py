"""Phase 12 B3 (optional): BP96 Fig. 2a on Santa Fe series A, if obtainable (PREREGISTRATION 3, B3).
Source used: public GitHub repository shubhomoydas/ad_examples (commit 4f20369), datasets/SantaFeCompetition/SF_A.dat
(1,000 points, the competition's data set A), cloned read-only outside this repository; the data are not committed.
Series E: not found on PyPI / GitHub through the allowed routes -> not obtainable.

    PYTHONPATH=. python -m experiments.phase12_original.b3_santafe /path/to/SF_A.dat
"""
import json
import sys
from pathlib import Path

import numpy as np

from experiments.phase12_original import titration12 as T12

OUT = Path(__file__).resolve().parent / "results" / "conf" / "b3_status.json"


def main(path):
    y = np.array([float(v) for v in Path(path).read_text().replace("\n", ",").split(",") if v.strip()])
    res = {"series_A": {"source": "github.com/shubhomoydas/ad_examples@4f20369 datasets/SantaFeCompetition/SF_A.dat",
                        "n": int(len(y)), "first_values": [float(v) for v in y[:5]],
                        "BP96": {"F_calc": 10.23, "F_crit": 1.16, "df": [968, 944], "MW_Z": -24.60}},
           "series_E": {"status": "not obtainable (no copy found on PyPI or GitHub via the allowed routes)"}}
    for rd in ("R0", "R0-BP96F", "R0-BP96"):
        o = T12.vwk12(y, name=rd, details=True)
        res["series_A"][rd] = {k: o.get(k) for k in ("nonlinear", "d_opt_gt1", "F_calc_plain", "F", "df1", "df2", "F_crit",
                                                    "MW_Z", "kappa_nl", "d_nl", "r_lin", "r_nl", "test", "p")}
        res["series_A"][rd]["inv_F_calc"] = 1 / o["F_calc_plain"]
    res["status"] = "series A run; series E not obtainable"
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res["series_A"]["R0"], default=float))


if __name__ == "__main__":
    main(sys.argv[1])
