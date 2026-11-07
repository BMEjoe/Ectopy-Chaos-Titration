"""
Phase 11 data helper: the nstdb noise records Phase 9 uses (bw, ma, em; .hea + .dat), downloaded from
PhysioNet nstdb 1.0.0 into experiments/phase9_waveform/nstdb_data/ (gitignored) and SHA-256 checked
against PhysioNet's SHA256SUMS.txt (same procedure as phase10_final.data10.download).

    python -m experiments.phase11_robustness.data11 --download nstdb
    python -m experiments.phase11_robustness.data11 --verify          # all Phase 10/11 data
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib

from experiments.phase10_final import data10 as D

ROOT = pathlib.Path(__file__).resolve().parents[2]
NSTDB = ROOT / "experiments" / "phase9_waveform" / "nstdb_data"
URL = "https://physionet.org/files/nstdb/1.0.0/"


def download_nstdb():
    NSTDB.mkdir(parents=True, exist_ok=True)
    sums = D._get(URL + "SHA256SUMS.txt").decode()
    (NSTDB / "SHA256SUMS.txt").write_text(sums)
    table = {ln.split(None, 1)[1].strip(): ln.split(None, 1)[0] for ln in sums.splitlines() if ln.strip()}
    manifest = {}
    for rec in ("bw", "ma", "em"):
        for ext in ("hea", "dat"):
            name = f"{rec}.{ext}"
            p = NSTDB / name
            if not p.exists():
                p.write_bytes(D._get(URL + name))
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            if table.get(name) != h:
                raise SystemExit(f"SHA-256 mismatch for nstdb/{name}")
            manifest[name] = h
    (NSTDB / "MANIFEST.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def verify_all():
    """Re-hash every downloaded file against the stored SHA256SUMS.txt of its database."""
    out = {}
    dirs = [(D.CONF_DIR / db) for db in D.CONF_DBS] + [(D.DEV_DIR / db) for db in D.DEV_DBS] + [NSTDB]
    for d in dirs:
        if not (d / "SHA256SUMS.txt").exists():
            out[str(d.relative_to(ROOT))] = "missing"
            continue
        table = {ln.split(None, 1)[1].strip(): ln.split(None, 1)[0]
                 for ln in (d / "SHA256SUMS.txt").read_text().splitlines() if ln.strip()}
        files = [p for p in d.iterdir() if p.name not in ("SHA256SUMS.txt", "MANIFEST.json")]
        bad = [p.name for p in files if table.get(p.name) != hashlib.sha256(p.read_bytes()).hexdigest()]
        out[str(d.relative_to(ROOT))] = {"n_files": len(files), "mismatch": bad}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", choices=["nstdb"])
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    if a.download == "nstdb":
        print("nstdb", len(download_nstdb()), "files; SHA-256 OK")
    if a.verify:
        print(json.dumps(verify_all(), indent=1))


if __name__ == "__main__":
    main()
