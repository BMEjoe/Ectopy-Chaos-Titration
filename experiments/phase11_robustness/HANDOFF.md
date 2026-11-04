# Phase 11 handoff (living document; update after every step)

Read this first when resuming. Phase 11 = registered SECONDARY robustness analyses for the manuscript
"Ectopic Beats Masquerade as Chaos in Heart Rate ...". Never merge into main, never open a pull request.
Commit and push after every step.

**Branches.** Work on `phase11-robustness` (created from main = Phase 10 head `1147f76`; the user asked for this
name). The session's designated branch `claude/beautiful-hawking-i1ommp` is kept identical (every push goes to both:
`git push origin HEAD:phase11-robustness HEAD:claude/beautiful-hawking-i1ommp`). The runner guard checks the
CURRENT branch's origin copy, so run from `phase11-robustness`.

## Task in brief (user, Phase 11)
- Integrity rule: never tune / re-select / re-define to rescue a claim; report weakening results prominently.
- Ground rules: frozen Phase 10 files; new code only in `experiments/phase11_robustness/`; final_pipeline.py opt-in
  only (not touched so far); seeds start at 11000; PREREGISTRATION.md committed + pushed before any analysis on
  confirmatory data (runner guard); amendments dated.
- Tier 1: A1 sham-handling control; A2 single-beat / HRT (measurement, generator, transplant V0-V7, synthetic,
  extended masking, beat type); A3 confounding (GEE + covariates, conditional logit); A4 original-study metrics;
  A5 titration readings R0-R7 (+ R8 only with the BP96 PDF), verification V1-V6 per reading, claims matrix.
- Tier 2: B1 seeds; B2 K3 at 1/128 s (Phase 9 TEST groups, degeneracy diagnostic); B3 variance-preserving
  spike-ins; B4 second base window.
- Tier 3: C1 paper_tables.py (+ manuscript check if the draft PDF exists); C2 figures F1-F7.
- Report: docs/PHASE11_ROBUSTNESS_RESULTS.md.

## Environment (rebuilt 2026-11-03 in a fresh container)
- Recipe `experiments/phase10_final/requirements-phase10.txt`, executed step by step (neurokit2 first, then
  pandas 3.0.6 over it, then statsmodels). Python **3.13.18** (Phase 10: 3.13.14), numpy 2.1.3, scipy 1.18.1,
  pandas 3.0.6, statsmodels 0.15.0, wfdb 4.3.1, numba 0.68.0, scikit-learn 1.9.1, matplotlib 3.11.2, patsy 1.0.3,
  wrapt 2.5.2 (Phase 10: 2.5.0). `/root/venv313/bin/python`, `OMP_NUM_THREADS=1`.
- pytest baseline (2026-11-03): **462 passed, 1 failed** (pre-existing
  `test_upo_stability.py::test_hybrid_period1_stability_matches_analytic_henon_multipliers[1-prl_norm]`). Same as Phase 10.
- Machine: 4 CPUs, 15 GB RAM.

## Data (gitignored; re-downloaded 2026-11-03, SHA-256 OK)
- `python -m experiments.phase10_final.run10 --phase conf --download` (guarded by the Phase 10 prereg): nsr2db 54,
  chf2db 29.
- `python -m experiments.phase10_final.data10 --download dev`: mitdb 48, nsrdb 18, chfdb 15.
- `python -m experiments.phase11_robustness.data11 --download nstdb` (new helper; Phase 9 had no script): bw, ma,
  em (.hea, .dat) into `experiments/phase9_waveform/nstdb_data/`.
- `python -m experiments.phase11_robustness.data11 --verify`: every file re-hashed, 0 mismatches.

## Private inputs
- `private_inputs/` added to .gitignore (commit 1, pushed before anything else).
- At session start **no file was present** in `private_inputs/` and nothing was attached:
  - BP96 / PM97 PDFs missing -> A5.1 (definitions from BP96) and R8 not done; A4 PM97 metric not extracted.
  - manuscript_draft.pdf missing -> C1 manuscript check not done (paper_tables.py is still built).

## Sanity check (2026-11-03)
`analysis10 --phase conf` rerun on the stored `results/conf/*.jsonl`: the regenerated `analysis.json` is
**identical** to the committed one (0 differing values; recursive comparison, tolerance 1e-9 relative). Headline:
P2 raw NSR 51.7 %, CHF 74.2 %, diff +22.4 % [12.8, 31.7]; masked diff -4.9 % [-13.7, +4.3]; P3 titration OR 2.31
[1.67, 3.19], CHF OR 0.77 [0.46, 1.28]; K3 8/988, U = 0.0146.

## Status
| Step | State | Commit |
|---|---|---|
| branch + .gitignore private_inputs/ pushed | done | 1st Phase 11 commit |
| env rebuilt, pytest 462 + 1 known failure | done | |
| data re-downloaded, SHA-256 OK | done | |
| sanity check analysis10 conf identical | done | |
| code: titration11 (+ exactness test, 20 series), hrt11 (+ tests), harness11 (+ tests), run11, data11 | done | |
| pilot (dev data nsrdb/chfdb, pilot seeds 11900+) for runtime budget: 0 errors; Tier 1 8.2 h + Tier 2 9.4 h = 17.6 h (no thinning) | done | |
| PREREGISTRATION.md committed + pushed (before any Phase 11 run on confirmatory data) | done | (this commit) |
| Tier 1 run started 2026-11-04 (`run_tier1.log`); prep, hrt, covar, sham_win, transplant done | running | 7514b6e |
| Tier 1 run: rdg_verify, rdg_synth, rdg_win done; sham_seg, rdg_seg running | running | |
| Tier 1 analysis: A1 windows, A2 (all), A3, A4, A5 windows/synthetic done (analysis11.json); A1 segments, A5 R1/R4 segments pending | partial | (this commit) |
| Tier 2 run + analysis | todo | |
| Tier 3 tables / figures | todo | |
| report docs/PHASE11_ROBUSTNESS_RESULTS.md | draft (Tier 1 windows filled; rest pending) | 5677e0a |
| exploratory dose-response (Amendment 1, `explore11.py`) | written; runs after Tier 2 | 1cf78a6 |
| paper_tables.py test run on partial JSON | ok (Table VI = Phase 10 §5.1) | |

## Resume commands
All runs are resumable (finished task ids skipped). Run from branch `phase11-robustness` (the guard checks
`origin/phase11-robustness`).
```bash
cd /home/user/arrhythmia-detection-pipeline
nohup env OMP_NUM_THREADS=1 PYTHONPATH=. /root/venv313/bin/python -m experiments.phase11_robustness.run11 --part tier1 --workers 4 >> experiments/phase11_robustness/results/conf/run_tier1.log 2>&1 &
nohup env OMP_NUM_THREADS=1 PYTHONPATH=. /root/venv313/bin/python -m experiments.phase11_robustness.run11 --part tier2 --workers 4 >> experiments/phase11_robustness/results/conf/run_tier2.log 2>&1 &
# optional, last: --part b3k1
PYTHONPATH=. /root/venv313/bin/python -m experiments.phase11_robustness.analysis11 --phase conf [--only A1,A2_hrt,...]
PYTHONPATH=. /root/venv313/bin/python -m experiments.phase11_robustness.paper_tables
PYTHONPATH=. /root/venv313/bin/python -m experiments.phase11_robustness.figures11
```
Progress: `tail experiments/phase11_robustness/results/conf/run_tier1.log`; line counts of `results/conf/*.jsonl`.

## Incidents
- 2026-11-04: container restart during Tier 1 (rdg_synth at 2,352/4,000). venv, data and results survived. Every
  JSONL checked for partial last lines (none found); the runner was restarted with the same command (resumable).

## Decisions
- D1 (branch): `phase11-robustness` as asked; mirrored to the session branch.
- D2 (pilot on development data): the runtime pilot runs on nsrdb / chfdb (same code paths) and pilot seeds
  11900-11901, so no confirmatory Phase 11 output exists before the preregistration is pushed.
- D3 (nstdb download): Phase 9 had no download script; `data11.py` does the same SHA-256-checked download.
- D4 (HRT definitions): filter rules confirmed from Cygankiewicz et al. 2003 (who state they follow the Schmidt
  group's published algorithm) and the Watanabe 2003 review; Bauer 2008 full text (JACC) returned 403; details in
  METHODS.md.

## Interim findings to carry into the report (not yet final)
- A5: verified readings R0, R3, R4, R6, R7 reproduce (i)-(iv) at window level -> no against-H result. Unverified:
  R1 (fails V2-V6), R2 (fails V3, V4, V6), R5 (fails V1, V4, V5). R5 (nested F) fails (ii) and REVERSES (iii) masked
  (-17.4 % [-31.7, -2.8]); report prominently even though R5 fails verification. R7 decisions identical to R0 on all
  windows: positivity is set by the base VWK test, the scheme changes only NL magnitude.

## Notes
- Never `pkill -f <pattern>` with a pattern that also appears in the calling shell command.
