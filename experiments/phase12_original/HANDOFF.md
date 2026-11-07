# Phase 12 handoff (living document; update after every step)

Phase 12 = check the titration implementation against the original papers [BP96] and [PM97], test the
static-transform discrepancy and the omitted surrogate step, and replicate PM97 on nsrdb / chfdb.
Never merge into main, never open a pull request. Commit and push after every step.

**Branches.** Work on `phase12-original-methods` (created from `phase11-robustness` at `0f2b05d`). The session's
designated branch `claude/beautiful-hawking-i1ommp` is kept identical:
`git push origin HEAD:phase12-original-methods HEAD:claude/beautiful-hawking-i1ommp`. The runner guard checks the
CURRENT branch's origin copy, so run from `phase12-original-methods`.

## Inputs
- The user attached the two papers in the session (2026-11-06). Copied to `private_inputs/` (git-ignored, never
  committed): `barahona_poon_1996.pdf` (3-page scan, no text layer; read from 220-dpi renders, page halves and
  figure legends cropped), `poon_merrill_1997.pdf` (text layer).
- `manuscript_draft.pdf`: NOT provided -> the manuscript check (Tier 2.3) is skipped.
- Quoting limit in committed files: at most about 15 words of either paper.

## Status
| Step | State | Commit |
|---|---|---|
| branch `phase12-original-methods`, papers in `private_inputs/` | done | (this commit) |
| environment: venv present (Python 3.13.18; every pinned package matches the recipe); pytest 496 passed + 1 known failure (462 + 34 Phase 11) | done | f1cdee0 |
| data re-downloaded with the guarded commands; SHA-256 OK (nsr2db 54, chf2db 29, mitdb 48, nsrdb 18, chfdb 15, nstdb) | done | |
| analysis11 sanity rerun on stored results | running | |
| Tier 1A fidelity audit (FIDELITY.md): 4 decision-relevant deviations | done | 65e195f |
| variant module titration12 + exactness test (22 passed) | done | f1cdee0 |
| pilot (synthetic 12900+ and mitdb only): projected 0.19 h, no thinning | done | (this commit) |
| PREREGISTRATION.md committed + pushed before any run on nsr2db/chf2db/nsrdb/chfdb | done | (this commit) |
| Tier 1B B1-B4 | todo | |
| Tier 1C PM97 replication | todo | |
| Tier 2 tables, figures, numbers ledger | todo | |
| report docs/PHASE12_ORIGINAL_METHODS_RESULTS.md | todo | |

## Resume commands
Run from branch `phase12-original-methods` (the guard checks `origin/phase12-original-methods`); resumable.
```bash
cd /home/user/arrhythmia-detection-pipeline
nohup env OMP_NUM_THREADS=1 PYTHONPATH=. /root/venv313/bin/python -m experiments.phase12_original.run12 --part tier1b --workers 4 >> experiments/phase12_original/results/conf/run_tier1b.log 2>&1 &
nohup env OMP_NUM_THREADS=1 PYTHONPATH=. /root/venv313/bin/python -m experiments.phase12_original.run12 --part tier1c --workers 4 >> experiments/phase12_original/results/conf/run_tier1c.log 2>&1 &
PYTHONPATH=. /root/venv313/bin/python -m experiments.phase12_original.analysis12
```

## Decisions
- D1 (branch): `phase12-original-methods` as asked; mirrored to the session branch.
- D2 (environment): the venv built in Phase 11 from `requirements-phase10.txt` was still present; every pinned version was
  checked against `pip freeze` (all match) instead of reinstalling the same versions.
- D3 (readings): R0-BP96 = plain F_calc + N = series length + chi-square gate -> Mann-Whitney (FIDELITY F7, F12, F16);
  R0-BP96F = without the gate; R0-S-FT / R0-S-AAFT = R0 + BP96 surrogate step (F18). R0 stays primary.
- D4 (pilot): synthetic pilot seeds and mitdb only, so nothing from nsr2db/chf2db/nsrdb/chfdb is touched before the
  preregistration is pushed.
- D5 (Mackey-Glass tau = 100): the BP96 scan describes its Mackey-Glass series as having an implicit dimension of seven
  (tau ~ 100 in Farmer 1982), so tau = 100 is primary and tau = 17 secondary. A one-seed smoke test (disclosed in the
  preregistration) showed tau = 100 undetected noise-free; parameters were not changed.

## Notes
- Never `pkill -f <pattern>` with a pattern that also appears in the calling shell command.
- Keep the session active while long jobs run (Phase 11 lost ~3 h to container restarts while idle).
