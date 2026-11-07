# Phase 12 handoff (living document; update after every step)

**STATUS: COMPLETE (2026-11-07).** Report: `docs/PHASE12_ORIGINAL_METHODS_RESULTS.md`; numbers ledger:
`docs/PHASE12_NUMBERS.md`. Nothing is running.

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
| branch `phase12-original-methods`, papers in `private_inputs/` | done | 8493f90 |
| environment: venv present (Python 3.13.18; every pinned package matches the recipe); pytest 496 passed + 1 known failure (462 + 34 Phase 11) | done | f1cdee0 |
| data re-downloaded with the guarded commands; SHA-256 OK (nsr2db 54, chf2db 29, mitdb 48, nsrdb 18, chfdb 15, nstdb) | done | |
| analysis11 sanity rerun on stored results: analysis11.json byte-identical to the committed file | done | (this commit) |
| Tier 1A fidelity audit (FIDELITY.md): 4 decision-relevant deviations | done | 65e195f |
| variant module titration12 + exactness test (22 passed) | done | f1cdee0 |
| pilot (synthetic 12900+ and mitdb only): projected 0.19 h, no thinning | done | 8d681ec |
| PREREGISTRATION.md committed + pushed before any run on nsr2db/chf2db/nsrdb/chfdb | done | 8d681ec |
| Tier 1B run (b1, b2, b4a, b4v, b4s, b4t, b4w; 0 errors) + B3 Santa Fe A; analysis12 Tier 1B sections | done | (this commit) |
| Tier 1C PM97 replication (pm97, pm97nl; 0 errors; 33 subjects, none excluded) + analysis12 PM97 | done | (this commit) |
| Tier 2: paper_tables12 (Tables 12-1..12-11), figures12 (F12-1..3), numbers ledger (43 entries) | done | c53d999 |
| report docs/PHASE12_ORIGINAL_METHODS_RESULTS.md | **complete** | (this commit) |

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

- D7 (post hoc, labelled): PM97 'descriptive_handling' and 'exploratory_clean_segments' were added after the PM97
  outputs were seen; labelled in analysis12.json and in the report; nothing registered was changed.
- D8 (manuscript check): no manuscript draft was provided, so manuscript_check.md / manuscript_text_updates.md were not
  produced; the report's section 6 gives replacement sentences by content.
- D6 (B3): Santa Fe series A obtained from the public GitHub repository shubhomoydas/ad_examples (commit 4f20369,
  datasets/SantaFeCompetition/SF_A.dat, 1,000 points), sparse-cloned read-only to /home/user/shubhomoydas/ad_examples
  (outside this repository; data not committed). Series E: not found via PyPI / GitHub -> not obtainable.

## Findings to carry into the report
- Sanity: analysis11 rerun identical. R0 recomputed = Phase 10 stored decisions (988/988 raw, masked, edited) and
  Phase 11 transplant decisions (1756/1756).
- B1: R0 passes (d_opt 20/20; F_calc median 6,885 vs BP96 9,387; loss point median 35 % = lower edge of 35-65 %;
  BP96 ~50 %). R0-BP96 FAILS (loss 30 %). R0-BP96F passes (35 %).
- B2: 4/9 primary misses under every reading (ecological 34 vs 50; GOPY 42 vs 25; Rossler 25 vs 75; Lorenz 72 vs 50);
  Mackey-Glass tau = 100 never detected (0 vs 15: 'consistent' only by the +/-15 rule); MG tau = 17: 42 vs 15.
- B3: Santa Fe A R0 F_calc 10.51 (BP96 10.23), Z -22.1 (-24.6), F crit 1.17 (1.16).
- B4: all readings pass V1-V6; claims (i), (iii)-(v) hold under all; (ii) fails under R0-S-AAFT (+4.2 [-5.0, 13.8]).
  N5 rule -> 'no verified surrogate reading removes N5' (Phase 10 N5 100/100 under all readings; AAFT removes most
  AR(1) static-transform positives). Transplant V1-V0: R0 +64.2, BP96 +27.8, S-AAFT +29.8, S-FT +61.0.
  S2 isolated 5 %: R0 4/100 but R0-BP96F 39, R0-BP96 46.

- PM97 (nsrdb/chfdb; 30 midnight blocks, 3 lowest-mean-HR fallbacks): raw R0 NSR 52.6 % vs CHF 69.2 %
  (+16.5 [-1.9, 33.3]): PM97 pattern NOT reproduced (healthy not ~100 %, CHF not lower). R0-BP96 raw: -9.6 [-27, 7].
  AFTER ectopy handling CHF is LOWER: masked -45.6 [-57.0, -33.7], edited -26.4 [-42.8, -8.2], Wu -27.9 [-44.1, -10.3]
  -> registered 'against my claims' criterion (PM97 decrease survives handling) is MET on these databases.
  Within-subject ectopy OR 2.35 [1.78, 3.10]. Post hoc (labelled): ectopy-free segments NSR 53 % vs CHF 19 % (8 CHF
  subjects) -34 [-56, -6]. CHF masked-arm analysable fraction 34 % (NSR 98 %).

## Notes
- Never `pkill -f <pattern>` with a pattern that also appears in the calling shell command.
- Keep the session active while long jobs run (Phase 11 lost ~3 h to container restarts while idle).
