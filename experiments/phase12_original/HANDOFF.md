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
| environment check, pytest, data re-download + SHA-256, analysis11 sanity check | todo | |
| Tier 1A fidelity audit (FIDELITY.md) | todo | |
| variant module + exactness test | todo | |
| pilot (synthetic + mitdb only) and PREREGISTRATION.md pushed | todo | |
| Tier 1B B1-B4 | todo | |
| Tier 1C PM97 replication | todo | |
| Tier 2 tables, figures, numbers ledger | todo | |
| report docs/PHASE12_ORIGINAL_METHODS_RESULTS.md | todo | |

## Resume commands
(filled in as runners are written)

## Decisions
- D1 (branch): `phase12-original-methods` as asked; mirrored to the session branch.

## Notes
- Never `pkill -f <pattern>` with a pattern that also appears in the calling shell command.
- Keep the session active while long jobs run (Phase 11 lost ~3 h to container restarts while idle).
