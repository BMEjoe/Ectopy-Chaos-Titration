# Phase 11: registered robustness analyses for "Ectopic Beats Masquerade as Chaos in Heart Rate"

**Status: DRAFT, in progress.** Tier 1 results are final. Still pending:
Tier 2 (B1-B4), the exploratory dose-response, and Tier 3 (tables, figures). Every number below is read from
`experiments/phase11_robustness/results/conf/analysis/analysis11.json`; nothing is hand-entered.

Branch `phase11-robustness` (mirrored to `claude/beautiful-hawking-i1ommp`). The preregistration
(`experiments/phase11_robustness/PREREGISTRATION.md`) was pushed in commit `2b0e459` before any Phase 11 analysis
touched confirmatory data. One dated amendment, `1cf78a6`, adds an exploratory analysis and changes nothing
registered. Phase 10 code and thresholds are frozen; `final_pipeline.py` is unchanged. CIs are 95% subject-cluster
bootstrap CIs (10,000 resamples) unless stated otherwise.

## Summary

| Analysis | Prediction under H | Result (95 % CI) | Verdict |
|---|---|---|---|
| A1 sham control, windows, masked | D > 0, Q < 0.5 | D = 0.356 [0.248, 0.469]; Q = 0.36 [0.26, 0.47] | supports, but a third of the masked-arm drop is data loss |
| A1 sham control, windows, edited | D > 0, Q < 0.5 | D = 0.672 [0.582, 0.752]; Q = 0.01 [0.00, 0.02] | supports (weak control, see text) |
| A1 sham control, segments, masked / edited / Wu | D > 0, Q < 0.5 | D = 0.554 [0.488, 0.615], Q = 0.20 / D = 0.591 [0.518, 0.658], Q = 0.007 / D = 0.602 [0.531, 0.668], Q = 0.007 | supports (all five cells: **A1 supports**) |
| A2 transplant: one ventricular-type beat into 439 clean real windows | positive rate rises | 22.1 % to 86.3 %, Δ = +64.2 pp [59.2, 69.3] | **supports** |
| A2 heart-rate turbulence as mechanism | (not predicted) | HRT alone Δ = −0.7 pp [−1.7, +0.2]; HRT added to the beat changes nothing | HRT is not the mechanism in real windows |
| A2 synthetic single beat (N1 linear RR) | rises | 0/100 to 89/100 | supports; contrasts with Phase 10's 5/100 at 5 % density |
| A2 extended masking (20 HRT intervals) | removes positives (and/or) | standard − extended = −1.5 pp [−4.7, +1.4] | nothing removed |
| A3 within-subject burden OR (conditional logit + covariates, windows) | OR > 1 | 3.40 [2.62, 4.41] (bootstrap [2.35, 5.97]) | **supports** |
| A3 CHF effect adjusted for burden + covariates | (not predicted) | OR 0.88 [0.46, 1.69] | no group effect |
| A3 segments, within-subject OR | OR > 1 | 2.66 [2.47, 2.87] | supports |
| A3 K1 (secondary), within-subject OR | (secondary) | 1.42 [0.85, 2.38] (16 positives) | **not significant** |
| A4 original-study metric (≥ 1 positive segment per subject) | ≈ 100 % in both groups on raw RR | NSR 54/54, CHF 29/29 | reproduced but uninformative |
| A5 verified readings (R0, R3, R4, R6, R7) reproduce claims (i)-(iv), windows | all reproduce | all reproduce | **supports** |
| A5 unverified readings (R1, R2, R5) | (flagged) | R1, R2 reproduce; **R5 fails (ii) and reverses (iii)** | reported, see A5 |
| A5 segments (R0, R1, R4) | reproduce (ii), (iii) | raw diff +22.4 / +20.6 / +28.3 pp (all CIs > 0); masked and Wu CIs include 0 | supports (R4 Wu +6.7 [−1.0, +14.5]) |
| B1 seeds (titration salts 1, 2; K3 and K1 seeds) | small ranges, no conclusion changes | 0 of 16 headline conclusions change; titration decisions 0 flips (by construction, see text); K3 8 → 11 / 10 detections (6 stable), K1 16 → 19 / 16 (10-11 stable) | supports (registered rule); **K3 / K1 window-level detections are seed-sensitive** |
| B2 K3 at 1/128 s, false positives | ≤ 7/100 in every FP condition | **3 of 155 conditions fail** (KTz P = 120 / 200 / 250, no ectopy: 14, 10, 10 per 100); pooled 55/10,092 vs 1/10,092 at 1/360 s | **weakens K3 claims** |
| B2 degeneracy of the 8 real K3 detections | only nsr029 w9 degenerate | only nsr029 w9 degenerate (surrogate SD 0); other 7 SD 0.04-0.08 | supports |
| B3 variance-preserving spike-ins | (descriptive bounds) | f_min unchanged in 15 of 16 rows at π < 0.05 and 14 of 16 at π < 0.20. Logistic r = 3.88 is weaker at π < 0.05 (0.5 to 0.7) but stronger at π < 0.20 (0.5 to 0.3); vdP (8, 3.3) is stronger at π < 0.20 (0.3 to 0.2) | bounds robust to the design |
| B4 second base window | (descriptive bounds) | *pending* | *pending* |
| C1 manuscript-number checks | | "7 of 8" run claim unsupported (4 windows); other items, see C1 | correct the manuscript |

## Deviations, missing inputs, incidents

- **Private inputs were absent.** `private_inputs/` was empty at session start.
  - No Barahona & Poon 1996 PDF: A5.1 (definitions from BP96) was not done and reading R8 was not run.
  - No Poon & Merrill 1997 PDF: their own metric was not extracted for A4; only the generic subject-level metrics were
    computed.
  - No manuscript draft: the C1 manuscript check was not done. `paper_tables.py` was still built, and the
    known-issue list was checked against the code and the Phase 10 outputs.
- **Clock time is unavailable.** The headers carry no start time, so A3 uses elapsed hours since the first beat.
- **Environment.** Python 3.13.18 (Phase 10: 3.13.14) and wrapt 2.5.2 (2.5.0); everything else matches the
  Phase 10 recipe. pytest gives 462 passed and 1 known failure, the same as Phase 10.
- **Sanity check.** Rerunning `analysis10 --phase conf` on the stored Phase 10 results reproduced `analysis.json`
  exactly (0 differing values).
- **Incident.** A container restart killed the Tier 1 runner during `rdg_synth`. No JSONL file had a partial line,
  and the resumable runner was restarted with the same command.
- **Amendment 1** (2026-11-04, exploratory only; written before running): an ectopy dose-response on N1. It is
  reported in its own section and has no verdict.

---

## A1. Sham-handling control

**Question.** Does masking or editing remove titration positives only because it throws data away?

**Design.**
- Each ectopic window got up to 5 sham masks: circular shifts of its own K3 mask that avoid every ectopy-related
  interval and the 20 heart-rate-turbulence (HRT) intervals after each abnormal beat.
- The primary statistic, per unit, compares the true arm with its shams (D = F_true − F_sham), with Q = F_sham / F_true.
- Units are windows that are raw positive, usable in the true arm, and have at least one usable sham.

**Windows, masked arm (pooled; 126 units, 56 subjects).**
- The true mask made 55.6 % of raw positives negative; the sham masks made 19.9 % negative.
- D = 0.356 [0.248, 0.469]; Q = 0.36 [0.26, 0.47].
- NSR: D = 0.331 [0.207, 0.465], Q = 0.39 [0.27, **0.53**].
- CHF: D = 0.436 [0.259, 0.640], Q = 0.27 [0.14, 0.40].
- **Verdict: supports.**

**Windows, edited arm (pooled; 258 units, 68 subjects).** D = 0.672 [0.582, 0.752]; Q = 0.009 [0.001, 0.021].
Supports.

**What weakens or limits this.**
1. **About a third of the masked-arm effect is data loss.** Removing the same number of regression rows at random
   positions, with the ectopy left in, turns 20 % of raw-positive windows negative. The masked-arm drop is therefore
   not all "ectopy removed". In NSR the upper confidence limit of Q (0.53) passes the registered 0.5 line. The point
   estimate, which is what the registered criterion uses, does not.
2. **The edited-arm sham is a weak control by construction.** Interpolating normal intervals changes a series very
   little, so F_sham ≈ 0 is expected. Its large D shows that editing *normal* beats does not remove positives. It does
   not show that the edited arm's drop is free of information loss in the way the masked-arm sham does.
3. **High-burden windows are excluded.** In 121 ectopic windows no sham could avoid the forbidden zone, and 132 more
   had fewer than 5 shams. A1 therefore speaks for low- to moderate-burden windows.
4. **Secondary (clean windows with donor masks).** Masking alone made 9.6 % [5.6, 14.3] of raw-positive clean windows
   negative, and editing 7.1 % [4.5, 9.8]. Among raw-negative ectopic windows, the true mask created positives in
   3.5 % and the sham in 4.5 %.

**Segments** (12-min segments, 2 shams each; pooled):

| Arm | Units (subjects) | F_true | F_sham | D [CI] | Q [CI] | Verdict |
|---|---|---|---|---|---|---|
| masked | 3,124 (78) | 0.689 | 0.135 | 0.554 [0.488, 0.615] | 0.20 [0.17, 0.23] | supports |
| edited | 3,308 (80) | 0.595 | 0.004 | 0.591 [0.518, 0.658] | 0.007 [0.004, 0.011] | supports |
| Wu | 3,308 (80) | 0.606 | 0.004 | 0.602 [0.531, 0.668] | 0.007 [0.004, 0.012] | supports |

- By group, masked: NSR D = 0.531 [0.446, 0.608], Q = 0.20; CHF D = 0.598 [0.496, 0.692], Q = 0.19.
- Edited and Wu cells look the same in NSR and CHF.
- **A1 summary verdict: supports.** All five pooled cells support.

Limits at segment level:
- **High-burden segments are excluded.** 1,349 segments had no admissible sham, and 1,372 had fewer than 2.
- **Masking alone removes positives.** Donor masks made 17.8 % [15.8, 20.0] of raw-positive clean segments negative
  in the masked arm (edited 4.6 %, Wu 7.3 %).
- **Wu-style deletion of normal intervals is nearly inert** (Q ≈ 0.007), like editing. Only the masked-arm sham
  removes enough regression rows to probe data loss. Each masked interval removes every row whose lag span contains
  it, whereas deleting one interval removes one interval.

## A2. The single-beat result: a beat, not its turbulence

**HRT in the confirmatory data (descriptive).** Subjects with at least 5 usable VPBs:

| Group | n | TO (%) median [IQR] | TS (ms/RR) median [IQR] | Abnormal TO / TS |
|---|---|---|---|---|
| NSR | 18 | −1.55 [−2.96, −0.91] | 6.71 [3.84, 9.81] | 0 % / 6 % |
| CHF | 24 | +0.52 [+0.06, +1.12] | 1.26 [0.46, 1.77] | 75 % / 88 % |

- CHF turbulence is blunted, as expected.
- Usable VPBs: NSR 2,111 of 8,463 V beats; CHF 4,450 of 86,096.
- Six details of the HRT definition could not be checked against a primary source; they are listed in METHODS.md.

**Transplant into 439 clean real windows** (69 subjects). Each window gets one event at a seeded position. The draws
are paired and the titration seed is shared, so the variants differ only in the inserted event. Pooled rates:

| Variant | Positive | Δ vs V0 (pp) [CI] | McNemar (only-V0 / only-variant) |
|---|---|---|---|
| V0 untouched | 22.1 % | – | – |
| V1 one ventricular-type beat, no HRT | **86.3 %** | **+64.2 [59.2, 69.3]** | 2 / 284 |
| V2 V1 + group's empirical HRT | 86.3 % | +64.2 [59.2, 69.3] | 2 / 284 |
| V3 V1 + normal HRT (TO −2.5 %, TS 13) | 86.8 % | +64.7 [59.5, 69.9] | 2 / 286 |
| V4 V1 + blunted HRT (TO 0 %, TS 2.5) | 86.8 % | +64.7 [59.7, 69.8] | 2 / 286 |
| V5 one atrial-type beat | 53.8 % | +31.7 [26.4, 37.0] | 5 / 144 |
| V6 one couplet | 90.7 % | +68.6 [63.4, 73.7] | 0 / 301 |
| V7 empirical HRT, no premature beat | 21.4 % | −0.7 [−1.7, +0.2] | 4 / 1 |

By group, V1: NSR +65.5 [60.3, 70.9]; CHF +58.8 [44.6, 73.3].

**Reading.**
- A single ventricular-type beat in 512 intervals is enough to make titration fire in most clean real windows.
- Turbulence is not the mechanism: HRT alone does nothing, and adding normal or blunted HRT to the beat changes
  nothing.
- The couplet gives the largest rise, as predicted, though its CI overlaps V1's.
- **Verdict: supports H.**

**Two discrepancies, reported as found.**
1. **Transplanted and real beats differ.** Real burden-1 windows are positive at 47 % [34, 60] with a single real
   V beat (25/53) and 76 % [64, 85] with a single real A beat (48/63). The transplanted ventricular-type beat (86 %)
   overshoots the real V rate, and the transplanted atrial-type beat (54 %) undershoots the real A rate. So the beat
   model reproduces the direction of the 22 % to 63 % jump but not its size by beat type. Plausible causes include
   prematurity, annotation quality and edge position, but none was tested; the cause is open.
2. **HRT does matter in some synthetic regimes.** The beat model on synthetic series (100 windows each):

   | Condition | V0 | V1 | V2 | V5 | V6 | V7 (HRT only) |
   |---|---|---|---|---|---|---|
   | N1 linear RR | 0 | **89** | 89 | 89 | 89 | 0 |
   | vdP regime 26 | 100 | 100 | 100 | 100 | 100 | 100 |
   | vdP 27 | 0 | 93 | 93 | 93 | 93 | 3 |
   | vdP 28 | 0 | 26 | 58 | 17 | 91 | **89** |
   | vdP 29 | 0 | 44 | 70 | 20 | 92 | **85** |
   | vdP 30 | 0 | 0 | 0 | 0 | 0 | 0 |
   | vdP 31 | 0 | 1 | 56 | 3 | 71 | **92** |

   In three low-variability vdP regimes, the HRT profile alone makes titration fire, whereas in real windows it does
   not. HRT can act as a trigger when intrinsic variability is small; it is not what drives the real-data effect.
   (vdP regime 26 is the forced quasi-periodic regime that Phase 10 already found positive without ectopy.)

**The synthetic single-beat result conflicts with Phase 10.** One isolated beat in N1 gives 89/100 positives
(Δ = +89 [82, 95]). Phase 10 reported 5/100 and 0/100 for the same beat model at 5 % and 10 % density (S2). Titration
therefore reacts to a single high-leverage event more than to many. The exploratory dose-response (below) maps the
curve between these points.

**Extended masking.** In 134 burden 1-4 windows analysable under both masks, masking the 20 post-beat intervals as well
changed the positive rate by standard − extended = −1.5 pp [−4.7, +1.4] (McNemar 1 / 3). It removed no extra
positives, consistent with the finding that HRT is not the mechanism.

## A3. Confounding

**Covariates.** All from NN intervals only: mean NN, log SDNN, log RMSSD, half-shift, log(1 + detrend ratio), and
elapsed-hour bins. One window had an undefined covariate and was dropped (987 windows used).

**Windows, titration raw (exposure x = log2(1 + burden)):**

| Model | Burden OR [CI] | CHF OR [CI] |
|---|---|---|
| GEE, x + CHF (Phase 10 model) | 2.31 [1.67, 3.19] | 0.77 [0.46, 1.28] |
| GEE, x + CHF + covariates | 3.17 [2.04, 4.91] | **0.88 [0.46, 1.69]** |
| **Conditional logit by subject, x + covariates (PRIMARY)** | **3.40 [2.62, 4.41]**; bootstrap [2.35, 5.97] | – |
| Conditional logit, x only | 3.29 [2.55, 4.24] | – |
| Within NSR, conditional logit + covariates | 8.35 [5.15, 13.55] | – |
| Within CHF, conditional logit + covariates | 2.17 [1.63, 2.88] | – |

- 77 of 83 subjects (915 windows) had both outcomes and inform the conditional logit.
- Covariates do not explain the burden effect; adjustment strengthens it.
- **Verdict: supports.**
- The burden effect is much steeper within NSR than within CHF. This fits CHF windows already being positive more
  often at low burden.

**Segments** (8,987 segments, 6 dropped): conditional logit + covariates OR 2.66 [2.47, 2.87]; GEE + covariates
2.06 [1.53, 2.78]; CHF effect 0.95 [0.61, 1.49].

**Secondary outcomes (windows).**
- LLE alone: within-subject OR 1.69 [1.41, 2.04].
- **K1 (LLE AND UPO): within-subject OR 1.42 [0.85, 2.38], not significant.** Only 16 windows are positive, the NSR
  models are not estimable (4 positives), and the GEE OR (1.48 [1.10, 1.99]) does not survive the within-subject
  design. The Phase 10 statement that K1's positives "track burden" is weaker than stated (see the manuscript section).

## A4. The original studies' own metrics (descriptive)

Fraction of subjects with at least one positive 12-min segment:

| Arm | NSR | CHF |
|---|---|---|
| raw | 54/54 (100 %) | 29/29 (100 %) |
| masked | 53/54 (98 %) | 22/27 (81 %) |
| edited | 54/54 (100 %) | 28/29 (97 %) |
| Wu | 54/54 (100 %) | 29/29 (100 %) |

- Subjects with at least 50 % of segments positive, raw: NSR 50 % [37, 63]; CHF 86 % [69, 95]. After Wu
  elimination: NSR 17 %, CHF 14 %.
- Wu-style detection rate per 3-h block, raw: NSR 46-60 % and CHF 71-78 %. Masked: NSR 17-35 % and CHF 19-27 %.
- Mean NL among positives, raw: NSR 0.37, CHF 0.59. Wu: NSR 0.19, CHF 0.33.
- As predicted, "chaos in every subject" is reproducible here and uninformative: it holds in both groups, raw and
  after Wu elimination.
- On raw RR, the direction is opposite to "decreased cardiac chaos in CHF": CHF has more positive segments. Poon &
  Merrill's exact metric could not be checked (no PDF).

## A5. Sensitivity to the implementation readings

**Verification (V1-V6, Phase 10 criteria):**

| Reading | Passes | Fails |
|---|---|---|
| R0 primary (frozen) | V1-V6 | – |
| R1 linear memory ≤ 6 | V1 | V2-V6: positive on limit cycle, torus, period-3 logistic and periodic Mackey–Glass (20/20 each) |
| R2 linear memory ≤ 20 | V1, V2, V5 | V3, V4 (8/18 periodic logistic NL = 0), V6 (periodic MG NL > 0) |
| R3 linear memory ≤ 50 | V1-V6 | – |
| R4 Mann–Whitney test | V1-V6 | – |
| R5 nested partial F-test | V2, V3, V6 | V1 (NL(3.7) = 1.36), V4, V5 (period-3 logistic positive 12/20) |
| R6 ε as a variance | V1-V6 | – |
| R7 noise increments of 1 % | V1-V6 | – |

**Claims matrix, windows (988).** Claims: (i) burden OR > 1; (ii) raw CHF − NSR > 0; (iii) post-handling CHF − NSR CI
includes 0, for both masked and edited; (iv) couplets positive (≥ 50/100) while N1 is not (≤ 5/100).

| Reading | Verified | (i) OR [CI] | (ii) raw diff [CI] | (iii) masked diff [CI] | (iv) E3 / N1 | Raw rate |
|---|---|---|---|---|---|---|
| R0 | yes | 2.31 [1.67, 3.19] ✓ | +23.5 [12.6, 34.3] ✓ | +1.1 [−10.5, +14.3] ✓ | 100 / 0 ✓ | 46.1 % |
| R1 | **no** | 2.88 ✓ | +26.3 [15.3, 37.4] ✓ | −1.3 [−11.3, +9.0] ✓ | ✓ | 54.9 % |
| R2 | **no** | 2.44 ✓ | +24.2 [13.0, 35.1] ✓ | −5.0 [−14.0, +4.1] ✓ | ✓ | 50.2 % |
| R3 | yes | 2.33 ✓ | +25.1 [14.2, 35.9] ✓ | −2.1 [−12.8, +10.0] ✓ | ✓ | 47.6 % |
| R4 | yes | 1.92 [1.58, 2.33] ✓ | +20.3 [9.6, 31.1] ✓ | −0.6 [−12.1, +12.1] ✓ | ✓ | 39.2 % |
| R5 | **no** | 3.27 [1.85, 5.78] ✓ | **+5.6 [−2.3, +13.4] ✗** | **−17.4 [−31.7, −2.8] reversed** | ✓ | 81.4 % |
| R6 | yes | 3.22 [2.19, 4.74] ✓ | +24.6 [13.8, 35.5] ✓ | +4.8 [−8.4, +19.5] ✓ | ✓ | 47.1 % |
| R7 | yes | = R0 | = R0 | = R0 | = R0 | 46.1 % |

The edited-arm CHF − NSR CI includes 0 in every reading. R5 is the closest to excluding it: −11.0 [−22.6, +0.4]. Full numbers are in `analysis11.json` → `A5.matrix`.

**Verdict.** Every reading that passes V1-V6 reproduces (i)-(iv) at window level, so there is **no against-H result**.

**What weakens it (reported with the same weight).**
- **R5, the nested partial F-test, changes the conclusions.**
  - Titration becomes positive in 81 % of raw windows.
  - The raw group difference disappears (+5.6 pp, CI includes 0).
  - After masking, CHF has *fewer* positives than NSR (−17.4 pp [−31.7, −2.8]).
  - R5 fails verification: NL at logistic r = 3.7 is twice the published value, and it is positive on a period-3
    orbit. The R5 result therefore shows the claims hold only for implementations that behave like the published
    method on the standard systems; an implementation choice that inflates sensitivity can change the group
    comparison.
- **R1 and R2 (short linear memory) fail verification** because they read periodic signals as nonlinear. They
  still reproduce every claim on the heart-rate data.
- **R7 is identical to R0 on every window.** Positivity is decided by the noise-free VWK test, and the noise schedule
  only changes NL's magnitude. R7 is therefore not an independent check of the binary claims.

**Segments** (CHF − NSR, difference of subject mean rates; R0 = Phase 10 stored results):

| Reading | (ii) raw [CI] | (iii) masked [CI] | (iii) Wu [CI] | Holds |
|---|---|---|---|---|
| R0 | +22.4 [+12.8, +31.7] | −4.9 [−13.7, +4.3] | +2.9 [−5.3, +11.4] | yes |
| R1 (unverified) | +20.6 [+12.1, +28.7] | +0.4 [−9.6, +10.5] | +1.0 [−8.6, +10.8] | yes |
| R4 Mann–Whitney | +28.3 [+17.8, +38.9] | −0.0 [−7.5, +7.6] | +6.7 [−1.0, +14.5] | yes |

- Claims (ii) and (iii) hold at segment level for all three readings.
- R4 strengthens the raw difference: NSR 36 % vs CHF 64 % of segments positive.
- Under R4, the post-Wu difference is the one that comes closest to excluding 0 (p = 0.08).

---

## B2. K3 at 1/128 s (Phase 9 TEST groups, new seeds)

**Data.** 12,290 windows, re-quantized from 1/360 s to 1/128 s (beat times rounded).
- 56 windows failed to generate ("fewer than n intervals", fast coupled-vdP regimes), the same generator limit
  Phase 9 TEST hit 63 times. They are left out.
- 12,234 windows were analysed: 10,092 false-positive (FP) windows in 155 conditions and 2,142 power windows.

**False positives (registered pass rule: ≤ 7 detections per 100 windows in every FP condition).**

| Resolution | Pooled FP | Conditions > 7/100 | Pass |
|---|---|---|---|
| 1/360 s (native) | 1 / 10,092 | none | yes |
| **1/128 s (the resolution of the real data)** | **55 / 10,092** | **KTz P = 120, 200, 250 (no ectopy, real noise channel): 14, 10, 10 per 100** | **no** |

**Verdict: weakens the K3 claims.** At the resolution of nsr2db / chf2db, K3 is not size-controlled in every
non-chaotic condition.

**Where the false positives come from (description, not a re-definition of the rule).**
- **They are confined to KTz.** 53 of the 55 are KTz windows, whose RR is the pacing input plus the real noise
  channel. Their NN SDNN is about 4.7 ms (median), less than one quantization step at 1/128 s (7.8 ms).
- **43 of the 55 are degenerate.** The surrogate prediction-error SD at the m attaining z_max is < 1e-6, and z_max
  reaches ~10^15. With intervals taking one or two values, the IAAFT surrogates are near-identical and the z-score
  divides by about zero.
- **Outside KTz there are 2 false positives in 6,992 windows** at 1/128 s, against 1 at 1/360 s.
- 2,182 FP windows are not analysable at 1/128 s. They count as not detected, as in Phase 9.

**Real detections (degeneracy diagnostic).** Of the 8 Phase 10 real K3 detections, only nsr029 window 9 is degenerate
(SDNN 15 ms, NN MAD 7.8 ms, surrogate SD 0, z_max 2 × 10^11), as Phase 10 reported. The other seven have SDNN 39-202 ms,
surrogate SD 0.039-0.080 and z_max 10.0-16.1, so the registered "weakens" criterion for real detections is not met.
Real windows with very low variability are exposed to the same artefact; nsr029 w9 is the one example.

**Power** (detections / windows; 1/128 s vs 1/360 s):

| Group | Family | 1/128 s | 1/360 s |
|---|---|---|---|
| chaos_primary | coupled vdP | 66/312 | 66/312 |
| chaos_primary | phase-reset | 39/630 | 25/630 |
| chaos_primary | KTz | 23/330 | 0/330 |
| chaos_second | coupled vdP | 40/120 | 34/120 |
| chaos_second | Mackey–Glass | 67/210 | 82/210 |
| chaos_second | phase-reset | 16/210 | 15/210 |
| chaos_second | KTz | 0/330 | 0/330 |
| **pooled** | primary / second | **10.1 % / 14.1 %** | **7.2 % / 15.1 %** |

The apparent power gain at 1/128 s for chaos_primary is partly the artefact above. The 23 KTz "detections" are on
series whose RR carries no chaos (Phase 9: KTz chaos is in the T wave only), and 22 of the 251 power detections at
1/128 s are degenerate.

## B3. Variance-preserving spike-ins (K3)

**Design.** x'_k = μ_NN + √(1 − f)(x_k − μ_NN) + √f σ_NN c_k, applied to every interval. The base windows,
chaotic components, f grid and 1/128 s quantization are the same as Phase 10's additive design.
- 1,328 tasks and 9,296 series, with 0 errors and 0 out-of-range intervals.
- Bounds use Phase 10's real K3 bound U = 0.0146, as registered.
- `primary_phase10_recomputed` reproduces Phase 10 Table VI exactly.

Detected / 83 at f = 0.05 … 0.9, and f_min (design formula) at π < 0.05 / 0.20:

| Family, parameter (λ) | Additive (Phase 10) | f_min | Variance-preserving | f_min | Wilson-lower f_min (0.20) |
|---|---|---|---|---|---|
| Hénon a = 1.22 (0.303) | 0 0 3 6 32 58 76 | 0.5 / 0.5 | 0 1 3 6 26 53 66 | 0.5 / 0.5 | 0.5 |
| Hénon a = 1.40 (0.419) | 0 1 4 9 29 58 76 | 0.5 / 0.3 | 0 0 3 8 27 53 71 | 0.5 / 0.3 | 0.5 |
| Hénon a = 1.14 (0.244) | 0 0 0 2 2 9 8 | none / 0.7 | 0 1 1 1 3 10 13 | none / 0.7 | 0.9 |
| logistic r = 3.88 (0.464) | 1 1 1 5 28 59 77 | 0.5 / 0.5 | 0 1 3 7 24 52 65 | **0.7** / **0.3** | 0.5 |
| logistic r = 4.00 (0.693) | 0 0 0 1 9 25 73 | 0.7 / 0.5 | 0 0 0 2 7 26 62 | 0.7 / 0.5 | 0.7 |
| coupled vdP (8, 3.3) (0.075) | 1 2 5 13 33 48 71 | 0.5 / 0.3 | 0 2 7 11 33 48 60 | 0.5 / **0.2** | 0.3 |
| coupled vdP (6, 3.3) (0.079) | 0 0 0 3 10 13 23 | none / 0.5 | 0 1 0 1 8 11 18 | none / 0.5 | 0.7 |
| coupled vdP (2, 2.7) (0.085) | 0 0 0 2 6 20 63 | 0.9 / 0.7 | 0 0 0 0 3 21 47 | 0.9 / 0.7 | 0.7 |
| coupled vdP (6, 4.0) (0.101) | 0 1 2 6 13 15 27 | 0.9 / 0.5 | 1 0 1 1 12 15 27 | 0.9 / 0.5 | 0.5 |
| weak maps (Hénon 1.08, logistic 3.58 / 3.65), all phase-reset | ≤ 5 at any f | none / none | ≤ 3 at any f | none / none | none |

**Reading (descriptive; no registered prediction).**
- The exclusion bounds do not depend on the spike-in design.
- Of 16 rows, the π < 0.05 f_min changes in one (logistic r = 3.88, 0.5 to 0.7, weaker).
- The π < 0.20 f_min changes in two: logistic r = 3.88 from 0.5 to 0.3, and vdP (8, 3.3) from 0.3 to 0.2, both
  stronger.
- Detection at f = 0.9 is lower under variance preservation (e.g. logistic r = 4.00: 62 vs 73), because the base
  window's own variability is shrunk rather than kept.
- The Wilson-lower variant is more conservative, typically one f step higher.

## B1. Seeds

**Registered rule.** A headline conclusion changes if its CI moves to the other side of its null. There are 16
checks: P3 burden OR and CHF OR in windows; P2 CHF − NSR raw, masked and Wu in segments; P2 masked change by group;
and LLE burden OR. Each is checked at two new seeds.
**Result: 0 of 16 change. Verdict: supports.**

**Titration (salts 1 and 2): no decision changes, by construction.**
- 0 of 988 windows flip (raw, masked, edited) and 0 of 8,993 segments flip (raw, masked, Wu). Every titration headline
  number is identical to Phase 10's.
- **The new salts were applied.** NL values differ in 453 of 455 raw-positive windows (mean |ΔNL| 0.035, max 0.31).
- **The decision cannot change.** Titration declares a series positive when the noise-free VWK test is significant.
  The noise realisations, which are the only thing the seed controls, change only the size of NL. (R7 gave the same
  result for the noise schedule.)
- So B1 cannot test the seed robustness of titration's yes/no claims. That robustness is guaranteed, and it also
  means "titration-positive" in this paper is the same as "the noise-free VWK test favours a nonlinear model". The
  manuscript should say so.

**K3 (masked_growth_seed_entropy 11101, 11102): the count is stable in size; the set of windows is not.**

| Seed | Detected / 988 | Flips vs Phase 10 | Phase 10 detections recovered | New windows |
|---|---|---|---|---|
| Phase 10 | 8 | – | – | – |
| 11101 | 11 | 5 | 7 of 8 (not nsr029 w9) | chf216 w7, w12; nsr035 w5; nsr048 w3 |
| 11102 | 10 | 6 | 6 of 8 (not nsr029 w9, chf216 w2) | chf225 w4; nsr014 w11; nsr048 w3, w12 |

- Six windows are detected under all three seeds: chf211 w7, chf218 w12, nsr049 w3/w5/w6, nsr052 w6. Fifteen windows
  are detected under at least one seed.
- The degenerate nsr029 w9 is not reproduced under either new seed.
- **Descriptive, not registered** (`B1_K3_bounds`). With 11 or 10 detections, Phase 10's real bound U rises from
  0.0146 to 0.0191 / 0.0171, and the primary-design f_min at π < 0.05 gets weaker:
  - logistic r = 4.00: 0.7 to 0.9;
  - logistic r = 3.88: 0.5 to 0.7;
  - Hénon a = 1.40: 0.5 to 0.7 (seed 11101 only);
  - coupled vdP (6, 4.0): 0.9 to no bound.
  The π < 0.20 bounds are unchanged.

**K1 and its components (seed pairs (11201, 11211), (11202, 11212)):**

| | Phase 10 | Seed pair 1 | Seed pair 2 | Flips (1 / 2) |
|---|---|---|---|---|
| K1 (LLE AND UPO) | 16 | 19 | 16 | 13 / 12 (about 10-11 windows stable) |
| LLE | 401 | 393 | 403 | 38 / 36 |
| UPO | 31 | 34 | 28 | 19 / 25 |
| LLE burden OR (GEE) | 1.51 [1.31, 1.75] | 1.51 [1.30, 1.76] | 1.58 [1.35, 1.85] | – |

The K1 count is stable, but about a third of its detections, and most of the UPO detections, depend on the seed.

## B4 (Tier 2): pending

## Exploratory: ectopy dose-response on N1 (Amendment 1): pending

## C1. Manuscript numbers (paper_tables.py)

`paper_tables.py` builds Tables I-VII and the supplementary tables from the analysis JSON only (*to be run*). There
is no manuscript draft, so the known issues were checked against code and data:

- **"7 of 8 windows" with an atrial run of three: unsupported.**
  - Only 4 confirmatory windows contain a single annotated supraventricular run of length ≥ 3: nsr010 w4 (5),
    nsr036 w7 (9), nsr046 w7 (3) and chf206 w4 (513, a whole window annotated A with three V).
  - Exactly one window has a run of exactly three.
  - All four are raw positive.
- **K3 count "9" vs "8".** At 512 intervals Phase 10 gives 8/988; the 9 probably comes from 9/320 at 1,024 intervals.
  The manuscript should give 8/988 with its length.
- **Detrend threshold.** The code gate is 0.7 (`PHASE6_DETREND`); "0.07" is a typo.
- **Table VI shifted strings.** The generated Table VI (`paper_tables/table_VI.md`) is identical, row by row, to Phase 10 report §5.1 (detected counts, replacement, f_min). A shift in the manuscript's copy cannot be confirmed without the draft; it should be replaced by the generated table.

## C2. Figures: pending

---

## What this changes in the manuscript (draft)

The manuscript draft was not available, so the passages are identified by content.

**Methods.**
- Add the sham-handling control:
  > "To separate the effect of removing ectopy from the effect of removing data, each ectopic window was also
  > processed with up to five sham masks (circular shifts of its own mask that avoid ectopic intervals and the 20
  > intervals after each abnormal beat)."
- Report that 121 high-burden windows had no admissible sham.
- State the implementation choices covered by A5 (linear memory, test, ε, noise schedule) and that only readings
  passing the standard-system verification are interpreted.
- State that clock time was unavailable and elapsed time was used.
- K3: give the resolution of every false-positive figure (1/360 s for Phase 9 TEST). Add the degeneracy diagnostic: the
  surrogate prediction-error SD at the m attaining z_max, flagged when < 1e-6.

**Results.**
- Add the transplant: "Inserting one ventricular-type premature beat into each of 439 annotated-clean windows raised
  the titration positive rate from 22.1 % to 86.3 % (Δ = +64.2 pp [59.2, 69.3]); adding heart-rate turbulence changed
  nothing, and turbulence alone did not raise it (−0.7 pp [−1.7, +0.2])."
- Report the sham result with its qualification. Removing the same number of rows at random, with the ectopy left in,
  reproduces about a third of the masked-arm reduction in windows (Q = 0.36) and a fifth in segments (Q = 0.20).
  D stays clearly positive in all five cells.
- Replace the GEE-only burden statement with the within-subject estimate: OR 3.40 [2.62, 4.41], covariate-adjusted.
- Report the CHF effect adjusted for burden and covariates: OR 0.88 [0.46, 1.69].

**Discussion / Limitations, Sec. VII B (the single-beat paragraph).**
- The paragraph should say the single-beat effect is reproduced experimentally, not just observed. It should also
  say the mechanism is the premature beat itself, a single high-leverage event in the polynomial autoregression,
  not heart-rate turbulence.
- It must qualify the Phase 10 statement that isolated synthetic beats rarely trigger titration. That holds at 5-10 %
  density (5/100, 0/100) but not for a single beat in 512 intervals (89/100).
- It should note the unexplained beat-type discrepancy: real single V beats give 47 % and real A beats 76 %, while the
  synthetic beat models give 86 % (V-type) and 54 % (A-type).
- Add: claims (ii)-(iii) depend on an implementation that passes the standard-system verification. With the nested
  partial F-test, which fails verification, the raw group difference disappears and masked CHF falls below NSR.

**Discussion / Limitations, K3 (wherever its specificity is stated).**
- "0/10,100 synthetic null and non-chaotic windows" holds at 1/360 s only.
- At the 1/128 s resolution of nsr2db / chf2db, K3 gave 55/10,092 false positives. It exceeded the registered
  7/100 in three KTz conditions (14, 10, 10 per 100), and 43 of the 55 were degenerate.
- Suggested text: "At 1/128 s, K3 is not size-controlled for rhythms whose beat-to-beat variability is below the
  quantization step (SDNN ≈ 5 ms); the IAAFT surrogates then degenerate. Seven of the eight real detections have
  SDNN 39-202 ms and non-degenerate surrogates; the eighth (nsr029, window 9) is degenerate."
- The K3 power gain at 1/128 s (chaos_primary 10.1 % vs 7.2 %) should not be cited as power. It includes 23
  detections on KTz series whose RR carries no chaos.

**Claims to weaken.**
- K1's burden association: the within-subject OR is 1.42 [0.85, 2.38], not significant. Say "K1's 16 positives
  occurred more often in high-burden windows (GEE), but the association was not significant within subjects."
- "Masking removed 49-66 % of positives": add that part of the reduction is generic data loss (sham Q = 0.36 in
  windows, 0.20 in segments). Editing and Wu-style elimination are not explained by data loss (Q < 0.01), but their
  shams are weak controls.
- K3 "most specific detector, 0/10,100": add the resolution, and the 1/128 s result above.
- "Seven of eight windows ... atrial run of three": remove or correct (4 windows).
- Detrend threshold 0.07 → 0.7; K3 count 8/988 at 512 intervals.

**Supplementary items.**
- Sham-control table (A1).
- Transplant table (A2).
- HRT table.
- A3 model table.
- A4 subject-level metrics.
- A5 verification and claims matrices.
- Seeds, K3 resolution and spike-in tables (Tier 2).

## Regeneration

Run from the repository root on branch `phase11-robustness`, with the environment from
`experiments/phase10_final/requirements-phase10.txt` and `OMP_NUM_THREADS=1`.

```bash
PYTHONPATH=. python -m experiments.phase10_final.run10 --phase conf --download
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.data11 --download nstdb
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.data11 --verify
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.run11 --part tier1 --workers 4
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.run11 --part tier2 --workers 4
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.explore11 --workers 4
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.analysis11 --phase conf
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.explore11 --analyse
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.paper_tables
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.figures11
```
