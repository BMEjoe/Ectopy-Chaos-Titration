# Phase 11 preregistration: robustness analyses for the ectopy / noise-titration manuscript

**Status of the data (read first).** The confirmatory databases (nsr2db, chf2db) and every Phase 10 output on them
(`experiments/phase10_final/results/conf/`, `docs/PHASE10_FINAL_RESULTS.md`) **have already been seen**. The
analyses below are therefore **registered secondary robustness analyses**. They are fixed here before their own
outputs exist, but they **cannot produce new confirmatory claims**; they can only strengthen, qualify or weaken the
Phase 10 claims. The runner (`run11.py`, every confirmatory part) refuses to start unless this file is committed,
unmodified and pushed (Phase 3 guard). Later changes are dated amendments at the end, each stating what had been
seen.

Before this file, Phase 11 code was run only on DEVELOPMENT data (nsrdb, chfdb) and on pilot seeds 11900-11901, for
the runtime budget (Section 9). No Phase 11 output on confirmatory data exists.

Methods listed for the runner guard:
- method: `TIT`
- method: `K3`
- method: `K1`
- method: `HRT`

**Integrity rule.** No detector, threshold, titration default or Phase 10 definition is changed. Nothing below is
tuned, re-selected or re-defined after its output is seen. A result that weakens the paper is reported with the
same prominence as one that supports it.

**The hypothesis under test (H).** "Titration positives in heart rate are driven largely by ectopy." For every
analysis: the prediction under H, and the result that counts against H.

## 1. Conventions

- **Data, windows, segments, labels, masks, arms**: exactly as Phase 10 (`PREREGISTRATION.md` there; code
  `harness10.py`, `run10.py`). 988 windows of 512 intervals (83 subjects; 549 with >= 1 ectopy-related interval,
  439 "clean" = burden 0 and no masked interval), 8,993 analysable 12-min segments. Phase 10 results are taken from
  `results/conf/real.jsonl` and `seg.jsonl` (verified reproducible: the Phase 11 sanity check regenerated
  `analysis.json` identically, and titration of stored windows reproduces the stored NL exactly).
- **Titration**: `titration10` (frozen, verified) through `harness10.titration` (seed = data CRC + salt) unless
  stated. Variant readings through `titration11.py`, which at the default reading reproduces `titration10.titrate`
  exactly (test on 20 seeded series, raw and masked).
- **"Usable"** arm result: masked arm analysable (>= 300 regression rows); edited arm eligible (NN fraction >= 0.80);
  Wu arm eligible (>= 400 NN intervals).
- **Positive** = NL > 0 (Phase 10).
- **CIs**: subject-cluster bootstrap, 10,000 resamples, percentile 95 % (resampling subjects within group for
  group-specific estimates and group differences; over all subjects for pooled estimates), bootstrap RNG seed
  20261104 unless stated. Window/segment-level proportions are ratio estimators (sum over resampled units).
  Wilson 95 % intervals for simple proportions; exact McNemar = two-sided binomial test on discordant pairs.
- **Seeds**: Phase 11 synthetic seeds 11000-11099 (power groups 11000-11029 / 11000-11009); RNG entropy 20261103 for
  every Phase 11 draw (`harness11.ss`); never an earlier seed. Pilot seeds 11900+ are never used in a reported part.
- **GEE**: Phase 10 model and Amendment 2 fallback (exchangeable; independence if non-finite), statsmodels 0.15.0;
  not fitted with < 10 positives or negatives.

## 2. A1: sham-handling control

**Question.** Do masking / editing / Wu elimination remove titration positives only because they throw data away?

**Procedure** (`run11 --part sham_win, sham_seg`).
- Windows with >= 1 ectopy-related interval (549): 5 sham masks per window, each a circular shift of the window's
  true K3 interval mask by an offset o ~ U{1, ..., 511} (seeded per window), accepted only if no sham-masked interval
  is a true ectopy-related interval or one of the 20 intervals after any abnormal beat (heart rate turbulence zone;
  `harness11.forbidden`); offsets distinct; up to 200 draws per sham, else "sham unavailable" (fewer shams).
  Arms: (a) masked arm with the sham mask (the real ectopy stays in the series and its rows are used);
  (b) edited arm with the sham mask (linear interpolation at the sham positions; eligibility by the sham mask's NN
  fraction). Titration seed as the Phase 10 arms (data CRC, salt 0), i.e. the same noise draws as the true arm.
- Clean windows (439): 5 donor masks per window, each the true mask of a different, randomly chosen (seeded) window
  of the same group with >= 1 masked interval, applied at the same interval positions; masked and edited arms.
- 12-min segments: the same with 2 shams per segment (circular shifts within the segment's in-range interval series;
  forbidden zone computed on the segment's beats), arms masked, edited and Wu-style (sham positions deleted, remaining
  intervals joined, >= 400 required). Clean segments (no masked interval): 2 donor masks from segments of the same
  group (applied at the same positions; truncated, or padded with unmasked intervals, to the recipient's length).
- Outputs: `results/conf/sham_win.jsonl`, `sham_seg.jsonl`.

**Primary statistic** (per level: windows, segments; per arm: masked, edited, [segments] Wu; by group and pooled).
Units = ectopic units that are raw positive (Phase 10 raw arm), usable in the true arm (Phase 10 stored result) and
with >= 1 usable sham. Per unit u: R_u = 1 if the true arm is negative; S_u = mean over u's usable shams of 1(sham arm
negative). F_true = mean R_u, F_sham = mean S_u, **D = F_true - F_sham** (paired), with its cluster-bootstrap CI;
also Q = F_sham / F_true (point and CI).
**Secondary:** clean units: fraction of raw-positive clean units made negative by donor-mask handling (data loss in
the absence of ectopy), per arm; sham-unavailable counts; positives created by sham handling among raw negatives.

**Prediction under H:** D > 0 with CI excluding 0, and Q < 0.5, in every pooled (level, arm) cell; clean-unit
removal fraction small.
**Against H:** in a pooled (level, arm) cell, CI(D) includes 0, or Q >= 0.5 (point estimate). Each of the five pooled
cells (windows masked, edited; segments masked, edited, Wu) gets its own verdict; the A1 summary verdict is
"supports" only if all five support.

## 3. A2: the single-beat result

### A2.1 Heart rate turbulence in the confirmatory data (DESCRIPTIVE; not a detector parameter)
`hrt11.py`, `run11 --part hrt`, full recordings.
- Usable isolated VPB (symbol V at beat v): beats v-6..v-1 and v+1..v+21 labelled normal; the 5 sinus intervals before
  the coupling interval and the 20 after the compensatory pause all in [300, 2000] ms; consecutive sinus intervals
  differ by <= 200 ms; every sinus interval within 20 % of the reference (mean of the 5 preceding sinus intervals);
  coupling interval <= 0.8 x reference; compensatory pause >= 1.2 x reference (Schmidt-group filter rules, METHODS.md).
- TO (%) per VPB = 100 [(RR1 + RR2) - (RR-2 + RR-1)] / (RR-2 + RR-1); subject TO = mean over the subject's VPBs.
  TS (ms per RR interval) = maximum least-squares slope over any 5 consecutive intervals within RR1..RR15 of the
  subject's averaged tachogram (consensus); TS20 within RR1..RR20 (Schmidt 1999) secondary.
- Reported per group: subjects with >= 1 and >= 5 usable VPBs; median [IQR] TO and TS over subjects with >= 5 usable
  VPBs (and >= 1, secondary); fraction with TO >= 0 % / TS <= 2.5 ms/RR; the group-averaged normalised tachogram
  h_j = RR_j / mean(RR-2, RR-1) - 1 (j = 1..20), averaged per subject, then over subjects with >= 1 usable VPB (equal
  subject weight). Compared descriptively with Grimm et al. (2003) healthy medians TS ~ 13 ms/RR, TO ~ -2.5 % and the
  cut-offs TS <= 2.5 ms/RR, TO >= 0 %. No test.

### A2.2 HRT generator (opt-in; `hrt11.insert_event`)
After the compensatory interval, the next 20 sinus intervals are multiplied by (1 + h_j). Profiles: HRT-empirical
(the group's averaged normalised tachogram from A2.1); HRT-normal (h = -a exp(-(j-1)/1.5) + b exp(-((j-7)/4)^2),
a, b >= 0 solved so that TO = -2.5 % and TS = 13 ms/RR on a constant baseline equal to the event's pre-VPB interval
mean(x_{j-2}, x_{j-1})); HRT-blunted (same form, TO = 0 %, TS = 2.5 ms/RR). Tests (`tests/test_hrt11.py`, passing):
TO and TS measured by `hrt11`'s own code on generated profiles and on a generated tachogram hit the targets (+/- 0.05).

### A2.3 Transplant into real clean windows (KEY EXPERIMENT)
`run11 --part transplant`. Every confirmatory clean window (burden 0, no masked interval). One event per window at a
seeded position j (premature interval index) with 30 <= j and j + 22 <= 512 - 31 (the whole event, including the 20
HRT intervals, >= 30 intervals from both ends); draws c, c' ~ U(0.6, 0.8), d ~ U(1.0, 1.1) per window, shared by all
variants (paired design). Beat times re-quantized to 1/128 s for V1-V7.
- V0 untouched; V1 one ventricular-type beat (x_j -> c x_j, compensatory x_{j+1} -> x_j + x_{j+1} - c x_j, i.e. full
  compensatory pause, Phase 6 model), no HRT; V2 V1 + HRT-empirical of the window's group; V3 V1 + HRT-normal; V4 V1 +
  HRT-blunted; V5 one atrial-type beat (x_j -> c x_j, x_{j+1} -> d x_{j+1}; Phase 6 E5); V6 one couplet (x_j -> c x_j,
  x_{j+1} -> c' x_j, x_{j+2} -> x_j + x_{j+1} + x_{j+2} - (c + c') x_j); V7 HRT-empirical with no premature beat.
- Titration raw arm, with ONE seed per window shared by V0-V7 (SeedSequence([20261103, 3, window key, 2])), so that
  the variants differ only in the series. (V0 is therefore re-titrated with this seed, not the Phase 10 seed.)
- **Statistic** per variant v = V1..V7, by group and pooled: Delta_v = rate(v) - rate(V0) (paired, window level),
  subject-cluster bootstrap CI; exact McNemar p. **Primary:** V1 and V2 pooled.

### A2.4 Synthetic single beats
`run11 --part synth1`. N1 (Phase 5 linear RR) and the six non-chaotic coupled-vdP regimes (Phase 10 Q2e generator,
1/360 s), 512 intervals, seeds 11000-11099 (100 windows per condition); the same V0-V7 events (HRT-empirical = the
NSR profile; re-quantized to 1/360 s). Statistic per (condition, variant): rate, Wilson CI, Delta_v vs V0 with a
window bootstrap CI (10,000) and exact McNemar.

### A2.5 Extended masking
`run11 --part extmask`. Real windows with burden 1-4: masked arm with the K3 mask extended to cover the 20 intervals
after every abnormal beat. Units: windows analysable under both masks. Statistic: positive rate (standard) - positive
rate (extended), all units and among raw positives, cluster-bootstrap CI, exact McNemar.

### A2.6 By beat type
Burden-1 windows by the annotation symbol of the single abnormal beat: V; supraventricular (A, a, J, S); other or
unknown (every other non-normal beat symbol, including Q and ?). Plus burden-0 windows whose only masked intervals come
from the 0.25-2.5 s interval rule, and the clean reference. Raw titration positive rate (Phase 10) with Wilson CIs.

**Prediction under H (A2):** inserting one ectopic event raises the clean-window positive rate (Delta > 0, CI
excluding 0, for V1 and/or V2), most for the couplet (V6); and/or extended masking removes additional positives.
**Against H:** no transplanted variant raises the positive rate (every CI includes 0) AND extended masking removes
nothing (CI includes 0). Then a single real beat does not cause the 22 % -> 63 % jump, and the association must come
from confounding (A3) or annotation problems; this would be reported plainly and the manuscript's single-beat
statement revised.

## 4. A3: confounding (window and segment level)

**Covariates** (`harness11.nn_covariates`, from NN intervals only = intervals not in the K3 mask): mean NN (s); SDNN of
NN intervals (log); RMSSD of adjacent NN pairs (log); half shift |mean NN 2nd half - mean NN 1st half| / SD NN; K1
detrend ratio (`final_pipeline.linear_trend_ratio` of the NN series; log(1 + ratio)); elapsed hours since the first
beat in 4 bins [0, 6), [6, 12), [12, 18), >= 18 (reference [0, 6)). **Clock time is unavailable** (no header start
time). Continuous covariates are standardised (z-scores over the model sample); units with an undefined covariate are
dropped (count reported).
**Exposure:** x = log2(1 + burden). **Outcome:** titration raw positive.
**Models:**
1. GEE (Phase 10 + Amendment 2 fallback): outcome ~ x + CHF + covariates.
2. **Conditional logistic regression stratified by subject** (`statsmodels ConditionalLogit`): outcome ~ x + covariates
   (PRIMARY for A3), and outcome ~ x alone. Wald 95 % CI; secondary: subject-cluster bootstrap percentile CI (2,000
   resamples). Only subjects with both outcomes contribute (count reported).
3. Models 1-2 within NSR only and within CHF only.
4. CHF effect from model 1 (adjusted for burden + covariates).
- Segment level: the same (segment burden, segment covariates, elapsed hours of the segment start).
- Secondary outcomes (window level): LLE alone and K1 (same models; sparse-outcome rule).
**Prediction under H:** within-subject burden OR > 1 with CI excluding 1 at window and segment level.
**Against H:** the within-subject burden OR (model 2 with covariates, windows, pooled) has a CI that includes 1.

## 5. A4: the original studies' own metrics (descriptive)

From Phase 10 segments, by group and arm (raw, masked, edited, Wu): fraction of subjects with >= 1 positive segment and
with >= 10 %, 25 %, 50 % of their usable segments positive (Wilson CIs); Wu et al.'s detection rate per 3-hour block
(elapsed time 0-3 h, ..., 21-24 h; mean subject block rate, cluster-bootstrap CI); mean NL among positive segments.
Poon & Merrill (1997) PDF: **not available** at preregistration -> their exact metric, segment length and preprocessing
are not extracted (if the PDF appears, the extraction is added by amendment before A4 is analysed).
**Prediction under H:** on raw RR almost every subject in both groups has >= 1 positive segment (so a per-subject
"chaos in all healthy subjects" statement is reproducible but uninformative); the per-subject fractions fall after
ectopy handling, more in CHF. No "against" criterion (descriptive).

## 6. A5: sensitivity of titration to the implementation readings

**A5.1** Barahona & Poon 1996 PDF: **not available** at preregistration -> definitions not extracted; METHODS.md
choices 1-3 remain "not addressed by the primary source"; **R8 is not run**. (If the PDF appears before the A5 runs
start, the extraction and R8 are added by amendment.)
**Readings** (`titration11.READINGS`): R0 primary (frozen); R1 linear memory <= 6; R2 <= 20; R3 <= 50; R4 Mann-Whitney
test on squared residuals; R5 nested partial F-test; R6 eps as a variance (EPS_POWER = 1); R7 noise increments of
1 % of sigma_y instead of bisection. Masked-arm row span = max(6, linear memory).
**Runs:** (a) verification V1-V6 (`verify_titration` systems, seeds and Phase 10 pass criteria, METHODS.md 2.5) for every
reading; (b) all 988 windows raw / masked / edited for R1-R7 (R0 = Phase 10 stored results); (c) synthetic N1, E3
couplets 5 %, E4 runs, N5 static warp, vdP (5.45, 5.6), 100 windows each, seeds 11000-11099, raw arm, R0-R7; (d) R1 and
R4: all segments raw / masked / Wu. Failing readings are kept and flagged.
**Matrix** (rows readings; columns claims):
(i) window GEE burden OR (outcome raw, model x + CHF) > 1 with CI excluding 1;
(ii) raw CHF - NSR > 0 (difference of group mean subject rates, CI excluding 0): windows for all readings; segments
for R0, R1, R4;
(iii) post-handling CHF - NSR CI includes 0: windows masked AND edited arms; segments masked AND Wu for R0, R1, R4;
(iv) couplets produce false positives on N1: E3 couplets >= 50/100 positive while N1 <= 5/100.
**Prediction under H:** every reading that passes V1-V6 reproduces (i)-(iv).
**Against H:** a reading that passes V1-V6 but fails (i) or (iii) at window level. ("Reverses" = OR CI entirely below
1, or a post-handling difference whose CI excludes 0; reported separately.)

## 7. Tier 2

**B1 seeds.** Titration salts 1 and 2 (raw, masked, edited on all windows; raw, masked, Wu on all segments); K3 with
`masked_growth_seed_entropy` 11101 and 11102 on all windows; K1 with (`lle_chaos_test_seed_entropy`, `random_seed`) =
(11201, 11211) and (11202, 11212) on all windows. Reported: each headline number (P2 raw / masked / Wu group differences,
P2 masked change by group, P3 OR and CHF OR, burden-bin rates, K3 / K1 / LLE / UPO counts) at salt 0 and each new seed,
its range, and the number (%) of unit decisions that flip vs Phase 10. Prediction: ranges small, no conclusion changes.
Weakens: any headline conclusion (CI excluding / including its null) changes under a new seed.

**B2 K3 at 1/128 s.** Phase 9 TEST groups (`run_test9.tasks` with new seeds: FP groups null_input, null_real,
noncha_input, noncha_real 11000-11099; power groups chaos_primary and chaos_second 11000-11029 / DEV MG 11000-11009),
same generators (TEST noise channel). K3 (frozen) on the native 1/360 s intervals and on intervals re-quantized to
1/128 s (beat times rounded). Pass rule (Phase 9): <= 7 detections per 100 windows in every FP condition, at 1/128 s
(primary) and 1/360 s. Power: detection rate of chaos_primary / chaos_second at both resolutions. **Degeneracy
diagnostic** (`harness11.k3_diag`, identical computation to the frozen K3 with the surrogate error SDs exposed):
degenerate iff the surrogate prediction-error SD at the m attaining z_max is < 1e-6. For the 8 real K3 detections:
SDNN, z_max, the diagnostic. Weakens K3 claims: any FP condition > 7/100 at 1/128 s, or degenerate real detections
beyond nsr029 window 9.

**B3 variance-preserving spike-ins.** x'_k = mu_NN + sqrt(1 - f)(x_k - mu_NN) + sqrt(f) sigma_NN c_k (every interval),
same base windows, chaotic components (same keys), f grid and 1/128 s quantization as Phase 10; K3. pi_upper and f_min
(pi < 0.05, 0.20; Phase 10 primary formula with Phase 10's real K3 bound U) beside the primary design; out-of-range
interval counts. K1 only if budget remains (Section 9).

**B4 second base window** (only if budget remains): j' = 1 + ((s + 6) mod 12) (next available cyclically, never the
Phase 10 base window); additive design; K3. Pooled with Phase 10 (166 base windows per condition): p_emp, Wilson lower
bound, pi_upper and f_min (primary and Wilson-lower variants).

## 8. Tier 3

C1 `paper_tables.py`: Tables I-VII (Phase 10 mapping: I detectors and parameters; II titration verification; III
synthetic false positives (Q2e); IV P2 segments by arm; V Q3 burden / GEE; VI Q1 K3 exclusion; VII Q1 K1 exclusion)
and supplementary tables A1-A5, B1-B3, built from the analysis JSON only (markdown + CSV). Manuscript check only if
`private_inputs/manuscript_draft.pdf` exists. C2 figures F1-F7 (vector PDF + 300-dpi PNG, colour-blind-safe palette,
3.37 in / 6.69 in widths).

## 9. Runtime budget (pilot)

Pilot (`run11 --pilot --part all`, `budget11.py`): ~2 % of each part's tasks on DEVELOPMENT data (nsrdb, chfdb:
every 50th window / spike task, every 16th subject for segment parts) and pilot seeds 11900-11901 (synthetic parts;
B2 1 seed per condition, every 3rd task). CPU seconds per task were measured with 4 workers busy (plus one analysis
process, so the projection is conservative). Projection = mean x confirmatory task count / 4 workers. 0 task errors.

| part | tier | pilot tasks | mean CPU s / task (4 workers busy) | confirmatory tasks | projected wall (4 workers) |
|---|---|---|---|---|---|
| prep | 1 | 1 | 36.4 | 1 | 0.00 h |
| hrt | 1 | 33 | 0.4 | 83 | 0.00 h |
| covar | 1 | 33 | 0.5 | 83 | 0.00 h |
| sham_win | 1 | 8 | 4.5 | 988 | 0.31 h |
| transplant | 1 | 5 | 5.4 | 439 | 0.16 h |
| synth1 | 1 | 14 | 12.7 | 700 | 0.62 h |
| extmask | 1 | 2 | 1.0 | 232 | 0.02 h |
| rdg_verify | 1 | 72 | 1.8 | 1680 | 0.21 h |
| rdg_synth | 1 | 80 | 2.6 | 4000 | 0.71 h |
| rdg_win | 1 | 56 | 2.5 | 6916 | 1.19 h |
| sham_seg | 1 | 3 | 336.5 | 83 | 1.94 h |
| rdg_seg | 1 | 6 | 263.0 | 166 | 3.03 h |
| seed_win | 2 | 16 | 1.8 | 1976 | 0.24 h |
| k3seed | 2 | 16 | 2.4 | 1976 | 0.33 h |
| b2 | 2 | 63 | 3.4 | 12290 | 2.92 h |
| b3 | 2 | 11 | 15.1 | 1328 | 1.40 h |
| seed_seg | 2 | 6 | 197.7 | 166 | 2.28 h |
| k1seed | 2 | 16 | 7.2 | 1976 | 0.99 h |
| b4 | 2 | 11 | 13.9 | 1328 | 1.28 h |
| **total** | | | | | **Tier 1 8.2 h + Tier 2 9.4 h = 17.6 h** |

- The projected total (17.6 h) is below 24 h: **no thinning**. B3 with K1 (`b3k1`, ~7 s per K1 call x 9,296 series
  ~ 4.5 h) is NOT part of the budget; it runs last, only if every other Tier 1-2 part is complete.
- B4 is in the budget and runs last among the budgeted parts.

**Thinning rule** (if the projected total exceeds 24 h, applied in this order until it does not): (1) drop B4;
(2) drop K1 from B3, then B3's f = 0.05 and 0.7; (3) B2 FP seeds 100 -> 50 per condition (pass rule <= 3 per 50), power
groups 30 -> 15; (4) B1 segment salt 2; (5) A5 segments R4. **Never thinned: the real-data windows and segments of
any analysis that runs.** Order of the run: Tier 1 (prep, hrt, covar, sham_win, transplant, synth1, extmask,
rdg_verify, rdg_synth, rdg_win, sham_seg, rdg_seg), then Tier 2 (seed_win, k3seed, b2, b3, seed_seg, k1seed, b4).
Any thinning is recorded as an amendment.

## 10. Outputs

`results/conf/<part>.jsonl` (one line per task), `results/conf/analysis/analysis11.json` (`analysis11.py`), tables
(`paper_tables.py`), figures (`figures11.py`), report `docs/PHASE11_ROBUSTNESS_RESULTS.md`.

## Amendments
