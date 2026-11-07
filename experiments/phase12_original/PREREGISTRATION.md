# Phase 12 preregistration: implementation fidelity to [BP96] / [PM97], static transforms and surrogates, PM97 replication

Written 2026-11-07, before any Phase 12 run on nsr2db, chf2db, nsrdb or chfdb. The runner (`run12.py`) refuses to run
on these databases unless this file is committed, unmodified and pushed. Later changes go only in dated amendments at
the end.

method: `TIT12`

## 0. Honesty statement
- **Every database here has been analysed before.** nsr2db / chf2db were the Phase 10–11 confirmatory sets; nsrdb /
  chfdb were Phase 9–10 development data. The Phase 10–11 results are known, including titration's group pattern on
  nsr2db / chf2db and the transplant result. Phase 12 therefore makes **registered secondary analyses**: every
  procedure, statistic and decision rule below is fixed before its own outputs exist.
- **Already done before writing this:**
  - FIDELITY.md, the Tier 1A audit, came from reading the papers and the code. It ran nothing on data.
  - The pytest baseline was run, and analysis11 was re-run on stored results as a sanity check.
  - `titration12.py` was written, with its exactness test (22 passed).
  - A **smoke test of the B2 generators**: one seed (12000) of each system, noise-free, R0 decision. All were
    detected except Mackey–Glass τ = 100, which was not detected noise-free. Parameters and sampling were fixed in
    `systems12.py` before that test and have not been changed since.
  - A **pilot** on synthetic pilot seeds 12900+ and mitdb only (`results/pilot/`), for timing. Its outputs were seen.
- **No PM97 outcome has been seen.** The nsrdb / chfdb headers (clock start times) have not been read in Phase 12.
- **Integrity rule.** Nothing is tuned, re-selected or re-defined to rescue a claim. A failed check is reported as
  prominently as a passed one.

## 1. Conventions
- **Code.** Frozen Phase 10–11 files are imported, never edited. New code lives in `experiments/phase12_original/`.
  - R0 = `titration10` defaults: κ ≤ 6, d ≤ 3, linear memory ≤ 83, common target rows n ≥ 83, α = 0.01, F_var.
  - "Positive" for a reading is its **noise-free decision**. Titration's NL > 0 reduces to this (Phase 11 B1: 0 of
    988 windows changed with the noise seed). **R0 on the confirmatory windows is taken from the Phase 10 stored
    results**, as in Phase 11 A5, and recomputed only for the surrogate gate. The agreement between the two is
    reported.
- **Seeds.** Synthetic seeds 12000–12099 (B1 and B2: 12000–12019); pilot seeds 12900+. Entropy 20261107. Every
  random draw comes from `SeedSequence([20261107, part code, unit key, …])` as coded in `run12.py`.
- **CIs.**
  - Subject-cluster bootstrap with 10,000 resamples (percentile), stratified by group, unless stated otherwise.
  - Rates over seeds use Wilson 95 % CIs.
  - GEE uses the Phase 10 model with the Amendment 2 fallback. Conditional logit uses a Wald CI.

## 2. Readings (FIDELITY.md section 4; `titration12.READINGS`)
R0 is the primary in every analysis. The others are opt-in.
- **R0-BP96F** (BP96 F12 and F7/F13):
  - F_calc = ε_lin²/ε_nl² = RSS_lin/RSS_nl, with no degrees-of-freedom adjustment.
  - N = number of points in the series, both in C(r) = ln ε + r/N and in df = (N − r_lin, N − r_nl).
  - Decision: C_nl < C_lin and F_calc above the upper 1 % point of F(df_lin, df_nl).
- **R0-BP96** (adds F16) = R0-BP96F plus a Gaussianity gate:
  - A χ² goodness-of-fit test of the standardised residuals of both best models against N(0, 1), with
    k = ⌈2n^0.4⌉ equiprobable bins and df = k − 3.
  - If either test rejects at 1 % (99 % cutoff), the F-test is replaced by the one-sided Mann–Whitney rank-sum
    test of squared residuals (linear larger) at 1 %.
  - BP96 does not specify the χ² binning or the ranked quantity; these are disclosed choices.
- **R0-S-FT and R0-S-AAFT** (adds F18): the R0 decision AND, for each of K = 10 surrogates s:
  F_s = [min(ε²_lin,s, ε²_nl,s)/df_s] / [ε²_nl,orig/df_orig] > upper 1 % point of F(df_s, df_orig).
  - df uses R0's convention (N = rows). df_s = N − r of the surrogate model attaining the minimum.
  - The worst case over the 10 surrogates decides.
  - Surrogates are fitted on the same rows as the original. FT = phase randomised, with DC (and Nyquist) kept. AAFT
    follows Theiler et al. 1992.
  - Surrogate source series: raw arm, the raw series; edited arm, the edited series; masked arm, the edited
    (interpolated) series fitted on the masked arm's rows.
  - **Surrogate readings are not titrated.** For V1–V6 their NL = the stored Phase 11 R0 NL of the same series and
    seeds if the surrogate requirement holds noise-free, else 0.
- **R0-BP96 and R0-BP96F are titrated for V1–V6** with the reading's test inside the bisection, using Phase 11's
  verification seeds.
- **Not changed** (not specified in BP96): the κ and d ranges, the target rows and the minimum length.

## 3. Tier 1B
### B1. BP96 Fig. 1 ecological model (exact benchmark)
- **Model.** The Fig. 1 legend's model, λ = 118, analysing y.
- **Initial conditions.** x₀ ~ U(10, 100), y₀ ~ U(0.1, 1.0), 1,000 iterations discarded. A draw is replaced (at most
  20 draws) if the orbit is non-finite or has SD < 1e-6.
- **Series.** 1,000 points, seeds 12000–12019.
- **Noise.** One correlated-noise series per seed: an FT surrogate of the clean series, scaled to unit SD. The noisy
  series is y + (p/100)·SD(y)·ξ for p = 0, 5, …, 60 %.
- **Reported per level and reading (R0, R0-BP96F, R0-BP96):**
  - the decision, d_opt > 1, the selected (κ, d) and r of both models;
  - F_calc (plain), the reading's F, df, the F critical value and 1/F_calc;
  - the one-sided Mann–Whitney Z (squared residuals; negative = nonlinear better);
  - the χ² p-values.
- **Detection-loss point** of a seed: the smallest p at which the decision is negative. If none up to 60 %, it is
  counted as 65 % (censored).
- **Pass rule** (each reading; R0 primary):
  1. d_opt > 1 at p = 0 in ≥ 19 of 20 seeds;
  2. median F_calc at p = 0 ≥ 100;
  3. median detection-loss point within 50 ± 15 % (35–65 %).

  All three must hold. A fail is reported prominently as evidence **against the fidelity of that reading**. BP96's
  own values are F_calc ≈ 9,387 against F(0.01, 915, 915) = 1.17, Z = −37.13 and about 50 %.

### B2. BP96 Table 1 (approximate benchmark)
- **Systems, 1,000 points, seeds 12000–12019** (`systems12.py`, fixed before the smoke test). BP96 gives no
  parameters or sampling, so these are standard choices:

  | System | Parameters | Sampling | Why |
  |---|---|---|---|
  | logistic | r = 4 | every iteration | standard fully chaotic map |
  | Hénon | (1.4, 0.3), x | every iteration | standard |
  | Ikeda | u = 0.9, Re z | every iteration | the value named in the task |
  | ecological | the Fig. 1 model, λ = 118 | every iteration | as in B1 |
  | non-chaotic fractal | Grebogi–Ott–Pelikan–Yorke attractor, σ = 1.5, ω = golden mean, x | every iteration | BP96 Table 1 ref. 21 |
  | Rössler | (0.2, 0.2, 5.7), x | Δt = 0.5 (RK4, step 0.01) | about 12 samples per mean cycle |
  | Duffing | double well, δ = 0.25, γ = 0.3, ω = 1, x | 12 samples per forcing period | standard chaotic setting |
  | Lorenz | (10, 28, 8/3), x | Δt = 0.075 | the Phase 10 verification sampling, from the 2001 paper |
  | Mackey–Glass | τ = 100 primary, τ = 17 secondary | Δt = 12.5 | Phase 10 verification generator. τ = 100 because BP96 (p. 217) describes its Mackey–Glass series as having an implicit dimension of seven, which is τ ≈ 100 in Farmer 1982 |

  Series D is a data set that cannot be generated, so it is not run.
- **Statistic.** Per seed and reading (R0, R0-BP96F, R0-BP96): the largest correlated-noise % still detected, by
  bisection on p in [0, 200] % with 8 halvings (resolution 0.8 pp). One noise shape per seed, as in B1. Not detected
  noise-free gives 0; still detected at 200 % gives 200 (censored).
- **Reported:** median and IQR beside the Table 1 value.
- **Consistency rule.** A system is consistent if the median is within ±15 pp of Table 1. Misses are counted and
  marked. B2 is described as approximate in every case.

### B3. Fig. 2 experimental series (optional)
Santa Fe series A and E are used only if they can be obtained from PyPI or GitHub through the allowed network.
GitHub code search goes through `gh api search/code`; files come from github.com or raw.githubusercontent.com. No
other sites are used and no block is worked around.
- If obtained: series A = the first 1,000 points of the laser data; series E = the 1,874 points described in BP96
  Fig. 2b if the segment boundaries are identifiable, otherwise not run.
- Reported: F_calc and 1/F_calc under R0 and R0-BP96, beside 10.23 and 1.08.
- Otherwise: recorded as "not obtainable".

### B4. Static transforms and the surrogate step
**B4.1 Static transforms.**
- **Series.** Gaussian AR(1) with φ = 0.5 and 0.9; n = 1,000 and 512; seeds 12000–12099, 100 series per cell;
  standardised to z.
- **Transforms:**
  - identity (control);
  - z³;
  - exp(z);
  - sign(z)|z|^0.5;
  - tanh(2z);
  - N5 = exp(0.5 z), the Phase 5 warp.
- **Readings:** all five. Reported: positive rate with Wilson CI per (φ, n, transform, reading).
- **Rules:**
  - "R0 flags static transforms" if any non-identity cell has an R0 Wilson lower bound > 5 %.
  - The surrogate step "removes" a cell's positives if that reading's rate is ≤ 5/100 where R0 is ≥ 50/100.
  - **Expectation, stated in advance:** FT surrogates test a linear Gaussian null, so they may keep
    static-transform positives; AAFT surrogates test a monotone static transform of a linear Gaussian process, so
    they may remove them.

**B4.2 Runs under all readings.** R0 on the windows is the stored Phase 10 result.
- (a) Verification V1–V6, under the Phase 10 criteria.
- (b) Phase 10 synthetic conditions: N1–N6, E3 couplets 5 %, E4 runs and S2 isolated 5 %. 512 intervals, seeds
  12000–12099, raw.
- (c) Phase 11 transplant variants V0, V1, V5 and V6 in all 439 clean confirmatory windows. The same event, position
  and draws as Phase 11; times re-quantized to the record's resolution.
- (d) All 988 confirmatory windows, raw / masked / edited (Phase 10 arm definitions).

**B4.3 Claims matrix** (rows = readings):
- (i) window GEE burden OR (raw outcome, x + CHF) > 1, CI excluding 1;
- (ii) raw CHF − NSR > 0 (difference of group mean subject rates, CI excluding 0);
- (iii) masked AND edited CHF − NSR CIs include 0;
- (iv) E3 couplets ≥ 50/100 positive while N1 ≤ 5/100;
- (v) transplant V1 − V0 > 0 (paired, subject-cluster bootstrap CI excluding 0).

**Against my claims:** a reading that passes V1–V6 fails (i), (iii) or (v). (ii) and (iv) are reported separately. A
"reversal" (an OR CI below 1, or a post-handling CI excluding 0) is reported separately.

**Reported plainly:**
- (1) Does R0-S-FT or R0-S-AAFT remove the N5 and static-transform positives (B4.1 rule; B4.2 (b) N5)?
- (2) Does it remove ectopy-driven positives? It "removes" them if V1 − V0 includes 0, or E3 < 50/100.
- (3) The N5 sentence for the manuscript, chosen by this rule:
  - If no surrogate reading that passes V1–V6 removes N5 → "R0, like titration as published in 2001, flags static
    transforms, and so does BP96's surrogate step (FT and AAFT)".
  - If R0-S-AAFT removes N5 and passes V1–V6 but R0-S-FT does not → "R0 flags static transforms; amplitude-adjusted
    surrogates would remove them, phase-randomised surrogates (BP96's description) would not".
  - Any other pattern is described as observed.
- (4) Whether any headline claim changes under a reading that passes V1–V6.

## 4. Tier 1C: PM97 replication on nsrdb / chfdb
- **Data.** nsrdb (18 records) and chfdb (15). Phase 10 dev download, beat annotations as in `data10.beat_series`.
  A record with fewer than 20,000 beat intervals is excluded and listed.
- **Block (20,000 consecutive intervals).**
  - If the header has a clock start time and at least 20,000 intervals follow the first beat at or after midnight,
    the block is those intervals.
  - Otherwise, the block with the lowest mean heart rate among blocks starting at multiples of 1,000 beats. Mean HR
    = 60 / mean of the block's intervals within [0.25, 3.0] s; ties go to the earliest block.
  - Recorded: the rule used, the start and end beat indices and times, and the base time.
- **Preprocessing (PM97).** The block's intervals as annotated, with no range filter, linearly detrended by least
  squares on beat index (mean restored). Nothing else, so ectopic beats stay in.
- **Segments.** 40 contiguous 500-interval segments (primary) and 10 × 2,000 (secondary).
- **Arms per segment:**
  - **raw** (PM97);
  - **rawf**: intervals within [0.25, 3.0] s, the Phase 10 range rule. Secondary; it separates range filtering from
    ectopy handling;
  - **masked**: the K3 rule's rows on rawf, at least 300 rows;
  - **edited**: Phase 6 editing of rawf, NN fraction at least 0.8;
  - **Wu**: masked intervals of rawf deleted, the rest joined, at least 400 left.
  - Masks use the Phase 10 labels: a non-normal or out-of-range beat makes the interval ectopy-related.
- **Test.** The noise-free decision under R0 (primary), R0-BP96 (secondary) and R0-BP96F (component). No titration.
  Secondary: R0 titration NL of every raw 500-beat segment, 10 realisations, one seed per segment.
- **Outcome.** Per subject and arm: detection rate = positive segments / analysable segments. Primary: raw, R0,
  500-beat, all subjects.
- **Subjects:**
  - Primary: all.
  - Secondary, PM97-like selection: within the block, rank subjects by (a) the count of non-normal plus unknown
    beats and (b) the SD of the 40 segment mean RR (in-range intervals). Sum the ranks (ties get the average rank) and
    take the 8 lowest nsrdb and the 11 lowest chfdb subjects; remaining ties go by record name. The selected IDs are
    listed.
- **Statistics per arm and reading:**
  - group mean and median;
  - CHF − NSR difference of means with its bootstrap CI (subjects resampled within groups);
  - two-sided Mann–Whitney p;
  - the number of CHF subjects whose rate is 20–80 pp below the NSR mean, with relative reductions also reported.
- **PM97 pattern reproduced** (primary cell) iff all three hold, each reported:
  1. NSR mean ≥ 0.90;
  2. the CHF − NSR CI lies entirely below 0;
  3. at least 2/3 of CHF subjects are 20–80 pp below the NSR mean.
- **Ectopy:**
  - the same outputs for masked, edited and Wu (and rawf);
  - the within-subject association: conditional logit of the raw R0 500-beat decision on log2(1 + segment ectopy
    count), stratified by subject, all subjects.
  - **Against my claims:** after ectopy handling (masked, edited AND Wu), a CHF − NSR CI entirely below 0 on these
    databases, i.e. PM97's decrease survives ectopy handling.
- **Fig. 3 reading (descriptive).**
  - PM97 does not define the model behind "the first 80 nonlinear coefficients". Registered reading: the full
    κ = 7, d = 3 model (M = 120), fitted by least squares on each standardised raw 500-beat segment that R0 calls
    positive; take its first 80 nonlinear coefficients in Eq. 1 order.
  - Per subject with at least 2 positive segments: the SD of each coefficient across those segments, averaged over
    the 80.
  - Reported: group mean ± s.e., the CHF − NSR difference with a bootstrap CI, and the Mann–Whitney p.
- **If PM97's pattern is not reproduced,** the candidate reasons are listed without claiming any as the cause:
  subject selection, the night window, the detrending choice, chfdb severity, and the implementation readings.

## 5. Tier 2
- `paper_tables12.py` builds every new table from `results/conf/analysis/analysis12.json` only: fidelity,
  B1/B2 benchmarks, the B4 matrix and the PM97 tables.
- `figures12.py` draws B1 (the Fig. 1d analogue), B4 and the PM97 figures.
- `docs/PHASE12_NUMBERS.md` is the numbers ledger for Phases 10–12. Any Phase 10/11 number that Phase 12 changes is
  flagged.
- No manuscript draft was provided, so the manuscript check is not done.

## 6. Runtime budget (pilot: `results/pilot/budget.md`, `budget12.py`)
The pilot used synthetic pilot seeds and mitdb windows and blocks, with 3 workers plus one analysis process. The
projection is mean CPU per task × confirmatory tasks / 4.
- Total **0.19 h**. B2's Mackey–Glass generation (about 6.5 s per series) is not in the pilot, adding about 0.05 h.
  Even a 10× error stays below the **12 h** limit, so nothing is thinned.
- **Thinning order if the total exceeds 12 h:**
  1. R0-S-AAFT on b4w and b4t;
  2. pm97nl, every second segment;
  3. B2 seeds 20 → 10;
  4. the B4.1 n = 512 cells.
- **Never thinned:** the PM97 decisions, the 988 windows under R0-BP96 and R0-S-FT, and B1.

## Amendments
