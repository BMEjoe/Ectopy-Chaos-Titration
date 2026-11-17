# Revision analyses M1 and M2: plan (EXPLORATORY, post hoc)

Written 2026-11-17, before any run of the analyses below. The runner (`run.py`) refuses to start unless this file is
committed, unmodified and pushed (the Phase 3 guard, as in Phases 10-12). Later changes go only in dated amendments at
the end, each stating what had been seen when it was written.

method: `REVM1M2`

## 0. Status and honesty statement

- **Everything here is exploratory.** These analyses answer a mock referee report (concerns M1 and M2). They were not
  registered before the data were analysed. Every database and every synthetic design has been analysed before
  (Phases 10-12), and those results are known and published in the manuscript draft. Every output of this plan is
  labelled "exploratory, post hoc" wherever it is reported.
- **Seen before writing this plan:**
  - the manuscript draft (FINAL_DRAFT_4), BP96, PM97 and PB01 (PMC34630, fetched 2026-11-17), the Phase 10-12
    reports, FIDELITY.md and the Phase 12 numbers ledger;
  - stored NL(r = 3.7) of the logistic map under the Phase 11 readings R1-R7 (0.633-0.671; nested F 1.36) and the
    Phase 12 readings (R0-BP96F 0.673, R0-BP96 0.667);
  - timings only: one VWK test on white noise (11-17 ms) and one series of each Phase 12 benchmark generator.
- **Not seen:** any output of a switch defined below except the four existing readings (R0, R0-BP96F, R0-BP96,
  R0-S-FT), whose Phase 12 results are published.
- **Integrity rule.** No configuration is chosen, dropped or re-defined because of its effect on the ectopy claims.
  R0-ORIG is built by the fixed rule in Section 3.5 from the benchmarks alone, before any headline or M2 run under it.
  Every configuration that is run is reported, including those that weaken the paper's claims.

## 1. Conventions

- **Code.** New code only in `analyses/revision_m1_m2/`; outputs only in `results/revision_m1_m2/`. Frozen Phase 10-12
  modules are imported, never edited. No preregistration, locked script or stored result is modified.
- **Configurable test** (`vwkcfg.py`). One function with explicit switches (Section 3.2). With every switch at its
  default it must reproduce `titration10.vwk_test` / `titrate` exactly. With the Phase 12 switch sets it must reproduce
  `titration12.vwk12` (R0-BP96F, R0-BP96, R0-S-FT) exactly.
- **Regression check first** (`tests/test_regression.py`, pytest). Before any new run, the test must pass on:
  (a) 30 seeded random series, raw and masked, against `titration10` and `titration12` (all output keys equal);
  (b) stored outputs: Phase 12 B1 details (seeds 12000-12002, 13 noise levels, R0 / R0-BP96F / R0-BP96: decision, F,
  p, r_lin, r_nl, C_lin, C_nl), Phase 12 B2 (`max_noise` for 3 systems x 2 seeds x 3 readings), Phase 10 verification
  (logistic r = 3.7, 3.575, 3.565 and 5 scan points: NL and NL_i), Phase 10 synthetic windows (20 tasks: TIT_raw
  positive, p, r_lin, r_nl, NL), Phase 11 exploratory dose-response (k = 1, 26: 10 seeds each, decision), Phase 10
  confirmatory windows (40 windows: TIT_raw / masked / edited positive, p, n_rows, NL), Phase 12 PM97 (2 subjects, 500
  beats, R0 / R0-BP96F / R0-BP96 decisions) and Phase 12 b4w (20 windows, R0-BP96 raw decision). Equality is exact
  (floats compared bit for bit after the stored JSON round trip).
- **"Positive"** for a configuration is its noise-free decision. Titration (NL) is computed only for the logistic
  benchmark (Phase 11 B1: NL > 0 iff the noise-free test is positive; checked again here in the regression test).
- **Seeds.** Entropy for new streams: 20261117. Every new draw comes from `SeedSequence([20261117, part code, key...])`.
  Existing designs keep their own seeds so that R0 reproduces the stored results (listed per analysis).
- **Uncertainty.** Proportions over seeds: Wilson 95 % CIs. Human data: subject-cluster bootstrap, 10,000 resamples,
  percentile, stratified by group where two groups are compared (Phase 10-12 code). Medians and means over seeds:
  bootstrap over seeds, 10,000 resamples. GEE: Phase 10 model with the Amendment 2 fallback, robust Wald CI.
  Conditional logit: Wald CI (Phase 11 A3).
- **Workers.** 4 processes, one BLAS thread each.

## 2. M1.1: deviations table

`deviations.csv` and `deviations.md`, one row per way R0 differs or might differ from BP96, PM97 or PB01 (titration).
Columns: id, aspect, what the original says (paraphrase, short quote, page), R0's choice, the alternative closest to
the original, confidence about what the original did (high / medium / low), decision-relevant (yes / no / by
construction no), switch tested (id or "none" with the reason). Sources: the three papers, FIDELITY.md, Phase 10
METHODS.md, the code. At minimum: N in C(r); F statistic and df; sidedness; linear memory and search range; degree cap
and higher degrees in the figures; embedding dimension; definition of eps and standardisation; matched maximum term
counts; significance level; noise type, repeats and step in titration; parametric vs surrogate test (FT, AAFT, IAAFT);
handling of the first kappa points; the decision rule (d_opt > 1 vs the 99 % test alone); initial conditions and
sampling of simulated systems.

## 3. M1.2: ablation on the benchmarks

### 3.1 Hypothesis (exploratory)
R0 is less sensitive than BP96 on BP96's own benchmarks because of one or more of its deviations. At least one switch
toward the original moves the benchmarks toward the published values. No prediction is made about which.

### 3.2 Configurations
Each switch changes one thing in R0. Defaults in brackets.
- **S1 N = length.** N in C(r) and in df = N - r is the series length [number of regression rows].
- **S2 plain F.** F = RSS_lin / RSS_nl (BP96 Fig. 1 legend) [(RSS_lin/df_lin)/(RSS_nl/df_nl)].
- **S3 chi-square gate.** If a chi-square test (99 %) rejects Gaussian residuals of either best model, the one-sided
  Mann-Whitney test of squared residuals replaces the F-test (Phase 12 R0-BP96 rule) [F-test always].
- **S3b Mann-Whitney only.** The one-sided Mann-Whitney test of squared residuals always (BP96's "alternatively";
  Phase 11 R4). S3 and S3b are alternatives of one deviation.
- **S4 surrogate step (FT).** R0 plus BP96's confirmatory surrogate requirement with 10 phase-randomised surrogates,
  worst case decides (Phase 12 R0-S-FT). Not titrated: NL = R0's NL if the requirement holds noise-free, else 0
  (Phase 12 convention).
- **S5a degree up to 8.** For each kappa = 1..6 the largest d <= 8 with M = (kappa+d)!/(kappa! d!) <= 84 (kappa = 1, 2:
  d = 8; 3: 6; 4: 4; 5, 6: 3), so BP96's kappa = 2, d = 8 (Fig. 1c) and kappa = 3, d = 5 (Fig. 2a) models are reachable
  and the linear and nonlinear searches stay matched at 84 terms [d <= 3].
- **S5b larger embedding.** kappa <= 7, d <= 3 (M = 120), linear memory <= 119 (matched), rows n >= 119. PM97 Fig. 3
  uses "the first 80 nonlinear coefficients", more than R0's largest model has (77). S5a and S5b are alternatives of
  one deviation.
- **S6 zero-padded start.** Every point is a regression target (BP96 Eq. 2 sums over n = 1..N); lags before the start
  are set to 0 (the mean of the standardised series); N = series length [common rows n >= 83].
- **S7 test-only decision.** Positive iff the test rejects at 1 %, without requiring C_nl < C_lin (BP96 Fig. 1d and
  Table 1 measure detection loss with "a 99 % cutoff") [C_nl < C_lin and the test].
- **P1 optimal sampling of flows** (benchmark protocol, Table 1 flows only: Rossler, Duffing, Lorenz, Mackey-Glass).
  BP96 (p. 216) chooses the delay "to maximize the difference between" eps_lin and eps_nl. Sampling step = m x the
  Phase 12 step, m in {0.25, 0.5, 1, 2, 3, 4} (rounded to the integration step). For each flow, m is chosen once on
  design seeds 12900-12904, noise-free, under R0, as the m with the largest median (eps_lin - eps_nl); ties go to the
  smaller m. P1 changes no test and is not part of R0-ORIG.
- **Existing readings:** R0-BP96F (S1 + S2), R0-BP96 (S1 + S2 + S3), and R0-TEXT (S1 + S2 + S3 + S4: everything BP96
  states explicitly).
- **R0-ORIG** (Section 3.5) and its leave-one-out variants.
- Masked rows for any configuration: target and every lag up to the configuration's longest lag free of
  ectopy-related intervals; at least 300 rows (R0 rule). Not changed: alpha = 0.01, one-sided tests, eps = normalised
  RMS, standardisation (decision-invariant), bisection titration with 10 realisations (PB01's scheme).

### 3.3 Benchmarks (every configuration)
- **B-T1, BP96 Table 1.** The 9 systems of Phase 12 B2 (logistic, Henon, Ikeda, ecological, Grebogi non-chaotic
  attractor, Rossler, Duffing, Lorenz, Mackey-Glass tau = 100; tau = 17 secondary), Phase 12 generators, parameters,
  sampling and noise (one FT-surrogate noise shape per seed). Seeds 12000-12099 (100 per system; 12000-12019 are
  Phase 12's). Statistic per seed: the largest correlated-noise % still detected (bisection on [0, 200] %, 8 halvings;
  0 if not detected noise-free; 200 censored). Per system: median, bootstrap CI. Summary: MAD = mean over the 9 systems
  of |median - Table 1 value| (pp), and the number within +/- 15 pp. With P1, the flows use the chosen m.
- **B-ECO, BP96 Fig. 1d.** Ecological model, seeds 12000-12099, correlated noise p = 0, 5, ..., 60 %. Loss point per
  seed = smallest p with a negative decision (65 if none, censored). Summary: median loss point (bootstrap CI);
  detection rate per level (Wilson CI); Phase 12 pass rule (d_opt > 1 at p = 0 in >= 95 % of seeds, median F_calc
  >= 100, median loss point 35-65 %).
- **B-LOG, PB01 Fig. 2.** NL at r = 3.7, logistic map, 100 initial conditions x0 ~ U(0.1, 0.9) (1,000 transient
  iterations, 1,000 points; seeds 0-99), titrated with the configuration's test, 10 noise realisations. Summary: mean
  NL, bootstrap CI. Secondary: the 101-value scan r = 3.50-4.00 of manuscript Fig. 1 (Phase 10 series and noise
  seeds), NL and the V4 criteria. Note: R0 was chosen partly on this value (not independent).
- **B-PM, PM97.** nsrdb, Phase 12 PM97 protocol (20,000-beat night block, linear detrending, 40 x 500 beats), raw
  arm, 18 subjects. Summary: mean subject positive rate, cluster-bootstrap CI. Target: "nearly 100 %". Secondary:
  chfdb, 2,000-beat segments, CHF - NSR.

### 3.4 Gap metric
Published targets: T1 = Table 1 values (MAD target 0); ECO = 50 %; LOG = 0.75; PM = 1.00. Gap of configuration c on
benchmark b: g_b(c) = |value_b(c) - target_b| (T1: MAD). Closure_b(c) = (g_b(R0) - g_b(c)) / g_b(R0) (1 = gap closed,
0 = no change, < 0 = moves away). Composite closure = mean of the four closures. CIs by paired bootstrap (the same
resampled seeds or subjects for R0 and c). "Gap closed" on a benchmark: T1 >= 8 of 9 systems within +/- 15 pp; ECO
median loss point 45-55 %; LOG |NL - 0.75| <= 0.05; PM rate >= 0.90.

### 3.5 Rule for R0-ORIG (fixed now)
1. Every switch S1-S7 (S3/S3b and S5a/S5b as two-option groups) is run alone on all four benchmarks.
2. A switch is "toward the original" iff its composite closure (point estimate) is > 0. In a two-option group, only
   the option with the larger positive composite closure is eligible.
3. R0-ORIG = R0 with every eligible switch. If none is eligible, R0-ORIG = R0, and this is reported.
4. R0-ORIG is run on all benchmarks, with each of its switches removed in turn (leave-one-out), to attribute the
   change. The switch (or combination) with the largest composite closure is reported as explaining most of the gap.
   If no configuration closes the gap on a benchmark (3.4), this is said plainly.
5. The resulting composition is recorded in a dated amendment, committed and pushed before any Section 4 or M2 run
   under R0-ORIG.

## 4. M1.3: headline results under R0, R0-BP96 and R0-ORIG

R0 is recomputed and compared with the stored Phase 10-12 decisions (agreement reported). For "other readings" the
claims check uses the stored Phase 12 results (R0-BP96F, R0-S-FT, R0-S-AAFT) where they exist.
- **H-a synthetic nulls with ectopy (Table III).** Phase 10 synthetic windows, 512 intervals, seeds 10100-10199: N1-N6,
  SETAR, the ten ectopy conditions, and the six non-chaotic van der Pol regimes x {none, isolated 5 %, couplets 10 %,
  bigeminy}; raw, edited and masked arms. Counts with Wilson CIs. The isolated-beat dose (k = 1-51) comes from M2.
- **H-b natural dose-response.** The 988 windows, raw arm: positive rate by burden (0, 1, 2-4, 5-15, >= 16) with
  Wilson CIs; GEE OR per doubling of (1 + burden), adjusted for CHF (Phase 10 Eq. 4); unadjusted CHF OR; conditional
  logit OR with the Phase 11 A3 covariates (mean NN, SDNN, RMSSD, half-window shift, trend ratio, elapsed time) and
  without; segment-level GEE (burden, CHF).
- **H-c single-beat transplant.** The 439 clean windows, Phase 11 positions and draws, V0 (none), V1 (ventricular),
  V5 (atrial), V6 (couplet), re-quantized as Phase 12 b4t. V - V0 in pp, paired, subject-cluster bootstrap CI.
- **H-d group comparison.** CHF - NSR difference of mean subject rates (cluster-bootstrap CI, permutation p): windows
  (raw, masked, edited, Wu-style; the window Wu arm is new: >= 400 intervals left) and 12-min segments (raw, masked,
  edited, Wu-style; Phase 10 rules); paired change raw -> arm by group.
- **H-e PM97 replication.** nsrdb and chfdb, Phase 12 protocol, 500 and 2,000 beats; raw, rawf, masked, edited, Wu;
  group rates, CHF - NSR (bootstrap CI), Mann-Whitney p, PM97 pattern criteria; within-subject conditional logit of
  the raw decision on log2(1 + segment ectopy).

## 5. M1.4: claims check

`claims_check.csv` / `.md`: one row per quantitative claim in the abstract, lead paragraph and conclusions, with its
value under R0, R0-BP96, R0-ORIG and any stored reading, whether it holds under every reading, and replacement
wording where it does not. A directional claim holds under a reading if its CI excludes the null in the claimed
direction (Phase 12 rules (i)-(v)); a claim of "no difference" holds if the CI includes 0; a numeric range holds if the
reading's estimate lies in it. Claims about K3 or the spike-in bounds do not depend on the titration reading and are
marked "not applicable".
- **"In heart rate, titration mainly detects ectopic beats."** Primary metric: the attributable fraction of raw
  window positives, AF = 1 - r0 / r, where r is the positive rate of the 988 windows and r0 that of the burden-0
  windows (subject-cluster bootstrap CI). Secondary: the GEE-based AF = 1 - mean predicted rate at burden 0 / observed
  rate (1,000 cluster-bootstrap refits). Rule: "mainly" is supported under a reading if AF > 0.50, and strongly
  supported if the CI's lower bound is > 0.50. Wording that every reading supports is proposed from the smallest AF and
  the smallest single-beat effect across readings.

## 6. M2: why a single premature beat flips the test, and the inverted U

### 6.1 Hypotheses (exploratory)
- **H2a, few equations.** At small k the positives are driven by the few regression equations whose target or lag
  vector contains the premature interval or the pause. They carry most of the drop in RSS from the best linear to the
  best nonlinear model, have the largest Cook's distances, and removing them (but not an equal number of other
  equations) removes most positives.
- **H2b, downslope by linear absorption.** At high density the short-long alternation becomes linear structure
  (negative lag-1 autocorrelation) that the best linear model absorbs, so the nonlinear model's relative improvement
  shrinks. Periodic placement produces the downslope at lower k than random placement.
- **H2c, beat type.** Any isolated large outlier (spike, missed or extra beat) also produces positives at small k; the
  paired short-long structure is not required. No direction is predicted for the downslope by type.

### 6.2 Design
- Base windows: N1 linear RR (`run11.synth_series("null", "N1_linear_rr", seed)`), 512 intervals, seeds 11000-11199
  (200 per cell; 11000-11099 are the Phase 11 exploratory dose-response windows). Times re-quantized to 1/360 s after
  insertion (k > 0), as Phase 11.
- k in {0, 1, 2, 4, 8, 16, 26, 51}.
- Beat types at the same positions and prematurity draws within a (seed, k, placement):
  - **ventricular** (default): x_i -> c x_i, x_{i+1} -> x_{i+1} + (1 - c) x_i, c ~ U(0.6, 0.8);
  - **atrial** (resetting): x_i -> c x_i, x_{i+1} -> d x_{i+1}, d ~ U(1.0, 1.1);
  - **spike**: x_i -> x_i + s (1 - c) x_i, s = +1 or -1 with equal probability (one additive outlier);
  - **missed beat**: x_i and x_{i+1} merged into one interval (series shortens by k);
  - **extra beat**: x_i split into u x_i and (1 - u) x_i, u ~ U(0.3, 0.7) (series lengthens by k).
- Placement: **random** (Phase 5 rule: indices 1..N-2, at least 3 apart; the ventricular draws equal Phase 11's
  `explore11` for seeds 11000-11099, SeedSequence([20261103, 99, k, seed, 1])); **periodic** for k in {16, 26, 51}:
  spacing s = floor((N - 2)/k), offset ~ U{1..s}, positions offset + j s.
- Other draws: SeedSequence([20261117, 60 + type code, placement code, k, seed]).
- Readings: R0 and R0-ORIG (all analyses); R0-BP96 (positive rate only).
- Ectopy-affected intervals: ventricular, atrial, extra: {i, i+1}; spike, missed: {i}. Affected equations (rows):
  target or any lag 1..6 (the nonlinear memory) on an affected interval; secondary: lags 1..kappa_nl.

### 6.3 Analyses
1. **Model anatomy.** In every window: kappa_lin, r_lin, kappa_nl, d_nl, r_nl and the terms of the best nonlinear
   model (lag sets). Frequency of each term among positive windows by k. Terms are flagged as powers of one lag or
   products of adjacent lags (which see a premature interval and its pause together), and by their "ectopic
   concentration" (share of the term's sum of squares from affected rows).
2. **Where the improvement comes from.** Delta = RSS_lin - RSS_nl = sum over rows of (e_lin^2 - e_nl^2), split into
   affected and other rows; share from affected rows and its ratio to the share of affected rows (enrichment), by k.
3. **Influence.** (a) Refit with affected rows removed (model selection and test redone; >= 300 rows else not
   analysable); (b) sham: the same number of rows removed at k random non-ectopic pseudo-positions (same rule);
   (c) leverage and Cook's distance of every row under the best nonlinear and best linear model; refit after removing
   the top m rows by Cook's distance (nonlinear model), m in {1, 2, 4, 8, 16}; share of the top 8 that are affected.
   For k = 1 also the R0 masked rule (lags up to 83).
4. **Downslope.** By k: lag-1 and lag-2 autocorrelation; kappa_lin, eps_lin^2, eps_nl^2; relative improvement
   1 - RSS_nl/RSS_lin; normalised RSS of each model on affected and other rows (RSS_group / SS_group); mean max |z| of
   the affected intervals (dilution check); C(r) curves of the linear and best nonlinear family for representative
   windows (seeds 11000-11002, every k, ventricular, random and periodic).

### 6.4 Decision rules (fixed now)
- **"Driven by a few high-leverage equations"**: ventricular, random, R0, k in {1, 2, 4, 8} pooled, raw-positive
  windows analysable after removal. Stated if >= 80 % become negative after removing the affected equations and
  <= 20 % after the sham removal. "Largely driven" if 50-80 %; otherwise "not driven". The same rule for each beat type
  and for R0-ORIG. Cook's: the smallest m at which >= 80 % of positives become negative, by k.
- **H2b supported** (ventricular, random, R0) if all hold: (a) mean lag-1 autocorrelation falls from k = 8 to k = 51
  (bootstrap CI of the difference excludes 0) and is monotone over k = 8, 16, 26, 51; (b) the linear model's normalised
  RSS on affected rows falls from k = 8 to 51 (CI excludes 0); (c) the relative improvement falls from k = 8 to 51
  (CI excludes 0). "Periodic gives the downslope at lower k" if the periodic positive rate at k = 16 is below the
  random one (Newcombe CI of the difference excludes 0). Alternatives reported with the same diagnostics: dilution
  (max |z| of affected intervals falls with k), nonlinear-model deterioration (its affected-row normalised RSS rises
  with k), change of selected model size. The simplest explanation consistent with all diagnostics is stated; if
  none fits, the alternatives are listed.

## 7. Pilot, runtime and thinning
- Pilot: every part on a small subset with pilot seeds (synthetic 12950-12952; data parts on mitdb with the Phase 12
  pilot rule), 4 workers; per-task CPU is projected to the full design and printed. Pilot outputs are kept in
  `results/revision_m1_m2/pilot/` and never analysed.
- Budget 12 h wall. If exceeded, thinning in this order: (1) S4 and R0-TEXT on B-T1 to 50 seeds; (2) the logistic scan
  only for R0, R0-BP96 and R0-ORIG; (3) M2 periodic placement for spike, missed and extra beats. Never thinned: the
  regression check, B-PM, the headline analyses under R0-ORIG, M2 ventricular.

## 8. Deliverables
`results/revision_m1_m2/`: REPORT.md (plan commit hash and timestamp at the top), deviations.csv/.md,
ablation_benchmarks.csv, headline_by_reading.csv, claims_check.csv/.md, m2_term_selection.csv, m2_ssr_share.csv,
m2_influence.csv, m2_downslope.csv, figures (PNG 300 dpi + PDF, Phase 11 style), manuscript_text.md, raw JSONL,
environment record. One command regenerates everything: `python analyses/revision_m1_m2/run_all.py`.

## Amendments

### Amendment 1 (2026-11-17, 16:05 UTC): verification controls for every configuration (added before any R0-ORIG run)
- **Seen when written:** the ablation was about 30 % complete. I looked at interim, incomplete summaries of the
  ecological loss point, NL(3.7) and the nsrdb rate by configuration (20-35 of 100 seeds, 2-9 of 33 subjects). The
  zero-padded start (S6) moved the ecological loss point to the censoring limit (65 %), past the published ~50 %.
  No R0-ORIG, headline or M2 output under any new configuration existed.
- **Reason:** a configuration chosen for sensitivity on BP96's benchmarks must also be checked for false positives.
  The plan had no size check for the new configurations.
- **Added (part `ver`):** the Phase 10 verification controls V5 (white noise, AR(2), limit cycle, 2-torus, logistic
  period 3; 20 series each) and V6 (Lorenz r = 28 and 160, Henon, Mackey-Glass tau = 17, 23, 30 and 14, 15, 16), same
  series and noise seeds as Phase 10, for R0, every single switch, R0-BP96F, R0-BP96, R0-TEXT and R0-ORIG. V1-V4 come
  from the logistic scan already in the plan. Criteria are Phase 10's (METHODS.md 2.5). A configuration that fails
  them is reported as such; R0-ORIG is still built by the 3.5 rule (no change) and still run, with the failure stated
  next to every result under it.
- Regression check extended: `verify_control` reproduces the stored Phase 10 control NL under R0 (test added, passed).

### Amendment 2 (2026-11-17, after the single-switch ablation; before any R0-ORIG run)
- **Seen when written:** the complete single-switch ablation (13 configurations, 22,542 tasks, 0 errors;
  `ablation_benchmarks.csv`, `ablation_table1_by_system.csv`, `p1_choice.json`).
- **R0-ORIG by the 3.5 rule** (`r0orig.json`): composite closures S1 +0.10, S2 +0.14, S3 -0.09, S3b +0.08, S4 -0.08,
  S5a +0.09, S5b +0.03, S6 -0.30, S7 0.00. Eligible: S1, S2, S3b, S5a, S5b; S5a beats S5b in its group.
  **R0-ORIG = S1 + S2 + S3b + S5a** (N = series length; plain F ratio replaced by the Mann-Whitney test as the
  decision test, since S3b replaces the F-test; degree up to 8 with M <= 84). Leave-one-out variants:
  R0-ORIG minus S1, S2, S3b, S5a. Note: with S3b in R0-ORIG the decision test is Mann-Whitney, so S2 only affects
  the reported F (and the surrogate F convention), not the decision.
- **P1** (sampling of flows): m = 4 (Rossler, Duffing, Mackey-Glass tau = 100), 3 (Lorenz), 1 (Mackey-Glass tau = 17);
  three of four choices are at the edge of the registered grid. P1 is a benchmark protocol, not part of R0-ORIG.
- Correction: Amendment 1's heading time should read 16:01 UTC (commit 4e137ea).
- Nothing else changes. Run order: head, m2, ver, abl2.
