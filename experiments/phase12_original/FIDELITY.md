# Phase 12 Tier 1A: fidelity of the titration implementation to [BP96] and [PM97]

**Sources** (read from the PDFs in the git-ignored `private_inputs/`; paraphrased, quotes kept to a few words):
- [BP96] Barahona M, Poon C-S. Nature 381, 215–217 (1996). This is a scanned PDF with no text layer. The pages were
  rendered at 220 dpi and read column by column, with the Fig. 1b and 2a panels zoomed.
- [PM97] Poon C-S, Merrill CK. Nature 389, 492–495 (1997). Read from its text layer.

**Code audited:** `experiments/phase10_final/titration10.py` (frozen, unchanged). The registered primary reading
**R0** is this code at its defaults (κ ≤ 6, d ≤ 3, linear memory ≤ 83, α = 0.01, F_var test).

**Verdicts:**
- **matches**: the code does what the source says.
- **deviates**: the code differs from the source. "Decision-relevant" means the difference can change a yes/no
  result.
- **not specified**: the source does not fix the point, and the code's choice is disclosed.

## 1. BP96: line by line

| # | Point | Source location | What the source says (paraphrased) | What the code does (`titration10.py`, line) | Verdict |
|---|---|---|---|---|---|
| F1 | Model form and term order | p. 215, Eq. 1 | Closed-loop Volterra–Wiener–Korenberg autoregression on y_{n−1..n−κ}. Terms in order: constant; linear lags 1..κ; then y_{n−1}², y_{n−1}y_{n−2}, …, up to y_{n−κ}^d. M = (κ+d)!/(d!κ!). | `design` (54–61): constant column, then for degree 1..d the lag products in `combinations_with_replacement` order. | **matches**. Eq. 1 shows only the first terms of each degree; the lexicographic order within a degree agrees with them. |
| F2 | Truncated models | p. 216, Eq. 3; Fig. 1b (open squares = truncated, filled = complete) | C(r) is evaluated for r = 1..M leading terms of a {κ, d} expansion. | `leading_rss` (64–89) gives the residual sum for every leading r; C is computed for every r (105, 113). | **matches** |
| F3 | Constant term | p. 215, Eq. 1 (a₀); Figs. 1b and 2a start at C = 0 for r = 1 | a₀ is a term and counts in r. The constant-only model has ε = 1. | Column of ones (57); r counts from 1. | **matches** |
| F4 | Coefficient estimation | p. 215 | Korenberg's recursive Gram–Schmidt estimation from linear and nonlinear autocorrelations of the series. | Modified Gram–Schmidt least squares on the target rows, re-orthogonalised; collinear terms skipped (64–89). | **matches**. Both give least-squares estimates; BP96 does not specify how the autocorrelation estimates treat the series edges. |
| F5 | Error normalisation | p. 215, Eq. 2; p. 216 ("error standard deviations") | ε² = Σ(ŷ−y)² / Σ(y−ȳ)²: a normalised residual variance, so ε is RMS-like. | ε = √(RSS/SST), with SST about the mean of the target rows (99, 154). C uses 0.5·log(RSS/SST) = log ε (92, 105). | **matches**. Phase 10 "choice 1" is now sourced. |
| F6 | Criterion | p. 216, Eq. 3; PM97 p. 494 (Fig. 2 legend) | C(r) = log ε(r) + r/N. BP96 leaves the log base unstated. PM97 writes log_e, and BP96's own C-values in Figs. 1b and 2a only agree with its F_calc values under a natural log. | `np.log`, natural log (105, 113). | **matches** |
| F7 | N in C(r) | p. 215, Eq. 2 (sums over n = 1..N of the series); p. 217, Fig. 2 legend (df); PM97 p. 494, Fig. 2 legend (N is the length of the series) | N = number of points in the series. | N = number of regression rows (100). All models share rows n ≥ 83 (146), so a 1,000-point series gives N = 917 and a 512-point window gives N = 429. | **deviates, decision-relevant**. The penalty r/N is 9 % (1,000 points) to 19 % (512) larger than in BP96. This affects the d_opt choice at the margin. |
| F8 | Linear search | p. 216; Fig. 1b; Fig. 2a | The best linear model is the κ_lin with minimum C at d = 1. The linear curves span the same r range as the nonlinear ones: about 85 terms in Fig. 1b, 56 in Fig. 2a. | κ_lin ≤ 83, i.e. up to 84 terms, the size of the largest nonlinear model (41, 101–107). | **matches**. Phase 10 "choice 3" is now sourced. |
| F9 | Nonlinear search range | p. 216 ("increasing κ and d > 1"); Fig. 1b (κ = 6, d = 1–3, complete models at r = 7, 28, 84); Fig. 2a (κ = 3, d up to 5, r = 56); Fig. 1c (a κ = 2, d = 8 model) | κ and d are increased, with no fixed limits. The examples use d up to 5 (Fig. 2a) and 8 (Fig. 1c). | κ = 1..6 and d ≤ 3 (35, 40, 108–122): the defaults of the 2006 Poon-lab routine, not BP96. | **not specified**. R0's range covers BP96's Fig. 1 setting (κ = 6, d ≤ 3) but cannot pick d > 3. Recorded as a limitation; left unchanged (the 12 h budget and the frozen R0). |
| F10 | Nonlinear candidates | p. 216 | The nonlinear family is d > 1. | The best nonlinear model is chosen among truncations containing at least one nonlinear term (114). | **matches** |
| F11 | Decision rule | p. 216 | Nonlinear determinism is indicated if d_opt > 1, i.e. the overall C-minimum is nonlinear. An F-test, or Mann–Whitney, corroborates. In the Fig. 1–2 examples all criteria agree ("three criteria"). | Nonlinear iff C_nl < C_lin AND p < α (155, 174). | **matches**: d_opt > 1 plus the corroborating test, required jointly. |
| F12 | F statistic | p. 216, Fig. 1 legend | F_calc = (ε_lin)² / (ε_nl)²: the plain ratio of the normalised error variances. | `F_var` (156–159): (RSS_lin/df_lin) / (RSS_nl/df_nl) = F_calc × df_nl/df_lin. | **deviates, decision-relevant near the threshold**. R0's F is smaller by df_nl/df_lin (e.g. 944/968 = 0.975), so R0 is slightly more conservative. The variance-ratio form is confirmed (Phase 10 "choice 2"). The nested partial F-test (Phase 11 R5) is **not** BP96's test. |
| F13 | F degrees of freedom | p. 217, Fig. 2a: F(0.01, 968, 944) for 1,000 points; Fig. 2b: F(0.06, 1855, 1783) for 1,874 points; p. 216, Fig. 1: F(0.01, 915, 915) | df = N − r for each model, with N = series length. Fig. 2a's optima are r_lin ≈ 32 and r_nl = 56, giving 968 and 944. In Fig. 1 both optima sit near r = 84–85, so 915, 915 is N − r to within one. | df = N − r with N = number of rows (157). | **deviates through F7** (N) |
| F14 | Significance level | p. 216, Fig. 1 legend (F(0.01, …), α = 0.01, Z vs −2.33); p. 217, Fig. 2a; Table 1 legend (99 % cutoff) | 1 %. The decision paragraph itself states no level; the figure legends and Table 1 do. | ALPHA = 0.01 (42). | **matches** |
| F15 | F-test sidedness | p. 216, Fig. 1 legend (F_calc ≫ critical value) | One-sided: the linear error larger than the nonlinear. | `stats.f.sf` upper tail (159). | **matches** |
| F16 | Gaussianity check | p. 216 | The F-test is for models with Gaussian residuals; Gaussianity was checked by a χ² test with a 99 % cutoff. | None. | **deviates, decision-relevant** when residuals are non-Gaussian, as RR residuals with ectopic spikes are. |
| F17 | Mann–Whitney alternative | p. 216; Fig. 1 legend (Z = −37.13 ≪ −2.33); Fig. 2a (Z = −24.60) | A nonparametric rank-sum test that does not need Gaussian residuals; one-sided at α = 0.01. Which residual quantity is ranked is not stated. | Only as an option (Phase 11 R4): one-sided Mann–Whitney on squared residuals, "greater" (167–170). | **sidedness matches; the ranked quantity is not specified**. Not part of R0. |
| F18 | Surrogate step | p. 216 (left and right columns) | Best linear and nonlinear models are also fitted to surrogates with the same autocorrelation and power spectrum. The original's nonlinear error must beat both surrogate errors. Surrogates play "only a confirmatory role". | None. | **deviates (step omitted)**. Decision-relevant if the surrogate comparison would fail where the F-test passes. Tested in Tier 1B (B4). |
| F19 | Preprocessing | p. 217 (Fig. 1 example analysed without reprocessing) | No preprocessing. | Standardises y before fitting (143–144). | **matches in effect**: every leading-term span is invariant to affine changes of y, so decisions are unchanged. |
| F20 | First κ points | p. 215, Eq. 2 | Sums run over n = 1..N. How lags before the series start are handled is not stated. | Common target rows n ≥ 83 for every model (146). | **not specified** |
| F21 | Minimum length | — | Not stated. BP96 examples are about 1,000 points. | At least 3·M = 252 rows (at least 300 for masked series) (149). | **not specified**: an implementation guard. |
| F22 | Noise titration (NL) | absent from BP96 and PM97 | NL comes from Poon & Barahona 2001, not from these papers. | `titrate` (178–220). | **outside the scope of BP96 and PM97**. Phase 11 (B1, R7) showed titration's positive/negative decision equals the noise-free VWK decision. |

## 2. PM97

| # | Point | Source location | What PM97 says (paraphrased) | Relation to this project | Verdict |
|---|---|---|---|---|---|
| P1 | Subjects and data | p. 493; Acknowledgements p. 495 | 8 healthy subjects and 11 CHF patients, aged 22–71. Continuous Holter heartbeat intervals at 128–250 samples/s from "an established database" (refs. 17, 20, 22). Data provided by A. L. Goldberger and J. E. Mietus. | PM97 does not name the databases. nsrdb (128 Hz) and chfdb (250 Hz) fit the stated sampling rates and providers, but that remains an **inference**. | as stated |
| P2 | Medication | p. 493 | None of the subjects took medication affecting heart rate. | Not documented for nsrdb / chfdb. | cannot be matched |
| P3 | Night window | p. 493; Fig. 1 legend (during sleep, 20,000 beats) | Nocturnal hours only, 20,000 beats per subject. | Tier 1C uses the night block, or a registered fallback when there is no clock time. | replicated by rule |
| P4 | Subject selection | p. 493 | Selected for a stable mean heart rate and a limited number of ectopic and undetected beats. | Tier 1C: all subjects (primary) plus a rank-based PM97-like selection (secondary). | replicated by rule |
| P5 | Preprocessing | p. 493 | Linear-regression detrending, otherwise no preprocessing, so ectopic beats stay in. | Phases 10–11 did not detrend and analysed 512-interval windows or 12-min segments. | **design differs from PM97** → Tier 1C |
| P6 | Segments | p. 493; p. 494 (Fig. 2c) | 40 contiguous 500-beat segments; repeated with 2,000-beat segments (their number is not stated). | Tier 1C: 40 × 500 (primary) and 10 × 2,000 (secondary). | replicated |
| P7 | Test | p. 494, Fig. 2 legend | The BP96 test. Selection by the F-ratio test of residuals at 1 %. C(r) = log_e ε + r/N with N = length of the series. The nonlinear model must beat the best linear model "with a similar number of" terms. | Same issues as F6, F7, F12. | see BP96 rows |
| P8 | No titration | pp. 493–494 | "Chaos" means the test fired; no noise titration is used. | Phases 10–11 called a segment positive when NL > 0, which equals the noise-free VWK decision (Phase 11 B1). | equivalent at the decision level |
| P9 | κ, d | — | Not stated. | R0 defaults. | not specified |
| P10 | Results | p. 494; Fig. 2b–c | Healthy: chaos in nearly all segments. CHF: detection markedly lower (by 20–80 %) in all but one patient (P1). The 2,000-beat segments show the same pattern. | Tier 1C compares. | — |
| P11 | Coefficient variability | p. 494; Fig. 3 | Leading nonlinear coefficient estimates (mean ± s.d.) across chaotic segments. Group mean of the s.d. of the "first 80 nonlinear coefficient estimates" is higher in CHF, in all CHF patients but P1. | The model behind the coefficients is not defined. Tier 1C registers one reading. | ambiguous |
| P12 | Diagnostic errors | p. 495 | One CHF patient failed both criteria; type I error 9 %, type II 0 %. | — | — |

**Secondary-source error.** Wu et al. 2009 described its 12-min (about 800-beat) segmentation as following PM97
(Phase 10 METHODS §2). PM97 actually used 500- and 2,000-beat segments.

## 3. Your summary checked against the PDFs

Every numbered point matches the PDFs except for these refinements:
- **Point 3 (C(r)).** BP96 does not state the log base or define N in Eq. 3. PM97's Fig. 2 legend gives log_e and N
  = series length, and BP96's figures are consistent with both.
- **Point 4 (search).** Correct for Fig. 1b. BP96 also varies d beyond 3 (d ≤ 5 in Fig. 2a; a κ = 2, d = 8 model in
  Fig. 1c), so R0's d ≤ 3 comes from the 2006 routine, not BP96.
- **Point 5 (decision).** The 1 % level is in the figure legends and Table 1, not in the decision paragraph. The
  degrees-of-freedom reading N − r is correct and is confirmed by Fig. 2b as well (1,874 points: 1855 and 1783).
- **Point 8 (spikes).** The spike sentence sits in the ecological-model paragraph (p. 217). "Not tested on linear
  signals" is correct in the sense that no such test is reported.
- **Point 9 (Fig. 1).** The Fig. 1d dashed line is the 99 % critical value (1/1.17 = 0.85). The Fig. 1d points cross
  it between 50 % and 60 % noise.
- **PM97 decision.** PM97 compares with the best linear model "with a similar number of" polynomial terms.
- **PM97 data source.** PM97 does not name nsrdb or chfdb (point P1).

## 4. Consequences

Decision-relevant deviations of R0 from BP96: **F7/F13** (N = rows, not series length), **F12** (df-adjusted F, not
the plain ratio), **F16** (no χ² gate, no Mann–Whitney fallback) and **F18** (no surrogate step). R0 remains the
registered primary. Opt-in readings (PREREGISTRATION.md):
- **R0-BP96**: F12 + F7/F13 + F16. The plain ratio F_calc = ε_lin²/ε_nl², with N = series length in C(r) and in the
  df. If a χ² test (99 %) rejects Gaussian residuals of either best model, the one-sided Mann–Whitney test at 1 % is
  used instead of the F-test.
- **R0-BP96F**: F12 + F7/F13 only, to show the χ² gate's separate effect.
- **R0-S-FT** and **R0-S-AAFT**: R0 plus the F18 surrogate requirement.

Not changed (not specified in BP96): the κ and d ranges (F9), the common target rows (F20), the minimum length (F21).
