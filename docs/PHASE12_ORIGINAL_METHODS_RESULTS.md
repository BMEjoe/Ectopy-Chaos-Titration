# Phase 12: the implementation checked against the original papers, and PM97 replicated

Branch `phase12-original-methods`, mirrored to `claude/beautiful-hawking-i1ommp`. The preregistration
(`experiments/phase12_original/PREREGISTRATION.md`) was pushed in commit `8d681ec`, before any Phase 12 run on nsr2db,
chf2db, nsrdb or chfdb; there are no amendments.

Every number below is read from `experiments/phase12_original/results/conf/analysis/analysis12.json`. Ready-to-paste
numbers with source keys are in `docs/PHASE12_NUMBERS.md` (IDs such as **R13** refer to it). Unless stated otherwise:
CIs are 95 % subject-cluster bootstrap CIs (10,000 resamples), rates over seeds use Wilson CIs, and pp = percentage
points.

## 1. Bottom line

**Is the implementation faithful?** Yes, with four listed deviations.
- R0 matches BP96 on the model, term order, error normalisation, criterion (natural log), linear search range,
  decision rule, significance level and test sidedness.
- It deviates on four decision-relevant points:
  1. N in the criterion and the degrees of freedom (rows, not series length);
  2. the F statistic (adjusted for degrees of freedom, not the plain ratio of error variances);
  3. no χ² Gaussianity check, so no Mann–Whitney fallback;
  4. no surrogate step.
- **Benchmarks.** R0 reproduces BP96's Santa Fe laser result closely (F_calc 10.51 vs 10.23). It passes the registered
  ecological-model benchmark, but only at the rule's lower edge: detection is lost at 35 % correlated noise, against
  about 50 % in BP96.
- **Faithful reading.** R0-BP96, which removes deviations 1–3, **fails** that benchmark (lost at 30 %).
- **Table 1 (approximate).** 4 of 9 systems fall outside ±15 pp under every reading.

**Do the headline claims survive under the faithful reading?** Yes for the ectopy claims; no for one CHF claim when
it is generalised.
- **Every reading passes V1–V6** (R0-BP96F, R0-BP96, R0-S-FT, R0-S-AAFT).
- Under every reading, the burden effect, the absence of a group difference after ectopy handling, couplets flagged
  while linear RR is not, and the single-beat transplant effect all hold. **No registered "against H" result.**
- **The single-beat effect is smaller under BP96's own test conventions:** +27.8 pp (R0-BP96) against +64.2 (R0).
- **Isolated 5 % ectopy is flagged far more often under BP96's conventions:** 39–46/100 against R0's 4/100.
- **The raw CHF > NSR difference disappears** with amplitude-adjusted surrogates: +4.2 pp [−5.0, +13.8].
- **The surrogate step does not remove the N5 counterexample** (100/100 under every reading). R0, like titration as
  published in 2001, flags static transforms, and so does BP96's surrogate step.

**PM97 replication (nsrdb / chfdb): the registered "against my claims" criterion is met.**
- **PM97's raw pattern is not reproduced.**
  - Healthy subjects are positive in 52.6 % of segments, not the ~100 % PM97 reported.
  - CHF is not lower: +16.5 pp [−1.9, +33.3].
- **After ectopy handling, CHF is clearly lower:**
  - masked −45.6 pp [−57.0, −33.7];
  - edited −26.4 [−42.8, −8.2];
  - Wu −27.9 [−44.1, −10.3].
- **In ectopy-free segments** (exploratory, post hoc), CHF is also lower: −34.4 [−56.3, −6.1].
- **On these databases, a CHF deficit in VWK detection appears once ectopy is removed.** On nsr2db / chf2db (Phases
  10–11) it does not.
- **Consequence for the manuscript.** It must not say that PM97's CHF decrease is an ectopy artefact in general. What
  the data support is narrower:
  - PM97's raw procedure, which leaves ectopy in, does not reproduce its own pattern on these databases.
  - Ectopy drives detection within subjects (OR 2.35 [1.78, 3.10]).
  - Whether the remaining CHF deficit reflects less chaos, lower variability or data loss is not resolved here (§5).

## 2. Summary

| Analysis | Registered rule | Result (CI) | Verdict |
|---|---|---|---|
| Sanity: analysis11 rerun | identical | byte-identical JSON | pass |
| R0 recomputed vs stored | — | 988/988 windows (raw, masked, edited) and 1,756/1,756 transplant decisions agree | pass |
| Tier 1A fidelity | — | 22 points: 12 match, 1 partly matches (Mann–Whitney), 5 rows deviate (4 distinct decision-relevant deviations; F13 follows from F7), 3 not specified, 1 outside scope (titration) | faithful with listed deviations |
| B1 ecological model, R0 | d_opt > 1 ≥ 19/20; median F_calc ≥ 100; loss point 35–65 % | 20/20; 6,885 [6,327, 8,067]; 35 % [30, 40] | **pass (at the lower edge)** |
| B1, R0-BP96 | same | 20/20; 6,885; **30 % [25, 31]** | **FAIL** |
| B1, R0-BP96F | same | 20/20; 6,885; 35 % | pass |
| B2 Table 1 (approximate) | median within ±15 pp | 4 of 9 miss under every reading: ecological 34 (50), GOPY 42 (25), Rössler 25 (75), Lorenz 72 (50) | **4 misses** |
| B2 Mackey–Glass τ = 100 | same | never detected noise-free (0/20); "within 15" only because Table 1 says 15 | not reproduced |
| B3 Santa Fe A (R0) | — | F_calc 10.51 (BP96 10.23); Z −22.1 (−24.6); F crit 1.17 (1.16). Series E not obtainable | close agreement |
| B4 V1–V6 | all readings | R0-BP96F, R0-BP96, R0-S-FT, R0-S-AAFT all pass | pass |
| B4 static transforms, R0 | flags if any Wilson lower bound > 5 % | flags cube and exp (φ = 0.5 and 0.9) and AR(1)-N5 (φ = 0.9) | **R0 flags static transforms** |
| B4 surrogate removal | ≤ 5/100 where R0 ≥ 50/100 | FT removes none of 4 cells (90–95/100 remain); AAFT 1 of 4 (4–11/100 remain); Phase 10 N5 100/100 under all | **not removed** |
| B4 claims (i)–(v), 5 readings | against H if a verified reading fails (i), (iii) or (v) | (i), (iii), (iv), (v) hold under all; (ii) fails under R0-S-AAFT (+4.2 [−5.0, +13.8]) | **no against-H result**; (ii) reading-dependent |
| B4 single-beat transplant V1 − V0 | (v) CI > 0 | R0 +64.2 [59.2, 69.3]; R0-BP96F +47.4; R0-BP96 +27.8 [21.9, 33.6]; R0-S-FT +61.0; R0-S-AAFT +29.8 [24.1, 35.6] | holds; **size depends on reading** |
| PM97 raw pattern (500-beat, R0, all 33 subjects) | NSR ≥ 90 %; CHF − NSR CI < 0; ≥ 2/3 of CHF 20–80 pp below | NSR 52.6 %; +16.5 [−1.9, +33.3]; 2/15 | **not reproduced** |
| PM97 raw, R0-BP96 | same | NSR 74.3 %; −9.6 [−27.0, +7.0]; 6/15 | not reproduced |
| PM97 after ectopy handling | against my claims if masked, edited and Wu CIs are all < 0 | −45.6 / −26.4 / −27.9 pp, all CIs < 0 | **AGAINST my claims (met)** |
| PM97 within-subject ectopy | — | conditional logit OR 2.35 [1.78, 3.10] per doubling | ectopy drives detection |
| PM97 Fig. 3 reading | descriptive | coefficient SD: CHF 5.0 ± 2.9 vs NSR 19.6 ± 9.2 (s.e.); difference −14.6 [−34.1, +1.8] | opposite direction to PM97, not significant |

## 3. Tier 1A: fidelity (`experiments/phase12_original/FIDELITY.md`, Table 12-1)

**Confirmed against the PDFs:**
- **Model (BP96 Eq. 1, p. 215).** Term order, the constant term and M = (κ+d)!/(d!κ!).
- **Error (Eq. 2).** ε² is the normalised residual variance; BP96 calls ε the error standard deviation.
- **Criterion.** C(r) = log ε + r/N (Eq. 3, p. 216). The log is natural: PM97's Fig. 2 legend writes log_e, and BP96's
  C-values agree with its own F_calc.
- **Search.** The linear search spans the same term range as the nonlinear one (Fig. 1b about 85 terms; Fig. 2a 56).
- **Decision.** d_opt > 1, corroborated by the F-test or Mann–Whitney, one-sided at 1 % (Fig. 1–2 legends, Table 1).
- **The three Phase 10 implementation choices are now sourced:**
  - ε as RMS: BP96 Eq. 2.
  - The variance-ratio F-test (not the nested partial F-test of Phase 11 R5): BP96 Fig. 1 legend.
  - A linear memory up to the size of the largest nonlinear model: Fig. 1b, 2a.

**Deviations:**
- **N (F7, F13).** BP96's Fig. 2 df (968 and 944 for 1,000 points; 1855 and 1783 for 1,874) and PM97's Fig. 2 legend
  set N to the series length. The code uses the number of regression rows.
- **F statistic (F12).** BP96 uses the plain ratio ε_lin²/ε_nl². The code divides each RSS by its df, which shrinks F
  by df_nl/df_lin.
- **Gaussianity check (F16).** BP96 checks Gaussian residuals with χ² at 99 % and offers Mann–Whitney as the
  alternative. The code has neither.
- **Surrogates (F18).** BP96 fits surrogates with the same autocorrelation and spectrum, but in a confirmatory role
  only. The code omits them.

**Not specified by BP96 (left unchanged):**
- **The κ and d ranges.** BP96 uses d up to 5 (Fig. 2a) and a κ = 2, d = 8 model (Fig. 1c), so R0's d ≤ 3 comes from
  the 2006 routine, not BP96.
- **The target rows** and **the minimum length**.

**Corrections to the task's summary** (FIDELITY.md §3):
- **1 % level.** It appears in the figure legends and Table 1, not in BP96's decision paragraph.
- **PM97's databases.** PM97 does not name them; nsrdb / chfdb remain an inference.
- **PM97's comparison.** PM97 compares with the best linear model of a similar number of terms.
- **Secondary-source error.** Wu et al. 2009 attributed its 12-min segmentation to PM97, but PM97 used 500- and
  2,000-beat segments.

## 4. Tier 1B: tests against the original papers' numbers

**B1. Ecological model, BP96 Fig. 1** (Table 12-2, F12-1; 20 seeds, 1,000 points)
- **Noise-free.** R0 selects κ = 6, d = 3 in all 20 seeds, with d_opt > 1 in 20/20. Median F_calc is 6,885 [6,327,
  8,067] against BP96's 9,387, the critical value 1.175 against 1.17, and Mann–Whitney Z −36.1 against −37.13.
- **With noise.** Detection is lost at 35 % [30, 40] correlated noise; BP96 reports about 50 %. 1/F_calc rises like
  BP96's Fig. 1d up to about 20 % noise, then faster, crossing the 1 % line between 30 % and 35 %.
- **Rule.** R0 passes at the lower edge of the registered window. R0-BP96 fails (30 % [25, 31]), because its
  Mann–Whitney fallback loses power once the noisy residuals are judged non-Gaussian. R0-BP96F passes.
- The difference from BP96's ~50 % is reported, not explained. Candidates include how BP96 generated its coloured
  noise and its unreported search range.

**B2. Table 1** (Table 12-3; approximate: BP96 gives no parameters or sampling)
- Medians within ±15 pp under R0: logistic 61 (70), Hénon 72 (70), Ikeda 65 (65) and Duffing 44 (40).
- Misses under every reading:
  - ecological model 34 (50);
  - Grebogi's non-chaotic attractor 42 (25), i.e. *more* noise-tolerant than in BP96;
  - Rössler 25 (75);
  - Lorenz 72 (50).
- **Mackey–Glass τ = 100** (chosen because the scan describes an implicit dimension of seven) is **never detected**,
  not even noise-free. With τ = 17 (secondary), 42 % (15).
- Sampling choices plausibly drive the flow systems (BP96 notes the optimal delay matters). These are reported as
  approximate disagreements, not as evidence about the code.

**B3. Santa Fe series A** (1,000 points; public GitHub copy, not committed)
- R0: F_calc 10.51 against BP96's 10.23; F critical value 1.17 (1.16); Z −22.1 (−24.6).
- Best linear model: 30 terms (BP96 about 32). Best nonlinear: κ = 6, d = 3 truncated to 74 terms (BP96: κ = 3,
  d = 5, 56 terms).
- Series E: not obtainable.

**B4. Static transforms and the surrogate step** (Tables 12-4 to 12-7, F12-2)
- **R0 flags static transforms of linear Gaussian data.** At φ = 0.9 (n = 1,000): x³ 97/100, exp 100/100, AR(1)-N5
  25/100. At φ = 0.5: 11–33/100. sign·√|x| and tanh are rarely flagged under R0.
  - Under R0-BP96, tanh is flagged 99/100 (φ = 0.9, n = 1,000): the χ² gate sends these saturated residuals to the
    Mann–Whitney test.
- **The surrogate step:**
  - **FT surrogates** (BP96's description) keep 90–95/100 of x³ and exp positives.
  - **AAFT surrogates** reduce them to 4–11/100, meeting the registered removal threshold (≤ 5/100) in only 1 of 4
    cells.
  - **Phase 10's N5 condition** (warped linear RR) stays **100/100 under all five readings**. A candidate reason, not
    tested, is the known spectral bias of AAFT surrogates for strongly correlated series.
- **Ectopy-driven positives are not removed:**
  - couplets 100/100 and runs 91/100 under every reading;
  - transplant V1 − V0 under R0-S-FT +61.0 pp, under R0-S-AAFT +29.8 [24.1, 35.6].
- **Claims matrix** (Table 12-6). Every reading passes V1–V6. Claims (i) burden OR > 1, (iii) no post-handling group
  difference, (iv) couplets vs N1 and (v) the single-beat effect hold under all five. Claim (ii), raw CHF > NSR in
  windows, holds under R0, R0-BP96F, R0-BP96 and R0-S-FT, but not under R0-S-AAFT (+4.2 pp [−5.0, +13.8]).
- **Isolated ectopy (S2, 5 %, N1)** is flagged 4/100 under R0 but 39/100 under R0-BP96F and 46/100 under R0-BP96.
  Phase 10's "isolated synthetic ectopy rarely triggers titration" holds only for R0's conventions.

**Deviations from the preregistration in Tier 1B:** none.
- The B3 data source (a public GitHub repository, read-only, outside this repository) is recorded in HANDOFF D6.
- The B2 smoke test was disclosed in the preregistration.

## 5. Tier 1C: PM97 replication (Tables 12-8 to 12-11, F12-3)

**Data.**
- 33 subjects: nsrdb 18 and chfdb 15. All have clock start times. 29 blocks start at midnight; 4 use the registered
  lowest-heart-rate fallback because fewer than 20,000 beats followed midnight.
- Beats used: Table 12-8.
- Preprocessing as PM97: linear detrending only, ectopy left in.

**Raw (PM97's procedure), 500-beat segments, R0:**
- Healthy 52.6 % of segments (median 53.8 %); CHF 69.2 % (median 77.5 %).
- CHF − NSR +16.5 pp [−1.9, +33.3]; Mann–Whitney p = 0.034.
- Only 2 of 15 CHF subjects are 20–80 pp below the healthy mean.
- **None of the three registered pattern components holds.**
- Under R0-BP96: healthy 74.3 %, CHF 64.7 %, −9.6 [−27.0, +7.0], 6/15. Not reproduced.
- 2,000-beat segments: healthy 81.7 %, CHF 88.0 %, +6.3 [−9.8, +20.3].
- The PM97-like subset (8 + 11 subjects, listed in Table 12-8): +8.0 [−17.4, +31.4].

**After ectopy handling, R0, 500-beat:**

| Arm | Healthy | CHF | CHF − NSR |
|---|---|---|---|
| masked | 53.3 % | 7.7 % | −45.6 [−57.0, −33.7] |
| edited | 53.3 % | 27.0 % | −26.4 [−42.8, −8.2] |
| Wu | 53.2 % | 25.3 % | −27.9 [−44.1, −10.3] |

- The 2,000-beat segments, the R0-BP96 reading and the PM97-like subset give the same direction.
- **The registered criterion "PM97's decrease survives ectopy handling" is met. This counts against my claims.**
- Within subjects, ectopy raises detection: conditional logit OR 2.35 [1.78, 3.10] per doubling of a segment's
  ectopic beats (32 informative subjects).

**Labelled descriptive and post-hoc analyses** (added after the outputs were seen):
- **Ectopy is far more common in CHF.** 80 % of CHF segments contain ectopy (median 5 ectopic beats per segment),
  against 3 % of healthy segments.
- **Masking costs CHF most of its data.** Only 34 % of CHF segments remain analysable after masking (98 % in healthy
  subjects).
- **Exploratory:** in ectopy-free segments (no non-normal beat, no out-of-range interval; 18 healthy and 8 CHF
  subjects), CHF is still lower: 19 % against 53 %, −34.4 pp [−56.3, −6.1]. Under R0-BP96: −54.1.

**Fig. 3 reading** (registered: full κ = 7, d = 3 model, first 80 nonlinear coefficients): the variability is
*lower* in CHF (5.0 ± 2.9 against 19.6 ± 9.2), the opposite direction to PM97, and not significant (p = 0.58).

**Titration NL** (secondary): median NL among positive raw segments is 0.17 in healthy subjects and 0.60 in CHF. The
titration decisions agree with the noise-free decisions in 1,320/1,320 segments.

**Candidate reasons the raw pattern differs from PM97** (none is claimed as the cause):
- subject selection (PM97 selected for few ectopic beats; the rank-based selection does not reproduce it either);
- the night window: blocks start at midnight here, PM97's exact window is unknown;
- detrending of the full 20,000-beat block versus per segment;
- chfdb severity;
- PM97's unreported κ and d;
- the implementation readings (R0-BP96 moves raw healthy detection to 74 % and reverses the sign of the group
  difference).

**The CHF deficit after handling has not been shown to be "decreased chaos".** Candidates:
- lower variability;
- the much larger data loss in CHF (Phase 11 A1 showed masking alone removes positives);
- the exploratory ectopy-free comparison, which comes from 8 CHF subjects with few ectopic beats.

**Deviations from the preregistration in Tier 1C:** the descriptive and exploratory additions above, labelled in
the JSON and in this report. Nothing registered was changed.

## 6. What this changes in the manuscript

The manuscript draft was not provided, so no `manuscript_check.md` was produced. Passages are identified by content;
numbers are in `docs/PHASE12_NUMBERS.md`.

**Methods**

- *Replace the sentence saying BP96 and PM97 were not accessible:*

  > "The implementation was checked line by line against Barahona and Poon (1996) and Poon and Merrill (1997) (Table
  > S12-1). It matches the published model, error definition, criterion and decision rule, and differs in four
  > decision-relevant details: the penalty and degrees of freedom use the number of regression rows rather than the
  > series length, the F statistic is adjusted for degrees of freedom, there is no χ² Gaussianity check with
  > Mann–Whitney fallback, and the confirmatory surrogate step is omitted (M1)."

- *Implementation choices, now sourced:*

  > "Following BP96 (Eq. 2), ε is the normalised RMS one-step error; the criterion is C(r) = ln ε(r) + r/N (BP96
  > Eq. 3; natural log as in PM97); the test is the one-sided variance-ratio F-test of the best linear against the
  > best nonlinear model at 1 % (BP96 Fig. 1 legend), not a nested F-test; and linear models may use as many terms as
  > the largest nonlinear model, as in BP96 Fig. 1b."

- *Add:*

  > "A titration-positive result is identical to a positive noise-free BP96 test: the noise schedule and seed only
  > change the size of NL (M3; Phase 11)."

- *Add the readings:*

  > "We repeated every headline analysis under the original paper's own conventions (R0-BP96) and with BP96's
  > surrogate step (phase-randomised or amplitude-adjusted surrogates); all readings passed the standard-system
  > verification (B6)."

**Results**

- *Single beat:* keep R12 and add R13:

  > "Under the original test's conventions the single-beat effect was smaller but still clear (+27.8 pp [+21.9,
  > +33.6])."

- *Benchmarks:*

  > "Our implementation reproduced BP96's Santa Fe laser result (F_calc 10.51 vs 10.23; B5) and detected nonlinearity
  > in BP96's ecological model in all 20 seeds (B1), but lost detection at 35 % correlated noise, against about 50 % in
  > BP96 (B2). Against BP96's approximate Table 1 it agreed within 15 points for 5 of 9 systems (B4)."

- *Isolated ectopy:* qualify the Phase 10 statement with S6:

  > "Isolated ectopy at 5 % density was rarely flagged under our conventions (4/100) but often under BP96's (39–46/100)."

**N5 counterexample (Results or Discussion)**

- *Registered sentence (S3):*

  > "Our implementation, like titration as published in 2001, flags static monotone transforms of linear data (N5:
  > 100/100), and so does BP96's surrogate step: neither phase-randomised nor amplitude-adjusted surrogates removed
  > the N5 positives (S1)."

- Add that BP96 reported such transforms as not flagged. Our implementation does flag them; amplitude-adjusted
  surrogates reduce but do not remove the AR(1) cases (S2).

**Introduction and Discussion: PM97**

- *Replace any description of PM97's data or test with:*

  > "Poon and Merrill analysed 20,000 nocturnal beats from 8 healthy subjects and 11 CHF patients, linearly detrended
  > but otherwise unprocessed (ectopic beats retained), in 40 segments of 500 beats, using the BP96 test at 1 %
  > without noise titration; they reported detection in nearly all healthy segments and a 20–80 % lower rate in all
  > but one CHF patient."

- *Replication, required by the integrity rule:*

  > "On nsrdb and chfdb, the databases whose sampling rates and providers match PM97's description, PM97's procedure
  > did not reproduce its pattern: healthy segments were positive in 52.6 % and CHF segments were not less often
  > positive (+16.5 pp [−1.9, +33.3]; P1–P3). Ectopy raised detection within subjects (OR 2.35 [1.78, 3.10]; P8). After
  > ectopy was masked, edited or removed, however, CHF segments were positive markedly less often than healthy ones
  > (−26 to −46 pp; P5–P7), also in ectopy-free segments (exploratory; P9). On nsr2db and chf2db no such difference
  > remained after handling (R4–R6). Whether the CHF deficit on chfdb reflects decreased chaos, lower variability or
  > the larger data loss of handling (P10) cannot be decided here."

- **Weaken the claim** "the decrease of cardiac chaos in CHF is not supported." Limit it to nsr2db / chf2db and to
  raw data. Say explicitly that on chfdb, after ectopy handling, a CHF deficit in detection is present.

**K3 and robustness qualifiers (Phase 11, unchanged)**
- K1: "0/10,100" holds at 1/360 s only; at 1/128 s, 55/10,092.
- K2: the K3 counts are seed-sensitive (8, 11, 10).
- K4: K1's burden association is not significant within subjects.
- The sham control (R14), and claim (ii) being reading-dependent (S5), belong in the limitations.

**Supplementary items**
- Tables 12-1 to 12-11: `experiments/phase12_original/results/conf/analysis/paper_tables12/`.
- Figures F12-1 to F12-3: `…/figures12/`.
- FIDELITY.md as a supplementary note, with paraphrases only.

## 7. Regeneration

Run from the repository root on branch `phase12-original-methods`, with the environment from
`experiments/phase10_final/requirements-phase10.txt` and `OMP_NUM_THREADS=1`. The papers are expected in the
git-ignored `private_inputs/`. The runner refuses confirmatory data unless the preregistration is committed and
pushed.

```bash
PYTHONPATH=. python -m experiments.phase10_final.run10 --phase conf --download
```

```bash
PYTHONPATH=. python -m experiments.phase10_final.data10 --download dev
```

```bash
PYTHONPATH=. python -m experiments.phase11_robustness.data11 --verify
```

```bash
PYTHONPATH=. python -m pytest -q experiments/phase12_original/tests
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.run12 --pilot --part all --workers 4
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.budget12
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.run12 --part tier1b --workers 4
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.run12 --part tier1c --workers 4
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.b3_santafe /path/to/SF_A.dat
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.analysis12
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.paper_tables12
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.figures12
```

```bash
PYTHONPATH=. python -m experiments.phase12_original.numbers12
```

`SF_A.dat` comes from github.com/shubhomoydas/ad_examples, commit 4f20369, `datasets/SantaFeCompetition/`. Total
runtime with 4 workers is about 30 minutes.
