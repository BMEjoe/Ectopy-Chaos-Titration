# Phase 11 methods, sources and implementation notes

Phase 11 adds registered secondary robustness analyses (PREREGISTRATION.md). No Phase 10 file, detector, threshold
or titration default is changed; `final_pipeline.py` is unchanged (0 diff lines), so the ground-rule check is not
required.

## 1. Code
| file | role |
|---|---|
| `titration11.py` | titration variant readings R0-R7 (reuses `titration10` functions; identical at R0: `tests/test_titration11.py`, 20 seeded series, raw and masked) |
| `hrt11.py` | heart rate turbulence measurement (TO, TS), HRT profiles and single-event insertion (`tests/test_hrt11.py`) |
| `harness11.py` | sham / donor / extended masks, NN covariates, K3 with degeneracy diagnostic (identical decision to the frozen K3: `tests/test_harness11.py`), variance-preserving spike-in |
| `run11.py` | guarded, resumable runner; `--pilot` on development data |
| `analysis11.py` | all Phase 11 statistics |
| `data11.py` | nstdb download + SHA-256 verification of all data |
| `paper_tables.py`, `figures11.py` | Tier 3 |

## 2. Heart rate turbulence (A2.1-A2.2)

### 2.1 Sources (accessed 2026-11-03)
| key | source | access | used for |
|---|---|---|---|
| [Sch99] | Schmidt G et al. Heart-rate turbulence after ventricular premature beats as a predictor of mortality after acute myocardial infarction. Lancet 353:1390 (1999) | **not accessed** (paywall); its definitions are taken from the restatements below | original TO / TS definitions; abnormal cut-offs TO >= 0 %, TS <= 2.5 ms/RR |
| [Bau08] | Bauer A et al. Heart rate turbulence: standards of measurement, physiological interpretation, and clinical use. ISHNE consensus. J Am Coll Cardiol 52:1353-1365 (2008) | **full text not accessible** (sciencedirect returned 403); abstract and bibliographic record only | consensus definitions (the user's specification) |
| [Wat03] | Watanabe MA. Heart rate turbulence: a review. Indian Pacing Electrophysiol J 3:10-22 (2003), PMC1555628 | open (PMC) | TO formula with the compensatory pause as beat 0 (two sinus intervals before the coupling interval vs the first two after the pause); TS = maximum slope of 5-interval regressions within the first 20 sinus intervals of the averaged tachogram; recent Schmidt-lab studies require >= 3-5 sinus intervals before the VPC and relaxed the post-pause minimum from 20 to 15; normal TO < 0 and TS > 2.5 |
| [Cyg03] | Cygankiewicz I et al. Ann Noninvasive Electrocardiol 8:289 (2003), PMC7313232 | open (PMC) | the filter rules "adapted from the web page popularizing the noncommercial use of HRT" (Schmidt group): exclude sinus intervals < 300 ms or > 2000 ms, > 200 ms difference to the preceding sinus interval, > 20 % difference to the reference = mean of the 5 last sinus intervals before the VPB; prematurity >= 20 %; post-extrasystolic interval >= 20 % longer than the normal interval; TO per VPB then averaged; TS on the averaged local tachogram |

### 2.2 Implementation (`hrt11.tachograms`, `subject_hrt`)
- Usable isolated VPB: symbol V at beat v; beats v-6 .. v-1 and v+1 .. v+21 labelled normal (Phase 10 label rule);
  5 sinus intervals before the coupling interval, 20 after the compensatory pause; filter rules of [Cyg03] applied to all
  25 sinus intervals (the step rule within the pre- and within the post-sequence).
- TO per VPB with RR-2, RR-1 = the two sinus intervals before the coupling interval and RR1, RR2 the first two after the
  compensatory pause; subject TO = mean over VPBs ([Cyg03], [Bau08]).
- TS on the subject's averaged tachogram: maximum least-squares slope over 5 consecutive intervals within RR1..RR15
  ([Bau08] as specified by the user; [Wat03] reports the Schmidt-lab relaxation to 15) and, secondary, within RR1..RR20
  ([Sch99] as restated by [Wat03]).
- Normalised tachogram h_j = RR_j / mean(RR-2, RR-1) - 1.

### 2.3 Details that could NOT be verified against [Sch99] / [Bau08] (both inaccessible)
1. Whether the 20 % reference rule and the 200 ms step rule apply to the post-pause intervals beyond RR15 (applied here
   to all 20, which is stricter).
2. The exact number of pre-VPB sinus intervals required by the consensus (5 here, from the reference definition of
   [Cyg03]; [Wat03] reports 3 or 5 in Schmidt-lab studies).
3. The minimum number of VPBs for a reliable subject value (>= 5 used for the primary per-subject report, >= 1 secondary;
   no source read here states a consensus minimum).
4. Whether TO in [Bau08] is computed per VPB and averaged (as [Cyg03] states) or on the averaged tachogram (as [Wat03]
   describes for the original 1999 method). Per-VPB averaging is used.
5. The reference values "TS ~ 13 ms/RR, TO ~ -2.5 % (healthy medians, Grimm et al. 2003)" are used as given in the task;
   Grimm et al. 2003 was not read. [Wat03] cites TS 26 +/- 2 (mean +/- sem) in 10 normal subjects (a different study).
6. Annotation resolution: nsr2db / chf2db beat times are at 1/128 s (7.8 ms); TS values below ~2 ms/RR are at the
   quantisation scale of a single interval, although averaging over VPBs reduces it.

### 2.4 HRT generator (`hrt11.insert_event`, `parametric_profile`)
- Profile h_j multiplies the 20 intervals after the compensatory interval (j = 1..20).
- HRT-normal / HRT-blunted: h = -a exp(-(j-1)/1.5) + b exp(-((j-7)/4)^2), (a, b) >= 0 solved by least squares so that
  TO and TS measured by `hrt11`'s own code on a constant baseline equal to the event's pre-VPB interval
  mean(x_{j-2}, x_{j-1}) hit (-2.5 %, 13 ms/RR) and (0 %, 2.5 ms/RR). The profile is recomputed for every insertion,
  so TS is in ms per RR interval at the window's own heart rate. Tests: targets met within 0.05 at baselines
  0.6-1.2 s, and re-measured on a generated tachogram with `tachograms()`.
- HRT-empirical: the group's subject-weighted mean normalised tachogram (A2.1), applied as measured.
- Ventricular, couplet and atrial events: the Phase 6 model (coupling factors U(0.6, 0.8); full compensatory pause for
  ventricular beats; atrial post-ectopic factor U(1.0, 1.1)).

## 3. Sham masks (A1)
- Circular shift of the true K3 interval mask (`np.roll`), so the sham has the true mask's run structure and count.
- Forbidden zone: the true ectopy-related intervals and intervals b .. b + 20 after every abnormal beat b (interval b
  runs from beat b to b + 1; the K3 mask already contains b - 1 and b), i.e. the 20 post-ectopic intervals where HRT
  occurs.
- Titration seed of a sham arm = the Phase 10 seed of the corresponding true arm (data CRC, salt 0), so true and sham
  arms use the same noise realisations.

## 4. Titration readings (A5)
- Barahona & Poon 1996 (Nature 381:215) was **not available** (no file in `private_inputs/`): Phase 10 METHODS.md 2.3
  choices 1-3 remain **not addressed** by the primary source. R8 ("BP96 as published") is not run.
- R1-R3 change only the linear memory limit; the masked arm's row span follows the longest candidate lag
  (max(6, kappa_lin)).
- R7 uses `titration10.titrate(scheme='increment')`'s logic (steps of 1 % of sigma_y up to 16 sigma_y).

## 5. Covariates (A3)
- NN = intervals not in the K3 mask. RMSSD uses only successive differences of two adjacent NN intervals. Half shift
  splits the window by interval index. The K1 detrend ratio is `final_pipeline.linear_trend_ratio` (|fitted linear
  change| / residual SD; K1's gate is 0.7) on the NN series.
- Elapsed hours since the record's first beat (no header start times: clock time is unavailable).

## 6. K3 degeneracy diagnostic (B2)
`harness11.k3_diag` repeats `final_pipeline.masked_growth_chaos_test` step by step (same RNG stream, same calls) and
also returns the surrogate prediction-error mean and SD per embedding dimension. Degenerate iff the SD at the m attaining
z_max is < 1e-6 (the frozen code floors the SD at 1e-12, which produces z ~ 1e11 when every surrogate gives the same
error).
