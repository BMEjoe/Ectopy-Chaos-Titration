# Claims check (EXPLORATORY, post hoc)

One row per quantitative claim in the abstract, lead paragraph and conclusions (plus the Results statements they rest on). Values per reading; 'holds' applies the stated criterion to every reading listed. Phase 12 values are stored results on that phase's seeds. R0-ORIG-minus-S3b (= S1 + S2 + S5a, the configuration closest to the published benchmarks) was added by PLAN.md Amendment 3 after the planned runs.

## A0 (abstract; conclusions (C8 'most often caused by ectopic beats'))
**Claim.** In heart rate, titration mainly detects ectopic beats

- R0: 52 [43, 62] %
- R0-BP96: 28 [19, 36] %
- R0-ORIG: 36 [27, 46] %
- R0-ORIG-minus-S3b: 24 [17, 32] %
- R0-BP96F (Phase 12): 26 [18, 34] %
- R0-S-FT (Phase 12): 64 [54, 73] %
- R0-S-AAFT (Phase 12): 41 [28, 54] %
- **Criterion:** attributable fraction AF > 0.50 under every reading (PLAN.md 5). **Holds under every reading:** no
- **Suggested wording:** In heart rate, ectopic beats are a major source of titration positives: they account for about 24-64 % of positive windows, depending on the test's conventions, and burden, not group, predicts a positive window

## L1 (lead paragraph)
**Claim.** One ectopic beat nearly tripled the chance of titration reporting chaos

- R0: 2.88x
- R0-BP96: 1.54x
- R0-ORIG: 1.59x
- R0-ORIG-minus-S3b: 1.59x
- R0-BP96F (Phase 12): 1.59x
- R0-S-FT (Phase 12): 3.62x
- R0-S-AAFT (Phase 12): 2.17x
- **Criterion:** rate ratio (burden 1 vs 0) >= 2.5 under every reading. **Holds under every reading:** no
- **Suggested wording:** One ectopic beat raised the chance of titration reporting chaos 1.5- to 3.6-fold, depending on the test's conventions

## V-B (Sec. V B; Discussion)
**Claim.** With burden in the model there was no separate effect of CHF (OR 0.77); ectopy accounts for the group difference at the window level

- R0: 0.77 [0.46, 1.28]
- R0-BP96: 0.85 [0.55, 1.31]
- R0-ORIG: 0.67 [0.40, 1.14]
- R0-ORIG-minus-S3b: 0.91 [0.52, 1.59]
- **Criterion:** adjusted CHF OR CI includes 1 under every reading. **Holds under every reading:** yes

## IV-A (Sec. IV A; Table III)
**Claim.** Isolated beats (5 %, 10 %) and atrial-type beats rarely triggered titration

- R0: isolated 5 %: 5/100; 10 %: 0/100; atrial 5 %: 2/100
- R0-BP96: isolated 5 %: 57/100; 10 %: 0/100; atrial 5 %: 88/100
- R0-ORIG: isolated 5 %: 82/100; 10 %: 9/100; atrial 5 %: 90/100
- R0-ORIG-minus-S3b: isolated 5 %: 63/100; 10 %: 10/100; atrial 5 %: 26/100
- **Criterion:** isolated 5 % and atrial 5 % <= 10/100 under every reading. **Holds under every reading:** no
- **Suggested wording:** Under my implementation, isolated beats at 5-10 % density rarely triggered titration; under the original test's conventions, 5 % isolated and atrial-type beats made most windows positive

## A1 (abstract)
**Claim.** For all 100 linear stochastic RR series with ectopic couplets, every titration version tested positive

- R0: 100/100
- R0-BP96: 100/100
- R0-ORIG: 100/100
- R0-ORIG-minus-S3b: 100/100
- R0-BP96F (Phase 12, other seeds): 100/100
- R0-S-FT (Phase 12, other seeds): 100/100
- R0-S-AAFT (Phase 12, other seeds): 100/100
- **Criterion:** couplets 5 %: 100/100 under every reading. **Holds under every reading:** yes

## A2 (abstract; lead)
**Claim.** In 988 real windows, positive rates increased from 22 % (no ectopic beats) to 63 % (1 ectopic beat)

- R0: 22 % -> 63 %
- R0-BP96: 37 % -> 57 %
- R0-ORIG: 31 % -> 50 %
- R0-ORIG-minus-S3b: 50 % -> 79 %
- **Criterion:** rate at burden 1 above burden 0 (Wilson CIs disjoint); the quoted numbers are R0's. **Holds under every reading:** yes
- **Suggested wording:** positive rates rose from 22-50 % (no ectopic beats) to 50-79 % (1 ectopic beat), depending on the test's conventions

## A3 (abstract; conclusions)
**Claim.** One premature beat inserted into clean windows added 28-64 percentage points, depending on conventions

- R0: 64.2 [59.2, 69.3] pp
- R0-BP96: 27.8 [21.9, 33.6] pp
- R0-ORIG: 33.9 [28.7, 39.3] pp
- R0-ORIG-minus-S3b: 45.6 [39.1, 52.2] pp
- R0-BP96F (Phase 12): 47.4 [40.9, 54.0] pp
- R0-S-FT (Phase 12): 61.0 [56.0, 66.2] pp
- R0-S-AAFT (Phase 12): 29.8 [24.1, 35.6] pp
- **Criterion:** V1 - V0 CI > 0 under every reading and the range covers every reading's estimate. **Holds under every reading:** yes

## A4 (abstract; lead; conclusions)
**Claim.** Raw CHF segments were positive 22.4 points more often than healthy ones, but -4.9 to +5.4 points (not significant) after ectopy handling

- R0: raw 22.4 [12.8, 31.7]; after handling -4.9 [-13.7, 4.3], 5.4 [-3.9, 14.7], 2.9 [-5.3, 11.4]
- R0-BP96: raw 22.3 [12.7, 31.8]; after handling -7.1 [-15.9, 1.8], 2.2 [-6.9, 11.5], 1.1 [-7.5, 9.8]
- R0-ORIG: raw 26.7 [16.2, 36.8]; after handling -3.9 [-12.5, 5.1], 7.6 [-2.0, 17.4], 6.4 [-2.5, 15.4]
- R0-ORIG-minus-S3b: raw 15.9 [7.9, 23.6]; after handling -12.0 [-23.6, -0.0], -2.6 [-13.5, 8.3], -2.2 [-12.8, 8.5]
- **Criterion:** raw CI > 0 and every handled-arm CI includes 0, under every reading. **Holds under every reading:** no
- **Suggested wording:** After editing or Wu-style deletion, CHF and healthy segments did not differ under any reading; after masking they did not differ under R0, R0-BP96 and R0-ORIG, but CHF was lower under the configuration closest to the published benchmarks (see values)

## V-A (Sec. V A; Discussion)
**Claim.** Raw window-level CHF excess +23.5 (R0), +14.2 (original conventions), +4.2 (AAFT surrogates)

- R0: 23.5 [12.6, 34.3] pp
- R0-BP96: 14.2 [4.6, 23.7] pp
- R0-ORIG: 18.0 [6.6, 29.5] pp
- R0-ORIG-minus-S3b: 17.9 [8.4, 27.2] pp
- R0-BP96F (Phase 12): 16.8 [7.0, 26.5] pp
- R0-S-FT (Phase 12): 23.2 [12.1, 34.6] pp
- R0-S-AAFT (Phase 12): 4.2 [-5.0, 13.8] pp
- **Criterion:** raw window CHF - NSR CI > 0 under every reading (known to fail for AAFT). **Holds under every reading:** no
- **Suggested wording:** keep the manuscript's wording: the raw excess is reading-dependent and disappears with amplitude-adjusted surrogates

## A5 (abstract; lead; conclusions)
**Claim.** On 2 older databases, Poon and Merrill's procedure did not reproduce their lower positive rate in CHF on raw intervals; it appeared only after ectopy handling (26-46 points lower)

- R0: NSR raw 52.6 %; raw diff 16.5 [-1.9, 33.3]; handled -45.6 [-57.0, -33.7], -26.4 [-42.8, -8.2], -27.9 [-44.1, -10.3]
- R0-BP96: NSR raw 74.3 %; raw diff -9.6 [-27.0, 7.0]; handled -62.1 [-74.8, -46.6], -50.4 [-64.7, -34.6], -48.1 [-62.1, -32.8]
- R0-ORIG: NSR raw 69.6 %; raw diff -3.4 [-22.3, 14.3]; handled -60.0 [-72.4, -45.7], -43.6 [-60.6, -24.6], -44.4 [-60.3, -26.7]
- R0-ORIG-minus-S3b: NSR raw 81.8 %; raw diff -7.5 [-22.9, 6.8]; handled -65.2 [-77.1, -51.2], -46.2 [-63.9, -26.6], -45.5 [-62.2, -27.0]
- **Criterion:** raw PM97 pattern not reproduced and every handled-arm CI < 0, under every reading. **Holds under every reading:** yes
- **Suggested wording:** holds under every reading; after handling CHF was 26-65 points lower across readings

## A7 (abstract; Discussion)
**Claim.** Titration also detects some static transforms of linear noise (N5 100/100)

- R0: 100/100
- R0-BP96: 100/100
- R0-ORIG: 100/100
- R0-ORIG-minus-S3b: 100/100
- R0-S-FT (Phase 12): 100/100
- R0-S-AAFT (Phase 12): 100/100
- **Criterion:** N5 >= 50/100 under every reading. **Holds under every reading:** yes

## C3 (conclusions; Sec. V B)
**Claim.** Positives increased steeply with ectopic burden (GEE OR 2.31); the same within subjects (OR 3.40)

- R0: GEE 2.31 [1.67, 3.19]; within-subject 3.40 [2.62, 4.41]
- R0-BP96: GEE 1.55 [1.34, 1.80]; within-subject 1.75 [1.45, 2.11]
- R0-ORIG: GEE 2.15 [1.73, 2.68]; within-subject 2.60 [2.03, 3.32]
- R0-ORIG-minus-S3b: GEE 2.39 [1.60, 3.56]; within-subject 3.11 [2.30, 4.20]
- **Criterion:** both OR CIs > 1 under every reading. **Holds under every reading:** yes

## C2 (conclusions; Sec. IV A)
**Claim.** 1-8 isolated premature beats made 91-99 % of linear windows positive; 26 and 51 beats 5 % and 0 %

- R0: k=0: 0 %, k=1: 90 %, k=2: 96 %, k=4: 100 %, k=8: 92 %, k=16: 39 %, k=26: 4 %, k=51: 0 %
- R0-BP96: k=0: 0 %, k=1: 90 %, k=2: 98 %, k=4: 100 %, k=8: 100 %, k=16: 90 %, k=26: 48 %, k=51: 0 %
- R0-ORIG: k=0: 0 %, k=1: 90 %, k=2: 98 %, k=4: 100 %, k=8: 100 %, k=16: 96 %, k=26: 78 %, k=51: 9 %
- **Criterion:** descriptive (exploratory dose-response). **Holds under every reading:** descriptive

