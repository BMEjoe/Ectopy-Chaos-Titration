# Draft text for the revision (M1, M2)

All numbers come from `results/revision_m1_m2/` (exploratory analyses; plan committed before any run). Section and
table numbers follow FINAL_DRAFT_4; new items are marked "new". Replace "Ref. [72]" with the archive reference of the
revised repository.

---

## Methods (new subsections, after Sec. II I)

### J. Ablation of implementation details (exploratory)

A referee pointed out that my test might be less sensitive than the original one. If so, large outliers such as
ectopic beats would make up a larger share of its positives. To find out which of my implementation details matter,
I listed every way my implementation differs, or might differ, from the original papers [9-11] (Table S-D). Nine
can be switched one at a time: N in the criterion (series length instead of fitted rows), the plain ratio of error
variances in the F test, a χ² check with a Mann-Whitney fallback, the Mann-Whitney test alone, the surrogate step,
polynomial degree up to 8, embedding dimension up to 7, a zero-padded start in which every point is a target, and a
decision by the test alone. For the flows of Table 1 in [9], I also chose the sampling step that maximizes the
difference between the linear and nonlinear errors, as [9] describes.

I measured four benchmarks against their published values: the noise tolerances in Table 1 of [9] (100 seeds per
system), the correlated-noise level at which the ecological model of [9] is lost (100 seeds), the noise limit of the
logistic map at r = 3.7 [11] (100 initial conditions), and the positive rate of healthy nsrdb segments under the
protocol of [10] (18 subjects). For each benchmark I computed the fraction of my implementation's gap to the
published value that a switch closes, with paired bootstrap intervals. A rule fixed before the runs combined every
switch that closed the gap on average into one configuration, which I call R0-ORIG. I checked every configuration on
the verification controls of Sec. II E, and I repeated the headline analyses under R0, R0-BP96 and R0-ORIG. The plan
was committed before any run (Ref. [72]).

### K. Influence of single equations (exploratory)

To find out why one premature beat can flip the test, I inserted k = 0 to 51 beats into 200 linear windows (N1) per
condition. Besides ventricular-type beats, I used atrial-type beats, a single additive spike, a missed beat (two
intervals merged) and an extra detected beat (one interval split in two), placed at random or periodically. In every
window I recorded the selected models. I split the drop in residual sum of squares from the best linear to the best
nonlinear model into the part from equations that contain an affected interval (as the target or as one of the six
nonlinear lags) and the rest. I computed the leverage and Cook's distance of every equation. I then refitted after
removing only the affected equations (about eight per beat), after removing as many equations around random
positions without ectopy (sham), and after removing the m equations with the largest Cook's distance (m = 1 to 16).

---

## Results

### M1 (new paragraphs, end of Sec. III C)

Which of my implementation details make the test less sensitive than the original? Two of the original's
conventions explain part of the difference: N as the series length in the criterion, and the plain ratio of error
variances. Together they raised the positive rate of healthy nsrdb segments from 52.6% to 79.0% (95% CI 72.8-84.7)
and the noise limit of the logistic map at r = 3.7 from 0.66 to 0.71 (mean of 100 initial conditions). Allowing
polynomial degree up to 8 raised them further, to 81.8% and 0.74. With the Mann-Whitney test as well (R0-ORIG), the
noise limit reached 0.751, the published value. No configuration moved the ecological model toward the published
50%: my implementation lost it at 35%, and the Mann-Whitney test lost it at 30%. No configuration agreed with Table
1 of [9] for more than 6 of its 9 systems, and two of my four misses there are more noise-tolerant than published.
Choosing the sampling of the flows as [9] describes did not help. A zero-padded start moved the ecological benchmark
past 60%, but it also made every periodic control test positive, which [9] reports did not happen. So the
conventions of the original test explain about half of the gap on two benchmarks and none of it on the other two
(Table S-A, Fig. S-A). All configurations except the zero-padded one passed the verification of Sec. II E.

### M1 (new paragraph, end of Sec. V B)

Do the ectopy results depend on the sensitivity of my test? Under every reading I tested, the positive rate rose with
ectopic burden (odds ratio per doubling 1.55-2.52; within subjects 1.75-3.40), one inserted premature beat added 28
to 64 points (34 under R0-ORIG), couplets made 100 of 100 linear windows positive, and CHF had no effect once burden
was in the model. After editing or deletion of ectopic beats the groups did not differ under any reading. After
masking they did not differ under R0, R0-BP96 or R0-ORIG, but under the configuration that came closest to the
original benchmarks, CHF segments were positive 12.0 points less often (95% CI 0.02-23.6). What changed most was the
baseline. Under the original test's conventions, 31-50% of windows without ectopic beats were positive, against 22%
under mine. So the share of positive windows that I can attribute to ectopy fell from 52% (95% CI 43-62) under my
conventions to 24-41% under the others (64% with phase-randomized surrogates). Ectopy is a major source of titration
positives in heart rate under every reading. It is the main source only under my conventions (Table S-B).

### M2 (new subsection in Sec. IV, after Sec. IV A)

**Why one premature beat matters.** One ventricular-type beat made 179 of 181 linear windows positive when it fell
inside the fitted equations, and none of 19 when it fell before the first fitted interval, where it can only enter the
long linear lags. Inside, the beat touches about eight of the 429 regression equations: those in which the premature
interval or the pause is the target or one of the six nonlinear lags. These equations were 1.9% of the rows but
carried 57% (median) of the drop in residual sum of squares from the best linear to the best nonlinear model. Under
the nonlinear model their leverage was 0.999 (other equations 0.036), and their largest Cook's distance about 900
(other equations 0.03). Squares and products of the first lags are close to zero except at the beat, so the
polynomial can fit these few equations almost exactly. The linear model cannot, because its coefficients are set by
hundreds of sinus equations. Removing only the affected equations made every positive window negative, for 1 to 8
beats and for every beat type I tested. Removing as many other equations made 3% negative. Removing the 8 equations
with the largest Cook's distance made 97% of single-beat windows negative. So the test is driven by a few
high-leverage equations. This was also true under R0-ORIG (Fig. S-B).

**Why the inverted U.** With more beats, the positive rate fell to 39% at 16 beats, 4.5% at 26 (5%) and 0% at 51 (10%).
The polynomial runs out of terms to fit every affected equation: their leverage fell from 0.98 at 8 beats to 0.41 at
26, and the nonlinear model's normalized error on them rose from 0.34 to 0.48 between 8 and 51 beats. The linear
model also took up part of the alternation of short and long intervals: the lag-1 autocorrelation fell from 0.11 to
-0.32, and its error on the affected equations fell from 0.52 to 0.45. In units of the window's standard deviation,
the beats also became less extreme. With periodic placement, the linear model with up to 83 lags described the beats
almost exactly, and no window was positive from 16 beats on. Under the original test's conventions the downslope
moved to higher densities: 26 beats still made 78% (R0-ORIG) and 48% (R0-BP96) of windows positive. The shape also
depends on the beat. A missed beat (one long interval) lost positives from 2 beats on. An extra detected beat (one
interval split in two) kept about 90% of windows positive up to 51 beats. I have not determined why.

### Corrections to existing text

- Sec. IV A, last sentences: replace "I have not determined what feature of a premature beat is picked up by the
  polynomial model (Sec. S12 of the supplementary material)" with "Sec. IV C (new) shows that a few high-leverage
  equations around each beat decide the test."
- Sec. IV A: "For the densities of Table III, isolated beats rarely triggered it" holds only under my conventions.
  Add: "Under the original test's conventions, isolated beats at 5% made 57-82 of 100 windows positive (Table S-B)."
- Sec. II E, "It failed only for 12 periodic parameters of the logistic map": add "and for one periodic van der Pol
  regime, (ρ, ω) = (6, 5.6), which was positive without noise in all 100 windows but had NL = 0".

---

## Alternative abstract sentences

1. "They do: in heart rate, ectopic beats are a major source of titration positives. One premature beat inserted into
   clean windows added 28-64 percentage points, depending on the test's conventions, and ectopic burden, not group,
   predicted a positive window."
2. "They do. A single premature beat made 179 of 181 linear series positive, because a few regression equations around
   the beat dominate the test; in 988 real windows, ectopy accounted for a quarter to a half of positive windows,
   depending on the test's conventions."
3. "They do: under every version of the test I examined, positive rates rose with the number of ectopic beats, and
   one inserted premature beat added 28-64 percentage points. Under the original test's conventions, a quarter of
   positive windows were attributable to ectopy; under mine, half."

## Lead paragraph (replacement for one sentence)

"One ectopic beat raised the chance of titration reporting chaos 1.5- to 2.9-fold, depending on the test's
conventions."

## Conclusions (replacement for the last sentence)

"A positive titration result for heart rate is evidence of nonlinearity, not of chaos. Ectopic beats explain a large
share of it, a quarter to a half of positive windows depending on the conventions, and they explain the excess of
positives in heart failure."

---

## Limitations (updated bullets for Sec. VII G)

- My implementation of titration differs from the original test in four details. Switching them toward the original
  explains about half of the difference on two benchmarks (logistic noise limit, healthy nsrdb rate) and none on the
  other two (ecological model, Table 1). The remaining difference is unexplained; it may lie in details the original
  papers do not report, such as system parameters, sampling and noise generation.
- Under the original test's conventions, more windows without ectopy test positive, so the share of positives
  attributable to ectopy is smaller (24-41%, against 52% under mine). The ectopy effects keep their direction and
  significance under every reading. Under the configuration closest to the original benchmarks, masked CHF segments
  were positive less often than healthy ones, as on the older databases.
- The single-beat effect rests on a few high-leverage regression equations. Any test that compares residual
  variances of fitted models, as the F test does, inherits this sensitivity to isolated outliers.
- The ablation and the influence analysis are exploratory. They were planned after a referee comment, and the
  configuration that best matched the benchmarks was identified after the planned runs.
- The mechanism analysis uses one linear null, one window length and stylized beats.

---

## Supplement table and figure captions (new)

**Table S-D.** Deviations of my titration implementation (R0) from Barahona and Poon (1996), Poon and Merrill (1997)
and Poon and Barahona (2001). For each point: what the original says, my choice, the closest alternative, my
confidence about what the original did, and the switch that tests it. Composite closure is the mean, over four
benchmarks, of the fraction of R0's gap to the published value that the switch closes.

**Table S-A.** Benchmarks by configuration: mean absolute deviation from Table 1 of [9] over 9 systems (100 seeds
each), median correlated-noise level at which the ecological model is lost (100 seeds), mean noise limit of the
logistic map at r = 3.7 (100 initial conditions), and the positive rate of healthy nsrdb segments under the protocol
of [10] (18 subjects). 95% intervals from bootstrap resampling of seeds or subjects.

**Table S-B.** Headline results under R0, R0-BP96 and R0-ORIG: positive rate by ectopic burden, GEE and within-subject
odds ratios, attributable fraction of positive windows, single-beat transplant, CHF minus healthy differences for raw
and handled intervals in windows and 12-min segments, and the replication of Poon and Merrill. 95% intervals as in
the main text.

**Table S-C.** Synthetic ectopy and inserted-beat dose response under each reading (100 or 200 linear windows per
row).

**Figure S-A.** Fraction of the benchmark gap closed by each configuration (1 = published value reached; negative =
further away), with 95% bootstrap intervals.

**Figure S-B.** Positive rate versus the number of inserted beats for five beat types under R0 and R0-ORIG (200 linear
windows per point), with the rate after removing the equations affected by the beats (dashed) and after removing as
many other equations (dotted).

**Figure S-C.** Mechanism of the inverted U: lag-1 autocorrelation, relative improvement of the nonlinear model,
normalized residual sum of squares of both models on affected and other equations, selected linear order and size of
the beats, and representative C(r) curves, for random and periodic placement.

**Figure S-D.** Share of positive windows whose best nonlinear model contains each polynomial term, by number of
inserted ventricular beats.
