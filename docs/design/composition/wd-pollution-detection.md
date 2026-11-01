# The white-dwarf pollution detection model (GEN.508)

Status: built 2026-10-23 by lane M. Code `src/planetgen/physics/wd_pollution.py`, data `src/planetgen/physics/data/wd_pollution_sample.json`, tests `src/tests/test_wd_pollution.py`. It is the selection model `wd-pollution-detection` of the tuning harness (GEN.503, `docs/design/composition/tunable-gaps.md` section 4) and serves the redox builder B5 (GEN.513) for rows a3 and a4.

## 1. What is measured, and what is reported

A polluted white dwarf with oxygen, magnesium, silicon and iron measured gives the oxygen fugacity of the body that fell on it: dIW = 2 log10(x_FeO / x_Fe,metal) with x_Fe,metal = 0.85 ([R1152], equations 3 and 4; [R1152]). The abundances have errors, and the papers propagate them by Monte Carlo draws; a draw with a negative FeO mole fraction is unphysical. The rule they report by ([R1152] section 2.3, applied again in [R1153]):

| Negative draws | Reported | Stars |
|---|---|---|
| under 16.5 percent | median and both bounds (`both`) | 16 |
| 16.5 to 50 percent | median and an upper bound (`upper`) | 6 |
| 50 percent or more | an upper limit only (`limit`) | 2 |

So a body with little FeO next to the abundance error is reported with a large error, or not at all. The 22 medians (mean -1.248, median -1.256, sample sd 0.447, recomputed from the two tables) are a sample cut by this rule and scattered by the measurement, and the tuned width of rows a3 and a4 must be compared with that sample and not with the bodies.

The earlier worker report (`tunable-gaps-sources/group3-redox-wd-radio.md`) counted 7 of 24 without a lower bound; the two tables give 8 (six with a median and an upper bound, two upper limits). The 24 rows, with the effective temperature of each star, are in the data file; Table 1 of [R1152] holds 16 stars and Table 4 of [R1153] holds 8 (temperatures from its Table 2).

## 2. The model

1. The abundance error on x_FeO is absolute: the measured fraction is `x_meas = x_true + sigma z0`, `z0` a standard normal, and `sigma` is log-normal across stars. A star with little FeO is therefore the one that is flagged.
2. The two parameters of `ln sigma` are fitted by maximum likelihood to the sample. The 16 stars with both bounds give `sigma` as half the spread of x_FeO between the 1-sigma bounds. A star with a median and an upper bound has at least 16.5 percent negative draws, so its `sigma` is at least `x / 0.974` (the normal quantile of 16.5 percent). A star with a limit only has no median to give a bound: it takes the largest bound of the others, a convention that says only that its error is large. Result: ln sigma has mean -2.18 and standard deviation 1.01 (a typical error of 0.11 in x_FeO, a factor 2.7 per sigma). Fitted to the 16 stars with both bounds alone the parameters are -2.75 and 0.52: that fit leaves out the stars whose error was large, so it is the optimistic bracket (`BOUNDED_ERROR_MODEL`).
3. The reported value is `2 log10(x_meas / 0.85)` with the class the rule gives from `x_meas / sigma`. A `limit` has no value, and the selection model `reported_median` keeps `both` and `upper` stars.
4. `expected_sample(mean, sd)` gives the sample a population N(mean, sd) shows, by quadrature over the body, the error of the star and the realisation of the abundances: no random numbers and no seed.

Validity: effective temperature 10,120 to 21,200 K; DA, DB and DZ atmospheres; bodies whose oxygen is measured.

## 3. What the model shows (COMPUTED, `test_wd_pollution.py`)

| Population N(-1.25, 0.2) | Both | Upper | Limit | Mean of reported | Sd of reported |
|---|---|---|---|---|---|
| Published sample (24 stars) | 0.667 | 0.250 | 0.083 | -1.248 | 0.447 |
| Whole-sample error model | 0.719 | 0.177 | 0.104 | -1.230 | 0.654 |
| Bounded-stars error model | 0.916 | 0.066 | 0.018 | -1.302 | 0.471 |

- The whole-sample model recovers the published classes inside their 2-sigma Wilson intervals and the mean inside 2 standard errors; that is the Done test. The bounded-stars model cannot: it gives 1.8 percent limits where 2 of 24 were seen and 7 percent upper bounds where 6 of 24 were.
- The spread of the reported medians is explained by the measurement. At zero intrinsic width the two models give 0.419 (bounded) and 0.625 (whole sample) against the published 0.447 (standard error about 0.07). The published value lies between the two, so the intrinsic width of the bodies is **not identified by the sample**: the bounded model needs about 0.15, the whole-sample model a width of zero (it overshoots). A width of 0.45, the recipe value of rows a3 and a4, is too wide in both. The crude intrinsic estimate of the worker report (0.22, from 13 well-constrained stars) is the same order.
- The centre: a population whose reported mean is -1.248 has a true centre of -1.27 (whole-sample model) to -1.20 (bounded model) at a width of 0.15, so the centre -1.25 survives the selection to within the standard error 0.10.

## 4. What it does not do

- The shipped rows a3 and a4 were retuned with this model (GEN.513, `docs/design/composition/wd-redox-tuning.md`): the width 0.45 became 0.22 and the four targets now name this model. The conversion between the core-formation dIW of the white dwarfs and the surface magma-ocean scale of row a1 is not tested.
- The settling of elements out of the atmosphere ([R1152] section 3.2 and Figs. 8 to 11) is not modelled: a DB observed long after accretion has waned can show a more reduced body than fell on it (the papers' time-dependent model, diffusion times, is the reference). The model's case is a body observed in its build-up or steady-state phase; a tuned width is an upper bound where the sample may hold stars observed later.
- The detection of pollution itself (which white dwarfs are polluted, 27 to about 50 percent, [R2415]) gates "rocky debris exists" and is not part of this selection: every simulated body is one that fell.
- The convention for the two limit stars' error and the Gaussian form of the abundance error are model choices; the fit changes if either is dropped (the two brackets above).

## 5. Use

    sampler = wd_pollution.observed_sampler(lambda i: draw.gauss(mu, sd))
    items = tuning_harness.simulate(target, sampler, 2000)        # target: selection_model "wd-pollution-detection", field "diw"
    tuning_harness.compare(target, items)

The harness imports the module the first time a target names the model (`tuning_harness.selection_model`).
