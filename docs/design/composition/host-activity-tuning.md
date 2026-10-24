# Host-star activity: the rotation track, row p6 (GEN.523, partly built)

Status: partly built 2026-10-24 by lane M (builder B15). Row p6 of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); tests in `src/tests/test_tuned_draws.py`.

## 1. p6: the rotation track

The row picked one of three percentile tracks (10th, 50th, 90th) of the initial rotation distribution of Tu 2015 with weights 0.25, 0.5 and 0.25. No paper gives those weights (they are the probabilities of quartile bins). The row now draws the percentile `u` uniform from 0.02 to 0.98, the range of the tracks Johnstone 2021 offers at every integer percentile, and the saturation time follows from `u`: the logarithm of the saturation time is linear in `u` between the published tracks, 5.7 Myr at the 10th, 23 Myr at the 50th and 226 Myr at the 90th percentile, extended by the same slope to the 2nd (4.3 Myr) and 98th (357 Myr; Tu 2015 states saturation lasting from about 10 to about 300 Myr). A uniform percentile needs no tuning, since the percentile is the observed 150 Myr distribution; the 200-fold spread of X-ray luminosity at 100 Myr that this draws is real. The drawn percentile is stored as `rotation_percentile` for the rows after it.

## 2. p4: checked, not changed

Row p4 already carries the mass-dependent saturation time (0.162 Gyr at 0.55, 0.6 at 0.25 and 2.3 at 0.15 solar masses; Jackson 2012, Medina 2022) from the earlier correction. A test pins the three anchors and their order.

## 3. What is not built

- **The mass dependence of the tracks.** Tu 2015 covers 0.9 to 1.1 solar masses only. Johnstone's tracks (MIT-licensed code, not downloaded) cover 0.1 to 1.25 and would give the saturation time as a function of mass and percentile; until they are fitted the one set of three times is used for every mass.
- **The consistency tests.** The simulated saturated fraction against age and mass, compared with Medina's 2.4 +- 0.3 Gyr, the rotation-period fractions (74, 19 and 8 percent) and the age of rapid rotators (2.0 +- 1.2 Gyr), needs a stellar-population model with the volume-complete and X-ray selection of those samples; it is the remainder of GEN.523 and is filed to the TODO thread.
- The recipes are not yet read by the generator, so no generated value changes and no epoch is taken.
