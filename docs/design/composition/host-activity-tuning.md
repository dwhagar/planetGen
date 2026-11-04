# Host-star activity: the rotation track, row p6 (GEN.523, partly built)

Status: partly built 2026-10-24 by lane M (builder B15). Row p6 of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); tests in `src/tests/test_tuned_draws.py`.

## 1. p6: the rotation track

The row picked one of three percentile tracks (10th, 50th, 90th) of the initial rotation distribution of Tu 2015 with weights 0.25, 0.5 and 0.25. No paper gives those weights (they are the probabilities of quartile bins). The row now draws the percentile `u` uniform from 0.02 to 0.98, the range of the tracks [R2183] offers at every integer percentile, and the saturation time follows from `u`: the logarithm of the saturation time is linear in `u` between the published tracks, 5.7 Myr at the 10th, 23 Myr at the 50th and 226 Myr at the 90th percentile, extended by the same slope to the 2nd (4.3 Myr) and 98th (357 Myr; Tu 2015 states saturation lasting from about 10 to about 300 Myr). A uniform percentile needs no tuning, since the percentile is the observed 150 Myr distribution; the 200-fold spread of X-ray luminosity at 100 Myr that this draws is real. The drawn percentile is stored as `rotation_percentile` for the rows after it.

## 2. p4: checked, not changed

Row p4 already carries the mass-dependent saturation time (0.162 Gyr at 0.55, 0.6 at 0.25 and 2.3 at 0.15 solar masses; [R2118, R3017]) from the earlier correction. A test pins the three anchors and their order.

## 3. What is not built

- **The mass dependence of the tracks.** Tu 2015 covers 0.9 to 1.1 solar masses only. Johnstone's tracks (MIT-licensed code, not downloaded) cover 0.1 to 1.25 and would give the saturation time as a function of mass and percentile; until they are fitted the one set of three times was used for every mass until the saturation time of section 4 replaced it (GEN.493 wiring).
- **The consistency tests.** The simulated saturated fraction against age and mass, compared with Medina's 2.4 +- 0.3 Gyr, the rotation-period fractions (74, 19 and 8 percent) and the age of rapid rotators (2.0 +- 1.2 Gyr), needs a stellar-population model with the volume-complete and X-ray selection of those samples; it is the remainder of GEN.523 and is filed to the TODO thread.
- The recipes are not yet read by the generator, so no generated value changes and no epoch is taken.

## 4. The mass dependence of the saturation time (GEN.493, lock-free step)

`physics/saturation_time.py` closes the first gap of section 3 without the model's code: Fig. 10 of [R2183] ([R2183], READ) plots the age at which a star falls below the X-ray saturation threshold against mass for the slow, medium and fast rotators, the observed 5th, 50th and 95th percentiles of the rotation distribution at 150 Myr. The table is that figure read by pixel position at thirteen masses (registry RS1, T3, L0).

- At one solar mass the curves give 5.2, 22.8 and 297 Myr, Tu et al.'s 5.7, 23 and 226 Myr at the 10th, 50th and 90th percentiles are the same numbers where they overlap, so the solar-mass times row p6 already used stay right.
- Lower mass stars stay saturated for gigayears: the median is 4.5 Gyr at 0.1, 1.5 Gyr at 0.3, 0.49 Gyr at 0.5 and 0.15 Gyr at 0.7 solar masses; the mass-function weighted mean over 0.1 to 0.3 solar masses is about 2.7 Gyr, against [R3017]'s measured 2.4 +- 0.3 Gyr for the transition (kinematic ages of 219 stars of those masses; [R3017], READ), a check the figure was not tuned to.
- The percentile of an observed distribution is uniform by construction, so row p6's uniform draw of the percentile (2nd to 98th) needs no weights of its own: that is the distribution, not a design choice. The "no paper weights the tracks" limit of the row goes with the wiring step.
- The activity table's `t_sat` (`ACTIVITY_TABLE`, design defaults good to 0.5 dex) agrees at 0.3 solar masses (2.0 against 1.5 Gyr), 0.45 (1.0 against 0.63) and 0.6 (0.4 against 0.3), and is 4 times high at 1 solar mass (0.10 Gyr against the 23 Myr median; [R2118]'s 63 to 200 Myr is the median of the cluster samples and sits between them: 10^7.8 to 10^8.3 yr, with the L_X/L_bol saturation at 10^-3.1 for late K to 10^-4.3 for early F). The wiring step (below) takes the star's own mass and percentile from the module; beyond 1.2 solar masses the module holds the 1.2 value, and the table's saturation ages are gone.

Not built: the model itself (the spin-down torque and the core-envelope coupling) and the decayed level after saturation (the exponent `b` of the table).

**Wired (GEN.493, merge train 7).** `activity.generate_activity` draws a star's initial-rotation percentile uniform from 2 to 98 (`PERCENTILE_RANGE`, the one draw that replaces the log-normal scatter on `t_sat`, so the run's stream is unchanged) and takes the saturation age from `saturation_time` at the star's mass and that percentile. `ACTIVITY_TABLE` is now `DECAY_EXPONENTS` (the exponent only); `coronal_xuv_at`, `reference_fluence_j` and `runaway_history` read the median age (`saturation_age_gy`). The old exponents were calibrated with the old ages, so from 0.45 solar masses up they are retuned (registry AB1, T9, L3) to hold each mass's `L_X/L_bol` at 5 Gyr where the table put it: the Sun still comes out at 1.3e27 erg/s today, with an exponent of 1.45 inside the Skumanich band. What moves: a star's saturated stretch and its fluence now follow its own spin (a slow solar-mass rotator leaves saturation at 5 Myr, a fast one at 300), so the water a Venus-orbit Earth loses to the runaway greenhouse around the median Sun is two oceans, not seven. Row p6 of the atmosphere recipes reads the same module through the recipe function `saturation_time_myr`; it draws its own percentile and does not yet read the star's, which waits for the recipes to be wired to the generator.
