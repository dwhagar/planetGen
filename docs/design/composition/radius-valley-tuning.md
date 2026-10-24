# The radius valley and the envelope: the tuning of rows e1 and e2 (GEN.509)

Status: built 2026-10-24 (lane "Source and tune atmosphere draws"). Module `src/planetgen/physics/radius_valley.py`, two targets in `physics/data/tuning_targets.json` (`e1-valley-slope`, `e1-valley-radius`), rows e1 and e2 of `scripts/atmosphere_recipes_data.py`, tests in `src/tests/test_radius_valley.py`. It uses the Kepler recovery model (GEN.504, `kepler-recovery.md`) and the Lopez and Fortney 2014 radius law of `bulk_structure.py`. Decision D2 of `tunable-gaps.md` (XUV energy-limited loss with variable efficiency, a core-mass distribution as a third knob, Lopez and Fortney as the radius law) is followed.

## 1. What is tuned and what to

Rows e1 (the planet keeps its envelope) and e2 (it lost it) weigh `min(1, f_lost / f_env0)` against its complement. Two constants matter: the starting envelope fraction `f_env0` and the efficiency of the loss. Nothing measures them for one planet; the population does, through the valley in the period-radius plane (Petigura et al. 2022, arXiv:2201.10020, READ in the worker report): slope d log R / d log P = -0.11 +- 0.02, and a valley radius of 1.84 +- 0.03 R_E (6371 km) at 10 days.

`radius_valley.survey` draws 4000 planets (host 0.8 to 1.3 solar masses, period log-uniform 2 to 100 d, core mass log-normal around 3.72 R_E-masses, the peak of Rogers and Owen 2021, a starting envelope log-normal around `f_env0`), applies the loss of the row, the Lopez and Fortney radius, and weights each planet by the chance Kepler would have found it. `valley_of` finds the minimum of the smoothed radius histogram in five period bins from 2 to 64 d and fits a line. One survey is one sample of the harness; 24 surveys are compared with each target (the harness's 2 sigma rule, sigma combining the published error and the spread of the surveys).

## 2. The result: one ridge, not a point

The two knobs are degenerate, as the tunable-gaps note predicted: a stronger loss and a larger starting envelope give the same valley. The surveys put both targets inside 2 sigma along `eta = 10 f_env0` for `f_env0` from 0.035 to 0.08 (slope -0.13 to -0.14, radius 1.82 to 1.85). Below 0.035 the slope steepens past -0.16; above 0.09 the envelopes of the starting draw reach giants. The row therefore draws `f_env0` log-uniform in 0.035 to 0.08 and sets its efficiency to 10 f_env0 (times `(15 km/s / v_esc)^2`, Owen and Wu 2017), instead of two independent constants; an off-ridge pair does not reproduce the valley (a test pins that with a no-loss and a ten-fold loss).

## 3. What it does not do

- It cannot choose the mechanism. Photoevaporation and core-powered mass loss both fit; separating them needs about 5000 planets (Rogers 2021). The simulated slope, -0.13 to -0.14, is 1 to 1.5 sigma steeper than -0.11 (other published slopes, -0.09 and -0.096, lie further away); the core-powered prediction is -0.11. The harness reports the difference; no knob is bent to remove it.
- The core-mass width (0.4 dex), the envelope width (0.6 dex) and the host sample are design choices, labelled stand-ins; the valley moves little with them in the ranges tried but they were not tuned.
- Mass-dependence (e-link: the valley against stellar mass, slope 0.18) is not tuned: the survey draws one host range. It belongs to a later builder with stellar-mass bins.
- Row e5 (envelope metal fraction) stays partly verified under GEN.509's other half; no JWST population mean exists.
- The recipes are not read by the generator yet, so no generated value changes and no epoch is taken. Row e1's loss uses the planet's total mass and core radius as the row did (`mass_earth`, `core_radius_earth`); the survey uses the core mass. The two differ by the envelope fraction, a few percent.
