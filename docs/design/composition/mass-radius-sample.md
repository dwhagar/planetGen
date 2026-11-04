# The mass-radius sample selection model (GEN.505)

Status: built 2026-10-23 by lane M. Code `src/planetgen/physics/mass_radius_sample.py`, fit script `scripts/fit_mass_radius_sample.py`, data `src/planetgen/physics/data/mass_radius_sample.json`, tests `src/tests/test_mass_radius_sample.py`. It is the selection model `mass-radius-sample` of the tuning harness (GEN.503, `docs/design/composition/tunable-gaps.md` section 4) for the water-rich share (B2, GEN.510) and the iron branch (B4, GEN.512).

## 1. What it answers

A density distribution of small planets ([R2774]; the density modes 0.94, 0.47 and 0.24 are recorded in `tunable-gaps.md` section 3.1) is measured on planets with mass better than 25 percent and radius better than 8 percent. Comparing the generator's drawn planets with those numbers needs the same cut: which of the drawn planets would have been measured that precisely.

## 2. The model

For a transiting planet the cut depends on two errors:

- the radial-velocity semi-amplitude `K = 0.08946 M P^(-1/3) M*^(-2/3)` m/s (M in Earth masses, P in years, M* in solar masses, circular edge-on orbit) and its error `sigma_K`, so `sigma_M/M = sigma_K/K`;
- the transit depth `(R/R*)^2` and its error `sigma_depth`, so `sigma_R/R = sigma_depth / (2 depth)`.

`ln sigma_K` is linear in `ln K` and `ln sigma_depth` in `ln depth`, each with a log-normal scatter, the two scatters correlated. The lines are fitted (`scripts/fit_mass_radius_sample.py`) to the NASA Exoplanet Archive's table `pscomppars` (retrieved 2026-10-23): the transiting planets with a published mass (provenance `Mass`, not a minimum mass) and both errors. A planet passes when `sigma_M/M < 0.25` and `sigma_R/R < 0.08`. Only the fitted numbers are in the repository; the archive rows are not.

| Sub-sample | Planets | Pass | sigma_K line (slope, scatter) | sigma_depth line (slope, scatter) | Correlation |
|---|---|---|---|---|---|
| `small`: radius under 4 Earth radii | 599 | 311 (52 percent) | 0.89, 0.81 | 0.84, 0.85 | 0.40 |
| `m-dwarf`: the same, hosts under 4000 K | 100 | 76 (76 percent) | 0.98, 0.76 | 0.71, 0.57 | 0.37 |

`pass_probability(item, subsample)` is the chance with no random numbers (a bivariate normal integral). `observe(item, subsample)` draws the two errors and returns the measurement, with the measured mass, radius and density; `observed_sampler` wraps a sampler for `tuning_harness.simulate`, and the selection model keeps the items that pass.

## 3. What it shows (COMPUTED)

- The slopes are close to 1: the relative error of a published mass is nearly independent of the size of the signal. Observers follow a planet up until it is precise, so the cut is a weak filter on mass, radius and period. Around a 0.3 solar-mass host at a 4 day period, an Earth-size planet passes with 0.44, a 5 Earth-mass planet of 1.8 Earth radii with 0.54 and a 20 Earth-mass planet of 3.5 Earth radii with 0.63 (`small`; the M-dwarf sub-sample gives 0.72, 0.78 and 0.81). A density distribution is tilted mildly toward the heavier and deeper-transit planets and is not cut hard.
- Applied to the Gaussian summary of the archive's planets (`population_pass_fraction`) the model gives 0.47 for `small` against the 0.52 counted, and 0.74 for `m-dwarf` against 0.76. The `small` gap is the Gaussian summary of a skewed population; the test tolerance is 0.06.

## 4. What it does not do

- It conditions on a planet that transits and has a published mass and radius. Whether it transits, is detected and is chosen for follow-up is not here (GEN.504 is the Kepler detection; there is no TESS model). The selection is therefore a relative weight inside the measured planets, not a count.
- The host's brightness is not an input: a generated system has no apparent magnitude, and the fit has the brightness inside the scatter. A host cooler than 4000 K has its own fit.
- The two errors are independent of the planet's other properties and Gaussian in the log; the archive mixes instruments and methods (transit-timing masses are in the sample).
- Validity: radius under 4 Earth radii; a larger planet is an extrapolation and the target's validity range flags it.
- The values depend on the archive's content on the retrieval date. Re-running the script after the archive grows changes the file and is a table version change.

## 5. The weak cut, as a finding (GEN.665)

The cut is kept as built, not strengthened. A stronger cut would need the host's apparent magnitude and the follow-up choice, and a generated system has neither, so the slopes near 1 of section 3 are the honest result: the cut tilts a density distribution mildly and does not remove light or shallow planets. A target tuned through it therefore sees almost the planets' own distribution, and its power is that of the published sample alone. Section 3 says how mild the tilt is; the tuning notes of B2 and B4 (GEN.510, GEN.512) read it that way.
