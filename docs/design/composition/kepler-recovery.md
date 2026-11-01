# The Kepler transit-plus-recovery model (GEN.504)

Status: built (lane M, 2026-10-23). Module `planetgen/physics/kepler_recovery.py`, data `physics/data/kepler_recovery.json`, fit `scripts/fit_kepler_recovery.py`, tests `src/tests/test_kepler_recovery.py`. Selection-model id `kepler-recovery` in the tuning harness (GEN.503).

## 1. What it is for

The features tuned to Kepler planets (the radius valley, the period-radius plane, the metallicity null of the California-Kepler Survey) are published for the planets Kepler could see. A transit survey sees a planet when the orbit is edge-on, the signal beats the star's noise over enough transits, and the pipeline recovers it. The model gives the chance for one drawn planet, so the harness compares the planets Kepler would have found with the published ones and not every drawn planet.

## 2. The three factors

1. **Geometry.** A circular orbit transits when the impact parameter is under 1, which is `R*/a` of the orientations, `a/R*` from Kepler's third law.
2. **Signal.** `MES = g(b) k^2 sqrt(N) / CDPP(T)`: `k = Rp/R*`; `g(b) = 1.170 - 0.089 b^2 - 0.514 b^4` is the effective depth over `k^2` (limb darkening and aperture correction), fitted to the pipeline's expected MES; `N = 0.962 x 1458.9 d / P` valid transits; `T` the first-to-fourth contact duration of a circular orbit (it reproduces the injected durations to 0.2 percent); `CDPP(T)` the host's combined differential photometric precision, log-linear in duration between 1.5, 3, 6, 12 and 15 hours. A generated host has none, so it takes `ln CDPP = a0 + a1 Kp + a2 Kp^2 + a3 ln R* + scatter z` from the DR25 stars, with the Kepler magnitude from the targets' quantiles. One normal draw `z` is shared by all durations (the residuals at 1.5 and 15 hours correlate at 0.79).
3. **Recovery.** `P = 0.943 GammaCDF(MES; 26.7, 0.312) S((N - 2.88) / 0.257) S((ln 31 h - ln T) / 0.363)`: the form of [R0594] three-transit gate and a cut on durations beyond about 31 hours, fitted by maximum likelihood to the 146,202 usable DR25 injections (44,772 recovered). The curve is fitted to the *modelled* MES, so the error of the modelled MES is part of the curve.

## 3. How well it reproduces the injections

| Check | Observed | Model |
|---|---|---|
| All usable injections, the stars' own CDPP | 44,772 recovered | 44,765 |
| Each MES bin (5 to 50) | see `calibration` in the data file | within 0.02 absolute in every bin with over 2,000 injections |
| Generated-host noise (Kp and R* kept, CDPP drawn from the noise model) | 44,764 | 42,328 (5.4 percent low) |
| Bright stars (Kp 5 to 11) | 14.2 percent | 22.3 percent (over-predicted) |

The generated-host path is a few percent low because the CDPP of real stars has a heavier tail than a log-normal; the bright-star over-prediction is variability and crowding the noise model does not carry. Both are limits, not tuned away.

## 4. The item

A mapping with `period_d`, `radius_earth`, `host_radius_solar`, `host_mass_solar`; optionally `kepmag` and `cdpp_z` (else `observed_sampler` draws them) and `impact_parameter` (given: the planet's orientation is realised, a planet with `b >= 1` has probability 0 and there is no geometric factor; absent: the geometric factor times the recovery averaged over `b`). `Target.field` picks any number of the item.

## 5. What it does not do

- It is the detection step (TPS and DV), not the vetting: the Robovetter removes more, and the injection files for it are not used here.
- Circular orbits; no eccentricity, no multiple-planet interaction, and no transit timing variations.
- It conditions on a Kepler target star. The targets are a flux-limited sample of the Kepler field (Kp 5 to 18); the `kepmag` quantiles are theirs.
- Validity: periods 0.5 to 500 days and radius ratios 0.005 to 0.07 (Earth-size to about 7 Earth radii around the Sun). Outside, the curve is an extrapolation the target's validity range flags.
- Grazing transits (`b >= 1`) are not in the injections; the model gives them zero.
- Stellar mass-dependent noise beyond radius is not modelled, so an M dwarf host gets the noise of its magnitude and radius only.

## 6. Licence

Both tables are public data of the NASA Exoplanet Archive (see `DATA-LICENSE.md`); only fitted numbers and quantiles are kept.
