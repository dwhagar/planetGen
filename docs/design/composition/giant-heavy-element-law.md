# The heavy-element mass of giant planets (GEN.460 input, read this session)

Research note, no behaviour change. GEN.460's text says "Chachan's coefficients are not yet in the item and must be read
first". They are now read. Source states: READ is a number found in the paper itself.

## 1. [R0710] ([R0710], 30 Sep 2025), READ

Sample: 147 warm giants (equilibrium temperature under 1000 K) with measured masses and radii, from 20 Earth masses to 20
Jupiter masses, from the Exoplanet Archive; thermal evolution models with updated hydrogen-helium equations of state and an
atmospheric boundary condition that varies with the envelope metallicity. More than three times the sample of Thorngren
et al. 2016.

| Quantity | Value | Notes |
|---|---|---|
| Law | `M_Z = M_core + f_Z (M_p - M_core)` (their equation 6) | `M_Z` is the heavy-element mass (everything but hydrogen and helium), bulk |
| `M_core` | 14.7 (+1.8 -1.6) Earth masses | set mainly by the lower-mass planets, where `M_Z` is about `M_core`; 15.9 +- 1.8 without 12 low-mass hot planets (Teq over 900 K, up to 0.3 Jupiter masses) |
| `f_Z` | 0.09 +- 0.01 | the heavy-element share of the gas accreted beyond `M_core`; about 7 times the solar 0.014 |
| Astrophysical scatter | `sigma_mult = 0.66 +- 0.08`, a fraction of the predicted `M_Z` (0.67 +0.08 -0.07 when normalised by the stellar metallicity) | applied as a spread proportional to `M_Z`; the paper does not give the distribution's shape beyond "scatter as a fraction of the model value" |
| Below `M_core` | the model assumes `Z_p = 1` up to `M_core`, so `M_Z = M_p` there (the straight formula would exceed the mass) | stated in the paper's text under equation 6 |
| Bulk metallicity | `Z_p = M_Z / M_p`; declines linearly with `M_p` for `M_p << 150` Earth masses and flattens at `f_Z` | |
| Stellar normalisation | `Z_p / Z_star = a / M_p + b`, `a = 2.6 +- 0.3` (M_p in Jupiter masses), `b = 3.3 +- 0.5` | enrichment of 29.3, 5.9 and 3.6 for 0.1, 1 and 10 Jupiter masses |
| Against a power law | the model beats `M_Z = alpha M_p^beta` with a Bayes factor of 921 | so [R4352]'s form is the older, poorer fit |
| Classical core accretion | `Z_p = 1` at 10 and `0.5` at 20 Earth masses is inconsistent with the population; `Z_p` reaches 0.5 at about 32 +- 3 Earth masses | |

Not read in the paper (and so not used): the fit to individual planets, the hot-planet sample (inflation changes the inferred
`M_Z`), planets above 20 Jupiter masses.

## 2. Against [R4352] ([R4352], read in an earlier lane, part C of the body composition research)

`M_Z = 57.9 (M / M_J)^0.61` Earth masses, scatter 0.26 dex (factor 1.82), 47 giants from 20 Earth masses to 20 Jupiter
masses. At 1 Jupiter mass: Thorngren 57.9, Chachan `14.7 + 0.09 (317.8 - 14.7)` = 42.0 Earth masses (COMPUTED); at 0.1
Jupiter mass (31.8 Earth masses): Thorngren 14.2, Chachan 16.2; at 10 Jupiter masses (3178): Thorngren 236, Chachan 299.
The ratio of 1.4 at Jupiter's mass agrees with the open-gaps report's "1.4 to 1.6 times lower between 0.3 and 3 Jupiter
masses". Under the source rule (two verified fits that disagree: the newest wins), **[R0710] is the default and
[R4352] stays selectable and the other end of the bracket**.

## 3. How the heavy-element mass reaches the interior draw (GEN.642)

`physics/giant_interior.py` takes `heavy_mass_earth` as an input. GEN.460 draws it: `M_Z = min(M, M_core + f_Z (M - M_core))`
times a draw with the paper's fractional spread, kept inside `0 < M_Z < M` and above the envelope's own metallicity floor
(by rejection, a stand-in until the distribution's shape is read from the paper's likelihood). For the envelope the
atmosphere's `Z` is `f Z_p` with `f` about 0.2 for Jupiter and Saturn ([R4350], read in the Atmosphere
lane's notes); `GiantSubject.envelope_z` carries it.
