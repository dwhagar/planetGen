# The comet admission model (GEN.507)

Status: built (lane M, 2026-10-23). Module `planetgen/physics/comet_admission.py`, data `physics/data/comet_admission.json`, fit `scripts/fit_comet_admission.py`, tests `src/tests/test_comet_admission.py`. Selection-model id `comet-admission` in the tuning harness (GEN.503).

## 1. What it is for

The comet volatile ratios (GEN.495, `comet_population`) and their heliocentric-distance term are published for the comets whose coma was measured, and those are the bright ones that come near the Sun. Tuning the ratios of every drawn comet to that sample would mix two populations. The model says which drawn comets would be in the sample and at what distance each would be measured, so the harness compares like with like.

## 2. The two parts

**Admission.** The share of known comets that have a measured composition, by group (`jupiter_family`, `halley_type`, `long_period`; a parabolic comet, whose `period_class` is None, is long period) and perihelion bin (0.3 to 1, 1 to 2, 2 to 8 AU). The database of [R3782] (909 rows of 96 comets) is joined to the orbits of the 2,422 known comets with perihelion from 0.3 to 8 AU (JPL Small-Body Database); 93 of them have a composition. The probability is the Jeffreys estimate `(k + 1/2) / (n + 1)`, so a cell with no measured comet (Halley-type comets beyond 2 AU) is not zero.

| Group | 0.3 to 1 AU | 1 to 2 AU | 2 to 8 AU |
|---|---|---|---|
| jupiter_family | 4 of 111 | 45 of 252 | 18 of 550 |
| halley_type | 4 of 18 | 1 of 60 | 0 of 32 |
| long_period | 5 of 413 | 7 of 375 | 9 of 611 |

**Distance.** Each composition is measured at a heliocentric distance `Rh`. For parent molecules and for radicals and dust separately, `ln Rh = intercept + slope ln q + scatter eps` is fitted to the rows (parents: 223 rows of 39 comets, slope 0.58, scatter 0.25; radicals: 686 rows of 86 comets, slope 0.34, scatter 0.32). A comet is never measured inside its perihelion, so the draw is truncated at `Rh >= q`: one uniform draw through the inverse CDF of the truncated normal.

## 3. How a target uses it

The item is a mapping with `period_class` (or None) and `perihelion_distance_au`. `observed_sampler(sampler, kind)` adds `rh_au`; the selection model returns `admission_probability(item)`; `Target.field` picks `rh_au` or any other number of the item.

## 4. What it does not do

- It conditions on a known comet. Comets nobody found are not in the fit, so the shares are shares of found comets.
- The nucleus size and the activity of a drawn comet are not inputs: the orbit table has no nucleus size for most comets, and the group and perihelion distance carry most of the signal.
- The admission shares are of a mixed database (many observers, 1980 to 2024); no time dependence is fitted.
- The fit is valid for perihelion 0.3 to 8 AU. A perihelion outside takes the nearest bin, which the harness flags as an extrapolation through the target's validity range.
- Interstellar comets (2I, 3I) are outside the fit; GEN.520 owns them.

## 5. Licence

The database is CC BY 4.0 [R4972]; only fitted numbers and counts are kept. See `DATA-LICENSE.md`.
