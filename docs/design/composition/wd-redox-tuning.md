# Redox of rocky debris: the retune of rows a3 and a4 (GEN.513)

Status: built 2026-10-23 by lane M. Rows a3 (building blocks of asteroids, C and D types) and a4 (debris of white dwarfs) of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`), their four targets in `physics/data/tuning_targets.json`, and the detection model of GEN.508 (`docs/design/composition/wd-pollution-detection.md`). Tests in `src/tests/test_wd_pollution.py` and `src/tests/test_tuning_harness.py`.

## 1. What was tuned

Both rows drew `dIW ~ N(mu, 0.45)` with the centre -1.26 (a3) or -1.25 (a4). The 0.45 is the spread of the 22 *reported* medians of [R1152] and 2023, which is the spread of the bodies plus the abundance error of each star, cut by the papers' reporting rule. Drawing the bodies as wide as the reported medians counts the measurement twice.

The rows now draw `N(-1.25, 0.22)`. The centre is the mean of the 22 values (standard error 0.10) for both rows; the detection model shows it survives the selection to within that error. The width 0.22 is the intrinsic spread the worker report computed from the 13 well-constrained stars (the least selection-biased subset).

## 2. What the targets now compare

The four targets (a3 and a4, mean and sd) name the selection model `wd-pollution-detection` and the field `diw`: a drawn body is measured with a log-normal abundance error, classified by the papers' rule and kept if a median is reported, and the mean and the standard deviation of the kept medians are compared with the published -1.25 (sigma 0.10) and 0.45 (sigma 0.069).

| Width of the draw | Reported sd, error from the stars with both bounds | Reported sd, error from the whole sample |
|---|---|---|
| 0 | 0.419 | 0.625 |
| 0.22 (shipped) | 0.48 | 0.66 |
| 0.38 | 0.583 (the edge of 2 sigma) | 0.711 |
| 0.45 (before) | 0.632 | 0.741 |

The published spread is 0.447 (sigma 0.069), so the accepted width is 0 to about 0.38 with the first error model, and the second model over-disperses at every width, including zero (3.3 sigma at 0.05). The comparison therefore names the first model for the spread; the mean passes with both. This is the finding of GEN.508 in its final form: the sample does not identify the intrinsic width, and 0.45 was too wide. 0.22 is inside the accepted range and agrees with the independent estimate.

## 3. What it does not do

- The scale conversion between the white-dwarf core-formation dIW and the surface magma-ocean dIW of row a1 is not tested. Row a1 keeps its floor and slope (flagged stand-ins); a3 and a4 are on the white-dwarf scale.
- The settling of elements out of a white dwarf atmosphere is not modelled; a tuned width is an upper bound for a sample that holds stars observed after the pollution waned.
- The recipes are not yet read by the generator, so this changes no generated value and takes no generator epoch.

## 4. The power of the targets (GEN.665)

The sample behind the spread has 7 stars with both bounds (16 with both bounds in the error fit, of which the redox values come from 22), not 8 of 24. The tests now say what the four targets can and cannot refuse (`test_wd_pollution.py`, last test): a centre 0.3 dex off and a spread of 0.9 or 0.22 against 0.447 are refused (more than 2 sigma of 0.10 and 0.069), the earlier width 0.45 is refused, and a width of 0.01 is not. The width 0.22 is kept: it is the intrinsic spread of the well-constrained stars and the only value of the accepted range 0 to about 0.38 with an independent estimate behind it. The data cannot tell it from a width of zero, and a better target needs more stars with both bounds.
