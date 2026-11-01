# Giant planets: the retune of rows f5 and f7 (GEN.516)

Status: built 2026-10-23 by lane M (builder B8, `docs/design/composition/tunable-gaps.md` section 4). Rows f5 (heavy-element enrichment of giants), f6 (its bins) and f7 (C/O) of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); seven targets in `physics/data/tuning_targets.json`; tests in `src/tests/test_tuned_draws.py`.

## 1. f5: the enrichment follows the later, larger fit

The row drew `log10(Zp/Z*) ~ N(log10(9.7 M^-0.45), 0.3)`: the [R4352] relation (47 giants) and a recalled design width. The 147 warm giants of [R0710] give `Zp/Z* = 2.6/M + 3.3` (M in Jupiter masses; errors 0.3 and 0.5) from 0.03 to 12 M_J and say the earlier relation is 1.4 to 1.6 times too high from 0.3 to 3 M_J. The later, larger sample wins (Boss's rule for two verified sources). The width is the intrinsic spread of the heavy-element mass, 0.26 dex (band 0.24 to 0.28), not 0.3.

| M (M_J) | Earlier relation | 147-giant fit (target) | Shipped draw (median) |
|---|---|---|---|
| 0.1 | 27.3 | 29.3 +/- 3.0 | 29.1 |
| 1 | 9.7 | 5.9 +/- 0.6 | 5.9 |
| 10 | 3.4 | 3.56 +/- 0.5 | 3.6 |

The targets are the three medians and the spread of the logarithm (0.26 +/- 0.022). The row is bulk enrichment relative to the star; the atmospheric metallicity (the factor 0.2 between bulk and atmosphere rests on Jupiter and Saturn only) is a separate Solar System constraint and is not tuned here. With the later source the row moves from partly verified to verified.

## 2. f7: the C/O core is as wide as the planets

The core was `log10(C/O) ~ N(log10 0.57, 0.15)` truncated to 0.2 to 0.8, a width taken from one planet (WASP-77Ab). LExACoM (69 giants with a retrieved C/O) has 12 below 0.35, 17 percent, where the old core gave 9 percent. The core width is now 0.21 dex, which gives 18 percent. The targets are the median 0.57, the share below 0.35 (12 of 69), the share above 0.8 (6 of 69) and the spread 0.27; every one passes within 2 sigma (the spread at 1.5 sigma: the core is truncated at 0.2, below which the table has planets the row cannot draw, and the table's spread carries the same method systematic of 0.2 to 0.3 as its median).

The tail weight stays the bracket 0.043 to 0.087 of the earlier correction, now a tuned target (6 of 69, Wilson 0.040 to 0.177).

## 3. What it does not do

- The retrievals have no completeness model; the draw is compared with the table on mass only (0.03 to 15 M_J).
- C/O above 0.8 is carried by hot and ultra-hot planets, so it mixes chemistry, retrieval systematics and formation.
- The recipes are not yet read by the generator, so this changes no generated value and takes no generator epoch.
