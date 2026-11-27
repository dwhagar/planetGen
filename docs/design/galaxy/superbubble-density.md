# Superbubbles shape the stars: the bubble field as a density modulation

**Status:** GEN.391 (the table of effects) is decided and GEN.392 (the handle, section 8) is built; the density terms (GEN.393), the scatters (GEN.394) and the sector fill (GEN.395) are not. Boss asked for it on 2026-10-16 at 9:27 PM EDT; GEN.391 turned this plan into the checked design. Boss decided the three open choices on 2026-10-16 at 9:38 and 9:39 PM EDT (section 7).

Boss (2026-10-16, 9:27 PM EDT): "instead of drawing the bubbles, we need a way to store the regions and have them adjust the star density function for the space they occupy, which means finding a way to shift the probability to favor one part of a sector over another. TODO should think about how to do that so when sectors get built they follow the bubbles without having to draw them."

## 1. What exists today (survey of main, 2026-10-16, 9:33 PM EDT)

- `galaxy/superbubble_field.py` (GEN.254) is already derived from the seed and stores nothing: 250 pc columns, each with its own seeded Poisson draw of bubbles (radius log-uniform 50 to 500 pc, centre height about the young layer, ambient gas), cached per column with `lru_cache`. There is no function that cuts a box through bubbles.
- Nothing in it changes a star's count or position. Bubbles only remove molecular clouds (`nebula_field.cloud_kept`), set the cell's gas and ambient medium in the profile (`galaxy/profile.py:302-307`), and are drawn on the Galaxy Map as three faint great circles (`galaxymap3d.js:2604-2639`).
- The stellar density is a function of the shape alone (`galaxy/density.py`: `relative_density`, `population_densities`, `population_densities_array`, `total_density_array`; three copies of the term sum, tested against each other). `ComponentField` (clusters and satellites) is the established way to carry seeded extra structure on the shape (`shape.components`, attached with `shape._replace` at `generation/run_plan.py:313` and `db/store.py:7508`, hashed by identity).
- The sector fill takes its count from a point value at the cell centre (`profile_at`, `profile.py:291-316`, then `generation/run_sector.py:654`) and places stars uniformly in the cell (`galaxy/sector.py:874-880`). The bright-star scatter thins by true density over a certified bound (`generation/object_first.py:98-132` `band_majorants`, `generation/bright_stars.py:433-457`); the phenomena scatter has no certified bound and uses bin-centre values with uniform positions inside the wedge (`generation/phenomenon_scatter.py:312-338`).
- A remnant in a bubble is sized by the cavity's interior density (GEN.298: `SystemConfig.AMBIENT_DENSITY_CM3`, rule in `intergalactic-and-halo-contents.md` section 5); the remnant rate keeps its factor, and the modulation below acts on the young stars only, so the remnants and the stars do not count the bubble twice (GEN.391 checks).
- Planned and different: `region_density_overrides` (DB.38, ADM.94, 8.6) store one 0.1 to 10 multiplier per whole sector region at fill; a bubble is a continuous seeded factor, so the two are complementary and share only the multiplicative read.

## 2. What a bubble does to the stars (young only, Boss 9:38 PM EDT; strengths checked in GEN.391)

A bubble is a gas structure, so it acts on young populations and leaves old ones alone: stars are collisionless and a 100 Myr-old star does not notice a cavity. The effect is a factor on the population weights, not a cut of the whole density.

| Population | Effect | Source named from memory |
|---|---|---|
| Young (under the bubble's age) | Inside: new star formation is suppressed once the gas is swept out. On the shell: swept-up gas collapses, so a share of new stars forms in a thin shell (collect and collapse). | [R4916] for the Local Bubble |
| Intermediate | Half-way between young and old, by age range | same |
| Old, bulge, halo, clusters, satellites | none | none |

The bubble's own OB association belongs to the young-cluster work (GEN.244, components); the field carries only the shell enhancement and the interior suppression, so the association is not counted twice. How strong each term is, how thick the shell is and how the strength scales with the bubble's radius are numbers GEN.391 takes from the papers, and where they are missing it asks Boss for research.

## 3. The modulation

For a point `x` and a bubble of centre `c` and radius `R`, depth `d = |x - c| / R`. Each young or intermediate population `p` gets a factor `m_p(d)`: below 1 inside, above 1 in a shell about `d = 1`, equal to 1 far outside, smooth in `d` (no step), so a bubble that cuts through part of a sector shifts the probability inside that sector toward its rim and away from its inside. Overlapping bubbles multiply. The density the rest of the code sees is `population_densities(x)[p] * m_p(x)`.

- **Bound.** The scatter needs `density <= bound`. A suppression needs no change to the bound, only a looser one; a shell boost needs the bound multiplied by `max(1, m_max)`, with `m_max` certified over the band: the product of `1 + A` of the bubbles that reach the band, or one global figure. `bright_stars._spots` multiplies its `wanted` by the exact `m(point)` and its `check=True` assertion tests the whole product.
- **Count and position in a cell.** The cell's expected count uses the cell mean of the factor (a short sub-sampled quadrature, as `ComponentField` does for the cube), not the centre value. Positions inside the cell follow the factor by rejection against the cell's maximum, so a sector that straddles a rim holds more young stars on the rim side.
- **Phenomena pass.** It thins each candidate by `m(x) / m_max` at its position, because its bin-centre counts lose the within-wedge shape.
- **No recursion.** The bubbles' own placement reads the gas tracer (`nebula_field.gas_factor`), which reads the young density; the modulation must not feed back into that read. The placement reads a bubble-free density (`with_bubbles=False` or a `_stellar_density` underneath), and `density.py` cannot import `superbubble_field` at module level (it imports `nebula_field`, which imports `density`).
- **Carry.** A `BubbleField` handle carries the galaxy seed on the shape like `ComponentField`; `profile_at` receives the seed through it, which also removes the trap that `run_galaxy.py:152` and `:688` do not pass `galaxy_seed` today.

## 4. Stored or derived

The regions stay derived from the seed (no fixed content, Boss 4:11 and 4:14 PM EDT): the same seed gives the same bubbles on every worker and every run, nothing is written, and an in-memory index per column serves lookups. A stored copy would drift from the seed and would need a migration for nothing. The edit-everything work (8.6) can later store an admin's override of one bubble as a region override (DB.38); it is not part of this chain. DB.85 does store a row for each bubble (`galaxy_containers`), written at a plan from the same `all_bubbles` call that the field answers with, so a stored row and the field cannot disagree for one seed and shape; the rows exist to name a bubble (an object ID) and to answer which containers hold a point and what a bubble holds, and the density lookups still read the field.

## 5. Touch points

| File | What changes | Who edits it now |
|---|---|---|
| `galaxy/density.py` | one call per term-sum copy (`population_densities`, `population_densities_array`, `total_density_array`) to a new `galaxy/density_modulation.py`; the factor never enters the term sums themselves | Core 3 (GEN.319), Core 4 (GEN.134); GEN.379 (8.3) rewrites the same three functions into components, so the factor stays one call site that GEN.379 keeps |
| `galaxy/superbubble_field.py` | `box_cut`, the bound, the bubble-free read | nobody |
| `generation/object_first.py` and `generation/bright_stars.py` | fold `m_max` into `band_majorants`, multiply `wanted` | Core 3 (GEN.322, GEN.323) |
| `generation/phenomenon_scatter.py` | thin by `m(x) / m_max` | Core 3 (GEN.323) |
| `galaxy/profile.py`, `galaxy/sector.py`, `generation/run_sector.py`, `run_galaxy.py` | cell-mean count, in-cell positions, the seed passed on | Core 3, Core 4 (GEN.134) |
| `generation/run_plan.py`, `db/store.py` | attach the field to the shape like `ComponentField` | Core 1 owns migrations; no migration here |
| `galaxy/nebula_field.py` | read the bubble-free density | Core 3 (GEN.319) |

The chain waits until Core 3's GEN.322, GEN.323 and GEN.319 have merged, so no lane edits a file under another's feet; GEN.391 (design) can be done before and names the interfaces.

## 6. The drawn layer

MAP.235 (the layer, toggle and key row) and MAP.236 (clearer bounds, more transparent) draw the bubbles. Once the stars follow the bubbles, a drawn bubble is a debug and explanation layer, not the way the bubbles are shown. Decided (Boss, 9:38 PM EDT): keep it as an optional layer that is off by default and keep MAP.236, at lower priority, for when it is on.

## 7. Decided (Boss, 2026-10-16, 9:38 and 9:39 PM EDT)

1. **Young only.** A bubble changes the young and intermediate star density only; old stars, halo, bulge, clusters and satellites are unchanged.
2. **Optional overlay.** The drawn layer (MAP.235) stays, off by default; MAP.236 stays at lower priority.
3. **Stay in 8.1.** The chain (GEN.391 to GEN.395, TEST.206, DOC.61) stays in 8.1 as planned: design in phase 3, build in phase 5. It does not hold the 8.1 headline GEN.301 but it does hold the 8.1 stamp.

## 8. The table of effects (GEN.391, decided by Boss on 2026-11-26 at 7:42 PM EST) and the handle (GEN.392)

Boss chose option 1 on card A, option 1 on card B and option 2 on card C. The numbers below come out of the code (`galaxy/bubble_field.py`) and the tuning constants named; every stand-in has a row in section 6.339 of the stand-in list (GEN.725, rows SB1 to SB5).

| Population | Inside the cavity | On the shell | Outside | Where the number comes from |
|---|---|---|---|---|
| Young (0 to 100 Myr) | `interior_factor`: the stars younger than the bubble are missing, down to a floor `(gas inside / disk gas)^1.5` (0.0003 to 0.002); a bubble of the full 40 Myr is 0.6, a small young one nearly 1 | `shell_factor`: 13.1 times as many (the shell is 0.06 of the radius thick, its gas `1 / (3 x 0.06)` = 5.6 times denser, star formation that to the power 1.5) | 1 | card A (physical), card B (derived shell), [R4916] (young stars on the Local Bubble's surface), [R4347] (14 to 30 percent of massive stars form triggered) |
| Intermediate (0.1 to 3 Gyr) | 1 | 1 | 1 | Boss, 2026-10-16: a gas structure acts on young stars only |
| Old disk, bulge, halo, clusters, satellites | 1 | 1 | 1 | Boss, 2026-10-16 |

**Age from radius (card C).** A bubble's age is Weaver's law ([R5374], eq. 24 of [R5373]) solved for the age, `age = (R / 0.76)^(5/3) (rho / L)^(1/3)`, at one driving power for all bubbles (`SUPERBUBBLE_DRIVING_POWER_ERG_S` 6.3e36 erg/s, `SUPERBUBBLE_DRIVING_DENSITY_CM3` 0.5, mean mass 1.25 hydrogen masses) and capped at the 40 Myr the driving power lasts ([R5373] section 2.1). No real bubble is the anchor. The ages and inside factors:

| Radius (pc) | 50 | 75 | 100 | 150 | 250 | 400 | 500 |
|---|---|---|---|---|---|---|---|
| Age (Myr) | 1.2 | 2.4 | 3.9 | 7.6 | 17.7 | 38.8 | 40 |
| Inside factor (gas 0.005 per cm3) | 0.988 | 0.976 | 0.961 | 0.924 | 0.823 | 0.612 | 0.600 |

Half of the field's bubbles are under 160 pc (the radius is log-uniform, 50 to 500 pc), so most cavities are young and shallow; only those over 400 pc reach 0.6. The size law of the field is GEN.254's, not this chain's; the observed number of superbubbles by size is a power law (slope 1.4 to 2.1 in the luminosity function, [R5373] Table 2), so real small bubbles are likelier than the log-uniform draw has them. That is for GEN.1843.

**The factor.** `young_factor(bubble, d)` with `d` the distance over the radius: the interior factor under `d` = 0.94, rising across `0.94 +/- 0.015` to the shell's peak, flat to `d` = 1, falling across `1 +/- 0.015` to 1 (smoothstep edges, `SUPERBUBBLE_EDGE_SHARE`). It is a mix of 1, the interior factor and the peak, so it never goes beyond the peak: the peak is the certified maximum of one bubble. Overlaps multiply.

**Finding for the census.** The ball's gas goes into the shell, but star formation goes as density to the 1.5, so the shell forms more stars than the swept ball would have: the young stars of a full-age bubble and its shell are about 2.4 times the unmodulated count (0.6 of the 0.83 of the volume plus 13.1 of the 0.17). The check is the census share of young stars on a shell against 14 to 30 percent ([R4347]); if it misses, the shell's thickness moves, not the law.

**The handle.** `GalaxyShape.bubbles` is a `BubbleField` holding the galaxy's seed (`None` for a shape with no seed), set where the shape is read (`db.store.get_galaxy_shape`) and planned (`generation.plan_outline.build`). It compares and hashes by seed. `bubbles_in_box(shape, low, high)` is the box cut (every bubble whose smooth edge reaches the box, by column and draw order); `young_factor(shape, point)` is the product at a point; `max_young_factor(shape, low, high)` is the certified bound (the peak to the power of the number of bubbles that reach the box, 1 for none); `density_without_bubbles(shape, point)` is the stellar density read with no bubble in it. `superbubble_field.bubble_free(shape)` strips the handle before the bubbles' own placement reads the shape, so the factor of GEN.393 cannot feed back into where the bubbles are, and a shape with and without the handle share one cache entry. Nothing reads the factor yet: the density is unchanged until GEN.393.
