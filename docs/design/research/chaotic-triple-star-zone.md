# Design note: the chaotic zone of a triple star (GEN.1437, part 1)

Written by Lane 21, 2026-11-11, from section 5.3 and item A5 of `/mnt/project-files/reports/flavor-text-and-sci-fi-anomalies-2026-11-08.md`. Design only: nothing here is built. Every number below that is not tied to a reference is a proposal and a stand-in to be registered (GEN.725) when it is built.

## 1. What the fiction and the physics agree on

A planet around one star of a bound triple (the "three suns" of the 3 Body Problem) really can have an orbit that is chaotic without being lost: it wanders, its eccentricity and its flux swing in irregular steps, and it can be ejected after a long quiet time. The formal idea is a split between Lagrange stability (nothing escapes) and Lyapunov stability (nearby orbits stay nearby); the region between them is bounded and chaotic ([R1855], Hayashi and others 2023, who study it for hierarchical triples; I have read the title only). For a planet around one star of a pair, the empirical edge of the safe region is the Holman and Wiegert limit ([R1990], equation 1) and, outside it, the chaotic band whose width a mass ratio sets (the resonance-overlap criterion of Wisdom 1980, [R4727], applied to planets by [R3226] and [R3175]; for a star as the perturber the same idea holds but the width is not fitted here).

## 2. What the generator already has

- **Stars:** a triple that fails the stability criterion is not kept. `physics/stability.py` (`is_triple_stable`, `auto` picks Equation 4, Mardling and Aarseth [R2868] or Eggleton and Kiseleva) and, for four stars, `physics/hierarchy_stability.py` and `physics/quadruple_stability.py` with a margin on the axis ratio (`level_margin`, `axis_margin`). So no generated triple is deep in the unstable regime. The band just inside the margin is where an orbit is marginal.
- **Planets:** the planet's outer limit near a companion is the Holman and Wiegert S-type limit (`physics/orbits.holman_wiegert_critical_semimajor_axis`, with the P-type for circumbinary, `holman_wiegert_circumbinary_a_crit_au`); the verifier `physics/clearance.py` checks it with the same fit; the disk is cut by the companion (`physics/disk_truncation.py`, GEN.1451).
- **Light:** `physics/radiation_field.py` sums the light of every star at a body, and `insolation.companion_swing` and `flux_variation` store the swing of a planet's starlight for one companion on an eccentric orbit. A third star is not in that swing yet (the swing takes the "other" luminosity as one sum, `other = sum(luminosity of the stars beside)`).
- **Loss:** a planet ejected by scattering has a kick and a destination (`physics/rogue_kick.py`, `rogue_age.py`), so "lost after an era" can end in a rogue the galaxy already draws.

## 3. The design: a derived chaotic zone, not a new draw

A planet is in a **chaotic zone** when its orbit is bounded but close enough to a critical limit that the companion stars' pull is not a small perturbation. It is an output of the finished orbits, like the class: nothing is rolled, no game rate.

1. **Position in the zone.** For a planet round star A with the nearest companion B (or a pair B, C seen as one body when they are a close pair): `x = a_planet / a_crit`, where `a_crit` is the Holman and Wiegert limit for the mass ratio and the eccentricity of B's orbit (for a P-type orbit the ratio is inverted: `a_crit / a_planet`). The zone is `X_LOW <= x <= 1`; the generator keeps no orbit beyond 1 (a stand-in `X_LOW = 0.6` to start, tier from inferred physics, tweak L3, falsified by the N-body samples below). A planet near the margin of a close triple is also in the zone when the triple's own margin (`hierarchy_stability`) is below `TRIPLE_ZONE_MARGIN` of the criterion (stand-in 1.2).
2. **The flux swing from every star.** The swing is computed from the three-dimensional geometry of all stars over one outer orbit (`radiation_field` already gives the field at a position; the missing step is the time series: sample the planet and the stars over the longest orbit and keep the minimum and maximum of the bolometric flux). For a zone planet the swing is the largest of the sampled values, and the existing flag `flux_variation` is read as it is.
3. **The era length.** An orbit in the zone has a time over which it keeps its character: the Lyapunov time, of the order of a few outer orbits near the stability edge ([R2868] gives the escape time scalings; I have not read its exponent). Stored as `chaotic_era_yr` (a number of years, or empty outside the zone): `ERA_ORBITS` times the longest orbital period of the triple, with `ERA_ORBITS` falling from 100 at `x = X_LOW` to 3 at `x = 1` (stand-in, tier from inferred physics, tweak L3). It is not a prediction of one orbit, it says how often the sky and the seasons of that planet reset.
4. **What the planet can still be.** A zone planet keeps its class, derived from the finished body, but its temperature class is judged on the mean and its eras on the swing: a habitable-looking mean with a swing over some factor is a world whose life lines name the eras (a gated flavor line, no class letter).
5. **The end of the planet.** A zone planet older than a few eras is unlikely to be there: the check is that the planet survives its age (`planet_survival`, the age-dependent ceiling GEN.325) with the era as the survival time, and a planet that fails is the ejected rogue the system's `rogue_kick` makes.

## 4. How we will check it

A REBOUND sample like the 2+2 samples of GEN.1155 (`scripts/build_quadruple_stability.py`'s checks): triples at x from 0.4 to 1.2 for planets of 1 and 10 Earth masses, 10^4 outer orbits, the fraction ejected, and the time of the first excursion past a factor of 2 in the flux. This sets `X_LOW`, the era law and the survival rule, and replaces the stand-ins. The tests: a zone planet's `x` is inside 0.6 to 1; the era falls as `x` rises; no zone planet is flagged by a lone star or a wide pair; the swing is at least the existing `flux_variation` for the same system.

## 5. Chain items this needs (not built here; TODO numbers them)

1. **Pure functions** `chaotic_zone.position(...)` and `chaotic_zone.era_years(...)` with tests on the Holman and Wiegert cases (Kepler-16 b, Alpha Centauri) and a three-star case. Lock-free.
2. **N-body calibration** of `X_LOW`, the era law and the survival rule, as a report and the three stand-in rows. Lock-free.
3. **The time-series flux swing** over three stars (extends `insolation.companion_swing` beyond one companion and `radiation_field` to a sampled orbit). Lock: changes stored `flux_variation` of triple systems.
4. **Stored fields** `chaotic_zone` and `chaotic_era_yr` on planets, with columns, sheet rows and the unit registry entry (an epoch and a migration).
5. **The oddity row** for the zone (a gate on the stored `chaotic_zone`, no extra rate), the gated flavor lines for the three-sun eras, and the plausibility limit (a zone planet keeps its mean temperature class in range).
