**Headline:** the mass of Jupiter and of the Sun are written once, so the moons, rings, giants and rogue-planet classes no longer disagree by a part in a million.

### Changed

- GEN.1152: moon_composition's `JUPITER_MASS_KG` (1.89813e27), the star mass of ring_origin's Saturn reference (1.989e30), giant_interior's `EARTH_MASS_JUPITER` (317.83), stellar_abundances' `JUPITER_MASS_SOL` and the rogue-planet mass bins of `tuning` (317.8) now read `physics/constants.py` (`JUPITER_MASS_TO_KG`, `SOLAR_MASS_TO_KG`, and the new `JUPITER_MASS_EARTH`, their ratio of the two conversion factors). A test fails on any of the rounded copies written again. The values move by 2.5e-6 (Jupiter in kg), 2e-4 (the Sun in kg) and 1e-4 (the rogue bins' edge), so a number computed from them can move in its last digits and a draw that lands within that of a threshold can fall on the other side; the pinned golden galaxy is unchanged.
generation-output: changed
generation-affects: fill
