**Headline:** A sea of so little water that Earth's salt would saturate it is now a rock-salt sea of real depth, a surface ocean with a water activity of 0.753, and not a puddle of a brine too salty for life.

### Changed

- GEN.1730 (every galaxy): surface water under `hydrosphere.BRINE_MAX_WATER_FRACTION` (2e-5 of the planet's mass, where Earth's salt per planet mass just saturates) is a halite-saturated sea where it is liquid above the sodium chloride eutectic of 252 K: water activity 0.753 (`NACL_SATURATION_ACTIVITY`, a standard value) and the ice melts at 252 K, where it was drawn from 0.30 to 0.65 and a perchlorate eutectic of 198 to 236 K. Such a sea 10 m deep or more (`SEA_MIN_DEPTH_KM`) is a `surface ocean`, not a `brine seep`: the oxygen law counts it as a sunlit sea, and the life law's solvent term at 0.753 is no longer near zero. Colder than 252 K, or shallower, nothing changes: the perchlorate brine and the seep stay. On the 15 desert worlds measured (5 to 20 percent of the surface wet, water mass 4.4e-6 to 2.9e-5) the sea is 0.6 to 2.8 km deep, and 12 of 15 hold more salt than the water dissolves.
- The Class M tier's own rule for a desert world's small sea (`abundance_tier.desert_sea`, an ordinary salt sea of water activity 0.97 to 0.99) is gone: the same law now makes every such sea, with water activity 0.753. `surface_liquid_world.NACL_SATURATION_ACTIVITY` is now `hydrosphere.NACL_SATURATION_ACTIVITY`.

generation-output: changed
generation-affects: fill
