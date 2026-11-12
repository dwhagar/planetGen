**Headline:** the ammonia of an ice-covered ocean is now concentrated in the liquid left under its ice shell, so more icy ocean worlds hold the strong ammonia-water solution that ammonia-water life and class X ask for.

### Changed

- `physics/surface_liquid_world.py` (GEN.1563, slice 1): `_solutes` multiplies the ocean's ammonia share of the water by `(liquid + ice shell) / liquid` (`_liquid_concentration`) and holds the stored percent at the eutectic (`AMMONIA_WATER_MASS_PERCENT`). The ice shell froze out of the ocean and holds almost no ammonia, so the liquid keeps it; the high-pressure ice below keeps its own (counting it in raised nearly every ammonia ocean to the eutectic). Registry LQ16.
- Measured on 180 seeds: 76 bodies hold ammonia in an ocean; those at or over 0.1 of the water go from 2 to 18 (classes E 1 to 12, C 0 to 5, W 1 to 1). Class X needs an airless ocean at 190 to 240 K: none in the sample either way (the airless ammonia oceans are colder, 83 to 134 K), so class X stays rare; ammonia-water life (class 8) is reached more often by icy ocean worlds with an air.

generation-output: changed
generation-affects: fill
