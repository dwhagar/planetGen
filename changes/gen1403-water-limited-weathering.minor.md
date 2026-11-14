**Headline:** Rock weathering on a world now counts only the land its sea waters, so dry worlds and worlds with a small sea keep more carbon dioxide and run warmer.

### Changed

- The continental silicate weathering that holds a world's carbon dioxide (`surface_pressure.carbonate_pco2_bar`, the law of Lehmer et al.) scaled with the land share of the surface. It now scales with `surface_pressure.wet_land_fraction`: the land times the share of the Earth's ocean cover the world has, `min(1, ocean_cover / 0.71)`. A world with the Earth's sea or more is unchanged; a dry world's land, which has no sea to rain on it, weathers nothing and keeps all the carbon dioxide its volcanoes gave it. Stand-ins WW1 and WW2 (educated guess from physics, L1) in section 6.287.
- This is the lever 4 of the Class M rarity report (the GEN.909 option 3). The temperate share on the same 200 systems before and after is in `reports/gen1403-water-limited-weathering-2026-11-11.md`.
- Changes the carbon dioxide, and with it the temperature and class, of every world with land and less than the Earth's ocean cover, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
