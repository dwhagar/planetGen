**Headline:** Rock weathering on a world now counts only the land its sea waters, so dry worlds and worlds with a small sea keep more carbon dioxide and run warmer.

### Changed

- GEN.1403 (Boss, 7:25 PM EST 11-17, "Build it"; the form is Lane 14's PR 6062 rebuilt on the runoff weathering law): `surface_air._weathering_land` scales the land with `surface_pressure.wet_land_fraction`: the land times the share of the Earth's ocean cover the world has, `min(1, ocean_cover / 0.71)`. A world with the Earth's sea or more is unchanged; a dry world's land, which has no sea to rain on it, weathers less and keeps more of the carbon dioxide its volcanoes gave it. Stand-ins WW1 and WW2 (educated guess from physics, L1) in section 6.309.
- Changes the carbon dioxide, and with it the temperature and class, of every world with land and less than the Earth's ocean cover, so it is a lock item.

generation-output: changed
generation-affects: fill
