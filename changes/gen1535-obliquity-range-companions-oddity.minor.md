**Headline:** A planet's tilt now feels the pull of its star's companion stars, each planet keeps the range its tilt wanders over, and a planet whose tilt wanders chaotically is named for it on its page.

### Changed

- GEN.1535 (GEN.1454 part 3, a lock item): the companion stars of a planet's star now drive a nodal mode in its orbit (`planets.companion_modes`, built on `obliquity_variation.companion_mode`) from the plane each companion moves in, which `StarSystem` now works out from two places of it a tenth of its period apart and keeps on `radiation_field.Companion.normal`. A planet in a wide or close binary, alone or among others, is set against these modes with its siblings' when `settle_obliquity` looks for a chaotic zone.
- `planets` and `moons` and `rogue_planets` gain `obliquity_min_deg` and `obliquity_max_deg` (migration 0244): the lowest and highest tilt a free-spinning planet's axis reaches over the age of its star, for the seasonal swing to read; NULL for a moon, a rogue planet, a locked or resonant body and rows made before the migration. The property sheet shows them.
- Oddities W6 (planet) and W61 (giant), "Wandering tilt": a planet that `settle_obliquity` put in a chaotic zone is named, with no extra draw (rate 1, registry OB7). Registry rows OB7 and OB8 added; OB5 and OB6 updated.
- Not done here: the Cassini state of a locked body from the real precession constant. The fluid Darwin-Radau figure gives a Mercury-like body a precession constant of 5 arcseconds a year (the real one must be some hundreds of times larger, because its figure is frozen in), so its tilt would come out at 4 to 7 degrees instead of 0.03; the log-uniform stand-in stays until the figure (J2 and C22) is stored.

generation-output: changed
generation-affects: fill
