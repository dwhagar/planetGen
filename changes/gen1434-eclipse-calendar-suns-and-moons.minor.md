**Headline:** A moon's page now says how often its planet hides its star, and a planet's page how often its moons hide the other suns of its sky.

### Added

- GEN.1434: the eclipse calendar counts the rest of the sky. `sky_clock.companion_eclipse` gives the eclipses of another sun of a planet's sky (GEN.1443) by each of the planet's moons: the sun stands still at its latitude off the orbit plane, the moon's track rises and falls with its tilt, and the mean interval, the longest eclipse and the deepest cover follow (`touching_share` takes the source's offset, `overlap_share` the hidden part of a partly covered disk). A planet's sheet gains "Eclipses of its second sun by <moon>" rows (`sun_eclipses_1` and on, the shortest interval first, at most four); a moon's sheet gains "Eclipses of its star by its planet" (`moon_eclipse` with the planet as the blocker).
- GEN.725: rows SK4 (the other suns), SK5 (the planet as blocker) and SK6 (the line-ups the calendar leaves out: sibling moons, a planet by a sibling planet, another sun by an inner planet) of section 6.272.

Derived on demand from stored masses, radii and orbits; nothing is stored and no number already shown moves.

generation-output: unchanged
generation-affects: none
