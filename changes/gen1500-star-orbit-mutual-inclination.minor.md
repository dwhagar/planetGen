**Headline:** The orbits of the stars of a binary or triple now line up when the stars are close and point anywhere when they are far apart.

### Changed

- Every companion's orbit used to be drawn with a random tilt (the whole sphere, at every separation). It is now drawn from two populations (`physics/star_orbit_alignment.py`): aligned, a small prograde tilt (a Rayleigh law of 5 degrees), and random, the sphere as before. The share that is aligned is all of them at 30 AU and inside, none at 1,000 AU and beyond, and falls with the logarithm of the separation between. The separation is the companion's own axis for a star's first companion, whose tilt is the tilt against the plane of that star's planets, and the axis of the pair it orbits for the inner pair of a triple, whose tilt is against that outer orbit.
- The tilt is drawn as the angle to the reference plane and turned about that plane's normal into the stored inclination and ascending node, so no column changes: the mutual inclination is still derived from the stored values (`stability.mutual_inclination`), and the triple stability cut, the Kozai-Lidov verdict and the quadruple verdicts read it as before.
- The stand-ins (the 30 and 1,000 AU edges, the 5 degree scale) are registry section 6.258 (AL1 to AL4). A first-companion tilt, which sets how a companion's plane sits against its star's planets, is drawn after the axis, so the draws of a multiple shift from the first companion on.
- Readers of the new tilt (the planet Kozai screen of GEN.1447, the forced inclination of GEN.1492, the precession of GEN.1467) are later slices and are not in this change.

generation-output: changed
generation-affects: scatter, fill
