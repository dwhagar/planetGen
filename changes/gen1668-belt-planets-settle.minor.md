**Headline:** A belt member made a planet gets the same finishing pass as any planet: its tilt range, the check for a chaotic tilt, its moons' tilts and the tides of its orbit.

### Changed

- Belt members promoted to planets (GEN.1510) are made after the first settle pass, so they missed it. They now go through `StarSystem._settle_promoted`: the sibling planets are registered as tide sources, the planet and its moons go through the same per-body settle as any planet (`_settle_body`, the orbit's tides, envelope and end stages), and `settle_obliquity` (GEN.1454, GEN.1535, GEN.1536) and `settle_moon_tilts` (GEN.1556) judge them against the modes of all the star's planets (new `only=` argument). A free-spinning promoted planet now has a tilt range (and, where its tilt lay in a chaotic zone, a new tilt); a locked one takes its Cassini state. The secular exchange and the Kozai pumping are not run again.
- The property sheet no longer says a promoted planet has no tilt range.

generation-output: changed
generation-affects: fill
