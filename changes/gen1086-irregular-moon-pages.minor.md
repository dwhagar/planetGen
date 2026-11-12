**Headline:** A planet that caught small bodies now has a page listing its irregular moons and a page for each one, with its size, makeup, orbit and how it was caught.

### Changed

- GEN.1086 (slice 1): a planet that keeps irregular moons shows a link on its page to `/object/planet:<ID>/irregular-moons`, the list of its bodies by designation with mass, diameter, type, makeup, distance, share of the planet's Hill radius, eccentricity, inclination, direction (forward or backward) and the window it was caught through, and what was not kept one by one. Each body has its own page, `/object/planet:<ID>/irregular/<rank>`, with a Sensor view, its orbit about the planet (closest and farthest distance, the period the planet's mass gives, the six elements) and links to the next heavier and lighter body. A body is called `<planet ID>.<rank hex>`, or by the belt's designation when it was caught from an inner belt.
- New routes `GET /api/objects/<ref>/irregular-moons` and `.../irregular-moons/<rank>` (`describe/irregular_moon_info.py`); the API version stays.
- Not yet: the irregular orbits on the System Map and the 3D view, and the clearance ledger taking the irregular zone as occupied space (the rest of GEN.1086).

generation-output: unchanged
