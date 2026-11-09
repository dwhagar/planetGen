**Headline:** The moon-history test accepts a dust ring as the origin of the ring of a planet that lost a moon, as the ring code has always allowed.

### Changed

- `tests/test_moon_history.py::test_every_planet_with_moons_ends_with_moons_it_can_keep` asserted that a planet with fallen moons and a ring had a ring of origin `disrupted_moon` or `captured_body`. The origin is the largest event of the ring (`ring_origin.ring_origin`), and the steady dust ring of the small moons near the planet (GEN.1233, `moon_erosion`) is one of the events: a fallen moon's ring spreads thin with the time since and a faint dust ring can outweigh it. The Train 21 merge changed what the 16 seeds draw, and one planet now has a dusty ring of 1e8 kg after a fallen moon. The test accepts the three origins (`rings.ORIGINS`) and requires a `moon_erosion` ring to be a dusty one. The ring code is unchanged.

generation-affects: none
