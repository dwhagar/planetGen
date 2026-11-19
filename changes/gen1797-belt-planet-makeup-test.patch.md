**Headline:** The promoted-belt-planet test compares a body's makeup with its member's after the body's own hydrogen air is taken out of the mass.

### Changed

- GEN.1797 (test only): `test_a_promoted_planet_has_the_orbit_and_makeup_of_its_member` held each mass fraction of a promoted belt planet to its member's share to 1e-6, but a cold small body can hold a secondary hydrogen air (GEN.1226) of a few millionths of its mass, and the other fractions make room for it (rock 0.75 x (1 - 2.9e-6) = 0.7499978). The seeds the test reads draw a body with such an air now that the architecture draws its own eccentricities; the expected fraction is now the member's share times (1 - the envelope fraction), at 1e-9, a tighter match than before.
