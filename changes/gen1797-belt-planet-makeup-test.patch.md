**Headline:** The promoted-belt-planet test compares a body's makeup with its member's after the body's own hydrogen air is taken out of the mass.

### Changed

- GEN.1797 (test only): `test_a_promoted_planet_has_the_orbit_and_makeup_of_its_member` held each mass fraction of a promoted belt planet to its member's share to 1e-6, but a cold small body can hold a secondary hydrogen air (GEN.1226) of a few millionths of its mass, and the other fractions make room for it (rock 0.75 x (1 - 2.9e-6) = 0.7499978). The seeds the test reads draw a body with such an air now that the architecture draws its own eccentricities; the expected fraction is now the member's share times (1 - the envelope fraction), at 1e-9, a tighter match than before.

- GEN.1797: the eccentricity a planet's slot of the architecture drew (`plan_eccentricity`) is now a field of the planet's build record (`BuildRecord.plan_eccentricity`, runtime state never stored) and not an attribute of the planet, which `test_bughunt_serialization_fuzz` found unlisted in `SERIALIZABLE_FIELDS`. Nothing is stored, so no column, migration or golden digest changes, and no generated value changes.

- GEN.1797 (test only): `test_a_nodes_planets_are_drawn_from_a_stream_of_its_own` compared the regenerated node planets with every stored planet, but the stored list holds the belts' promoted members (GEN.1510, `belt_rank` set), which come with the belts after the draw it regenerates; the sample the architecture's new draws give has ten of them. The stored side leaves them out.
