**Headline:** The young-star planet test checks that no life is drawn, not which class letter a planet gets.

### Changed
- TEST.462: `test_planets_respect_their_stars_age_and_history` failed now and then because it required that no planet of a star under 0.1 Gy carry a class letter from `HABITABLE_PLANET_CLASSES`. Since the class cut-over (GEN.646) the letter is derived from the finished conditions and says nothing of life (GEN.915), so a class V world around a 0.08 Gy star is right as long as it holds no life. The test now asserts that these planets and their moons have no life class. The test draws unseeded systems, which is why the failure showed only in some runs and orders.
