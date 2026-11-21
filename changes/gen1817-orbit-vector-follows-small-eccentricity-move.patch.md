**Headline:** A planet whose forced eccentricity moves by less than the settle tolerance now has its position and velocity made on the new eccentricity, so the orbit they describe is the one the eccentricity names.

### Fixed

- `planets.finish_orbit_oddities` returned early for an eccentricity move under `ECCENTRICITY_SETTLE_TOLERANCE` (1e-4), after the eccentricity had already changed, so the position, the velocity, the moons' host eccentricity and the light stayed on the old ellipse. A chain planet held by oddity W5 (Efriggex III of seed 21, 12th system of `test_every_generated_orbit_is_the_ellipse_it_was_drawn_on`) showed it: eccentricity 0.010461925 against 0.010461783 from `orbit_from_vector`, 1.4e-7 against the test's 1e-7. A move under the tolerance now runs `follow_orbit_change` (position, host eccentricity, exposure) and leaves the air and the class as they were; no move changes nothing (GEN.1817).

generation-output: changed
generation-affects: fill
