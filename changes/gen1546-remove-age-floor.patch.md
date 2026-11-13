**Headline:** A star's age is no longer raised to suit the planets it hosts, so planets are made on the age the star was born with and generation is about a quarter faster.

### Changed

- GEN.1546 (with GEN.1545): the age floor is removed. `Star.adjust_age_for_planets`, `adjust_pair_age_for_planets`, the binary proxy's `adjust_age_for_planets`, `planet_age_requirement_gy`, `age_ceiling_gy`, `physics/age_needs.py` and `age_needs.json`, the tuning constants only they read (`MAX_PLANET_AGE_ADJUSTMENT_FACTOR`, `UNREACHABLE_PLANET_AGE_MIN_LIFESPAN_RATIO`, `WHITE_DWARF_AGE_ADDITION_GY`), the post-age trim and belt rebuild, and the second planet round under the settled age (`system.planets_regenerated_after_age`) are gone. A pair's stars share the age they are built with. The habitable-request tilt in `sample_living_star` and the life timeline caps stay. Measured by the star-last research on 100 default seeds: the same planet and rocky counts, bodies with air 219 to 216, class mix nearly unchanged, run time down 28 percent, M-dwarf median age 6.8 to 4.2 Gy, no pile of ages at 13.8 Gy.

generation-output: changed
generation-affects: fill
