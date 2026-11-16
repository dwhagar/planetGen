**Headline:** A restored planet or moon keeps the eccentricity a resonance holds it on, and a restored star keeps the mass it was born with and the age its phase ends, so a body or system that is settled or edited again starts from what it was made with.

### Added

- DB.132: `body_history` gains `held_eccentricity` (`_resonance_e`, the floor `settle_orbit` keeps under the tide: a planet of a resonant chain, a moon of a first-order pair or held in the evection resonance) and `forced_eccentricity` (the eccentricity a chain forces on a planet, `oddity_quantities["forced_eccentricity"]`). `body_history.attach` gives them back to a restored body, so the tide of a restored planet settled again cannot wear it below its floor.
- DB.133: `stars` gains `initial_mass_sol` (the birth mass the population model gave the star) and `phase_end_age_gy` (the age its drawn phase ends; NULL beside a birth mass reads back as infinite, a white dwarf's). A restored star keeps what `remnant_orbits` (the orbits a white dwarf's planets formed on, the engulfment radius) and `place_terms.born_mass_solar` (the companion floor) read.
- Migration 0268 adds the four columns; a galaxy made earlier has NULL in them, which reads as "not kept" and nothing is backfilled. The columns are left out of the fingerprint.

generation-output: unchanged
generation-affects: none
