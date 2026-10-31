**Headline:** The old stars of the disc are as metal-poor as their longer age range says: the thick disc's 8 to 13 billion year old stars pull the old population's mean [Fe/H] down by about 0.125 dex.

### Changed
- GEN.1025 follow-up (train 7): the old population's ages now run 3 to 13 Gyr, and the age-metallicity relation's mix over that window averages -0.255 against -0.130 over 3 to 10 Gyr. `tuning.METALLICITY_BY_POPULATION["old"]["offset"]` was -0.2 and is -0.325, so the population's [Fe/H] is still the relation's mix over its window and the age tilt (`age_metallicity.tilted_age`) leaves the population's mean age where it was (the mean age had fallen from 8.04 to 7.14 Gyr). Stand-in tier T6, tweak L2 (tuned to the relation). `test_star_matrix` now reads the main-sequence lifespan on the track its alpha branch stands for (GEN.1024), as the generator does.

generation-output: changed
