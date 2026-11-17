**Headline:** `test_the_rings_under_one_solar_mass_are_the_near_masses_and_reach_less_far` judges the rings on the sample the sky really draws, not on 40 stars.

### Changed
- The test took the 40-star sample of `few_samples`. A ring's reach is where 0.5 percent of its stars still show, which is a fifth of a star in 40, so the reach was the brightest star the stream happened to draw and the heavier ring did not always reach farther (38.1 pc for the 0.3 ring against 36.9 for the 0.7 one, since the star table read the birth [Fe/H], GEN.1623). It now builds the rings with the production `RING_SAMPLES` (400 per population): 11, 21, 42, 49 and 68 pc, all assertions unchanged, no bound widened.

generation-output: unchanged
