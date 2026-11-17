**Headline:** A system build is about a tenth faster: the tide sources no longer work out the same orbit means again for every body and every pass.

### Changed
- PERF.276: `tide_sources.orbit_mean_inverse_cube` and `orbit_mean_inverse_sixth` (a mean over the points of the source's orbit) are cached; they are pure in their three numbers and a build asked for the same pairs again and again (20,095 calls for 1,364 different pairs in three G2V systems, 9.3 s of 90.6 profiled seconds). `planet_sources` fell from 9.3 s to 0.8 s over those three systems, the build from 90.6 s to 83.1 s profiled (45.5 s to 42.1 s plain for four systems). A digest of every body's mass, distance, eccentricity, day, tilt and temperature, and each moon's, is the same before and after.

generation-output: unchanged
