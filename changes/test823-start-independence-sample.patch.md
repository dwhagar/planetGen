**Headline:** The air solver's start-independence test counts 24 generated systems, not 12, and the train gate now runs it when generation or physics code changes.

### Fixed
- The test of how many air bodies have a start-dependent answer (limit 40 percent) sampled only 19 to 22 air bodies, so one or two bodies swung the share by 5 to 10 points. It went 7 of 22 (32 percent) to 8 of 19 (42 percent) at Release 4717 (GEN.1493 slice 3 cut the veneer: cold outer bodies of about 100 K, all single-answer, now hold under 10 kPa and leave the count). On 40 seeds it is 24 percent before and 28 percent on main. The sample is now 24 seeds (13 of 40, 32 percent on main) and the limit is unchanged.
- The train gate runs `tests/test_surface_air_start_independence.py` whenever a file under `generation/` or `physics/` changes (it was left out of the sample tests, so the Train 56 and 57 lists did not hold it).
