**Headline:** A test that never ran now no longer hides behind its twin: the duplicate pulsar-planet isotope test is deleted.

### Fixed
- **TEST.732: `test_isotope_record.py` defined `test_a_pulsar_planet_with_no_orbit_distance_takes_the_snow_lines_ice_and_still_gets_a_record` twice with the same body, so Python kept the second and the first never ran.** The first copy is deleted. A scan of every `def test_` in `src/tests` (module level and inside classes) finds no other name defined twice in one scope. Tests only; no generated output moves.

generation-affects: none
