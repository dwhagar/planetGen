**Headline:** The class census script draws systems again: it read a pair attribute that the pair-orbit work removed, so every system it drew was an error.

### Fixed
- GEN.1879 (bug): `scripts/class_census.py` wrote `system.wide_binary` into each system's record, and GEN.640 slice b removed that attribute (a wide pair's orbit is its node's, read as `system.pair_orbit`), so every system of a census failed (25 of 25 in the first minute of a default census). The record now carries `pair_orbit` (separation, eccentricity, period, periapsis and apoapsis of the pair, or none for a single star) and no `wide_binary` key; no table of the report read the old key. `tests/test_class_census_script.py` draws a forced wide pair and a single star through `draw_one`, and its fixed-seed list test now counts the `mdwarf` sample (300 systems) that the script gained.

generation-affects: none
