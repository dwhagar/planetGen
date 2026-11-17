**Headline:** the census audit fills its sectors with four workers at a time, so the whole-sphere tests no longer run past the hang guard.

### Fixed

- TEST.783: `tests.audit.representation.run` filled every seeded sector of every region one after another, and `test_census_component_wholesphere` (30 sectors in each region) ran past the hang guard (1,500 s in a gate). It now fills the sectors with forked workers (`workers`, up to 4 by default, 1 for the old one-by-one run). Each sector is bound to its own seed, so the census and the star counts of the fill are the same; a new test in `test_census_representation.py` compares a run with 1 worker and one with 3 (every sector's counts and the fill's star and sector counts). `test_census_representation` uses the same `run` and gets the same speedup.
