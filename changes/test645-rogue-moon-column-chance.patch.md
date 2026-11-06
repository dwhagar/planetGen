**Headline:** The column round-trip test no longer fails now and then on whether the test galaxy's rogue planets hold a moon.

### Fixed
- TEST.645: `test_every_column_is_written` listed `moons.rogue_planet_id` as never written in the test galaxy, but whether a rogue planet that galaxy's donors threw out holds a moon depends on the draws of the build, so the test failed in some combined runs ("written after all"). The column is now listed as chance, like the other columns the draw decides; `test_parent_link` and `test_rogue_body_storage` still write it and read it back.
