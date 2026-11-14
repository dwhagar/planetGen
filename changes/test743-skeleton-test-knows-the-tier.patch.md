**Headline:** The test of the stored-skeleton serialization knows that the Class M abundance tier is a setting for the systems made in a galaxy, not a number the Galaxy Map reads.

### Fixed
- `test_galaxy_gen::test_galaxy_density_shape_serializes_the_stored_skeleton` expected every `GalaxyShape` field in the dict the Galaxy Map reads, and failed with `KeyError: 'abundance_tier'` since GEN.1704 added the tier. `density.shape_with_terms` leaves the tier out on purpose, with the other fields that are objects or settings for the generator and not numbers for the map's density cloud, and the dict still carries the tier through `pins` (`class_m`) when the plan pinned one. The test now lists `abundance_tier` with them and asserts it is absent. No product code changes. TEST.743.

generation-output: unchanged
generation-affects: none
