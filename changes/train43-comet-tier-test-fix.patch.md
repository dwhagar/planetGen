**Headline:** A test that reads the stored galaxy shape back now skips the comet tier pin the way it skips the Class M tier pin.

### Fixed
- **`test_galaxy_density_shape_serializes_the_stored_skeleton` read a field the map's shape does not carry (GEN.1775).** `GalaxyShape.comet_tier` is a pin kept for the systems made here and is dropped from the map's JSON like `abundance_tier`; the test now skips it too. Test only.

generation-affects: none

generation-output: unchanged
