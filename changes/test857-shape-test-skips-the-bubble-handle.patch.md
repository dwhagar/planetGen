**Headline:** The galaxy shape test no longer expects the superbubbles' handle in the Galaxy Map's shape dict.

### Fixed

- TEST.857: `test_galaxy_density_shape_serializes_the_stored_skeleton` loops the fields of `GalaxyShape` and read `shape_dict["bubbles"]`, which GEN.392 added as the superbubbles' handle. The product is right: `density.shape_with_terms` already leaves `bubbles` out of the dict the map reads (the handle is an object for the density model and is not in the stored row), so the test lists `bubbles` with the other object fields.

generation-output: unchanged
