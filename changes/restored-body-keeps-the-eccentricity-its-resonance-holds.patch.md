**Headline:** A saved body that a resonance holds on an eccentricity keeps it when the galaxy is loaded again.

### Fixed
- GEN.1786: since the build record (GEN.1651) the generators keep the eccentricity a resonance holds a body on in `body.build.resonance_e`, but `body_history.numbers_of` still read an old `_resonance_e` attribute, so a generated body never had its `held_eccentricity` saved, and `body_history.attach` wrote the old attribute back, which `planets.settle_orbit` no longer reads, so a restored resonant body lost the floor under its eccentricity when it was settled again. Both now use `build.resonance_e`. Only the saved `held_eccentricity` column and a restored body settled again are touched; a body made fresh comes out as before.
- `test_restored_resonance_and_star_birth` searched its seeds for the old `_orbit_birth` attribute, so none of 80 suited it; it and `test_body_history` now read `build.orbit_birth` and `build.resonance_e`. The assertions are as they were.

generation-affects: none

generation-output: unchanged
