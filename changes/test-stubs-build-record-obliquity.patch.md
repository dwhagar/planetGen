**Headline:** The stand-in moon of the Cassini-tilt test again carries a build record, so main is green again after GEN.1651.

### Fixed
- `test_settle_obliquity` made a stand-in moon without a `build` record, which `settle_moon_tilts` now reads for the moon's frame axis (GEN.1651); the stand-in carries a `BuildRecord()`.

generation-output: unchanged
