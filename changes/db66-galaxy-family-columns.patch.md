**Headline:** The galaxy shape row can now hold a galaxy's family and the numbers of its bar, thick disk and spheroid.

### Added
- `galaxy_shape` gains `family` and fourteen component columns (the bar's axes and angle, the thick disk's amplitude, the Sersic spheroid, the irregular family's clump share) in migration 0117. Every existing galaxy becomes a `spiral` with no spheroid and no clump share, so no stored density changes. Nothing reads the columns yet (GEN.379 does).
