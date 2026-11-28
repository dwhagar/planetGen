**Headline:** Tests now check that every galaxy the generator can draw falls inside the sourced ranges for its family.

### Added
- **TEST.204: family plausibility tests.** `test_family_plausibility.py` draws 2,500 seeded spiral galaxies across four stellar mass bins (3e9 to 2e11) and fails if a parameter falls outside its range (stellar mass, disk length against the size-mass law, bulge-to-disk by stage, bar length, arm count, pitch angle, bulge-to-total, bulge radius), if the size-mass slope or the median size residual is off, if the halo or the black hole drifts from `median_halo_mass` or the bulge relation, or if the Milky Way or Andromeda cannot be reached. A dial sweep pins each lever at its floor, middle and ceiling and checks the result is finite or refused, stays inside the published bands, and that the GEN.625 warnings fire only where they should. Only the spiral family is built (GEN.380 to 385 and GEN.389 are open), so only it is tested.
- The thin disk scale height has no row yet: with a mass pinned, 10 of 2,500 seeds (0.4 percent) draw a height above the 1,155 pc ceiling in `DERIVED_LIMITS`. No tolerance was added; GEN.1902 fixes the draw and adds the row.

generation-affects: none
generation-output: unchanged
