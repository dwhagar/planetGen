**Headline:** Three tests of a galaxy's central black hole no longer fail by chance on a small draw.

### Fixed
- **`test_quasar_physics_follows_from_mass_and_eddington_ratio` failed once in a gate with a mass of 80,504 solar masses (TEST.286).** The test asserted that a nucleus weighs at least a tenth of the supermassive minimum (1e5), but the galaxy's own black hole follows its bulge or stellar mass with a log-normal scatter (GEN.402), so 38 percent of the draws are under the minimum and about one in 400 under a tenth of it; 50 unseeded draws failed one run in eight. The bound was never a property of the draw. The test now asks only that the mass is a finite positive number, and the two single-draw tests with the same bound (`test_supermassive_black_hole_sits_still_and_glows_faintly`, `test_a_quiescent_nucleus_is_stored_as_a_supermassive_black_hole_at_the_center`, each failing about one run in 500) are changed the same way.

TODO items resolved: TEST.286. Added: none.
