**Headline:** The roundness engine now has strengths for the other ices of an icy body (high-pressure water ices, methane hydrate, dry ice, ammonia dihydrate) and a way to name a body's ice type from its composition.

### Changed
- **GEN.1177.** `physics/data/roundness_strengths.json` gains an `ice_types` section: ice II and VI as strong as ice Ih, V a bit lower (0.5), III near zero, methane clathrate hydrate 25 times ice (20 to 30), solid CO2 0.3 MPa (0.1 to 0.5), ammonia dihydrate 2 times ice at 110 K, 1 at 135 K and 0.3 at 170 K. The new `physics/ice_types.py` reads it: `strength_ratio`, `strength_pa`, `water_phase` (the phase of ice at a central pressure) and `ice_type_of` (names the ice type of a body from the components of `small_body_composition`). Nine registry rows (6.129, IT1 to IT9). Nothing reads the types yet, so generated worlds do not change.
generation-output: unchanged
