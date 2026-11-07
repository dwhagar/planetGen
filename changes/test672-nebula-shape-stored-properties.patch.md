**Headline:** The legacy-nebula shape test no longer fails on a nebula whose species text a save rewrites; it checks the shape its stored properties give.

### Changed

- **Nebula shape test follows the stored species (TEST.672).** A nebula's shape is seeded from its properties, `dominant_species` among them, and since GEN.1276 a save rewrites that text from the numbers after the shape is stored. `test_a_nebula_saved_before_shapes_gets_the_shape_a_fresh_save_would` compared the shape of a row with no stored shape against the generator's own, which can never match; it failed on every run, with a different nebula each time because the test builds a random one. It now checks the shape equals the one the stored properties give and is the same on a second read (renamed `..._gets_the_shape_its_stored_properties_give`). The `query.nebula_shape` docstring says the same; no behavior changes.

generation-output: unchanged
