**Headline:** A nebula's shape is drawn from its class and its numbers only, so a nebula whose shape was never stored gets back the shape the generator drew, and every nebula's shape changes once.

### Changed

- **The nebula shape seed no longer holds the species text (GEN.1332).** `shape_for_nebula` hashed the class letter, four numbers and `dominant_species`. Since GEN.1276 a save rewrites that species text from the nebula's composition, after the shape is stored, so a row with no stored shape (saved before shapes were kept) drew a different shape from the one a fresh save stores. The seed is now the class letter and the radius, density, temperature and extinction only (`shape_for_nebula` loses its sixth argument; the generator and `query.nebula_shape` pass five). Every nebula gets another shape, drawn the same way as before: a lock item, with the epoch and the golden digests added at the train. `test_a_nebula_saved_before_shapes_gets_the_shape_a_fresh_save_would` is back to comparing with the generator's shape and runs three nebulae; a new test holds the seed to the numbers.

generation-output: changed
