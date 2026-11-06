**Headline:** The galaxy's molecular and dark clouds now store their composition too, each with its own gas, dust and phases.

### Added
- `nebula_field` gains the 28 composition columns of `nebulae` without the progenitor mass (GEN.1277, part h of GEN.1268), FLOAT, written by `planetgen plan` in the pass that draws each cloud. The design note (docs/design/composition/nebula-composition.md section 4) weighs stored floats against packed integers: about 18 MB more for 163,000 clouds, searchable without decoding.
- Migration `0215_nebula_field_composition` adds the columns and fills the rows of an older galaxy from its seed and shape.
- A field cloud's composition is keyed on its cell and place (`nebula_columns.field_key`), so the `nebulae` row it becomes when a sector stores it reads the same numbers as its `nebula_field` row.

### Changed
- A stored class M or N nebula with a centre is now keyed on its place, not its object ID, and the composition is worked out inside the galaxy's own scale (`nebula_columns.in_galaxy`), in the generator, the plan and both migrations, so a backfill reproduces what a regeneration stores.

generation-output: unchanged
