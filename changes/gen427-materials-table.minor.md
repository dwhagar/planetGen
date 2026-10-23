**Headline:** One shared table now holds every material a body can be made of, with its density, phase points and equation-of-state constants and the source of each number.

### Added
- **The shared materials table (GEN.427).** `physics/materials.py` loads `physics/data/materials.json`: 158 rows, from part E's materials table, part B's SLB 2024 mantle phases and iron rows, part C's ice phases and the numbers read from Seager et al. 2007 (Tables 1 to 3), the Handbook of Mineralogy and PubChem. Every number carries its source and a read, computed, second-hand or not-reached code; a value with no source is empty, and part E's recalled cells are dropped. Each row has a class gate (`any`, or `planetary` for the high-pressure phases and the envelopes), a source and a verification state, and the rows are registered as the one shared `materials` table that other tables' fractions refer to. The 84 asteroid components, the field families, the comet list and the planet class composition texts all resolve to rows, and a test fails when one does not. Nothing reads it yet, so no draw changes.

### Changed
- **The licence is no longer a column.** The gas, element, isotope and decay-series files, the possibility-table row format, the dataset catalogue and the new materials table lose their licence field (Boss, 2026-10-22); the licence and source notes are in `DATA-LICENSE.md` in the repository root.

generation-output: unchanged
