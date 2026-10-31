**Headline:** The bright-star backfill path a run really takes now has a test that runs it directly.

### Changed
- **`backfill_bright_stars_unit` is tested against the serial backfill (TEST.532).** A run backfills in `scatter-backfill` units, but the older tests covered the serial `backfill_bright_stars_around`. A new test runs `expand_backfill_bright_stars` and each unit through `backfill_bright_stars_unit` on one database and the serial form on another, for the same new sector, and shows the sectors, stars (mass, luminosity, type, position, seed) and the unit's per-sector progress reports agree. The serial form stays as the test reference (it is already registered by design in the wiring register). No code or generated output changes.

generation-output: unchanged

TODO items resolved: TEST.532.
