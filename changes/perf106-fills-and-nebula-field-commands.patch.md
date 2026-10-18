**Headline:** The map luminosity fill, the habitability re-score and the nebula field run on the work queue as commands.

### Changed
- **`fill-map-luminosity` and `rescore-habitability` are commands (PERF.106).** `planetgen plan --fill-map-luminosity` and `--fill-habitability` queue a unit per range of ids (200,000 scattered rows, 500 systems) through the wrapper, so the ranges can run on several workers; each unit writes only its own rows, so the result is the same at any worker count. The bar counts ids, and the stage's row in the stage history is stored by the unit that completes it.
- **`clean-nebula-field` and `draw-nebula-field` are commands (PERF.106).** `--redo-nebula-field` checks the stored field clouds in one unit (the dry run is always done and logged first; with `apply` the same unit deletes the rejected ones and refreshes containment), and the nebula field table is drawn by a unit per slab of cells, here and in the phenomena scatter, replacing the older work queue the draw used. The table holds the same rows at one and four workers, as it did at one before.
- **`docs/design/redis-and-the-work-queue.md` describes them.**

TODO items resolved: none (slice five of PERF.106). Added: none.
