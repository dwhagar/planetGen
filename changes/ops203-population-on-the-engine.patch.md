**Headline:** The population and territories passes now run as pieces the job engine counts, so the Queue page shows real numbers for the scan, the dating and the territories instead of one quiet unit.

### Changed

- OPS.203: `run-population` and `run-territories` hand the engine their pieces (`population.model.population_pieces`, `territories_pieces`). Population: a forget for a rescan, one scan piece for each batch of planet ids, the drop, the dating of the civilizations, the map's technology rows, and the counts. Territories: the recompute, which counts the polities that have claimed and the owner rows written as its second-level bars (OPS.192).
- The engine commits after every piece through `ctx.commit()`; the model's pass code has no commit of its own. A scan piece the watermark has already passed (another pass took the lock first) leaves its planets alone. The old whole-unit drivers `population_unit` and `territories_unit` and the `on_stage` announcements (UX.83's three steps) are gone: the engine's counts are the bar. A bring-to-epoch redo and the tests run the same pieces through `model.run_population` and `model.run_territories`, one commit a piece.
- The pass lock is taken around each piece (the unit's beat says it waits on `mysql` while it waits) and let go after the piece's commit, so a pooled connection never carries it away and two passes take turns a piece at a time.

generation-output: unchanged
generation-affects: none
