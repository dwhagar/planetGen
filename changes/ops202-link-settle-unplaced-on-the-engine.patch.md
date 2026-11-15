**Headline:** The neighbour link, the sector-path settle and the `sector` run's sectors now run as pieces the job engine counts, so their units show real counts and honour a cancel or a pause between two batches.

### Changed

- OPS.202, first part: `link-neighbors` gives the engine one piece per batch of `store.LINK_BATCH_SECTORS` sectors (each batch still one transaction under the neighbour lock, retried on a deadlock as before), `settle-sector-paths` one piece per batch of `sector_paths.SETTLE_BATCH` sectors that sit together (`sector_paths.settle_batch`; each sector committed through `ctx.commit()`), and `generate-unplaced-sector` one piece per sector of a `sector` run. The old unit functions `link_neighbors_unit`, `save_sector_paths_unit` and `unplaced_sectors_unit` are gone.
- `settle_sectors` (for the API's edits and the tests) runs the same `settle_batch`; the process pool a crowded sector's bodies are integrated in is opened by the first batch and shut down when the unit ends (new `ctx.on_close`, which the engine runs last in first out however the unit ends).
- A stage's results are unchanged: `linked`, `saved` and `visited`, and the sectors a `sector` run reports.

generation-output: unchanged
generation-affects: none
