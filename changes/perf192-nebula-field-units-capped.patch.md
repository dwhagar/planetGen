**Headline:** A galaxy with more than 4,000 nebula field slabs no longer fails the nebula field stage; its units cover several slabs each.

### Changed
- **The nebula field stage cuts like the other fixed-size stages.** It queued one unit per 50 pc slab, so a galaxy wide enough for more than 4,000 slabs failed with "a stage queues at most 4000 units". Beyond that size a unit now draws a run of consecutive slabs (`units.coarsen`), the same cut at any worker count and the same clouds; below it nothing changes. The log still has a line per slab.
