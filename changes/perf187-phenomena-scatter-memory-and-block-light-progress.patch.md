**Headline:** The phenomena scatter no longer holds gigabytes in one process, and the block light totals at its end show their progress instead of going quiet.

### Changed
- **The block light totals are written as they are complete (PERF.187).** `db.block_lights.rebuild` kept every block of a unit of 729 sector shells in memory until the unit ended; on a 12 pc galaxy that was the scatter step's 6.4 GB (7.9 GB measured while the step ran). It now reads one shell at a time and writes a block as soon as the shells have passed it, so what it holds is the blocks of the current block shell of each level.
- **The nebula block totals write one level at a time (PERF.187).** `galaxy.nebula_blocks.iter_rows` gives the levels one by one and lets go of a level's cover once its rows are made; `db.nebula_blocks.rebuild` writes each level as it comes instead of holding the rows of all of them. `all_rows` still returns the dict for the callers that want it.
- **The block light totals report their progress (PERF.189, Boss).** After the hypervelocity stars' epoch is stamped, a scatter used to go quiet while the map's block light totals were made. They are now a step of their own ("Making the block light totals", bar over the units), with a log line before the scan, after it, before the levels above the units and about every 30 seconds (`tuning.PROGRESS_LOG_SECONDS`).
- **The smoke budget for the scatter step is lowered (PERF.187)** to the new measured peak.
- **The logging standard has the rule for work after a progress block (DOC.80):** it needs a step of its own, with a short checklist.
