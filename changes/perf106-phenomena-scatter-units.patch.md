**Headline:** The phenomenon scatter's shells run as commands on the work queue, like the star scatters.

### Changed
- **`scatter-phenomena` is a command (PERF.106).** Its shells (neutron stars, black holes, planetary nebulae, supernova remnants, hypervelocity stars, and the clusters' and satellites' own) are cut into units of consecutive shells, queued through the wrapper as a run of their own and followed with the same bar and per-shell lines as before; a unit runs alone from its JSON arguments (`planetgen command-run`). The same rows are written at any worker count, and on the test galaxy they are the rows the previous version wrote.
- **The units' speeds are recorded by the unit that completes the stage**, the drawing and the writing of the rows under their own kinds as before. The stage's row in the stage history is still stored by the process that started the run, because the stage also clears the earlier scatter, draws the special phenomena, builds the nebula field and the indexes around its units.
- **The old per-run task queue of the phenomenon scatter is gone** with no replacement.
- **`docs/design/redis-and-the-work-queue.md` describes it.**

TODO items resolved: none (slice three of PERF.106). Added: none.
