**Headline:** `planetgen sector` and `planetgen system` run through the work queue like every other generation command, so they show on the Queue page and take their turn.

### Changed
- **`sector` and `system` are runs of one command (PERF.119, first slice).** Each builds its step list from what was asked and hands it to the wrapper: the sectors of a `sector` run are cut into units that workers generate and save (each sector draws from its own seed, so a run gives the same sectors at any worker count), and a `system` run is one unit. A `system --output` run saves nothing, so it is a maintenance run and does not hold up generation. The old in-process pool for these two commands is gone; both now need Redis, as the other generation commands do.

TODO items resolved: none (PERF.119 slice 1 of 3). Added: none.
