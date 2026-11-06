**Headline:** A cache warm-up no longer takes every worker: a database check, a visitor's map nodes and the short jobs get workers while it runs.

### Changed
- The maintenance workers now have a queue for each kind of task (the warm-up's units, a database check's units, the map-node prefetch and rebuild jobs) and take one job from each in turn, so a check or a visitor's nodes no longer wait behind thousands of warm-up units. Measured on a 3-worker box with 300 warm-up units queued: a check and a node job waited 11.5 s before and 0.1 s after.
- The host's worker cap now keeps places for the short and long job lanes, so a warm-up or a generation run that has used up the host's workers never leaves a search or a sector list with no worker: with the cap used up, the short job waited for ever before and 1.7 s after.
- The workers for warm-up units, node jobs and slow page answers build the site's app once when they start instead of for every job: 30 warm-up units that each build the app took 69.2 s before and 9.8 s after, and a second short job 6.9 s before and 0.3 s after (the first job after a worker starts takes as long as before).
- The whole-map build runs at the lowered priority queue workers use and with one process per core less the places kept for short and long jobs (`tile_cache.build_workers` 0), so it no longer takes the cores a visitor's request needs.

generation-output: unchanged
