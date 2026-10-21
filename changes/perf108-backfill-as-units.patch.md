**Headline:** The bright-star backfill of `planetgen galaxy` runs on the work queue as a command.

### Changed
- **`scatter-backfill` is a command (PERF.108).** The backfill after a `galaxy` run finds the sectors to draw once, cuts them into chunks of at least 10,000 sectors (at most 2,000 chunks) and queues them through the wrapper, so the chunks can run on several workers; every sector is drawn from its own seed, so the stars are the same at one and at four workers and whichever way the sectors are cut. The bar counts the sectors visited, the per-shell lines read as before, and the stage's row in the stage history is stored by the unit that completes it. `docs/design/redis-and-the-work-queue.md` describes it.

TODO items resolved: none (the backfill of PERF.108). Added: none.
