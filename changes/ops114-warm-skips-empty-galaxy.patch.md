**Headline:** The cache warm-up and the nightly whole-map build now skip an empty galaxy cleanly instead of trying to work.

### Changed
- **OPS.114: a warm-up run right after a New galaxy wipe tried every task.** `planetgen warm` (from update.sh, after a change, or by hand) and the whole-map build (`planetgen.cli.tiles`, the nightly timer) now check first whether the database holds any galaxy (no sector and no galaxy plan row, or no tables yet). If not, each prints `galaxy is empty, nothing to warm up`, does nothing and exits 0. update.sh reports that step as skipped because the galaxy is empty, not as failed, so the Update button and the maintenance timers do not report a failure. Once a galaxy exists the warm-up runs as before.
