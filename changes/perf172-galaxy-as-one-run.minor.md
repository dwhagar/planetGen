**Headline:** `planetgen galaxy` is one run in the queue, so one estimate covers the whole galaxy.

### Changed
- **The galaxy is one run.** `planetgen galaxy` plans its sectors, then hands the wrapper one step list: the sectors, the link, the mass and bright-star scatters, the phenomena and their indexes and nebula field, the block light totals, the backfill, the population pass, the territories, the settle and the map's block totals. The estimate shows the whole list, the Queue page shows one run, and the same seed makes the same galaxy at one worker and at four. The process that starts the run only follows it; the work between a stage's units moved into the stage's command.
- **A stage that cannot start, or cannot end, fails the run and says which.** The empty-the-table step before a scatter and the store-the-record step after it used to run in the starting process; they now run in a worker, and a failure shows as "stage X could not start" or "the end of stage X failed" instead of a unit that looked done.
- **A link that fails now fails the run** (it only warned before), and a failed galaxy prints the `planetgen link --since` line to run it again.
- **A stage with nothing to do is stored as skipped** (no sector saved since the run started, nothing to backfill), so it adds nothing to the estimates.
- **`planetgen scatter` ends with the map's block light totals** as its own stage, and the phenomena scatter shows its indexes and nebula field as stages of their own.

### Added
- **Commands `finish-phenomena`, `refresh-block-lights` and `refresh-run-blocks`**, and the stages `phenomena-index`, `block-lights` and `run-blocks` (`planetgen commands` lists them).
