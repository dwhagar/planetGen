**Headline:** A unit whose failed attempt is waiting to be retried is no longer taken for dead and run twice, a second save of a cell waits and finds it charted, and what holds or restarts a unit is said in the job's log.

### Fixed
- A silent beat left by an attempt that raised is not a dead unit while RQ has its retry waiting in the queue: the follower no longer stops it, queues it a second time and holds the whole queue for it to run alone, which put two copies of one sector unit on the same cell at once (the `uq_star_systems_uid` duplicate in Boss's neighborhood run).
- Saving a sector takes the cell's named lock and looks at the chart ledger inside the transaction, so a second copy that finishes after the first returns the first one's sector instead of stopping on the object IDs the first gave its systems (DB.131's check ran before the sector was generated, and could not see a copy still running).

### Changed
- The follower's lines (a unit stuck, stopped, queued again, held to run alone, released) and a unit held after failing twice are written to the job's log as well as the run's event stream.
- The no-unit-finished line names the stage ("Stage 'Generate the neighborhood': ...") instead of the run's own stage number, which read as Stage 3 inside a job on its Stage 5.
