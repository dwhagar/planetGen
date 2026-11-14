**Headline:** The supervisor now judges the units of every live run, so a run whose starting process is gone is still watched.

### Changed
- **OPS.208, output unchanged.** The beat judge of the running units (`unitwatch.UnitWatch`) and the stuck action (`unitstuck.Stuck`: dump, stop, queue again once, then fail the run) run in the supervisor (`supervisor.RunWatches`), one turn a second while a run is live, for the stage in flight of each live run. `wrapper.run` and `units.follow_stage` only follow now: they lose their `watch` argument. A run whose starter was killed no longer goes unjudged. Judgments, log lines and the Queue page's unit states are the same.
generation-output: unchanged
generation-affects: none
