**Headline:** The process that follows a run now judges every running unit's heartbeat at each poll, says once in the run log and `error.log` when a unit goes quiet, dead or loops, and shows each unit's state on a new API route (OPS.138, heartbeat part 3).

### Added
- **`queue/unitwatch.py` `UnitWatch`.** At every poll `follow_stage` has the watch read the beats of the stage's running units (`unitbeat.read`), judge each with `beatjudge.judge` on the Redis server clock (with the beat the unit had at the check before, so a count that goes backward inside a phase is seen), and log **each change of state once**: a run log line (`unit-state`, with the stage and unit) and, for dead, quiet and suspect-looping units, a line in `error.log` too. A unit that works from the start is not announced; one that works again after a flag gets a line saying so; a second attempt of a unit that finished is a new change. `UnitWatch.agreed()` lists the units dead or quiet in two checks in a row (`beatjudge.agreed`), which is what the stuck action (OPS.140) takes; a unit whose count moves is never in it.
- **`GET /api/admin/work/<id>/units`** (admin): the state of each running unit as the follower last judged it (state, severity, reason, whether two checks agree, counter, rate, projected finish, seconds still, phase, what it waits on, done and estimate, restarts, resident MB, beat age, running seconds). The verdicts are kept in `runrecord.unitstate_key` (a hash to a stage, cleared when the stage ends and deleted with the run). The API version stays: a new route is additive.

### Changed
- `units.follow_stage` takes `watch`. Only the process that started the run (`wrapper.run`) makes one; `wrapper.follow`, which only watches, does not, so two processes never log the same change twice.
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
