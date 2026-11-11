**Headline:** A unit that is dead or stuck is now found by its heartbeat, its stack is dumped, it is stopped and queued again once, and a second time it fails the run with both stacks; a unit whose count moves is never stopped, however long it takes (OPS.140, heartbeat part 7).

### Added
- **`queue/unitstuck.py` `Stuck`.** For a unit that `UnitWatch.agreed()` lists (dead or quiet in two checks in a row): the follower sends its horse `SIGUSR1`, which has `faulthandler` write every thread's stack to the unit's stack file (`unitbeat.register_dump`, set up by `_do_unit`; the beat record carries `pid`, `host` and `stack_file`); stops the job with RQ's stop command; and, once the unit has left its horse, queues the same job again. RQ does not retry a job that was stopped on purpose, so this is done here: the job is taken out of the started or failed registry, the stale beat is cleared, and `redisqueue.requeue_job` queues it. A unit stuck again on its second attempt fails the run, with both stacks in the message. A stop that has not taken effect in `STOP_WAIT_SECONDS` (120 s) fails the run. What was done to each unit is kept in `runrecord.stuck_key`, so a restarted follower does not stop a unit twice.
- **The killed flag** (OPS.144, Lane 13's `runrecord.mark_killed`): the stuck action sets it on every stop, with the reason, and ends a run it fails with `finish_run(killed=True)`, so the Queue page offers Resume run.

### Changed
- **The unit time limit is only a backstop.** `units.unit_timeout` is now always the 6 hour ceiling (`redisqueue.MAX_TIMEOUT_SECONDS`), not 10 times the stage's expected time, and `_unit_limit` (OPS.134) is gone with it. A moving unit is never ended for taking longer than estimated (Boss, 5:16 PM EST on 11-09).
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
