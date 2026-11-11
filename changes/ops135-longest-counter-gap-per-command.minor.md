**Headline:** Every command now remembers the longest time its counter has stood still between two changes, so the heartbeat's "quiet" rule is 5 times what that command has really done and not a guess (OPS.135, heartbeat part 2).

### Added
- **The gap, measured in the horse.** `UnitBeat.longest_gap()` is the longest seconds between two changes of the unit's counter, from the unit's start to its end, less the time spent waiting on the sky lock, MySQL or the gate (waiting is not standing still). It is also in the beat record as `longest_gap`.
- **`queue/commandgaps.py`.** A finished unit notes its gap; the run keeps the largest per command (`runrecord.gaps_key`, a hash of command to seconds, expiring with the run). A run that ends `done` folds its gaps into the server's (`planetgen:command-gaps`, beside `stagetimes`): a longer gap replaces the stored one, a shorter one pulls it down by half, so one slow night does not loosen the rule for ever. A failed or stopped run adds nothing.
- **For the judge.** `commandgaps.history_of(conn, run_id)` is the `history_of` of the `UnitWatch` that `wrapper.run` makes: the larger of the server's gap and this run's so far, which `beatjudge.quiet_limit_seconds` turns into `max(5 min, 5 x the longest gap)`.

### Changed
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
