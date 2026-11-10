**Headline:** A running unit now tells the queue it is alive from a thread of its own, on the Redis server clock, so a long quiet stretch of work no longer looks like a dead unit (OPS.136, heartbeat part 2).

### Added
- **The beat record.** `queue/unitbeat.py` `UnitBeat` writes one JSON record per running unit into the stage's progress hash (field = unit): `done`, `estimate` (0 = no counter), `beat` (the Redis server time of the write, `conn.time()`, so a horse with a wrong clock cannot look dead or alive), `started`, `moved` (when `done` last changed; absent until it first does), `phase`, `waiting_on` (`sky lock`, `mysql`, or empty), `restarts` (times the current phase began again) and `resident_mb`. `unitbeat.read(conn, run_id, stage)` returns `{unit: beatjudge.Beat}` for the judge (OPS.137).
- **Tags.** `unitbeat.waiting(tag)` and `unitbeat.phase(text)` reach the running unit from code that is not handed the reporter. The sky draw lock (`DrawLock`) and the contention grant (`db/contention.py`) tag their waits, so a unit waiting for another unit's turn is `waiting`, not `quiet`.

### Changed
- **The write moved off the work loop.** The unit's reporter only sets numbers in its own process; a thread of the horse writes the record every 10 s (`BEAT_SECONDS`), and sooner, at most every 2 s, while the counter moves. A frozen work loop (a long GET_LOCK, a big query, a C call) keeps a fresh beat; a hung process stops it. The follower reads `done` and `estimate` from the record as it read the old two-number field. Redis writes for a fast reporter go from up to four a second (`REPORT_SECONDS` 0.25) to at most one every 2 s.
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
