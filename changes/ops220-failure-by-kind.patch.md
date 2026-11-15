**Headline:** a unit that fails because the database or Redis hiccuped is now run again by the engine after a growing wait, and the run says what kind of failure ended it.

### Added

- OPS.220: `queue/failures.py` names every failure one of four kinds (passing, killed, idle, lasting). A passing failure (MySQL deadlock 1213, lock wait 1205, connection lost 2006 or 2013, a closed connection, a Redis connection that failed) runs the unit again after a wait that doubles from 2 s to 5 min, spread at random, inside a per-run budget of 40 retries and 30 minutes of waiting. Every attempt is a row in the ledger, and the budget is one Redis hash per run.

### Changed

- OPS.220: the run's error text names the kind of failure. A killed or idle unit keeps the existing alone-with-a-bigger-claim and stack-dump paths. `rq.Retry(max=1)` and the in-transaction deadlock retries stay until OPS.219 (checkpoints), OPS.221 (held items) and the supervisor hand-out remove their last users.
