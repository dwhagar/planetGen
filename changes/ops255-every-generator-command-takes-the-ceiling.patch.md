**Headline:** every way a generation job starts now takes the 30-day job limit.

### Fixed

- OPS.255: OPS.255: every command of the generator in a Generate job (finish the last run, sector, system, warm, check-db, backfill and the rest, not only galaxy, plan and run-job) takes the 30-day RQ ceiling, so no way of starting a run on the queue is ended by a clock of a smaller galaxy's recorded time.
