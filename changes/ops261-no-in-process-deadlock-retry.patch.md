**Headline:** A deadlock or lost connection in a write unit now reaches the job engine, which retries it from its cursor; the in-process retry loop is gone.

### Changed
- `db.contention` keeps `begin_read_committed` and `run_once` (one try, rollback, lost-connection explained, error raised). The sleeping retry loop, its gate and its waits are removed (OPS.261); the engine's failure classes and per-run budget (OPS.219, OPS.220) retry from the cursor.

generation-output: unchanged
