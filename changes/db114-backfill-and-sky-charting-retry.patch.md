**Headline:** The bright-star backfills, the special phenomena and the sky charting no longer stop a run on a MySQL deadlock.

### Fixed

- **DB.114: every writer of the bright-star and phenomenon scatter tables now runs the Object ID step inside a transaction that is retried.** After DB.113 fixed the scatter units, the sector backfills (by luminosity and by mass), the special phenomena write and each batch of the sky charting still ran the same ID step with other writes in the same transaction and no retry, so a deadlock (MySQL error 1213) or lock wait timeout (1205) with other workers stopped the run. Each of them is now one whole transaction at READ COMMITTED (the sector locks, the draw, the stars with their IDs, the levels and the ledger) that `db.contention.run_with_retries` rolls back and writes again, with the same waits, logging and exclusive grant after the eighth failed try as DB.113. The stars written are the ones an undisturbed run writes, and a progress tick is not given twice.

generation-output: unchanged
