**Headline:** The habitability re-score no longer dies with a MySQL deadlock when four workers run its units at once.

### Fixed
- `rescore-habitability` ran each unit as one transaction at the default isolation level, so two workers re-scoring neighbouring id ranges held gap locks on the shared `system_interests` and planet index ranges and each waited on the other (MySQL error 1213), which stopped the stage. The unit now runs at READ COMMITTED, which takes no gap locks, and starts over on a deadlock or a lock wait timeout (`store.RETRYABLE_ERRORS`) on a connection of its own, up to `store.SECTOR_SAVE_ATTEMPTS` times (8) with a pause that grows with each try, then fails as before. The numbers are the same whenever a unit is run, so a retry changes nothing stored.

### Added
- `tests/test_perf228_rescore_deadlock.py`: four threads re-scoring neighbouring ranges of 16 saved systems for 12 rounds each finish with no error (it failed on the old unit with error 1213), a unit that meets one deadlock gives the same answer, a unit that keeps meeting one gives up after 8 attempts, and any other database error is not retried.

generation-output: unchanged
TODO items resolved: PERF.228
