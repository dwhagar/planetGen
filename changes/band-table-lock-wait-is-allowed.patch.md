**Headline:** The sampling-table lock wait is named in the no-sleep test, so main passes it again.

### Fixed

- PERF.281: `test_no_sleep_in_commands` (OPS.197) lists the table lock's poll in `band_table_cache._building` with its reason: the wait is declared with `unitbeat.waiting` and beats the unit on every poll, so it is not a silent stretch. A blocking `flock` would be silent to the unit's heartbeat, so the poll stays. No behavior changes.

generation-output: unchanged
generation-affects: none
