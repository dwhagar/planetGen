**Headline:** `scripts/test.py` checks MySQL, Redis and the test setup before a run, restarts a server that dies during it and carries on without redoing the tests that passed; `full` is no longer called the merge gate.

### Added
- **A check before every run (TEST.376).** `scripts/test.py` loads the `PLANETGEN_TEST_*` variables from `/etc/profile.d/planetgen-test.sh` when they are missing, starts MySQL and Redis when they do not answer, purges the MySQL binary logs when under 3 GB of disk is free, and checks `argon2-cffi` and `playwright` 1.56.0. What it cannot repair stops the run with `PREFLIGHT FAILED: ...` and exit 2. `scripts/test.py env` does only the check; `--no-preflight` skips it. Before this a database test with MySQL down was skipped, not failed, and a run without `PLANETGEN_TEST_MYSQL_USER` failed with `DatabaseLoginError` as the wrong account.
- **A watchdog during every run.** MySQL and Redis are asked every 15 seconds. A server that is gone is restarted, the phase's pytest is stopped and started again without the tests that already passed (the progress log, as `--resume` does), up to 5 times a phase; the closing line says how many restarts there were. A server that cannot be restarted ends the run with `FAIL` and exit 3, and `--resume` continues it.

### Changed
- **`full` is the end-of-version run, not the merge gate (TEST.372).** The help of `scripts/test.py` and docs/reference/testing.md now say an ordinary merge needs a green `changed` run (finished with `--resume` after an interruption), and that the whole suite is for TEST.261, SEC.33 and GEN.301.

TODO items resolved: TEST.372, TEST.376. Added: none.
