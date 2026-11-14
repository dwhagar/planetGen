**Headline:** Queued runs work again: the database step parameters are turned into a connection config by `MySQLConfig.from_args`, which the cleanup rounds called but never defined, so every queued galaxy run stopped at its first step.

### Fixed
- **`MySQLConfig.from_args` exists.** `queue/wrapper.py`, `queue/units.py` and `db/check.py` (layered-fixes rounds 2-4, PR 6454) call `store.MySQLConfig.from_args(named)` for a step's `database` parameters, but the class had no such method (`AttributeError: type object 'MySQLConfig' has no attribute 'from_args'`), so a queued run failed in the wrapper's database clock. The method builds the config exactly as the three call sites did before the cleanup (host, port, user and database named; the password is the process's settings'). Found by the gate on `test_work_units::test_a_galaxy_run_stores_each_of_its_stages_once_and_follows_its_sector_units`.

generation-output: unchanged
generation-affects: none
