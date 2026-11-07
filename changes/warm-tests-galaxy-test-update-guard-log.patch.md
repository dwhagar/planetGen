**Headline:** The tests of the warm-up pass again after it learned to skip an empty galaxy, and the update guard writes its one warning through the program's own logger.

### Changed

- `src/tests/test_warm_phases.py`, `test_warm_log.py`, `test_cache_warmup.py`, `test_maintenance.py` (found while checking OPS.114): the warm-up now asks the database whether it holds a galaxy before it warms anything and, when it holds none, says so and exits 0. Ten tests ran the warm-up against a stand-in database that is only a name (`x`, or a Redis key prefix), or against arguments with no database at all, and so failed with "Unknown database" or a missing `database`. They now run as if the database held a galaxy; `test_warm_empty_galaxy.py` still holds the empty case. The warm-up itself is unchanged.
- `src/planetgen/admin/updateguard.py`: the update guard made its own standard-library logger (`logging.getLogger`), which `test_no_own_log_formatting.py` refuses (every line goes through `planetgen.util.log`). Its one warning, about a run record it cannot read, now goes through `log.warning`. Same text, same level; it appears in the program's log format.

generation-output: unchanged
