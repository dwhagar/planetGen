**Headline:** `planetgen reset-jobs` starts again: it no longer looks for the worker-start counters that OPS.235 slice 1b removed.

### Fixed
- `cli/reset_jobs.py` named `workers.STARTING_PREFIX`, which slice 1b deleted, so importing the command raised `AttributeError` and `tests/test_reset_jobs.py` failed at collection. The key patterns are now the RQ keys and the runs' keys; its docstring no longer mentions the start counters or the admin's worker count.
- `tests/test_reset_jobs.py` and `tests/test_maintenance_lane.py` no longer set or clean the counters and the worker-count key.
