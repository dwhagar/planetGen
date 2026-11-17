**Headline:** the rich test galaxy can be built once for every machine that tests the same code, and loaded from a dump by the rest.

### Changed
- `tests/rich_galaxy_cache.py`: with `PLANETGEN_RICH_GALAXY_CACHE` set to a folder every machine can reach, the first process that needs the `rich_galaxy` of a code state builds it as before and leaves a gzipped `mysqldump` of its database there; the others load the dump in seconds instead of building (about half an hour). `test_db_column_round_trip`, `test_db_check_constraints` and `test_cellkey_motion` each build it, on a merge train on three machines from the same commit. The galaxy is the same galaxy with the same systems (nothing made smaller), a minus zero in a double is kept (`mysqldump` would write it as plus zero), and without the variable nothing changes.
- The dump is named by a key over every file of the package, the support modules of the build, the system file and `config.json`, the settings in `PLANETGEN_*` variables, the dataset files (names and sizes), and the versions of the database server, Python, numpy and numba; a different commit is a different key. A process that finds the build taken waits for it (up to 3 hours, then builds its own without caching); a lock older than 3 hours is taken over.

generation-output: unchanged
