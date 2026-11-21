**Headline:** The rich test galaxy is built once per machine and loaded by the other test files, which takes the three files that use it from about 50 minutes each to about 50 minutes for the first and under two minutes for the others.

### Changed

- The rich galaxy's dump cache (TEST.788) was off unless `PLANETGEN_RICH_GALAXY_CACHE` named a folder, so each of `test_cellkey_motion`, `test_db_check_constraints` and `test_db_column_round_trip` built its own copy of the galaxy, which GEN.1745 slice 4b made slow (62 minutes for the first on the Train 51 gate). It is on by default now, in `planetgen-rich-galaxy` in the machine's temporary folder; a folder every machine reaches shares it between machines, and the variable set empty turns it off. A dump nobody loaded for three days is removed when a new one is kept.
- The cache could never have hit under pytest: the key held the folders and databases each pytest process makes for itself (`PLANETGEN_LOG_DIR`, `PLANETGEN_MAINTENANCE_FILE`, `PLANETGEN_CONTROL_DATABASE` and the like), so no two processes had the same key. Variables that only say where something is (names ending `_DIR`, `_FILE`, `_DATABASE`) are out of the key; every setting the generator reads, every file of the package, the data and the support files stay in it, so a change to generation code is a new key and a fresh build.
- A process that loads a dump now makes the galaxy's cluster field, as the warm-up does; the dump does not carry it, and without it the pages that need it answered 503 (TEST.807).

generation-output: unchanged
