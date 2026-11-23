**Headline:** The test that checks the pages show what the command line wrote now plans the same galaxy every run, so how long its command-line step takes no longer changes from run to run.

### Fixed
- TEST.825: `test_web_fresh_after_cli.py` planned its toy galaxy with no seed, so every run drew a new galaxy and the one-sector command-line step took anywhere from 6 s to 223 s on one quiet box (same code, six runs each on two heads: 5.9 to 223.1 s and 6.2 to 131.0 s). `_plan_galaxy` now passes `GALAXY_SEED` (`bytes(range(16))`, the seed `test_cache_warmup.py` plans with, the only other test that plans such a galaxy and it already seeds), so each run builds the same galaxy: the same step, six runs, took 28.5 to 30.7 s. Every assertion is kept; no source file changes. The file alone on a quiet box: before 23m 10s and 31m 18s (red at line 236, the window TEST.824 fixed), after (with the seed and TEST.824 on main) 16m 29s with 4 passed.

generation-output: unchanged
