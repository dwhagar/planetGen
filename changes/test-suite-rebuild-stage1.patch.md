### Added
- **`scripts/test.py` picks a slice of the test suite** (TEST.139): `quick` (not database, browser or heavy tests), `area web db` (by part of the program), `browser`, and `full` (the merge gate, in two phases: every core for the non-browser tests, then 2 workers for the browser tests so software WebGL is not starved). It prints the start time and how long the same command took last time on this checkout.
- New test markers `heavy` (the 19 files over 100 s each) and `area_web`, `area_api`, `area_db`, `area_generation`, `area_physics`, `area_names`, `area_wiki`, `area_ops`, set from the test file's name (`src/tests/areas.py`); no test file moved.
- **Failures are kept and re-run alone** (`.test-failures.jsonl`): every failing test's full report with the worker and the machine's load, and when the run ends the failed tests run again alone with the verdict "passes alone" (a load failure) or "FAILS ALONE" (a real one). `quick` takes 2 min 33 s for 11,545 tests.

### Fixed
- **The `db` test marker missed 2,210 of the 3,681 tests that need MySQL**, so `-m "not db"` still ran them (they skipped without a server); tests that reach the server through `_mysql_server_available` are marked too.
