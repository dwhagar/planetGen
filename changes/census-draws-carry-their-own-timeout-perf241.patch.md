**Headline:** The three live census tests state their own time limit (3,600 s), because on 4 vCPUs each takes 886 to 1,036 s against the default 900 s.

### Changed

- `src/tests/test_prevalence_census.py`, `test_prevalence_life_census.py`, `test_prevalence_systems_census.py` (PERF.241): `@pytest.mark.timeout(3600)` on the live census test of each. They draw every cell with the number of draws their reference was measured with (300 per cell, 100 systems per star type), already on every core through a process pool, so a smaller draw would change what the ranges mean. Measured on the Train Manager's 4 vCPU machine after train 12: class census 1,020 s (23 min 39 s for the file), life census 886 s (14 min 46 s for the file), systems census 1,036 s (17 min 16 s for the file). The limit is the one `tests/hang_guard.py` holds a test with its own mark to; a hang still ends the test at 3,600 s.

generation-output: unchanged
