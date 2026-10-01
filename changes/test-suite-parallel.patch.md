### Changed
- **The test suite runs in parallel.** `pytest -n auto` (pytest-xdist, now
  in the `test` extra) runs one worker per core: the full suite went from
  17 min 21 s to 5 min 57 s on 4 cores, and CI and the weekly deep fuzz
  run use it. Each test process now gets its own throwaway control
  database instead of falling back to `planetgen_control`, so workers
  never share one and the tests never touch a real control schema.
  Plain `pytest` still runs serially.
