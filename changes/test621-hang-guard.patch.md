**Headline:** A test that hangs past its time limit is now stopped, named and reported, and the test run goes on instead of waiting for hours.

### Fixed
- The full non-browser run hung twice with every worker idle (1.5 h and 3 h). `pytest-timeout`'s signal fires once and raises in the test; a test that holds a process pool then waits for its children as the `with` block closes (`shutdown(wait=True)`), and a child that is still busy never ends, so the alarm was spent and nothing stopped the wait. Reproduced with a two-child pool that sleeps.
- New `tests/hang_guard.py` (loaded from `pytest.ini`): a second clock on a thread of its own for every test, over setup, call and teardown. Past the test's limit plus `hang_guard_grace` (120 s) it prints `HUNG TEST <node id> (worker ...)` and kills the worker's child processes, so the wait ends and the timeout is reported as the test's failure; 120 s later, if the test is still running, it ends the worker (`os._exit(70)`), which xdist reports as a crashed worker on that test and replaces.
- `scripts/test.py` stops a phase that prints nothing for two hours (`PLANETGEN_TEST_IDLE_LIMIT_S`, 0 turns it off), says `HUNG RUN`, and returns exit code 124.

### Added
- `tests/test_hang_guard.py` (a hung pool is killed and the next test runs; a test stuck on nothing a child can end loses its worker under xdist and the run goes on; no limit, no watch; the process-tree helpers) and three tests of the idle limit in `tests/test_test_runner.py`.

generation-output: unchanged
TODO items resolved: TEST.621
