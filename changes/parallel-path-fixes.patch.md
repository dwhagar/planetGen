### Fixed

- Generation with more than one worker no longer stalls for 50 seconds at a time. A worker waiting for the neighbour lock kept the rows it had written, and when the lock's holder needed one of them (a name collision renames an existing system) neither moved until MySQL's lock wait timeout. The holder now gives up after 3 seconds and starts its sector save over. A seven-ring test run at 2 workers went from 266 s to 18 s (PERF.21).
- A generation run no longer hangs forever on Python 3.12 when a worker process dies. The pool stops the other workers with SIGTERM, but a worker caught it in the middle of a task, ended only that task and then waited for another one, so Python 3.12's pool cleanup waited on it for good. A worker sent SIGTERM now exits as soon as its task has been reported (PERF.22).
- Cancelling or interrupting a parallel run (SIGTERM, Ctrl+C) can no longer be lost. When the signal landed inside a database call it became a database error the work queue logged and carried on past, and the run finished as if nothing had happened. The queue now remembers the signal, stops at the next step and ends the job as cancelled (TEST.73).
- The bright-star progress bar no longer ends at 101%: a progress report that arrives after its layer has finished is ignored, and the terminal, `progress.json` and the Generate page never show a share above 100% (PERF.23).
- Re-running `generate.py plan --bright-stars-down-to` after a run of it stopped part way no longer puts the band in twice: the re-run first removes what the stopped run left below the star-fill level, then draws the whole band (GEN.32).
- A bright-star placement test now works on Python 3.9 and 3.10 (TEST.76).

### Added

- CI runs the generation tests again with 2 and 4 worker processes (Python 3.12 and 3.13), and tests can patch the workers as well as their own process (`tests/worker_patches.py`, described in `docs/testing.md`) (TEST.74, PERF.21).
