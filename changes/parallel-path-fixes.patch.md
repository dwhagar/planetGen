### Fixed

- A generation run no longer hangs forever on Python 3.12 when a worker process dies. The pool stops the other workers with SIGTERM, but a worker caught it in the middle of a task, ended only that task and then waited for another one, so Python 3.12's pool cleanup waited on it for good. A worker sent SIGTERM now exits as soon as its task has been reported (PERF.22).
- The bright-star progress bar no longer ends at 101%: a progress report that arrives after its layer has finished is ignored, and the terminal, `progress.json` and the Generate page never show a share above 100% (PERF.23).
- A bright-star placement test now works on Python 3.9 and 3.10 (TEST.76).
