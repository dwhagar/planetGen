**Headline:** One file, error.log in the log folder, now records every error and crash from everything, while each error still goes to its own log too.

### Changed

- `error.log` (beside `warmup.log` and the others, `logpaths.error_log_path`) holds every line of severity `error` and `critical` that any logger writes, the program's or another library's (OPS.116). A line is in the display form with the fields `log=` (the log the error also went to: `warmup.log`, `tiles.log`, `orbits.log`, `update.log`, `apache` for the web app, else `console`) and `pid=`, and a traceback follows its line. It is made on the first error, mode 0660, by the one logging module (`util/log.py`); the setup scripts make it for the web server's user, the folder's rotation rule (`*.log`) covers it, and the macOS rotation file names it.
- Crashes: an uncaught exception in the main thread or another is written as `critical` with its traceback whether or not the debug log is on, and a fatal signal in a library (SIGSEGV and the like) leaves the interpreter's dump of every thread in the file. A Ctrl+C is not an error.
- Killed or crashed processes are written by their parent: a detached worker (exit status or signal), a work horse the out-of-memory killer or `kill -9` ended (`units.horse_killed`, for any job, not only a run's unit), a job a dead worker left running (`reap_dead_jobs`), a web job's step that exited non-zero and a one-off command that failed, timed out or did not start, each with the last error lines the program printed. A stop on purpose (SIGTERM, SIGINT, SIGHUP, an update's refusal) is not an error.
- The web app writes an unhandled exception with its route and traceback, and any other 5xx answer with its route, once; the deliberate 503 of a closed site or a cache being rebuilt is not an error.
- `update.sh` and `install.sh` write every error they print to the file too, name the command that made them stop with its exit status and the step, and a warm-up stage that does not finish is written there though the screen only warns.
- `docs/design/operations/logging-standard.md` has an "The error log" section, `docs/reference/logging.md` lists the file, and `log_dir`'s help names it. `src/tests/test_error_log.py` holds one test per source.

generation-output: unchanged
