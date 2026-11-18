**Headline:** the job engine's stage 4 kill tests run with real processes, and a unit no longer fails because its template had just been killed.

### Fixed
- OPS.236 (bug): a fork asked of a template that had just died (killed between two turns of the launcher) raised `BrokenPipeError`, which the launcher took for a failure of the unit and marked its row `failed`. `JobHost.fork` now says the template ended (`NotReady`), the host is dead, the launcher puts the row back to `ready` and the next turn starts a new template (`TemplateSet.host`).

### Added
- OPS.236: `tests/test_engine_kill_processes.py` kills real processes under a real template and a real control database: the template while its job runs (the job finishes, and the next job gets a new template), the template and its job together (handed out again once, finished once), a job the kernel kills (retried alone on a claim grown by one and a half), a job killed before it took its attempt, the supervisor killed with a job running (the job stops at its piece, the next supervisor hands it out again, every piece finished once) and a supervisor killed after its job wrote its end (nothing to recover).

generation-output: unchanged
generation-affects: none
