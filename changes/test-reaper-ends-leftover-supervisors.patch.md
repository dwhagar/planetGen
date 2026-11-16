**Headline:** The test runner ends the supervisor a test leaves behind, as it does the workers, so a check run in the next test is no longer served by it.

### Fixed

- `tests/stray_workers.py` ends `planetgen.cli.supervisor` processes on a private test Redis database with the workers (`parse`, `strays`, `reap_database`, `reap_all`). A supervisor stays up for its five idle minutes after the test that started it, and `queue_command_env` ended only workers, so `test_the_status_command_exits_0_for_a_missing_dataset_with_a_fallback_when_told_to` failed on main whenever an earlier test on the same Redis database had started one.
- `tests/test_stray_workers.py` and `tests/test_queue_command_env.py` prove it for supervisors.

generation-output: unchanged
generation-affects: none
