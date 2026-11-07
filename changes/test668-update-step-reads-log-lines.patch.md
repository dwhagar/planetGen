**Headline:** The two update-step tests that checked for a bare "error:" on stderr now read the scripts' log line, so they pass again.

### Changed

- **`test_update_database_step` reads the log line the scripts print (TEST.668).** Since the scripts write every message as a log record (`scripts/logrecord.sh`, for example `2026-11-07 03:32:27 UTC ERROR    command: ...`), an unreachable server and a refused migration no longer print a line that starts with `error:`. `test_an_unreachable_server_stops_the_step` and `test_a_failed_migration_stops_the_step_and_keeps_the_version` now look for that line (time, zone, `ERROR`, component, message). Only the tests change.

generation-output: unchanged
