**Headline:** The tests of update.sh's database step pass again: they now look for the failure line as the update scripts write it, with its time and severity.

### Changed

- `src/tests/test_update_database_step.py` (TEST.668): `test_an_unreachable_server_stops_the_step` and `test_a_failed_migration_stops_the_step_and_keeps_the_version` expected the failure to begin with `error:`. The database step was not wrong: since OPS.98 and OPS.100 every line a command writes to a pipe carries its time and severity, and the `error:` its message began with is the severity (`2026-11-07 03:25:40 UTC ERROR    command: (2003, "Can't connect to MySQL server ...")`). The step still stops with a non-zero exit, no traceback, the schema version unchanged, and the server's own reason on stderr. The two tests now look for that line (a time, a zone, `ERROR`, the component) and for the reason. They ran into no GRANT problem on a machine with a normal MySQL.

generation-output: unchanged
