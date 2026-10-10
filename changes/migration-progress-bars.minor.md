### Added
- The migration script shows a bar over its steps and, for a revision that works in batches (`alembic_runner.report_progress`), a second bar for that step, both with the time left; when output is not a terminal it prints a line per step and, for a step in batches, at most one line every 30 seconds and at each 10 percent. Migration statements wait up to an hour for a lock (DB.15).
