**Headline:** The engine's own wording and tests now follow the house rules for units, logging and the command list (TEST.795).

### Changed

- The pressure report, the stop-grace line, the database-signal line, the item-file size errors and the light-job time lines now write durations and sizes with the shared formatters (`format_duration_seconds`, `format_bytes`) instead of by hand.
- The template logs, at debug level, why a job's arguments could not be read before telling the supervisor.
- The help test no longer lists the removed `worker` command, and the update-order test knows the `retire_config_options` step.
