**Headline:** The population pass and the territories pass are now two commands, `planetgen population` and `planetgen territories`.

### Changed
- **`planetgen population` stops at the polities (PERF.125).** It scans the life worlds, dates the civilizations and founds the polities; `--territories-only` is gone. After a first run or a `--rescan`, run `planetgen territories` to work out who owns which system. A `galaxy` run does both as consecutive stages, and the installer's population question runs both.
- **Both passes are work-queue commands** (`run-population`, `run-territories`, one unit each), so they show on the Queue page like any other run.

### Added
- **`planetgen territories`.** Recomputes which polity owns which system from the polities stored, and nothing else.
- **Rules for a new command, held by a test (PERF.120).** The command inventory lists the seven rules every command follows, and `test_command_rules.py` checks the structural ones for every registered command.
