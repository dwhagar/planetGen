### Added
- The Galaxy Map's opening view (its first tiles and the galaxy stage) is built into the tile cache in the background at the end of every `update.sh` / `update.ps1`, so the first visit after an update no longer waits for the database. `python -m planetgen.cli.warm_map` does the same by hand, e.g. after clearing the tile cache (MAP.134).
