### Added

- Database tests (TEST.6-9, TEST.11-18): SQL portability and strict `sql_mode`, every released schema migrated and compared with a new database, crashed migrations re-run, every column written and read back, boundary values, collation collisions, CHECK constraints, failed sector saves, id-block edges, batch limits and search edge cases.

### Fixed

- Databases from v8-v37 can migrate again (steps re-run safely after a crash, and MariaDB column CHECKs are dropped correctly). Galaxy schema v50 brings migrated databases to exactly the new shape, keeps nebulae and asteroid fields when their sector is deleted, and stores the `--comets`/`--wide-binary` choices with a system's recipe.
- Saving one-off systems at the same time now names them Alpha/Beta like sector saves, and retries on a deadlock.
- Names that differ only by accent ("Vega"/"Véga") are treated as the same name everywhere, and the first holder keeps its own spelling.
- Very large batched saves are split to fit the server's packet limit.
- Search finds words longer than the full-text index's longest token.
- The Systems list works under MariaDB's strict GROUP BY mode.
- A failed sector save no longer leaves objects renamed.
- System and star names are limited to 200 characters, so their planets and moons always fit.
