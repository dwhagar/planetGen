### Added
- `planetgen check-db` and a "Check the database" section on the Generate page (DB.8): a read-only check of the schema version and models, table health, rows whose parent is gone, ids that would clash, impossible values, sector counts and version keys. It ends with a pass or fail line per check and exits 1 on damage and 2 when a check could not run.
