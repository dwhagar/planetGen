**Headline:** A database that has its tables but no recorded schema version is now refused with a clear message instead of being guessed at.

### Removed
- `store.detect_schema_version` and the `_VERSION_MARKERS` table (DB.4): the code no longer works out a version from the shape of the tables, and a new migration no longer needs a marker row.
- The "tables look like schema version N" line of `planetgen.cli.checkdb`; the recorded version, the Alembic revision and the model comparison still cover it.

### Changed
- A database with the galaxy tables but an empty or missing `schema_migrations` raises `SchemaVersionMissingError` from every connection, `migrate --status` and `migrate`, before any change is made. Restore it from a backup, or drop it and run `update.sh`.
