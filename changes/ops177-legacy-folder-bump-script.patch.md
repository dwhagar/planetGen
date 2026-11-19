**Headline:** A legacy/ folder with the one script that brings a pre-8.0 database up to the oldest structure 8.0 accepts (OPS.177, DB.126).

### Added

- `legacy/bump_database_to_8_0.py`: run by hand on a database you keep (back it up first), it runs the old revisions 0062 to 0130 (kept in `legacy/db_revisions/`, as they were before DB.98 replaced them with an empty baseline) and records each version, stopping at v130. Then `update.sh` takes it on. The error that a too-old database gets has named this script since DB.127; now it exists. Nothing in the install, update, migration or program calls it, and it does nothing else.
