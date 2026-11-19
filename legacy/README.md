# legacy/

Scripts that are never run by the program, the installer or an update. Each one does one thing, only when a person runs it on purpose (Boss's hard cut at 8.0, 2026-11-12: nothing before 8.0 is accepted by the program).

- `bump_database_to_8_0.py`: brings a database from v61 up to v130, the oldest structure 8.0 accepts, with the old Alembic revisions in `db_revisions/` (the revisions as they were before DB.98 replaced them with an empty baseline). Back the database up first; then run `update.sh`.

Other legacy scripts go here when they are needed.
