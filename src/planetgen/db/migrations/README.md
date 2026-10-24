# Database migrations (Alembic)

Schema changes after v130 are Alembic revisions in `versions/`. A revision's
id is its schema version as four digits (`0131`), so `store.SCHEMA_VERSION`
is the head revision's number. `schema_migrations` still gets one row per
version (the bookkeeping the rest of the code reads), written by
`planetgen.db.alembic_runner` after each revision runs.

To change the schema:

1. Edit `planetgen/db/schema.sql` to the new shape (fresh databases are
   created from it).
2. Add `versions/00NN_<name>.py` (copy `script.py.mako`'s layout) with
   `down_revision` set to the previous head and an `upgrade()` that makes an
   existing database match `schema.sql`, using `op.*` or `op.execute(...)`.
3. Run `python scripts/generate_db_models.py` to regenerate `planetgen/db/models.py`
   (`tests/test_db_models.py` fails while it disagrees with `schema.sql`).
4. Check in `tests/fixtures/old_schemas/schema_v<previous>.sql.gz`, as
   `tests/fixtures/old_schemas/README.md` describes, so the migration is
   tested against the schema it starts from.

There are no downgrades.

## Slow changes on a big table

A revision that touches a table with millions of rows uses `planetgen/db/migrations/helpers.py` (DB.18)
instead of `op.add_column`, `op.create_index`, `op.alter_column` or `op.batch_alter_table` (never on the big
tables: it copies the table). The helpers name the `ALGORITHM` (and `LOCK=NONE`), so the server refuses when it
cannot do the change online instead of copying silently, and each one does nothing when its work is already
there, so a stopped revision can be run again:

- `add_column_instant`, `drop_column_instant`: a change of the table's definition only.
- `add_index_online`, `drop_index_online`: `ALGORITHM=INPLACE, LOCK=NONE`; minutes on 10^8 rows, no progress.
- `batched_update`: pages of the primary key, a commit per page, progress on the migration's bar.
- `exists_table`, `exists_column`, `exists_index`, `index_columns`, `column_type`: the checks.
- `lock_wait_timeout` (a block), `require_no_long_transactions` (refuses to start an index build behind a
  transaction open for a minute or more) and a retry of a statement that timed out waiting for a lock.

**`sectors` is not instant (DB.84).** Its FULLTEXT index on `name` forces a full table rebuild (`ALGORITHM=COPY`) for any
column change, with room needed for a second copy of the table. Batch every change to `sectors` that a version needs into
one revision, say in that revision's docstring that it rebuilds `sectors`, and keep the index.
