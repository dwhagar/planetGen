# Database migrations (Alembic)

Schema changes after v61 are Alembic revisions in `versions/`. A revision's
id is its schema version as four digits (`0062`), so `store.SCHEMA_VERSION`
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
4. Add a marker for the new version to `_VERSION_MARKERS` in
   `planetgen/db/store.py` (a column or index the version introduced).
5. Check in `tests/fixtures/old_schemas/schema_v<previous>.sql.gz`, as
   `tests/fixtures/old_schemas/README.md` describes, so the migration is
   tested against the schema it starts from.

There are no downgrades.
