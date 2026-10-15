# planetgen/db/__init__.py

"""
The content database (MySQL): writing, reading, migrating and checking it.

`store.py` is the one place generated objects are converted and saved
(connections, schema version, inserts, `load_sector`); `models.py` describes
every table as SQLAlchemy `Table` objects, `alembic_runner.py` applies the
migrations in `db/migrations/` and `login.py` logs in to MySQL. Reads are
`query.py` (lists and display detail), `near.py` (what is near a place),
`corridor.py` and `stones.py` (what lies along a course and the stops in
unknown space), `sector_paths.py` (saved paths through a sector),
`block_stats.py` (the Galaxy Map's block aggregates), `countcache.py` (counts
too big for a request), `species_fit.py` (which species could live where)
and `stats.py` (statistics for the admin page); `render.py` renders a
stored system's text and `edits.py` writes an admin's edits back. `check.py`
is the database check, `fingerprint.py` a galaxy's digest. `docs/database-schema.md`
describes the tables.
"""
