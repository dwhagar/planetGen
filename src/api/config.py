# api/config.py

"""
Configuration for the read-only Flask API.

Intentionally tiny -- one setting today. `DB_PATH` is read from an
environment variable so a WSGI deployment (see `wsgi.py`) can point at a
specific `.db` file without editing code, falling back to the same
`stellarObjects._db.DEFAULT_DB_PATH` every other entry point in this project
(`sectorGen.py`, `systemGen.py`, `queryDb.py`) already defaults to.
"""

import os

from stellarObjects._db import DEFAULT_DB_PATH


class Config:
    # TODO: DB_PATH assumes a SQLite file path (see PLANETGEN_DB_PATH above
    # and stellarObjects._db.DEFAULT_DB_PATH). Part of the Phase 5 MySQL
    # migration (docs/TODO.md) is replacing this with a connection
    # string/host+credentials pair and a real secrets-handling story (env
    # vars at minimum) instead of a bare path -- mirrors the same
    # SQLite-path assumption in src/queryDb.py's `--db-path`/
    # `open_readonly`. See docs/TODO.md, "Phase 5 -- Web interface".
    DB_PATH = os.environ.get("PLANETGEN_DB_PATH", DEFAULT_DB_PATH)
