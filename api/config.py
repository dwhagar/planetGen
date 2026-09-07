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
    DB_PATH = os.environ.get("PLANETGEN_DB_PATH", DEFAULT_DB_PATH)
