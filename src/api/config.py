# api/config.py

"""
Configuration for the read-only Flask API.

Intentionally tiny -- one setting today. `MYSQL_CONFIG` is a
`stellarObjects._db.MySQLConfig`, itself built from the same
`PLANETGEN_MYSQL_*` environment variables every other entry point in this
project (`sectorGen.py`, `systemGen.py`, `queryDb.py`) reads, so a WSGI
deployment (see `wsgi.py`) points this API at a specific database without
editing code -- typically a read-only account's credentials (this API
never writes; see `queryDb.py`'s module docstring for the same
read-only-by-grant convention), set via the vhost's `SetEnv` directives or
the `gunicorn` service's environment file.
"""

from stellarObjects._db import MySQLConfig


class Config:
    MYSQL_CONFIG = MySQLConfig()
